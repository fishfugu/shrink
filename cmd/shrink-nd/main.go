package main

import (
	"encoding/csv"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mathext"
)

/*
Build & run:

  go mod init shrink_nd
  go get gonum.org/v1/gonum@latest
  go run main.go -mode=chi -n=3 -sigma=0.3 -rmax=0.95 -nGrid=2000

Empirical (CSV of nD points, each row is x1,...,xn):

  go run cmd/shrink-nd/main.go -mode=empirical -csv=examples/nd/data_3d.csv -csv-out=examples/nd/out_scaled.csv -autoscale

Flags explained below in main().

go run cmd/shrink-nd/main.go \
  -gen=examples/nd/data_3d.csv -gen-n=3 -gen-rows=200 -gen-sigma=0.25 -gen-seed=42 \
  -gen-mean=0.3,0.1,-0.2

go run cmd/shrink-nd/main.go -mode=empirical -csv=examples/nd/data_3d.csv -csv-out=examples/nd/scaled_3d.csv -autoscale

*/

// ---------- Flags ----------

var (
	mode      = flag.String("mode", "chi", "mode: 'empirical' (CSV points) or 'chi' (model weighting)")
	csvIn     = flag.String("csv", "", "input CSV for empirical mode; rows are x1,...,xn (numeric)")
	csvOut    = flag.String("csv-out", "", "output CSV path (empirical): original + r + d + scaled coords")
	nFlag     = flag.Int("n", 0, "dimension for chi mode (required in chi mode). In empirical mode, n is inferred")
	sigma     = flag.Float64("sigma", 0.25, "per-axis std σ for Chi(n,σ) model (chi mode)")
	rmax      = flag.Float64("rmax", 0.999, "upper integration limit for r in chi mode (0<rmax<1)")
	nGrid     = flag.Int("nGrid", 1000, "grid size for chi-mode integration")
	printRows = flag.Bool("printrows", true, "print per-row diagnostics")
	autoscale = flag.Bool("autoscale", false, "empirical: rescale all radii by max radius so that max r < 1")
	applyMode = flag.String("apply", "perpoint", "empirical apply mode: perpoint | global-js | global-chi")

	genOut   = flag.String("gen", "", "generate n-D Gaussian points to this CSV (instead of processing)")
	genN     = flag.Int("gen-n", 2, "dimensions for generated points")
	genRows  = flag.Int("gen-rows", 200, "number of points to generate")
	genSigma = flag.Float64("gen-sigma", 1.0, "per-axis std σ for generated points")
	genSeed  = flag.Int64("gen-seed", 42, "RNG seed (fixed default for reproducibility)")
	genMean  = flag.String("gen-mean", "", "comma-separated mean vector for generated points")
)

// ---------- n-ball volume & spherical caps ----------

// Vn returns the volume of an n-ball of radius R: V_n(R) = π^{n/2} / Γ(n/2+1) * R^n
func Vn(R float64, n int) float64 {
	lg, _ := math.Lgamma(float64(n)/2 + 1)
	return math.Exp((float64(n)/2)*math.Log(math.Pi)-lg) * math.Pow(R, float64(n))
}

// IncBetaReg is a thin wrapper around Gonum's regularized incomplete beta I_z(a,b).
func IncBetaReg(z, a, b float64) float64 {
	// Clamp to [0,1] for numerical safety.
	if z <= 0 {
		return 0
	}
	if z >= 1 {
		return 1
	}
	return mathext.RegIncBeta(a, b, z)
}

// CapVolume_n: volume of an n-D spherical cap of an n-ball (radius R, cap height h in [0,2R]).
// Uses the small-cap formula for h<=R and complementary reflection for h>R.
func CapVolume_n(R, h float64, n int) float64 {
	if h <= 0 {
		return 0
	}
	if h >= 2*R {
		return Vn(R, n)
	}
	smallCap := func(R, h float64, n int) float64 {
		// z = (2Rh - h^2)/R^2 in [0,1]
		z := (2*R*h - h*h) / (R * R)
		if z <= 0 {
			return 0
		}
		if z >= 1 {
			return Vn(R, n)
		}
		a := 0.5 * float64(n+1)
		b := 0.5
		return 0.5 * Vn(R, n) * IncBetaReg(z, a, b)
	}
	if h <= R {
		return smallCap(R, h, n)
	}
	// Large cap (h>R): volume = full - complementary small cap of height (2R - h)
	return Vn(R, n) - smallCap(R, 2*R-h, n)
}

// OverlapVolume_n: intersection volume of two n-balls with radii 1 and r, centers distance d apart.
func OverlapVolume_n(d, r float64, n int) float64 {
	if d >= 1+r {
		return 0
	}
	if d <= math.Abs(1-r) {
		// Smaller ball is wholly contained.
		return Vn(math.Min(1, r), n)
	}
	// Cutting hyperplane position from center of radius-1 ball
	x := (d*d - r*r + 1) / (2 * d)
	h1 := 1 - x
	h2 := r - (d*d+r*r-1)/(2*d)
	return CapVolume_n(1, h1, n) + CapVolume_n(r, h2, n)
}

// ---------- Root solving (Brent–Dekker) ----------

const (
	maxIter = 100
	tol     = 1e-12
	eps     = 1e-15
)

// SolveD_n finds d in (|1-r|,1+r) such that overlap volume is half of the r-ball volume.
func SolveD_n(n int, r float64) (float64, error) {
	if n < 2 {
		return 0, fmt.Errorf("n must be >= 2; got %d", n)
	}
	if !(r > 0 && r < 1) {
		return 0, fmt.Errorf("r must be in (0,1); got %g", r)
	}
	lo := 1 - r + 1e-15
	hi := 1 + r - 1e-15
	target := 0.5 * Vn(r, n)
	F := func(d float64) float64 { return OverlapVolume_n(d, r, n) - target }

	flo, fhi := F(lo), F(hi)
	if !(flo > 0 && fhi < 0) {
		lo2, hi2 := lo+1e-12, hi-1e-12
		flo, fhi = F(lo2), F(hi2)
		if !(flo > 0 && fhi < 0) {
			return 0, fmt.Errorf("failed to bracket root; n=%d r=%.9f (F(lo)=%.3e, F(hi)=%.3e)", n, r, flo, fhi)
		}
		lo, hi = lo2, hi2
	}
	root, err := brent(F, lo, hi, flo, fhi, tol)
	if err != nil {
		return 0, err
	}
	return root, nil
}

func brent(F func(float64) float64, a, b, fa, fb, tol float64) (float64, error) {
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if fa*fb > 0 {
		return 0, errors.New("brent: root not bracketed")
	}
	c, fc := a, fa
	d, e := b-a, b-a
	for i := 0; i < maxIter; i++ {
		if fb*fc > 0 {
			c, fc = a, fa
			d, e = b-a, b-a
		}
		if math.Abs(fc) < math.Abs(fb) {
			a, b, c = b, c, b
			fa, fb, fc = fb, fc, fb
		}
		tol1 := 2*eps*math.Abs(b) + 0.5*tol
		xm := 0.5 * (c - b)
		if math.Abs(xm) <= tol1 || fb == 0 {
			return b, nil
		}
		var p, q float64
		if math.Abs(e) >= tol1 && math.Abs(fa) > math.Abs(fb) {
			s := fb / fa
			if a != c && fa != fc {
				r := fb / fc
				t := fa / fc
				p = s * (2*xm*r*(r-t) - (b-a)*(t-1))
				q = (r - 1) * (t - 1) * (s - 1)
			} else {
				p = 2 * xm * s
				q = 1 - s
			}
			if p > 0 {
				q = -q
			}
			p = math.Abs(p)
			if 2*p < math.Min(3*xm*q-math.Abs(tol1*q), math.Abs(e*q)) {
				e, d = d, p/q
			} else {
				d = xm
				e = d
			}
		} else {
			d = xm
			e = d
		}
		a, fa = b, fb
		if math.Abs(d) > tol1 {
			b += d
		} else {
			b += math.Copysign(tol1, xm)
		}
		fb = F(b)
	}
	return b, errors.New("brent: maximum iterations exceeded")
}

// ---------- Chi(n,σ) PDF for model weighting ----------

// ChiPDF returns the χ_n PDF for radius r with per-axis std σ (R = ||X||, X~N(0,σ^2 I_n)).
func ChiPDF(n int, sigma, r float64) float64 {
	if r <= 0 || sigma <= 0 {
		return 0
	}
	// f(r) = 1/(2^{n/2-1} Γ(n/2) σ^n) * r^{n-1} * exp(-r^2/(2σ^2))
	lnGamma, _ := math.Lgamma(float64(n) / 2)
	logC := -((float64(n)/2 - 1) * math.Log(2)) - lnGamma - float64(n)*math.Log(sigma)
	return math.Exp(logC + (float64(n)-1)*math.Log(r) - (r*r)/(2*sigma*sigma))
}

// ---------- CSV I/O (empirical) ----------

func readNDCSV(path string) (data [][]float64, n int, err error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, 0, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.FieldsPerRecord = -1

	var first = true
	var cols []int // numeric columns indices
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, 0, err
		}
		if len(rec) == 0 {
			continue
		}
		if first {
			first = false
			// Try to auto-detect numeric columns by attempting to parse all cells.
			cand := make([]int, 0, len(rec))
			okAll := true
			for j := range rec {
				if _, err := strconv.ParseFloat(strings.TrimSpace(rec[j]), 64); err != nil {
					okAll = false
					break
				}
				cand = append(cand, j)
			}
			if okAll && len(cand) >= 2 {
				cols = cand
			} else {
				// Treat first row as header; the next non-empty line must be numeric
				continue
			}
		}
		if len(cols) == 0 {
			// Establish numeric cols from first numeric row
			cand := make([]int, 0, len(rec))
			okAll := true
			for j := range rec {
				if _, err := strconv.ParseFloat(strings.TrimSpace(rec[j]), 64); err != nil {
					okAll = false
					break
				}
				cand = append(cand, j)
			}
			if !okAll || len(cand) < 2 {
				continue
			}
			cols = cand
		}
		vec := make([]float64, 0, len(cols))
		ok := true
		for _, j := range cols {
			x, err := strconv.ParseFloat(strings.TrimSpace(rec[j]), 64)
			if err != nil {
				ok = false
				break
			}
			vec = append(vec, x)
		}
		if ok {
			data = append(data, vec)
		}
	}
	if len(data) == 0 {
		return nil, 0, errors.New("no numeric rows found")
	}
	return data, len(data[0]), nil
}

func norms(data [][]float64) []float64 {
	rs := make([]float64, len(data))
	for i, v := range data {
		rs[i] = l2(v)
	}
	return rs
}

func l2(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += x * x
	}
	return math.Sqrt(s)
}

// ---------- Drivers ----------
func empiricalDriver(pathIn, pathOut string, doPrint, doAutoscale bool) error {
	points, n, err := readNDCSV(pathIn)
	if err != nil {
		return err
	}
	if n < 2 {
		return fmt.Errorf("need at least 2 dimensions; inferred n=%d", n)
	}
	rs := norms(points)

	// Optional rescale to keep r in (0,1) (scales both points and radii)
	scale := 1.0
	if doAutoscale {
		var rmax float64
		for _, r := range rs {
			if r > rmax {
				rmax = r
			}
		}
		if rmax <= 0 {
			return errors.New("all-zero data after parsing")
		}
		scale = 0.999 / rmax
		for i := range points {
			for j := range points[i] {
				points[i][j] *= scale
			}
		}
		for i := range rs {
			rs[i] *= scale
		}
	}

	// Output CSV writer (if requested)
	var out *csv.Writer
	var outFile *os.File
	if pathOut != "" {
		f, err := os.Create(pathOut)
		if err != nil {
			return err
		}
		defer f.Close()
		outFile = f
		out = csv.NewWriter(f)
		defer out.Flush()

		// header: x1..xn, r, d, x1_scaled..xn_scaled
		hdr := make([]string, 0, 2*n+2)
		for j := 1; j <= n; j++ {
			hdr = append(hdr, fmt.Sprintf("x_%d", j))
		}
		hdr = append(hdr, "r", "d")
		for j := 1; j <= n; j++ {
			hdr = append(hdr, fmt.Sprintf("x_%d_scaled", j))
		}
		out.Write(hdr)
	}

	switch strings.ToLower(*applyMode) {
	case "perpoint":
		var sumD float64
		var used int
		for i, r := range rs {
			if !(r > 0 && r < 1) {
				continue
			}
			d, err := SolveD_n(n, r)
			if err != nil {
				return fmt.Errorf("row %d: %v", i, err)
			}
			if doPrint {
				fmt.Printf("r = %.6f -> d (shrink) = %.12f\n", r, d)
			}
			sumD += d
			used++

			if out != nil {
				row := make([]string, 0, 2*n+2)
				for j := 0; j < n; j++ {
					row = append(row, fmt.Sprintf("%.12f", points[i][j]))
				}
				row = append(row, fmt.Sprintf("%.12f", r), fmt.Sprintf("%.12f", d))
				for j := 0; j < n; j++ {
					row = append(row, fmt.Sprintf("%.12f", d*points[i][j]))
				}
				out.Write(row)
			}
		}
		if used == 0 {
			if outFile != nil {
				os.Remove(outFile.Name())
			}
			return errors.New("no rows with r in (0,1) after filtering/scaling")
		}
		fmt.Printf("\nSingle shrink factor (empirical mean of d(r)): %.12f\n", sumD/float64(used))
		if outFile != nil {
			fmt.Printf("Wrote %s (%d rows)\n", outFile.Name(), used)
		}
		return nil

	case "global-chi":
		// Compute one geometric factor = empirical mean of d(r), then apply uniformly
		var sumD float64
		var used int
		var dlist []float64
		for i, r := range rs {
			if !(r > 0 && r < 1) {
				continue
			}
			d, err := SolveD_n(n, r)
			if err != nil {
				return fmt.Errorf("row %d: %v", i, err)
			}
			dlist = append(dlist, d)
			sumD += d
			used++
		}
		if used == 0 {
			return errors.New("no rows with r in (0,1) after filtering/scaling")
		}
		dGlobal := sumD / float64(used)
		fmt.Printf("Global (empirical-chi) shrink factor: %.12f\n", dGlobal)

		// Print rows (optional): show each r with the *global* factor
		if doPrint {
			for _, r := range rs {
				if r > 0 && r < 1 {
					fmt.Printf("r = %.6f -> d_global (shrink) = %.12f\n", r, dGlobal)
				}
			}
		}
		if out != nil {
			for i, r := range rs {
				if !(r > 0 && r < 1) {
					continue
				}
				row := make([]string, 0, 2*n+2)
				for j := 0; j < n; j++ {
					row = append(row, fmt.Sprintf("%.12f", points[i][j]))
				}
				row = append(row, fmt.Sprintf("%.12f", r), fmt.Sprintf("%.12f", dGlobal))
				for j := 0; j < n; j++ {
					row = append(row, fmt.Sprintf("%.12f", dGlobal*points[i][j]))
				}
				out.Write(row)
			}
			fmt.Printf("Wrote %s (%d rows)\n", outFile.Name(), used)
		}
		return nil

	case "global-js":
		// James–Stein-style global factor, applied uniformly to all rows.
		xbar := meanVec(points, n)
		rbar := norm2(xbar) // ||xbar||
		if rbar == 0 {
			return errors.New("sample mean is zero; global-JS factor undefined (would be 0).")
		}
		s2 := sigma2HatIso(points, n, xbar) // per-axis variance
		N := float64(len(points))
		// JS factor for the mean: s = max(0, 1 - ((p-2)*σ^2/n)/||xbar||^2)
		// Here: p=n, "n" in the formula is sample size; per-axis var of mean is σ^2/N
		// So  ((p-2)*σ^2/N) / ||xbar||^2:
		num := (float64(n) - 2) * (s2 / N)
		den := rbar * rbar
		sJS := 1.0 - num/den
		if sJS < 0 {
			sJS = 0
		}
		fmt.Printf("Global (James–Stein-style) shrink factor: %.12f  [||xbar||=%.6f, sigma^2=%.6f, N=%d]\n",
			sJS, rbar, s2, len(points))

		// Print per-row using the same global factor (optional)
		if doPrint {
			for _, r := range rs {
				if r > 0 && r < 1 {
					fmt.Printf("r = %.6f -> d_global (shrink) = %.12f\n", r, sJS)
				}
			}
		}
		// Output scaled data
		if out != nil {
			used := 0
			for i, r := range rs {
				if !(r > 0 && r < 1) {
					continue
				}
				row := make([]string, 0, 2*n+2)
				for j := 0; j < n; j++ {
					row = append(row, fmt.Sprintf("%.12f", points[i][j]))
				}
				row = append(row, fmt.Sprintf("%.12f", r), fmt.Sprintf("%.12f", sJS))
				for j := 0; j < n; j++ {
					row = append(row, fmt.Sprintf("%.12f", sJS*points[i][j]))
				}
				out.Write(row)
				used++
			}
			fmt.Printf("Wrote %s (%d rows)\n", outFile.Name(), used)
		}
		return nil
	default:
		return fmt.Errorf("unknown -apply mode: %s (use perpoint | global-js | global-chi)", *applyMode)
	}
}

func chiDriver(n int, sigma, rmax float64, nGrid int, doPrint bool) error {
	if n < 2 {
		return fmt.Errorf("-n (dimension) must be >= 2; got %d", n)
	}
	if sigma <= 0 {
		return fmt.Errorf("-sigma must be > 0; got %g", sigma)
	}
	if !(rmax > 0 && rmax < 1) {
		return fmt.Errorf("-rmax must be in (0,1); got %g", rmax)
	}
	if nGrid < 2 {
		nGrid = 2
	}
	dr := rmax / float64(nGrid-1)

	var num, den float64
	for i := 0; i < nGrid; i++ {
		r := dr * float64(i)
		if r <= 0 {
			continue
		}
		w := ChiPDF(n, sigma, r)
		d, err := SolveD_n(n, r)
		if err != nil {
			return err
		}
		if doPrint {
			fmt.Printf("r = %.6f -> d (shrink) = %.12f\n", r, d)
		}
		num += w * d
		den += w
	}
	if den == 0 {
		return errors.New("zero weight in chi integration (check sigma/rmax)")
	}
	fmt.Printf("\nSingle shrink factor (Chi-weighted, n=%d, σ=%.3f, rmax=%.3f): %.12f\n",
		n, sigma, rmax, num/den)
	return nil
}

// helper to generate CSV:
func generateGaussianCSV(path string, n, rows int, sigma float64, seed int64) error {
	if n < 1 || rows < 1 {
		return fmt.Errorf("invalid gen-n=%d or gen-rows=%d", n, rows)
	}
	if sigma <= 0 {
		return fmt.Errorf("sigma must be >0; got %g", sigma)
	}
	rng := rand.New(rand.NewSource(seed))
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	defer w.Flush()

	// header
	hdr := make([]string, n)
	for j := range hdr {
		hdr[j] = fmt.Sprintf("x_%d", j+1)
	}
	w.Write(hdr)

	meanVec, err := parseMeanVector(*genMean, n)
	if err != nil {
		return err
	}

	for i := 0; i < rows; i++ {
		row := make([]string, n)
		for j := 0; j < n; j++ {
			x := meanVec[j] + rng.NormFloat64()*sigma
			row[j] = fmt.Sprintf("%.12f", x)
		}
		w.Write(row)
	}
	fmt.Printf("Generated %d×%d Gaussian points to %s (σ=%.3f, seed=%d)\n", rows, n, path, sigma, seed)
	return nil
}

func parseMeanVector(s string, n int) ([]float64, error) {
	parts := strings.Split(s, ",")
	if len(parts) != n {
		return nil, fmt.Errorf("mean vector length %d != gen-n %d", len(parts), n)
	}
	mean := make([]float64, n)
	for i, p := range parts {
		val, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err != nil {
			return nil, fmt.Errorf("bad mean value at index %d: %v", i, err)
		}
		if math.IsNaN(val) || math.IsInf(val, 0) {
			return nil, fmt.Errorf("invalid mean value at index %d: %v", i, val)
		}
		mean[i] = val
	}
	return mean, nil
}

// mean vector of points (n dims)
func meanVec(points [][]float64, n int) []float64 {
	m := make([]float64, n)
	if len(points) == 0 {
		return m
	}
	for _, v := range points {
		for j := 0; j < n; j++ {
			m[j] += v[j]
		}
	}
	inv := 1.0 / float64(len(points))
	for j := 0; j < n; j++ {
		m[j] *= inv
	}
	return m
}

func norm2(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += x * x
	}
	return math.Sqrt(s)
}

// Unbiased per-axis variance estimator for isotropic Gaussian:
// sigma^2_hat = (1/(n*(N-1))) * sum_i ||x_i - xbar||^2
func sigma2HatIso(points [][]float64, n int, xbar []float64) float64 {
	N := len(points)
	if N <= 1 {
		return 0
	}
	var sum float64
	for _, v := range points {
		for j := 0; j < n; j++ {
			d := v[j] - xbar[j]
			sum += d * d
		}
	}
	return sum / (float64(n) * float64(N-1))
}

// ---------- Main ----------

func main() {
	flag.Parse()

	// If generation requested, do it and exit
	if *genOut != "" {
		if strings.TrimSpace(*genMean) == "" {
			fmt.Fprintln(os.Stderr, "generate mode requires -gen-mean with exactly gen-n comma-separated values, e.g. -gen-mean=0.5,0.5 for gen-n=2")
			os.Exit(1)
		}
		if err := generateGaussianCSV(*genOut, *genN, *genRows, *genSigma, *genSeed); err != nil {
			fmt.Fprintf(os.Stderr, "generate error: %v\n", err)
			os.Exit(1)
		}
		return
	}

	switch strings.ToLower(*mode) {
	case "empirical":
		if *csvIn == "" {
			fmt.Fprintln(os.Stderr, "empirical mode requires -csv=path/to/input.csv (rows are x1,...,xn)")
			os.Exit(1)
		}
		if err := empiricalDriver(*csvIn, *csvOut, *printRows, *autoscale); err != nil {
			fmt.Fprintf(os.Stderr, "empirical error: %v\n", err)
			os.Exit(1)
		}
	case "chi":
		if *nFlag == 0 {
			fmt.Fprintln(os.Stderr, "chi mode requires -n (dimension)")
			os.Exit(1)
		}
		if err := chiDriver(*nFlag, *sigma, *rmax, *nGrid, *printRows); err != nil {
			fmt.Fprintf(os.Stderr, "chi error: %v\n", err)
			os.Exit(1)
		}
	default:
		fmt.Fprintln(os.Stderr, "unknown -mode; use 'empirical' or 'chi'")
		os.Exit(1)
	}
}
