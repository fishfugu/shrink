package main

import (
	"encoding/csv"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"
)

// go run likelihood-weighted-shrinkage-factor.go																	# Rayleigh default (sigma=0.25)
// go run likelihood-weighted-shrinkage-factor.go -sigma=0.3 -rmax=0.8     											# Rayleigh with custom params
// go run likelihood-weighted-shrinkage-factor.go -mode=empirical -csv=data_r_theta.csv -csv-out=out_rtheta.csv		# Empirical XY -> scaled CSV output
// go run likelihood-weighted-shrinkage-factor.go -mode=empirical -csv=data_xy.csv -csv-out=out_xy.csv				# Empirical, xy -> scaled output

const (
	maxIter = 100
	tol     = 1e-12
	eps     = 1e-15
)

// ---------- Flags ----------

var (
	mode    = flag.String("mode", "rayleigh", "weighting mode: 'empirical' or 'rayleigh'")
	csvIn   = flag.String("csv", "", "CSV file: either columns x_m,y_m or a single column r_m in (0,1). Required in empirical mode.")
	csvOut  = flag.String("csv-out", "", "Output CSV path (only used in empirical mode). No default; if omitted, no file is written.")
	sigma   = flag.Float64("sigma", 0.25, "Rayleigh sigma (per-axis std) for rayleigh mode (scaled units)")
	rmax    = flag.Float64("rmax", 0.999, "Upper integration limit for r_m in rayleigh mode")
	nGrid   = flag.Int("n", 1000, "Number of grid points for rayleigh mode")
	theta   = flag.Float64("theta", 45.0, "Angle in degrees for reporting/constructing M=(d cosθ, d sinθ) if only r_m is provided")
	println = flag.Bool("printrows", true, "Print per-row results (may be verbose for large n)")
)

// ---------- Core geometry: d(r_m) with r_o=1 ----------

func SolveD(rm float64, theta *float64) (d float64, Mx, My float64, err error) {
	if rm == 0.0 {
		if theta != nil {
			return 1.0, math.Cos(*theta), math.Sin(*theta), nil
		}
		return 1.0, 0, 0, nil
	}
	if rm < 0.0 || rm >= 1.0 {
		return 0, 0, 0, fmt.Errorf("require 0 < r_m < 1 for the half-area split (r_o=1). r_m=%.12f", rm)
	}
	lo := 1.0 - rm + 1e-15
	hi := 1.0 + rm - 1e-15

	target := 0.5 * math.Pi * rm * rm
	F := func(x float64) float64 { return areaRO1(x, rm) - target }

	flo, fhi := F(lo), F(hi)
	if math.IsNaN(flo) || math.IsNaN(fhi) {
		return 0, 0, 0, errors.New("NaN in initial bracket")
	}
	if !(flo > 0 && fhi < 0) {
		lo2, hi2 := lo+1e-10, hi-1e-10
		flo2, fhi2 := F(lo2), F(hi2)
		if !(flo2 > 0 && fhi2 < 0) {
			return 0, 0, 0, errors.New("failed to bracket root; check r_m")
		}
		lo, hi, flo, fhi = lo2, hi2, flo2, fhi2
	}

	root, err := brent(F, lo, hi, flo, fhi, tol)
	if err != nil {
		return 0, 0, 0, err
	}
	d = root
	if theta != nil {
		ct, st := math.Cos(*theta), math.Sin(*theta)
		Mx, My = d*ct, d*st
	}
	return d, Mx, My, nil
}

// Intersection area for r_o = 1:
func areaRO1(d, rm float64) float64 {
	if !(math.Abs(1.0-rm) < d && d < 1.0+rm) {
		return 0
	}
	aco := func(x float64) float64 {
		if x < -1 {
			x = -1
		} else if x > 1 {
			x = 1
		}
		return math.Acos(x)
	}
	alpha := aco((d*d + 1 - rm*rm) / (2 * d))
	beta := aco((d*d + rm*rm - 1) / (2 * d * rm))
	gamma := aco((1 + rm*rm - d*d) / (2 * rm))
	return alpha + rm*rm*beta - rm*math.Sin(gamma)
}

// Brent–Dekker root finder on [a,b].
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

// ---------- Rayleigh weights (if modelling) ----------

func rayleighPDF(r, sigma float64) float64 {
	if r <= 0 {
		return 0
	}
	return (r / (sigma * sigma)) * math.Exp(-(r*r)/(2*sigma*sigma))
}

// ---------- CSV readers ----------

type inKind int

const (
	kindUnknown inKind = iota
	kindRadii
	kindXY
)

type inputData struct {
	kind inKind
	rms  []float64
	xs   []float64
	ys   []float64
}

// readMixedCSV: detects header. If it finds x_m,y_m (or x,y), returns XY; else tries r_m column; else first column as radii.
func readMixedCSV(path string) (*inputData, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.FieldsPerRecord = -1

	var (
		first      = true
		hasHeader  = false
		colX, colY = -1, -1
		colR       = -1
		xs, ys     []float64
		rms        []float64
	)

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if len(rec) == 0 {
			continue
		}
		if first {
			first = false
			// try header
			lc := make([]string, len(rec))
			for i := range rec {
				lc[i] = strings.ToLower(strings.TrimSpace(rec[i]))
			}
			for i, c := range lc {
				switch c {
				case "x_m", "x":
					colX = i
				case "y_m", "y":
					colY = i
				case "r_m", "r":
					if colR == -1 {
						colR = i
					}
				}
			}
			if colX >= 0 && colY >= 0 {
				hasHeader = true
				continue
			}
			if colR >= 0 {
				hasHeader = true
				continue
			}
			// no header—fall through and try to parse this first row as data
		}

		// If we have XY columns
		if colX >= 0 && colY >= 0 {
			if colX >= len(rec) || colY >= len(rec) {
				continue
			}
			x, err1 := strconv.ParseFloat(strings.TrimSpace(rec[colX]), 64)
			y, err2 := strconv.ParseFloat(strings.TrimSpace(rec[colY]), 64)
			if err1 == nil && err2 == nil {
				xs = append(xs, x)
				ys = append(ys, y)
			}
			continue
		}

		// If we have an R column
		if colR >= 0 {
			if colR >= len(rec) {
				continue
			}
			val := strings.TrimSpace(rec[colR])
			if v, err := strconv.ParseFloat(val, 64); err == nil {
				rms = append(rms, v)
			}
			continue
		}

		// No header: try interpret as XY (first two cols) or single radii col
		if len(rec) >= 2 {
			x, err1 := strconv.ParseFloat(strings.TrimSpace(rec[0]), 64)
			y, err2 := strconv.ParseFloat(strings.TrimSpace(rec[1]), 64)
			if err1 == nil && err2 == nil {
				xs = append(xs, x)
				ys = append(ys, y)
				continue
			}
		}
		// else try first col as r_m
		if v, err := strconv.ParseFloat(strings.TrimSpace(rec[0]), 64); err == nil {
			rms = append(rms, v)
		}
	}

	out := &inputData{}
	if len(xs) > 0 && len(ys) == len(xs) {
		out.kind = kindXY
		out.xs = xs
		out.ys = ys
		return out, nil
	}
	if len(rms) > 0 {
		out.kind = kindRadii
		out.rms = rms
		return out, nil
	}
	if hasHeader {
		return nil, errors.New("no valid rows under detected header")
	}
	return nil, errors.New("could not interpret CSV (need x_m,y_m or r_m)")
}

// ---------- Weighted shrinkage calculators ----------

func empiricalFromXY(xs, ys []float64, printRows bool, csvOutPath string) (float64, error) {
	if len(xs) == 0 || len(ys) != len(xs) {
		return 0, errors.New("empty or mismatched XY data")
	}
	var sum float64
	var n int

	var out *csv.Writer
	var outFile *os.File
	if csvOutPath != "" {
		f, err := os.Create(csvOutPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		outFile = f
		out = csv.NewWriter(f)
		defer out.Flush()
		out.Write([]string{"x_m", "y_m", "r_m", "d", "theta_rad", "Mx", "My", "x_scaled", "y_scaled"})
	}

	for i := range xs {
		x, y := xs[i], ys[i]
		rm := math.Hypot(x, y)
		if !(rm > 0 && rm < 1) {
			// skip outside domain
			continue
		}
		th := math.Atan2(y, x)
		d, Mx, My, err := SolveD(rm, &th)
		if err != nil {
			return 0, err
		}
		if printRows {
			fmt.Printf("r_m = %.6f -> d (scaling factor) = %.12f -> M = (%.12f, %.12f)\n", rm, d, Mx, My)
		}
		sum += d
		n++

		if out != nil {
			xScaled := d * x
			yScaled := d * y
			out.Write([]string{
				fmt.Sprintf("%.12f", x),
				fmt.Sprintf("%.12f", y),
				fmt.Sprintf("%.12f", rm),
				fmt.Sprintf("%.12f", d),
				fmt.Sprintf("%.12f", th),
				fmt.Sprintf("%.12f", Mx),
				fmt.Sprintf("%.12f", My),
				fmt.Sprintf("%.12f", xScaled),
				fmt.Sprintf("%.12f", yScaled),
			})
		}
	}

	if n == 0 {
		if outFile != nil {
			os.Remove(outFile.Name())
		}
		return 0, errors.New("no valid rows in (0,1) after filtering XY")
	}
	return sum / float64(n), nil
}

func empiricalFromRadii(rms []float64, thetaRad float64, printRows bool, csvOutPath string) (float64, error) {
	if len(rms) == 0 {
		return 0, errors.New("no radii provided")
	}
	var sum float64
	var n int

	var out *csv.Writer
	var outFile *os.File
	if csvOutPath != "" {
		f, err := os.Create(csvOutPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		outFile = f
		out = csv.NewWriter(f)
		defer out.Flush()
		// We don't have original x_m,y_m—emit constructed scaled coordinates on the θ ray.
		out.Write([]string{"r_m", "d", "theta_rad", "Mx", "My", "x_scaled", "y_scaled"})
	}

	for _, rm := range rms {
		if !(rm > 0 && rm < 1) {
			continue
		}
		d, Mx, My, err := SolveD(rm, &thetaRad)
		if err != nil {
			return 0, err
		}
		if printRows {
			fmt.Printf("r_m = %.6f -> d (scaling factor) = %.12f -> M = (%.12f, %.12f)\n", rm, d, Mx, My)
		}
		sum += d
		n++

		if out != nil {
			// Without original (x_m,y_m), we can only output the scaled point on the θ ray
			xScaled := d * math.Cos(thetaRad)
			yScaled := d * math.Sin(thetaRad)
			out.Write([]string{
				fmt.Sprintf("%.12f", rm),
				fmt.Sprintf("%.12f", d),
				fmt.Sprintf("%.12f", thetaRad),
				fmt.Sprintf("%.12f", Mx),
				fmt.Sprintf("%.12f", My),
				fmt.Sprintf("%.12f", xScaled),
				fmt.Sprintf("%.12f", yScaled),
			})
		}
	}

	if n == 0 {
		if outFile != nil {
			os.Remove(outFile.Name())
		}
		return 0, errors.New("no r_m values in (0,1) after filtering")
	}
	return sum / float64(n), nil
}

// Rayleigh expected shrink: quadrature over (0, rmax] with weights f_R(r; sigma).
func rayleighShrinkage(sigma, rmax float64, n int, thetaRad float64, printRows bool) (float64, error) {
	if sigma <= 0 {
		return 0, errors.New("sigma must be > 0")
	}
	if rmax <= 0 || rmax >= 1 {
		return 0, errors.New("rmax must be in (0,1)")
	}
	if n < 2 {
		n = 2
	}
	dr := rmax / float64(n-1)
	var num, den float64
	for i := 0; i < n; i++ {
		rm := dr * float64(i)
		if rm == 0 {
			continue
		}
		w := rayleighPDF(rm, sigma)
		d, Mx, My, err := SolveD(rm, &thetaRad)
		if err != nil {
			return 0, err
		}
		if printRows {
			fmt.Printf("r_m = %.6f -> d (scaling factor) = %.12f -> M = (%.12f, %.12f)\n", rm, d, Mx, My)
		}
		num += w * d
		den += w
	}
	if den == 0 {
		return 0, errors.New("zero total Rayleigh weight (check sigma/rmax)")
	}
	return num / den, nil
}

// ---------- Main ----------

func main() {
	flag.Parse()
	thetaRad := (*theta) * math.Pi / 180.0

	switch strings.ToLower(*mode) {
	case "empirical":
		if *csvIn == "" {
			fmt.Fprintln(os.Stderr, "empirical mode requires -csv=path/to/input.csv")
			os.Exit(1)
		}
		in, err := readMixedCSV(*csvIn)
		if err != nil {
			fmt.Fprintf(os.Stderr, "read error: %v\n", err)
			os.Exit(1)
		}

		switch in.kind {
		case kindXY:
			tilde, err := empiricalFromXY(in.xs, in.ys, *println, *csvOut)
			if err != nil {
				fmt.Fprintf(os.Stderr, "empirical XY error: %v\n", err)
				os.Exit(1)
			}
			fmt.Printf("\nSingle shrink factor (empirical mean of d(r_m)): %.12f\n", tilde)
			if *csvOut != "" {
				fmt.Printf("Wrote %s\n", *csvOut)
			}
		case kindRadii:
			tilde, err := empiricalFromRadii(in.rms, thetaRad, *println, *csvOut)
			if err != nil {
				fmt.Fprintf(os.Stderr, "empirical radii error: %v\n", err)
				os.Exit(1)
			}
			fmt.Printf("\nSingle shrink factor (empirical mean of d(r_m)): %.12f\n", tilde)
			if *csvOut != "" {
				fmt.Printf("Wrote %s\n", *csvOut)
			}
		default:
			fmt.Fprintln(os.Stderr, "could not determine input kind (need x_m,y_m or r_m)")
			os.Exit(1)
		}

	case "rayleigh":
		tilde, err := rayleighShrinkage(*sigma, *rmax, *nGrid, thetaRad, *println)
		if err != nil {
			fmt.Fprintf(os.Stderr, "rayleigh error: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("\nSingle shrink factor (Rayleigh PDF-weighted, σ=%.3f, rmax=%.3f): %.12f\n", *sigma, *rmax, tilde)

	default:
		fmt.Fprintln(os.Stderr, "unknown -mode; use 'empirical' or 'rayleigh'")
		os.Exit(1)
	}
}
