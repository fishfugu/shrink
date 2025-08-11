package main

import (
	"encoding/csv"
	"errors"
	"fmt"
	"math"
	"os"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

const (
	maxIter = 100
	tol     = 1e-12
	eps     = 1e-15
)

// ---------- Solver (same logic you had) ----------

// SolveD returns d given r_m (with r_o = 1) by solving A(d) = 0.5 * pi * r_m^2.
// If theta is provided (radians), it also returns the inferred centre M = d * (cos theta, sin theta).
func SolveD(rm float64, theta *float64) (d float64, Mx, My float64, err error) {
	if rm == 0.0 {
		if theta != nil {
			return 1.0, math.Cos(*theta), math.Sin(*theta), nil
		}
		return 1.0, 0, 0, nil
	}
	if rm < 0.0 || rm >= 1.0 {
		return 0, 0, 0, fmt.Errorf("require 0 < r_m < 1 for a clean half-area split w.r.t r_o=1. r_m = %.12f", rm)
	}
	// Bracket where overlap goes from full (d <= 1 - r_m) to zero (d >= 1 + r_m)
	lo := 1.0 - rm + 1e-15 // just inside intersection regime
	hi := 1.0 + rm - 1e-15 // just inside intersection regime

	target := 0.5 * math.Pi * rm * rm
	F := func(x float64) float64 { return areaRO1(x, rm) - target }

	flo, fhi := F(lo), F(hi)
	if math.IsNaN(flo) || math.IsNaN(fhi) {
		return 0, 0, 0, errors.New("NaN in initial bracket")
	}
	// Expect flo > 0 (area near full) and fhi < 0 (area near 0). If not, nudge.
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

// Intersection area with r_o = 1, r_m = rm, centres at distance d:
// A(d) = alpha + rm^2 * beta - rm * sin(gamma)
func areaRO1(d, rm float64) float64 {
	// Domain guard for intersection
	if !(math.Abs(1.0-rm) < d && d < 1.0+rm) {
		// Outside intersection: return 0 or full as you prefer. Here we use 0.
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

// Brent–Dekker root finder on [a,b] with F(a)=fa>0 and F(b)=fb<0 (or vice versa).
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
			// Inverse quadratic interpolation or secant
			s := fb / fa
			if a != c && fa != fc {
				// Inverse quadratic
				r := fb / fc
				t := fa / fc
				p = s * (2*xm*r*(r-t) - (b-a)*(t-1))
				q = (r - 1) * (t - 1) * (s - 1)
			} else {
				// Secant
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

// ---------- Visualisation + CSV ----------

func main() {
	// measurement at theta = 45 degrees
	th := math.Pi / 4

	// rm = 0 should give us the original measurement exactly
	dOrig, mxOrig, myOrig, err := SolveD(0, &th)
	if err != nil {
		panic(err)
	}
	fmt.Println()
	fmt.Printf("Original M should be at a distance of 1. d = %.12f\n", dOrig)
	if dOrig != 1.0 {
		panic("d should be 1. If it's not, could be a rounding error?")
	}
	fmt.Printf("Original measurement = (%.12f, %.12f)\n", mxOrig, myOrig)
	fmt.Println()

	// Prepare CSV writers
	drmFile, err := os.Create("drm.csv")
	if err != nil {
		panic(err)
	}
	defer drmFile.Close()
	drmW := csv.NewWriter(drmFile)
	defer drmW.Flush()
	_ = drmW.Write([]string{"r_m", "d"}) // header

	locusFile, err := os.Create("locus.csv")
	if err != nil {
		panic(err)
	}
	defer locusFile.Close()
	locusW := csv.NewWriter(locusFile)
	defer locusW.Flush()
	_ = locusW.Write([]string{"r_m", "theta_rad", "mx", "my"}) // header

	// For plotting
	n := 1000
	pts := make(plotter.XYs, 0, n)
	ptsRay := make(plotter.XYs, 0, n)

	// Print loop in same format AND write CSV rows
	for i := 0; i < n; i++ {
		rm := float64(i) / 1000.0 // 0.000 to 0.999 step 0.001
		d, mx, my, err := SolveD(rm, &th)
		if err != nil {
			panic(err)
		}
		fmt.Printf("r_m = %.6f -> d (scaling factor) = %.12f -> M = (%.12f, %.12f)\n", rm, d, mx, my)

		// CSV rows
		_ = drmW.Write([]string{
			fmt.Sprintf("%.6f", rm),
			fmt.Sprintf("%.12f", d),
		})
		_ = locusW.Write([]string{
			fmt.Sprintf("%.6f", rm),
			fmt.Sprintf("%.12f", th),
			fmt.Sprintf("%.12f", mx),
			fmt.Sprintf("%.12f", my),
		})

		// Collect for plots
		pts = append(pts, plotter.XY{X: rm, Y: d})
		ptsRay = append(ptsRay, plotter.XY{X: mx, Y: my})
	}
	drmW.Flush()
	locusW.Flush()

	fmt.Println("-----")

	// Plot shrink factor d(r_m)
	p1 := plot.New()
	p1.Title.Text = "Shrink Factor d as a Function of r_m (r_o = 1)"
	p1.X.Label.Text = "r_m"
	p1.Y.Label.Text = "d (scaled distance of M)"
	l1, err := plotter.NewLine(pts)
	if err != nil {
		panic(err)
	}
	plotutil.AddLines(p1, l1)
	if err := p1.Save(6*vg.Inch, 4*vg.Inch, "drm.png"); err != nil {
		panic(err)
	}
	fmt.Println("Wrote drm.png and drm.csv")

	// Plot locus of M for θ = 45°
	p2 := plot.New()
	p2.Title.Text = "Locus of Inferred Centre M for θ = 45°"
	p2.X.Label.Text = "Mx"
	p2.Y.Label.Text = "My"
	l2, err := plotter.NewLine(ptsRay)
	if err != nil {
		panic(err)
	}
	plotutil.AddLines(p2, l2)
	if err := p2.Save(6*vg.Inch, 6*vg.Inch, "locus.png"); err != nil {
		panic(err)
	}
	fmt.Println("Wrote locus.png and locus.csv")
}
