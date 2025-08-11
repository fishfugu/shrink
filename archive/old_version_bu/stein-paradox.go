package main_old

import (
	"errors"
	"fmt"
	"math"
)

const (
	maxIter = 100
	tol     = 1e-12
	eps     = 1e-15
)

// SolveRO computes r_o and the shrink factor s = r_o/d for a given (x_m, y_m, r_m).
// It scales so that d=1, solves for r'_o by Brent’s method, then returns r_o = d * r'_o.
func SolveRO(xm, ym, rm float64) (ro, shrink float64, err error) {
	d := math.Hypot(xm, ym)
	if d <= 0 {
		return 0, 0, errors.New("d = sqrt(x_m^2 + y_m^2) must be > 0")
	}
	if rm <= 0 {
		return 0, 0, errors.New("r_m must be > 0")
	}

	// Scale: d' = 1
	rmS := rm / d

	// Domain for intersection (and thus half-area solvability): r_o' in (1 - r_m', 1 + r_m')
	lo := math.Max(1.0-rmS+1e-15, 1e-15) // guard positivity
	hi := 1.0 + rmS - 1e-15
	if !(lo < hi) {
		return 0, 0, errors.New("invalid bracket; check r_m relative to d")
	}

	// Target: A(r_o') - 0.5 * pi * r_m'^2 = 0
	target := 0.5 * math.Pi * rmS * rmS
	F := func(roS float64) float64 {
		return areaScaled(roS, rmS) - target
	}

	fLo := F(lo)
	fHi := F(hi)
	// At lo, area ~ 0 -> negative; at hi, area ~ pi r_m^2 -> positive
	if math.IsNaN(fLo) || math.IsNaN(fHi) {
		return 0, 0, errors.New("NaN encountered in initial bracket")
	}
	if fLo > 0 && fHi > 0 || fLo < 0 && fHi < 0 {
		// Very rare due to numeric extremes; try a tiny inward nudge
		lo2 := lo + 1e-10
		hi2 := hi - 1e-10
		fLo2, fHi2 := F(lo2), F(hi2)
		if !(fLo2 <= 0 && fHi2 >= 0) {
			return 0, 0, errors.New("failed to bracket the root")
		}
		lo, hi, fLo, fHi = lo2, hi2, fLo2, fHi2
	}

	// Brent–Dekker root finder
	roS, err := brent(F, lo, hi, fLo, fHi, tol)
	if err != nil {
		return 0, 0, err
	}

	// Unscale
	ro = roS * d
	shrink = roS // because shrink = r_o/d = r'_o
	return ro, shrink, nil
}

// areaScaled computes the intersection area A for d=1, using the compact formula:
// A = r_o^2 * alpha + r_m^2 * beta - r_o r_m * sin(gamma),
// where alpha = arccos((1 + r_o^2 - r_m^2)/(2 r_o)), etc.
func areaScaled(ro, rm float64) float64 {
	// Guard domain: |ro - rm| < 1 < ro + rm
	if ro <= 0 || rm <= 0 || !(math.Abs(ro-rm) < 1 && 1 < ro+rm) {
		// Return 0 or full area depending on containment; here 0 is fine for our bracket usage.
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
	// Angles
	alpha := aco((1 + ro*ro - rm*rm) / (2 * ro))
	beta := aco((1 + rm*rm - ro*ro) / (2 * rm))
	gamma := aco((ro*ro + rm*rm - 1) / (2 * ro * rm))

	A := ro*ro*alpha + rm*rm*beta - ro*rm*math.Sin(gamma)
	return A
}

// brent finds a root of F in [a,b] given F(a)=fa <= 0 <= F(b)=fb or vice versa.
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
	xm := 0.0
	var s, p, q float64

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
		xm = 0.5 * (c - b)

		if math.Abs(xm) <= tol1 || fb == 0 {
			return b, nil
		}

		if math.Abs(e) >= tol1 && math.Abs(fa) > math.Abs(fb) {
			// Attempt inverse quadratic interpolation or secant
			if a != c && fa != fc {
				// Inverse quadratic interpolation
				s = fb / fc
				p = s * (2 * xm) // placeholder, replaced below
				// Using quadratic through (a,fa), (b,fb), (c,fc):
				// p/q form from Brent’s method
				q = (fa/fc - 1)
				r := (fb/fc - 1)
				t := (fa/fb - 1)
				p = r * (2*xm*q*(q-t) - (b-a)*(t-1))
				q = (q - 1) * (r - 1) * (t - 1)
			} else {
				// Secant
				s = fb / fa
				p = 2 * xm * s
				q = 1 - s
			}
			if p > 0 {
				q = -q
			}
			p = math.Abs(p)

			// Acceptance criteria
			if 2*p < math.Min(3*xm*q-math.Abs(tol1*q), math.Abs(e*q)) {
				e = d
				d = p / q
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

func main() {
	// Example: measurement at (4,4), rm = 1.0
	xm, ym, rm := 4.0, 4.0, 1.0
	ro, shrink, err := SolveRO(xm, ym, rm)
	if err != nil {
		panic(err)
	}
	fmt.Printf("r_o = %.12f\n", ro)
	fmt.Printf("shrink factor (r_o/d) = %.12f\n", shrink)
}
