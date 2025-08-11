# Stein Paradox Experiments — Geometric & James–Stein Shrinkage

This repository contains two related CLIs exploring shrinkage estimators:

- **[`shrink-nd`](./cmd/shrink-nd)** — the main tool.  
  An $n$-dimensional geometric shrinkage rule based on **half-overlap** of $n$-balls,  
  plus two global variants (James–Stein-style and χ-weighted).
- **[`lwsf`](./cmd/lwsf)** — the original 2-D prototype (*likelihood-weighted shrinkage factor*).

A generator for reproducible Gaussian $n$-D datasets is included.

---

## 📂 Directory layout

```

.
├── cmd/
│   ├── lwsf/         # 2D prototype CLI
│   └── shrink-nd/    # nD generalization CLI
├── examples/
│   ├── 2d/           # sample CSVs for lwsf
│   └── nd/           # sample CSVs and outputs for shrink-nd
├── docs/             # figures, background notes
├── archive/          # older code kept for reference
└── README.md

````

---

## 🔧 Install

```bash
go mod tidy
# The repo is a single Go module; commands live under ./cmd/...
````

> **Dependencies:**
> [Gonum](https://github.com/gonum/gonum) is pulled automatically via `go.mod`
> (used for the regularised incomplete beta function).

---

## 🚀 Quick start

### `shrink-nd` (main CLI)

#### 1. Chi-model (no CSV)

```bash
go run ./cmd/shrink-nd -mode=chi -n=3 -sigma=0.30 -rmax=0.95 -nGrid=2000
```

#### 2. Empirical (CSV of \$n\$-D points)

```bash
# Per-point shrinkage (each row gets its own d(r))
go run ./cmd/shrink-nd -mode=empirical -csv=examples/nd/data_3d.csv \
  -csv-out=examples/nd/scaled_3d_perpoint.csv -autoscale -apply=perpoint

# Global geometric factor (empirical mean of d(r)), applied to all rows
go run ./cmd/shrink-nd -mode=empirical -csv=examples/nd/data_3d.csv \
  -csv-out=examples/nd/scaled_3d_global_chi.csv -autoscale -apply=global-chi

# Global James–Stein-style factor, applied to all rows
go run ./cmd/shrink-nd -mode=empirical -csv=examples/nd/data_3d.csv \
  -csv-out=examples/nd/scaled_3d_global_js.csv -autoscale -apply=global-js
```

#### 3. Generate synthetic \$n\$-D data

```bash
# Generate 3D Gaussian around mean (0.3, 0.1, -0.2) with σ=0.25
go run ./cmd/shrink-nd \
  -gen=examples/nd/data_3d.csv -gen-n=3 -gen-rows=200 \
  -gen-sigma=0.25 -gen-seed=42 -gen-mean=0.3,0.1,-0.2
```

---

### `lwsf` — Prototype (2-D)

The original 2-D shrinkage rule:

* Computes \$d(r\_m)\$ so that a **unit circle** \$C\_o\$ divides the measurement circle \$C\_m\$ into **two equal areas**.
* Modes:

  * **Empirical** — shrink from measured radii in CSV.
  * **Rayleigh** — integrate over a Rayleigh PDF with given \$\sigma\$.

Example Rayleigh mode:

```bash
go run ./cmd/lwsf -mode=rayleigh -sigma=0.25 -rmax=0.8 -n=1000 -theta=45 -printrows
```

Example empirical mode:

```bash
go run ./cmd/lwsf -mode=empirical -csv=examples/2d/data_r_theta.csv \
  -csv-out=examples/2d/out_xy.csv -theta=45 -printrows
```

[Example CSVs here](./examples/2d)

---

## 📐 Theory

### Geometric shrinkage rule

Work in **scaled coordinates** with a unit reference \$n\$-ball \$B\_1\$ at the origin.

For a point with radius \$r\in(0,1)\$, define \$d(r)\in(|1-r|,1+r)\$ by the **half-volume condition**:

$$
\mathrm{Vol}\big(B_1 \cap B_r(d)\big) \;=\; \frac12\,\mathrm{Vol}(B_r).
$$

The intersection is the sum of **two spherical caps**. Each cap volume is

$$
V_{\text{cap}}^{(n)}(R,h) =
\frac12\, V_n(R)\; I_z\!\Big(\frac{n+1}{2},\frac12\Big),
\quad
z=\frac{2Rh-h^2}{R^2},
$$

where
\$V\_n(R)=\dfrac{\pi^{n/2}}{\Gamma(\frac{n}{2}+1)}R^n\$
and \$I\_z\$ is the **regularised incomplete beta**.

The half-volume equation is solved for \$d\$ via **Brent–Dekker** root finding.

---

## ⚙️ `shrink-nd` flags

<details>
<summary><strong>Mode & I/O</strong></summary>

* `-mode` — `"empirical"` or `"chi"` (default `"chi"`).
* `-csv` — input CSV of rows `x1,...,xn` (empirical mode).
* `-csv-out` — output CSV path (empirical mode).

</details>

<details>
<summary><strong>Empirical options</strong></summary>

* `-autoscale` — rescale dataset so `max radius < 1`.
* `-apply` — `perpoint` | `global-chi` | `global-js`.
* `-printrows` — print per-row diagnostics.

</details>

<details>
<summary><strong>Chi-model options</strong></summary>

* `-n` — dimension (required).
* `-sigma` — per-axis std for isotropic Gaussian.
* `-rmax` — integrate over \$r\in(0,r\_{\max}]\subset(0,1)\$.
* `-nGrid` — number of grid points.

</details>

<details>
<summary><strong>Generator options</strong></summary>

* `-gen` — write generated points and exit.
* `-gen-n` — dimension.
* `-gen-rows` — number of rows.
* `-gen-sigma` — per-axis std.
* `-gen-seed` — RNG seed.
* `-gen-mean` — comma-separated mean vector of length `gen-n`.

</details>

---

## 📝 Output format (empirical mode)

```
x_1, ..., x_n, r, d, x_1_scaled, ..., x_n_scaled
```

* `r` — Euclidean norm after optional autoscaling.
* `d` — shrinkage factor (per-row or global, depending on `-apply`).
* `x_j_scaled = d * x_j`.

---

## 📂 Examples

* [examples/nd/](./examples/nd) — generated 3D data and scaled outputs for each apply mode.
* [examples/2d/](./examples/2d) — earlier 2D inputs/outputs for `lwsf`.

---

## ❗ Troubleshooting

* **“failed to bracket root”** — ensure radii are in \$(0,1)\$. Use `-autoscale` if needed.
* **`global-js` very small or zero** — if \$|\bar{x}|\$ is close to the origin relative to noise, the JS factor may clip to zero.

---

## 📄 License

See [LICENSE](./LICENSE).

---

## 📚 Notes

* Long-form background and derivations: [docs/](./docs)
* Older experiments: [archive/](./archive)

```
