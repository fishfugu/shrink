# Shrinkage Factor Calculator

This program calculates the shrinkage factor $d(r_m)$ for a set of scaled radii $r_m \in (0,1)$ such that a unit circle $C_o$ divides the measurement circle $C_m$ into two equal areas.

It supports both:

* **Empirical mode** — using radii from a CSV file.
* **Rayleigh mode** — integrating over a Rayleigh distribution with given $\sigma$.

The program can:

1. Print per-row results for each $r_m$, including:

   ```
   r_m = 0.781000 -> d (scaling factor) = 0.890083405769 -> M = (0.629384012041, 0.629384012041)
   ```
2. Output a calculated **single shrink factor** — the weighted or empirical mean of $d(r_m)$.
3. (New) Apply this shrink factor to a 2D dataset from a CSV file, producing scaled coordinates.

---

## Usage

### Rayleigh mode

```bash
go run main.go -mode=rayleigh -sigma=0.25 -rmax=0.8 -n=1000 -theta=45 -printrows
```

* `-sigma`: Rayleigh $\sigma$ (per-axis std).
* `-rmax`: Upper limit for $r_m$ integration.
* `-n`: Number of integration points.
* `-theta`: Angle in degrees for reporting $M = (d\cos\theta, d\sin\theta)$.
* `-printrows`: Whether to print per-row results.

---

### Empirical mode with CSV

```bash
go run main.go -mode=empirical -csv=input.csv -csv-out=output.csv -theta=45 -printrows
```

* `-csv`: Path to input CSV containing $r_m$ in (0,1).
* `-csv-out`: Path to output CSV with scaled coordinates.
* `-theta`: Angle for reporting $M$.
* `-printrows`: Whether to print per-row results.

**CSV format (empirical mode)**:

* Must contain a header row with a column named `r_m` or have $r_m$ as the first column.
* Optional columns `x_m` and `y_m` will be scaled using the computed shrink factor:

  * `x_scaled = x_m * shrink_factor`
  * `y_scaled = y_m * shrink_factor`

---

## Output

Example empirical mode output:

```
Loaded 5 r_m values from input.csv

r_m = 0.781000 -> d (scaling factor) = 0.890083405769 -> M = (0.629384012041, 0.629384012041)
...
Single shrink factor (empirical mean of d(r_m)): 0.959712041447
Scaled data written to output.csv
```

Example Rayleigh mode output:

```
r_m = 0.100000 -> d (scaling factor) = 0.995004165278 -> M = (0.703562363973, 0.703562363973)
...
Single shrink factor (Rayleigh PDF-weighted, σ=0.250, rmax=0.800): 0.962385219704
```

---

## Theory

The calculation of $d(r_m)$ solves:

$$
A(d) = \frac{1}{2} \pi r_m^2
$$

where $A(d)$ is the area of intersection of two circles:

* $C_o$: radius 1, centred at origin.
* $C_m$: radius $r_m$, centre distance $d$ from origin.

The solution is found via Brent–Dekker root-finding on:

$$
F(d) = A(d) - \frac{1}{2} \pi r_m^2
$$

---

## Example datasets

You can generate artificial datasets:

```bash
# Rayleigh distributed radii
go run main.go -mode=rayleigh -sigma=0.3 -rmax=0.8 -n=50 -printrows

# Using real measured radii
go run main.go -mode=empirical -csv=my_data.csv -csv-out=my_scaled_data.csv
```
