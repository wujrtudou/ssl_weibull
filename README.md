# SSL mixed-missingness simulation code

This directory implements the four simulation studies in Section 5 of the manuscript.
Both the recorded-mechanism and latent-mechanism estimators are fitted by **direct numerical maximisation** of the appropriate observed-data likelihood using `nlminb`; no ECM iterations are used in the Monte Carlo runs.

## Final simulation design

1. **Parameter recovery** (`scripts/01_parameter_recovery.R`)
   - `p = 3`
   - baseline: `n=1000, Delta=2, gamma=0.30, xi1=3, alpha=0.10`
   - one-factor-at-a-time variation:
     - `n = 500, 1000`
     - `Delta = 1, 2, 3`
     - `gamma = 0.15, 0.30, 0.45`
     - `xi1 = 1, 3, 5`
     - `alpha = 0, 0.10, 0.30`
   - 10 distinct configurations, 500 replications each.

2. **Efficiency along the theorem path** (`scripts/02_efficiency.R`)
   - `p = 1, 5`, `Delta = 1, 2, 3`
   - `gamma = 0.30`, `xi1 = 3`
   - `alpha = 0, 0.10, 0.20, 0.30, 0.40`
   - `n = 1000` for the full grid (30 configurations)
   - `n = 500` additionally for `Delta=2` (10 configurations)
   - 500 replications each.

3. **Fixed overall missingness** (`scripts/03_fixed_missingness.R`)
   - `Gamma = 0.30, 0.50`, `p = 1, 5`, `alpha = 0, 0.10, 0.20`
   - `n=1000, Delta=2, xi1=3`
   - `gamma(alpha) = (Gamma-alpha)/(1-alpha)` and `xi0` recalibrated for every alpha.
   - 12 configurations, 500 replications each.

4. **Latent-mechanism robustness** (`scripts/04_latent_robustness.R`)
   - `p=5`, `n=500,1000`, `Delta=1,2`
   - `alpha=0,0.10,0.20,0.30`, `gamma=0.30`, `xi1=3`
   - 16 configurations, 500 replications each.
   - compares CC, recorded mixed, latent mixed, oracle latent, MCAR-only and MAR-only fits.
   - latent alpha is box-constrained on `[0,1)` so the MLE can attain exactly zero.

Total: 78 configurations and 39,000 Monte Carlo replications.

## Requirements

The core simulation uses base R only (`stats`, `parallel`). No contributed package is required.

## Run

From this directory:

```r
source("run_all.R")
```

or run one study, e.g.

```r
source("scripts/04_latent_robustness.R")
run_latent_robustness(B = 500, n_cores = 8)
```

Each configuration is checkpointed to `results/raw/<study>/`. Re-running a script resumes completed work unless `overwrite=TRUE`.

For a quick code-path check before the full Monte Carlo run:

```r
source("scripts/99_smoke_test.R")
run_smoke_test()
```

## Numerical details

- Gaussian model: two classes, common covariance matrix.
- Mixture proportion: `pi1 = plogis(eta_pi)`.
- Positive entropy slope: `xi1 = exp(eta_xi)`.
- Positive-definite common covariance: lower-triangular Cholesky factor with log-diagonal.
- `xi0` is calibrated by one-dimensional numerical integration over the Gaussian discriminant score followed by `uniroot`.
- Entropy is recomputed from the current mixture parameter at every likelihood evaluation.
- Recorded mixed likelihood uses the observed mutually-exclusive `m1`, `m2` indicators and `alpha_hat = mean(m1)`.
- Latent mixed likelihood uses only `m = m1 + m2` and directly optimises alpha with a box constraint.
- Multiple starts are retained by maximum observed-data log likelihood.

## Output

Raw per-replication RDS files are written under `results/raw/`. Study-level CSV summaries are written under `results/summary/`.

The optional functions in `R/06_theory_are.R` numerically approximate the leading-order ARE from expected observed information. This is deliberately separated from the Monte Carlo runners because it is substantially more computationally expensive than the simulation summaries.
