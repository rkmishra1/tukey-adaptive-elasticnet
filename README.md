# Tukey Adaptive Elasticnet Simulations

This repository contains R code for the computational study of the Tukey adaptive Elasticnet estimator (`Tukey-AdEnet`) fitted by proximal AdaGrad and tuned by robust BIC.

The code follows the manuscript setup:

- Linear model `y = X beta + e`.
- AR(1) Gaussian design with correlations `rho = 0.30, 0.60, 0.80`.
- Active set size `s = 3 * floor(p / 9)`.
- Three scenarios: clean data, response contamination, and combined response/design contamination.
- Three dimensional regimes: `zeta_1_2`, `zeta_2_3`, and `zeta_5_6`.
- RBIC selection over a two-dimensional `(lambda1, lambda2)` grid for `Tukey-AdEnet`.

## Files

- `R/tukey_adenet.R`: Tukey loss, gradient, proximal AdaGrad estimator, adaptive weights, and RBIC tuning.
- `R/simulate_data.R`: data-generating mechanisms and manuscript simulation grid.
- `R/metrics.R`: MSPE and variable-selection metrics.
- `R/competitors.R`: wrappers for AdL, AdEnet, LAD-Lasso, Tukey-AdL, S-LTS, R-LARS, and Tukey-AdEnet.
- `scripts/run_simulation.R`: command-line simulation runner.
- `scripts/install_packages.R`: installs comparison-method packages from CRAN.
- `scripts/smoke_test.R`: small quick check that the estimator and RBIC tuning run.
- `scripts/make_boxplots.R`: creates MSPE boxplots from a raw results CSV.

## Dependencies

```sh
Rscript scripts/install_packages.R
```

The full seven-method comparison needs `glmnet`, `rqPen`, `robustHD`, and `robustbase`. If a package is missing, `scripts/run_simulation.R` records that method as skipped by default. For manuscript runs, use `--missing_action=stop` so missing packages stop the job immediately.

## Quick Check

```sh
Rscript scripts/smoke_test.R
```

## Full Manuscript Grid

The full study uses 200 replications across all 81 configurations. This is computationally expensive because each replication fits a two-dimensional RBIC grid.

```sh
Rscript scripts/run_simulation.R --reps=200 --missing_action=stop --output_dir=results
```

For a smaller pilot run:

```sh
Rscript scripts/run_simulation.R \
  --reps=5 \
  --scenarios=response_design \
  --regimes=zeta_2_3 \
  --n_values=500 \
  --rhos=0.30 \
  --n_lambda1=8 \
  --lambda2_factors=0,0.1,0.5 \
  --max_iter=500 \
  --output_dir=results
```

To run only a subset of methods:

```sh
Rscript scripts/run_simulation.R --methods=AdL,AdEnet,Tukey-AdL,Tukey-AdEnet --reps=10
```

The runner writes two CSV files:

- `comparison_raw_*.csv`: one row per method and replication.
- `comparison_summary_*.csv`: manuscript-style averages by method and configuration.

## Boxplots

```sh
Rscript scripts/make_boxplots.R --raw=results/comparison_raw_YYYYMMDD_HHMMSS.csv --output_dir=figures
```

## Main Output Columns

- `C`: number of true zero coefficients estimated as zero.
- `IC`: number of true nonzero coefficients incorrectly estimated as zero.
- `MSPE`: `(beta_hat - beta_true)' Sigma (beta_hat - beta_true)` with AR(1) `Sigma`.
- `lambda1`, `lambda2`: RBIC-selected tuning parameters.
- `criterion`: BIC/RBIC value for the selected model, depending on the method.
- `converged`: whether proximal AdaGrad met the stopping tolerance.

## Notes

The estimator uses the Tukey biweight score

```text
psi(u) = u * (1 - (u / d)^2)^2 * I(|u| <= d)
```

and the coordinate-wise proximal AdaGrad update

```text
u_j    = beta_j - eta_j * grad_j
beta_j = sign(u_j) max(|u_j| - eta_j lambda1 w_j, 0) / (1 + eta_j lambda2)
```

If the `robustbase` package is installed and `p < n`, the initialization uses `lmrob`. Otherwise, the code falls back to ridge initialization so the high-dimensional settings remain runnable without extra package requirements.
