# tukeyAdEnet

> Robust sparse regression via Tukey's biweight loss + adaptive elastic net penalty, fitted by proximal AdaGrad and tuned by robust BIC.

![R](https://img.shields.io/badge/R-%3E%3D3.6-276DC3?logo=r&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-blue)

**tukeyAdEnet** is an R package for robust, sparse linear regression. It
replaces the squared-error loss of classical penalized regression (lasso,
elastic net, adaptive elastic net) with **Tukey's redescending biweight
loss**, so that observations with standardised residuals beyond a tuning
constant `d` have *exactly zero* influence on the fit — rather than merely
bounded influence, as with Huber-type losses. The estimator retains the
variable-selection properties of the adaptive elastic net penalty and is
fitted with a coordinate-wise proximal AdaGrad algorithm. Regularisation
parameters are chosen automatically by minimising a **robust BIC** over a
two-dimensional `(lambda1, lambda2)` grid, avoiding the need for
cross-validation (which is itself sensitive to outliers landing in a
held-out fold).

## Installation

```r
# from GitHub, until the package is on CRAN
# install.packages("remotes")
remotes::install_github("rkmishra1/tukey-adaptive-elasticnet")
```

## Usage

```r
library(tukeyAdEnet)

set.seed(1)
n <- 150; p <- 12
x <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -3, 2, rep(0, p - 3))
y <- as.numeric(x %*% beta_true + rnorm(n))

# contaminate 10% of the responses with gross outliers
idx <- sample(seq_len(n), size = round(0.1 * n))
y[idx] <- y[idx] + rnorm(length(idx), mean = 30, sd = 5)

# fit at a single (lambda1, lambda2) pair
fit <- tukeyAdEnet(x, y, lambda1 = 0.3, lambda2 = 0.05)
coef(fit)

# or tune both parameters automatically via robust BIC
tuned <- tukeyAdEnetRBIC(x, y)
coef(tuned)
predict(tuned, newx = x[1:5, ])
plot(tuned)
```

See `vignette("tukeyAdEnet")` for a full walk-through, including a
comparison against ordinary least squares under contamination.

## Method

The estimator minimises

```
sum_i rho_d((y_i - x_i'beta) / sigma) + lambda1 * sum_j w_j |beta_j| + (lambda2 / 2) * sum_j beta_j^2
```

where `rho_d` is Tukey's biweight loss, `sigma` is a robust scale estimate
(normalised MAD of residuals from an initial fit), and `w_j = 1 / |beta_init_j|`
are adaptive weights from an initial robust (or ridge, when `p >= n`) fit.
The non-convex, non-smooth objective is minimised by a coordinate-wise
proximal AdaGrad update:

```
u_j  = beta_j - eta_j * grad_j
beta_j = sign(u_j) * max(|u_j| - eta_j * lambda1 * w_j, 0) / (1 + eta_j * lambda2)
```

`lambda1` and `lambda2` are selected by minimising a robust BIC (the
biweight loss plus a `log(n)` penalty on the number of non-zero
coefficients) over a two-dimensional grid.

## Simulation study & manuscript reproduction

The [`simulation-study/`](simulation-study/) directory contains the full
simulation study, real-data benchmarks (against `AdL`, `AdEnet`,
`LAD-Lasso`, `S-LTS`, `R-LARS`, `Tukey-AdL`), and figures accompanying the
manuscript that introduces this estimator. That material is not part of
the installable package (it depends on several additional packages used
only for comparison) — see [`simulation-study/README.md`](simulation-study/README.md)
for details and instructions to reproduce it.

## License

MIT © Ramakrushna Mishra
