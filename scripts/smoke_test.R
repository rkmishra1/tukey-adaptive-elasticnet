#!/usr/bin/env Rscript

source("R/tukey_adenet.R")
source("R/simulate_data.R")
source("R/metrics.R")
source("R/competitors.R")

set.seed(1)

dat <- generate_simulation_data(
  n = 80,
  p = 20,
  rho = 0.30,
  scenario = "response"
)

init <- initial_beta(dat$x, dat$y)
sigma <- mad_sigma(dat$y - as.numeric(dat$x %*% init$beta))
grid <- make_tuning_grid(
  x = dat$x,
  y = dat$y,
  beta_init = init$beta,
  sigma = sigma,
  weights = init$weights,
  n_lambda1 = 4,
  lambda2_factors = c(0, 0.1)
)

tuned <- tune_tukey_adenet_rbic(
  x = dat$x,
  y = dat$y,
  grid = grid,
  beta_init = init$beta,
  weights = init$weights,
  sigma = sigma,
  max_iter = 100,
  tol = 1e-5
)

sel <- selection_metrics(tuned$fit$beta, dat$beta)
mspe <- mspe_beta(tuned$fit$beta, dat$beta, dat$rho)

print(data.frame(
  method = "Tukey-AdEnet",
  C = unname(sel["C"]),
  IC = unname(sel["IC"]),
  MSPE = mspe,
  lambda1 = tuned$fit$lambda1,
  lambda2 = tuned$fit$lambda2,
  rbic = tuned$rbic,
  converged = tuned$fit$converged
))

stopifnot(length(tuned$fit$beta) == dat$p)
stopifnot(is.finite(mspe))
stopifnot(is.finite(tuned$rbic))

opts <- list(
  glmnet_nlambda = 8,
  adenet_alphas = c(0.3, 0.7),
  n_lambda1 = 4,
  lambda1_min_ratio = 0.02,
  lambda2_factors = c(0, 0.1),
  eta = 0.5,
  tol = 1e-5,
  max_iter = 100,
  zero_tol = 1e-8,
  lad_nlambda = 8,
  lad_max_iter = 500,
  sparse_lts_lambda = c(0.2, 0.1),
  rlars_s_max = NULL
)

for (method in c("AdL", "AdEnet", "Tukey-AdL", "Tukey-AdEnet")) {
  fit <- fit_competitor(method, dat$x, dat$y, opts)
  stopifnot(length(fit$beta) == dat$p)
  stopifnot(all(is.finite(fit$beta)))
}

cat("Competitor wrapper smoke checks passed for installed-package methods.\n")
