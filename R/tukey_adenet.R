# Tukey adaptive Elasticnet fitted by proximal AdaGrad.

tukey_rho <- function(u, d = 4.685) {
  inside <- abs(u) <= d
  out <- rep(d^2 / 6, length(u))
  v <- u[inside] / d
  out[inside] <- (d^2 / 6) * (1 - (1 - v^2)^3)
  out
}

tukey_psi <- function(u, d = 4.685) {
  inside <- abs(u) <= d
  out <- numeric(length(u))
  v <- u[inside] / d
  out[inside] <- u[inside] * (1 - v^2)^2
  out
}

tukey_loss <- function(x, y, beta, sigma, d = 4.685, average = FALSE) {
  r <- as.numeric(y - x %*% beta)
  value <- sum(tukey_rho(r / sigma, d = d))
  if (average) value / length(y) else value
}

tukey_gradient <- function(x, y, beta, sigma, d = 4.685, average = TRUE) {
  r <- as.numeric(y - x %*% beta)
  psi <- tukey_psi(r / sigma, d = d)
  grad <- -as.numeric(crossprod(x, psi)) / sigma
  if (average) grad / length(y) else grad
}

soft_threshold <- function(z, threshold) {
  sign(z) * pmax(abs(z) - threshold, 0)
}

mad_sigma <- function(residuals) {
  sigma <- stats::mad(as.numeric(residuals), center = 0, constant = 1.4826)
  if (!is.finite(sigma) || sigma <= 1e-4) {
    sigma <- stats::sd(as.numeric(residuals))
  }
  if (!is.finite(sigma) || sigma <= 1e-4) sigma <- 1e-4
  sigma
}

ridge_initial_beta <- function(x, y, lambda = NULL) {
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(lambda)) lambda <- if (p >= n) 1 else 1e-6

  if (p <= n) {
    solve(crossprod(x) + diag(lambda, p), crossprod(x, y))
  } else {
    as.numeric(crossprod(x, solve(tcrossprod(x) + diag(lambda, n), y)))
  }
}

initial_beta <- function(x, y, gamma = 1, eps = 1e-4, weight_cap = 1e6) {
  p <- ncol(x)

  beta0 <- NULL
  if (p < nrow(x) && requireNamespace("robustbase", quietly = TRUE)) {
    beta0 <- tryCatch(
      {
        data <- data.frame(y = as.numeric(y), x)
        stats::coef(robustbase::lmrob(y ~ . - 1, data = data, setting = "KS2014"))
      },
      error = function(e) NULL
    )
  }

  if (is.null(beta0) || any(!is.finite(beta0))) {
    beta0 <- as.numeric(ridge_initial_beta(x, y))
  }

  beta0 <- as.numeric(beta0)
  weights <- 1 / (abs(beta0) + eps)^gamma
  weights <- pmin(weights, weight_cap)

  list(beta = beta0, weights = weights)
}

tukey_adenet_fit <- function(x,
                             y,
                             lambda1,
                             lambda2,
                             weights = NULL,
                             beta_init = NULL,
                             sigma = NULL,
                             d = 4.685,
                             eta = 0.5,
                             adagrad_eps = 1e-8,
                             tol = 1e-6,
                             max_iter = 2000,
                             zero_tol = 1e-8,
                             verbose = FALSE) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(beta_init) || is.null(weights)) {
    init <- initial_beta(x, y)
    if (is.null(beta_init)) beta_init <- init$beta
    if (is.null(weights)) weights <- init$weights
  }

  beta <- as.numeric(beta_init)
  weights <- as.numeric(weights)

  if (is.null(sigma)) {
    sigma <- mad_sigma(y - as.numeric(x %*% beta))
  }

  g_accum <- numeric(p)
  loss_prev <- tukey_loss(x, y, beta, sigma = sigma, d = d, average = TRUE)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    grad <- tukey_gradient(x, y, beta, sigma = sigma, d = d, average = TRUE)
    g_accum <- g_accum + grad^2
    rates <- eta / (sqrt(g_accum) + adagrad_eps)

    u <- beta - rates * grad
    beta_new <- soft_threshold(u, rates * lambda1 * weights) / (1 + rates * lambda2)

    loss_new <- tukey_loss(x, y, beta_new, sigma = sigma, d = d, average = TRUE)
    delta_loss <- abs(loss_new - loss_prev)
    delta_beta <- sqrt(sum((beta_new - beta)^2)) / (sqrt(sum(beta^2)) + 1e-12)

    beta <- beta_new
    loss_prev <- loss_new

    if (verbose && (iter == 1 || iter %% 50 == 0)) {
      message(
        "iter=", iter,
        " loss=", signif(loss_new, 6),
        " delta_loss=", signif(delta_loss, 4),
        " delta_beta=", signif(delta_beta, 4)
      )
    }

    if (delta_loss < tol || delta_beta < tol) {
      converged <- TRUE
      break
    }
  }

  beta_hat <- (1 + lambda2 / n) * beta
  beta_hat[abs(beta_hat) < zero_tol] <- 0

  structure(
    list(
      beta = beta_hat,
      beta_raw = beta,
      sigma = sigma,
      weights = weights,
      lambda1 = lambda1,
      lambda2 = lambda2,
      loss = tukey_loss(x, y, beta_hat, sigma = sigma, d = d, average = FALSE),
      df = sum(abs(beta_hat) > zero_tol),
      iterations = iter,
      converged = converged
    ),
    class = "tukey_adenet_fit"
  )
}

rbic_tukey_adenet <- function(x, y, fit, d = 4.685, zero_tol = 1e-8) {
  n <- nrow(x)
  loss <- tukey_loss(x, y, fit$beta, sigma = fit$sigma, d = d, average = FALSE)
  df <- sum(abs(fit$beta) > zero_tol)
  loss + log(n) * df
}

lambda1_max <- function(x, y, beta_init, sigma, weights, d = 4.685) {
  grad0 <- tukey_gradient(x, y, beta_init, sigma = sigma, d = d, average = TRUE)
  max(abs(grad0) / pmax(weights, 1e-12))
}

make_tuning_grid <- function(x,
                             y,
                             beta_init,
                             sigma,
                             weights,
                             d = 4.685,
                             n_lambda1 = 20,
                             lambda1_min_ratio = 0.02,
                             lambda2_factors = c(0, 0.01, 0.05, 0.1, 0.5, 1)) {
  lam_max <- lambda1_max(x, y, beta_init, sigma, weights, d = d)
  if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1

  lambda1 <- exp(seq(log(lam_max), log(lam_max * lambda1_min_ratio), length.out = n_lambda1))
  lambda2 <- lam_max * lambda2_factors

  expand.grid(lambda1 = lambda1, lambda2 = lambda2, KEEP.OUT.ATTRS = FALSE)
}

tune_tukey_adenet_rbic <- function(x,
                                   y,
                                   grid = NULL,
                                   beta_init = NULL,
                                   weights = NULL,
                                   sigma = NULL,
                                   d = 4.685,
                                   eta = 0.5,
                                   tol = 1e-6,
                                   max_iter = 2000,
                                   zero_tol = 1e-8,
                                   verbose = FALSE) {
  x <- as.matrix(x)
  y <- as.numeric(y)

  if (is.null(beta_init) || is.null(weights)) {
    init <- initial_beta(x, y)
    if (is.null(beta_init)) beta_init <- init$beta
    if (is.null(weights)) weights <- init$weights
  }

  if (is.null(sigma)) {
    sigma <- mad_sigma(y - as.numeric(x %*% beta_init))
  }

  if (is.null(grid)) {
    grid <- make_tuning_grid(
      x = x,
      y = y,
      beta_init = beta_init,
      sigma = sigma,
      weights = weights,
      d = d
    )
  }

  fits <- vector("list", nrow(grid))
  rbic <- rep(Inf, nrow(grid))

  for (k in seq_len(nrow(grid))) {
    fits[[k]] <- tukey_adenet_fit(
      x = x,
      y = y,
      lambda1 = grid$lambda1[k],
      lambda2 = grid$lambda2[k],
      beta_init = beta_init,
      weights = weights,
      sigma = sigma,
      d = d,
      eta = eta,
      tol = tol,
      max_iter = max_iter,
      zero_tol = zero_tol,
      verbose = FALSE
    )
    rbic[k] <- rbic_tukey_adenet(x, y, fits[[k]], d = d, zero_tol = zero_tol)

    if (verbose) {
      message(
        "grid ", k, "/", nrow(grid),
        " lambda1=", signif(grid$lambda1[k], 4),
        " lambda2=", signif(grid$lambda2[k], 4),
        " rbic=", signif(rbic[k], 6)
      )
    }
  }

  best <- which.min(rbic)
  list(
    fit = fits[[best]],
    rbic = rbic[best],
    best_index = best,
    grid = cbind(grid, rbic = rbic),
    sigma = sigma,
    beta_init = beta_init,
    weights = weights
  )
}
