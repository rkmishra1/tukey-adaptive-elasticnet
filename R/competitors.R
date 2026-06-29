# Competitor estimators used in the computational study.

available_methods <- function() {
  c("AdL", "AdEnet", "HAdL", "Tukey-AdL", "S-LTS", "R-LARS", "Tukey-AdEnet")
}

required_package <- function(method) {
  switch(
    method,
    "AdL" = "glmnet",
    "AdEnet" = "glmnet",
    "HAdL" = "hqreg",
    "S-LTS" = "robustHD",
    "R-LARS" = "robustHD",
    NULL
  )
}

method_is_available <- function(method) {
  pkg <- required_package(method)
  is.null(pkg) || requireNamespace(pkg, quietly = TRUE)
}

missing_package_message <- function(method) {
  pkg <- required_package(method)
  if (is.null(pkg)) return(NULL)
  paste0("Method ", method, " requires the R package '", pkg, "'.")
}

as_beta_vector <- function(coefs, p, zero_tol = 1e-8) {
  if (is.list(coefs) && !is.data.frame(coefs)) coefs <- coefs[[1]]
  coefs <- as.matrix(coefs)

  if (nrow(coefs) == p || nrow(coefs) == p + 1) {
    beta <- coefs[, ncol(coefs)]
  } else if (ncol(coefs) == p || ncol(coefs) == p + 1) {
    beta <- coefs[nrow(coefs), ]
  } else {
    beta <- as.numeric(coefs)
  }

  beta <- as.numeric(beta)
  if (length(beta) == p + 1) beta <- beta[-1]
  if (length(beta) != p) {
    stop("Could not extract a coefficient vector of length ", p, ".", call. = FALSE)
  }

  beta[abs(beta) < zero_tol] <- 0
  beta
}

bic_gaussian <- function(x, y, beta, zero_tol = 1e-8) {
  n <- nrow(x)
  rss <- sum((as.numeric(y) - as.numeric(x %*% beta))^2)
  rss <- max(rss, .Machine$double.eps)
  df <- sum(abs(beta) > zero_tol)
  n * log(rss / n) + log(n) * df
}

fit_glmnet_bic <- function(x,
                           y,
                           alpha,
                           penalty_factor = NULL,
                           nlambda = 100,
                           zero_tol = 1e-8) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required.", call. = FALSE)
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  p <- ncol(x)
  if (is.null(penalty_factor)) penalty_factor <- rep(1, p)

  fit <- glmnet::glmnet(
    x = x,
    y = y,
    family = "gaussian",
    alpha = alpha,
    nlambda = nlambda,
    penalty.factor = penalty_factor,
    intercept = FALSE,
    standardize = FALSE
  )

  betas <- as.matrix(stats::coef(fit, s = fit$lambda))[-1, , drop = FALSE]
  bics <- apply(betas, 2, function(beta) bic_gaussian(x, y, beta, zero_tol = zero_tol))
  best <- which.min(bics)
  beta <- as.numeric(betas[, best])
  beta[abs(beta) < zero_tol] <- 0

  list(
    beta = beta,
    criterion = bics[best],
    lambda = fit$lambda[best],
    lambda1 = NA_real_,
    lambda2 = NA_real_,
    df = sum(abs(beta) > zero_tol),
    iterations = NA_integer_,
    converged = TRUE
  )
}

fit_adaptive_lasso <- function(x,
                               y,
                               gamma = 1,
                               eps = 1e-4,
                               nlambda = 100,
                               zero_tol = 1e-8) {
  beta0 <- as.numeric(ridge_initial_beta(as.matrix(x), as.numeric(y)))
  weights <- pmin(1 / (abs(beta0) + eps)^gamma, 1e6)
  out <- fit_glmnet_bic(
    x = x,
    y = y,
    alpha = 1,
    penalty_factor = weights,
    nlambda = nlambda,
    zero_tol = zero_tol
  )
  out$method <- "AdL"
  out
}

fit_adaptive_enet <- function(x,
                              y,
                              alphas = c(0.1, 0.3, 0.5, 0.7, 0.9),
                              gamma = 1,
                              eps = 1e-4,
                              nlambda = 100,
                              zero_tol = 1e-8) {
  beta0 <- as.numeric(ridge_initial_beta(as.matrix(x), as.numeric(y)))
  weights <- pmin(1 / (abs(beta0) + eps)^gamma, 1e6)

  fits <- lapply(alphas, function(alpha) {
    out <- fit_glmnet_bic(
      x = x,
      y = y,
      alpha = alpha,
      penalty_factor = weights,
      nlambda = nlambda,
      zero_tol = zero_tol
    )
    out$alpha <- alpha
    out
  })

  best <- which.min(vapply(fits, `[[`, numeric(1), "criterion"))
  out <- fits[[best]]
  out$method <- "AdEnet"
  out
}

fit_huber_adl <- function(x,
                          y,
                          nlambda = 100,
                          zero_tol = 1e-8) {
  if (!requireNamespace("hqreg", quietly = TRUE)) {
    stop("Package 'hqreg' is required.", call. = FALSE)
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  p <- ncol(x)
  
  # Obtain robust adaptive weights
  init <- initial_beta(x, y)
  weights <- init$weights

  fit <- hqreg::hqreg(
    X = x,
    y = y,
    method = "huber",
    alpha = 1,
    nlambda = nlambda,
    penalty.factor = weights
  )

  betas <- fit$beta[-1, , drop = FALSE]
  
  # Compute BIC
  n <- nrow(x)
  bics <- apply(betas, 2, function(b) {
    r <- y - as.numeric(x %*% b)
    g <- fit$gamma
    h_loss <- ifelse(abs(r) <= g, 0.5 * r^2, g * (abs(r) - 0.5 * g))
    df <- sum(abs(b) > zero_tol)
    2 * sum(h_loss) + log(n) * df
  })

  best <- which.min(bics)
  beta <- as.numeric(betas[, best])
  beta[abs(beta) < zero_tol] <- 0

  list(
    method = "HAdL",
    beta = beta,
    criterion = bics[best],
    lambda = fit$lambda[best],
    lambda1 = fit$lambda[best],
    lambda2 = 0,
    df = sum(abs(beta) > zero_tol),
    iterations = NA_integer_,
    converged = TRUE
  )
}

fit_tukey_adl <- function(x,
                          y,
                          n_lambda1 = 20,
                          lambda1_min_ratio = 0.02,
                          eta = 0.5,
                          tol = 1e-6,
                          max_iter = 2000,
                          zero_tol = 1e-8) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  init <- initial_beta(x, y)
  sigma <- mad_sigma(y - as.numeric(x %*% init$beta))
  grid <- make_tuning_grid(
    x = x,
    y = y,
    beta_init = init$beta,
    sigma = sigma,
    weights = init$weights,
    n_lambda1 = n_lambda1,
    lambda1_min_ratio = lambda1_min_ratio,
    lambda2_factors = 0
  )

  tuned <- tune_tukey_adenet_rbic(
    x = x,
    y = y,
    grid = grid,
    beta_init = init$beta,
    weights = init$weights,
    sigma = sigma,
    eta = eta,
    tol = tol,
    max_iter = max_iter,
    zero_tol = zero_tol
  )

  list(
    method = "Tukey-AdL",
    beta = tuned$fit$beta,
    criterion = tuned$rbic,
    lambda = tuned$fit$lambda1,
    lambda1 = tuned$fit$lambda1,
    lambda2 = 0,
    df = tuned$fit$df,
    iterations = tuned$fit$iterations,
    converged = tuned$fit$converged
  )
}

fit_sparse_lts <- function(x,
                           y,
                           lambda = seq(0.20, 0.05, by = -0.05),
                           mode = "fraction",
                           zero_tol = 1e-8) {
  if (!requireNamespace("robustHD", quietly = TRUE)) {
    stop("Package 'robustHD' is required.", call. = FALSE)
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  fit <- robustHD::sparseLTS(
    x = x,
    y = y,
    lambda = lambda,
    mode = mode,
    intercept = FALSE,
    crit = "BIC"
  )
  beta <- as_beta_vector(stats::coef(fit, fit = "reweighted", zeros = TRUE), ncol(x), zero_tol)
  bic <- suppressWarnings(stats::BIC(fit))

  list(
    method = "S-LTS",
    beta = beta,
    criterion = suppressWarnings(min(bic, na.rm = TRUE)),
    lambda = if (!is.null(fit$lambda) && length(bic)) fit$lambda[which.min(bic)] else NA_real_,
    lambda1 = NA_real_,
    lambda2 = NA_real_,
    df = sum(abs(beta) > zero_tol),
    iterations = NA_integer_,
    converged = TRUE
  )
}

fit_robust_lars <- function(x,
                            y,
                            s_max = NULL,
                            zero_tol = 1e-8) {
  if (!requireNamespace("robustHD", quietly = TRUE)) {
    stop("Package 'robustHD' is required.", call. = FALSE)
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  if (is.null(s_max)) s_max <- min(ncol(x), max(1, floor(nrow(x) / 2) - 1))

  fit <- robustHD::rlars(
    x = x,
    y = y,
    sMax = s_max,
    crit = "BIC"
  )
  beta <- as_beta_vector(stats::coef(fit, zeros = TRUE), ncol(x), zero_tol)
  bic <- suppressWarnings(stats::BIC(fit))

  list(
    method = "R-LARS",
    beta = beta,
    criterion = suppressWarnings(min(bic, na.rm = TRUE)),
    lambda = NA_real_,
    lambda1 = NA_real_,
    lambda2 = NA_real_,
    df = sum(abs(beta) > zero_tol),
    iterations = NA_integer_,
    converged = TRUE
  )
}

fit_tukey_adenet_method <- function(x,
                                    y,
                                    n_lambda1 = 20,
                                    lambda1_min_ratio = 0.02,
                                    lambda2_factors = c(0, 0.01, 0.05, 0.1, 0.5, 1),
                                    eta = 0.5,
                                    tol = 1e-6,
                                    max_iter = 2000,
                                    zero_tol = 1e-8) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  init <- initial_beta(x, y)
  sigma <- mad_sigma(y - as.numeric(x %*% init$beta))
  grid <- make_tuning_grid(
    x = x,
    y = y,
    beta_init = init$beta,
    sigma = sigma,
    weights = init$weights,
    n_lambda1 = n_lambda1,
    lambda1_min_ratio = lambda1_min_ratio,
    lambda2_factors = lambda2_factors
  )

  tuned <- tune_tukey_adenet_rbic(
    x = x,
    y = y,
    grid = grid,
    beta_init = init$beta,
    weights = init$weights,
    sigma = sigma,
    eta = eta,
    tol = tol,
    max_iter = max_iter,
    zero_tol = zero_tol
  )

  list(
    method = "Tukey-AdEnet",
    beta = tuned$fit$beta,
    criterion = tuned$rbic,
    lambda = NA_real_,
    lambda1 = tuned$fit$lambda1,
    lambda2 = tuned$fit$lambda2,
    df = tuned$fit$df,
    iterations = tuned$fit$iterations,
    converged = tuned$fit$converged
  )
}

fit_competitor <- function(method, x, y, opts) {
  switch(
    method,
    "AdL" = fit_adaptive_lasso(x, y, nlambda = opts$glmnet_nlambda, zero_tol = opts$zero_tol),
    "AdEnet" = fit_adaptive_enet(
      x, y,
      alphas = opts$adenet_alphas,
      nlambda = opts$glmnet_nlambda,
      zero_tol = opts$zero_tol
    ),
    "HAdL" = fit_huber_adl(
      x, y,
      nlambda = opts$glmnet_nlambda,
      zero_tol = opts$zero_tol
    ),
    "Tukey-AdL" = fit_tukey_adl(
      x, y,
      n_lambda1 = opts$n_lambda1,
      lambda1_min_ratio = opts$lambda1_min_ratio,
      eta = opts$eta,
      tol = opts$tol,
      max_iter = opts$max_iter,
      zero_tol = opts$zero_tol
    ),
    "S-LTS" = fit_sparse_lts(
      x, y,
      lambda = opts$sparse_lts_lambda,
      zero_tol = opts$zero_tol
    ),
    "R-LARS" = fit_robust_lars(
      x, y,
      s_max = opts$rlars_s_max,
      zero_tol = opts$zero_tol
    ),
    "Tukey-AdEnet" = fit_tukey_adenet_method(
      x, y,
      n_lambda1 = opts$n_lambda1,
      lambda1_min_ratio = opts$lambda1_min_ratio,
      lambda2_factors = opts$lambda2_factors,
      eta = opts$eta,
      tol = opts$tol,
      max_iter = opts$max_iter,
      zero_tol = opts$zero_tol
    ),
    stop("Unknown method: ", method, call. = FALSE)
  )
}
