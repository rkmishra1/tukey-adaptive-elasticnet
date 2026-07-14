#' @keywords internal
#' @noRd
soft_threshold <- function(z, threshold) {
  sign(z) * pmax(abs(z) - threshold, 0)
}

#' @keywords internal
#' @noRd
mad_sigma <- function(residuals) {
  sigma <- stats::mad(as.numeric(residuals), center = 0, constant = 1.4826)
  if (!is.finite(sigma) || sigma <= 1e-4) {
    sigma <- stats::sd(as.numeric(residuals))
  }
  if (!is.finite(sigma) || sigma <= 1e-4) sigma <- 1e-4
  sigma
}

#' @keywords internal
#' @noRd
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

#' @keywords internal
#' @noRd
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
