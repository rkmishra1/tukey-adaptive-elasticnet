#' @keywords internal
#' @noRd
rbic_tukey_adenet <- function(x, y, fit, d = 4.685, zero_tol = 1e-8) {
  n <- nrow(x)
  loss <- tukey_loss(x, y, fit$beta, sigma = fit$sigma, d = d, average = FALSE)
  df <- sum(abs(fit$beta) > zero_tol)
  loss + log(n) * df
}

#' @keywords internal
#' @noRd
lambda1_max <- function(x, y, beta_init, sigma, weights, d = 4.685) {
  grad0 <- tukey_gradient(x, y, beta_init, sigma = sigma, d = d, average = TRUE)
  max(abs(grad0) / pmax(weights, 1e-12))
}

#' @keywords internal
#' @noRd
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

#' Tune a Tukey adaptive elastic net by robust BIC
#'
#' Fits [tukeyAdEnet()] over a two-dimensional `(lambda1, lambda2)` grid and
#' selects the pair that minimises a robust BIC (biweight loss plus a
#' `log(n)` complexity penalty on the number of non-zero coefficients).
#' This replaces cross-validation, which is itself sensitive to outliers in
#' the held-out folds.
#'
#' By default the grid spans `n_lambda1` values of `lambda1` on a log scale
#' from the smallest value that zeroes out every coefficient
#' (`lambda1_max`) down to `lambda1_min_ratio * lambda1_max`, crossed with
#' `lambda2 = lambda1_max * lambda2_factors`. Supply `grid` directly to use
#' a custom set of `(lambda1, lambda2)` pairs instead.
#'
#' @inheritParams tukeyAdEnet
#' @param grid Optional data frame with numeric columns `lambda1` and
#'   `lambda2` giving the tuning grid to search. If `NULL` (the default), a
#'   grid is constructed automatically (see Details).
#' @param n_lambda1 Number of `lambda1` values in the automatically
#'   constructed grid. Ignored if `grid` is supplied.
#' @param lambda1_min_ratio Ratio of the smallest to the largest `lambda1`
#'   in the automatically constructed grid. Ignored if `grid` is supplied.
#' @param lambda2_factors Numeric vector of multipliers applied to
#'   `lambda1_max` to obtain the `lambda2` grid values. Ignored if `grid`
#'   is supplied.
#'
#' @return An object of class `"tukeyAdEnetRBIC"`, a list with components:
#'   \describe{
#'     \item{`fit`}{The selected [tukeyAdEnet()] fit (class `"tukeyAdEnet"`).}
#'     \item{`rbic`}{The robust BIC value of the selected fit.}
#'     \item{`best_index`}{Row index of the selected pair in `grid`.}
#'     \item{`grid`}{The searched grid, with an appended `rbic` column.}
#'     \item{`sigma`}{The robust scale shared across all grid fits.}
#'     \item{`beta_init`, `weights`}{The initial coefficients and adaptive
#'       penalty weights shared across all grid fits.}
#'     \item{`call`}{The matched call.}
#'   }
#'
#' @seealso [tukeyAdEnet()] for fitting at a single `(lambda1, lambda2)`
#'   pair.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, -2, 1.5, rep(0, p - 3))
#' y <- as.numeric(x %*% beta + rnorm(n))
#' \donttest{
#' tuned <- tukeyAdEnetRBIC(x, y, n_lambda1 = 8, lambda2_factors = c(0, 0.1))
#' coef(tuned)
#' }
#'
#' @export
tukeyAdEnetRBIC <- function(x,
                             y,
                             grid = NULL,
                             weights = NULL,
                             beta_init = NULL,
                             sigma = NULL,
                             d = 4.685,
                             n_lambda1 = 20,
                             lambda1_min_ratio = 0.02,
                             lambda2_factors = c(0, 0.01, 0.05, 0.1, 0.5, 1),
                             eta = 0.5,
                             adagrad_eps = 1e-8,
                             tol = 1e-6,
                             max_iter = 2000,
                             zero_tol = 1e-8,
                             verbose = FALSE) {
  cl <- match.call()
  x <- validate_x(x)
  y <- validate_y(y, nrow(x))

  init <- resolve_init(x, y, weights = weights, beta_init = beta_init)
  weights <- init$weights
  beta_init <- init$beta

  if (is.null(sigma)) {
    sigma <- mad_sigma(y - as.numeric(x %*% beta_init))
  }
  validate_scalar_pos(sigma, "sigma")

  if (is.null(grid)) {
    grid <- make_tuning_grid(
      x = x, y = y, beta_init = beta_init, sigma = sigma, weights = weights,
      d = d, n_lambda1 = n_lambda1, lambda1_min_ratio = lambda1_min_ratio,
      lambda2_factors = lambda2_factors
    )
  } else {
    if (!all(c("lambda1", "lambda2") %in% names(grid))) {
      stop("`grid` must have columns `lambda1` and `lambda2`.", call. = FALSE)
    }
  }

  n_grid <- nrow(grid)
  fits <- vector("list", n_grid)
  rbic <- rep(Inf, n_grid)

  for (k in seq_len(n_grid)) {
    fits[[k]] <- tukey_adenet_fit_core(
      x = x, y = y,
      lambda1 = grid$lambda1[k], lambda2 = grid$lambda2[k],
      weights = weights, beta_init = beta_init, sigma = sigma,
      d = d, eta = eta, adagrad_eps = adagrad_eps,
      tol = tol, max_iter = max_iter, zero_tol = zero_tol,
      verbose = FALSE
    )
    rbic[k] <- rbic_tukey_adenet(x, y, fits[[k]], d = d, zero_tol = zero_tol)

    if (verbose) {
      message(
        "grid ", k, "/", n_grid,
        " lambda1=", signif(grid$lambda1[k], 4),
        " lambda2=", signif(grid$lambda2[k], 4),
        " rbic=", signif(rbic[k], 6)
      )
    }
  }

  best <- which.min(rbic)
  best_fit <- fits[[best]]
  names(best_fit$beta) <- colnames(x)
  names(best_fit$beta_raw) <- colnames(x)
  names(best_fit$weights) <- colnames(x)
  best_fit$d <- d
  best_fit$call <- cl
  best_fit <- structure(best_fit, class = "tukeyAdEnet")

  structure(
    list(
      fit = best_fit,
      rbic = rbic[best],
      best_index = best,
      grid = cbind(grid, rbic = rbic),
      sigma = sigma,
      beta_init = beta_init,
      weights = weights,
      call = cl
    ),
    class = "tukeyAdEnetRBIC"
  )
}
