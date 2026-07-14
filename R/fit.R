#' @keywords internal
#' @noRd
tukey_adenet_fit_core <- function(x,
                                   y,
                                   lambda1,
                                   lambda2,
                                   weights,
                                   beta_init,
                                   sigma,
                                   d,
                                   eta,
                                   adagrad_eps,
                                   tol,
                                   max_iter,
                                   zero_tol,
                                   verbose) {
  n <- nrow(x)
  p <- ncol(x)

  beta <- as.numeric(beta_init)
  weights <- as.numeric(weights)

  g_accum <- numeric(p)
  loss_prev <- tukey_loss(x, y, beta, sigma = sigma, d = d, average = TRUE)
  converged <- FALSE
  iter <- 0L

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
  )
}

#' Fit a Tukey adaptive elastic net regression
#'
#' Fits a robust, sparse linear model by minimising Tukey's redescending
#' biweight loss subject to an adaptive elastic net penalty
#' (\eqn{\lambda_1 \sum_j \hat{w}_j |\beta_j| + \tfrac{\lambda_2}{2} \sum_j \beta_j^2}),
#' using a coordinate-wise proximal AdaGrad algorithm. Unlike the classical
#' squared-error loss, the biweight loss gives observations with
#' standardised residuals beyond `d` exactly zero influence on the fit,
#' so the estimator has a positive breakdown point under contamination in
#' the response and/or the design.
#'
#' `x` should not include an intercept column; the model is fitted on
#' (implicitly or explicitly) centred data, consistent with
#' \code{y = x \%*\% beta + error}. Standardising the columns of `x` prior
#' to fitting is recommended, since the elastic net penalty is not
#' scale-invariant.
#'
#' @param x Numeric predictor matrix (or an object coercible to one) with
#'   `n` rows and `p` columns. Must not contain an intercept column.
#' @param y Numeric response vector of length `n`.
#' @param lambda1 Non-negative adaptive-lasso tuning parameter (L1 penalty
#'   scale). Use [tukeyAdEnetRBIC()] to select this automatically.
#' @param lambda2 Non-negative ridge tuning parameter (L2 penalty scale).
#'   Defaults to `0`, i.e. an adaptive-lasso fit.
#' @param weights Optional numeric vector of length `p` giving the adaptive
#'   penalty weights \eqn{\hat{w}_j}. If `NULL`, weights are computed as
#'   \eqn{1 / |\tilde{\beta}_j|} from a robust (or ridge, if `p >= n`)
#'   initial fit.
#' @param beta_init Optional numeric vector of length `p` giving starting
#'   coefficient values for the AdaGrad iterations. If `NULL`, the same
#'   initial fit used for `weights` is reused.
#' @param sigma Optional positive robust scale estimate used to standardise
#'   residuals inside the biweight loss. If `NULL`, it is estimated by the
#'   (normalised) median absolute deviation of the residuals from the
#'   initial fit.
#' @param d Positive tuning constant for Tukey's biweight (see
#'   [tukeyRho()]). Defaults to `4.685` (approx. 95\% Gaussian efficiency).
#' @param eta Base AdaGrad step size.
#' @param adagrad_eps Small constant added to the AdaGrad denominator for
#'   numerical stability.
#' @param tol Convergence tolerance on the relative change in the
#'   coefficient vector and on the absolute change in the (averaged)
#'   objective value.
#' @param max_iter Maximum number of AdaGrad iterations.
#' @param zero_tol Coefficients with absolute value below this threshold
#'   are set to exact zero in the returned fit.
#' @param verbose If `TRUE`, print iteration progress.
#'
#' @return An object of class `"tukeyAdEnet"`, a list with components:
#'   \describe{
#'     \item{`beta`}{Named numeric vector of fitted (elastic-net-rescaled)
#'       coefficients, length `p`.}
#'     \item{`beta_raw`}{Coefficients before the `(1 + lambda2 / n)`
#'       rescaling.}
#'     \item{`sigma`}{The robust scale used in the biweight loss.}
#'     \item{`weights`}{The adaptive penalty weights used.}
#'     \item{`lambda1`, `lambda2`}{The tuning parameters used.}
#'     \item{`loss`}{Total (unaveraged) biweight loss at `beta`.}
#'     \item{`df`}{Number of non-zero coefficients.}
#'     \item{`iterations`}{Number of AdaGrad iterations performed.}
#'     \item{`converged`}{Logical; whether the convergence tolerance was
#'       reached before `max_iter`.}
#'     \item{`d`}{The Tukey tuning constant used.}
#'     \item{`call`}{The matched call.}
#'   }
#'
#' @seealso [tukeyAdEnetRBIC()] for automatic selection of `lambda1` and
#'   `lambda2` via a robust BIC grid search; [predict.tukeyAdEnet()],
#'   [coef.tukeyAdEnet()].
#'
#' @references
#' Zou, H. and Zhang, H. H. (2009). On the adaptive elastic-net with a
#' diverging number of parameters. *Annals of Statistics*, 37(4),
#' 1733-1751. \doi{10.1214/09-AOS699}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, -2, 1.5, rep(0, p - 3))
#' y <- as.numeric(x %*% beta + rnorm(n))
#' fit <- tukeyAdEnet(x, y, lambda1 = 0.2, lambda2 = 0.05)
#' coef(fit)
#'
#' @export
tukeyAdEnet <- function(x,
                         y,
                         lambda1,
                         lambda2 = 0,
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
  cl <- match.call()
  x <- validate_x(x)
  y <- validate_y(y, nrow(x))
  validate_scalar_nonneg(lambda1, "lambda1")
  validate_scalar_nonneg(lambda2, "lambda2")

  init <- resolve_init(x, y, weights = weights, beta_init = beta_init)
  weights <- init$weights
  beta_init <- init$beta

  if (is.null(sigma)) {
    sigma <- mad_sigma(y - as.numeric(x %*% beta_init))
  }
  validate_scalar_pos(sigma, "sigma")

  fit <- tukey_adenet_fit_core(
    x = x, y = y,
    lambda1 = lambda1, lambda2 = lambda2,
    weights = weights, beta_init = beta_init, sigma = sigma,
    d = d, eta = eta, adagrad_eps = adagrad_eps,
    tol = tol, max_iter = max_iter, zero_tol = zero_tol,
    verbose = verbose
  )

  names(fit$beta) <- colnames(x)
  names(fit$beta_raw) <- colnames(x)
  names(fit$weights) <- colnames(x)
  fit$d <- d
  fit$call <- cl

  structure(fit, class = "tukeyAdEnet")
}

#' @keywords internal
#' @noRd
resolve_init <- function(x, y, weights = NULL, beta_init = NULL) {
  if (is.null(weights) || is.null(beta_init)) {
    init <- initial_beta(x, y)
    if (is.null(beta_init)) beta_init <- init$beta
    if (is.null(weights)) weights <- init$weights
  }
  list(beta = as.numeric(beta_init), weights = as.numeric(weights))
}

#' @keywords internal
#' @noRd
validate_x <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (anyNA(x)) stop("`x` must not contain missing values.", call. = FALSE)
  if (nrow(x) < 2 || ncol(x) < 1) {
    stop("`x` must have at least 2 rows and 1 column.", call. = FALSE)
  }
  x
}

#' @keywords internal
#' @noRd
validate_y <- function(y, n) {
  y <- as.numeric(y)
  if (length(y) != n) {
    stop("`y` must have length equal to `nrow(x)`.", call. = FALSE)
  }
  if (anyNA(y)) stop("`y` must not contain missing values.", call. = FALSE)
  y
}

#' @keywords internal
#' @noRd
validate_scalar_nonneg <- function(v, name) {
  if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 0) {
    stop(sprintf("`%s` must be a single non-negative number.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
validate_scalar_pos <- function(v, name) {
  if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v <= 0) {
    stop(sprintf("`%s` must be a single positive number.", name), call. = FALSE)
  }
  invisible(TRUE)
}
