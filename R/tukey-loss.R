#' Tukey's biweight loss and influence functions
#'
#' `tukeyRho()` evaluates Tukey's (bisquare) biweight rho function, and
#' `tukeyPsi()` evaluates its derivative (the influence function). Both are
#' redescending: `tukeyPsi(u, d)` returns exactly zero once `abs(u) > d`, so
#' observations with standardised residuals beyond the tuning constant `d`
#' have no influence on the fit.
#'
#' @param u Numeric vector of (typically standardised) residuals.
#' @param d Positive tuning constant controlling the redescending cutoff.
#'   The default `4.685` gives approximately 95\% efficiency at the Gaussian
#'   model for location/regression M-estimation.
#'
#' @return A numeric vector of the same length as `u`.
#'
#' @references
#' Maronna, R. A., Martin, R. D., Yohai, V. J., & Salibian-Barrera, M. (2019).
#' *Robust Statistics: Theory and Methods (with R)*. Wiley.
#'
#' @examples
#' u <- seq(-8, 8, length.out = 401)
#' plot(u, tukeyRho(u), type = "l", main = "Tukey biweight loss")
#' plot(u, tukeyPsi(u), type = "l", main = "Tukey biweight influence")
#'
#' @export
tukeyRho <- function(u, d = 4.685) {
  inside <- abs(u) <= d
  out <- rep(d^2 / 6, length(u))
  v <- u[inside] / d
  out[inside] <- (d^2 / 6) * (1 - (1 - v^2)^3)
  out
}

#' @rdname tukeyRho
#' @export
tukeyPsi <- function(u, d = 4.685) {
  inside <- abs(u) <= d
  out <- numeric(length(u))
  v <- u[inside] / d
  out[inside] <- u[inside] * (1 - v^2)^2
  out
}

#' @keywords internal
#' @noRd
tukey_loss <- function(x, y, beta, sigma, d = 4.685, average = FALSE) {
  r <- as.numeric(y - x %*% beta)
  value <- sum(tukeyRho(r / sigma, d = d))
  if (average) value / length(y) else value
}

#' @keywords internal
#' @noRd
tukey_gradient <- function(x, y, beta, sigma, d = 4.685, average = TRUE) {
  r <- as.numeric(y - x %*% beta)
  psi <- tukeyPsi(r / sigma, d = d)
  grad <- -as.numeric(crossprod(x, psi)) / sigma
  if (average) grad / length(y) else grad
}
