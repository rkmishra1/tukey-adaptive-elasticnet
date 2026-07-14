#' Inspect a fitted Tukey adaptive elastic net
#'
#' `coef()` extracts the fitted coefficient vector, and `print()` /
#' `summary()` display a short description of the fit.
#'
#' @param object,x A `"tukeyAdEnet"` object as returned by [tukeyAdEnet()]
#'   or `tukeyAdEnetRBIC(...)$fit`.
#' @param ... Currently ignored.
#'
#' @return `coef()` returns the named numeric coefficient vector. `print()`
#'   and `summary()` return `object`, invisibly.
#'
#' @name tukeyAdEnet-methods
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, -2, 1.5, rep(0, p - 3))
#' y <- as.numeric(x %*% beta + rnorm(n))
#' fit <- tukeyAdEnet(x, y, lambda1 = 0.2, lambda2 = 0.05)
#' fit
#' summary(fit)
#' coef(fit)
NULL

#' @rdname tukeyAdEnet-methods
#' @export
coef.tukeyAdEnet <- function(object, ...) {
  object$beta
}

#' @rdname tukeyAdEnet-methods
#' @export
print.tukeyAdEnet <- function(x, ...) {
  cat("Tukey adaptive elastic net fit\n")
  cat(sprintf(
    "  lambda1 = %.4g, lambda2 = %.4g, d = %.4g\n",
    x$lambda1, x$lambda2, x$d
  ))
  cat(sprintf(
    "  non-zero coefficients: %d / %d\n", x$df, length(x$beta)
  ))
  cat(sprintf(
    "  converged: %s (%d iterations)\n", x$converged, x$iterations
  ))
  invisible(x)
}

#' @rdname tukeyAdEnet-methods
#' @export
summary.tukeyAdEnet <- function(object, ...) {
  nz <- object$beta[abs(object$beta) > 0]
  print(object)
  if (length(nz)) {
    cat("\nNon-zero coefficients:\n")
    print(nz)
  } else {
    cat("\nNo non-zero coefficients.\n")
  }
  invisible(object)
}

#' Predict from a fitted Tukey adaptive elastic net
#'
#' @param object A `"tukeyAdEnet"` or `"tukeyAdEnetRBIC"` object, as
#'   returned by [tukeyAdEnet()] or [tukeyAdEnetRBIC()].
#' @param newx Numeric predictor matrix (or an object coercible to one) with
#'   the same number and order of columns used to fit `object`.
#' @param ... Currently ignored.
#'
#' @return A numeric vector of predicted responses, `newx %*% coef(object)`.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, -2, 1.5, rep(0, p - 3))
#' y <- as.numeric(x %*% beta + rnorm(n))
#' fit <- tukeyAdEnet(x, y, lambda1 = 0.2, lambda2 = 0.05)
#' predict(fit, newx = x[1:5, ])
#'
#' @export
predict.tukeyAdEnet <- function(object, newx, ...) {
  if (missing(newx)) stop("`newx` must be supplied.", call. = FALSE)
  newx <- validate_x(newx)
  if (ncol(newx) != length(object$beta)) {
    stop(
      "`newx` must have ", length(object$beta), " columns, matching the fitted model.",
      call. = FALSE
    )
  }
  as.numeric(newx %*% object$beta)
}

#' Inspect a robust-BIC-tuned Tukey adaptive elastic net
#'
#' `coef()` and `predict()` delegate to the selected fit (`object$fit`).
#' `print()` summarises the tuning grid and the selected fit.
#'
#' @param object,x A `"tukeyAdEnetRBIC"` object as returned by
#'   [tukeyAdEnetRBIC()].
#' @param newx Numeric predictor matrix (or an object coercible to one) with
#'   the same number and order of columns used to fit `object`.
#' @param ... Currently ignored.
#'
#' @return `coef()` returns the named numeric coefficient vector of the
#'   selected fit. `predict()` returns a numeric vector of predicted
#'   responses. `print()` returns `object`, invisibly.
#'
#' @name tukeyAdEnetRBIC-methods
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, -2, 1.5, rep(0, p - 3))
#' y <- as.numeric(x %*% beta + rnorm(n))
#' \donttest{
#' tuned <- tukeyAdEnetRBIC(x, y, n_lambda1 = 8, lambda2_factors = c(0, 0.1))
#' tuned
#' coef(tuned)
#' predict(tuned, newx = x[1:5, ])
#' }
NULL

#' @rdname tukeyAdEnetRBIC-methods
#' @export
coef.tukeyAdEnetRBIC <- function(object, ...) {
  coef(object$fit)
}

#' @rdname tukeyAdEnetRBIC-methods
#' @export
predict.tukeyAdEnetRBIC <- function(object, newx, ...) {
  predict(object$fit, newx = newx, ...)
}

#' @rdname tukeyAdEnetRBIC-methods
#' @export
print.tukeyAdEnetRBIC <- function(x, ...) {
  cat("Tukey adaptive elastic net, tuned by robust BIC\n")
  cat(sprintf("  grid size: %d combinations\n", nrow(x$grid)))
  cat(sprintf("  selected rbic: %.4g\n", x$rbic))
  cat("  selected fit:\n")
  print(x$fit)
  invisible(x)
}

#' Plot the robust BIC tuning path
#'
#' Plots the robust BIC criterion against `lambda1`, one line per distinct
#' value of `lambda2`, with the selected `(lambda1, lambda2)` pair marked.
#'
#' @param x A `"tukeyAdEnetRBIC"` object as returned by [tukeyAdEnetRBIC()].
#' @param y Unused; included for S3 method consistency.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return The RBIC grid (invisibly), as a data frame.
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, -2, 1.5, rep(0, p - 3))
#' y <- as.numeric(x %*% beta + rnorm(n))
#' \donttest{
#' tuned <- tukeyAdEnetRBIC(x, y, n_lambda1 = 8, lambda2_factors = c(0, 0.1))
#' plot(tuned)
#' }
#'
#' @export
plot.tukeyAdEnetRBIC <- function(x, y, ...) {
  grid <- x$grid
  lambda2_levels <- sort(unique(grid$lambda2))
  cols <- grDevices::hcl.colors(max(length(lambda2_levels), 2), palette = "Dark 3")

  ord <- order(grid$lambda1)
  grid <- grid[ord, , drop = FALSE]

  plot(
    NA, NA,
    xlim = range(grid$lambda1), ylim = range(grid$rbic),
    log = "x", xlab = expression(lambda[1]), ylab = "Robust BIC",
    main = "Tukey-AdEnet robust BIC tuning path", ...
  )
  for (i in seq_along(lambda2_levels)) {
    sub <- grid[grid$lambda2 == lambda2_levels[i], , drop = FALSE]
    sub <- sub[order(sub$lambda1), , drop = FALSE]
    graphics::lines(sub$lambda1, sub$rbic, col = cols[i], lwd = 2)
  }
  graphics::points(x$fit$lambda1, x$rbic, pch = 19, col = "black", cex = 1.3)
  graphics::legend(
    "topright",
    legend = signif(lambda2_levels, 3),
    col = cols[seq_along(lambda2_levels)],
    lwd = 2, title = expression(lambda[2]), bty = "n", cex = 0.8
  )

  invisible(x$grid)
}
