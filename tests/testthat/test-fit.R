sim_data <- function(n = 80, p = 8, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  colnames(x) <- paste0("v", seq_len(p))
  beta <- c(2, -2, 1.5, rep(0, p - 3))
  y <- as.numeric(x %*% beta + stats::rnorm(n))
  list(x = x, y = y, beta = beta)
}

test_that("tukeyAdEnet fits and returns expected structure", {
  d <- sim_data()
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 0.1, lambda2 = 0.05, max_iter = 300)

  expect_s3_class(fit, "tukeyAdEnet")
  expect_length(fit$beta, ncol(d$x))
  expect_equal(names(fit$beta), colnames(d$x))
  expect_true(fit$df <= ncol(d$x))
  expect_true(is.logical(fit$converged))
})

test_that("large lambda1 shrinks all coefficients to zero", {
  d <- sim_data()
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 1e6, lambda2 = 0, max_iter = 200)
  expect_true(all(fit$beta == 0))
  expect_equal(fit$df, 0)
})

test_that("lambda1 = 0, lambda2 = 0 recovers something close to an unpenalised fit", {
  d <- sim_data(n = 300, p = 3, seed = 2)
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 0, lambda2 = 0, max_iter = 2000, tol = 1e-8)
  expect_equal(as.numeric(fit$beta), d$beta, tolerance = 0.5)
})

test_that("coef.tukeyAdEnet returns the beta vector", {
  d <- sim_data()
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 0.1, max_iter = 200)
  expect_identical(coef(fit), fit$beta)
})

test_that("predict.tukeyAdEnet computes newx %*% beta", {
  d <- sim_data()
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 0.1, max_iter = 200)
  pred <- predict(fit, newx = d$x)
  expect_equal(pred, as.numeric(d$x %*% fit$beta))
})

test_that("predict.tukeyAdEnet errors on mismatched newx dimensions", {
  d <- sim_data()
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 0.1, max_iter = 200)
  expect_error(predict(fit, newx = d$x[, 1:3]), "columns")
})

test_that("input validation catches malformed inputs", {
  d <- sim_data()
  expect_error(tukeyAdEnet(d$x, d$y[1:5], lambda1 = 0.1), "length")
  expect_error(tukeyAdEnet(d$x, d$y, lambda1 = -1), "non-negative")

  x_na <- d$x
  x_na[1, 1] <- NA
  expect_error(tukeyAdEnet(x_na, d$y, lambda1 = 0.1), "missing")
})

test_that("works when p > n using ridge initial fit", {
  set.seed(3)
  n <- 20; p <- 30
  x <- matrix(stats::rnorm(n * p), n, p)
  beta <- c(3, -3, rep(0, p - 2))
  y <- as.numeric(x %*% beta + stats::rnorm(n))
  fit <- tukeyAdEnet(x, y, lambda1 = 0.3, lambda2 = 0.1, max_iter = 300)
  expect_length(fit$beta, p)
  expect_true(all(is.finite(fit$beta)))
})

test_that("print and summary methods run without error", {
  d <- sim_data()
  fit <- tukeyAdEnet(d$x, d$y, lambda1 = 0.1, max_iter = 200)
  expect_output(print(fit), "Tukey adaptive elastic net fit")
  expect_output(summary(fit), "coefficients")
})
