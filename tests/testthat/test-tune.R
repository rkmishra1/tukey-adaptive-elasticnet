sim_data <- function(n = 80, p = 8, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  colnames(x) <- paste0("v", seq_len(p))
  beta <- c(2, -2, 1.5, rep(0, p - 3))
  y <- as.numeric(x %*% beta + stats::rnorm(n))
  list(x = x, y = y, beta = beta)
}

test_that("tukeyAdEnetRBIC selects a grid point and returns valid structure", {
  d <- sim_data()
  tuned <- tukeyAdEnetRBIC(
    d$x, d$y,
    n_lambda1 = 5, lambda2_factors = c(0, 0.1), max_iter = 300
  )

  expect_s3_class(tuned, "tukeyAdEnetRBIC")
  expect_s3_class(tuned$fit, "tukeyAdEnet")
  expect_equal(nrow(tuned$grid), 10)
  expect_true(all(c("lambda1", "lambda2", "rbic") %in% names(tuned$grid)))
  expect_equal(tuned$rbic, min(tuned$grid$rbic))
  expect_equal(tuned$fit$lambda1, tuned$grid$lambda1[tuned$best_index])
})

test_that("custom grid is respected", {
  d <- sim_data()
  grid <- expand.grid(lambda1 = c(0.05, 0.2), lambda2 = c(0, 0.1))
  tuned <- tukeyAdEnetRBIC(d$x, d$y, grid = grid, max_iter = 200)
  expect_equal(nrow(tuned$grid), nrow(grid))
})

test_that("malformed custom grid errors informatively", {
  d <- sim_data()
  expect_error(
    tukeyAdEnetRBIC(d$x, d$y, grid = data.frame(a = 1, b = 2)),
    "lambda1"
  )
})

test_that("coef and predict delegate to the selected fit", {
  d <- sim_data()
  tuned <- tukeyAdEnetRBIC(
    d$x, d$y,
    n_lambda1 = 4, lambda2_factors = 0, max_iter = 200
  )
  expect_equal(coef(tuned), coef(tuned$fit))
  expect_equal(predict(tuned, newx = d$x), predict(tuned$fit, newx = d$x))
})

test_that("print.tukeyAdEnetRBIC runs without error", {
  d <- sim_data()
  tuned <- tukeyAdEnetRBIC(
    d$x, d$y,
    n_lambda1 = 4, lambda2_factors = 0, max_iter = 200
  )
  expect_output(print(tuned), "robust BIC")
})

test_that("plot.tukeyAdEnetRBIC runs without error", {
  d <- sim_data()
  tuned <- tukeyAdEnetRBIC(
    d$x, d$y,
    n_lambda1 = 4, lambda2_factors = c(0, 0.1), max_iter = 200
  )
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  expect_silent(plot(tuned))
})
