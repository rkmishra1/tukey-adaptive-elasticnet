test_that("tukeyRho and tukeyPsi are symmetric and bounded", {
  u <- seq(-10, 10, length.out = 201)
  d <- 4.685

  rho <- tukeyRho(u, d = d)
  psi <- tukeyPsi(u, d = d)

  expect_equal(rho, tukeyRho(-u, d = d))
  expect_equal(psi, -tukeyPsi(-u, d = d))

  expect_equal(tukeyRho(0, d = d), 0)
  expect_equal(tukeyPsi(0, d = d), 0)

  expect_true(all(rho <= d^2 / 6 + 1e-10))
  expect_true(all(rho >= 0))
})

test_that("tukeyPsi redescends exactly to zero beyond d", {
  d <- 4.685
  beyond <- c(d + 1e-6, d + 1, 100)
  expect_equal(tukeyPsi(beyond, d = d), c(0, 0, 0))
  expect_equal(tukeyPsi(-beyond, d = d), c(0, 0, 0))
})

test_that("tukeyRho is flat beyond d", {
  d <- 4.685
  expect_equal(tukeyRho(d + 1, d = d), tukeyRho(d + 5, d = d))
})
