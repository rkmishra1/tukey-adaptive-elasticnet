# Data generation for the computational study.

generate_ar1_design <- function(n, p, rho) {
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  x <- z
  if (p > 1) {
    innov_scale <- sqrt(1 - rho^2)
    for (j in 2:p) {
      x[, j] <- rho * x[, j - 1] + innov_scale * z[, j]
    }
  }
  x
}

ar1_quadratic_form <- function(v, rho) {
  p <- length(v)
  if (p == 1) return(v^2)
  idx <- seq_len(p)
  sigma <- rho^abs(outer(idx, idx, "-"))
  as.numeric(crossprod(v, sigma %*% v))
}

generate_beta <- function(p) {
  q <- floor(p / 9)
  s <- 3 * q
  beta <- numeric(p)
  active <- seq_len(s)
  signs <- sample(c(-1, 1), s, replace = TRUE)
  beta[active] <- signs * stats::runif(s, min = 1, max = 3)
  list(beta = beta, active = active, s = s)
}

generate_simulation_data <- function(n,
                                     p,
                                     rho,
                                     scenario = c("clean", "response", "response_design"),
                                     eps = 0.10,
                                     sigma_e = 7,
                                     mu_shift = 15,
                                     mu_x = 10) {
  scenario <- match.arg(scenario)

  beta_info <- generate_beta(p)
  x_clean <- generate_ar1_design(n, p, rho)
  contaminated <- rep(FALSE, n)

  if (scenario == "clean") {
    errors <- stats::rnorm(n, mean = 0, sd = sigma_e)
    x_obs <- x_clean
  } else {
    contaminated <- stats::runif(n) < eps
    errors <- stats::rnorm(n, mean = 0, sd = sigma_e)
    errors[contaminated] <- stats::rnorm(sum(contaminated), mean = mu_shift, sd = sigma_e)

    x_obs <- x_clean
    if (scenario == "response_design" && any(contaminated)) {
      x_obs[contaminated, ] <- matrix(
        stats::rnorm(sum(contaminated) * p, mean = mu_x, sd = 1),
        nrow = sum(contaminated),
        ncol = p
      )
    }
  }

  y <- as.numeric(x_clean %*% beta_info$beta + errors)

  list(
    x = x_obs,
    y = y,
    beta = beta_info$beta,
    active = beta_info$active,
    contaminated = contaminated,
    scenario = scenario,
    n = n,
    p = p,
    rho = rho
  )
}

simulation_grid <- function() {
  regimes <- list(
    zeta_1_2 = data.frame(regime = "zeta_1_2", n = c(500, 1000, 1500), p = c(84, 121, 149)),
    zeta_2_3 = data.frame(regime = "zeta_2_3", n = c(500, 1000, 1500), p = c(246, 394, 519)),
    zeta_5_6 = data.frame(regime = "zeta_5_6", n = c(500, 1000, 1500), p = c(675, 1200, 1800))
  )

  base <- do.call(rbind, regimes)
  scenarios <- c("clean", "response", "response_design")
  rhos <- c(0.30, 0.60, 0.80)

  out <- merge(
    merge(data.frame(scenario = scenarios), base, all = TRUE),
    data.frame(rho = rhos),
    all = TRUE
  )
  out[order(out$scenario, out$regime, out$n, out$rho), ]
}
