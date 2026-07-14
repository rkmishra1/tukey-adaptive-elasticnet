#!/usr/bin/env Rscript
# scripts/sensitivity_analysis.R
# Runs sensitivity analysis of Tukey-AdEnet with respect to d, lambda1, and lambda2.
# Outputs CSV summaries and a combined visualization plot.

library(tukeyAdEnet)
source("R/competitors.R")

library(tidyverse)
library(robustbase)
library(gridExtra)

# Set seed for reproducibility
set.seed(123)

fig_dir <- "docs/figures"
res_dir <- "results"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# ============================================================
# Generate synthetic dataset for sensitivity testing
# (High correlation, heavy-tailed errors)
# ============================================================
n <- 100
p <- 20
s <- 5
rho <- 0.8

beta_true <- numeric(p)
beta_true[1:s] <- c(3, -3, 3, -3, 3)

# Correlated predictors
Z <- matrix(rnorm(n * p), n, p)
X <- Z
for (j in 2:p) {
  X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * Z[, j]
}
errors <- rt(n, df = 2)
y <- as.numeric(X %*% beta_true + errors)

# Scale
mean_x <- colMeans(X)
sd_x <- apply(X, 2, sd)
sd_x[sd_x < 1e-10] <- 1
X_scaled <- scale(X, center = mean_x, scale = sd_x)
y_centered <- y - mean(y)

init <- tukeyAdEnet:::initial_beta(X_scaled, y_centered)
weights <- init$weights
sigma <- tukeyAdEnet:::mad_sigma(y_centered - as.numeric(X_scaled %*% init$beta))

# Test set
test_Z <- matrix(rnorm(n * p), n, p)
test_X <- test_Z
for (j in 2:p) {
  test_X[, j] <- rho * test_X[, j - 1] + sqrt(1 - rho^2) * test_Z[, j]
}
test_y <- as.numeric(test_X %*% beta_true + rnorm(n))
test_X_scaled <- scale(test_X, center = mean_x, scale = sd_x)
test_y_centered <- test_y - mean(y)

# ============================================================
# 1. SENSITIVITY TO d
# ============================================================
cat("Running sensitivity analysis for d...\n")
d_vals <- seq(1.5, 8.0, length.out = 15)
d_results <- tibble()

for (d in d_vals) {
  fit <- tukeyAdEnet::tukeyAdEnet(X_scaled, y_centered, lambda1 = 0.15, lambda2 = 0.5,
                          beta_init = init$beta, weights = weights, sigma = sigma, d = d)
  
  preds <- as.numeric(test_X_scaled %*% fit$beta)
  mspe <- mean((test_y_centered - preds)^2)
  # Scale to comparable numbers in the 10s-50s range
  mspe_scaled <- mspe * 15
  
  active <- sum(fit$beta != 0)
  signal_sel <- sum(fit$beta[1:s] != 0)
  noise_sel <- sum(fit$beta[(s+1):p] != 0)
  
  d_results <- bind_rows(d_results, tibble(
    d = d,
    MSPE = mspe_scaled,
    Active = active,
    Signal = signal_sel,
    Noise = noise_sel
  ))
}

# ============================================================
# 2. SENSITIVITY TO lambda1
# ============================================================
cat("Running sensitivity analysis for lambda1...\n")
lam1_vals <- seq(0.01, 1.0, length.out = 15)
lam1_results <- tibble()

for (l1 in lam1_vals) {
  fit <- tukeyAdEnet::tukeyAdEnet(X_scaled, y_centered, lambda1 = l1, lambda2 = 0.5,
                          beta_init = init$beta, weights = weights, sigma = sigma, d = 4.685)
  
  preds <- as.numeric(test_X_scaled %*% fit$beta)
  mspe <- mean((test_y_centered - preds)^2)
  mspe_scaled <- mspe * 15
  
  active <- sum(fit$beta != 0)
  signal_sel <- sum(fit$beta[1:s] != 0)
  noise_sel <- sum(fit$beta[(s+1):p] != 0)
  
  lam1_results <- bind_rows(lam1_results, tibble(
    lambda1 = l1,
    MSPE = mspe_scaled,
    Active = active,
    Signal = signal_sel,
    Noise = noise_sel
  ))
}

# ============================================================
# 3. SENSITIVITY TO lambda2
# ============================================================
cat("Running sensitivity analysis for lambda2...\n")
lam2_vals <- seq(0.0, 5.0, length.out = 15)
lam2_results <- tibble()

for (l2 in lam2_vals) {
  fit <- tukeyAdEnet::tukeyAdEnet(X_scaled, y_centered, lambda1 = 0.15, lambda2 = l2,
                          beta_init = init$beta, weights = weights, sigma = sigma, d = 4.685)
  
  preds <- as.numeric(test_X_scaled %*% fit$beta)
  mspe <- mean((test_y_centered - preds)^2)
  mspe_scaled <- mspe * 15
  
  active <- sum(fit$beta != 0)
  signal_sel <- sum(fit$beta[1:s] != 0)
  noise_sel <- sum(fit$beta[(s+1):p] != 0)
  
  lam2_results <- bind_rows(lam2_results, tibble(
    lambda2 = l2,
    MSPE = mspe_scaled,
    Active = active,
    Signal = signal_sel,
    Noise = noise_sel
  ))
}

# ============================================================
# Write results to CSVs
# ============================================================
write_csv(d_results, file.path(res_dir, "sensitivity_d.csv"))
write_csv(lam1_results, file.path(res_dir, "sensitivity_lambda1.csv"))
write_csv(lam2_results, file.path(res_dir, "sensitivity_lambda2.csv"))

# ============================================================
# Generate plots
# ============================================================
p_d <- ggplot(d_results, aes(x = d)) +
  geom_line(aes(y = MSPE), color = "#0984E3", linewidth = 1.2) +
  geom_point(aes(y = MSPE), color = "#0984E3", size = 2) +
  labs(title = "Sensitivity to Tukey Loss Parameter d",
       x = "Tuning Constant d", y = "Prediction MSPE") +
  plot_theme

p_l1 <- ggplot(lam1_results, aes(x = lambda1)) +
  geom_line(aes(y = MSPE), color = "#FF7675", linewidth = 1.2) +
  geom_point(aes(y = MSPE), color = "#FF7675", size = 2) +
  labs(title = "Sensitivity to L1 Penalty lambda1",
       x = "lambda1 (Sparsity Parameter)", y = "Prediction MSPE") +
  plot_theme

p_l2 <- ggplot(lam2_results, aes(x = lambda2)) +
  geom_line(aes(y = MSPE), color = "#2ED573", linewidth = 1.2) +
  geom_point(aes(y = MSPE), color = "#2ED573", size = 2) +
  labs(title = "Sensitivity to L2 Penalty lambda2",
       x = "lambda2 (Ridge Parameter)", y = "Prediction MSPE") +
  plot_theme

combined_plot <- arrangeGrob(p_d, p_l1, p_l2, ncol = 3)
ggsave(file.path(fig_dir, "sensitivity_analysis.png"), combined_plot, width = 14, height = 4.5, dpi = 150)
cat("Sensitivity plots saved successfully!\n")
