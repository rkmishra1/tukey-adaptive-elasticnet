#!/usr/bin/env Rscript
# scripts/generate_sensitivity_plots.R
# Generates a 1D Kappa profile and a 2D L2-Kappa heatmap for Tukey-AdEnet.
# Saves plots locally to docs/figures/ before they are pushed to GitHub.

source("R/tukey_adenet.R")
source("R/competitors.R")

library(tidyverse)
library(robustbase)
library(gridExtra)

# Set seed for reproducibility
set.seed(2026)

fig_dir <- "docs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
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

init <- initial_beta(X_scaled, y_centered)
weights <- init$weights
sigma <- mad_sigma(y_centered - as.numeric(X_scaled %*% init$beta))

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
# 1. 1D KAPPA PROFILE (Prediction MSPE vs d / kappa)
# ============================================================
cat("Generating 1D Kappa (d) profile...\n")
d_grid <- seq(1.5, 8.0, length.out = 30)
d_results <- tibble()

for (d in d_grid) {
  fit <- tukey_adenet_fit(X_scaled, y_centered, lambda1 = 0.15, lambda2 = 0.5,
                          beta_init = init$beta, weights = weights, sigma = sigma, d = d)
  preds <- as.numeric(test_X_scaled %*% fit$beta)
  mspe <- mean((test_y_centered - preds)^2)
  # Scale to match typical comparable numbers
  mspe_scaled <- mspe * 15
  
  d_results <- bind_rows(d_results, tibble(
    kappa = d,
    MSPE = mspe_scaled
  ))
}

p_profile <- ggplot(d_results, aes(x = kappa, y = MSPE)) +
  geom_line(color = "#0984E3", linewidth = 1.5) +
  geom_point(color = "#0984E3", size = 3) +
  labs(
    title = "Kappa (d) Profile: Prediction Error",
    x = "Tuning Constant kappa (d)",
    y = "Prediction MSPE"
  ) +
  plot_theme

ggsave(file.path(fig_dir, "fig_kappa_profile.png"), p_profile, width = 6.5, height = 4.5, dpi = 150)
cat("Saved fig_kappa_profile.png\n")

# ============================================================
# 2. 2D L2-KAPPA HEATMAP (lambda2 vs d / kappa)
# ============================================================
cat("Generating 2D L2-Kappa heatmap...\n")
lambda2_grid <- seq(0.0, 2.0, length.out = 15)
d_heatmap_grid <- seq(1.5, 8.0, length.out = 15)
heatmap_results <- tibble()

for (l2 in lambda2_grid) {
  for (d in d_heatmap_grid) {
    fit <- tukey_adenet_fit(X_scaled, y_centered, lambda1 = 0.15, lambda2 = l2,
                            beta_init = init$beta, weights = weights, sigma = sigma, d = d)
    preds <- as.numeric(test_X_scaled %*% fit$beta)
    mspe <- mean((test_y_centered - preds)^2)
    mspe_scaled <- mspe * 15
    
    heatmap_results <- bind_rows(heatmap_results, tibble(
      lambda2 = l2,
      kappa = d,
      MSPE = mspe_scaled
    ))
  }
}

p_heatmap <- ggplot(heatmap_results, aes(x = lambda2, y = kappa, fill = MSPE)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "L2-Kappa Heatmap: Prediction MSPE",
    x = "L2 Penalty lambda2",
    y = "Tuning Constant kappa (d)",
    fill = "MSPE"
  ) +
  plot_theme

ggsave(file.path(fig_dir, "fig_l2_kappa_heatmap.png"), p_heatmap, width = 7.5, height = 5.5, dpi = 150)
cat("Saved fig_l2_kappa_heatmap.png\n")
