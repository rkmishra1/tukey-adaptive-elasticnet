#!/usr/bin/env Rscript
# scripts/real_data_analysis.R
# Real data analysis comparing Tukey-AdEnet with competitors on three public datasets

# Set paths and source package files
source("R/tukey_adenet.R")
source("R/simulate_data.R")
source("R/metrics.R")
source("R/competitors.R")

library(tidyverse)
library(robustbase)
library(robustHD)
library(rqPen)

# Set random seed for reproducibility
set.seed(42)

# Directory to save figures
fig_dir <- "docs/figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Directory to save raw results
res_dir <- "results"
if (!dir.exists(res_dir)) {
  dir.create(res_dir, recursive = TRUE)
}

# Competitor options (tuning grid and tolerances)
opts <- list(
  glmnet_nlambda = 50,
  adenet_alphas = c(0.1, 0.3, 0.5, 0.7, 0.9),
  n_lambda1 = 20,
  lambda1_min_ratio = 0.02,
  lambda2_factors = c(0, 0.01, 0.05, 0.1, 0.5, 1),
  eta = 0.5,
  tol = 1e-5,
  max_iter = 1000,
  zero_tol = 1e-5,
  lad_nlambda = 30,
  lad_max_iter = 1000,
  sparse_lts_lambda = seq(0.20, 0.05, by = -0.05),
  rlars_s_max = NULL
)

# Helper function to split, scale, and center training/testing sets
prepare_train_test <- function(X, y, train_prop = 0.70) {
  n <- nrow(X)
  train_idx <- sample(1:n, size = floor(train_prop * n))
  
  X_train <- X[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  X_test  <- X[-train_idx, , drop = FALSE]
  y_test  <- y[-train_idx]
  
  # Compute mean and SD of X_train for standardisation
  mean_x <- colMeans(X_train)
  sd_x <- apply(X_train, 2, sd)
  sd_x[sd_x < 1e-10] <- 1  # prevent division by zero
  
  X_train_scaled <- scale(X_train, center = mean_x, scale = sd_x)
  X_test_scaled  <- scale(X_test, center = mean_x, scale = sd_x)
  
  # Center y
  mean_y <- mean(y_train)
  y_train_centered <- y_train - mean_y
  y_test_centered  <- y_test - mean_y
  
  list(
    X_train = X_train_scaled,
    y_train = y_train_centered,
    X_test  = X_test_scaled,
    y_test  = y_test_centered,
    mean_y  = mean_y
  )
}

# Core function to run benchmarks on a dataset
run_dataset_benchmark <- function(X_orig, y_orig, p_noise, B = 20, dataset_name = "") {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("Benchmarking Dataset: %s (n = %d, p_orig = %d, p_noise = %d)\n", 
              dataset_name, nrow(X_orig), ncol(X_orig), p_noise))
  cat(sprintf("==================================================\n"))
  
  methods <- c("AdL", "AdEnet", "LAD-Lasso", "S-LTS", "R-LARS", "Tukey-AdL", "Tukey-AdEnet")
  n_methods <- length(methods)
  
  # Container for results
  results_list <- list()
  
  p_orig <- ncol(X_orig)
  p_total <- p_orig + p_noise
  n_samples <- nrow(X_orig)
  
  for (b in 1:B) {
    cat(sprintf("  Replication %d/%d...\n", b, B))
    
    # 1. Generate noise variables
    X_noise <- matrix(rnorm(n_samples * p_noise), nrow = n_samples, ncol = p_noise)
    colnames(X_noise) <- paste0("Noise", 1:p_noise)
    X_full <- cbind(X_orig, X_noise)
    
    # 2. Split and scale
    split_data <- prepare_train_test(X_full, y_orig, train_prop = 0.70)
    
    for (m in methods) {
      # Fit model with error handling
      fit_res <- tryCatch({
        fit_competitor(m, split_data$X_train, split_data$y_train, opts)
      }, error = function(e) {
        cat(sprintf("    Method %s failed in rep %d: %s\n", m, b, e$message))
        NULL
      })
      
      if (!is.null(fit_res) && length(fit_res$beta) == p_total) {
        beta_hat <- fit_res$beta
        
        # Calculate test predictions
        # Note: y_test is centered, so prediction is X_test %*% beta_hat
        preds <- as.numeric(split_data$X_test %*% beta_hat)
        residuals <- split_data$y_test - preds
        
        # Metrics
        mspe <- mean(residuals^2)
        med_spe <- median(residuals^2)
        
        # Active set selection
        active_vars <- sum(beta_hat != 0)
        signal_sel  <- sum(beta_hat[1:p_orig] != 0)
        noise_sel   <- sum(beta_hat[(p_orig+1):p_total] != 0)
        
        results_list[[length(results_list) + 1]] <- tibble(
          Dataset = dataset_name,
          Rep = b,
          Method = m,
          MSPE = mspe,
          MedSPE = med_spe,
          Active = active_vars,
          SignalSelected = signal_sel,
          NoiseSelected = noise_sel
        )
      }
    }
  }
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # Summarise results
  summary_df <- results_df %>%
    group_by(Method) %>%
    summarise(
      MSPE_mean = mean(MSPE, na.rm = TRUE),
      MSPE_se   = sd(MSPE, na.rm = TRUE) / sqrt(n()),
      MedSPE_mean = mean(MedSPE, na.rm = TRUE),
      MedSPE_se   = sd(MedSPE, na.rm = TRUE) / sqrt(n()),
      Active_mean = mean(Active, na.rm = TRUE),
      Signal_mean = mean(SignalSelected, na.rm = TRUE),
      Noise_mean  = mean(NoiseSelected, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_df)
  return(list(raw = results_df, summary = summary_df))
}

# ==================================================
# DATASET 1: TopGear (robustHD)
# ==================================================
data(TopGear, package = "robustHD")
# Clean and subset to key numeric predictors, and omit missing values
tg_cols <- c("Price", "Displacement", "BHP", "Torque", "Acceleration", "TopSpeed", "Weight", "Length", "Width", "Height")
tg_clean <- na.omit(TopGear[, tg_cols])

X_tg <- as.matrix(tg_clean[, -1])
# Predict Price on log scale due to extreme range
y_tg <- log(tg_clean$Price)

tg_results <- run_dataset_benchmark(X_tg, y_tg, p_noise = 20, B = 20, dataset_name = "TopGear")

# ==================================================
# DATASET 2: pulpfiber (robustbase)
# ==================================================
data(pulpfiber, package = "robustbase")
X_pf <- as.matrix(pulpfiber[, 1:4])
y_pf <- pulpfiber$Y1  # Breaking length

pf_results <- run_dataset_benchmark(X_pf, y_pf, p_noise = 15, B = 20, dataset_name = "pulpfiber")

# ==================================================
# DATASET 3: toxicity (robustbase)
# ==================================================
data(toxicity, package = "robustbase")
X_tx <- as.matrix(toxicity[, -1])
y_tx <- toxicity$toxicity

tx_results <- run_dataset_benchmark(X_tx, y_tx, p_noise = 15, B = 20, dataset_name = "toxicity")

# ==================================================
# Save Results and Generate Figures
# ==================================================
# Save raw results
all_raw <- bind_rows(tg_results$raw, pf_results$raw, tx_results$raw)
write_csv(all_raw, file.path(res_dir, "real_data_raw_results.csv"))

# Save summary tables
all_summaries <- bind_rows(
  tg_results$summary %>% mutate(Dataset = "TopGear"),
  pf_results$summary %>% mutate(Dataset = "pulpfiber"),
  tx_results$summary %>% mutate(Dataset = "toxicity")
)
write_csv(all_summaries, file.path(res_dir, "real_data_summary_results.csv"))

cat("\nGenerating and saving plots...\n")

# Reorder methods for logical comparison
method_order <- c("AdL", "AdEnet", "LAD-Lasso", "S-LTS", "R-LARS", "Tukey-AdL", "Tukey-AdEnet")
all_raw$Method <- factor(all_raw$Method, levels = method_order)

# Set theme
plot_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Function to plot results for a dataset
plot_dataset_metrics <- function(dataset_name, dataset_res, output_path) {
  # We will plot prediction MedSPE (lower is better) and Noise Variables Selected (lower is better)
  summary_data <- dataset_res$summary
  summary_data$Method <- factor(summary_data$Method, levels = method_order)
  
  # 1. MedSPE plot with SE bars
  p_mse <- ggplot(summary_data, aes(x = Method, y = MedSPE_mean, fill = Method)) +
    geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
    geom_errorbar(aes(ymin = MedSPE_mean - MedSPE_se, ymax = MedSPE_mean + MedSPE_se), width = 0.2) +
    scale_fill_manual(values = c(rep("#BDC5C7", 5), "#FF7675", "#0984E3")) +
    labs(
      title = sprintf("%s: Out-of-Sample Prediction Accuracy", dataset_name),
      x = "Method",
      y = "Median Squared Prediction Error (MedSPE)"
    ) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
    
  # 2. Variable selection plot: original vs noise variables
  # Reshape data for stacked bar plot
  select_data <- summary_data %>%
    select(Method, Signal_mean, Noise_mean) %>%
    pivot_longer(cols = c(Signal_mean, Noise_mean), names_to = "VariableType", values_to = "Count") %>%
    mutate(VariableType = factor(VariableType, levels = c("Noise_mean", "Signal_mean"), 
                                 labels = c("Noise Variables (False Positives)", "Signal Variables (True Positives)")))
  
  p_vars <- ggplot(select_data, aes(x = Method, y = Count, fill = VariableType)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("#FF7675", "#55E6C1")) +
    labs(
      title = sprintf("%s: Variable Selection Performance", dataset_name),
      x = "Method",
      y = "Average Number of Selected Variables",
      fill = ""
    ) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
    
  # Combine plots side-by-side
  library(gridExtra)
  combined_plot <- gridExtra::arrangeGrob(p_mse, p_vars, ncol = 2)
  ggsave(output_path, combined_plot, width = 11, height = 5, dpi = 150)
  cat(sprintf("  Saved %s\n", output_path))
}

plot_dataset_metrics("TopGear", tg_results, file.path(fig_dir, "real_data_topgear.png"))
plot_dataset_metrics("pulpfiber", pf_results, file.path(fig_dir, "real_data_pulpfiber.png"))
plot_dataset_metrics("toxicity", tx_results, file.path(fig_dir, "real_data_toxicity.png"))

cat("\nReal data analysis benchmark complete!\n")
