#!/usr/bin/env Rscript

source("R/tukey_adenet.R")
source("R/simulate_data.R")
source("R/metrics.R")
source("R/competitors.R")

parse_args <- function(args) {
  opts <- list(
    reps = 200,
    seed = 20260623,
    scenarios = "clean,response,response_design",
    regimes = "zeta_1_2,zeta_2_3,zeta_5_6",
    rhos = "0.30,0.60,0.80",
    n_values = "500,1000,1500",
    methods = "AdL,AdEnet,LAD-Lasso,Tukey-AdL,S-LTS,R-LARS,Tukey-AdEnet",
    missing_action = "skip",
    n_lambda1 = 20,
    lambda1_min_ratio = 0.02,
    lambda2_factors = "0,0.01,0.05,0.1,0.5,1",
    eta = 0.5,
    max_iter = 2000,
    tol = 1e-6,
    zero_tol = 1e-8,
    glmnet_nlambda = 100,
    lad_nlambda = 100,
    lad_max_iter = 5000,
    adenet_alphas = "0.1,0.3,0.5,0.7,0.9",
    sparse_lts_lambda = "0.20,0.15,0.10,0.05",
    rlars_s_max = NA,
    output_dir = "results",
    verbose = FALSE
  )

  for (arg in args) {
    key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(key_value) != 2) next
    key <- key_value[1]
    value <- key_value[2]
    if (!key %in% names(opts)) next
    opts[[key]] <- value
  }

  integer_keys <- c("reps", "seed", "n_lambda1", "max_iter", "glmnet_nlambda", "lad_nlambda", "lad_max_iter")
  numeric_keys <- c("lambda1_min_ratio", "eta", "tol", "zero_tol", "rlars_s_max")
  for (key in integer_keys) opts[[key]] <- as.integer(opts[[key]])
  for (key in numeric_keys) opts[[key]] <- as.numeric(opts[[key]])

  opts$verbose <- as.logical(opts$verbose)
  opts$scenarios <- strsplit(opts$scenarios, ",", fixed = TRUE)[[1]]
  opts$regimes <- strsplit(opts$regimes, ",", fixed = TRUE)[[1]]
  opts$rhos <- as.numeric(strsplit(opts$rhos, ",", fixed = TRUE)[[1]])
  opts$n_values <- as.integer(strsplit(opts$n_values, ",", fixed = TRUE)[[1]])
  opts$methods <- strsplit(opts$methods, ",", fixed = TRUE)[[1]]
  opts$lambda2_factors <- as.numeric(strsplit(opts$lambda2_factors, ",", fixed = TRUE)[[1]])
  opts$adenet_alphas <- as.numeric(strsplit(opts$adenet_alphas, ",", fixed = TRUE)[[1]])
  opts$sparse_lts_lambda <- as.numeric(strsplit(opts$sparse_lts_lambda, ",", fixed = TRUE)[[1]])
  if (is.na(opts$rlars_s_max)) opts$rlars_s_max <- NULL

  unknown_methods <- setdiff(opts$methods, available_methods())
  if (length(unknown_methods)) {
    stop("Unknown methods: ", paste(unknown_methods, collapse = ", "), call. = FALSE)
  }
  if (!opts$missing_action %in% c("skip", "stop")) {
    stop("--missing_action must be 'skip' or 'stop'.", call. = FALSE)
  }

  opts
}

result_row <- function(config, dat, rep_id, method, fit = NULL, status = "ok", error = "") {
  if (status == "ok") {
    beta_hat <- fit$beta
    sel <- selection_metrics(beta_hat, dat$beta)
    mspe <- mspe_beta(beta_hat, dat$beta, rho = config$rho)
    lambda <- fit$lambda
    lambda1 <- fit$lambda1
    lambda2 <- fit$lambda2
    criterion <- fit$criterion
    df <- fit$df
    iterations <- fit$iterations
    converged <- fit$converged
  } else {
    sel <- c(C = NA_real_, IC = NA_real_)
    mspe <- NA_real_
    lambda <- NA_real_
    lambda1 <- NA_real_
    lambda2 <- NA_real_
    criterion <- NA_real_
    df <- NA_real_
    iterations <- NA_integer_
    converged <- NA
  }

  data.frame(
    scenario = config$scenario,
    regime = config$regime,
    n = config$n,
    p = config$p,
    rho = config$rho,
    replication = rep_id,
    method = method,
    active = length(dat$active),
    C = unname(sel["C"]),
    IC = unname(sel["IC"]),
    MSPE = mspe,
    lambda = lambda,
    lambda1 = lambda1,
    lambda2 = lambda2,
    criterion = criterion,
    df = df,
    iterations = iterations,
    converged = converged,
    status = status,
    error = error,
    stringsAsFactors = FALSE
  )
}

run_method <- function(config, dat, rep_id, method, opts) {
  if (!method_is_available(method)) {
    msg <- missing_package_message(method)
    if (opts$missing_action == "stop") stop(msg, call. = FALSE)
    return(result_row(config, dat, rep_id, method, status = "skipped", error = msg))
  }

  fit <- tryCatch(
    fit_competitor(method, dat$x, dat$y, opts),
    error = function(e) {
      if (opts$missing_action == "stop") stop(e)
      e
    }
  )

  if (inherits(fit, "error")) {
    return(result_row(config, dat, rep_id, method, status = "failed", error = conditionMessage(fit)))
  }

  result_row(config, dat, rep_id, method, fit = fit)
}

run_one_replication <- function(config, rep_id, opts) {
  dat <- generate_simulation_data(
    n = config$n,
    p = config$p,
    rho = config$rho,
    scenario = config$scenario
  )

  rows <- lapply(opts$methods, function(method) run_method(config, dat, rep_id, method, opts))
  do.call(rbind, rows)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  set.seed(opts$seed)

  grid <- simulation_grid()
  grid <- subset(
    grid,
    scenario %in% opts$scenarios &
      regime %in% opts$regimes &
      rho %in% opts$rhos &
      n %in% opts$n_values
  )

  if (!dir.exists(opts$output_dir)) dir.create(opts$output_dir, recursive = TRUE)

  raw <- vector("list", nrow(grid) * opts$reps)
  idx <- 1

  for (i in seq_len(nrow(grid))) {
    config <- grid[i, ]
    message(
      "Configuration ", i, "/", nrow(grid), ": ",
      config$scenario, ", ", config$regime,
      ", n=", config$n, ", p=", config$p, ", rho=", config$rho
    )

    for (rep_id in seq_len(opts$reps)) {
      raw[[idx]] <- run_one_replication(config, rep_id, opts)
      if (opts$verbose) {
        ok <- raw[[idx]][raw[[idx]]$status == "ok", , drop = FALSE]
        msg <- if (nrow(ok)) paste(paste(ok$method, signif(ok$MSPE, 5), sep = "="), collapse = ", ") else "no completed methods"
        message("  replication ", rep_id, "/", opts$reps, ": ", msg)
      }
      idx <- idx + 1
    }
  }

  raw_results <- do.call(rbind, raw)
  summary_results <- summarise_comparison_results(raw_results)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  raw_path <- file.path(opts$output_dir, paste0("comparison_raw_", stamp, ".csv"))
  summary_path <- file.path(opts$output_dir, paste0("comparison_summary_", stamp, ".csv"))

  utils::write.csv(raw_results, raw_path, row.names = FALSE)
  utils::write.csv(summary_results, summary_path, row.names = FALSE)

  message("Wrote raw results to: ", raw_path)
  message("Wrote summary results to: ", summary_path)
}

main()
