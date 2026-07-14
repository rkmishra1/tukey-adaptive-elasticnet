selection_metrics <- function(beta_hat, beta_true, zero_tol = 1e-8) {
  true_zero <- abs(beta_true) <= zero_tol
  selected_zero <- abs(beta_hat) <= zero_tol

  c(
    C = sum(true_zero & selected_zero),
    IC = sum(!true_zero & selected_zero)
  )
}

mspe_beta <- function(beta_hat, beta_true, rho) {
  ar1_quadratic_form(beta_hat - beta_true, rho = rho)
}

summarise_comparison_results <- function(raw_results) {
  key <- c("scenario", "regime", "n", "p", "rho", "method")
  split_results <- split(raw_results, raw_results[key], drop = TRUE)

  summaries <- lapply(split_results, function(dat) {
    if (!"status" %in% names(dat)) dat$status <- "ok"
    if (!"error" %in% names(dat)) dat$error <- ""
    if (!"lambda" %in% names(dat)) dat$lambda <- NA_real_
    if (!"lambda1" %in% names(dat)) dat$lambda1 <- NA_real_
    if (!"lambda2" %in% names(dat)) dat$lambda2 <- NA_real_

    dat_ok <- dat[dat$status == "ok" & is.finite(dat$MSPE), , drop = FALSE]
    if (!nrow(dat_ok)) {
      first_error <- dat$error[nzchar(dat$error)][1]
      if (is.na(first_error)) first_error <- ""
      return(data.frame(
        scenario = dat$scenario[1],
        regime = dat$regime[1],
        n = dat$n[1],
        p = dat$p[1],
        rho = dat$rho[1],
        method = dat$method[1],
        requested_reps = nrow(dat),
        completed_reps = 0,
        active = dat$active[1],
        C_mean = NA_real_,
        IC_mean = NA_real_,
        MSPE_mean = NA_real_,
        MSPE_cv = NA_real_,
        lambda_median = NA_real_,
        lambda1_median = NA_real_,
        lambda2_median = NA_real_,
        convergence_rate = NA_real_,
        status = dat$status[1],
        error = first_error,
        stringsAsFactors = FALSE
      ))
    }

    mspe_mean <- mean(dat_ok$MSPE)
    data.frame(
      scenario = dat_ok$scenario[1],
      regime = dat_ok$regime[1],
      n = dat_ok$n[1],
      p = dat_ok$p[1],
      rho = dat_ok$rho[1],
      method = dat_ok$method[1],
      requested_reps = nrow(dat),
      completed_reps = nrow(dat_ok),
      active = dat_ok$active[1],
      C_mean = mean(dat_ok$C),
      IC_mean = mean(dat_ok$IC),
      MSPE_mean = mspe_mean,
      MSPE_cv = stats::sd(dat_ok$MSPE) / mspe_mean,
      lambda_median = stats::median(dat_ok$lambda, na.rm = TRUE),
      lambda1_median = stats::median(dat_ok$lambda1, na.rm = TRUE),
      lambda2_median = stats::median(dat_ok$lambda2, na.rm = TRUE),
      convergence_rate = mean(dat_ok$converged, na.rm = TRUE),
      status = "ok",
      error = "",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summaries)
}

summarise_tukey_adenet_results <- summarise_comparison_results
