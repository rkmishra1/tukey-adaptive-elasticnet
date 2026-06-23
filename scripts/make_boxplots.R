#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    raw = NULL,
    output_dir = "figures"
  )

  for (arg in args) {
    key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(key_value) != 2) next
    if (key_value[1] %in% names(opts)) opts[[key_value[1]]] <- key_value[2]
  }

  if (is.null(opts$raw)) {
    stop("Please provide --raw=path/to/tukey_adenet_raw_*.csv", call. = FALSE)
  }
  opts
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dat <- utils::read.csv(opts$raw)
  if (!dir.exists(opts$output_dir)) dir.create(opts$output_dir, recursive = TRUE)

  groups <- unique(dat[c("scenario", "regime", "n", "p", "rho")])

  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    sub <- dat[
      dat$scenario == g$scenario &
        dat$regime == g$regime &
        dat$n == g$n &
        dat$p == g$p &
        dat$rho == g$rho,
    ]

    file <- file.path(
      opts$output_dir,
      paste0(
        "MSPE_",
        safe_name(g$scenario), "_",
        safe_name(g$regime), "_",
        "n", g$n, "_p", g$p, "_rho", format(g$rho, nsmall = 2),
        ".png"
      )
    )

    grDevices::png(file, width = 1400, height = 900, res = 160)
    graphics::boxplot(
      MSPE ~ method,
      data = sub,
      main = paste0(
        "MSPE: ", g$scenario,
        ", ", g$regime,
        ", n=", g$n,
        ", p=", g$p,
        ", rho=", format(g$rho, nsmall = 2)
      ),
      ylab = "MSPE",
      xlab = "",
      col = "gray85",
      border = "gray30"
    )
    graphics::stripchart(
      MSPE ~ method,
      data = sub,
      vertical = TRUE,
      method = "jitter",
      pch = 16,
      cex = 0.7,
      col = grDevices::adjustcolor("steelblue", alpha.f = 0.45),
      add = TRUE
    )
    grDevices::dev.off()

    message("Wrote: ", file)
  }
}

main()
