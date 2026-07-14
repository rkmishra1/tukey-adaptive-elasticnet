#!/usr/bin/env Rscript

repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  repos <- c(CRAN = "https://cloud.r-project.org")
}

pkgs <- c("glmnet", "rqPen", "robustHD", "robustbase")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (!length(missing)) {
  message("All required comparison packages are already installed.")
} else {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = repos)
}
