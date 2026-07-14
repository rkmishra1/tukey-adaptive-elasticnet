## Test environments

* local Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* `checking for future file timestamps ... NOTE / unable to verify current
  time`: caused by no internet access to the CRAN time-check endpoint in
  this build environment; not expected on machines with normal outbound
  access.

## Downstream dependencies

This is a new submission; there are no reverse dependencies.

## Notes for reviewers

* `robustbase` is used only conditionally (`requireNamespace()`) to build a
  robust initial fit when `p < n`; the package falls back to a ridge
  initial fit otherwise, so `robustbase` is listed under `Suggests`, not
  `Imports`.
