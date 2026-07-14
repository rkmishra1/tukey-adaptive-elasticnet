#' tukeyAdEnet: Robust Adaptive Elastic Net Regression via Tukey's Biweight Loss
#'
#' Fits sparse linear regression models that are robust to outliers in the
#' response and in the design matrix by minimising Tukey's redescending
#' biweight loss subject to an adaptive elastic net penalty. Models are
#' fitted with a coordinate-wise proximal AdaGrad algorithm
#' ([tukeyAdEnet()]), and regularisation parameters are selected
#' automatically via a robust BIC grid search ([tukeyAdEnetRBIC()]).
#'
#' @keywords internal
#' @importFrom stats coef predict
"_PACKAGE"
