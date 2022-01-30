

#' Simulated data
#'
# binHeckman was generated from bivariate normal error with correlation 0.5
#' There are 1000 observations and 12 predictors.
#' beta <- c(-2.78,0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.7, 0.7, 0.7)
#' gamma <- c(1.9,0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.7, 0.7, 0.7, 1).

#' Variable x12 is used for exclusion restriction in the selection equation

#'
#' @usage data(binHeckman)
#'
#' @format A data frame with 1000 rows and 15 variables
#' \describe{
#'   \item{yobs}{observed outcome variable}
#'   \item{uu}{selection indicator}
#'   \item{ystar}{underlying outcome variable. It is not needed for the modelling}
#'   \item{x1-x12}{predictor variables}
#'
#' }
#'
#' @keywords datasets
#'
#' @references Ogundimu EO (2022) On Lasso and Adaptive Lasso for non-random
#' sample in credit scoring
#'
#' @source
#'
#' @examples
#' data(binHeckman)
#'

#' \donttest{
#' selection <- uu~ X1+ X2 + X3+ X4+ X5+ X6+ X7+ X8 + X9 + X10 +X11+X12
#' outcome <- yobs ~ X1+ X2 + X3+ X4+ X5+ X6+ X7+ X8 + X9 + X10 +X11
#'pp <- HeckSelect(selection, outcome, data=binHeckman, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")}

#'
"binHeckman"
