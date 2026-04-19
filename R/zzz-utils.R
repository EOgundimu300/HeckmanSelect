# -----------------------------------------------------------------------
# Shared utilities for HeckSelect and RidgeBPSS
#
# Small general-purpose helpers used by both fitters.  Defined once
# here so the two fitter files can stay focused on model-specific
# logic.
#
# This file defines (all internal, none exported):
#   - %||% : null-coalesce operator
#   - dbivnorm : bivariate normal density
#   - coerce_binary01 : safe binary response coercion
#   - validate_penalty_arg : argument validator for lambda / lambda.rho
#   - soft_threshold : elementwise lasso soft-thresholding operator
#   - backtransform_coefs : standardised-to-raw coefficient back-transform
#
# The `zzz-` filename prefix is a convention marking this as a
# utility file that loads after other R/ files alphabetically.
# -----------------------------------------------------------------------


# =======================================================================
# NULL-COALESCE OPERATOR
# =======================================================================
#
# Returns the left operand if non-NULL, otherwise the right.  Mirrors
# the behaviour of rlang::`%||%` so users do not need to depend on
# rlang for this one operator.

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# =======================================================================
# BIVARIATE NORMAL DENSITY
# =======================================================================
#
# Used by HeckSelect (via loglik_amh-style computations where the
# joint density is needed as a building block) and by RidgeBPSS (via
# gradlikDn / hesslikDn for the BPSS gradient and Hessian).  A
# minimal self-contained implementation that avoids taking on
# mvtnorm as a dependency.

#' @keywords internal
#' @noRd
dbivnorm <- function(x1, x2, rho) {
  r2 <- 1 - rho^2
  exp(-(x1^2 - 2 * rho * x1 * x2 + x2^2) / (2 * r2)) / (2 * pi * sqrt(r2))
}


# =======================================================================
# BINARY RESPONSE COERCION
# =======================================================================
#
# Safely coerces a response vector to 0/1 integers.  Accepts numeric,
# logical, character, or factor inputs; validates numeric inputs
# BEFORE coercing (as.integer(0.7) silently becomes 0, which would
# pass a post-coercion check but corrupt the data).

#' @keywords internal
#' @noRd
coerce_binary01 <- function(y, name) {
  if (is.character(y)) y <- factor(y)
  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      stop(name, " must be binary (exactly 2 factor levels), got ",
           nlevels(y), " levels: ", paste(levels(y), collapse = ", "))
    }
    # Standard R convention: first level = 0, second level = 1.
    # This matches glm(family = binomial) behaviour.
    y <- as.integer(y == levels(y)[2L])
  } else if (is.logical(y)) {
    y <- as.integer(y)
  } else if (is.numeric(y)) {
    # Validate BEFORE coercing -- as.integer(0.7) gives 0, which would
    # pass a post-coercion check but silently corrupt the data.
    if (!all(is.na(y) | y %in% c(0, 1))) {
      stop(name, " must contain only 0/1 values (and possibly NA for outcome).")
    }
    y <- as.integer(y)
  } else {
    stop(name, " must be numeric, logical, character, or factor.")
  }
  y
}


# =======================================================================
# PENALTY ARGUMENT VALIDATION
# =======================================================================
#
# Validates that lambda / lambda.rho arguments supplied by the user
# are either NULL (automatic tuning), or non-negative finite numeric
# scalars / vectors.  Negative penalties have no statistical meaning
# here and non-finite values would poison the coordinate descent.

#' @keywords internal
#' @noRd
validate_penalty_arg <- function(x, name) {
  if (is.null(x)) return(invisible(NULL))
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0))
    stop(name, " must be NULL or a non-negative finite numeric value.")
  invisible(x)
}


# =======================================================================
# SOFT-THRESHOLDING OPERATOR
# =======================================================================
#
# Elementwise soft-threshold used by the coordinate-descent updates
# in HeckSelect's LSA surrogate:
#   S(z, lambda) = sign(z) * max(|z| - lambda, 0)
# Vectorised over z; lambda may be scalar or vector of matching length.
# This is the unique minimiser of 0.5*(theta - z)^2 + lambda*|theta|.

#' @keywords internal
#' @noRd
soft_threshold <- function(z, lam) {
  sign(z) * pmax(abs(z) - lam, 0)
}


# =======================================================================
# BACK-TRANSFORMATION TO ORIGINAL COVARIATE SCALE
# =======================================================================
#
# When standardize = TRUE, fitters work internally on centred and
# scaled columns.  This function converts the resulting standardised-
# scale coefficients back to the ORIGINAL (raw) covariate scale for
# user-facing output.
#
# Transformation: if column j was standardised with centre mu_j and
# scale sigma_j, then
#   beta_j(raw)  = beta_j(std) / sigma_j    for standardised j
#   beta_j(raw)  = beta_j(std)              otherwise
#   beta_0(raw)  = beta_0(std) - sum_j beta_j(raw) * mu_j
# (sum over standardised j).
#
# The linear predictor X * beta is unchanged; only the coefficient
# parameterisation differs.  Predictions, likelihoods, fitted values
# are all invariant.

#' @keywords internal
#' @noRd
backtransform_coefs <- function(coefs_std, scaling,
                                 equation = c("S", "O")) {
  equation <- match.arg(equation)
  if (is.null(scaling)) return(coefs_std)

  idx    <- if (equation == "S") scaling$idx_S    else scaling$idx_O
  center <- if (equation == "S") scaling$center_S else scaling$center_O
  scale  <- if (equation == "S") scaling$scale_S  else scaling$scale_O

  if (length(idx) == 0L) return(coefs_std)

  # Coefficients in positions idx were fit on (x - mu)/sigma.
  # Divide by sigma to put them back on raw-x scale.
  coefs_out <- coefs_std
  coefs_out[idx] <- coefs_std[idx] / scale

  # Intercept absorbs the centring: -sum(beta_j_raw * mu_j) over
  # standardised j.
  coefs_out[1L] <- coefs_std[1L] - sum(coefs_out[idx] * center)

  coefs_out
}