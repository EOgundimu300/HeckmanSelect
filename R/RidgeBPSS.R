# -----------------------------------------------------------------------
# RidgeBPSS: Ridge-Penalized Bivariate Probit Model with Sample Selection
#
# Main user-facing function: RidgeBPSS()
# Shared bootstrap infrastructure: see R/bpss-shared.R
# Shared utilities: see R/zzz-utils.R
# -----------------------------------------------------------------------

#' @importFrom pbivnorm pbivnorm
#' @importFrom trust trust
#' @importFrom stats glm.fit binomial pnorm dnorm model.frame model.matrix
#'   terms delete.response na.pass model.response coef predict logLik
#'   setNames sd .getXlevels
NULL

# =======================================================================
# CORE BPSS LIKELIHOOD (unpenalized negative log-likelihood)
# Parameter vector: theta = (gammaS, betaO, eta)
# =======================================================================

#' @keywords internal
#' @noRd
loglikDn <- function(beta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS; ibetaO <- seq(tail(ibetaS,1)+1, length=NXO); irho <- tail(ibetaO,1)+1
  betaS <- beta[ibetaS]; betaO <- beta[ibetaO]; rho <- tanh(beta[irho])
  nObs <- length(YS)
  rho_rep <- rep(rho, nObs); k1 <- drop(XO %*% betaO); k2 <- drop(XS %*% betaS)
  # Floor pnorm(-k2) to prevent 0*log(0) = NaN when YS=1 and k2 is large
  pphink2 <- pmax(pnorm(-k2), 1e-16)
  pphi2n <- pmin(pmax(pbivnorm(k2,-k1,-rho_rep),1e-16),1-1e-16)
  pphi2p <- pmin(pmax(pbivnorm(k2, k1, rho_rep),1e-16),1-1e-16)
  loglik <- (1-YS)*log(pphink2) + YS*(1-YO)*log(pphi2n) + YS*YO*log(pphi2p)
  -sum(loglik)
}

#' @keywords internal
#' @noRd
gradlikDn <- function(beta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS; ibetaO <- seq(tail(ibetaS,1)+1,length=NXO); irho <- tail(ibetaO,1)+1
  betaS <- beta[ibetaS]; betaO <- beta[ibetaO]; rho <- tanh(beta[irho])
  r <- sqrt(1-rho^2)
  k1 <- drop(XO%*%betaO); k2 <- drop(XS%*%betaS); phik1 <- dnorm(k1); phik2 <- dnorm(k2)
  omg1 <- (k2-rho*k1)/r; omg2 <- (-k1+rho*k2)/r
  pphiomg1 <- pnorm(omg1); pphinomg2 <- pnorm(-omg2)
  pphink2 <- pmax(pnorm(-k2), 1e-16)
  pphi2n <- pmin(pmax(pbivnorm(k2,-k1,-rho),1e-16),1-1e-16)
  pphi2p <- pmin(pmax(pbivnorm(k2, k1, rho),1e-16),1-1e-16)
  phi2 <- dbivnorm(k2, k1, rho)
  z1n <- (phik1*pphiomg1)/pphi2n; z1p <- (phik1*pphiomg1)/pphi2p
  z2n <- (phik2*pnorm(omg2))/pphi2n; z2p <- (phik2*pphinomg2)/pphi2p
  vn <- phi2/pphi2n; vp <- phi2/pphi2p
  db <- YS*(1-YO)*(-z1n)+YS*YO*z1p
  dg <- (1-YS)*(-phik2/pphink2)+YS*(1-YO)*z2n+YS*YO*z2p
  drhostar <- YS*(1-YO)*(-r^2*vn)+YS*YO*(r^2*vp)
  score <- rep(NA,irho)
  score[ibetaS] <- drop(crossprod(XS,dg)); score[ibetaO] <- drop(crossprod(XO,db)); score[irho] <- sum(drhostar)
  -score
}

#' @keywords internal
#' @noRd
hesslikDn <- function(beta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS; ibetaO <- seq(tail(ibetaS,1)+1,length=NXO); irho <- tail(ibetaO,1)+1; nParam <- irho
  betaS <- beta[ibetaS]; betaO <- beta[ibetaO]; rho <- tanh(beta[irho])
  r <- sqrt(1-rho^2); r1 <- 1/r
  k1 <- drop(XO%*%betaO); k2 <- drop(XS%*%betaS); phik1 <- dnorm(k1); phik2 <- dnorm(k2)
  pphink2 <- pmax(pnorm(-k2), 1e-16)
  omg1 <- (k2-rho*k1)/r; omg2 <- (-k1+rho*k2)/r
  pphi2n <- pmin(pmax(pbivnorm(k2,-k1,-rho),1e-16),1-1e-16)
  pphi2p <- pmin(pmax(pbivnorm(k2, k1, rho),1e-16),1-1e-16)
  phi2 <- dbivnorm(k2, k1, rho)
  z1n <- (phik1*pnorm(omg1))/pphi2n; z1p <- (phik1*pnorm(omg1))/pphi2p
  z2n <- (phik2*pnorm(omg2))/pphi2n; z2p <- (phik2*pnorm(-omg2))/pphi2p
  vn <- phi2/pphi2n; vp <- phi2/pphi2p
  dz1nb <- -k1*z1n-rho*vn+z1n^2; dz1ng <- vn-z1n*z2n; dz1nr <- r^2*vn*(omg2*r1+z1n)
  dz1pb <- -k1*z1p-rho*vp-z1p^2; dz1pg <- vp-z1p*z2p; dz1pr <- r^2*vp*(omg2*r1-z1p)
  dz2ng <- -k2*z2n+rho*vn-z2n^2; dz2nr <- r^2*vn*(omg1*r1+z2n)
  dz2pg <- -k2*z2p-rho*vp-z2p^2; dz2pr <- r^2*vp*(-omg1*r1-z2p)
  w_bb <- YS*(1-YO)*(-dz1nb)+YS*YO*dz1pb; ddlb <- crossprod(XO*as.vector(w_bb),XO)
  w_gg0 <- (1-YS)*((k2*phik2/pphink2)-(phik2/pphink2)^2)
  w_gg <- w_gg0+YS*(1-YO)*dz2ng+YS*YO*dz2pg; ddlg <- crossprod(XS*as.vector(w_gg),XS)
  w_bg <- YS*(1-YO)*(-dz1ng)+YS*YO*dz1pg
  ddlbg <- crossprod(XO*as.vector(w_bg),XS); ddlgb <- crossprod(XS*as.vector(w_bg),XO)
  w_aa <- YS*(1-YO)*(-r^2*vn*(-rho-omg1*omg2+r^2*vn))+YS*YO*(r^2*vp*(-rho-omg1*omg2-r^2*vp))
  ddlrhostar <- sum(w_aa)
  w_ba <- YS*(1-YO)*(-dz1nr)+YS*YO*dz1pr; ddllba <- drop(crossprod(XO,w_ba))
  w_ga <- YS*(1-YO)*dz2nr+YS*YO*dz2pr; ddllga <- drop(crossprod(XS,w_ga))
  hessD <- matrix(NA,nParam,nParam)
  hessD[ibetaO,ibetaO] <- ddlb; hessD[ibetaS,ibetaS] <- ddlg
  hessD[ibetaS,ibetaO] <- ddlgb; hessD[ibetaO,ibetaS] <- ddlbg; hessD[irho,irho] <- ddlrhostar
  hessD[ibetaO,irho] <- ddllba; hessD[ibetaS,irho] <- ddllga
  hessD[irho,ibetaO] <- ddllba; hessD[irho,ibetaS] <- ddllga
  -hessD
}


# =======================================================================
# PENALTY WEIGHTS (internal)
# =======================================================================

#' @keywords internal
#' @noRd
make_penalty_weights <- function(NXS, NXO, lambda, lambda.rho) {
  npar <- NXS + NXO + 1L
  w    <- numeric(npar)
  if (NXS > 1L) w[2:NXS] <- lambda
  if (NXO > 1L) w[(NXS + 2L):(NXS + NXO)] <- lambda
  w[npar] <- lambda.rho
  w
}

# =======================================================================
# DATA PREPARATION (internal)
# =======================================================================

#' @keywords internal
#' @noRd
prepare_bpss_data <- function(selection, outcome, data, standardize = TRUE) {

  # --- Check intercepts are present (BPSS requires them) ---
  stopifnot("Selection formula must include an intercept" =
              attr(terms(selection), "intercept") == 1)
  stopifnot("Outcome formula must include an intercept" =
              attr(terms(outcome), "intercept") == 1)

  # --- Extract model frames ---
  mf_sel <- model.frame(selection, data, na.action = na.pass)
  mf_out <- model.frame(outcome,   data, na.action = na.pass)

  # --- Safely coerce binary responses (validate BEFORE coercing) ---
  YS <- coerce_binary01(model.response(mf_sel), "Selection response (S)")
  YO <- coerce_binary01(model.response(mf_out), "Outcome response (Y_obs)")

  # --- YS must not contain NA ---
  if (anyNA(YS)) stop("Selection response (S) cannot contain missing values.")

  # --- Accepted applicants must have observed outcomes ---
  if (any(is.na(YO[YS == 1L]))) {
    stop("Outcome contains NA among accepted observations (YS = 1). ",
         "All accepted applicants must have observed Y_obs.")
  }

  # --- Mask YO for rejected applicants ONLY ---
  YO[YS == 0L] <- 0L

  # --- Build design matrices ---
  XS <- model.matrix(selection, mf_sel)
  XO <- model.matrix(outcome,   mf_out)

  if (anyNA(XS)) stop("Selection design matrix contains NA values.")
  if (anyNA(XO)) stop("Outcome design matrix contains NA values.")

  # --- Always store factor levels, contrasts, column names for predict() ---
  # This ensures predict() works correctly with factor variables in newdata,
  # even when newdata contains only a subset of factor levels.
  scaling <- list(
    xlev_S      = .getXlevels(terms(selection), mf_sel),
    xlev_O      = .getXlevels(terms(outcome),   mf_out),
    contrasts_S = attr(XS, "contrasts"),
    contrasts_O = attr(XO, "contrasts"),
    colnames_S  = colnames(XS),
    colnames_O  = colnames(XO),
    center_S = NULL, scale_S = NULL, idx_S = integer(0),
    center_O = NULL, scale_O = NULL, idx_O = integer(0),
    standardize = isTRUE(standardize)
  )

  # --- Standardisation rule (glmnet convention) ---
  # When standardize = TRUE (the default):
  #   - Every non-intercept column of the design matrix is centred to mean 0
  #     and scaled by the n-denominator standard deviation
  #        sigma_j = sqrt( sum_i (x_ij - xbar_j)^2 / n ),
  #     matching glmnet's internal convention (NOT base R's n-1 sd()).
  #   - Applies regardless of whether the underlying predictor is numeric,
  #     factor-derived dummy, logical, or interaction.
  #   - The intercept column is never standardised.
  #   - Columns with near-zero spread (e.g. an all-constant dummy in a
  #     subset) are left unscaled (scale = 1) to avoid division by zero.
  #
  # When standardize = FALSE:
  #   - The design matrix is used as-is.
  #
  # In either case, the fitted coefficients returned to the user in
  # $coefficients are back-transformed to the ORIGINAL covariate scale
  # (see backtransform_coefs()).  The internal standardised-scale
  # coefficients remain available in $std_coefficients.  This matches
  # glmnet's user-facing contract: standardisation is an internal
  # optimisation detail, not something the user has to undo.

  if (standardize) {
    # Standardise ALL non-intercept columns uniformly, using the glmnet
    # convention: divide the centred sum of squares by n (not n-1).
    # This corresponds to sigma_j = sqrt(sum_i (x_ij - xbar_j)^2 / n),
    # which is what glmnet uses internally.  The difference from base
    # R's sd() is a factor of sqrt(n/(n-1)) and vanishes for large n.
    idx_S <- if (ncol(XS) > 1L) 2:ncol(XS) else integer(0)
    idx_O <- if (ncol(XO) > 1L) 2:ncol(XO) else integer(0)

    n <- nrow(XS)

    glmnet_standardize <- function(X, cols) {
      Xsub <- X[, cols, drop = FALSE]
      centre <- colMeans(Xsub)
      Xc <- sweep(Xsub, 2L, centre, "-")
      scale <- sqrt(colSums(Xc^2) / nrow(Xc))   # n-denominator
      scale[scale < 1e-10] <- 1
      Xscaled <- sweep(Xc, 2L, scale, "/")
      list(Xscaled = Xscaled, centre = centre, scale = scale)
    }

    if (length(idx_S) > 0L) {
      out <- glmnet_standardize(XS, idx_S)
      XS[, idx_S] <- out$Xscaled
      scaling$center_S <- out$centre
      scaling$scale_S  <- out$scale
    }

    if (length(idx_O) > 0L) {
      out <- glmnet_standardize(XO, idx_O)
      XO[, idx_O] <- out$Xscaled
      scaling$center_O <- out$centre
      scaling$scale_O  <- out$scale
    }

    scaling$idx_S <- idx_S
    scaling$idx_O <- idx_O

    # Clean attributes left by scale()
    attr(XS, "dimnames") <- list(NULL, colnames(XS))
    attr(XO, "dimnames") <- list(NULL, colnames(XO))
  }

  list(YS = YS, YO = YO, XS = XS, XO = XO,
       n = length(YS), NXS = ncol(XS), NXO = ncol(XO),
       scaling = scaling, formula.sel = selection, formula.out = outcome)
}


# =======================================================================
# STARTING VALUES (internal)
# =======================================================================

#' @keywords internal
#' @noRd
make_start_values <- function(prep) {
  gamma_start <- tryCatch({
    coefs <- suppressWarnings(
      glm.fit(x = prep$XS, y = prep$YS, family = binomial(link = "probit"))
    )$coefficients
    coefs[!is.finite(coefs)] <- 0; coefs
  }, error = function(e) rep(0, prep$NXS))

  acc <- (prep$YS == 1L)
  beta_start <- tryCatch({
    if (sum(acc) > prep$NXO + 5L) {
      coefs <- suppressWarnings(
        glm.fit(x = prep$XO[acc, , drop = FALSE], y = prep$YO[acc],
                family = binomial(link = "probit"))
      )$coefficients
      coefs[!is.finite(coefs)] <- 0; coefs
    } else { rep(0, prep$NXO) }
  }, error = function(e) rep(0, prep$NXO))

  c(gamma_start, beta_start, 0)
}


# =======================================================================
# IC COMPUTATION (internal)
#
# IMPORTANT: renamed from `compute_ic` to `compute_ic_ridge` 
# =======================================================================

#' @keywords internal
#' @noRd
compute_ic_ridge <- function(theta_hat, prep, pen_wts, crit = "bic") {
  npar <- length(theta_hat); n <- prep$n
  neg_loglik <- loglikDn(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO)

  if (all(pen_wts == 0)) {
    df <- npar
  } else {
    H_unpen <- hesslikDn(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO)
    H_unpen <- 0.5 * (H_unpen + t(H_unpen))
    H_pen   <- H_unpen + 2 * diag(pen_wts)

    df <- tryCatch({
      # df = tr[H_unpen * H_pen^{-1}] = tr[H_pen^{-1} * H_unpen]
      # Using solve(A, B) avoids explicit matrix inversion
      sum(diag(solve(H_pen, H_unpen)))
    }, error = function(e) {
      warning("penalized Hessian could not be inverted; ",
              "falling back to df = npar for IC computation.")
      npar
    })
    df <- max(0, min(df, npar))
  }

  ic <- switch(crit,
    bic = 2 * neg_loglik + df * log(n),
    aic = 2 * neg_loglik + 2 * df,
    gcv = { -2 * (-neg_loglik) / (n * (1 - df / n)^2) },
    stop("Unknown criterion: ", crit)
  )
  list(ic = ic, df = df, neg_loglik = neg_loglik)
}

# =======================================================================
# SINGLE FIT AT FIXED (lambda, lambda.rho)
# =======================================================================

#' @keywords internal
#' @noRd
fit_single <- function(prep, lambda, lambda.rho, start = NULL,
                       crit = "bic",
                       rinit = 1, rmax = 100, iterlim = 500) {

  pen_wts <- make_penalty_weights(prep$NXS, prep$NXO, lambda, lambda.rho)
  npar <- prep$NXS + prep$NXO + 1L

  if (is.null(start)) { theta0 <- make_start_values(prep) }
  else                 { theta0 <- start }
  stopifnot(length(theta0) == npar)

  # Counter for diagnostic purposes: how many times did the objfun fallback
  # trigger during this fit?  Non-zero values flag that the optimizer traversed
  # regions where the objective was non-finite and had to be rescued.
  fallback_hits <- 0L

  objfun <- function(theta) {
    val  <- loglikDn(theta,  prep$YS, prep$XS, prep$YO, prep$XO)
    grad <- gradlikDn(theta, prep$YS, prep$XS, prep$YO, prep$XO)
    hess <- hesslikDn(theta, prep$YS, prep$XS, prep$YO, prep$XO)
    val  <- val  + sum(pen_wts * theta^2)
    grad <- grad + 2 * pen_wts * theta
    diag(hess) <- diag(hess) + 2 * pen_wts
    # Defensive check: if any quantity is non-finite (e.g. from extreme eta
    # pushing rho near ±1), return a large value with zero gradient and
    # identity Hessian so the trust region shrinks rather than crashes.
    if (!is.finite(val) || any(!is.finite(grad)) || any(!is.finite(hess))) {
      fallback_hits <<- fallback_hits + 1L
      val  <- 1e20
      grad <- rep(0, length(theta))
      hess <- diag(length(theta))
    }
    list(value = val, gradient = grad, hessian = hess)
  }

  fit <- tryCatch(
    trust(objfun, parinit = theta0, rinit = rinit, rmax = rmax, iterlim = iterlim),
    error = function(e) {
      list(argument = theta0, converged = FALSE, iterations = 0,
           hessian = matrix(NA_real_, npar, npar))
    }
  )

  # --- Post-fit sanity check: reject fallback pseudo-convergence ---
  # The fallback in objfun() returns val = 1e20 with grad = 0 and hess = I
  # whenever any of val/grad/hess is non-finite.  A trust region method can
  # mistake this for a stationary point.  Crucially, the fallback can trigger
  # when only the gradient or Hessian is non-finite even if the objective
  # value itself is finite (e.g., r = sqrt(1-rho^2) near zero produces finite
  # probabilities but divergent derivative terms via divisions by r).
  # We therefore re-evaluate ALL THREE true penalized quantities at the
  # proposed solution; if any is non-finite (or the objective sits near the
  # fallback value), the fit is not a genuine optimum.
  #
  # Side-effect: if the sanity check passes, we keep the freshly recomputed
  # penalized Hessian (th_recomputed) and use it downstream for cond_number
  # and storage, rather than fit$hessian from trust().  This guarantees we
  # never store the identity-matrix fallback Hessian on a converged fit.
  th_recomputed <- NULL
  if (isTRUE(fit$converged)) {
    sane <- tryCatch({
      theta_test <- fit$argument
      tv  <- loglikDn(theta_test,  prep$YS, prep$XS, prep$YO, prep$XO) +
               sum(pen_wts * theta_test^2)
      tg  <- gradlikDn(theta_test, prep$YS, prep$XS, prep$YO, prep$XO) +
               2 * pen_wts * theta_test
      th  <- hesslikDn(theta_test, prep$YS, prep$XS, prep$YO, prep$XO)
      diag(th) <- diag(th) + 2 * pen_wts
      ok <- is.finite(tv) && tv < 1e19 &&
              all(is.finite(tg)) && all(is.finite(th))
      if (ok) th_recomputed <<- th
      ok
    }, error = function(e) FALSE)
    if (!isTRUE(sane)) {
      fit$converged <- FALSE
    }
  }

  # --- Early return on failure ---
  # Critically, we do NOT expose theta0 as if it were a fitted parameter
  # vector.  $theta, $gamma, $beta, $eta, $rho are all set to NA so that
  # any caller inspecting a failed fit cannot mistake starting values
  # for estimates.  The starting vector is preserved separately as
  # $start_values for debugging purposes.
  if (!isTRUE(fit$converged)) {
    na_gamma <- setNames(rep(NA_real_, prep$NXS), colnames(prep$XS))
    na_beta  <- setNames(rep(NA_real_, prep$NXO), colnames(prep$XO))
    return(list(
      theta = rep(NA_real_, npar),
      gamma = na_gamma,
      beta  = na_beta,
      eta = NA_real_, rho = NA_real_,
      lambda = lambda, lambda.rho = lambda.rho,
      ic = Inf, df = NA_real_, neg_loglik = NA_real_, logLik = NA_real_,
      converged = FALSE, iterations = if(is.numeric(fit$iterations)) fit$iterations else 0L,
      pen_wts = pen_wts, hessian = matrix(NA_real_, npar, npar),
      cond_number = NA_real_, crit = crit,
      fallback_hits = fallback_hits,
      start_values = theta0
    ))
  }

  theta_hat <- fit$argument
  eta_hat <- theta_hat[npar]; rho_hat <- tanh(eta_hat)
  gamma_hat <- setNames(theta_hat[1:prep$NXS], colnames(prep$XS))
  beta_hat  <- setNames(theta_hat[(prep$NXS+1):(prep$NXS+prep$NXO)], colnames(prep$XO))

  ic_result <- tryCatch(
    compute_ic_ridge(theta_hat, prep, pen_wts, crit),
    error = function(e) list(ic = Inf, df = NA_real_, neg_loglik = NA_real_)
  )

  # Prefer the freshly recomputed penalized Hessian (guaranteed finite, never
  # the fallback I_n).  Fall back to trust()'s stored Hessian only if the
  # recomputation path was not taken — which cannot happen for a converged
  # fit past the sanity check, but we keep the guard for robustness.
  H_Q <- if (!is.null(th_recomputed)) th_recomputed else fit$hessian
  cond_number <- tryCatch({
    H_sym <- 0.5 * (H_Q + t(H_Q))
    sv <- svd(H_sym, nu = 0, nv = 0)$d
    max(sv) / max(tail(sv, 1), 1e-15)
  }, error = function(e) NA_real_)

  list(theta = theta_hat, gamma = gamma_hat, beta = beta_hat,
       eta = eta_hat, rho = rho_hat,
       lambda = lambda, lambda.rho = lambda.rho,
       ic = ic_result$ic, df = ic_result$df,
       neg_loglik = ic_result$neg_loglik, logLik = -ic_result$neg_loglik,
       converged = TRUE, iterations = fit$iterations,
       pen_wts = pen_wts, hessian = H_Q,
       cond_number = cond_number, crit = crit,
       fallback_hits = fallback_hits)
}


# =======================================================================
# GRID SEARCH WITH WARM STARTS (internal)
# =======================================================================

#' @keywords internal
#' @noRd
tune_grid <- function(prep, lambda.grid, lambda.rho.grid,
                      start = NULL, crit = "bic",
                      rinit = 1, rmax = 100, iterlim = 500) {

  # Sort DESCENDING: large lambda (well-conditioned) first for stable warm starts
  lambda.grid     <- sort(lambda.grid, decreasing = TRUE)
  lambda.rho.grid <- sort(lambda.rho.grid, decreasing = TRUE)

  best_ic <- Inf; best_fit <- NULL
  grid_results <- list(); counter <- 0L

  if (is.null(start)) start <- make_start_values(prep)

  for (lr in lambda.rho.grid) {
    warm_start <- start
    for (l in lambda.grid) {
      counter <- counter + 1L
      fit <- fit_single(prep, l, lr, start = warm_start, crit = crit,
                        rinit = rinit, rmax = rmax, iterlim = iterlim)
      if (fit$converged) warm_start <- fit$theta
      if (fit$converged && is.finite(fit$ic) && fit$ic < best_ic) {
        best_ic <- fit$ic; best_fit <- fit
      }
      grid_results[[counter]] <- list(
        lambda = l, lambda.rho = lr, ic = fit$ic,
        df = fit$df, rho = fit$rho, converged = fit$converged
      )
    }
  }
  grid_df <- do.call(rbind, lapply(grid_results, as.data.frame))
  list(best_fit = best_fit, grid = grid_df)
}


# =======================================================================
# MAIN USER-FACING FUNCTION
# =======================================================================

#' Ridge-Penalized Bivariate Probit Model with Sample Selection
#'
#' Fits a ridge-penalized bivariate probit with sample selection model
#' using trust-region optimisation of the penalized negative
#' log-likelihood.  Supports automatic tuning of the penalty
#' parameters via BIC, AIC, or GCV.
#'
#' @param selection a formula for the selection equation, including an
#'   intercept.
#' @param outcome a formula for the outcome equation, including an
#'   intercept.
#' @param data a data frame containing all variables referenced in the
#'   two formulas, including the binary selection and outcome responses.
#' @param lambda numeric scalar, vector, or NULL.  If NULL (default),
#'   the tuning grid `control$lambda.grid` is used with IC tuning.
#'   If scalar, that value is used.  If a vector, those values are used
#'   as the grid.
#' @param lambda.rho numeric scalar, vector, or NULL.  Same semantics
#'   as `lambda`, applied to eta = atanh(rho).  Default NULL.
#' @param standardize logical; if TRUE (default), design matrix columns
#'   are centred and scaled internally (glmnet convention) and
#'   coefficients are returned on the raw scale.
#' @param crit information criterion for tuning: "bic" (default), "aic",
#'   or "gcv".
#' @param start optional starting values of length `NXS + NXO + 1`.
#' @param control an optional list of control arguments:
#'   \describe{
#'     \item{rinit, rmax, iterlim}{trust-region settings (defaults 1,
#'           100, 500).}
#'     \item{lambda.grid, lambda.rho.grid}{tuning grids when the
#'           corresponding argument is NULL.}
#'   }
#' @param ... currently ignored.
#'
#' @return an object of class "RidgeBPSS" with components including
#'   `coefficients`, `std_coefficients`, `rho`, `lambda`, `lambda.rho`,
#'   `logLik`, `ic`, `df`, `converged`, `fitted.values`, `grid`, and
#'   `fit_spec`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500
#' x1 <- rnorm(n); x2 <- rnorm(n); Z <- rnorm(n); V <- rnorm(n)
#' S <- as.integer(0.2 + 0.5*x1 - 0.4*x2 - 0.5*Z + rnorm(n) > 0)
#' Y <- as.integer(0.1 + 0.6*x1 - 0.3*x2 + rnorm(n) > 0)
#' Y[S == 0] <- 0L
#' dat <- data.frame(S = S, Y_obs = Y, x1 = x1, x2 = x2, Z = Z, V = V)
#' fit <- RidgeBPSS(S ~ x1 + x2 + Z, Y_obs ~ x1 + x2 + V, dat,
#'                    lambda = 0.1, lambda.rho = 0)
#' print(fit)
#' summary(fit)
#' }
#'
#' @seealso [HeckSelect()] for the lasso / adaptive-lasso variant,
#'   [bootValidate()] for bootstrap optimism correction.
#' @export
RidgeBPSS <- function(selection, outcome,
                      data = sys.frame(sys.parent()),
                      lambda = NULL, lambda.rho = NULL,
                      standardize = TRUE,
                      crit = c("bic", "aic", "gcv"),
                      start = NULL,
                      control = list(), ...) {

  cl <- match.call()

  # Capture the ORIGINAL user-supplied fitting specification BEFORE any
  # argument defaulting.  This is used by bootValidate.RidgeBPSS() to
  # replay the original fitting procedure on each bootstrap sample.
  # Storing the INPUT (not the selected) values is essential so that a
  # user who auto-tuned gets bootstrap refits that also auto-tune, and
  # a user who pinned (lambda, lambda.rho) gets refits that also pin.
  fit_spec <- list(
    lambda_input     = lambda,          # NULL, scalar, or vector
    lambda_rho_input = lambda.rho,      # NULL, scalar, or vector
    control_input    = control,         # user's control list (un-defaulted)
    crit             = match.arg(crit),
    standardize      = standardize
  )

  crit <- fit_spec$crit

  validate_penalty_arg(lambda, "lambda")
  validate_penalty_arg(lambda.rho, "lambda.rho")

  con <- list(rinit = 1, rmax = 100, iterlim = 500,
              lambda.grid     = c(0, 10^seq(-4, 2, length.out = 10)),
              lambda.rho.grid = c(0, 10^seq(-4, 2, length.out = 10)))
  con[names(control)] <- control
  stopifnot("lambda.grid must be non-negative" = all(con$lambda.grid >= 0))
  stopifnot("lambda.rho.grid must be non-negative" = all(con$lambda.rho.grid >= 0))

  prep <- prepare_bpss_data(selection, outcome, data, standardize)

  if (!is.null(start)) {
    expected <- prep$NXS + prep$NXO + 1L
    if (length(start) != expected) stop("start must have length ", expected)
  }

  if (!is.null(lambda) && !is.null(lambda.rho)) {
    if (length(lambda) != 1L) stop("For a single fit, lambda must be scalar.")
    if (length(lambda.rho) != 1L) stop("For a single fit, lambda.rho must be scalar.")
  }

  # --- Dispatch ---
  if (is.null(lambda) && is.null(lambda.rho)) {
    tr <- tune_grid(prep, con$lambda.grid, con$lambda.rho.grid,
                    start, crit, con$rinit, con$rmax, con$iterlim)
    best <- tr$best_fit; grid <- tr$grid
  } else if (is.null(lambda)) {
    tr <- tune_grid(prep, con$lambda.grid, lambda.rho,
                    start, crit, con$rinit, con$rmax, con$iterlim)
    best <- tr$best_fit; grid <- tr$grid
  } else if (is.null(lambda.rho)) {
    tr <- tune_grid(prep, lambda, con$lambda.rho.grid,
                    start, crit, con$rinit, con$rmax, con$iterlim)
    best <- tr$best_fit; grid <- tr$grid
  } else {
    best <- fit_single(prep, lambda, lambda.rho, start, crit,
                       con$rinit, con$rmax, con$iterlim)
    grid <- data.frame(lambda = lambda, lambda.rho = lambda.rho,
                       ic = best$ic, df = best$df,
                       rho = best$rho, converged = best$converged)
  }

  # Helper to build the "failed fit" return object.  Used both when the
  # tuning grid produced no converged fit AND when a single fixed-penalty
  # fit_single() call returned converged = FALSE.  Critically, this
  # function does NOT expose the starting values of theta as if they
  # were fitted coefficients: $coefficients, $linear.predictors, and
  # $fitted.values are NA.  A copy of the start vector (if available)
  # is preserved under $start_values for debugging.
  build_failed_object <- function(grid_df, best_obj = NULL) {
    structure(list(
      call = cl, converged = FALSE,
      coefficients = list(selection = NA, outcome = NA),
      std_coefficients = list(selection = NA, outcome = NA),
      rho = NA, eta = NA,
      lambda = if (!is.null(best_obj)) best_obj$lambda else NA,
      lambda.rho = if (!is.null(best_obj)) best_obj$lambda.rho else NA,
      ic = Inf, df = NA, logLik = NA, crit = crit,
      iterations = if (!is.null(best_obj)) best_obj$iterations else NA,
      fitted.values = NA,
      linear.predictors = list(outcome = NA, selection = NA),
      hessian = NA, condition = NA,
      formula = list(selection = selection, outcome = outcome),
      scaling = prep$scaling, nobs = prep$n, npar = NA,
      theta = NA, pen_wts = NA, grid = grid_df, prep = prep,
      start_values = if (!is.null(best_obj)) {
        if (!is.null(best_obj$start_values)) best_obj$start_values
        else best_obj$theta
      } else NULL,
      fallback_hits = if (!is.null(best_obj)) best_obj$fallback_hits else NA_integer_,
      fit_spec = fit_spec
    ), class = "RidgeBPSS")
  }

  # --- Guard against total failure (tuning grid produced nothing) ---
  if (is.null(best)) {
    warning("No converged fit found across the tuning grid.")
    return(build_failed_object(grid_df = grid, best_obj = NULL))
  }

  # --- Guard against failed single fit ---
  # If the user requested a single fixed-penalty fit and it did not
  # converge, do NOT build a fitted object from theta0.  Return the
  # same NA skeleton as the tuning-grid failure case.
  if (!isTRUE(best$converged)) {
    warning("Single fit did not converge; returning NA-populated failed-fit object.")
    return(build_failed_object(grid_df = grid, best_obj = best))
  }

  k1 <- drop(prep$XO %*% best$beta)
  k2 <- drop(prep$XS %*% best$gamma)
  rho_hat <- best$rho
  prob_s <- pnorm(k2)
  prob_joint <- pbivnorm(k1, k2, rho_hat)
  prob_y_given_s1 <- ifelse(prob_s > 1e-8,
                            pmin(pmax(prob_joint / prob_s, 0), 1), NA_real_)

  # --- Back-transform coefficients to the ORIGINAL covariate scale ---
  # This ensures user-facing $coefficients are directly comparable to the
  # true data-generating coefficients (e.g. in Monte Carlo bias/RMSE).
  # The internal $theta remains on the standardised scale for reproducibility
  # and for predict(), which re-applies the stored centres/scales to newdata.
  gamma_raw <- backtransform_coefs(best$gamma, prep$scaling, equation = "S")
  beta_raw  <- backtransform_coefs(best$beta,  prep$scaling, equation = "O")

  structure(list(
    call = cl,
    coefficients = list(selection = gamma_raw, outcome = beta_raw),
    std_coefficients = list(selection = best$gamma, outcome = best$beta),
    eta = best$eta, rho = best$rho,
    lambda = best$lambda, lambda.rho = best$lambda.rho,
    logLik = best$logLik, ic = best$ic, df = best$df, crit = crit,
    converged = best$converged, iterations = best$iterations,
    fitted.values = prob_y_given_s1,
    linear.predictors = list(outcome = k1, selection = k2),
    hessian = best$hessian, condition = best$cond_number,
    formula = list(selection = selection, outcome = outcome),
    scaling = prep$scaling, nobs = prep$n,
    npar = length(best$theta), theta = best$theta,
    pen_wts = best$pen_wts, grid = grid, prep = prep,
    fallback_hits = best$fallback_hits,
    fit_spec = fit_spec
  ), class = "RidgeBPSS")
}


# =======================================================================
# S3 METHODS
# =======================================================================

#' Print a RidgeBPSS fit
#' @param x a RidgeBPSS fit.
#' @param ... ignored.
#' @return invisibly returns the fit.
#' @export
print.RidgeBPSS <- function(x, ...) {
  cat("Ridge-penalized Bivariate Probit with Sample Selection\n\nCall:\n")
  print(x$call); cat("\n")

  # Helper: print the fallback_hits diagnostic if available and non-zero.
  # This is most informative on failed fits, so we print it regardless
  # of convergence status.
  print_fallback <- function() {
    if (!is.null(x$fallback_hits) && !is.na(x$fallback_hits) && x$fallback_hits > 0L) {
      cat(sprintf("Fallback hits:      %d (non-finite objfun evaluations during optimisation)\n",
                  x$fallback_hits))
    }
  }

  if (!isTRUE(x$converged)) {
    cat("Model did NOT converge.\n")
    print_fallback()
    return(invisible(x))
  }
  cat(sprintf("rho (selectivity):  %.4f  (eta = %.4f)\n", x$rho, x$eta))
  cat(sprintf("lambda:             %.6f\n", x$lambda))
  cat(sprintf("lambda.rho:         %.6f\n", x$lambda.rho))
  cat(sprintf("Log-likelihood:     %.2f\n", x$logLik))
  cat(sprintf("%s:                %.2f  (df = %.1f)\n", toupper(x$crit), x$ic, x$df))
  cat(sprintf("Converged:          %s  (%d iterations)\n", x$converged, x$iterations))
  cat(sprintf("Condition number:   %.1f\n", x$condition))
  cat(sprintf("N observations:     %d\n", x$nobs))
  print_fallback()
  cat("\n")
  invisible(x)
}

#' Summarise a RidgeBPSS fit
#' @param object a RidgeBPSS fit.
#' @param ... ignored.
#' @return invisibly returns the fit.
#' @export
summary.RidgeBPSS <- function(object, ...) {
  if (!isTRUE(object$converged)) { cat("Model did not converge.\n"); return(invisible(object)) }
  cat("Ridge-penalized BPSS\n====================\n\nCall:\n"); print(object$call); cat("\n")
  cat(sprintf("Observations:  %d\n", object$nobs))
  cat(sprintf("Log-lik:       %.4f\n", object$logLik))
  cat(sprintf("%s:           %.4f  (df = %.2f)\n", toupper(object$crit), object$ic, object$df))
  cat(sprintf("rho:           %.4f  (eta = %.4f)\n", object$rho, object$eta))
  cat(sprintf("lambda:        %.6f\n", object$lambda))
  cat(sprintf("lambda.rho:    %.6f\n", object$lambda.rho))
  cat(sprintf("Condition:     %.1f\n\n", object$condition))
  if (isTRUE(object$scaling$standardize)) {
    cat("NOTE: Internal standardisation was applied (glmnet convention).\n")
    cat("      Coefficients shown below are on the ORIGINAL covariate scale\n")
    cat("      (back-transformed).  Standardised-scale coefficients are\n")
    cat("      available in $std_coefficients.\n\n")
  } else {
    cat("NOTE: Fitted without internal standardisation (standardize = FALSE).\n\n")
  }
  cat("Selection equation coefficients:\n")
  print(round(object$coefficients$selection, 5))
  cat("\nOutcome equation coefficients:\n")
  print(round(object$coefficients$outcome, 5))
  cat(sprintf("\nDependence parameter:\n  eta = %.5f    rho = tanh(eta) = %.5f\n",
              object$eta, object$rho))
  invisible(object)
}

#' Extract coefficients from a RidgeBPSS fit
#' @param object a RidgeBPSS fit.
#' @param type \code{"all"} (default), \code{"selection"},
#'   \code{"outcome"}, or \code{"rho"}.
#' @param ... ignored.
#' @return numeric vector or list, depending on `type`.
#' @export
coef.RidgeBPSS <- function(object, type = c("all", "selection", "outcome", "rho"), ...) {
  type <- match.arg(type)
  switch(type,
    all       = list(selection = object$coefficients$selection,
                     outcome = object$coefficients$outcome,
                     rho = object$rho, eta = object$eta),
    selection = object$coefficients$selection,
    outcome   = object$coefficients$outcome,
    rho       = c(rho = object$rho, eta = object$eta))
}

#' Predict from a RidgeBPSS fit
#' @param object a RidgeBPSS fit.
#' @param newdata optional new data frame; if NULL, returns in-sample
#'   predictions.
#' @param type one of \code{"conditional"}, \code{"marginal"}, or
#'   \code{"selection"}.
#' @param ... ignored.
#' @return numeric vector of predicted probabilities.
#' @export
predict.RidgeBPSS <- function(object, newdata = NULL,
                               type = c("conditional", "marginal", "selection"), ...) {
  type <- match.arg(type)

  # --- Failed-fit handling ---
  # For the in-sample case (newdata = NULL), the stored fitted.values,
  # linear.predictors, etc. are already NA on a failed fit, so we just
  # return them (with a warning).  For newdata, however, we cannot build
  # predictions from NA coefficients, so stop early with a clear message
  # rather than warning and then failing obscurely later.
  if (!isTRUE(object$converged)) {
    if (!is.null(newdata)) {
      stop("Model did not converge; prediction on newdata is unavailable. ",
           "Check $converged and refit with different starting values, ",
           "penalty, or trust-region settings.")
    }
    warning("Model did not converge; returning NA predictions.")
  }

  if (is.null(newdata)) {
    return(switch(type,
      conditional = object$fitted.values,
      marginal = pnorm(object$linear.predictors$outcome),
      selection = pnorm(object$linear.predictors$selection)))
  }
  tt_sel <- delete.response(terms(object$formula$selection))
  tt_out <- delete.response(terms(object$formula$outcome))

  # Use stored factor levels and contrasts from training data.
  # This ensures correct dummy-variable structure even when newdata
  # contains only a subset of factor levels (e.g. scoring a single row).
  sc <- object$scaling
  mf_sel <- model.frame(tt_sel, newdata, na.action = na.pass,
                         xlev = sc$xlev_S)
  mf_out <- model.frame(tt_out, newdata, na.action = na.pass,
                         xlev = sc$xlev_O)
  XS_new <- model.matrix(tt_sel, mf_sel, contrasts.arg = sc$contrasts_S)
  XO_new <- model.matrix(tt_out, mf_out, contrasts.arg = sc$contrasts_O)

  # --- Column-name alignment check ---
  if (!identical(colnames(XS_new), names(object$coefficients$selection)))
    stop("Selection design matrix columns in newdata do not match training: ",
         paste(colnames(XS_new), collapse=","), " vs ",
         paste(names(object$coefficients$selection), collapse=","))
  if (!identical(colnames(XO_new), names(object$coefficients$outcome)))
    stop("Outcome design matrix columns in newdata do not match training: ",
         paste(colnames(XO_new), collapse=","), " vs ",
         paste(names(object$coefficients$outcome), collapse=","))

  # Coefficients in object$coefficients are on the ORIGINAL covariate scale
  # (they were back-transformed during model construction).  Therefore we
  # must apply them to the RAW (unstandardised) newdata design matrix —
  # do NOT re-apply centres/scales.
  k1 <- drop(XO_new %*% object$coefficients$outcome)
  k2 <- drop(XS_new %*% object$coefficients$selection)
  rho_hat <- object$rho
  switch(type,
    conditional = { ps <- pnorm(k2); pj <- pbivnorm(k1, k2, rho_hat)
      ifelse(ps > 1e-8, pmin(pmax(pj/ps, 0), 1), NA_real_) },
    marginal = pnorm(k1), selection = pnorm(k2))
}

#' Log-likelihood of a RidgeBPSS fit
#' @param object a RidgeBPSS fit.
#' @param ... ignored.
#' @return an object of class "logLik".
#' @export
logLik.RidgeBPSS <- function(object, ...) {
  val <- object$logLik; attr(val, "df") <- object$df
  attr(val, "nobs") <- object$nobs; class(val) <- "logLik"; val
}


# =======================================================================
# BOOTSTRAP METHOD
# =======================================================================

#' @rdname bootValidate
#' @param data the original training data frame used to fit `object`.
#' @param mboot number of bootstrap samples (default 200).
#' @param seed integer for set.seed, or NULL for no reseeding.
#' @param verbose logical; print progress dots if TRUE.
#' @param retune logical; if TRUE (default), each refit replays the
#'   original fit_spec; if FALSE, each refit pins to the selected
#'   lambda / lambda.rho from the original fit.
#' @param control optional list of control overrides merged into
#'   `fit_spec$control_input`.
#' @export
bootValidate.RidgeBPSS <- function(object, data, mboot = 200L,
                                     seed = NULL, verbose = FALSE,
                                     retune = TRUE, control = NULL,
                                     ...) {
  if (!isTRUE(object$converged))
    stop("Cannot validate a non-converged fit; check object$converged.")
  if (is.null(object$fit_spec))
    stop("Object lacks $fit_spec; cannot replay fitting procedure. ",
         "This usually means the fit was created by an older version ",
         "of RidgeBPSS; refit and try again.")
  n_orig <- nrow(data)
  if (n_orig != length(object$prep$YS))
    stop("'data' has ", n_orig, " rows but original fit had ",
         length(object$prep$YS), " observations.")

  # Extras via ... are ignored with a warning.  See the double-control
  # bug story in HeckSelect for context -- this guard prevents the
  # same silent-failure class from recurring in RidgeBPSS.
  if (...length() > 0L) {
    unused <- ...names()
    warning("Unused arguments to bootValidate.RidgeBPSS: ",
             paste(sprintf("'%s'", unused), collapse = ", "),
             ".  Refit configuration comes from object$fit_spec and ",
             "the explicit 'control' argument; any other arguments ",
             "are ignored.  If you wanted to override the control ",
             "list, use control = list(...) explicitly.")
  }

  fs <- object$fit_spec
  selection <- object$formula$selection
  outcome   <- object$formula$outcome

  # Build the control list for refits: start from fs$control_input,
  # merge any user-supplied overrides.
  if (!is.null(control) && !is.list(control))
    stop("'control' must be a list or NULL.")
  control_boot <- fs$control_input
  if (!is.null(control)) {
    control_boot[names(control)] <- control
  }

  # Refit closure (single-argument, no ... forwarding -- no
  # possibility of double-argument conflicts).  retune = TRUE replays
  # fs input values; retune = FALSE pins selected (lambda, lambda.rho).
  refit_fn <- function(boot_dat) {
    if (retune) {
      RidgeBPSS(selection, outcome, data = boot_dat,
                 lambda     = fs$lambda_input,
                 lambda.rho = fs$lambda_rho_input,
                 crit       = fs$crit,
                 standardize = fs$standardize,
                 control    = control_boot)
    } else {
      RidgeBPSS(selection, outcome, data = boot_dat,
                 lambda     = object$lambda,
                 lambda.rho = object$lambda.rho,
                 crit       = fs$crit,
                 standardize = fs$standardize,
                 control    = control_boot)
    }
  }

  # Extractor closure: (y, p) on selected observations with finite p.
  extract_yp <- function(fit) {
    sel <- which(fit$prep$YS == 1L)
    if (length(sel) == 0L) return(list(y = numeric(0), p = numeric(0)))
    p <- fit$fitted.values[sel]
    y <- fit$prep$YO[sel]
    ok <- !is.na(p)
    list(y = y[ok], p = p[ok])
  }

  # Evaluation index: aligned once at the start (so apparent_y,
  # apparent_p, and predict_on all use the same observations).
  sel_orig <- which(object$prep$YS == 1L)
  p_orig_at_sel <- object$fitted.values[sel_orig]
  ok_app <- !is.na(p_orig_at_sel)
  eval_orig <- sel_orig[ok_app]

  apparent_y <- object$prep$YO[eval_orig]
  apparent_p <- object$fitted.values[eval_orig]

  predict_on_eval <- function(fit, newdata) {
    p_full <- suppressWarnings(predict(fit, newdata = newdata,
                                          type = "conditional"))
    p_full[eval_orig]
  }

  out <- .bootValidate_engine(
    refit_fn    = refit_fn,
    extract_yp  = extract_yp,
    predict_on  = predict_on_eval,
    data        = data,
    apparent_y  = apparent_y,
    apparent_p  = apparent_p,
    mboot       = mboot,
    seed        = seed,
    verbose     = verbose
  )

  # Attach RidgeBPSS-specific return fields
  out$coef       <- object$coefficients
  out$lambda     <- object$lambda
  out$lambda.rho <- object$lambda.rho
  out
}