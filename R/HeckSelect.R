# -----------------------------------------------------------------------
# HeckSelect: Lasso and Adaptive Lasso for Bivariate Probit Model with Sample Selection
#
# Main user-facing function: HeckSelect()
# Shared bootstrap infrastructure: see R/bpss-shared.R
# Shared utilities: see R/zzz-utils.R
# -----------------------------------------------------------------------

#' @importFrom pbivnorm pbivnorm
#' @importFrom stats glm.fit pnorm dnorm qnorm model.frame model.matrix
#'   terms delete.response na.pass model.response quantile .getXlevels
#'   setNames var sd coef predict as.formula
#' @importFrom graphics par plot lines abline axis legend matplot points
#'   mtext
#' @importFrom grDevices rainbow
#' @importFrom utils tail head
NULL


# =======================================================================
# NORMAL-DEPENDENCE LIKELIHOODS (internal)
# =======================================================================

#' @keywords internal
#' @noRd
loglik_normal <- function(theta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho <- tail(ibetaO, 1) + 1
  betaS <- theta[ibetaS]; betaO <- theta[ibetaO]; rho <- tanh(theta[irho])
  nObs <- length(YS)
  rho_rep <- rep(rho, nObs)
  k1 <- drop(XO %*% betaO); k2 <- drop(XS %*% betaS)
  pphink2 <- pmax(pnorm(-k2), 1e-16)
  pphi2n  <- pmin(pmax(pbivnorm(k2, -k1, -rho_rep), 1e-16), 1 - 1e-16)
  pphi2p  <- pmin(pmax(pbivnorm(k2,  k1,  rho_rep), 1e-16), 1 - 1e-16)
  ll <- (1 - YS) * log(pphink2) +
        YS * (1 - YO) * log(pphi2n) +
        YS * YO * log(pphi2p)
  -sum(ll)
}

#' @keywords internal
#' @noRd
gradlik_normal <- function(theta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho <- tail(ibetaO, 1) + 1
  betaS <- theta[ibetaS]; betaO <- theta[ibetaO]; rho <- tanh(theta[irho])
  r <- sqrt(1 - rho^2)
  k1 <- drop(XO %*% betaO); k2 <- drop(XS %*% betaS)
  phik1 <- dnorm(k1); phik2 <- dnorm(k2)
  omg1 <- (k2 - rho * k1) / r; omg2 <- (-k1 + rho * k2) / r
  pphiomg1 <- pnorm(omg1); pphinomg2 <- pnorm(-omg2)
  pphink2 <- pmax(pnorm(-k2), 1e-16)
  pphi2n <- pmin(pmax(pbivnorm(k2, -k1, -rho), 1e-16), 1 - 1e-16)
  pphi2p <- pmin(pmax(pbivnorm(k2,  k1,  rho), 1e-16), 1 - 1e-16)
  phi2 <- dbivnorm(k2, k1, rho)
  z1n <- (phik1 * pphiomg1) / pphi2n
  z1p <- (phik1 * pphiomg1) / pphi2p
  z2n <- (phik2 * pnorm(omg2)) / pphi2n
  z2p <- (phik2 * pphinomg2) / pphi2p
  vn <- phi2 / pphi2n; vp <- phi2 / pphi2p
  db <- YS * (1 - YO) * (-z1n) + YS * YO * z1p
  dg <- (1 - YS) * (-phik2 / pphink2) + YS * (1 - YO) * z2n + YS * YO * z2p
  drhostar <- YS * (1 - YO) * (-r^2 * vn) + YS * YO * (r^2 * vp)
  score <- rep(NA_real_, irho)
  score[ibetaS] <- drop(crossprod(XS, dg))
  score[ibetaO] <- drop(crossprod(XO, db))
  score[irho]   <- sum(drhostar)
  -score
}

#' @keywords internal
#' @noRd
hesslik_normal <- function(theta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho <- tail(ibetaO, 1) + 1; nParam <- irho
  betaS <- theta[ibetaS]; betaO <- theta[ibetaO]; rho <- tanh(theta[irho])
  r <- sqrt(1 - rho^2); r1 <- 1 / r
  k1 <- drop(XO %*% betaO); k2 <- drop(XS %*% betaS)
  phik1 <- dnorm(k1); phik2 <- dnorm(k2)
  pphink2 <- pmax(pnorm(-k2), 1e-16)
  omg1 <- (k2 - rho * k1) / r; omg2 <- (-k1 + rho * k2) / r
  pphi2n <- pmin(pmax(pbivnorm(k2, -k1, -rho), 1e-16), 1 - 1e-16)
  pphi2p <- pmin(pmax(pbivnorm(k2,  k1,  rho), 1e-16), 1 - 1e-16)
  phi2 <- dbivnorm(k2, k1, rho)
  z1n <- (phik1 * pnorm(omg1)) / pphi2n
  z1p <- (phik1 * pnorm(omg1)) / pphi2p
  z2n <- (phik2 * pnorm(omg2)) / pphi2n
  z2p <- (phik2 * pnorm(-omg2)) / pphi2p
  vn <- phi2 / pphi2n; vp <- phi2 / pphi2p
  dz1nb <- -k1 * z1n - rho * vn + z1n^2
  dz1ng <- vn - z1n * z2n
  dz1nr <- r^2 * vn * (omg2 * r1 + z1n)
  dz1pb <- -k1 * z1p - rho * vp - z1p^2
  dz1pg <- vp - z1p * z2p
  dz1pr <- r^2 * vp * (omg2 * r1 - z1p)
  dz2ng <- -k2 * z2n + rho * vn - z2n^2
  dz2nr <- r^2 * vn * (omg1 * r1 + z2n)
  dz2pg <- -k2 * z2p - rho * vp - z2p^2
  dz2pr <- r^2 * vp * (-omg1 * r1 - z2p)
  w_bb <- YS * (1 - YO) * (-dz1nb) + YS * YO * dz1pb
  ddlb <- crossprod(XO * as.vector(w_bb), XO)
  w_gg0 <- (1 - YS) * ((k2 * phik2 / pphink2) - (phik2 / pphink2)^2)
  w_gg  <- w_gg0 + YS * (1 - YO) * dz2ng + YS * YO * dz2pg
  ddlg <- crossprod(XS * as.vector(w_gg), XS)
  w_bg <- YS * (1 - YO) * (-dz1ng) + YS * YO * dz1pg
  ddlbg <- crossprod(XO * as.vector(w_bg), XS)
  ddlgb <- crossprod(XS * as.vector(w_bg), XO)
  w_aa <- YS * (1 - YO) * (-r^2 * vn * (-rho - omg1 * omg2 + r^2 * vn)) +
          YS * YO       * ( r^2 * vp * (-rho - omg1 * omg2 - r^2 * vp))
  ddlrhostar <- sum(w_aa)
  w_ba <- YS * (1 - YO) * (-dz1nr) + YS * YO * dz1pr
  ddllba <- drop(crossprod(XO, w_ba))
  w_ga <- YS * (1 - YO) * dz2nr + YS * YO * dz2pr
  ddllga <- drop(crossprod(XS, w_ga))
  hessD <- matrix(NA_real_, nParam, nParam)
  hessD[ibetaO, ibetaO] <- ddlb
  hessD[ibetaS, ibetaS] <- ddlg
  hessD[ibetaS, ibetaO] <- ddlgb
  hessD[ibetaO, ibetaS] <- ddlbg
  hessD[irho, irho]     <- ddlrhostar
  hessD[ibetaO, irho]   <- ddllba
  hessD[ibetaS, irho]   <- ddllga
  hessD[irho, ibetaO]   <- ddllba
  hessD[irho, ibetaS]   <- ddllga
  -hessD
}


#' @keywords internal
#' @noRd
.prob_floor <- function(p, eps = 1e-16) pmin(pmax(p, eps), 1 - eps)

# =======================================================================
# AMH COPULA LIKELIHOODS (internal)
# =======================================================================

#' @keywords internal
#' @noRd
.amh_primitives <- function(theta, YS, XS, YO, XO) {
  NXS <- ncol(XS); NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho   <- tail(ibetaO, 1) + 1
  betaS  <- theta[ibetaS]; betaO <- theta[ibetaO]; eta <- theta[irho]
  rho    <- tanh(eta); r2 <- 1 - rho^2
  k1     <- drop(XO %*% betaO); k2 <- drop(XS %*% betaS)
  u      <- pnorm(k2); v <- pnorm(k1)
  ubar   <- 1 - u; vbar <- 1 - v
  # Floors on u, v, ubar, vbar to guard reciprocals (never floor rho).
  u_f    <- pmax(u,    1e-16); v_f    <- pmax(v,    1e-16)
  ubar_f <- pmax(ubar, 1e-16); vbar_f <- pmax(vbar, 1e-16)
  D      <- 1 - rho * ubar * vbar
  D_f    <- pmax(D,    1e-16)        # D > 0 theoretically; floor for safety
  E      <- 1 - rho * ubar
  E_f    <- pmax(E,    1e-16)
  phi1   <- dnorm(k1); phi2 <- dnorm(k2)
  A1     <- rho * vbar / D_f
  A2     <- rho * ubar / D_f
  B      <- ubar * vbar / D_f
  list(
    NXS = NXS, NXO = NXO, ibetaS = ibetaS, ibetaO = ibetaO, irho = irho,
    rho = rho, r2 = r2, k1 = k1, k2 = k2,
    u = u, v = v, ubar = ubar, vbar = vbar,
    u_f = u_f, v_f = v_f, ubar_f = ubar_f, vbar_f = vbar_f,
    D = D, D_f = D_f, E = E, E_f = E_f,
    phi1 = phi1, phi2 = phi2,
    A1 = A1, A2 = A2, B = B
  )
}

#' @keywords internal
#' @noRd
loglik_amh <- function(theta, YS, XS, YO, XO) {
  P <- .amh_primitives(theta, YS, XS, YO, XO)
  i00 <- YS == 0
  i10 <- YS == 1 & YO == 0
  i11 <- YS == 1 & YO == 1
  # Cell log-probabilities, with floors.
  # P(S=0) = ubar;  P(S=1,Y=1) = uv/D;  P(S=1,Y=0) = u*vbar*E/D.
  p00 <- .prob_floor(P$ubar)
  p10 <- .prob_floor(P$u * P$vbar * P$E / P$D)
  p11 <- .prob_floor(P$u * P$v / P$D)
  ll <- numeric(length(YS))
  ll[i00] <- log(p00[i00])
  ll[i10] <- log(p10[i10])
  ll[i11] <- log(p11[i11])
  -sum(ll)
}

#' @keywords internal
#' @noRd
gradlik_amh <- function(theta, YS, XS, YO, XO) {
  P <- .amh_primitives(theta, YS, XS, YO, XO)
  i00 <- YS == 0
  i10 <- YS == 1 & YO == 0
  i11 <- YS == 1 & YO == 1
  nObs <- length(YS)

  # Per-observation contributions to gamma-coefficient score:
  #   S=0:   -phi2 / ubar
  #   S=10:  phi2 * (1/u + rho/E - A1)
  #   S=11:  phi2 * (1/u - A1)
  dg_scalar <- numeric(nObs)
  dg_scalar[i00] <- -P$phi2[i00] / P$ubar_f[i00]
  dg_scalar[i10] <-  P$phi2[i10] * (1 / P$u_f[i10] + P$rho / P$E_f[i10] - P$A1[i10])
  dg_scalar[i11] <-  P$phi2[i11] * (1 / P$u_f[i11] - P$A1[i11])

  # Per-observation contributions to beta-coefficient score:
  #   S=0:   0
  #   S=10:  phi1 * (-1/vbar - A2)
  #   S=11:  phi1 * (1/v - A2)
  db_scalar <- numeric(nObs)
  db_scalar[i10] <- P$phi1[i10] * (-1 / P$vbar_f[i10] - P$A2[i10])
  db_scalar[i11] <- P$phi1[i11] * ( 1 / P$v_f[i11]    - P$A2[i11])

  # Per-observation contributions to eta score:
  #   S=0:   0
  #   S=10:  -r2 * ubar * v / (D * E)
  #   S=11:   r2 * B
  de_scalar <- numeric(nObs)
  de_scalar[i10] <- -P$r2 * P$ubar[i10] * P$v[i10] / (P$D_f[i10] * P$E_f[i10])
  de_scalar[i11] <-  P$r2 * P$B[i11]

  score <- rep(NA_real_, P$irho)
  score[P$ibetaS] <- drop(crossprod(XS, dg_scalar))
  score[P$ibetaO] <- drop(crossprod(XO, db_scalar))
  score[P$irho]   <- sum(de_scalar)
  -score
}

#' @keywords internal
#' @noRd
hesslik_amh <- function(theta, YS, XS, YO, XO) {
  P <- .amh_primitives(theta, YS, XS, YO, XO)
  nObs <- length(YS); npar <- P$irho
  i00 <- YS == 0
  i10 <- YS == 1 & YO == 0
  i11 <- YS == 1 & YO == 1

  # Abbreviations
  u <- P$u_f; v <- P$v_f; ubar <- P$ubar_f; vbar <- P$vbar_f
  D <- P$D_f; E <- P$E_f
  rho <- P$rho; r2 <- P$r2
  phi1 <- P$phi1; phi2 <- P$phi2
  A1 <- P$A1; A2 <- P$A2; B <- P$B
  D2 <- D^2; E2 <- E^2

  # Scalar derivatives of score multipliers per observation type.
  # NOTE: these are contributions to Hess of LOG-likelihood.
  # -------------- S = 0: m_gamma = -phi2/ubar --------------
  # d(m_gamma)/dk2 = -(dphi2/dk2)/ubar - phi2*(-1)*(-(-phi2))/ubar^2
  #                = -(-k2*phi2)/ubar - phi2^2/ubar^2
  #                =  k2*phi2/ubar - (phi2/ubar)^2
  M_gg_00 <-  P$k2[i00] * phi2[i00] / ubar[i00] - (phi2[i00] / ubar[i00])^2
  # All other M-scalars are zero for S=0.

  # -------------- S = 1, Y = 1 --------------
  # m_gamma = phi2 * (1/u - A1)
  # d/dk2: dphi2/dk2 * (1/u - A1) + phi2 * (-phi2/u^2 - dA1/dk2)
  dphi2_dk2_11 <- -P$k2[i11] * phi2[i11]
  dphi1_dk1_11 <- -P$k1[i11] * phi1[i11]
  bracket_g_11 <- (1 / u[i11] - A1[i11])
  M_gg_11 <- dphi2_dk2_11 * bracket_g_11 +
    phi2[i11] * (-phi2[i11] / u[i11]^2 - (-rho^2 * vbar[i11]^2 * phi2[i11] / D2[i11]))
  # d m_gamma / dk1 = phi2 * (0 - dA1/dk1) = phi2 * (rho * phi1 / D^2)
  M_gb_11 <- phi2[i11] * (rho * phi1[i11] / D2[i11])
  # d m_gamma / deta = phi2 * (0 - dA1/deta) = phi2 * (-r2 * vbar / D^2)
  M_ge_11 <- phi2[i11] * (-r2 * vbar[i11] / D2[i11])
  # m_beta = phi1 * (1/v - A2)
  # d m_beta / dk1 = dphi1/dk1*(1/v - A2) + phi1*(-phi1/v^2 - dA2/dk1)
  M_bb_11 <- dphi1_dk1_11 * (1 / v[i11] - A2[i11]) +
    phi1[i11] * (-phi1[i11] / v[i11]^2 - (-rho^2 * ubar[i11]^2 * phi1[i11] / D2[i11]))
  # d m_beta / deta = phi1 * (-dA2/deta) = phi1 * (-r2 * ubar / D^2)
  M_be_11 <- phi1[i11] * (-r2 * ubar[i11] / D2[i11])
  # m_eta = r2 * B
  # d m_eta / deta: r2 = 1 - rho^2, d r2/d eta = -2 rho * r2;
  #                 dB/deta = r2 * ubar^2 * vbar^2 / D^2
  # d m_eta / d eta = d(r2)/d eta * B + r2 * dB/d eta
  #                = -2 rho * r2 * B + r2 * r2 * ubar^2 * vbar^2 / D^2
  M_ee_11 <- -2 * rho * r2 * B[i11] + r2^2 * ubar[i11]^2 * vbar[i11]^2 / D2[i11]

  # -------------- S = 1, Y = 0 --------------
  # m_gamma = phi2 * (1/u + rho/E - A1)
  dphi2_dk2_10 <- -P$k2[i10] * phi2[i10]
  dphi1_dk1_10 <- -P$k1[i10] * phi1[i10]
  bracket_g_10 <- (1 / u[i10] + rho / E[i10] - A1[i10])
  # d/dk2 bracket = -phi2/u^2 + rho*(-dE/dk2)/E^2 - dA1/dk2
  #               = -phi2/u^2 + rho*(-rho*phi2)/E^2 + rho^2*vbar^2*phi2/D^2
  dbg_dk2_10 <- -phi2[i10] / u[i10]^2 - rho^2 * phi2[i10] / E2[i10] +
                rho^2 * vbar[i10]^2 * phi2[i10] / D2[i10]
  M_gg_10 <- dphi2_dk2_10 * bracket_g_10 + phi2[i10] * dbg_dk2_10
  # d/dk1 bracket = 0 + 0 - dA1/dk1 = rho*phi1/D^2
  M_gb_10 <- phi2[i10] * (rho * phi1[i10] / D2[i10])
  # d/deta bracket = 0 + (r2/E - rho*(-r2*ubar)/E^2) - dA1/deta
  #               = r2/E + r2*rho*ubar/E^2 - r2*vbar/D^2
  #               = r2*(E + rho*ubar)/E^2 - r2*vbar/D^2
  #               = r2*(1)/E^2 - r2*vbar/D^2        [since E + rho*ubar = 1]
  dbg_deta_10 <- r2 / E2[i10] - r2 * vbar[i10] / D2[i10]
  M_ge_10 <- phi2[i10] * dbg_deta_10
  # m_beta = phi1 * (-1/vbar - A2)
  # d/dk1: dphi1/dk1*(-1/vbar - A2) + phi1*(-(-phi1)/vbar^2*(-1) - dA2/dk1)
  # Careful: d(-1/vbar)/dk1 = -(-1)*(-(-phi1))/vbar^2 = -phi1/vbar^2
  # Actually: -1/vbar, vbar = 1-v, d(-1/vbar)/dv = -(-1)/vbar^2 * (d vbar/d v) = -(-1)(-1)/vbar^2 = -1/vbar^2
  # d(-1/vbar)/dk1 = d(-1/vbar)/dv * phi1 = -phi1/vbar^2
  M_bb_10 <- dphi1_dk1_10 * (-1 / vbar[i10] - A2[i10]) +
    phi1[i10] * (-phi1[i10] / vbar[i10]^2 - (-rho^2 * ubar[i10]^2 * phi1[i10] / D2[i10]))
  # d m_beta / deta = phi1 * (0 - dA2/deta) = phi1 * (-r2 * ubar / D^2)
  M_be_10 <- phi1[i10] * (-r2 * ubar[i10] / D2[i10])
  # m_eta = -r2 * ubar * v / (D * E)
  # This needs careful differentiation.
  # Let N = -r2 * ubar * v,  M = D * E, so m_eta = N / M.
  # d m_eta / d eta = (dN/deta * M - N * dM/deta) / M^2
  N_10 <- -r2 * ubar[i10] * v[i10]
  M_10 <- D[i10] * E[i10]
  dN_deta_10 <- 2 * rho * r2 * ubar[i10] * v[i10]      # d(-r2)/deta = 2 rho r2
  dM_deta_10 <- -r2 * ubar[i10] * vbar[i10] * E[i10] + D[i10] * (-r2 * ubar[i10])
  M_ee_10 <- (dN_deta_10 * M_10 - N_10 * dM_deta_10) / M_10^2

  # ----- Assemble six blocks -----
  # M_gg, M_bb, M_ee: diagonal-like contributions (per-obs scalars)
  # Work with n-vectors, zero where that obs does not contribute.
  M_gg <- numeric(nObs); M_gg[i00] <- M_gg_00; M_gg[i10] <- M_gg_10; M_gg[i11] <- M_gg_11
  M_gb <- numeric(nObs);                        M_gb[i10] <- M_gb_10; M_gb[i11] <- M_gb_11
  M_ge <- numeric(nObs);                        M_ge[i10] <- M_ge_10; M_ge[i11] <- M_ge_11
  M_bb <- numeric(nObs);                        M_bb[i10] <- M_bb_10; M_bb[i11] <- M_bb_11
  M_be <- numeric(nObs);                        M_be[i10] <- M_be_10; M_be[i11] <- M_be_11
  M_ee <- numeric(nObs);                        M_ee[i10] <- M_ee_10; M_ee[i11] <- M_ee_11

  # These are contributions to Hess of LOG-likelihood.
  # Build 2D blocks with diag-weighting trick.
  Hll_gg <- crossprod(XS * as.vector(M_gg), XS)            # p_S x p_S
  Hll_gb <- crossprod(XS * as.vector(M_gb), XO)            # p_S x p_O
  Hll_ge <- drop(crossprod(XS, M_ge))                       # p_S x 1
  Hll_bb <- crossprod(XO * as.vector(M_bb), XO)
  Hll_be <- drop(crossprod(XO, M_be))
  Hll_ee <- sum(M_ee)

  Hll <- matrix(0, npar, npar)
  Hll[P$ibetaS, P$ibetaS] <- Hll_gg
  Hll[P$ibetaS, P$ibetaO] <- Hll_gb
  Hll[P$ibetaO, P$ibetaS] <- t(Hll_gb)
  Hll[P$ibetaS, P$irho]   <- Hll_ge
  Hll[P$irho,   P$ibetaS] <- Hll_ge
  Hll[P$ibetaO, P$ibetaO] <- Hll_bb
  Hll[P$ibetaO, P$irho]   <- Hll_be
  Hll[P$irho,   P$ibetaO] <- Hll_be
  Hll[P$irho,   P$irho]   <- Hll_ee

  # Negative log-likelihood Hessian
  -Hll
}


# =======================================================================
# DATA PREP + PENALTY STRUCTURE (internal)
# =======================================================================

#' @keywords internal
#' @noRd
prepare_heckselect_data <- function(selection, outcome, data, standardize = TRUE) {

  stopifnot("Selection formula must include an intercept" =
              attr(terms(selection), "intercept") == 1)
  stopifnot("Outcome formula must include an intercept" =
              attr(terms(outcome), "intercept") == 1)

  mf_sel <- model.frame(selection, data, na.action = na.pass)
  mf_out <- model.frame(outcome,   data, na.action = na.pass)

  YS <- coerce_binary01(model.response(mf_sel), "Selection response (S)")
  YO <- coerce_binary01(model.response(mf_out), "Outcome response (Y_obs)")

  if (anyNA(YS)) stop("Selection response (S) cannot contain missing values.")
  if (any(is.na(YO[YS == 1L]))) {
    stop("Outcome contains NA among accepted observations (YS = 1). ",
         "All accepted applicants must have observed Y_obs.")
  }
  YO[YS == 0L] <- 0L

  XS <- model.matrix(selection, mf_sel)
  XO <- model.matrix(outcome,   mf_out)
  if (anyNA(XS)) stop("Selection design matrix contains NA values.")
  if (anyNA(XO)) stop("Outcome design matrix contains NA values.")

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

  if (standardize) {
    idx_S <- if (ncol(XS) > 1L) 2:ncol(XS) else integer(0)
    idx_O <- if (ncol(XO) > 1L) 2:ncol(XO) else integer(0)

    glmnet_standardize <- function(X, cols) {
      Xsub   <- X[, cols, drop = FALSE]
      centre <- colMeans(Xsub)
      Xc     <- sweep(Xsub, 2L, centre, "-")
      scale  <- sqrt(colSums(Xc^2) / nrow(Xc))  # n-denominator
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
    attr(XS, "dimnames") <- list(NULL, colnames(XS))
    attr(XO, "dimnames") <- list(NULL, colnames(XO))
  }

  list(YS = YS, YO = YO, XS = XS, XO = XO,
       n = length(YS), NXS = ncol(XS), NXO = ncol(XO),
       scaling = scaling, formula.sel = selection, formula.out = outcome)
}

#' @keywords internal
#' @noRd
make_heck_penalty_factors <- function(theta_hat, NXS, NXO, penalty,
                                       gamma_alasso = 1, eps_alasso = 1e-6,
                                       unpenalized_idx = integer(0)) {
  npar <- NXS + NXO + 1L
  w_main <- rep(1, npar)
  # Intercepts (positions 1 and NXS+1) and eta (position npar) get weight 0
  # against the main lambda; the eta-penalty is applied via lambda.rho.
  w_main[1L]        <- 0
  w_main[NXS + 1L]  <- 0
  w_main[npar]      <- 0

  w_eta <- 1

  if (penalty == "alasso") {
    th <- theta_hat
    # Adaptive weights on the non-intercept, non-eta coefficients.
    # penalized positions are 2:NXS (selection) and (NXS+2):(NXS+NXO) (outcome).
    pen_idx <- c(if (NXS > 1L) 2:NXS else integer(0),
                 if (NXO > 1L) (NXS + 2L):(NXS + NXO) else integer(0))
    for (j in pen_idx) {
      w_main[j] <- 1 / max(abs(th[j]), eps_alasso)^gamma_alasso
    }
    # Adaptive weight for eta
    w_eta <- 1 / max(abs(th[npar]), eps_alasso)^gamma_alasso
  }

  # User-specified unpenalized coordinates override everything else.
  if (length(unpenalized_idx) > 0L) {
    valid <- unpenalized_idx[unpenalized_idx >= 1L & unpenalized_idx <= npar]
    w_main[valid] <- 0
  }

  list(w_main = w_main, w_eta = w_eta)
}



# =======================================================================
# PILOT FITTING (internal)
# =======================================================================

#' @keywords internal
#' @noRd
fit_unpenalized_heck <- function(prep, model, control = list()) {
  # Resolve pilot preference.  NULL (default) = auto: prefer GJRM, fall back
  # to internal on failure.  Explicit "gjrm" or "internal" disables fallback.
  user_method <- control$pilot_method
  auto <- is.null(user_method)
  pilot_method <- if (auto) {
    if (requireNamespace("GJRM", quietly = TRUE)) "gjrm" else "internal"
  } else {
    match.arg(user_method, c("gjrm", "internal"))
  }

  if (pilot_method == "gjrm") {
    if (!requireNamespace("GJRM", quietly = TRUE)) {
      stop("GJRM package is required for pilot_method = 'gjrm'.  ",
           "Install it with install.packages('GJRM'), or set ",
           "control$pilot_method = 'internal' to use the built-in ",
           "trust-region pilot, or leave control$pilot_method unset ",
           "to let HeckSelect choose automatically.")
    }
    gjrm_fit <- .fit_pilot_gjrm(prep, model, control)
    if (isTRUE(gjrm_fit$converged)) return(gjrm_fit)
    # GJRM failed.  If the user explicitly requested GJRM, return the
    # failed object.  Otherwise, fall back to the internal pilot.
    if (!auto) return(gjrm_fit)
    message("HeckSelect: GJRM pilot did not converge on this dataset; ",
            "falling back to internal trust-region pilot. ",
            "(This is expected on rare-event data or highly correlated ",
            "predictors.  Force GJRM with control$pilot_method = 'gjrm' ",
            "to see the failure, or force internal with ",
            "control$pilot_method = 'internal' to skip GJRM entirely.)")
  }
  .fit_pilot_internal(prep, model, control)
}

#' @keywords internal
#' @noRd
.fit_pilot_gjrm <- function(prep, model, control) {
  # Build data frame with back-named columns for gjrm formulas.
  df <- data.frame(
    .S = prep$YS,
    .Y = prep$YO
  )
  # Attach selection-equation columns (skip intercept).
  if (prep$NXS > 1L) {
    for (j in 2:prep$NXS) {
      df[[paste0("s_", j - 1)]] <- prep$XS[, j]
    }
  }
  if (prep$NXO > 1L) {
    for (j in 2:prep$NXO) {
      df[[paste0("o_", j - 1)]] <- prep$XO[, j]
    }
  }
  s_rhs <- if (prep$NXS > 1L)
    paste(paste0("s_", 1:(prep$NXS - 1)), collapse = " + ") else "1"
  o_rhs <- if (prep$NXO > 1L)
    paste(paste0("o_", 1:(prep$NXO - 1)), collapse = " + ") else "1"
  f_sel <- as.formula(paste(".S ~", s_rhs))
  f_out <- as.formula(paste(".Y ~", o_rhs))

  bivd <- switch(model, "Normal" = "N", "AMH" = "AMH")

  # GJRM can emit benign convergence warnings (e.g. about smoothing
  # parameters) on perfectly fine fits; we must NOT treat those as
  # convergence failure.  Only hard errors kill the fit.  Actual
  # convergence is judged below by inspecting the fit object directly.
  #
  # GJRM API compatibility: the current CRAN-documented public API for
  # gjrm() (as of GJRM 0.2-6.8, July 2025) uses LOWERCASE argument names
  # `copula` and `model`.  Older versions of GJRM (and the internal
  # SemiParBIV() / copulaSampleSel() dispatchers underneath) used
  # capitalised `BivD` and `Model`.  We dispatch explicitly based on
  # what the installed version's gjrm() function actually exposes, by
  # inspecting its formal argument names -- this is robust to error
  # message text, locale differences, and future version drift.
  gjrm_args <- names(formals(GJRM::gjrm))
  use_lowercase <- all(c("copula", "model") %in% gjrm_args)
  use_uppercase <- all(c("BivD", "Model") %in% gjrm_args)

  if (!use_lowercase && !use_uppercase) {
    return(list(converged = FALSE,
                theta_hat = rep(NA_real_, prep$NXS + prep$NXO + 1L),
                H = NULL, logLik = NA_real_, method = "gjrm",
                message = paste("Installed GJRM::gjrm() does not expose",
                                 "expected arguments (neither",
                                 "'copula'/'model' nor 'BivD'/'Model');",
                                 "cannot dispatch safely.")))
  }

  # Build the call with the appropriate argument names.  We construct
  # the call via do.call() so the same argument list can be reused
  # with either naming convention without repeating code.
  gjrm_call_args <- if (use_lowercase) {
    list(formula = list(f_sel, f_out), data = df,
         copula = bivd, margins = c("probit", "probit"),
         model = "BSS")
  } else {
    list(formula = list(f_sel, f_out), data = df,
         BivD = bivd, margins = c("probit", "probit"),
         Model = "BSS")
  }

  fit <- tryCatch(
    suppressWarnings(do.call(GJRM::gjrm, gjrm_call_args)),
    error = function(e) NULL
  )

  # Convergence check: GJRM uses different convergence flags across versions
  # (fit$fit$conv, fit$fit$conv.sp).  The most reliable signal is
  # (a) fit not NULL, (b) logLik finite, (c) coefficient vector of the
  # correct length AND all finite.  The finite-coef check is essential:
  # in difficult data, GJRM can return a fit object with NaN in some
  # coordinates, which would silently corrupt the downstream LSA
  # surrogate centre.
  coefs_try <- tryCatch(coef(fit, part = "full"), error = function(e) NULL)
  gjrm_converged <- !is.null(fit) &&
                    is.finite(tryCatch(fit$logLik, error = function(e) NA_real_)) &&
                    !is.null(coefs_try) &&
                    length(coefs_try) == (prep$NXS + prep$NXO + 1L) &&
                    all(is.finite(coefs_try))
  if (!gjrm_converged) {
    reason <- if (is.null(fit)) "GJRM returned NULL (no fit object)."
              else if (!is.finite(tryCatch(fit$logLik, error = function(e) NA_real_)))
                "GJRM returned non-finite logLik."
              else if (is.null(coefs_try))
                "GJRM coef(part='full') extraction failed."
              else if (length(coefs_try) != (prep$NXS + prep$NXO + 1L))
                sprintf("GJRM returned %d coefficients, expected %d.",
                        length(coefs_try), prep$NXS + prep$NXO + 1L)
              else "GJRM returned non-finite coefficients."
    return(list(converged = FALSE,
                theta_hat = rep(NA_real_, prep$NXS + prep$NXO + 1L),
                H = NULL, logLik = NA_real_, method = "gjrm",
                message = reason))
  }

  coefs  <- coefs_try
  npar   <- prep$NXS + prep$NXO + 1L
  # GJRM returns coefficients in the order: selection eqn, outcome eqn, then
  # the dependence parameter.  Our theta ordering matches.
  theta_hat <- as.numeric(coefs)
  names(theta_hat) <- c(paste0("S:", colnames(prep$XS)),
                        paste0("O:", colnames(prep$XO)),
                        "eta")

  # Evaluate analytic value, gradient, and Hessian at the GJRM point.
  # ALL THREE must be finite AND the gradient must be small (stationary)
  # for us to accept GJRM's point as the LSA centre.  A finite logLik
  # and Hessian are necessary but not sufficient: GJRM may have stopped
  # at a point that is stationary under ITS penalized smoothing-spline
  # objective but NOT stationary under our analytic likelihood.  Using
  # such a non-stationary point as the LSA centre would make the
  # quadratic surrogate incorrect.
  loglik_fun <- if (model == "Normal") loglik_normal else loglik_amh
  grad_fun   <- if (model == "Normal") gradlik_normal else gradlik_amh
  hess_fun   <- if (model == "Normal") hesslik_normal else hesslik_amh
  H <- tryCatch(hess_fun(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO),
                 error = function(e) NULL)
  logLik_val <- tryCatch(
    -loglik_fun(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO),
    error = function(e) NA_real_)
  g_at_pilot <- tryCatch(
    grad_fun(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO),
    error = function(e) NULL)

  # Finiteness checks
  if (is.null(H) || any(!is.finite(H)) || !is.finite(logLik_val) ||
      is.null(g_at_pilot) || any(!is.finite(g_at_pilot))) {
    return(list(converged = FALSE,
                theta_hat = rep(NA_real_, npar),
                H = NULL, logLik = NA_real_, method = "gjrm",
                message = paste("GJRM converged but our analytic value,",
                                 "gradient, or Hessian at the returned",
                                 "point is non-finite; rejecting for stability.")))
  }

  # Stationarity check: require sup-norm of gradient to be small in both
  # absolute and relative terms.  For BPSS the likelihood scales as O(n),
  # so we use a scale-aware threshold that matches the post-fit sanity
  # check in .fit_pilot_internal().  A gradient of 1e-2 at n = 1000 is
  # well within typical sampling noise and does not necessarily indicate
  # non-stationarity; we keep the bar intentionally modest to avoid
  # rejecting GJRM pilots that are "good enough" for LSA centring, while
  # still catching genuinely non-stationary returns.
  grad_sup <- max(abs(g_at_pilot))
  grad_tol <- max(1e-2, 1e-4 * max(abs(logLik_val), 1))
  if (grad_sup > grad_tol) {
    return(list(converged = FALSE,
                theta_hat = rep(NA_real_, npar),
                H = NULL, logLik = NA_real_, method = "gjrm",
                message = sprintf(paste("GJRM converged but our analytic",
                                          "gradient at the returned point",
                                          "has sup-norm %.3g > tol %.3g;",
                                          "rejecting because the LSA centre",
                                          "would not be stationary."),
                                    grad_sup, grad_tol)))
  }

  list(converged = TRUE, theta_hat = theta_hat, H = H, logLik = logLik_val,
       method = "gjrm")
}

#' @keywords internal
#' @noRd
.fit_pilot_internal <- function(prep, model, control) {
  if (!requireNamespace("trust", quietly = TRUE)) {
    stop("trust package is required for pilot_method = 'internal'. ",
         "Install with install.packages('trust').")
  }
  loglik_fun <- if (model == "Normal") loglik_normal else loglik_amh
  grad_fun   <- if (model == "Normal") gradlik_normal else gradlik_amh
  hess_fun   <- if (model == "Normal") hesslik_normal else hesslik_amh

  # Starting values: probit on YS for gammaS, probit on YO (conditional on YS=1)
  # for betaO, eta = 0.
  gamma0 <- tryCatch(
    suppressWarnings(glm.fit(prep$XS, prep$YS,
                              family = binomial("probit")))$coefficients,
    error = function(e) rep(0, prep$NXS)
  )
  gamma0[!is.finite(gamma0)] <- 0
  accepted <- prep$YS == 1
  beta0 <- tryCatch(
    suppressWarnings(glm.fit(prep$XO[accepted, , drop = FALSE],
                              prep$YO[accepted],
                              family = binomial("probit")))$coefficients,
    error = function(e) rep(0, prep$NXO)
  )
  beta0[!is.finite(beta0)] <- 0
  theta0 <- c(gamma0, beta0, 0)

  objfun <- function(theta) {
    # Robust fallback pattern (matches RidgeBPSS): if ANY of value,
    # gradient, or Hessian is non-finite (either via an error or via
    # NaN/Inf entries), return a SINGLE coherent fallback bundle:
    #   value    = 1e20   (trust-region step will shrink)
    #   gradient = 0      (no preferred direction)
    #   hessian  = I      (identity, well-conditioned, no ill-guidance)
    #
    # The earlier element-wise patching (zero non-finite gradient
    # components, zero non-finite Hessian entries) was dangerous:
    # trust() would receive a mixed bundle where a finite value and
    # partially-valid derivatives could guide the optimiser toward a
    # corrupted step.  The full fallback is strictly safer.
    n_th <- length(theta)
    bad_bundle <- list(value = 1e20,
                       gradient = rep(0, n_th),
                       hessian  = diag(n_th))

    v <- tryCatch(loglik_fun(theta, prep$YS, prep$XS, prep$YO, prep$XO),
                  error = function(e) NA_real_)
    if (!is.finite(v)) return(bad_bundle)

    g <- tryCatch(grad_fun(theta, prep$YS, prep$XS, prep$YO, prep$XO),
                  error = function(e) rep(NA_real_, n_th))
    if (any(!is.finite(g))) return(bad_bundle)

    H <- tryCatch(hess_fun(theta, prep$YS, prep$XS, prep$YO, prep$XO),
                  error = function(e) matrix(NA_real_, n_th, n_th))
    if (any(!is.finite(H))) return(bad_bundle)

    list(value = v, gradient = g, hessian = H)
  }

  fit <- tryCatch(
    trust::trust(objfun, theta0, rinit = 1, rmax = 100, iterlim = 500,
                 blather = FALSE),
    error = function(e) list(converged = FALSE, argument = theta0, value = Inf)
  )

  # Accept-if-effectively-stationary fallback:
  # trust::trust can report converged = FALSE when its trust-region radius
  # collapses because the proposed steps all fail to improve a function
  # that is already numerically at its minimum (this happens for the
  # intercept-only model, where the likelihood surface is essentially
  # flat near the MLE on the eta axis).  We handle this with a fallback:
  # accept the returned point if the true gradient there is small in
  # either absolute or relative terms.
  #
  # Thresholds:
  #   abs  < 1e-4 : gradient sup-norm below 1e-4.  This is one order of
  #                  magnitude tighter than the old 1e-3 threshold; the
  #                  looser value was flagged as potentially too permissive
  #                  in review.  The intercept-only test case (Test 8)
  #                  returns a gradient ~1e-7 at the MLE, well within 1e-4.
  #   rel  < 1e-6 * max(|f|, 1) : scale-aware gradient check.  For BPSS
  #                  the NLL is O(n), so for n=1000 a relative threshold of
  #                  1e-6 corresponds to ~1e-3 absolute.  Tightened from
  #                  1e-5 in response to review.
  #
  # These thresholds are intentionally NOT taken all the way to the
  # trust-region default (1e-8) because in the intercept-only case trust()
  # flags non-convergence even at machine-epsilon gradients.  The fallback
  # is here specifically to rescue that case; the post-fit sanity check
  # below still validates via a true likelihood/gradient/Hessian
  # re-evaluation.
  if (!isTRUE(fit$converged)) {
    arg_try <- if (!is.null(fit$argument)) fit$argument else theta0
    vcheck <- tryCatch(loglik_fun(arg_try, prep$YS, prep$XS, prep$YO, prep$XO),
                        error = function(e) NA_real_)
    gcheck <- tryCatch(grad_fun(arg_try, prep$YS, prep$XS, prep$YO, prep$XO),
                        error = function(e) rep(NA_real_, length(theta0)))
    if (is.finite(vcheck) && all(is.finite(gcheck)) &&
        (max(abs(gcheck)) < 1e-4 ||
         max(abs(gcheck)) < 1e-6 * max(abs(vcheck), 1))) {
      fit <- list(converged = TRUE, argument = arg_try,
                  iterations = if (!is.null(fit$iterations)) fit$iterations else 0L,
                  value = vcheck)
    }
  }

  if (!isTRUE(fit$converged)) {
    return(list(converged = FALSE,
                theta_hat = rep(NA_real_, length(theta0)),
                H = NULL, logLik = NA_real_, method = "internal",
                message = "Internal trust-region pilot did not converge."))
  }

  # Post-fit sanity check.  During the trust iterations the objfun
  # returns sentinel values (1e20, zero gradient, identity Hessian) at
  # points where the likelihood evaluation fails.  If trust happened to
  # terminate at such a point, we must NOT accept the solution.  So we
  # re-evaluate the TRUE objective, gradient, and Hessian at the
  # returned argument and reject if any is non-finite, or if the
  # objective is still sitting near the sentinel.
  theta_hat <- fit$argument
  v_true <- tryCatch(loglik_fun(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO),
                      error = function(e) NA_real_)
  g_true <- tryCatch(grad_fun(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO),
                      error = function(e) rep(NA_real_, length(theta_hat)))
  H_true <- tryCatch(hess_fun(theta_hat, prep$YS, prep$XS, prep$YO, prep$XO),
                      error = function(e) matrix(NA_real_, length(theta_hat),
                                                    length(theta_hat)))
  if (!is.finite(v_true) || v_true >= 1e19 ||
      any(!is.finite(g_true)) || any(!is.finite(H_true))) {
    return(list(converged = FALSE,
                theta_hat = rep(NA_real_, length(theta0)),
                H = NULL, logLik = NA_real_, method = "internal",
                message = paste("Internal pilot converged nominally but the",
                                 "likelihood/gradient/Hessian at the returned",
                                 "point is not finite.  This usually indicates",
                                 "the optimiser terminated in a region where",
                                 "rho is too close to +/-1 or the bivariate",
                                 "CDF evaluation is unstable.")))
  }
  logLik_val <- -v_true
  H <- H_true
  names(theta_hat) <- c(paste0("S:", colnames(prep$XS)),
                        paste0("O:", colnames(prep$XO)),
                        "eta")
  list(converged = TRUE, theta_hat = theta_hat, H = H, logLik = logLik_val,
       method = "internal")
}



# =======================================================================
# LAMBDA GRID AND HESSIAN STABILISATION (internal)
# =======================================================================

#' @keywords internal
#' @noRd
stabilize_hessian <- function(H, eig_floor = 1e-8) {
  # Step 1: symmetrise.
  H <- 0.5 * (H + t(H))
  # Step 2: check if floor is even needed (common case: H is already PD).
  # Use a cheap lower bound via Gershgorin first; fall back to eigen only
  # when necessary.
  diag_min <- min(diag(H))
  if (diag_min > eig_floor) {
    # Try Cholesky; if it succeeds, H is already PD and we can skip eigen.
    ch <- tryCatch(chol(H), error = function(e) NULL)
    if (!is.null(ch)) return(H)
  }
  # Full eigen decomposition with floor.
  e <- eigen(H, symmetric = TRUE)
  # Floor negative / tiny eigenvalues.
  vals <- pmax(e$values, eig_floor)
  V <- e$vectors
  V %*% diag(vals, length(vals)) %*% t(V)
}

#' @keywords internal
#' @noRd
compute_lambda_max <- function(theta_hat, H, w_main, NXS, NXO,
                                 extra_zero_idx = integer(0)) {
  # Compute the KKT-boundary value of the main-penalty lambda.  By
  # default (extra_zero_idx = empty), this computes lambda_max at the
  # reference assumption that eta and all other w_main == 0 coords are
  # FREE to adjust in the trial point.  When extra_zero_idx specifies
  # additional coordinates to force to zero at the trial (e.g. eta,
  # when lambda.rho is large enough to zero it), the trial point has
  # theta = 0 on (penalized union extra_zero_idx), and only the
  # remaining coordinates adjust via the linear system.
  #
  # Using extra_zero_idx = length(theta_hat) (the eta slot) gives
  # lambda_max under the assumption that eta is zeroed by lambda.rho.
  # This value is used in the 2D dispatch as the RHS-endpoint bound
  # that complements the LHS-endpoint bound computed with eta free.
  #
  # CAVEAT: The 2D dispatch combines both endpoint values via
  #   lambda_max_envelope = max( lambda_max(eta free), lambda_max(eta = 0) )
  # to obtain a CONSERVATIVE upper-envelope for the 2D grid.  This is
  # not a fully joint KKT derivation -- it does not solve the coupled
  # system that would yield the exact KKT boundary at every intermediate
  # lambda.rho slice.  In practice the envelope appears to cover the
  # relevant KKT boundaries encountered during tuning, but users should
  # be aware it is a practical compromise rather than a proof.
  npar <- length(theta_hat)
  # penalized coords P have w_main > 0.
  P <- which(w_main > 0)
  if (length(P) == 0L) return(0)

  # Zeroed set at trial: penalized coords plus any user-supplied extra
  # zeros (usually eta).  extra_zero_idx must not overlap P.
  extra_zero_idx <- intersect(extra_zero_idx, seq_len(npar))
  zeroed <- union(P, extra_zero_idx)
  # Adjusting set U: everything not zeroed.
  U <- setdiff(seq_len(npar), zeroed)

  if (length(U) == 0L) {
    # No coords adjust; trial is exactly zero everywhere except we have
    # to use theta_trial - theta_hat = -theta_hat.
    g <- drop(H %*% (-theta_hat))
    return(max(abs(g[P]) / w_main[P]))
  }

  # Solve H_UU * diff_U = H_{U, zeroed} * theta_hat_{zeroed} so that
  # the stationarity condition on U holds at the trial point
  #   theta_trial = (theta_hat_U + diff_U, 0_{zeroed}).
  H_UU <- H[U, U, drop = FALSE]
  H_UZ <- H[U, zeroed, drop = FALSE]
  theta_Z <- theta_hat[zeroed]
  diff_U <- tryCatch(solve(H_UU, H_UZ %*% theta_Z), error = function(e) NULL)
  if (is.null(diff_U)) {
    diff_U <- solve(H_UU + diag(1e-6, length(U)), H_UZ %*% theta_Z)
  }
  diff_U <- as.numeric(diff_U)

  # Build full diff = theta_trial - theta_hat
  diff <- -theta_hat
  diff[U] <- diff_U

  # Gradient g_j = [H diff]_j for penalized j
  g <- drop(H %*% diff)
  max(abs(g[P]) / w_main[P])
}

#' @keywords internal
#' @noRd
compute_lambda_rho_max <- function(theta_hat, H, w_eta,
                                     extra_zero_idx = integer(0)) {
  # Compute the KKT-boundary value of lambda.rho.  By default
  # (extra_zero_idx empty) the trial point forces ONLY eta to zero and
  # allows every other coordinate to adjust via the linear system; this
  # corresponds to the LHS endpoint of the 2D grid (main lambda = 0).
  # For the RHS endpoint (main lambda = lambda_max, all penalized
  # regression coefs zero), pass the penalized regression-coord indices
  # in extra_zero_idx so they are also forced to zero at the trial.
  npar <- length(theta_hat)
  eta_idx <- npar
  extra_zero_idx <- setdiff(intersect(extra_zero_idx, seq_len(npar)), eta_idx)
  zeroed <- union(eta_idx, extra_zero_idx)
  U <- setdiff(seq_len(npar), zeroed)

  if (length(U) == 0L) {
    # Trial is all zero; diff = -theta_hat.
    g_eta <- drop(H[eta_idx, , drop = FALSE] %*% (-theta_hat))
    return(as.numeric(abs(g_eta)) / max(w_eta, .Machine$double.eps))
  }

  # Solve H_UU * diff_U = H_{U, zeroed} * theta_hat_{zeroed} so that
  # stationarity holds on U at the trial point
  #   theta_trial = (theta_hat_U + diff_U, 0_{zeroed}).
  H_UU <- H[U, U, drop = FALSE]
  H_UZ <- H[U, zeroed, drop = FALSE]
  theta_Z <- theta_hat[zeroed]
  diff_U <- tryCatch(solve(H_UU, H_UZ %*% theta_Z), error = function(e) NULL)
  if (is.null(diff_U)) {
    diff_U <- solve(H_UU + diag(1e-6, length(U)), H_UZ %*% theta_Z)
  }
  diff_U <- as.numeric(diff_U)

  # Build full diff = theta_trial - theta_hat
  diff <- -theta_hat
  diff[U] <- diff_U

  # Gradient on eta at the trial point
  g_eta <- drop(H[eta_idx, , drop = FALSE] %*% diff)
  as.numeric(abs(g_eta)) / max(w_eta, .Machine$double.eps)
}

#' @keywords internal
#' @noRd
build_lambda_grid <- function(lambda_max, n_lambda = 100, lambda_ratio = 1e-4) {
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    # Degenerate: return a trivial grid that will produce the unpenalized fit
    # in the single-point loop.
    return(c(0))
  }
  lambda_min <- lambda_ratio * lambda_max
  exp(seq(log(lambda_max), log(lambda_min), length.out = n_lambda))
}

# =======================================================================
# COORDINATE DESCENT AND PATH (internal)
# =======================================================================

#' @keywords internal
#' @noRd
cd_lsa <- function(theta_hat, H, lambda, lambda.rho, penalty_factors,
                   theta_start = NULL, control = list()) {
  max_iter <- if (!is.null(control$max_iter)) control$max_iter else 1000L
  tol      <- if (!is.null(control$tol))      control$tol      else 1e-8

  npar  <- length(theta_hat)
  NXS_plus_NXO <- npar - 1L
  irho  <- npar

  w_main <- penalty_factors$w_main
  w_eta  <- penalty_factors$w_eta

  # Per-coord L1 penalty magnitudes (the soft-threshold amount, PRE-divide).
  pen <- lambda * w_main
  # eta slot gets its own penalty:  lambda.rho * w_eta
  pen[irho] <- lambda.rho * w_eta

  H_diag <- diag(H)
  # Any non-positive H_diag entries would be pathological — warn and floor.
  if (any(H_diag <= 0)) {
    H_diag[H_diag <= 1e-12] <- 1e-12
  }

  theta <- if (is.null(theta_start)) theta_hat else theta_start
  # Gradient of quadratic part at starting point
  g <- drop(H %*% (theta - theta_hat))

  converged <- FALSE
  for (it in seq_len(max_iter)) {
    max_change <- 0
    for (j in seq_len(npar)) {
      theta_j_old <- theta[j]
      # z_j = theta_j - g_j / H_jj  is the unconstrained update.
      z_j <- theta_j_old - g[j] / H_diag[j]
      if (pen[j] == 0) {
        theta_j_new <- z_j
      } else {
        theta_j_new <- soft_threshold(z_j, pen[j] / H_diag[j])
      }
      d <- theta_j_new - theta_j_old
      if (d != 0) {
        # Update gradient: g += H[, j] * d
        g <- g + H[, j] * d
        theta[j] <- theta_j_new
        if (abs(d) > max_change) max_change <- abs(d)
      }
    }
    if (max_change < tol) { converged <- TRUE; break }
  }

  list(theta = theta, iterations = it, converged = converged, max_change = max_change)
}

#' @keywords internal
#' @noRd
fit_lasso_path <- function(theta_hat, H, lambda_grid, lambda.rho, penalty_factors,
                           control = list()) {
  n_lambda <- length(lambda_grid)
  npar <- length(theta_hat)
  theta_path <- matrix(NA_real_, nrow = npar, ncol = n_lambda)
  rownames(theta_path) <- names(theta_hat)
  iters     <- integer(n_lambda)
  converged <- logical(n_lambda)

  theta_current <- theta_hat   # warm start at the pilot MLE

  for (k in seq_len(n_lambda)) {
    res <- cd_lsa(theta_hat, H, lambda_grid[k], lambda.rho, penalty_factors,
                  theta_start = theta_current, control = control)
    theta_path[, k] <- res$theta
    iters[k]        <- res$iterations
    converged[k]    <- res$converged
    # Warm-start next lambda with current solution
    theta_current <- res$theta
  }
  list(theta_path = theta_path, iterations = iters, converged = converged,
       lambda_grid = lambda_grid)
}

#' @keywords internal
#' @noRd
compute_ic <- function(theta_path, prep, model, pen_factors, lambda_rho,
                        crit = "bic", tol_zero = 1e-8,
                        converged = NULL) {
  loglik_fun <- if (model == "Normal") loglik_normal else loglik_amh
  n_lambda <- ncol(theta_path)
  npar     <- nrow(theta_path)

  # If the caller supplied a convergence mask, path points that did NOT
  # converge during coordinate descent are immediately excluded: their
  # IC is set to Inf so they cannot win the argmin, and nll/df are NA.
  # This prevents a nonconverged partial iterate from being selected by
  # the tuning criterion.  Passing converged = NULL (default) reproduces
  # the pre-fix behaviour (all finite-likelihood points are eligible).
  if (is.null(converged)) converged <- rep(TRUE, n_lambda)
  if (length(converged) != n_lambda) {
    stop("converged mask length (", length(converged),
         ") does not match n_lambda (", n_lambda, ").")
  }

  # Classify coordinates.
  # w_main == 0 marks coordinates that are NOT penalized by the main
  # lambda.  This set includes intercepts (both equations) AND any
  # user-specified always-included covariates supplied via
  # unpenalized.S / unpenalized.O.  Eta (position npar) is handled
  # separately because its penalisation is governed by lambda.rho and
  # w_eta, not by w_main.
  unpen_main <- which(pen_factors$w_main == 0)
  # Exclude eta from unpen_main; it's position npar.
  unpen_main <- setdiff(unpen_main, npar)
  # Coords penalized by the main lambda (positive w_main):
  pen_main <- which(pen_factors$w_main > 0)
  # Eta is unpenalized IFF lambda_rho == 0 (regardless of alasso weight).
  eta_unpen <- (lambda_rho == 0)

  n    <- prep$n
  nll  <- numeric(n_lambda)
  df   <- numeric(n_lambda)
  ic   <- numeric(n_lambda)
  finite <- logical(n_lambda)

  for (k in seq_len(n_lambda)) {
    # Exclude path points that did not converge during CD.
    if (!isTRUE(converged[k])) {
      nll[k] <- NA_real_; df[k] <- NA_real_; ic[k] <- Inf
      finite[k] <- FALSE
      next
    }
    th <- theta_path[, k]
    v <- tryCatch(loglik_fun(th, prep$YS, prep$XS, prep$YO, prep$XO),
                  error = function(e) NA_real_)
    if (!is.finite(v) || v > 1e18) {
      nll[k] <- NA_real_; df[k] <- NA_real_; ic[k] <- Inf
      finite[k] <- FALSE
      next
    }
    finite[k] <- TRUE
    nll[k] <- v
    # Always-unpenalized (intercepts)
    d <- length(unpen_main)
    # penalized main coords: count non-zero
    d <- d + sum(abs(th[pen_main]) > tol_zero)
    # Eta
    if (eta_unpen) {
      d <- d + 1
    } else {
      if (abs(th[npar]) > tol_zero) d <- d + 1
    }
    df[k] <- d
    ic[k] <- switch(crit,
                    bic = 2 * v + log(n) * d,
                    aic = 2 * v + 2 * d,
                    # GCV: 2*NLL / (n*(1-df/n)^2) -- likelihood analog of the
                    # standard RSS/(n*(1-df/n)^2) form.  The factor of 2
                    # matches the convention used in RidgeBPSS so that GCV
                    # values are directly comparable across the two packages.
                    gcv = 2 * v / (n * (1 - d / n)^2))
  }
  list(nll = nll, df = df, ic = ic, finite = finite,
       best = if (any(finite)) which.min(ic) else NA_integer_)
}



# =======================================================================
# MAIN USER-FACING FUNCTION
# =======================================================================

#' (A)LASSO-Penalized Bivariate Probit Model with Sample Selection 
#'
#' Fits a penalized bivariate probit with sample selection model using
#' lasso (Tibshirani 1996) or adaptive lasso (Zou 2006) penalties.
#' Supports Normal-error and AMH copula dependence structures and
#' automatic tuning of the penalty parameters via BIC, AIC, or GCV.
#'
#' @param selection a formula for the selection equation, including an
#'   intercept.
#' @param outcome a formula for the outcome equation, including an
#'   intercept.
#' @param data a data frame containing all variables referenced in the
#'   two formulas, including the binary selection and outcome responses.
#' @param lambda numeric scalar, vector, or NULL.  If NULL (default),
#'   an automatic log-spaced grid is built and tuned by `crit`.  If a
#'   scalar, that single value is used.  If a vector, those values are
#'   used as the grid.
#' @param lambda.rho numeric scalar or NULL.  If 0 (default), the
#'   dependence parameter eta is unpenalized (Model 3a).  If positive,
#'   eta receives an L1 penalty with this weight (Model 3b).  If NULL,
#'   a second log-spaced grid is built and tuned jointly with lambda.
#' @param penalty one of \code{"lasso"} or \code{"alasso"}.
#' @param Model one of \code{"Normal"} or \code{"AMH"}.
#' @param standardize logical; if TRUE (default), design matrix columns
#'   are centred and scaled internally (glmnet convention) and
#'   coefficients are returned on the raw scale.
#' @param crit information criterion for tuning: "bic", "aic", or "gcv".
#' @param gamma_alasso positive scalar, adaptive-lasso exponent gamma.
#'   Default 1.  Used only when penalty = "alasso".
#' @param unpenalized.S character vector of DESIGN-MATRIX column names
#'   (not formula terms) in the selection equation to leave unpenalized.
#' @param unpenalized.O character vector of column names in the outcome
#'   equation to leave unpenalized.
#' @param control an optional list of control arguments:
#'   \describe{
#'     \item{n_lambda}{grid size when lambda = NULL (default 100)}
#'     \item{lambda_ratio}{ratio of lambda_min to lambda_max (default 1e-4)}
#'     \item{n_lambda_rho}{grid size when lambda.rho = NULL (default 20)}
#'     \item{tol}{coordinate descent convergence tolerance (default 1e-8)}
#'     \item{max_iter}{maximum CD sweeps (default 1000)}
#'     \item{pilot_method}{NULL (GJRM-first), "internal", or "gjrm"}
#'   }
#' @param ... currently ignored.
#'
#' @return an object of class "HeckSelect" with components including
#'   `coefficients`, `std_coefficients`, `rho`, `lambda`, `lambda.rho`,
#'   `logLik`, `ic`, `df`, `converged`, `fitted.values`, `path`, and
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
#' fit <- HeckSelect(S ~ x1 + x2 + Z, Y_obs ~ x1 + x2 + V, dat,
#'                    penalty = "lasso", Model = "Normal")
#' print(fit)
#' summary(fit)
#' }
#'
#' @seealso [RidgeBPSS()] for the ridge-penalized variant,
#'   [bootValidate()] for bootstrap optimism correction.
#' @export
HeckSelect <- function(selection, outcome,
                        data = sys.frame(sys.parent()),
                        lambda = NULL, lambda.rho = 0,
                        penalty = c("lasso", "alasso"),
                        Model   = c("Normal", "AMH"),
                        standardize = TRUE,
                        crit = c("bic", "aic", "gcv"),
                        gamma_alasso = 1,
                        unpenalized.S = character(),
                        unpenalized.O = character(),
                        control = list(), ...) {
  cl <- match.call()

  # Capture the ORIGINAL user-supplied fitting specification BEFORE any
  # argument defaulting.  This is used by bootValidate() to replay the
  # original fitting procedure on each bootstrap sample -- the bootstrap
  # optimism-correction logic is only valid if each refit follows the
  # SAME procedure as the original (e.g. auto-tune 2D grid if original
  # did, pin fixed values if original did).  Storing the INPUT not the
  # selected values is essential here.
  fit_spec <- list(
    lambda_input     = lambda,          # NULL, scalar, or vector
    lambda_rho_input = lambda.rho,      # 0 (default), scalar, or NULL
    control_input    = control,         # user's control list (un-defaulted)
    penalty          = match.arg(penalty),
    Model            = match.arg(Model),
    crit             = match.arg(crit),
    standardize      = standardize,
    gamma_alasso     = gamma_alasso,
    unpenalized.S    = unpenalized.S,
    unpenalized.O    = unpenalized.O
  )

  penalty <- fit_spec$penalty
  Model   <- fit_spec$Model
  crit    <- fit_spec$crit

  validate_penalty_arg(lambda,     "lambda")
  validate_penalty_arg(lambda.rho, "lambda.rho")

  # Control defaults
  con <- list(
    n_lambda      = 100L,
    lambda_ratio  = 1e-4,
    n_lambda_rho  = 20L,           # coarser grid for the 2D case
    lambda_rho_max_factor = 1,     # scale factor for the rho-axis in 2D
    tol           = 1e-8,
    max_iter      = 1000L,
    tol_zero      = 1e-8,
    eps_alasso    = 1e-6,
    pilot_method  = NULL
  )
  con[names(control)] <- control

  # --- Prepare data ---
  prep <- prepare_heckselect_data(selection, outcome, data, standardize)
  npar <- prep$NXS + prep$NXO + 1L

  # --- Translate user-specified unpenalized variable names to positions ---
  # unpenalized.S is a character vector of variable names in the selection
  # design matrix (e.g. "x1", "age").  Those variables will receive w_main
  # = 0 regardless of LASSO/ALASSO weighting, so they are always retained
  # in the model.  Same logic for unpenalized.O.  Intercepts are already
  # unpenalized; specifying "(Intercept)" is a no-op but not an error.
  unpen_idx <- integer(0)
  if (length(unpenalized.S) > 0L) {
    match_idx <- match(unpenalized.S, colnames(prep$XS))
    if (any(is.na(match_idx))) {
      stop("unpenalized.S contains variable names not in the selection ",
           "design matrix: ",
           paste(unpenalized.S[is.na(match_idx)], collapse = ", "),
           ".  Available names: ",
           paste(colnames(prep$XS), collapse = ", "))
    }
    unpen_idx <- c(unpen_idx, match_idx)
  }
  if (length(unpenalized.O) > 0L) {
    match_idx <- match(unpenalized.O, colnames(prep$XO))
    if (any(is.na(match_idx))) {
      stop("unpenalized.O contains variable names not in the outcome ",
           "design matrix: ",
           paste(unpenalized.O[is.na(match_idx)], collapse = ", "),
           ".  Available names: ",
           paste(colnames(prep$XO), collapse = ", "))
    }
    # Outcome positions in the concatenated theta vector start at NXS+1.
    unpen_idx <- c(unpen_idx, match_idx + prep$NXS)
  }
  unpen_idx <- unique(unpen_idx)

  # Helper to build NA-populated failed-fit object.
  build_failed_object <- function(reason, pilot = NULL, lambda_used = NA,
                                   lambda_rho_used = NA,
                                   path_grid = NULL, start_values = NULL) {
    structure(list(
      call = cl, converged = FALSE,
      coefficients = list(selection = NA, outcome = NA),
      std_coefficients = list(selection = NA, outcome = NA),
      eta = NA_real_, rho = NA_real_,
      lambda = lambda_used, lambda.rho = lambda_rho_used,
      logLik = NA_real_, ic = Inf, df = NA_real_, crit = crit,
      converged_path = NULL, iterations = NA_integer_,
      fitted.values = NA, linear.predictors = list(outcome = NA, selection = NA),
      path = path_grid,
      formula = list(selection = selection, outcome = outcome),
      scaling = prep$scaling, nobs = prep$n, npar = npar,
      penalty = penalty, Model = Model, standardize = isTRUE(standardize),
      pilot = pilot, message = reason, start_values = start_values,
      fit_spec = fit_spec
    ), class = "HeckSelect")
  }

  # --- Pilot fit ---
  pilot <- fit_unpenalized_heck(prep, Model, control = con)
  if (!isTRUE(pilot$converged)) {
    warning("Unpenalized pilot fit did not converge; returning failed-fit object.")
    return(build_failed_object(reason = paste0("pilot did not converge: ",
                                                  pilot$message %||% ""),
                                pilot = pilot))
  }

  # Stabilise the pilot Hessian (symmetrise + eigenvalue floor) so the LSA
  # surrogate is a proper convex quadratic.  See stabilize_hessian() header
  # for rationale.  All downstream routines (compute_lambda_max,
  # compute_lambda_rho_max, cd_lsa) use this stabilised H.
  pilot$H <- stabilize_hessian(pilot$H)

  # --- Penalty factors ---
  pen_factors <- make_heck_penalty_factors(pilot$theta_hat, prep$NXS, prep$NXO,
                                             penalty, gamma_alasso, con$eps_alasso,
                                             unpenalized_idx = unpen_idx)

  # --- Dispatch on lambda / lambda.rho ---
  if (is.null(lambda) && is.null(lambda.rho)) {
    # 2D grid.
    #
    # CALIBRATION STRATEGY: Because the KKT boundary for one penalty
    # depends on the value of the other, a single compute_lambda_max
    # computed only at lambda.rho = 0 is NOT guaranteed to bound the
    # true KKT boundary at arbitrary lambda.rho.  We use a CONSERVATIVE
    # UPPER-ENVELOPE calibration: compute each bound at both extremes
    # of the other axis (eta free vs eta = 0 for lambda_max; reg-coefs
    # free vs reg-coefs at 0 for lambda_rho_max) and take the max.
    #
    # This endpoint-envelope calibration is substantially better than a
    # single-endpoint bound -- it correctly handles the case where one
    # endpoint gives a larger boundary than the other -- and it appears
    # to cover the relevant KKT boundaries in practice.  We do NOT
    # claim it is a provably tight upper bound for every intermediate
    # lambda.rho slice; a full joint KKT derivation would require
    # solving a coupled system for each slice.  The endpoint envelope
    # is a practical compromise that gives a noticeably better 2D grid
    # than the single-endpoint heuristic used earlier.
    #
    # lambda_max_A: at lambda.rho = 0 slice (eta free in the trial)
    # lambda_max_B: at lambda.rho = lrho_max slice (eta = 0 in the trial)
    # lambda_rho_max_A: at lambda = 0 slice (reg coefs free)
    # lambda_rho_max_B: at lambda = lmax slice (reg coefs at 0)
    eta_idx <- length(pilot$theta_hat)
    pen_idx_main <- which(pen_factors$w_main > 0)
    lmax_A <- compute_lambda_max(pilot$theta_hat, pilot$H, pen_factors$w_main,
                                    prep$NXS, prep$NXO)
    lmax_B <- compute_lambda_max(pilot$theta_hat, pilot$H, pen_factors$w_main,
                                    prep$NXS, prep$NXO,
                                    extra_zero_idx = eta_idx)
    lmax   <- max(lmax_A, lmax_B)

    lrho_max_A <- compute_lambda_rho_max(pilot$theta_hat, pilot$H,
                                            pen_factors$w_eta)
    lrho_max_B <- compute_lambda_rho_max(pilot$theta_hat, pilot$H,
                                            pen_factors$w_eta,
                                            extra_zero_idx = pen_idx_main)
    lrho_max   <- max(lrho_max_A, lrho_max_B) * con$lambda_rho_max_factor

    lambda_grid <- build_lambda_grid(lmax, con$n_lambda, con$lambda_ratio)
    lrho_grid <- build_lambda_grid(lrho_max, con$n_lambda_rho, con$lambda_ratio)
    # Prepend 0 so that we include the lambda.rho = 0 slice (Model 3a).
    lrho_grid <- c(0, lrho_grid)

    best_overall <- NULL
    best_overall_ic <- Inf
    all_path_info <- vector("list", length(lrho_grid))

    for (i in seq_along(lrho_grid)) {
      lrho_i <- lrho_grid[i]
      path_i <- fit_lasso_path(pilot$theta_hat, pilot$H, lambda_grid, lrho_i,
                                pen_factors, control = con)
      ic_i <- compute_ic(path_i$theta_path, prep, Model, pen_factors, lrho_i,
                          crit, con$tol_zero, converged = path_i$converged)
      all_path_info[[i]] <- list(lambda.rho = lrho_i, ic = ic_i, path = path_i)
      if (!is.na(ic_i$best)) {
        if (ic_i$ic[ic_i$best] < best_overall_ic) {
          best_overall_ic <- ic_i$ic[ic_i$best]
          best_overall <- list(
            theta_std = path_i$theta_path[, ic_i$best],
            lambda = lambda_grid[ic_i$best],
            lambda.rho = lrho_i,
            iterations = path_i$iterations[ic_i$best],
            ic = best_overall_ic,
            df = ic_i$df[ic_i$best],
            nll = ic_i$nll[ic_i$best]
          )
        }
      }
    }
    if (is.null(best_overall)) {
      warning("No converged fit found across the 2D grid; returning failed object.")
      return(build_failed_object("No fit on 2D grid converged.", pilot = pilot,
                                  path_grid = all_path_info,
                                  start_values = pilot$theta_hat))
    }
    # Build a flat grid data frame for the $path slot.
    grid_df <- do.call(rbind, lapply(all_path_info, function(p) {
      data.frame(lambda.rho = p$lambda.rho, lambda = lambda_grid,
                  ic = p$ic$ic, df = p$ic$df, nll = p$ic$nll)
    }))
    path_result <- list(grid = grid_df, per_rho = all_path_info)

  } else if (is.null(lambda) && !is.null(lambda.rho)) {
    # 1D grid on main lambda, lambda.rho fixed scalar.  For the single
    # specified lambda.rho we can compute lambda_max exactly: if
    # lambda.rho is 0 then compute_lambda_max with eta free is correct;
    # otherwise we take the max of the two endpoint computations as a
    # safe upper bound (covers the case where eta is at or above its
    # own shrinkage boundary for the supplied lambda.rho).
    if (length(lambda.rho) != 1L) stop("lambda.rho must be a scalar or NULL.")
    eta_idx <- length(pilot$theta_hat)
    lmax_A <- compute_lambda_max(pilot$theta_hat, pilot$H, pen_factors$w_main,
                                    prep$NXS, prep$NXO)
    lmax_B <- compute_lambda_max(pilot$theta_hat, pilot$H, pen_factors$w_main,
                                    prep$NXS, prep$NXO,
                                    extra_zero_idx = eta_idx)
    lmax   <- max(lmax_A, lmax_B)
    lambda_grid <- build_lambda_grid(lmax, con$n_lambda, con$lambda_ratio)
    path <- fit_lasso_path(pilot$theta_hat, pilot$H, lambda_grid, lambda.rho,
                            pen_factors, control = con)
    ic_info <- compute_ic(path$theta_path, prep, Model, pen_factors, lambda.rho,
                           crit, con$tol_zero, converged = path$converged)
    if (is.na(ic_info$best)) {
      warning("No converged fit found on lambda grid; returning failed object.")
      return(build_failed_object("No fit on lambda grid converged.", pilot = pilot,
                                  start_values = pilot$theta_hat))
    }
    k <- ic_info$best
    best_overall <- list(
      theta_std = path$theta_path[, k], lambda = lambda_grid[k],
      lambda.rho = lambda.rho,
      iterations = path$iterations[k], ic = ic_info$ic[k],
      df = ic_info$df[k], nll = ic_info$nll[k]
    )
    grid_df <- data.frame(lambda.rho = lambda.rho, lambda = lambda_grid,
                           ic = ic_info$ic, df = ic_info$df, nll = ic_info$nll)
    path_result <- list(grid = grid_df, per_rho = list(list(lambda.rho = lambda.rho,
                                                              ic = ic_info, path = path)))

  } else if (!is.null(lambda) && is.null(lambda.rho)) {
    # 1D grid on lambda.rho, lambda fixed
    if (length(lambda) != 1L) stop("lambda must be a scalar or NULL.")
    # 1D grid on lambda.rho, lambda fixed scalar.  For the supplied
    # fixed lambda, regression coefs may be anywhere from their MLE
    # values (at lambda = 0) down to zero (at lambda >= lambda_max), so
    # we compute lrho_max at both endpoints and take the max for a
    # conservative upper-envelope bound that covers both extremes.
    pen_idx_main <- which(pen_factors$w_main > 0)
    lrho_max_A <- compute_lambda_rho_max(pilot$theta_hat, pilot$H,
                                            pen_factors$w_eta)
    lrho_max_B <- compute_lambda_rho_max(pilot$theta_hat, pilot$H,
                                            pen_factors$w_eta,
                                            extra_zero_idx = pen_idx_main)
    lrho_max   <- max(lrho_max_A, lrho_max_B) * con$lambda_rho_max_factor
    lrho_grid <- c(0, build_lambda_grid(lrho_max, con$n_lambda_rho, con$lambda_ratio))
    all_path_info <- vector("list", length(lrho_grid))
    best_overall <- NULL; best_overall_ic <- Inf
    for (i in seq_along(lrho_grid)) {
      lrho_i <- lrho_grid[i]
      path_i <- fit_lasso_path(pilot$theta_hat, pilot$H, lambda, lrho_i,
                                pen_factors, control = con)
      ic_i <- compute_ic(path_i$theta_path, prep, Model, pen_factors, lrho_i,
                          crit, con$tol_zero, converged = path_i$converged)
      all_path_info[[i]] <- list(lambda.rho = lrho_i, ic = ic_i, path = path_i)
      if (!is.na(ic_i$best) && ic_i$ic[ic_i$best] < best_overall_ic) {
        best_overall_ic <- ic_i$ic[ic_i$best]
        best_overall <- list(
          theta_std = path_i$theta_path[, ic_i$best],
          lambda = lambda, lambda.rho = lrho_i,
          iterations = path_i$iterations[ic_i$best],
          ic = best_overall_ic, df = ic_i$df[ic_i$best],
          nll = ic_i$nll[ic_i$best]
        )
      }
    }
    if (is.null(best_overall)) {
      warning("No converged fit found on lambda.rho grid; returning failed object.")
      return(build_failed_object("No fit on lambda.rho grid converged.",
                                  pilot = pilot, start_values = pilot$theta_hat))
    }
    grid_df <- do.call(rbind, lapply(all_path_info, function(p) {
      data.frame(lambda.rho = p$lambda.rho, lambda = lambda,
                  ic = p$ic$ic, df = p$ic$df, nll = p$ic$nll)
    }))
    path_result <- list(grid = grid_df, per_rho = all_path_info)

  } else {
    # Single (lambda, lambda.rho) fit
    if (length(lambda) != 1L) stop("For a single fit, lambda must be scalar.")
    if (length(lambda.rho) != 1L) stop("For a single fit, lambda.rho must be scalar.")
    res <- cd_lsa(pilot$theta_hat, pilot$H, lambda, lambda.rho, pen_factors,
                   theta_start = pilot$theta_hat, control = con)
    if (!isTRUE(res$converged)) {
      warning("Coordinate descent did not converge at requested (lambda, lambda.rho).")
      return(build_failed_object("Coordinate descent failed to converge.",
                                  pilot = pilot, lambda_used = lambda,
                                  lambda_rho_used = lambda.rho,
                                  start_values = pilot$theta_hat))
    }
    # Evaluate IC at this single point
    ic_info <- compute_ic(matrix(res$theta, ncol = 1L), prep, Model, pen_factors,
                           lambda.rho, crit, con$tol_zero)
    if (is.na(ic_info$best)) {
      warning("Single fit produced non-finite likelihood; returning failed object.")
      return(build_failed_object("Non-finite likelihood at single fit.",
                                  pilot = pilot, lambda_used = lambda,
                                  lambda_rho_used = lambda.rho,
                                  start_values = pilot$theta_hat))
    }
    best_overall <- list(theta_std = res$theta, lambda = lambda,
                         lambda.rho = lambda.rho, iterations = res$iterations,
                         ic = ic_info$ic[1L], df = ic_info$df[1L],
                         nll = ic_info$nll[1L])
    grid_df <- data.frame(lambda.rho = lambda.rho, lambda = lambda,
                           ic = ic_info$ic, df = ic_info$df, nll = ic_info$nll)
    path_result <- list(grid = grid_df, per_rho = NULL)
  }

  # --- Back-transform selected coefficients ---
  NXS <- prep$NXS; NXO <- prep$NXO
  th_std <- best_overall$theta_std
  gamma_std <- th_std[1:NXS];  names(gamma_std) <- colnames(prep$XS)
  beta_std  <- th_std[(NXS + 1):(NXS + NXO)]; names(beta_std) <- colnames(prep$XO)
  eta_std   <- th_std[npar]
  gamma_raw <- backtransform_coefs(gamma_std, prep$scaling, "S")
  beta_raw  <- backtransform_coefs(beta_std,  prep$scaling, "O")
  rho_hat   <- tanh(eta_std)

  # --- In-sample predictions on the original scale ---
  # Linear predictors via original-scale coefficients applied to original X.
  # Rebuild original-scale design matrices.
  mf_sel <- model.frame(selection, data, na.action = na.pass)
  mf_out <- model.frame(outcome,   data, na.action = na.pass)
  XS_raw <- model.matrix(selection, mf_sel)
  XO_raw <- model.matrix(outcome,   mf_out)

  k1 <- drop(XO_raw %*% beta_raw)
  k2 <- drop(XS_raw %*% gamma_raw)
  pS <- pnorm(k2)

  if (Model == "Normal") {
    pJ <- pbivnorm(k1, k2, rho_hat)
  } else {
    # AMH joint P(S=1, Y=1) = C(Phi(k2), Phi(k1); rho) in closed form.
    u <- pnorm(k2); v <- pnorm(k1)
    pJ <- u * v / (1 - rho_hat * (1 - u) * (1 - v))
  }
  fitted.values <- ifelse(pS > 1e-8, pmin(pmax(pJ / pS, 0), 1), NA_real_)

  structure(list(
    call = cl,
    coefficients = list(selection = gamma_raw, outcome = beta_raw),
    std_coefficients = list(selection = gamma_std, outcome = beta_std),
    eta = eta_std, rho = rho_hat,
    lambda = best_overall$lambda, lambda.rho = best_overall$lambda.rho,
    logLik = -best_overall$nll, ic = best_overall$ic, df = best_overall$df,
    crit = crit,
    converged = TRUE, iterations = best_overall$iterations,
    fitted.values = fitted.values,
    linear.predictors = list(outcome = k1, selection = k2),
    formula = list(selection = selection, outcome = outcome),
    scaling = prep$scaling, nobs = prep$n, npar = npar,
    penalty = penalty, Model = Model, standardize = isTRUE(standardize),
    pilot = pilot, path = path_result, prep = prep,
    gamma_alasso = gamma_alasso,
    unpenalized = list(selection = unpenalized.S, outcome = unpenalized.O),
    fit_spec = fit_spec
  ), class = "HeckSelect")
}

# =======================================================================
# S3 METHODS
# =======================================================================

#' Print a HeckSelect fit
#' @param x a HeckSelect fit.
#' @param ... ignored.
#' @return invisibly returns the fit.
#' @export
print.HeckSelect <- function(x, ...) {
  cat("\nHeckSelect fit\n")
  cat("Penalty:       ", x$penalty, "\n")
  cat("Model:         ", x$Model,   "\n")
  cat("Criterion:     ", toupper(x$crit), "\n")
  cat("Converged:     ", x$converged, "\n")
  if (!isTRUE(x$converged)) {
    cat("  (fit did not converge; coefficients and predictions are NA)\n")
    if (!is.null(x$message)) cat("  Reason:", x$message, "\n")
    invisible(x); return(invisible(x))
  }
  cat(sprintf("%s:           %.4f  (df = %.2f)\n", toupper(x$crit), x$ic, x$df))
  cat(sprintf("rho:           %.4f  (eta = %.4f)\n", x$rho, x$eta))
  cat(sprintf("lambda:        %.6g\n", x$lambda))
  cat(sprintf("lambda.rho:    %.6g\n", x$lambda.rho))
  # Identify coordinates that were actually subject to the main-lambda
  # penalty.  Intercepts (always unpenalized) and user-specified
  # unpenalized.S/unpenalized.O variables must be excluded from both the
  # numerator (non-zero count) and the denominator (total penalized).
  sel_nms <- names(x$std_coefficients$selection)
  out_nms <- names(x$std_coefficients$outcome)
  sel_unpen_user <- x$unpenalized$selection %||% character()
  out_unpen_user <- x$unpenalized$outcome   %||% character()
  sel_pen <- setdiff(sel_nms, c("(Intercept)", sel_unpen_user))
  out_pen <- setdiff(out_nms, c("(Intercept)", out_unpen_user))
  n_nz_sel <- sum(abs(x$std_coefficients$selection[sel_pen]) > 1e-8)
  n_nz_out <- sum(abs(x$std_coefficients$outcome[out_pen])   > 1e-8)
  cat(sprintf("Non-zero penalized coefs:  selection %d / %d, outcome %d / %d\n",
              n_nz_sel, length(sel_pen),
              n_nz_out, length(out_pen)))
  if (length(sel_unpen_user) > 0 || length(out_unpen_user) > 0) {
    cat("Unpenalized (user-specified): selection: ",
        if (length(sel_unpen_user) > 0)
          paste(sel_unpen_user, collapse = ", ") else "(none)",
        "; outcome: ",
        if (length(out_unpen_user) > 0)
          paste(out_unpen_user, collapse = ", ") else "(none)",
        "\n", sep = "")
  }
  invisible(x)
}

#' Summarise a HeckSelect fit
#' @param object a HeckSelect fit.
#' @param ... ignored.
#' @return invisibly returns the fit.
#' @export
summary.HeckSelect <- function(object, ...) {
  cat("\nHeckSelect fit summary\n")
  cat("Call:\n  "); print(object$call); cat("\n")
  cat("Penalty:       ", object$penalty, "\n")
  cat("Model:         ", object$Model,   "\n")
  cat("Criterion:     ", toupper(object$crit), "\n")
  cat("Converged:     ", object$converged, "\n")
  if (!isTRUE(object$converged)) {
    cat("\n(fit did not converge; no coefficients to report)\n")
    if (!is.null(object$message)) cat("Reason:", object$message, "\n")
    return(invisible(object))
  }
  cat(sprintf("%s:           %.4f  (df = %.2f)\n",
              toupper(object$crit), object$ic, object$df))
  cat(sprintf("rho:           %.4f  (eta = %.4f)\n", object$rho, object$eta))
  cat(sprintf("lambda:        %.6g\n", object$lambda))
  cat(sprintf("lambda.rho:    %.6g\n", object$lambda.rho))

  if (isTRUE(object$standardize)) {
    cat("\nNOTE: Internal standardisation was applied (glmnet convention).\n")
    cat("      Coefficients shown below are on the ORIGINAL covariate scale\n")
    cat("      (back-transformed).  $std_coefficients contains the coefficients\n")
    cat("      on the INTERNAL STANDARDISED scale -- these are NOT standard\n")
    cat("      errors.  (Standard errors are not computed for penalized fits;\n")
    cat("      use bootstrap if SEs are needed.)\n\n")
  } else {
    cat("\nNOTE: Fitted without internal standardisation (standardize = FALSE).\n\n")
  }
  cat("Selection equation coefficients:\n")
  print(round(object$coefficients$selection, 5))
  cat("\nOutcome equation coefficients:\n")
  print(round(object$coefficients$outcome, 5))
  cat(sprintf("\nDependence parameter:\n  eta = %.5f    rho = tanh(eta) = %.5f\n",
              object$eta, object$rho))
  invisible(object)
}

#' Extract coefficients from a HeckSelect fit
#' @param object a HeckSelect fit.
#' @param scale \code{"original"} (default) or \code{"standardized"}.
#' @param type \code{"all"} (default), \code{"selection"},
#'   \code{"outcome"}, or \code{"rho"}.
#' @param ... ignored.
#' @return numeric vector or list, depending on `type`.
#' @export
coef.HeckSelect <- function(object,
                             scale = c("original", "standardized"),
                             type  = c("all", "selection", "outcome", "rho"),
                             ...) {
  scale <- match.arg(scale)
  type  <- match.arg(type)
  coefs <- if (scale == "original") object$coefficients else object$std_coefficients
  switch(type,
    all       = list(selection = coefs$selection, outcome = coefs$outcome,
                       rho = unname(object$rho), eta = unname(object$eta)),
    selection = coefs$selection,
    outcome   = coefs$outcome,
    rho       = c(rho = unname(object$rho), eta = unname(object$eta)))
}

#' Log-likelihood of a HeckSelect fit
#' @param object a HeckSelect fit.
#' @param ... ignored.
#' @return an object of class "logLik".
#' @export
logLik.HeckSelect <- function(object, ...) {
  val <- object$logLik
  attr(val, "df")   <- object$df
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' Predict from a HeckSelect fit
#' @param object a HeckSelect fit.
#' @param newdata optional new data frame; if NULL, returns in-sample
#'   predictions.
#' @param type one of \code{"conditional"}, \code{"joint"},
#'   \code{"latent_outcome"}, or \code{"selection"}.
#' @param ... ignored.
#' @return numeric vector of predicted probabilities.
#' @export
predict.HeckSelect <- function(object, newdata = NULL,
                                 type = c("conditional", "joint",
                                          "latent_outcome", "selection",
                                          "marginal"), ...) {
  type <- match.arg(type)

  # "marginal" was ambiguous in earlier versions (returned Phi(k1), the
  # latent outcome margin).  The semantically correct observable marginal
  # in a sample-selection model is P(S=1, Y=1 | X), now type = "joint".
  # We error on "marginal" to force the user to pick one explicitly.
  if (type == "marginal") {
    stop('type = "marginal" is no longer supported because it was ambiguous. ',
         'Use type = "joint" for the observable marginal P(S = 1, Y = 1 | X), ',
         'or type = "latent_outcome" for the latent-equation probit margin ',
         'Phi(X_O * beta) that ignores selection.')
  }

  # Failed-fit handling: mirrors RidgeBPSS.
  if (!isTRUE(object$converged)) {
    if (!is.null(newdata)) {
      stop("Model did not converge; prediction on newdata is unavailable. ",
           "Check $converged and refit with different starting values or penalties.")
    }
    warning("Model did not converge; returning NA predictions.")
  }

  # Helper: compute joint probability P(S=1, Y=1 | X) for Normal or AMH.
  joint_prob <- function(k1, k2, rho_hat, model) {
    if (model == "Normal") {
      pbivnorm(k1, k2, rho_hat)
    } else {
      u <- pnorm(k2); v <- pnorm(k1)
      u * v / (1 - rho_hat * (1 - u) * (1 - v))
    }
  }

  if (is.null(newdata)) {
    k1 <- object$linear.predictors$outcome
    k2 <- object$linear.predictors$selection
    return(switch(type,
      conditional     = object$fitted.values,
      joint           = joint_prob(k1, k2, object$rho, object$Model),
      latent_outcome  = pnorm(k1),
      selection       = pnorm(k2)))
  }
  tt_sel <- delete.response(terms(object$formula$selection))
  tt_out <- delete.response(terms(object$formula$outcome))
  sc <- object$scaling
  mf_sel <- model.frame(tt_sel, newdata, na.action = na.pass, xlev = sc$xlev_S)
  mf_out <- model.frame(tt_out, newdata, na.action = na.pass, xlev = sc$xlev_O)
  XS_new <- model.matrix(tt_sel, mf_sel, contrasts.arg = sc$contrasts_S)
  XO_new <- model.matrix(tt_out, mf_out, contrasts.arg = sc$contrasts_O)

  if (!identical(colnames(XS_new), names(object$coefficients$selection)))
    stop("Selection design matrix columns in newdata do not match training: ",
         paste(colnames(XS_new), collapse = ","), " vs ",
         paste(names(object$coefficients$selection), collapse = ","))
  if (!identical(colnames(XO_new), names(object$coefficients$outcome)))
    stop("Outcome design matrix columns in newdata do not match training: ",
         paste(colnames(XO_new), collapse = ","), " vs ",
         paste(names(object$coefficients$outcome), collapse = ","))

  # Apply RAW coefficients to RAW newdata.  No re-standardisation.
  k1 <- drop(XO_new %*% object$coefficients$outcome)
  k2 <- drop(XS_new %*% object$coefficients$selection)
  rho_hat <- object$rho
  switch(type,
    conditional = {
      pS <- pnorm(k2)
      pJ <- joint_prob(k1, k2, rho_hat, object$Model)
      ifelse(pS > 1e-8, pmin(pmax(pJ / pS, 0), 1), NA_real_)
    },
    joint           = joint_prob(k1, k2, rho_hat, object$Model),
    latent_outcome  = pnorm(k1),
    selection       = pnorm(k2))
}

#' Plot method for HeckSelect fits
#' @param x a HeckSelect fit.
#' @param type one of \code{"coefficients"}, \code{"ic"}, \code{"nzero"},
#'   or \code{"rho"}.
#' @param facet one of \code{"separate"} (default) or \code{"combined"}.
#' @param xvar one of \code{"lambda"} or \code{"l1norm"}.
#' @param lambda.rho optional: for 2D fits, the slice to display.
#' @param ... passed to underlying graphics functions.
#' @return invisibly NULL.
#' @export
plot.HeckSelect <- function(x, type = c("coefficients", "ic", "nzero", "rho"),
                             facet = c("separate", "combined"),
                             xvar = c("lambda", "l1norm"),
                             lambda.rho = NULL,
                             ...) {
  type  <- match.arg(type)
  facet <- match.arg(facet)
  xvar  <- match.arg(xvar)

  # Resolve the displayed slice and slice-local annotation scalars via
  # the helper.  All markers and titles below MUST use these slice-local
  # values (lambda_slice, ic_slice, rho_slice), not the global x$lambda
  # / x$ic / x$rho, unless you want to reintroduce Reviewer 1's
  # cross-slice annotation bug.
  info <- .plot_slice_info(x, lambda.rho = lambda.rho)
  slice        <- info$slice
  k_slice      <- info$k_slice
  lambda_slice <- info$lambda_slice
  ic_slice     <- info$ic_slice
  rho_slice    <- info$rho_slice

  lg <- slice$path$lambda_grid
  path_mat <- slice$path$theta_path   # npar x n_lambda
  n_lambda <- length(lg)
  NXS <- x$prep$NXS; NXO <- x$prep$NXO
  eta_idx <- NXS + NXO + 1L

  # Identify regression-coefficient indices (non-intercept).  For display
  # in the coefficient path we want ALL regression coords (both penalized
  # and user-unpenalized) because the user still wants to see them.
  sel_idx <- if (NXS > 1L) 2:NXS                  else integer(0)
  out_idx <- if (NXO > 1L) (NXS + 2L):(NXS + NXO) else integer(0)
  reg_idx <- c(sel_idx, out_idx)

  # nzero accounting: count ONLY coords that were actually subject to
  # the main-lambda penalty.  This excludes intercepts (already excluded
  # via sel_idx / out_idx starting at 2 and NXS+2), AND excludes
  # user-specified unpenalized variables.  Eta is always excluded from
  # nzero; it has its own "rho" plot.  This matches the print() method's
  # accounting and matches the documented semantics of nzero.
  unpen_user_S <- x$unpenalized$selection %||% character()
  unpen_user_O <- x$unpenalized$outcome   %||% character()
  sel_pen_names <- setdiff(colnames(x$prep$XS)[-1L], unpen_user_S)
  out_pen_names <- setdiff(colnames(x$prep$XO)[-1L], unpen_user_O)
  sel_pen_idx <- match(sel_pen_names, colnames(x$prep$XS))
  out_pen_idx <- match(out_pen_names, colnames(x$prep$XO)) + NXS
  pen_idx_for_nzero <- c(sel_pen_idx, out_pen_idx)
  pen_idx_for_nzero <- pen_idx_for_nzero[!is.na(pen_idx_for_nzero)]

  nzero <- if (length(pen_idx_for_nzero) > 0L) {
    apply(path_mat[pen_idx_for_nzero, , drop = FALSE], 2,
          function(v) sum(abs(v) > 1e-8))
  } else {
    rep(0L, n_lambda)
  }

  # X-axis values: log(lambda) (standard) or L1 norm (glmnet alternative).
  xvals <- switch(xvar,
                   lambda = log(pmax(lg, .Machine$double.eps)),
                   l1norm = colSums(abs(path_mat[reg_idx, , drop = FALSE])))
  xlab_txt <- switch(xvar,
                      lambda = expression(log(lambda)),
                      l1norm = "L1 norm of coefficients")

  # Slice-local selected-lambda marker.
  x_selected <- switch(xvar,
                        lambda = log(pmax(lambda_slice, .Machine$double.eps)),
                        l1norm = xvals[k_slice])

  # Thin the top-axis nzero labels.  For a 100-point grid, drawing every
  # label produces unreadable overlap.  We draw labels only at grid
  # positions where nzero *changes*, plus the endpoints -- this gives
  # users the meaningful information (where does model size transition)
  # without clutter.  For small grids (<= 15 points) we label all.
  nzero_ticks <- if (n_lambda <= 15L) {
    seq_len(n_lambda)
  } else {
    # Endpoints plus any index where nzero differs from its left neighbour
    changes <- which(c(FALSE, diff(nzero) != 0))
    unique(c(1L, changes, n_lambda))
  }

  # Extract graphics-specific dots (lambda.rho is consumed above; any
  # other ... args go through to base plot/matplot/axis).
  dots <- list(...)
  if ("lambda.rho" %in% names(dots)) dots$lambda.rho <- NULL

  # Helper: call a plotting function with (args, dots) merged.
  do_plot <- function(fun, args) do.call(fun, c(args, dots))

  # Dispatch ------------------------------------------------------------
  if (type == "coefficients") {
    # Intercept-only fits (NXS == 1 AND NXO == 1) have no regression
    # coefficients to plot.  Error cleanly rather than hand a 0-column
    # matrix to matplot().  The other plot types (ic, nzero, rho)
    # remain meaningful for intercept-only fits and proceed normally.
    if (length(reg_idx) == 0L) {
      stop(paste("type = 'coefficients' requires at least one non-intercept",
                  "regression covariate.  This fit has intercept-only formulas",
                  "in both equations; try type = 'ic', 'nzero', or 'rho'",
                  "instead."))
    }
    if (facet == "separate" && length(sel_idx) > 0 && length(out_idx) > 0) {
      op <- par(mfrow = c(2, 1), mar = c(4.5, 4.3, 4, 1) + 0.1,
                oma = c(0, 0, 2, 0))
      on.exit(par(op))

      # Selection panel
      ymat_S <- path_mat[sel_idx, , drop = FALSE]
      do_plot(matplot, list(xvals, t(ymat_S), type = "l", lty = 1, lwd = 1.4,
                             xlab = xlab_txt,
                             ylab = "std. selection coef",
                             main = "Selection equation",
                             col = seq_along(sel_idx)))
      abline(v = x_selected, lty = 2, col = "grey50")
      abline(h = 0, lty = 3, col = "grey70")
      axis(3, at = xvals[nzero_ticks], labels = nzero[nzero_ticks],
           tick = FALSE, line = -0.5, cex.axis = 0.7, col.axis = "grey40")
      mtext("nzero (penalized regression coefs only)",
            side = 3, line = 1.4, cex = 0.7, col = "grey40")
      legend("topright", legend = rownames(path_mat)[sel_idx],
             col = seq_along(sel_idx), lty = 1, lwd = 1.4,
             bty = "n", cex = 0.7)

      # Outcome panel
      ymat_O <- path_mat[out_idx, , drop = FALSE]
      do_plot(matplot, list(xvals, t(ymat_O), type = "l", lty = 1, lwd = 1.4,
                             xlab = xlab_txt,
                             ylab = "std. outcome coef",
                             main = "Outcome equation",
                             col = seq_along(out_idx)))
      abline(v = x_selected, lty = 2, col = "grey50")
      abline(h = 0, lty = 3, col = "grey70")
      legend("topright", legend = rownames(path_mat)[out_idx],
             col = seq_along(out_idx), lty = 1, lwd = 1.4,
             bty = "n", cex = 0.7)

      mtext(sprintf("HeckSelect %s path (%s, lambda.rho = %.4g)",
                     toupper(x$penalty), x$Model, slice$lambda.rho),
            outer = TRUE, line = 0.3, cex = 1.05, font = 2)
    } else {
      # Combined or single-equation-only plot
      ymat <- path_mat[reg_idx, , drop = FALSE]
      do_plot(matplot, list(xvals, t(ymat), type = "l", lty = 1, lwd = 1.4,
                             xlab = xlab_txt,
                             ylab = "standardised coefficient",
                             main = sprintf(paste("HeckSelect %s path",
                                                    "(%s, lambda.rho = %.4g)"),
                                              toupper(x$penalty), x$Model,
                                              slice$lambda.rho),
                             col = seq_along(reg_idx)))
      abline(v = x_selected, lty = 2, col = "grey50")
      abline(h = 0, lty = 3, col = "grey70")
      axis(3, at = xvals[nzero_ticks], labels = nzero[nzero_ticks],
           tick = FALSE, line = -0.5, cex.axis = 0.7, col.axis = "grey40")
      mtext("nzero (penalized coefs)", side = 3, line = 1.4,
            cex = 0.7, col = "grey40")
      legend("topright", legend = rownames(path_mat)[reg_idx],
             col = seq_along(reg_idx), lty = 1, lwd = 1.4,
             bty = "n", cex = 0.7)
    }

  } else if (type == "ic") {
    ic_vals <- slice$ic$ic
    # Replace Inf (non-converged / non-finite) with NA for plotting
    ic_plot <- ifelse(is.finite(ic_vals), ic_vals, NA_real_)
    # Leave extra top-margin room for the nzero axis
    op <- par(mar = c(4.5, 4.3, 4.5, 1) + 0.1)
    on.exit(par(op))
    do_plot(plot, list(xvals, ic_plot, type = "b", pch = 19, cex = 0.6,
                        xlab = xlab_txt,
                        ylab = toupper(x$crit),
                        main = sprintf(paste("%s = %.3f at slice min",
                                               "(%s, lambda.rho = %.4g)"),
                                         toupper(x$crit), ic_slice, x$Model,
                                         slice$lambda.rho)))
    abline(v = x_selected, lty = 2, col = "red")
    axis(3, at = xvals[nzero_ticks], labels = nzero[nzero_ticks],
         tick = FALSE, line = -0.5, cex.axis = 0.7, col.axis = "grey40")
    mtext("nzero (penalized coefs)", side = 3, line = 1.3,
          cex = 0.7, col = "grey40")

  } else if (type == "nzero") {
    do_plot(plot, list(xvals, nzero, type = "s", lwd = 2,
                        xlab = xlab_txt,
                        ylab = "number of non-zero penalized coefficients",
                        main = sprintf("HeckSelect model size (%s, lambda.rho = %.4g)",
                                        x$Model, slice$lambda.rho)))
    abline(v = x_selected, lty = 2, col = "red")
    # Dot at slice-local selected point
    points(xvals[k_slice], nzero[k_slice], pch = 19, col = "red", cex = 1.2)

  } else if (type == "rho") {
    # rho along the path, reconstructed from theta_path[eta_idx, ]
    eta_path <- path_mat[eta_idx, ]
    rho_path <- tanh(eta_path)
    do_plot(plot, list(xvals, rho_path, type = "l", lwd = 2, ylim = c(-1, 1),
                        xlab = xlab_txt,
                        ylab = expression(hat(rho)),
                        main = sprintf("HeckSelect dependence path (%s, lambda.rho = %.4g)",
                                        x$Model, slice$lambda.rho)))
    abline(v = x_selected, lty = 2, col = "red")
    abline(h = 0,           lty = 3, col = "grey70")
    abline(h = c(-1, 1),    lty = 3, col = "grey85")
    # Slice-local selected rho (sits on the plotted curve by construction)
    points(x_selected, rho_slice, pch = 19, col = "red", cex = 1.2)
  }

  invisible(x)
}

#' @keywords internal
#' @noRd
.plot_slice_info <- function(x, lambda.rho = NULL) {
  if (!isTRUE(x$converged) || is.null(x$path) || is.null(x$path$per_rho)) {
    stop("No path information available.  This may be a single-point fit.")
  }
  per_rho <- x$path$per_rho
  available_lrho <- vapply(per_rho, function(p) p$lambda.rho, numeric(1))
  if (is.null(lambda.rho)) {
    lrho_target <- x$lambda.rho
  } else {
    if (length(lambda.rho) != 1L || !is.numeric(lambda.rho)) {
      stop("lambda.rho must be a single numeric value.")
    }
    lrho_target <- lambda.rho
  }
  matches <- vapply(available_lrho,
                     function(v) isTRUE(all.equal(v, lrho_target)),
                     logical(1))
  if (!any(matches)) {
    stop(sprintf(paste("Requested lambda.rho = %.6g not found in fitted path.",
                        "Available slices: %s"),
                  lrho_target,
                  paste(sprintf("%.6g", available_lrho), collapse = ", ")))
  }
  sel <- which(matches)[1L]
  slice <- per_rho[[sel]]

  # Slice-local best index.  If the slice has no converged path point
  # to annotate (every IC is Inf), we do NOT silently fall back to the
  # global x$lambda -- that would reintroduce the cross-slice marker
  # bug in a subtler form.  We error cleanly instead.
  k_slice <- slice$ic$best
  if (is.na(k_slice)) {
    stop(sprintf(paste("Requested lambda.rho = %.6g slice has no converged",
                        "path points; nothing to annotate.  This is an edge",
                        "case usually caused by overriding to a pathological",
                        "slice where every path point failed IC evaluation."),
                  lrho_target))
  }

  lg <- slice$path$lambda_grid
  NXS <- x$prep$NXS; NXO <- x$prep$NXO
  eta_idx <- NXS + NXO + 1L

  list(
    slice        = slice,
    sel          = sel,
    k_slice      = k_slice,
    lambda_slice = lg[k_slice],
    ic_slice     = slice$ic$ic[k_slice],
    rho_slice    = tanh(slice$path$theta_path[eta_idx, k_slice])
  )
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
bootValidate.HeckSelect <- function(object, data, mboot = 200L,
                                       seed = NULL, verbose = FALSE,
                                       retune = TRUE, control = NULL,
                                       ...) {
  if (!isTRUE(object$converged))
    stop("Cannot validate a non-converged fit; check object$converged.")
  if (is.null(object$fit_spec))
    stop("Object lacks $fit_spec; cannot replay fitting procedure. ",
         "This usually means the fit was created by an older version ",
         "of HeckSelect; refit and try again.")
  n_orig <- nrow(data)
  if (n_orig != length(object$prep$YS))
    stop("'data' has ", n_orig, " rows but original fit had ",
         length(object$prep$YS), " observations.")

  # WARN about extra ... args.  Previously ... was forwarded into the
  # refit closure, which caused a silent double-argument bug: if the
  # user passed control = X to bootValidate(), it collided with the
  # control = fs$control_input the closure supplies internally, R
  # raised "formal argument 'control' matched by multiple actual
  # arguments", and the tryCatch in .bootValidate_engine() swallowed
  # the error and returned NULL for every refit.  All "successful"
  # bootstrap runs turned into empty failure objects.
  #
  # The fix: ... is not forwarded to refits.  All refit-influencing
  # arguments must come from the fitted object's fit_spec or from the
  # explicit `control` argument below.  Any extra ... args are
  # diagnostic noise and we warn.
  if (...length() > 0L) {
    unused <- ...names()
    warning("Unused arguments to bootValidate.HeckSelect: ",
             paste(sprintf("'%s'", unused), collapse = ", "),
             ".  Refit configuration comes from object$fit_spec and ",
             "the explicit 'control' argument; any other arguments ",
             "are ignored.  If you wanted to override the control ",
             "list, use control = list(...) explicitly.")
  }

  # Replay the ORIGINAL user-supplied fitting specification.  Using the
  # input values (not the selected ones) is what makes the bootstrap
  # validly estimate optimism for the ENTIRE fitting procedure.  For
  # example, if the original call was HeckSelect(lambda = NULL,
  # lambda.rho = NULL) -- a 2D-tuned fit -- each bootstrap refit must
  # also do 2D tuning.  If the original pinned lambda.rho = 10 and
  # re-tuned lambda, each refit must too.  Replaying the SELECTED
  # values instead would collapse tuning variance and understate
  # optimism.
  fs <- object$fit_spec
  selection <- object$formula$selection
  outcome   <- object$formula$outcome

  # Build the control list for refits: start from fs$control_input,
  # and merge any user-supplied overrides.  control = NULL (default)
  # means "use fs$control_input unchanged".  Named entries in the
  # user's control override the corresponding fs entries; fs entries
  # the user doesn't override are preserved.
  if (!is.null(control) && !is.list(control))
    stop("'control' must be a list or NULL.")
  control_boot <- fs$control_input
  if (!is.null(control)) {
    control_boot[names(control)] <- control
  }

  # Refit closure.  Takes ONLY boot_dat (no ... forwarding, no
  # possibility of double-argument conflicts).  retune = TRUE replays
  # fs (input values).  retune = FALSE pins to the realised selected
  # (lambda, lambda.rho) from the original fit -- useful for
  # sensitivity analysis of the specific chosen model, but NOT the
  # standard bootstrap-optimism setup.
  refit_fn <- function(boot_dat) {
    if (retune) {
      HeckSelect(selection, outcome, data = boot_dat,
                  lambda     = fs$lambda_input,
                  lambda.rho = fs$lambda_rho_input,
                  Model = fs$Model, penalty = fs$penalty, crit = fs$crit,
                  gamma_alasso = fs$gamma_alasso,
                  standardize  = fs$standardize,
                  unpenalized.S = fs$unpenalized.S,
                  unpenalized.O = fs$unpenalized.O,
                  control = control_boot)
    } else {
      HeckSelect(selection, outcome, data = boot_dat,
                  lambda     = object$lambda,
                  lambda.rho = object$lambda.rho,
                  Model = fs$Model, penalty = fs$penalty, crit = fs$crit,
                  gamma_alasso = fs$gamma_alasso,
                  standardize  = fs$standardize,
                  unpenalized.S = fs$unpenalized.S,
                  unpenalized.O = fs$unpenalized.O,
                  control = control_boot)
    }
  }

  # Extractor closure: given a fitted object, return (y, p) pairs on
  # the observations where the outcome is observed AND fitted.values
  # is finite.
  extract_yp <- function(fit) {
    sel <- which(fit$prep$YS == 1L)
    if (length(sel) == 0L) return(list(y = numeric(0), p = numeric(0)))
    p <- fit$fitted.values[sel]
    y <- fit$prep$YO[sel]
    ok <- !is.na(p)
    list(y = y[ok], p = p[ok])
  }

  # CRITICAL: compute eval_orig ONCE, and use it consistently for the
  # apparent target, apparent predictions, and the predict_on closure.
  # An earlier version used sel_orig (length = N_selected) for predict
  # and apparent_y/p after ok_app (length = N_selected - N_NA_preds),
  # which misaligned pairing inside the engine when fitted.values had
  # any NAs.  Fixing alignment once here guarantees correctness.
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

  # Attach HeckSelect-specific return fields
  out$coef       <- object$coefficients
  out$lambda     <- object$lambda
  out$lambda.rho <- object$lambda.rho
  out
}


