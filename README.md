# HeckmanSelect

Regularized estimation for the bivariate probit model with 
sample-selection (BPSS).

## Overview

`HeckmanSelect` implements three complementary penalized estimators 
for the bivariate probit model with sample selection, suitable 
for reject-inference applications in credit scoring and related 
selection-bias problems in binary outcomes.

- **`HeckSelect()`** — lasso (Tibshirani 1996) and adaptive-lasso 
  (Zou 2006) penalization via a Least-Squares Approximation (LSA) 
  surrogate with coordinate descent and warm starts. Supports 
  Normal-error and Ali–Mikhail–Haq (AMH) copula dependence structures.
- **`RidgeBPSS()`** — ridge penalization via direct trust-region 
  optimization of the penalized negative log-likelihood, with 
  separate penalty parameters for regression coefficients and the 
  selectivity parameter.

Both fitters support automatic tuning of the penalty parameter(s) 
via BIC, AIC, or GCV, and return S3 objects with `print`, `summary`, 
`coef`, `predict`, and `logLik` methods.

Bootstrap optimism-corrected internal validation is provided via 
the `bootValidate()` S3 generic, which reports paired-optimism-
corrected AUC, AUPR (area under precision-recall curve), Brier 
score, ECE (expected calibration error), and MCE (maximum 
calibration error). Paired computation is used because AUC and 
AUPR can fail when a bootstrap sample misses one class, while 
Brier, ECE, and MCE remain computable — a naive 
`mean(training) - mean(test)` would mix different replicate sets 
per metric.

## Included data

The package bundles the **American Express credit card dataset** 
(Greene 1998): 13,444 applications with 56 variables covering 
demographics, credit history, spending patterns, and ZIP-code-level 
socioeconomic indicators. This is the standard benchmark dataset 
for reject-inference research. See `?AmEx` for full variable 
documentation.

## Installation

```r
# install.packages("devtools")
devtools::install_github("EOgundimu300/HeckmanSelect")
```

To install a specific tagged version:

```r
devtools::install_github("EOgundimu300/HeckmanSelect@v2.0.0")
```

## Quick example (synthetic data)

```r
library(HeckmanSelect)

set.seed(1)
n <- 500
x1 <- rnorm(n); x2 <- rnorm(n); Z <- rnorm(n); V <- rnorm(n)
S <- as.integer(0.2 + 0.5*x1 - 0.4*x2 - 0.5*Z + rnorm(n) > 0)
Y <- as.integer(0.1 + 0.6*x1 - 0.3*x2 + rnorm(n) > 0)
Y[S == 0] <- 0L
dat <- data.frame(S = S, Y_obs = Y, x1 = x1, x2 = x2, Z = Z, V = V)

# Adaptive-lasso fit
fit_hs <- HeckSelect(S ~ x1 + x2 + Z, Y_obs ~ x1 + x2 + V,
                      data = dat, penalty = "alasso", Model = "Normal")
summary(fit_hs)

# Ridge fit
fit_rb <- RidgeBPSS(S ~ x1 + x2 + Z, Y_obs ~ x1 + x2 + V,
                     data = dat, lambda = 0.1, lambda.rho = 0)
summary(fit_rb)

# Bootstrap validation (requires the PRROC package)
val <- bootValidate(fit_hs, data = dat, mboot = 200, seed = 1)
val$resu
```

## Example on AmEx

```r
library(HeckmanSelect)
data(AmEx)

outcome   <- DEFAULT  ~ AGE + ACADMOS + ADEPCNT + MAJORDRG + 
             MINORDRG + OWNRENT + INCOME + SELFEMPL
selection <- CARDHLDR ~ AGE + ACADMOS + ADEPCNT + MAJORDRG + 
             MINORDRG + OWNRENT + INCOME + SELFEMPL + BANKSAV + BANKCH

fit <- HeckSelect(selection, outcome, data = AmEx,
                   penalty = "alasso", Model = "AMH", crit = "bic")
summary(fit)

# Paired-optimism bootstrap validation
val <- bootValidate(fit, data = AmEx, mboot = 200, seed = 1)
val$resu
```

## Methodology

The methods implemented in this package are described in:

- Ogundimu, E. O. (2019). Prediction of default probability by using 
  statistical models for rare events. *Journal of the Royal 
  Statistical Society: Series A*, 182(4), 1143–1162.
- Ogundimu, E. O. (2022). On Lasso and adaptive Lasso for non-random 
  sample in credit scoring. *Statistical Modelling*, 22(6), 519–542.


## Citation

To cite `HeckmanSelect` in publications, use:

```r
citation("HeckmanSelect")
```

## License

MIT © Emmanuel Ogundimu. See `LICENSE` for details.

## Contributing

Issues and pull requests are welcome at 
[https://github.com/EOgundimu300/HeckmanSelect/issues](https://github.com/EOgundimu300/HeckmanSelect/issues).
