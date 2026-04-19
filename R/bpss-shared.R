# -----------------------------------------------------------------------
# Shared bootstrap infrastructure for HeckSelect and RidgeBPSS
#
# This file defines:
#   - bootValidate: S3 generic (exported)
#   - .calibration_metrics: ECE / MCE helper (internal)
#   - .prroc_metrics: per-(y, prob) metrics via PRROC (internal)
#   - .summarize_boot_metrics: paired-optimism summariser (internal)
#   - .bootValidate_engine: class-agnostic bootstrap engine (internal)
#
# The HeckSelect- and RidgeBPSS-specific methods live alongside their
# respective fitter functions in R/HeckSelect.R and R/RidgeBPSS.R.
# -----------------------------------------------------------------------

#' @importFrom stats quantile
NULL


# =======================================================================
# S3 GENERIC (exported)
# =======================================================================

#' Bootstrap internal validation for penalized BPSS fits
#'
#' S3 generic implementing Harrell's bootstrap-optimism correction for
#' apparent in-sample predictive performance.  Refits the model on
#' `mboot` bootstrap samples, evaluates each refit on both the
#' bootstrap sample (training-optimistic) and the original data
#' (test), and reports optimism-corrected performance as
#' `apparent - optimism`, where `optimism` is the mean of paired
#' differences (training - test) over replicates where both sides
#' are finite.
#'
#' Paired computation is essential because AUC and AUPR can fail
#' when a bootstrap sample misses one class, while Brier, ECE, and
#' MCE remain computable; `colMeans(training) - colMeans(test)`
#' would mix different replicate sets per metric and yield a biased
#' optimism estimate.
#'
#' Metrics reported: AUC (area under ROC), AUPR (area under
#' precision-recall curve -- note this is NOT "partial AUC"), Brier
#' score, ECE (expected calibration error), MCE (maximum calibration
#' error).  PRROC is used for AUC / AUPR computation; install the
#' PRROC package before calling.
#'
#' Methods are provided for objects of class \code{"HeckSelect"} and
#' \code{"RidgeBPSS"}.
#'
#' @param object a fitted model object of class \code{"HeckSelect"}
#'   or \code{"RidgeBPSS"}.
#' @param ... method-specific arguments; see
#'   \code{\link{bootValidate.HeckSelect}} and
#'   \code{\link{bootValidate.RidgeBPSS}}.
#'
#' @return a list with components:
#'   \describe{
#'     \item{\code{coef}, \code{lambda}, \code{lambda.rho}}{from the
#'       original fit (for reference).}
#'     \item{\code{nonconvergence}}{count of bootstrap replicates
#'       where every metric in the training matrix is NA.}
#'     \item{\code{resu}}{data frame with one row per metric, columns
#'       \code{index.orig} (apparent), \code{training} (bootstrap
#'       in-sample mean), \code{test} (original-data evaluation
#'       mean), \code{optimism} (paired mean of training - test),
#'       \code{index.corrected} (apparent - optimism), and
#'       \code{n} (metric-specific count of usable replicates).}
#'     \item{\code{*_r}, \code{*test_r}}{per-replicate metric vectors
#'       (useful for plotting bootstrap distributions).}
#'   }
#'
#' @references
#' Harrell, F. E. (2015).  \emph{Regression Modeling Strategies}, 2nd
#' ed.  Springer, Chapter 5.
#'
#' @examples
#' \dontrun{
#' # HeckSelect example
#' fit_hs <- HeckSelect(S ~ x1 + x2, Y_obs ~ x1 + x2, data = dat,
#'                       penalty = "lasso", Model = "Normal")
#' val_hs <- bootValidate(fit_hs, data = dat, mboot = 100, seed = 1)
#' val_hs$resu
#'
#' # RidgeBPSS example
#' fit_rb <- RidgeBPSS(S ~ x1 + x2, Y_obs ~ x1 + x2, data = dat,
#'                      lambda = 0.1, lambda.rho = 0)
#' val_rb <- bootValidate(fit_rb, data = dat, mboot = 100, seed = 1)
#' val_rb$resu
#' }
#'
#' @export
bootValidate <- function(object, ...) UseMethod("bootValidate")


# =======================================================================
# INTERNAL HELPERS
# =======================================================================

#' Expected Calibration Error (ECE) and Maximum Calibration Error (MCE)
#'
#' Internal helper: computes calibration error between predicted
#' probabilities and binary outcomes via quantile-based binning.
#' Quantile-based binning via `cut()` + `unique(quantile(...))`
#' handles ties (common when lasso zeroes many coefficients and
#' predictions quantise) and unequal sample / bin divisibility.
#' ECE is the count-weighted mean of per-bin absolute gaps.
#'
#' @param y numeric/integer/logical vector of binary outcomes (0/1).
#' @param prob numeric vector in [0, 1], same length as `y`.
#' @param g integer number of bins, default 10.
#' @param group optional grouping vector for per-subgroup calibration.
#'
#' @return list with `ece`, `mce`, `bins` data frame, and optional
#'   `by_group` list.
#'
#' @keywords internal
#' @noRd
.calibration_metrics <- function(y, prob, g = 10, group = NULL) {
  # --- Input validation ---
  if (!is.numeric(y) && !is.integer(y) && !is.logical(y))
    stop("'y' must be numeric, integer, or logical.")
  if (!is.numeric(prob))
    stop("'prob' must be numeric.")
  if (length(y) != length(prob))
    stop("'y' and 'prob' must have the same length.")
  if (!is.null(group) && length(group) != length(y))
    stop("'group' must have the same length as 'y' when supplied.")
  if (length(g) != 1L || !is.numeric(g) || g < 1L || g != as.integer(g))
    stop("'g' must be a positive integer.")

  # Single-bin helper.  Computes binned ECE/MCE on a (y, prob) pair
  # using quantile-based cut with robust tie handling.
  one_ece <- function(y_, prob_, g_) {
    # Drop NAs
    ok <- !is.na(y_) & !is.na(prob_)
    y_   <- y_[ok]
    prob_ <- prob_[ok]
    n <- length(prob_)
    if (n == 0L) {
      return(list(ece = NA_real_, mce = NA_real_,
                   bins = data.frame(bin = integer(0), count = integer(0),
                                     mean_prob = numeric(0),
                                     obs_rate  = numeric(0),
                                     abs_gap   = numeric(0))))
    }
    # Quantile-based breaks.  If ties collapse cutpoints we end up with
    # fewer bins than g, which is the correct behaviour (can't have
    # more bins than unique predicted probabilities).
    probs_seq <- seq(0, 1, length.out = g_ + 1L)
    brks <- unique(stats::quantile(prob_, probs = probs_seq, na.rm = TRUE,
                                     names = FALSE))
    if (length(brks) < 2L) {
      # All predictions identical: one bin.
      gap <- abs(mean(y_ == 1) - mean(prob_))
      return(list(ece = gap, mce = gap,
                   bins = data.frame(bin = 1L, count = n,
                                     mean_prob = mean(prob_),
                                     obs_rate  = mean(y_ == 1),
                                     abs_gap   = gap)))
    }
    bin <- cut(prob_, breaks = brks, include.lowest = TRUE, labels = FALSE)
    # Per-bin summaries
    counts     <- tabulate(bin, nbins = length(brks) - 1L)
    mean_prob <- tapply(prob_, bin, mean)
    obs_rate  <- tapply(y_ == 1, bin, mean)
    abs_gap   <- abs(obs_rate - mean_prob)
    # Drop any empty bins (possible if ties collapsed mid-range)
    nonempty <- counts > 0L
    counts     <- counts[nonempty]
    mean_prob <- as.numeric(mean_prob)[nonempty]
    obs_rate  <- as.numeric(obs_rate)[nonempty]
    abs_gap   <- as.numeric(abs_gap)[nonempty]
    # ECE: count-weighted mean of absolute gaps (standard definition).
    ece <- sum(counts * abs_gap) / sum(counts)
    mce <- max(abs_gap)
    bins_df <- data.frame(bin = seq_along(counts),
                           count = counts,
                           mean_prob = mean_prob,
                           obs_rate  = obs_rate,
                           abs_gap   = abs_gap)
    list(ece = ece, mce = mce, bins = bins_df)
  }

  # Overall
  out <- one_ece(y, prob, g)

  # Per-group (if requested)
  if (!is.null(group)) {
    grp_fac <- as.factor(group)
    by_group <- lapply(levels(grp_fac), function(lv) {
      idx <- which(grp_fac == lv & !is.na(grp_fac))
      one_ece(y[idx], prob[idx], g)
    })
    names(by_group) <- levels(grp_fac)
    out$by_group <- by_group
  }

  out
}


#' Per-(y, prob) metrics: AUC, AUPR, Brier, ECE, MCE
#'
#' Internal helper used by `.bootValidate_engine()` in three places
#' (apparent, training, test).  AUC and AUPR require PRROC and both
#' classes present in `y`; they return NA if either condition fails.
#' Brier, ECE, and MCE are always computable.
#'
#' @param y numeric/integer vector of binary outcomes (0/1).
#' @param prob numeric vector of predicted probabilities.
#'
#' @return list with `auc`, `aupr`, `brier`, `ece`, `mce`.
#'
#' @keywords internal
#' @noRd
.prroc_metrics <- function(y, prob) {
  fg <- prob[y == 1]
  bg <- prob[y == 0]
  # PRROC needs at least one obs of each class for AUC/AUPR.  If the
  # bootstrap happens to miss one class entirely, return NAs for the
  # discrimination metrics but still compute Brier/ECE/MCE.
  auc <- aupr <- NA_real_
  if (length(fg) > 0L && length(bg) > 0L) {
    roc <- tryCatch(PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg),
                     error = function(e) NULL)
    pr  <- tryCatch(PRROC::pr.curve(scores.class0  = fg, scores.class1 = bg),
                     error = function(e) NULL)
    if (!is.null(roc)) auc  <- roc$auc
    if (!is.null(pr))  aupr <- pr$auc.integral
  }
  brier <- mean((prob - y)^2)
  cal   <- .calibration_metrics(y, prob, g = 10L)
  list(auc = auc, aupr = aupr, brier = brier,
       ece = cal$ece, mce = cal$mce)
}


#' Summarise paired-optimism bootstrap metrics
#'
#' Internal helper: given apparent-metric vector + per-replicate
#' training/test metric matrices, compute Harrell's paired optimism
#' correction:
#'   `optimism[m] = mean(boot[, m] - test[, m])` over rows where
#'   both are finite.
#' Reports also the `resu` data frame and `refit_failures` count
#' (replicates where the refit itself produced nothing, i.e. every
#' entry in the training row is NA).
#'
#' Factored out of `.bootValidate_engine()` so the paired-optimism
#' math is directly unit-testable independently of the refit/predict
#' machinery.
#'
#' @param apparent named numeric vector of apparent (in-sample)
#'   metric values.
#' @param boot_metrics matrix of per-replicate training metrics,
#'   one row per bootstrap sample, one column per metric.
#' @param test_metrics matrix of per-replicate test metrics, same
#'   shape and column names as `boot_metrics`.
#'
#' @return list with `resu` (data frame) and `refit_failures`
#'   (integer).
#'
#' @keywords internal
#' @noRd
.summarize_boot_metrics <- function(apparent, boot_metrics, test_metrics) {
  metric_names <- colnames(boot_metrics)
  if (is.null(metric_names) ||
      !identical(metric_names, colnames(test_metrics)))
    stop(".summarize_boot_metrics: boot_metrics and test_metrics must ",
         "have identical column names.")
  n_metric <- length(metric_names)

  training <- numeric(n_metric)
  test     <- numeric(n_metric)
  optimism <- numeric(n_metric)
  n_used   <- integer(n_metric)

  for (m in seq_len(n_metric)) {
    ok <- is.finite(boot_metrics[, m]) & is.finite(test_metrics[, m])
    n_used[m] <- sum(ok)
    if (n_used[m] > 0L) {
      training[m] <- mean(boot_metrics[ok, m])
      test[m]     <- mean(test_metrics[ok, m])
      optimism[m] <- mean(boot_metrics[ok, m] - test_metrics[ok, m])
    } else {
      training[m] <- test[m] <- optimism[m] <- NA_real_
    }
  }

  resu <- data.frame(
    index.orig      = apparent,
    training        = training,
    test            = test,
    optimism        = optimism,
    index.corrected = apparent - optimism,
    n               = n_used
  )
  resu <- round(resu, 4)
  rownames(resu) <- metric_names

  # "Refit failures" (reported as nonconvergence) = rows where ALL
  # metrics are NA in boot_metrics.  This is a narrower and more
  # accurate definition than "any NA" -- a converged fit whose
  # bootstrap-in-sample evaluation produced a single NA (e.g. one
  # class missing from the selected observations) is not a refit
  # failure, just a metric-computability artefact, which the
  # per-metric n column already reports.
  refit_failures <- sum(apply(boot_metrics, 1,
                                function(r) all(is.na(r))))

  list(resu = resu, refit_failures = refit_failures)
}


#' Class-agnostic bootstrap-validation engine
#'
#' Internal helper: given closures for refit, (y, p) extraction, and
#' predict-on-newdata, run `mboot` bootstrap iterations, collect
#' training and test metrics, and summarise via
#' `.summarize_boot_metrics()`.
#'
#' This engine is called by both `bootValidate.HeckSelect()` and
#' `bootValidate.RidgeBPSS()`, each of which supplies
#' class-specific closures that know how to refit their respective
#' models and extract predictions.
#'
#' @param refit_fn single-argument closure: takes a bootstrap-sample
#'   data frame, returns a fitted object.
#' @param extract_yp closure: takes a fitted object, returns a list
#'   with elements `y` and `p` (observed outcomes and predicted
#'   probabilities on the fit's in-sample selected observations).
#' @param predict_on closure: takes a fitted object and a new data
#'   frame, returns a vector of predictions aligned with
#'   `apparent_y` and `apparent_p`.
#' @param data the original training data frame.
#' @param apparent_y observed outcomes for the apparent (original)
#'   fit, aligned with `apparent_p`.
#' @param apparent_p predicted probabilities from the apparent fit.
#' @param mboot number of bootstrap samples.
#' @param seed optional integer for set.seed.
#' @param verbose logical; print progress dots and surface refit
#'   errors via `message()` if TRUE.
#'
#' @return list with `nonconvergence`, `resu`, and per-replicate
#'   metric vectors.
#'
#' @keywords internal
#' @noRd
.bootValidate_engine <- function(refit_fn, extract_yp, predict_on,
                                   data, apparent_y, apparent_p,
                                   mboot = 200L, seed = NULL,
                                   verbose = FALSE) {
  if (!requireNamespace("PRROC", quietly = TRUE))
    stop("Package 'PRROC' is required for bootValidate(); ",
         "install.packages('PRROC').")
  if (!is.null(seed)) set.seed(seed)

  metric_names <- c("auc", "aupr", "brier", "ece", "mce")
  n_metric <- length(metric_names)

  apparent <- unlist(.prroc_metrics(apparent_y, apparent_p))
  # (Names of apparent match metric_names.)
  apparent <- apparent[metric_names]

  boot_metrics <- matrix(NA_real_, nrow = mboot, ncol = n_metric,
                          dimnames = list(NULL, metric_names))
  test_metrics <- matrix(NA_real_, nrow = mboot, ncol = n_metric,
                          dimnames = list(NULL, metric_names))

  n_orig <- nrow(data)
  for (i in seq_len(mboot)) {
    iboot <- sample.int(n_orig, replace = TRUE)
    boot_dat <- data[iboot, , drop = FALSE]

    # Refit.  refit_fn is a single-argument closure built by the
    # method (no ... forwarding); any refit-influencing configuration
    # is baked into the closure from the original fit_spec and the
    # method's explicit control argument.  If verbose, we log error
    # messages on refit failure (previously errors were silently
    # swallowed, which hid the double-control bug for months).
    res <- tryCatch(suppressWarnings(refit_fn(boot_dat)),
                     error = function(e) {
                       if (verbose) {
                         message("bootValidate refit ", i, " errored: ",
                                 conditionMessage(e))
                       }
                       NULL
                     })
    if (is.null(res) || !isTRUE(res$converged)) {
      if (verbose) cat("x")
      next
    }

    # Training-optimistic (bootstrap in-sample) evaluation
    yp_b <- extract_yp(res)
    if (length(yp_b$p) > 0L) {
      boot_metrics[i, ] <- unlist(
        .prroc_metrics(yp_b$y, yp_b$p)
      )[metric_names]
    }

    # Test: bootstrap fit evaluated on the ORIGINAL data
    p_test <- tryCatch(predict_on(res, data),
                        error = function(e) NULL)
    if (!is.null(p_test)) {
      ok_t <- !is.na(p_test)
      if (any(ok_t)) {
        test_metrics[i, ] <- unlist(
          .prroc_metrics(apparent_y[ok_t], p_test[ok_t])
        )[metric_names]
      }
    }

    if (verbose) cat(".")
  }
  if (verbose) cat("\n")

  # Summarise via the standalone helper.  Factoring the paired-
  # optimism math out of the engine makes it unit-testable
  # independently of the refit/predict machinery (see Test 56).
  summ <- .summarize_boot_metrics(apparent, boot_metrics, test_metrics)

  list(
    nonconvergence = summ$refit_failures,
    resu           = summ$resu,
    auc_r          = boot_metrics[, "auc"],
    auctest_r      = test_metrics[, "auc"],
    aupr_r         = boot_metrics[, "aupr"],
    auprtest_r     = test_metrics[, "aupr"],
    brier_r        = boot_metrics[, "brier"],
    briertest_r    = test_metrics[, "brier"],
    ece_r          = boot_metrics[, "ece"],
    ecetest_r      = test_metrics[, "ece"],
    mce_r          = boot_metrics[, "mce"],
    mcetest_r      = test_metrics[, "mce"]
  )
}
