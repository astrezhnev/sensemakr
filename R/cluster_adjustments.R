# Cluster-adjusted sensitivity statistics (clustered treatment assignment)
# ------------------------------------------------------------------------
# When treatment is assigned at the cluster level, the outcome-side partial
# R^2 in an individual-level regression is contaminated by within-cluster
# variation that a cluster-level confounder cannot explain. The correction
# factor is Pearson's partial correlation ratio
#
#     eta^2_{Y | D, X} = between-cluster residual variance / total residual variance,
#
# which is equal to the R^2 from a regression of the restricted-regression
# residuals on a complete set of cluster indicators. Dividing the
# outcome-treatment partial f by eta recovers the cluster-level robustness
# value and extreme robustness value from the unit-level regression, without
# having to fit the cluster-aggregated regression. See Strezhnev (2026),
# "Omitted variable bias sensitivity analysis with clustered treatment
# assignment."


# Internal: resolve a cluster column name to a vector aligned with model rows.
# Works for both lm and fixest (feols): lm keeps the model frame, feols does
# not, so we fall back to the original data and drop the unused rows (NA rows
# tracked in `na.action` for lm and `obs_selection$obsRemoved` for fixest).
.get_cluster <- function(model, cluster) {
  if (!(is.character(cluster) && length(cluster) == 1)) {
    stop("`cluster` must be a single column name (character).")
  }
  n <- length(stats::residuals(model))
  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (!is.null(mf) && cluster %in% names(mf)) return(mf[[cluster]])
  # feols records the calling environment in `call_env`; lm carries it on the
  # formula. Try both so the data argument resolves for either model type.
  envs <- list(model$call_env, environment(stats::formula(model)))
  data <- NULL
  for (env in envs) {
    if (!is.environment(env)) next
    data <- tryCatch(eval(stats::getCall(model)$data, envir = env),
                     error = function(e) NULL)
    if (is.data.frame(data) && cluster %in% names(data)) break
    data <- NULL
  }
  if (is.null(data) || !cluster %in% names(data)) {
    stop("Cluster variable '", cluster, "' not found in the model frame ",
         "or in the model's data.")
  }
  vec <- data[[cluster]]
  if (length(vec) != n) {
    if (!is.null(model$na.action)) {
      vec <- vec[-as.integer(model$na.action)]           # lm: positive indices
    } else if (!is.null(model$obs_selection$obsRemoved)) {
      vec <- vec[model$obs_selection$obsRemoved]          # fixest: negative indices
    }
  }
  if (length(vec) != n) {
    stop("Could not align cluster variable '", cluster, "' to the model's ",
         "estimation sample (length ", length(vec), " vs ", n, ").")
  }
  vec
}


# Internal: error unless treatment is constant within every cluster.
.check_cluster_treatment <- function(model, treatment, cluster_vec) {
  mm <- stats::model.matrix(model)
  if (!treatment %in% colnames(mm)) {
    stop("Treatment '", treatment, "' not found in the model design matrix.")
  }
  d <- mm[, treatment]
  ss_total <- sum((d - mean(d))^2)
  ss_within <- sum((d - stats::ave(d, cluster_vec))^2)
  if (ss_total <= 0 || ss_within > 1e-10 * ss_total) {
    stop("Treatment variable has a non-zero variance within cluster")
  }
  invisible(TRUE)
}


# Internal: between-cluster share of the variance of `v`, i.e. the R^2 of a
# regression of `v` on the cluster indicators. Equals 1 - SS_within/SS_total.
.eta2_of <- function(v, cluster_vec, w = NULL) {
  aux <- if (is.null(w)) {
    stats::lm(v ~ factor(cluster_vec))
  } else {
    stats::lm(v ~ factor(cluster_vec), weights = w)
  }
  summary(aux)$r.squared
}


# Internal: cluster means of `v`, weighted by `w` when the model carries weights.
.cluster_mean <- function(v, cluster_vec, w = NULL) {
  if (is.null(w)) return(stats::ave(v, cluster_vec))
  wsum <- rowsum(w, cluster_vec)
  idx  <- match(as.character(cluster_vec), rownames(wsum))
  (rowsum(w * v, cluster_vec) / wsum)[idx]
}


# Internal: warn when the model constrains a covariate's within-cluster and
# between-cluster slopes to be equal. A covariate is fine if it is constant within
# clusters, or if its cluster mean lies in the column span of the design (so the
# two slopes are free to differ). Only a covariate that varies within clusters
# *and* whose cluster mean is absent from the model triggers the restriction.
#
# The cluster mean is the *weighted* one when the model has weights: a weighted fit
# orthogonalizes against the weighted cluster space, so it is the weighted
# decomposition whose slopes are free to differ. Using the unweighted mean here
# would warn on a correctly weighted Mundlak specification and stay silent on an
# incorrectly weighted one. This matches the condition under which .compute_eta2()
# departs from the between-cluster share of the residual variance.
#
# The warning deliberately avoids the phrase "varies within clusters" used by
# .check_cluster_benchmarks(), so that callers muffling one do not swallow the other.
.check_cluster_model <- function(model, cluster_vec, tol = 1e-10) {
  mm  <- stats::model.matrix(model)
  nms <- setdiff(colnames(mm), "(Intercept)")
  if (!length(nms)) return(invisible(TRUE))
  qrX <- qr(mm)
  w   <- stats::weights(model)

  bad <- Filter(function(nm) {
    v      <- mm[, nm]
    ss_tot <- sum((v - mean(v))^2)
    m      <- .cluster_mean(v, cluster_vec, w)              # its cluster mean
    ss_w   <- sum((v - m)^2)
    if (ss_tot <= 0 || ss_w <= tol * ss_tot) return(FALSE)  # cluster-constant: fine
    ss_m   <- sum((m - mean(m))^2)
    # A covariate with no between-cluster component (e.g. a within-cluster
    # deviation regressor) constrains nothing: its cluster mean is the zero
    # vector, trivially spanned. Compare ss_m against the covariate's own scale,
    # not against zero, or floating-point noise in a numerically-zero cluster
    # mean gets compared against noise in its own projection.
    if (ss_m <= tol * ss_tot) return(FALSE)
    sum(qr.resid(qrX, m)^2) > tol * ss_m                    # mean not spanned
  }, nms)

  if (length(bad)) {
    warning("Covariate(s) ", paste0("'", bad, "'", collapse = ", "),
            " constrain their between-cluster and within-cluster slopes to be ",
            "equal: each one differs across units inside a cluster, while its ",
            "cluster mean is absent from the model. The treatment coefficient ",
            "therefore differs from the unrestricted estimate, and r2yz.dx is no ",
            "longer just the share of between-cluster outcome variance a ",
            "cluster-level confounder explains: it also includes the confounder's ",
            "influence on the fitted coefficient of these covariates. Add their ",
            "cluster means to the model.",
            call. = FALSE)
  }
  invisible(TRUE)
}


# Internal: eta^2_{Y | D, X}, the largest outcome-side partial R^2 that a
# cluster-level confounder can attain given the model's regressors. This is the
# partial R^2 of the outcome with the cluster indicators:
#
#     eta2 = 1 - SSR(regressors + cluster dummies) / SSR(regressors)
#
# By Frisch-Waugh-Lovell the numerator is the SSR of the within-cluster demeaned
# regression, so the G cluster dummies are never formed.
#
# This coincides with the between-cluster share of the residual variance --
# .eta2_of(residuals(model), cluster_vec) -- exactly when every fitted value is
# constant within clusters. It does not when a covariate varies within clusters
# and its cluster mean is absent from the model: that specification ties the
# covariate's within- and between-cluster slopes together, letting a cluster-level
# confounder reduce within-cluster residual variance by shifting the constrained
# coefficient. The between-cluster share then falls below the true ceiling, which
# would make the cluster-adjusted robustness values anti-conservative. See
# .check_cluster_model(), which warns on exactly that specification.
.compute_eta2 <- function(model, cluster_vec) {
  # model.frame() is empty for fixest, so recover the response this way
  y  <- stats::fitted(model) + stats::residuals(model)
  mm <- stats::model.matrix(model)
  w  <- stats::weights(model)
  if (is.null(w)) w <- rep(1, length(y))
  sw <- sqrt(w)

  # weighted within-cluster demeaning
  dm <- function(v) v - .cluster_mean(v, cluster_vec, w)

  yt <- dm(y)
  Xt <- apply(mm, 2, dm)
  # regressors constant within clusters demean to zero (the intercept among them)
  keep <- colSums((sw * Xt)^2) > 1e-10 * pmax(colSums((sw * mm)^2), 1)

  ssr_full <- if (any(keep)) {
    sum(qr.resid(qr(sw * Xt[, keep, drop = FALSE]), sw * yt)^2)
  } else {
    sum(w * yt^2)                     # all regressors cluster-constant
  }
  1 - ssr_full / sum(w * stats::residuals(model)^2)
}


# Internal: flatten `benchmark_covariates` (a character vector, or a named list
# of character vectors for group benchmarks) to the covariate names it refers to.
.benchmark_names <- function(benchmark_covariates) {
  if (is.list(benchmark_covariates)) unlist(benchmark_covariates, use.names = FALSE)
  else benchmark_covariates
}


#' Between-cluster share of the variance of benchmark covariates
#'
#' @description
#' Computes, for each benchmark covariate, the share of its variance that lies
#' between clusters. This is an \emph{unconditional} quantity: the covariate is not
#' residualized on the other regressors, so it is not \eqn{\eta^2_{W \mid X}}.
#' It equals 1 for a cluster-level benchmark. A benchmark for a cluster-level
#' confounder must itself be cluster-constant (\eqn{\eta^2_W = 1}). When it is
#' not, its within-cluster component is orthogonal to a cluster-assigned
#' treatment: it inflates the outcome-side bound while contributing nothing to
#' the treatment-side bound, so the resulting bounds do not describe a feasible
#' cluster-level confounder.
#'
#' Warns, rather than errors, because deliberately contrasting a unit-level with
#' a cluster-level benchmark is a legitimate analysis.
#'
#' @param model an \code{lm} or \code{fixest} model.
#' @param benchmark_covariates benchmark covariates, as passed to
#'   \code{\link{sensemakr}}: a character vector, or a named list of character
#'   vectors for group benchmarks.
#' @param cluster_vec the cluster variable, aligned with the model's rows.
#' @param tol tolerance on the between-cluster share. Default \code{1e-8}.
#' @return A named numeric vector giving the between-cluster share of each benchmark
#'   covariate's variance.
#' @keywords internal
#' @noRd
.check_cluster_benchmarks <- function(model, benchmark_covariates, cluster_vec,
                                      tol = 1e-8) {
  bench_names <- .benchmark_names(benchmark_covariates)
  mm <- stats::model.matrix(model)
  w  <- stats::weights(model)

  eta2_w <- vapply(bench_names, function(nm) {
    if (!nm %in% colnames(mm)) return(NA_real_)
    .eta2_of(mm[, nm], cluster_vec, w)
  }, numeric(1))

  for (nm in names(eta2_w)) {
    if (!is.na(eta2_w[nm]) && eta2_w[nm] < 1 - tol) {
      warning("Benchmark covariate '", nm, "' varies within clusters ",
              "(between-cluster share of its variance = ",
              formatC(eta2_w[nm], digits = 2, format = "f"), "). ",
              "Only its cluster mean can serve as a benchmark for a cluster-level ",
              "confounder; the within-cluster component is orthogonal to treatment ",
              "and inflates the outcome-side bound. Consider entering the covariate ",
              "in Mundlak form and benchmarking against its cluster mean.",
              call. = FALSE)
    }
  }
  eta2_w
}


# Internal: recover the benchmark name from a bound label ("1x female" -> "female",
# "1/2x female" -> "female"). Group benchmarks are labelled by the list name, which
# is not a covariate, so the lookup they feed yields NA.
.bound_label_benchmark <- function(bound_label) {
  sub("^[^ ]+x ", "", bound_label)
}


# Internal: put the outcome-side bounds on the cluster scale. The unit-level
# r2yz.dx is the cluster-level one deflated by eta2 (the benchmark's partial f2
# with the outcome deflates by exactly eta2, while the treatment-side multiplier
# is level-invariant for a cluster-constant benchmark), so the conversion is a
# division. A value above 1 means the benchmarked confounder is impossible at the
# cluster level: it would have to explain more than all the between-cluster
# residual variance of the outcome.
.add_cluster_bounds <- function(bounds, eta2) {
  if (is.null(bounds) || is.null(eta2)) return(bounds)
  bounds$r2yz.dx_cluster <- bounds$r2yz.dx / eta2
  bounds$feasible <- bounds$r2yz.dx_cluster <= 1
  if (any(!bounds$feasible)) {
    warning("Implied bound on r2yz.dx_cluster greater than 1 for ",
            sum(!bounds$feasible), " of ", nrow(bounds), " bound(s). ",
            "A cluster-level confounder cannot explain more than 100% of the ",
            "between-cluster residual variance of the outcome, so the benchmarked ",
            "confounder is infeasible. Try a lower kd and/or ky.",
            call. = FALSE)
  }
  bounds
}
