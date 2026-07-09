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


# Internal: eta^2_{Y | D, X}, the partial correlation ratio of the outcome.
# Equal to the R^2 of the restricted-regression residuals regressed on the
# cluster indicators (respecting model weights, if any).
.compute_eta2 <- function(model, cluster_vec) {
  .eta2_of(stats::residuals(model), cluster_vec, stats::weights(model))
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
#' Computes \eqn{\eta^2_{W \mid X}} for each benchmark covariate: the share of
#' its variance that lies between clusters. A benchmark for a cluster-level
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
#' @return A named numeric vector of \eqn{\eta^2_W}, one per benchmark covariate.
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
              "(eta^2_{W|X} = ", formatC(eta2_w[nm], digits = 2, format = "f"), "). ",
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
