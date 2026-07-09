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


# Internal: eta^2_{Y | D, X}, the partial correlation ratio of the outcome.
# Equal to the R^2 of the restricted-regression residuals regressed on the
# cluster indicators (respecting model weights, if any).
.compute_eta2 <- function(model, cluster_vec) {
  res <- stats::residuals(model)
  w   <- stats::weights(model)
  aux <- if (is.null(w)) {
    stats::lm(res ~ factor(cluster_vec))
  } else {
    stats::lm(res ~ factor(cluster_vec), weights = w)
  }
  summary(aux)$r.squared
}
