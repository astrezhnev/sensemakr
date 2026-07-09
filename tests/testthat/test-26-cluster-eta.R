context("cluster-adjusted sensitivity statistics")

# Unbalanced clustered design: treatment and covariate assigned at the cluster
# level, outcome has both a between-cluster and a within-cluster component.
make_cluster_data <- function(seed = 42, G = 50, beta = 0.4, within_sd = 1.5) {
  set.seed(seed)
  Ng <- sample(2:9, G, replace = TRUE)
  N  <- sum(Ng)
  cl <- factor(rep(seq_len(G), Ng))
  D  <- rbinom(G, 1, 0.5)[cl]
  X  <- rnorm(G)[cl]
  Y  <- beta * D + 0.6 * X + rnorm(G, sd = 1)[cl] + rnorm(N, sd = within_sd)
  data.frame(Y = Y, D = D, X = X, cl = cl)
}

# Cluster means, with the cluster sizes as regression weights. The weighted
# cluster-aggregated regression reproduces the unit-level treatment effect
# exactly under unbalanced cluster sizes.
aggregate_clusters <- function(dat) {
  agg <- aggregate(cbind(Y, D, X) ~ cl, data = dat, FUN = mean)
  agg$ng <- as.numeric(table(dat$cl))
  agg
}


test_that("eta2 equals the R2 of the restricted residuals on cluster indicators", {
  dat <- make_cluster_data()
  model <- lm(Y ~ D + X, data = dat)
  sens  <- sensemakr(model, treatment = "D", cluster = "cl")

  res <- residuals(model)
  # definition: between-cluster residual variance / total residual variance
  eta2_definition <- sum(ave(res, dat$cl)^2) / sum(res^2)
  # auxiliary regression: R2 of residuals on the cluster indicators
  eta2_auxiliary <- summary(lm(res ~ factor(dat$cl)))$r.squared

  expect_equal(sens$sensitivity_stats$eta2, eta2_definition)
  expect_equal(sens$sensitivity_stats$eta2, eta2_auxiliary)
  expect_true(sens$sensitivity_stats$eta2 > 0 && sens$sensitivity_stats$eta2 < 1)
})


test_that("cluster-adjusted RV and XRV match the weighted cluster-aggregated regression", {
  dat <- make_cluster_data()
  agg <- aggregate_clusters(dat)

  unit_model    <- lm(Y ~ D + X, data = dat)
  cluster_model <- lm(Y ~ D + X, data = agg, weights = ng)

  # the weighted cluster regression recovers the unit-level point estimate
  expect_equal(unname(coef(unit_model)["D"]), unname(coef(cluster_model)["D"]))

  unit_sens    <- sensemakr(unit_model, treatment = "D", cluster = "cl")
  cluster_sens <- sensemakr(cluster_model, treatment = "D")

  # robustness value: corrected unit-level == standard cluster-level
  expect_equal(as.numeric(unit_sens$sensitivity_stats$rv_q_cluster),
               as.numeric(cluster_sens$sensitivity_stats$rv_q))

  # extreme robustness value: corrected unit-level == standard cluster-level
  expect_equal(as.numeric(unit_sens$sensitivity_stats$xrv_q_cluster),
               as.numeric(extreme_robustness_value(cluster_model,
                                                   covariates = "D",
                                                   alpha = 1)))

  # the correction inflates the partial f, so the uncorrected unit-level RV is
  # strictly more conservative than the cluster-adjusted one
  expect_true(as.numeric(unit_sens$sensitivity_stats$rv_q) <
                as.numeric(unit_sens$sensitivity_stats$rv_q_cluster))
})


test_that("the cluster-adjusted RV lies on the equi-confounding line and zeroes the estimate", {
  dat   <- make_cluster_data()
  model <- lm(Y ~ D + X, data = dat)
  sens  <- sensemakr(model, treatment = "D", cluster = "cl")
  stats <- sens$sensitivity_stats

  rv   <- as.numeric(stats$rv_q_cluster)
  eta2 <- stats$eta2

  # a confounder equally strong at the cluster level sits at
  # (r2dz.x, r2yz.dx) = (rv, eta2 * rv) on the unit-level axes
  adjusted <- adjusted_estimate(estimate = stats$estimate,
                                se = stats$se,
                                dof = stats$dof,
                                r2dz.x = rv,
                                r2yz.dx = rv * eta2)
  expect_equal(adjusted, 0)
})


test_that("cluster argument errors when treatment varies within cluster", {
  dat <- make_cluster_data()
  set.seed(1)
  # individual-level treatment: varies within cluster
  dat$D_unit <- rbinom(nrow(dat), 1, 0.5)
  model <- lm(Y ~ D_unit + X, data = dat)

  expect_error(sensemakr(model, treatment = "D_unit", cluster = "cl"),
               "Treatment variable has a non-zero variance within cluster")
})


test_that("cluster = NULL reproduces the standard, unadjusted analysis", {
  dat   <- make_cluster_data()
  model <- lm(Y ~ D + X, data = dat)

  standard <- sensemakr(model, treatment = "D")
  clustered <- sensemakr(model, treatment = "D", cluster = "cl")

  # no cluster columns, and the unadjusted statistics are untouched
  expect_null(standard$cluster)
  expect_false("eta2" %in% names(standard$sensitivity_stats))
  expect_equal(as.numeric(standard$sensitivity_stats$rv_q),
               as.numeric(clustered$sensitivity_stats$rv_q))
  expect_equal(standard$sensitivity_stats$r2yd.x,
               clustered$sensitivity_stats$r2yd.x)
})


# Returns the plotting window (xlim, ylim) actually used by the contour plot,
# without leaving a device or file behind.
contour_window <- function(sens, ...) {
  path <- tempfile(fileext = ".png")
  png(path)
  on.exit({dev.off(); unlink(path)}, add = TRUE)
  plot(sens, ...)
  par("usr")
}

test_that("contour plot widens its limits only when the cluster RV falls off-canvas", {
  # a weak effect with substantial within-cluster noise keeps the
  # cluster-adjusted RV inside the default 0.4 window
  small <- sensemakr(lm(Y ~ D + X, data = make_cluster_data(seed = 1, G = 40,
                                                            beta = 0.15,
                                                            within_sd = 2.0)),
                     treatment = "D", cluster = "cl")
  expect_true(as.numeric(small$sensitivity_stats$rv_q_cluster) < 0.4)
  # the default window is left exactly as the unadjusted plot would have it
  expect_equal(contour_window(small), contour_window(small, lim = 0.4))

  # a large effect with little within-cluster noise pushes the cluster-adjusted
  # RV past the default window, which must then widen to contain the point
  big <- sensemakr(lm(Y ~ D + X, data = make_cluster_data(seed = 7, G = 40,
                                                          beta = 1.6,
                                                          within_sd = 0.25)),
                   treatment = "D", cluster = "cl")
  rv   <- as.numeric(big$sensitivity_stats$rv_q_cluster)
  eta2 <- big$sensitivity_stats$eta2
  expect_true(rv > 0.4)

  big_window <- contour_window(big)
  expect_true(rv        >= big_window[1] && rv        <= big_window[2])
  expect_true(rv * eta2 >= big_window[3] && rv * eta2 <= big_window[4])

  # an explicit lim is always respected, even when the point falls outside it
  user_window <- contour_window(big, lim = 0.45)
  expect_true(user_window[2] < 0.5)
})


test_that("print and summary report the cluster block only when clustering is enabled", {
  dat   <- make_cluster_data()
  model <- lm(Y ~ D + X, data = dat)

  clustered <- sensemakr(model, treatment = "D", cluster = "cl")
  standard  <- sensemakr(model, treatment = "D")

  clustered_print   <- paste(capture.output(print(clustered)), collapse = "\n")
  clustered_summary <- paste(capture.output(summary(clustered)), collapse = "\n")
  standard_print    <- paste(capture.output(print(standard)), collapse = "\n")
  standard_summary  <- paste(capture.output(summary(standard)), collapse = "\n")

  # the cluster block appears, names the cluster variable, and reports eta2,
  # the cluster-adjusted RV and the cluster-adjusted extreme RV
  expect_true(grepl("Cluster-Adjusted Statistics (cluster: cl)", clustered_print, fixed = TRUE))
  expect_true(grepl("Partial eta2 of outcome with cluster", clustered_print, fixed = TRUE))
  expect_true(grepl("Extreme Robustness Value", clustered_print, fixed = TRUE))
  expect_true(grepl(as.character(round(clustered$sensitivity_stats$eta2, 5)),
                    clustered_print, fixed = TRUE))

  # summary additionally carries the verbal interpretation
  expect_true(grepl("Verbal interpretation of cluster-adjusted sensitivity statistics",
                    clustered_summary, fixed = TRUE))
  expect_true(grepl("between clusters, not within them", clustered_summary, fixed = TRUE))

  # without a cluster argument, nothing cluster-related is printed
  expect_false(grepl("Cluster", standard_print, fixed = TRUE))
  expect_false(grepl("Cluster", standard_summary, fixed = TRUE))
  expect_false(grepl("eta2", standard_summary, fixed = TRUE))
})


test_that("cluster adjustment agrees between lm and fixest", {
  skip_if_not_installed("fixest")
  dat <- make_cluster_data()

  lm_sens <- sensemakr(lm(Y ~ D + X, data = dat),
                       treatment = "D", cluster = "cl")
  fe_sens <- sensemakr(fixest::feols(Y ~ D + X, data = dat),
                       treatment = "D", cluster = "cl")

  expect_equal(fe_sens$sensitivity_stats$eta2,
               lm_sens$sensitivity_stats$eta2)
  expect_equal(as.numeric(fe_sens$sensitivity_stats$rv_q_cluster),
               as.numeric(lm_sens$sensitivity_stats$rv_q_cluster))
  expect_equal(as.numeric(fe_sens$sensitivity_stats$xrv_q_cluster),
               as.numeric(lm_sens$sensitivity_stats$xrv_q_cluster))
})
