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


# A design whose covariates and benchmark are all cluster-constant, plus one
# benchmark that varies within cluster. The weighted cluster-aggregated
# regression is then the exact counterpart of the unit-level regression.
make_benchmark_data <- function(seed = 11, G = 60) {
  set.seed(seed)
  Ng <- sample(3:10, G, replace = TRUE)
  N  <- sum(Ng)
  cl <- factor(rep(seq_len(G), Ng))
  D  <- rbinom(G, 1, 0.5)[cl]
  Wg <- rnorm(G)[cl]                       # cluster-constant benchmark
  Xg <- rnorm(G)[cl]                       # cluster-constant covariate
  Wu <- Wg + rnorm(N, sd = 1.2)            # benchmark varying within cluster
  Y  <- 0.4 * D + 0.5 * Wg + 0.3 * Xg + rnorm(G, sd = 1)[cl] + rnorm(N, sd = 1.2)
  data.frame(Y = Y, D = D, Wg = Wg, Xg = Xg, Wu = Wu, cl = cl)
}


test_that("bounds gain the cluster-scale outcome bound, matching the cluster regression", {
  dat   <- make_benchmark_data()
  model <- lm(Y ~ D + Wg + Xg, data = dat)
  sens  <- sensemakr(model, treatment = "D", benchmark_covariates = "Wg",
                     kd = 1:3, cluster = "cl")
  eta2  <- sens$cluster$eta2

  # exact conversion, not an approximation
  expect_equal(sens$bounds$r2yz.dx_cluster, sens$bounds$r2yz.dx / eta2)
  expect_true(all(sens$bounds$feasible))
  # a cluster-constant benchmark has eta2_W == 1
  expect_equal(sens$bounds$eta2_benchmark, rep(1, 3))

  # and it equals the bound from the size-weighted cluster-aggregated regression
  agg <- aggregate(cbind(Y, D, Wg, Xg) ~ cl, data = dat, FUN = mean)
  agg$ng <- as.numeric(table(dat$cl))
  cluster_bounds <- ovb_bounds(lm(Y ~ D + Wg + Xg, data = agg, weights = ng),
                               treatment = "D", benchmark_covariates = "Wg", kd = 1:3)
  expect_equal(sens$bounds$r2yz.dx_cluster, cluster_bounds$r2yz.dx)
  # the treatment-side bound is level-invariant for a cluster-constant benchmark
  expect_equal(sens$bounds$r2dz.x, cluster_bounds$r2dz.x)
})


test_that("an infeasible cluster-level bound is flagged and warned about", {
  dat   <- make_benchmark_data()
  model <- lm(Y ~ D + Wg + Xg, data = dat)
  expect_warning(sens <- sensemakr(model, treatment = "D",
                                   benchmark_covariates = "Wg",
                                   kd = 8, cluster = "cl"),
                 "r2yz.dx_cluster greater than 1")
  expect_false(sens$bounds$feasible)
  expect_true(sens$bounds$r2yz.dx_cluster > 1)
})


test_that("a benchmark covariate varying within cluster triggers a warning", {
  dat   <- make_benchmark_data()
  model <- lm(Y ~ D + Wu + Xg, data = dat)

  expect_warning(sensemakr(model, treatment = "D", benchmark_covariates = "Wu",
                           kd = 1, cluster = "cl"),
                 "varies within clusters")

  # a cluster-constant benchmark is silent
  clean <- lm(Y ~ D + Wg + Xg, data = dat)
  expect_silent(sensemakr(clean, treatment = "D", benchmark_covariates = "Wg",
                          kd = 1, cluster = "cl"))
})


test_that("manual bounds and group benchmarks coexist with the cluster columns", {
  dat   <- make_benchmark_data()
  model <- lm(Y ~ D + Wg + Xg, data = dat)

  # manual bound plus a benchmark: the rows must bind together
  both <- sensemakr(model, treatment = "D", benchmark_covariates = "Wg", kd = 1,
                    r2dz.x = 0.1, r2yz.dx = 0.1, cluster = "cl")
  expect_equal(nrow(both$bounds), 2)
  expect_true(is.na(both$bounds$eta2_benchmark[1]))   # manual bound
  expect_equal(both$bounds$eta2_benchmark[2], 1)      # benchmark

  # group benchmarks take the cluster columns too; eta2_benchmark is NA for groups
  grouped <- sensemakr(model, treatment = "D",
                       benchmark_covariates = list(grp = c("Wg", "Xg")),
                       kd = 1, cluster = "cl")
  expect_true("r2yz.dx_cluster" %in% names(grouped$bounds))
  expect_true(is.na(grouped$bounds$eta2_benchmark))
})


test_that("bounds are untouched when no cluster is supplied", {
  dat   <- make_benchmark_data()
  model <- lm(Y ~ D + Wg + Xg, data = dat)
  sens  <- sensemakr(model, treatment = "D", benchmark_covariates = "Wg", kd = 1)
  expect_false("r2yz.dx_cluster" %in% names(sens$bounds))
  expect_false("feasible" %in% names(sens$bounds))
  expect_false("eta2_benchmark" %in% names(sens$bounds))
})


test_that("$cluster carries eta2, n_clusters, the column name and the cluster sizes", {
  dat   <- make_benchmark_data()
  sens  <- sensemakr(lm(Y ~ D + Wg + Xg, data = dat), treatment = "D", cluster = "cl")

  expect_equal(sens$cluster$n_clusters, nlevels(dat$cl))
  expect_equal(sens$cluster$cluster, "cl")
  expect_equal(sum(sens$cluster$cluster_sizes), nrow(dat))
  expect_equal(as.integer(sens$cluster$cluster_sizes), as.integer(table(dat$cl)))
})


test_that("a numeric cluster column works the same as a factor", {
  dat <- make_benchmark_data()
  dat$cl_num <- as.integer(dat$cl)   # e.g. an integer FIPS code

  by_factor  <- sensemakr(lm(Y ~ D + Wg + Xg, data = dat), treatment = "D", cluster = "cl")
  by_numeric <- sensemakr(lm(Y ~ D + Wg + Xg, data = dat), treatment = "D", cluster = "cl_num")

  expect_equal(by_numeric$cluster$eta2, by_factor$cluster$eta2)
  expect_equal(by_numeric$cluster$n_clusters, by_factor$cluster$n_clusters)
  expect_equal(as.numeric(by_numeric$sensitivity_stats$rv_q_cluster),
               as.numeric(by_factor$sensitivity_stats$rv_q_cluster))
})


test_that("ovb_minimal_reporting reports the cluster row only when clustering is on", {
  dat   <- make_benchmark_data()
  model <- lm(Y ~ D + Wg + Xg, data = dat)
  clustered <- sensemakr(model, treatment = "D", benchmark_covariates = "Wg",
                         kd = 1, cluster = "cl")
  standard  <- sensemakr(model, treatment = "D", benchmark_covariates = "Wg", kd = 1)

  for (fmt in c("latex", "html", "pure_html")) {
    tab_c <- ovb_minimal_reporting(clustered, format = fmt, verbose = FALSE)
    tab_s <- ovb_minimal_reporting(standard,  format = fmt, verbose = FALSE)
    expect_true(grepl("Cluster-adjusted", tab_c, fixed = TRUE))
    expect_true(grepl("60 clusters", tab_c, fixed = TRUE))
    expect_false(grepl("Cluster", tab_s, fixed = TRUE))
  }
})


test_that("extreme plot scales the outcome scenarios to the cluster ceiling", {
  dat   <- make_benchmark_data()
  sens  <- sensemakr(lm(Y ~ D + Wg + Xg, data = dat), treatment = "D", cluster = "cl")
  eta2  <- sens$cluster$eta2
  xrv   <- as.numeric(sens$sensitivity_stats$xrv_q_cluster)

  path <- tempfile(fileext = ".png")
  png(path); on.exit({dev.off(); unlink(path)}, add = TRUE)

  # scaled: the returned scenarios are fractions of the ceiling
  scaled <- ovb_extreme_plot(sens, r2yz.dx = c(1, 0.75, 0.5))
  expect_equal(scaled[[1]]$r2yz.dx[1], eta2)

  # unscaled: the old behaviour remains reachable
  unscaled <- ovb_extreme_plot(sens, r2yz.dx = c(1, 0.75, 0.5), scale_to_ceiling = FALSE)
  expect_equal(unscaled[[1]]$r2yz.dx[1], 1)

  # the 100%-of-ceiling curve crosses the threshold exactly at the cluster XRV
  at_xrv <- adjusted_estimate(estimate = sens$sensitivity_stats$estimate,
                              se = sens$sensitivity_stats$se,
                              dof = sens$sensitivity_stats$dof,
                              r2dz.x = xrv, r2yz.dx = eta2)
  expect_equal(at_xrv, 0)
})


test_that("extreme plot legend labels are rounded and overridable", {
  path <- tempfile(fileext = ".png")
  png(path); on.exit({dev.off(); unlink(path)}, add = TRUE)

  # a non-round r2yz.dx must not produce a full-precision legend; the plot
  # simply has to render, and custom labels must be accepted
  expect_silent(ovb_extreme_plot(estimate = 0.1, se = 0.05, dof = 100,
                                 r2dz.x = 0.1, r2yz.dx = c(0.318739766639781, 0.2)))
  expect_silent(ovb_extreme_plot(estimate = 0.1, se = 0.05, dof = 100,
                                 r2dz.x = 0.1, r2yz.dx = c(0.5, 0.25),
                                 legend.labels = c("full ceiling", "half ceiling")))
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


test_that("cluster bounds and the benchmark guard work for fixest as for lm", {
  skip_if_not_installed("fixest")
  dat <- make_benchmark_data()

  lm_sens <- sensemakr(lm(Y ~ D + Wg + Xg, data = dat), treatment = "D",
                       benchmark_covariates = "Wg", kd = 1:2, cluster = "cl")
  fe_sens <- sensemakr(fixest::feols(Y ~ D + Wg + Xg, data = dat), treatment = "D",
                       benchmark_covariates = "Wg", kd = 1:2, cluster = "cl")

  expect_equal(fe_sens$bounds$r2yz.dx_cluster, lm_sens$bounds$r2yz.dx_cluster)
  expect_equal(fe_sens$bounds$eta2_benchmark, lm_sens$bounds$eta2_benchmark)

  # the within-cluster benchmark guard fires for fixest too
  expect_warning(sensemakr(fixest::feols(Y ~ D + Wu + Xg, data = dat),
                           treatment = "D", benchmark_covariates = "Wu",
                           kd = 1, cluster = "cl"),
                 "varies within clusters")
})
