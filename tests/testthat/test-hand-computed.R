################################################################################
### Hand-computed values on data small enough to check by hand (#12)
################################################################################

test_that("expit and logit are inverses of each other", {
  expect_equal(expit(0), 0.5)
  expect_equal(logit(0.5), 0)
  expect_equal(expit(logit(c(0.1, 0.25, 0.5, 0.9))), c(0.1, 0.25, 0.5, 0.9))

  # log(0.75 / 0.25) = log(3)
  expect_equal(logit(0.75), log(3))
  expect_equal(expit(log(3)), 0.75)

  expect_equal(expit(-Inf), 0)
  expect_equal(expit(Inf), 1)
})

test_that("auc_mannwhitney matches the rank-based definition", {
  # Perfect separation.
  expect_equal(auc_mannwhitney(c(0.1, 0.3, 0.6, 0.9), c(0, 0, 1, 1)), 1)
  # Perfectly reversed.
  expect_equal(auc_mannwhitney(c(0.1, 0.3, 0.6, 0.9), c(1, 1, 0, 0)), 0)
  # All scores tied: every comparison is a coin flip.
  expect_equal(auc_mannwhitney(rep(1, 4), c(0, 0, 1, 1)), 0.5)

  # One inversion out of four positive-negative pairs.
  #   positives 0.4, 0.9 ; negatives 0.2, 0.6
  #   pairs won: (0.4>0.2), (0.9>0.2), (0.9>0.6) = 3 of 4
  expect_equal(auc_mannwhitney(c(0.2, 0.4, 0.6, 0.9), c(0, 1, 0, 1)), 0.75)

  # A single tie counts as half a win: pairs are (1 vs 1) and (1 vs 3).
  expect_equal(auc_mannwhitney(c(1, 1, 3), c(0, 1, 1)), 0.75)

  expect_true(is.na(auc_mannwhitney(c(1, 2, 3), c(1, 1, 1))))
  expect_true(is.na(auc_mannwhitney(c(1, 2, 3), c(0, 0, 0))))
})

test_that("auc_mannwhitney agrees with a brute-force pair count", {
  set.seed(11L)
  score <- stats::rnorm(60)
  label <- stats::rbinom(60, 1, 0.4)

  positives <- score[label == 1L]
  negatives <- score[label == 0L]
  comparisons <- outer(positives, negatives, function(p, n) {
    (p > n) + 0.5 * (p == n)
  })

  expect_equal(auc_mannwhitney(score, label), mean(comparisons))
})

test_that("get_pooled_covariance matches the pooled formula", {
  # Two groups, one covariate. Group A: 1, 3 (deviations -1, +1).
  # Group B: 10, 14 (deviations -2, +2). Sum of squares = 2 + 8 = 10,
  # residual df = 4 - 2 = 2, so the pooled variance is exactly 5.
  data <- data.frame(x = c(1, 3, 10, 14),
                     T = c("A", "A", "B", "B"),
                     stringsAsFactors = FALSE)

  pooled <- get_pooled_covariance(data, X_vars = "x", treatment_var = "T")

  expect_equal(unname(pooled$S[1, 1]), 5)
  expect_equal(unname(pooled$S_inv[1, 1]), 1 / 5)
  expect_identical(dimnames(pooled$S), list("x", "x"))
})

test_that("get_pooled_covariance matches stats::cov on a single group", {
  # With one treatment group the pooled covariance reduces to the ordinary
  # sample covariance.
  set.seed(3L)
  data <- data.frame(x1 = stats::rnorm(30), x2 = stats::rnorm(30))
  data$T <- "A"

  pooled <- get_pooled_covariance(data, X_vars = c("x1", "x2"),
                                  treatment_var = "T")
  expect_equal(unname(pooled$S),
               unname(stats::cov(as.matrix(data[, c("x1", "x2")]))))
})

test_that("mahalanobis_distance_matrix matches stats::mahalanobis", {
  set.seed(5L)
  S_inv <- solve(matrix(c(4, 1, 1, 3), 2))
  X_query <- matrix(stats::rnorm(6), 3, 2)
  X_reference <- matrix(stats::rnorm(8), 4, 2)

  distances <- mahalanobis_distance_matrix(X_query, X_reference, S_inv)

  expect_equal(dim(distances), c(3L, 4L))
  for (i in 1:3) {
    expected <- sqrt(stats::mahalanobis(X_reference, X_query[i, ],
                                        S_inv, inverted = TRUE))
    expect_equal(distances[i, ], expected)
  }
})

test_that("mahalanobis distance reduces to Euclidean under the identity", {
  X_query <- matrix(c(0, 0, 3, 4), 2, 2, byrow = TRUE)
  X_reference <- matrix(c(0, 0, 0, 1), 2, 2, byrow = TRUE)

  distances <- mahalanobis_distance_matrix(X_query, X_reference, diag(2))

  expect_equal(distances[1, ], c(0, 1))
  expect_equal(distances[2, ], c(5, sqrt(9 + 9)))
})

test_that("build_group_distance_matrices lines up with its own pieces", {
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")
  anchor_rows <- which(labels == fixture$anchor)

  built <- build_group_distance_matrices(
    fixture$data, X_vars = fixture$X_vars, treatment_var = "treatment",
    anchor_rows = anchor_rows, groups = fixture$groups
  )

  pooled <- get_pooled_covariance(fixture$data, fixture$X_vars, "treatment")
  expect_equal(built$S_inv, pooled$S_inv)

  for (g in fixture$groups) {
    expect_equal(built$group_rows[[g]], which(labels == g))
    expect_equal(dim(built$D[[g]]),
                 c(length(anchor_rows), sum(labels == g)))
    expect_true(all(built$D[[g]] >= 0))
  }

  # Spot-check one entry against a direct computation.
  g <- fixture$groups[1]
  X <- covariate_matrix(fixture$data, fixture$X_vars)
  expected <- sqrt(stats::mahalanobis(
    X[built$group_rows[[g]][3], , drop = FALSE],
    X[anchor_rows[2], ], pooled$S_inv, inverted = TRUE
  ))
  expect_equal(built$D[[g]][2, 3], unname(expected))
})

test_that("compute_smd_balance matches the SMD formula by hand", {
  # Anchor values 1, 2, 3, 4 (mean 2.5, sample variance 5/3).
  # Group values 3, 4, 5, 6 (mean 4.5, sample variance 5/3).
  # Pooled sd = sqrt(5/3), SMD = (2.5 - 4.5) / sqrt(5/3).
  data <- data.frame(x = c(1, 2, 3, 4, 3, 4, 5, 6),
                     T = rep(c("A", "B"), each = 4),
                     stringsAsFactors = FALSE)
  matched <- data.frame(anchor = 1:4, B = 5:8)

  balance <- compute_smd_balance(data, matched, X_vars = "x", groups = "B")

  expect_equal(balance$by_covariate$smd, -2 / sqrt(5 / 3))
  expect_equal(balance$by_covariate$abs_smd, 2 / sqrt(5 / 3))
  expect_true(balance$by_covariate$smd_defined)
  expect_equal(balance$summary$mean_abs_smd, 2 / sqrt(5 / 3))
  expect_equal(balance$summary$n_undefined, 0L)
})

test_that("calc_caliper_3way matches its definition", {
  # 0.6 * sqrt(sum of the across-group mean variances / 3).
  set.seed(8L)
  ps <- matrix(stats::runif(60), 30, 2)
  groups <- rep(c("A", "B", "C"), each = 10)

  caliper <- calc_caliper_3way(ps, groups)

  variance_by_group <- sapply(1:2, function(j) {
    tapply(ps[, j], groups, stats::var)
  })
  expected <- 0.6 * sqrt(sum(rowMeans(variance_by_group)) / 3)

  expect_equal(caliper, expected)
  expect_gt(caliper, 0)
})

test_that("calc_caliper_3way rejects the wrong shape", {
  expect_error(calc_caliper_3way(matrix(1:9, 3, 3), rep("A", 3)))
  expect_error(calc_caliper_3way(matrix(1:6, 3, 2), rep("A", 2)))
})

test_that("compute_effective_sample_size is exact for uniform weights", {
  # With every weight equal, Kish's ESS equals the sample size.
  data <- data.frame(T = rep(c("A", "B"), times = c(10, 6)),
                     stringsAsFactors = FALSE)

  ess <- compute_effective_sample_size(data, rep(2, 16), treatment_var = "T")

  expect_equal(ess$ess, ess$n)
  expect_equal(ess$n, c(10, 6))
  expect_equal(ess$max_over_mean, c(1, 1))
})
