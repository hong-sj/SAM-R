################################################################################
### Consistent SMD orientation across the two balance functions (#10)
################################################################################

test_that("both balance tables report anchor minus comparator", {
  # A covariate the anchor group scores higher on must come out positive in
  # both tables. Previously compute_weighted_balance returned group minus
  # anchor, so the same imbalance had opposite signs in the two reports.
  data <- data.frame(
    x = c(rep(c(10, 12), 20), rep(c(0, 2), 20)),
    T = rep(c("A", "B"), each = 40),
    stringsAsFactors = FALSE
  )
  matched <- data.frame(anchor = 1:40, B = 41:80)

  matched_balance <- compute_smd_balance(data, matched, "x", "B")
  weighted_balance <- compute_weighted_balance(
    data, rep(1, nrow(data)), X_vars = "x",
    treatment_var = "T", anchor_level = "A"
  )

  expect_gt(matched_balance$by_covariate$smd, 0)
  expect_gt(weighted_balance$by_covariate$smd, 0)

  # The two differ only by the variance convention -- compute_smd_balance uses
  # the sample variance, the weighted version a weighted population variance,
  # which for uniform weights is the sqrt((n-1)/n) factor below. The point here
  # is the orientation, so that ratio is stated rather than assumed away.
  n_per_arm <- 40
  expect_equal(
    matched_balance$by_covariate$smd,
    weighted_balance$by_covariate$smd * sqrt((n_per_arm - 1) / n_per_arm)
  )
})

test_that("the two by_covariate tables share their column names", {
  fixture <- sam_fixture()

  matched_balance <- compute_smd_balance(fixture$data, fixture$match$matched,
                                         fixture$X_vars, fixture$groups)
  weights <- compute_balancing_weights(
    fixture$data, method = "overlap", X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )
  weighted_balance <- compute_weighted_balance(
    fixture$data, weights$weights, X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )

  expect_identical(names(matched_balance$by_covariate),
                   names(weighted_balance$by_covariate))
  expect_identical(names(matched_balance$summary),
                   names(weighted_balance$summary))
})

test_that("abs_smd agrees with smd in both tables", {
  fixture <- sam_fixture()

  matched_balance <- compute_smd_balance(fixture$data, fixture$match$matched,
                                         fixture$X_vars, fixture$groups)
  expect_equal(matched_balance$by_covariate$abs_smd,
               abs(matched_balance$by_covariate$smd))

  weights <- compute_balancing_weights(
    fixture$data, method = "iptw", X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )
  weighted_balance <- compute_weighted_balance(
    fixture$data, weights$weights, X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )
  expect_equal(weighted_balance$by_covariate$abs_smd,
               abs(weighted_balance$by_covariate$smd))
})

test_that("unweighted weighted-balance reproduces the raw-cohort SMD", {
  # With uniform weights the weighted table must reduce to the plain
  # standardized mean difference of the full cohort, which pins the
  # orientation independently of the matched-cohort code path.
  data <- sam_small_data(n = 120L, n_groups = 2L)

  balance <- compute_weighted_balance(
    data, rep(1, nrow(data)), X_vars = "x1",
    treatment_var = "T", anchor_level = "A"
  )

  is_anchor <- data$T == "A"
  population_variance <- function(x) sum((x - mean(x))^2) / length(x)
  expected <- (mean(data$x1[is_anchor]) - mean(data$x1[!is_anchor])) /
    sqrt((population_variance(data$x1[is_anchor]) +
            population_variance(data$x1[!is_anchor])) / 2)

  expect_equal(balance$by_covariate$smd, expected)
})
