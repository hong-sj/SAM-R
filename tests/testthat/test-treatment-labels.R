################################################################################
### Consistent treatment-label handling (#1)
################################################################################

numeric_treatment_data <- function() {
  data <- sam_small_data(n = 150L, n_groups = 3L)
  # 0/1/2 rather than A/B/C: the level names are now valid `levels()` positions,
  # which is what made relevel() silently pick the wrong reference category.
  data$T <- as.numeric(factor(data$T)) - 1
  data
}

test_that("treatment_labels canonicalises every column type", {
  data <- data.frame(chr = c("a", "b"), num = c(0, 1),
                     fct = factor(c("a", "b"), levels = c("b", "a")))

  expect_identical(treatment_labels(data, "chr"), c("a", "b"))
  expect_identical(treatment_labels(data, "num"), c("0", "1"))
  expect_identical(treatment_labels(data, "fct"), c("a", "b"))

  expect_error(treatment_labels(data, "nope"), "Treatment column not found")
})

test_that("treatment_level rejects what cannot name a single level", {
  expect_identical(treatment_level(0), "0")
  expect_identical(treatment_level(factor("a")), "a")
  expect_null(treatment_level(NULL))

  expect_error(treatment_level(c(1, 2)), "single value")
  expect_error(treatment_level(NA), "must not be NA")
})

test_that("a numeric anchor_level names a level, not a position", {
  data <- numeric_treatment_data()
  X_vars <- c("x1", "x2", "x3")

  # Previously: relevel() read `ref = 0` as an index into levels() and raised
  # "'ref' must be an existing level" for a level that is plainly present.
  fit_numeric <- estimate_gps_multinom(data, X_vars = X_vars,
                                       treatment_var = "T", anchor_level = 0)
  expect_identical(colnames(fit_numeric$gps)[1], "0")

  # Level "2" sits at position 3, so a positional reading would have selected a
  # different reference category without saying so.
  fit_two <- estimate_gps_multinom(data, X_vars = X_vars,
                                   treatment_var = "T", anchor_level = 2)
  fit_two_chr <- estimate_gps_multinom(data, X_vars = X_vars,
                                       treatment_var = "T", anchor_level = "2")
  expect_identical(colnames(fit_two$gps)[1], "2")
  expect_equal(fit_two$gps, fit_two_chr$gps)
})

test_that("a numeric treatment column runs end to end", {
  data <- numeric_treatment_data()
  X_vars <- c("x1", "x2", "x3")

  fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = 0)
  search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                 anchor_level = 0, top_m = 5)

  # Previously: anchor_rows was empty and the whole pipeline returned NaN.
  expect_equal(length(search$anchor_rows), sum(data$T == 0))
  expect_setequal(search$groups, c("1", "2"))

  result <- sam_match(data, search, X_vars = X_vars, treatment_var = "T")
  expect_gt(result$matching_rate, 0)
  expect_true(all(treatment_labels(data, "T")[result$matched$anchor] == "0"))
  for (g in search$groups) {
    expect_true(all(treatment_labels(data, "T")[result$matched[[g]]] == g))
  }
})

test_that("an anchor level that selects no rows raises", {
  fixture <- sam_fixture()

  # Previously: anchor_rows of length 0, `groups` still holding every level,
  # an empty matched frame and a NaN matching rate -- with no warning.
  expect_error(
    gps_candidate_search(fixture$data, fixture$gps, treatment_var = "treatment",
                         anchor_level = "ZZZ", gps_space = "logit"),
    "not a column of `gps`"
  )
  expect_error(
    estimate_gps_multinom(fixture$data, X_vars = fixture$X_vars,
                          treatment_var = "treatment", anchor_level = "ZZZ"),
    "No rows found"
  )
})

test_that("the error names the levels that are actually present", {
  data <- sam_small_data(n = 60L, n_groups = 2L)
  expect_error(
    estimate_gps_multinom(data, X_vars = c("x1", "x2"), treatment_var = "T",
                          anchor_level = "Z"),
    "Levels present"
  )
})

test_that("a factor treatment keeps its own level order", {
  data <- sam_small_data(n = 150L, n_groups = 3L)
  data$T <- factor(data$T, levels = c("C", "B", "A"))
  X_vars <- c("x1", "x2", "x3")

  fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "B")
  expect_identical(colnames(fit$gps)[1], "B")
  expect_setequal(colnames(fit$gps), c("A", "B", "C"))

  search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                 anchor_level = "B", top_m = 5)
  expect_equal(length(search$anchor_rows), sum(data$T == "B"))
  expect_silent(sam_match(data, search, X_vars = X_vars, treatment_var = "T"))
})

test_that("the weighting module normalises anchor_level the same way", {
  data <- numeric_treatment_data()
  X_vars <- c("x1", "x2", "x3")

  numeric_anchor <- evaluate_comparator_weighting(
    data, method = "overlap", X_vars = X_vars,
    treatment_var = "T", anchor_level = 0
  )
  character_anchor <- evaluate_comparator_weighting(
    data, method = "overlap", X_vars = X_vars,
    treatment_var = "T", anchor_level = "0"
  )
  expect_equal(numeric_anchor$balance, character_anchor$balance)

  expect_error(
    compute_weighted_balance(data, rep(1, nrow(data)), X_vars = X_vars,
                             treatment_var = "T", anchor_level = "Z"),
    "No rows found"
  )
})
