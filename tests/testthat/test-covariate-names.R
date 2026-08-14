################################################################################
### Covariate names that are not syntactic R names (#11)
################################################################################

test_that("a covariate name with a space survives GPS fitting", {
  data <- sam_small_data(n = 120L, n_groups = 3L)
  names(data)[names(data) == "x1"] <- "age yrs"
  names(data)[names(data) == "x2"] <- "sbp-1"
  names(data)[names(data) == "x3"] <- "1st lactate"
  X_vars <- c("age yrs", "sbp-1", "1st lactate")

  # Previously: data.frame() mangled these to age.yrs / sbp.1 / X1st.lactate
  # while the formula still asked for the originals, so multinom() was handed a
  # formula referencing columns that were not in the frame.
  fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "A")

  expect_equal(nrow(fit$gps), nrow(data))
  expect_setequal(colnames(fit$gps), c("A", "B", "C"))
  expect_equal(unname(rowSums(fit$gps)), rep(1, nrow(data)), tolerance = 1e-8)
})

test_that("the rest of the pipeline handles the same names", {
  data <- sam_small_data(n = 120L, n_groups = 3L)
  names(data)[names(data) == "x1"] <- "age yrs"
  X_vars <- c("age yrs", "x2", "x3")

  fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "A")
  search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                 anchor_level = "A", top_m = 5)
  result <- sam_match(data, search, X_vars = X_vars, treatment_var = "T")

  expect_gt(result$matching_rate, 0)

  balance <- compute_smd_balance(data, result$matched, X_vars, search$groups)
  expect_true("age yrs" %in% balance$by_covariate$covariate)
  expect_true(all(is.finite(balance$summary$max_abs_smd)))
})

test_that("a covariate colliding with the response column is rejected", {
  data <- sam_small_data(n = 60L, n_groups = 2L)
  data$treatment_factor <- data$x1

  expect_error(
    estimate_gps_multinom(data, X_vars = c("x1", "treatment_factor"),
                          treatment_var = "T", anchor_level = "A"),
    "collides with the internal response column"
  )
})

test_that("quoting the formula leaves ordinary names untouched", {
  data <- sam_small_data(n = 120L, n_groups = 3L)
  X_vars <- c("x1", "x2", "x3")

  renamed <- data
  names(renamed)[names(renamed) == "x1"] <- "a b"

  plain <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                                 anchor_level = "A")
  quoted <- estimate_gps_multinom(renamed, X_vars = c("a b", "x2", "x3"),
                                  treatment_var = "T", anchor_level = "A")

  # Renaming a column must not change the fit.
  expect_equal(unname(plain$gps), unname(quoted$gps))
})
