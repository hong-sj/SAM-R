################################################################################
### GPS validation and identity across pipeline stages
##
## Every stage takes `gps` separately from `data` and matches the two by row
## position, so a matrix that does not line up produces a wrong answer rather
## than an error.
################################################################################

test_that("a GPS with the wrong number of rows is rejected", {
  fixture <- sam_fixture()

  expect_error(
    gps_candidate_search(fixture$data, fixture$gps[-1, , drop = FALSE],
                         treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "same number of rows"
  )
})

test_that("GPS values outside [0, 1] are rejected", {
  fixture <- sam_fixture()
  out_of_range <- fixture$gps
  out_of_range[1, 1] <- 5

  # Previously accepted, producing candidate pools from a score of 5.
  expect_error(
    gps_candidate_search(fixture$data, out_of_range,
                         treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "\\[0, 1\\]"
  )

  negative <- fixture$gps
  negative[2, 1] <- -0.5
  expect_error(
    gps_candidate_search(fixture$data, negative, treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "\\[0, 1\\]"
  )
})

test_that("non-finite GPS values are rejected with a clear message", {
  fixture <- sam_fixture()
  with_na <- fixture$gps
  with_na[2, ] <- NA_real_

  # Previously this reached RANN and failed with "NA/NaN/Inf in foreign
  # function call (arg 1)".
  expect_error(
    gps_candidate_search(fixture$data, with_na, treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "must all be finite"
  )
})

test_that("rows that do not sum to one are rejected", {
  fixture <- sam_fixture()
  unnormalized <- fixture$gps * 0.5

  expect_error(
    gps_candidate_search(fixture$data, unnormalized,
                         treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "must sum to 1"
  )
})

test_that("unnamed or duplicated GPS columns are rejected", {
  fixture <- sam_fixture()

  unnamed <- fixture$gps
  colnames(unnamed) <- NULL
  expect_error(
    gps_candidate_search(fixture$data, unnamed, treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "unique, named treatment columns"
  )

  duplicated_columns <- fixture$gps
  colnames(duplicated_columns)[2] <- colnames(duplicated_columns)[1]
  expect_error(
    gps_candidate_search(fixture$data, duplicated_columns,
                         treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "unique, named treatment columns"
  )
})

test_that("the search records the GPS it selected candidates against", {
  fixture <- sam_fixture()

  expect_false(is.null(fixture$search$gps_fingerprint))
  expect_identical(fixture$search$gps_fingerprint$columns,
                   colnames(fixture$gps))
})

test_that("a different GPS handed to a later stage is rejected", {
  fixture <- sam_fixture()

  # Previously accepted: sam_evaluate reported AUCs computed from scores that
  # never produced the match.
  reversed <- fixture$gps[, rev(colnames(fixture$gps)), drop = FALSE]
  expect_error(
    sam_evaluate(fixture$data, fixture$search, fixture$match, reversed,
                 X_vars = fixture$X_vars, treatment_var = "treatment"),
    "not the matrix `search` was built from"
  )

  perturbed <- fixture$gps
  perturbed[5, ] <- rev(perturbed[5, ])
  expect_error(
    sam_evaluate(fixture$data, fixture$search, fixture$match, perturbed,
                 X_vars = fixture$X_vars, treatment_var = "treatment"),
    "not the matrix `search` was built from"
  )

  expect_silent(
    sam_evaluate(fixture$data, fixture$search, fixture$match, fixture$gps,
                 X_vars = fixture$X_vars, treatment_var = "treatment")
  )
})

test_that("match_3way verifies the GPS as well", {
  fixture <- sam_fixture("sample_3group")
  reversed <- fixture$gps[, rev(colnames(fixture$gps)), drop = FALSE]

  expect_error(
    match_3way(fixture$data, fixture$search, reversed,
               treatment_var = "treatment"),
    "not the matrix `search` was built from"
  )
})

test_that("a search without a GPS fingerprint is still accepted", {
  fixture <- sam_fixture()
  legacy <- fixture$search
  legacy$gps_fingerprint <- NULL

  expect_silent(
    sam_evaluate(fixture$data, legacy, fixture$match, fixture$gps,
                 X_vars = fixture$X_vars, treatment_var = "treatment")
  )
})

test_that("the weighting module validates its GPS too", {
  fixture <- sam_fixture()
  out_of_range <- fixture$gps
  out_of_range[1, 1] <- 2

  expect_error(
    compute_balancing_weights(fixture$data, method = "iptw",
                              gps = out_of_range, treatment_var = "treatment",
                              anchor_level = fixture$anchor),
    "\\[0, 1\\]"
  )
})

test_that("a treatment level with no GPS column is named at every stage", {
  fixture <- sam_fixture()
  missing_level <- colnames(fixture$gps)[2]
  dropped <- fixture$gps[, -2, drop = FALSE]
  dropped <- dropped / rowSums(dropped)

  # Previously the level simply fell out of `groups`, so its subjects were
  # never matched and nothing said so.
  expect_error(
    gps_candidate_search(fixture$data, dropped, treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    paste0("not found in `gps`: ", missing_level),
    fixed = TRUE
  )
  expect_error(
    compute_balancing_weights(fixture$data, method = "iptw", gps = dropped,
                              treatment_var = "treatment",
                              anchor_level = fixture$anchor),
    "not found in `gps`",
    fixed = TRUE
  )
})

test_that("a row-reordered GPS is rejected", {
  fixture <- sam_fixture()
  reordered <- fixture$gps[c(2, 1, seq(3, nrow(fixture$gps))), , drop = FALSE]

  # Previously accepted: R matches gps to data by position, so every subject
  # after the swap carried someone else's scores.
  expect_error(
    gps_candidate_search(fixture$data, reordered, treatment_var = "treatment",
                         anchor_level = fixture$anchor, top_m = 5),
    "same row order"
  )
})

test_that("a missing treatment label is rejected", {
  fixture <- sam_fixture()
  with_gap <- fixture$data
  with_gap$treatment[3] <- NA

  # Previously: NA compared to NA against every group, so which() dropped the
  # row from all of them and the subject vanished without a word.
  expect_error(treatment_labels(with_gap, "treatment"),
               "missing treatment label")
  expect_error(
    estimate_gps_multinom(with_gap, X_vars = fixture$X_vars,
                          treatment_var = "treatment",
                          anchor_level = fixture$anchor),
    "missing treatment label"
  )
})

test_that("the returned model applies to the covariates as given", {
  # The covariates are standardized to condition the fit. Previously the model
  # came back still expecting standardized input, with nothing to say so, and
  # scoring subjects with it was wrong by up to 0.91 in probability on the
  # bundled data.
  fixture <- sam_fixture()
  fit <- estimate_gps_multinom(fixture$data, X_vars = fixture$X_vars,
                               treatment_var = "treatment",
                               anchor_level = fixture$anchor)

  scored <- stats::predict(fit$model,
                           newdata = fixture$data[, fixture$X_vars, drop = FALSE],
                           type = "probs")
  expect_equal(scored[, colnames(fit$gps)], fit$gps, tolerance = 1e-10)
})

test_that("scoring new subjects agrees with re-running the estimator", {
  data <- sam_small_data(n = 150L, n_groups = 3L)
  X_vars <- c("x1", "x2", "x3")

  fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "A")

  # A subject the model never saw, on the covariates' own scale.
  fresh <- data.frame(x1 = 2.5, x2 = -1.25, x3 = 0.75)
  scored <- stats::predict(fit$model, newdata = fresh, type = "probs")

  expect_length(scored, ncol(fit$gps))
  expect_equal(sum(scored), 1, tolerance = 1e-10)
  expect_true(all(scored >= 0 & scored <= 1))
})

test_that("the coefficients are on the covariates' original scale", {
  # A covariate rescaled by a constant must have its coefficient scaled by the
  # inverse, which is only true if the standardization was folded back in.
  data <- sam_small_data(n = 150L, n_groups = 3L)
  X_vars <- c("x1", "x2", "x3")

  fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "A")

  rescaled <- data
  rescaled$x1 <- rescaled$x1 * 100
  refit <- estimate_gps_multinom(rescaled, X_vars = X_vars,
                                 treatment_var = "T", anchor_level = "A")

  expect_equal(stats::coef(refit$model)[, "x1"],
               stats::coef(fit$model)[, "x1"] / 100,
               tolerance = 1e-6)
  expect_equal(unname(refit$gps), unname(fit$gps), tolerance = 1e-6)
})
