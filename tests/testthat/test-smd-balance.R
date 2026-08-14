################################################################################
### Guarded SMD denominator and unassessable covariates (#4)
################################################################################

test_that("a well-behaved cohort reports nothing undefined", {
  fixture <- sam_fixture()
  balance <- compute_smd_balance(fixture$data, fixture$match$matched,
                                 fixture$X_vars, fixture$groups)

  expect_true(all(balance$by_covariate$smd_defined))
  expect_true(all(balance$summary$n_undefined == 0L))
  expect_true(all(is.finite(balance$summary$mean_abs_smd)))
  expect_true(all(is.finite(balance$summary$max_abs_smd)))
})

test_that("a zero-variance covariate no longer destroys the summary", {
  fixture <- sam_fixture()
  before <- compute_smd_balance(fixture$data, fixture$match$matched,
                                fixture$X_vars, fixture$groups)

  degenerate <- fixture$data
  degenerate$constant <- 1
  after <- compute_smd_balance(degenerate, fixture$match$matched,
                               c(fixture$X_vars, "constant"), fixture$groups)

  # Previously: mean_abs_smd and max_abs_smd became NaN for every group,
  # taking nine perfectly assessable covariates down with the tenth.
  expect_true(all(is.finite(after$summary$mean_abs_smd)))
  expect_true(all(is.finite(after$summary$max_abs_smd)))

  # The assessable covariates are unaffected.
  expect_equal(after$summary$mean_abs_smd, before$summary$mean_abs_smd)
  expect_equal(after$summary$max_abs_smd, before$summary$max_abs_smd)

  # And the covariate that could not be assessed is visible rather than absent.
  expect_equal(after$summary$n_undefined, rep(1L, nrow(after$summary)))
  undefined <- after$by_covariate[!after$by_covariate$smd_defined, ]
  expect_equal(unique(undefined$covariate), "constant")
  expect_true(all(undefined$smd == 0))
})

test_that("differing means with zero variance do not produce an infinity", {
  # The other failure mode of the same expression: 1/0 rather than 0/0.
  data <- data.frame(
    x = c(rep(1, 4), rep(2, 4)),
    T = rep(c("A", "B"), each = 4),
    stringsAsFactors = FALSE
  )
  matched <- data.frame(anchor = 1:4, B = 5:8)

  balance <- compute_smd_balance(data, matched, "x", "B")

  expect_false(balance$by_covariate$smd_defined)
  expect_equal(balance$by_covariate$smd, 0)
  expect_true(is.na(balance$summary$mean_abs_smd))
  expect_equal(balance$summary$n_undefined, 1L)
})

test_that("compute_weighted_balance carries the same guard", {
  fixture <- sam_fixture()
  degenerate <- fixture$data
  degenerate$constant <- 1

  weights <- compute_balancing_weights(
    degenerate, method = "overlap", X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )
  balance <- compute_weighted_balance(
    degenerate, weights$weights, X_vars = c(fixture$X_vars, "constant"),
    treatment_var = "treatment", anchor_level = fixture$anchor
  )

  expect_true(all(is.finite(balance$summary$mean_abs_smd)))
  expect_equal(balance$summary$n_undefined, rep(1L, nrow(balance$summary)))
  expect_true("smd_defined" %in% names(balance$by_covariate))
})

test_that("both balance tables report the same summary columns", {
  fixture <- sam_fixture()

  matched_balance <- compute_smd_balance(fixture$data, fixture$match$matched,
                                         fixture$X_vars, fixture$groups)
  weights <- compute_balancing_weights(
    fixture$data, method = "matching", X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )
  weighted_balance <- compute_weighted_balance(
    fixture$data, weights$weights, X_vars = fixture$X_vars,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )

  expect_identical(names(matched_balance$summary),
                   names(weighted_balance$summary))
})

test_that("a covariate that is undefined in every group yields NA, not NaN", {
  data <- data.frame(
    x = rep(0, 8),
    T = rep(c("A", "B"), each = 4),
    stringsAsFactors = FALSE
  )
  matched <- data.frame(anchor = 1:4, B = 5:8)

  balance <- compute_smd_balance(data, matched, "x", "B")

  expect_true(is.na(balance$summary$mean_abs_smd))
  expect_true(is.na(balance$summary$max_abs_smd))
  expect_false(is.nan(balance$summary$mean_abs_smd))
})
