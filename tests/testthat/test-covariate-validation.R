################################################################################
### Covariate validation before the linear algebra (#2)
################################################################################

test_that("solve() does not raise on a non-finite matrix", {
  # This is the reason the check exists at all: the tryCatch() around solve()
  # never fires for missing data, so an is.finite() test on the result is what
  # actually catches it.
  #
  # Exactly what comes back is a LAPACK detail and differs between platforms --
  # all NA on some, a mix of NA and NaN on others -- so only the two properties
  # get_pooled_covariance() relies on are asserted.
  S <- matrix(c(1, 0.5, 0.5, 1), 2)
  S[1, 1] <- NA_real_

  expect_silent(inverse <- solve(S))
  expect_false(all(is.finite(inverse)))
})

test_that("a single missing covariate cell raises instead of emptying the match", {
  fixture <- sam_fixture()
  broken <- fixture$data
  broken[[fixture$X_vars[1]]][5] <- NA_real_

  # Previously: matching_rate fell from 0.132 to 0, with no error and no
  # warning, because the all-NA precision matrix made every distance NA and
  # every anchor look exhausted.
  expect_error(
    sam_match(broken, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment"),
    "missing or non-finite"
  )
  expect_gt(fixture$match$matching_rate, 0)
})

test_that("the message names the offending columns and their row counts", {
  fixture <- sam_fixture()
  broken <- fixture$data
  broken[[fixture$X_vars[1]]][c(2, 7, 11)] <- NA_real_
  broken[[fixture$X_vars[2]]][4] <- Inf

  expect_error(
    get_pooled_covariance(broken, fixture$X_vars, "treatment"),
    paste0(fixture$X_vars[1], " \\(3 rows\\)")
  )
  expect_error(
    get_pooled_covariance(broken, fixture$X_vars, "treatment"),
    paste0(fixture$X_vars[2], " \\(1 rows\\)")
  )
})

test_that("non-numeric covariates are rejected with an actionable message", {
  fixture <- sam_fixture()
  with_character <- fixture$data
  with_character$site <- "A"

  expect_error(
    sam_match(with_character, fixture$search,
              X_vars = c(fixture$X_vars, "site"), treatment_var = "treatment"),
    "not numeric"
  )
  expect_error(
    sam_match(with_character, fixture$search,
              X_vars = c(fixture$X_vars, "site"), treatment_var = "treatment"),
    "model.matrix"
  )

  # A factor is not silently reduced to its integer codes either.
  with_factor <- fixture$data
  with_factor$site <- factor(rep(c("A", "B"), length.out = nrow(with_factor)))
  expect_error(
    get_pooled_covariance(with_factor, c(fixture$X_vars, "site"), "treatment"),
    "not numeric"
  )

  # Logicals are accepted as 0/1 indicators.
  with_logical <- fixture$data
  with_logical$flag <- rep(c(TRUE, FALSE), length.out = nrow(with_logical))
  expect_silent(
    get_pooled_covariance(with_logical, c(fixture$X_vars, "flag"), "treatment")
  )
})

test_that("covariates used only for matching are validated too", {
  # sam_match deliberately allows a different covariate set from the one GPS
  # screening used, so these columns are not covered by anything upstream.
  fixture <- sam_fixture()
  broken <- fixture$data
  broken$extra <- c(NA_real_, rep(0, nrow(broken) - 1L))

  expect_error(
    sam_match(broken, fixture$search, X_vars = c(fixture$X_vars, "extra"),
              treatment_var = "treatment"),
    "extra \\(1 rows\\)"
  )
})

test_that("a missing covariate column is named", {
  fixture <- sam_fixture()

  expect_error(
    get_pooled_covariance(fixture$data, c(fixture$X_vars, "nope"), "treatment"),
    "not found in data: nope"
  )
})

test_that("non-positive residual degrees of freedom raise", {
  data <- data.frame(x1 = c(1, 2, 3), x2 = c(4, 5, 6),
                     T = c("a", "b", "c"), stringsAsFactors = FALSE)

  expect_error(
    get_pooled_covariance(data, c("x1", "x2"), "T"),
    "Residual degrees of freedom"
  )
})

test_that("a rank-deficient covariance still yields a finite inverse", {
  fixture <- sam_fixture()
  collinear <- fixture$data
  collinear$duplicate <- collinear[[fixture$X_vars[1]]]

  # Whether solve() raises on an exactly singular matrix, returns a non-finite
  # result, or returns finite values from a barely-invertible one is a LAPACK
  # detail that differs between platforms. What get_pooled_covariance()
  # guarantees on all of them is a finite inverse; when it reaches
  # MASS::ginv() to get one, it says so.
  warned <- NULL
  pooled <- withCallingHandlers(
    get_pooled_covariance(collinear, c(fixture$X_vars, "duplicate"),
                          "treatment"),
    warning = function(w) {
      warned <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  expect_true(all(is.finite(pooled$S_inv)))
  if (!is.null(warned)) {
    expect_match(warned, "numerically singular")
  }
})

test_that("estimate_gps_multinom validates its covariates as well", {
  fixture <- sam_fixture()
  broken <- fixture$data
  broken[[fixture$X_vars[1]]][1] <- NA_real_

  expect_error(
    estimate_gps_multinom(broken, X_vars = fixture$X_vars,
                          treatment_var = "treatment",
                          anchor_level = fixture$anchor),
    "missing or non-finite"
  )
})
