################################################################################
### A match that formed no sets stays usable downstream
##
## Producing no matched sets is a legitimate outcome -- match_3way() with a
## tight caliper produces one -- so the diagnostics have to describe it as an
## empty result rather than failing, or worse, describing it as something else.
################################################################################

empty_fixture <- function() {
  fixture <- sam_fixture()
  fixture$empty <- empty_matched_frame(fixture$groups)
  fixture
}

test_that("the empty frame carries the same column types as a populated one", {
  fixture <- sam_fixture()
  empty <- empty_matched_frame(fixture$groups)
  populated <- fixture$match$matched

  expect_identical(names(empty), names(populated))
  expect_identical(
    vapply(empty, function(x) class(x)[1], character(1)),
    vapply(populated, function(x) class(x)[1], character(1))
  )
  expect_equal(nrow(empty), 0L)
})

test_that("an empty match is an empty balance report, not unassessable covariates", {
  fixture <- empty_fixture()

  # Previously: every covariate came back with smd_defined = FALSE and
  # n_undefined = 10, which is what a constant covariate looks like. "Nothing
  # matched" and "this covariate has no variance" are different findings and
  # must not arrive through the same channel.
  balance <- compute_smd_balance(fixture$data, fixture$empty,
                                 X_vars = fixture$X_vars,
                                 groups = fixture$groups)

  expect_equal(nrow(balance$by_covariate), 0L)
  expect_equal(balance$summary$group, fixture$groups)
  expect_true(all(is.na(balance$summary$mean_abs_smd)))
  expect_true(all(is.na(balance$summary$max_abs_smd)))
  expect_true(all(balance$summary$n_undefined == 0L))

  # The column set is the same as for a real report, so downstream code that
  # reads it does not have to special-case the empty one.
  populated <- compute_smd_balance(fixture$data, fixture$match$matched,
                                   X_vars = fixture$X_vars,
                                   groups = fixture$groups)
  expect_identical(names(balance$by_covariate), names(populated$by_covariate))
  expect_identical(names(balance$summary), names(populated$summary))
})

test_that("an empty match gives an empty AUC report", {
  fixture <- empty_fixture()

  auc <- compute_pairwise_treatment_auc(fixture$gps, fixture$empty,
                                        groups = fixture$groups,
                                        anchor_level = fixture$anchor)

  expect_equal(nrow(auc$pairwise), 0L)
  expect_true(is.na(auc$mean_auc))

  populated <- compute_pairwise_treatment_auc(fixture$gps,
                                              fixture$match$matched,
                                              groups = fixture$groups,
                                              anchor_level = fixture$anchor)
  expect_identical(names(auc$pairwise), names(populated$pairwise))
})

test_that("sam_evaluate reports an empty match without failing", {
  fixture <- empty_fixture()
  empty_result <- fixture$match
  empty_result$matched <- fixture$empty

  report <- sam_evaluate(fixture$data, fixture$search, empty_result,
                         fixture$gps, X_vars = fixture$X_vars,
                         treatment_var = "treatment")

  expect_true(all(is.na(unlist(report$loss_distribution))))
  expect_true(all(is.na(unlist(report$dispersion_distribution))))
  expect_equal(nrow(report$smd_balance$by_covariate), 0L)
  expect_equal(nrow(report$treatment_discrimination_auc$pairwise), 0L)
})

test_that("a caliper that matches nothing is a result, not an error", {
  fixture <- sam_fixture("sample_3group")

  trio <- match_3way(fixture$data, fixture$search, fixture$gps,
                     treatment_var = "treatment", caliper = 1e-9)

  expect_equal(nrow(trio$matched), 0L)
  expect_identical(names(trio$matched),
                   matched_frame_columns(fixture$groups, "rassen_perimeter"))
  expect_equal(trio$matching_rate, 0)
  expect_setequal(trio$unmatched_anchor_rows, fixture$search$anchor_rows)

  # And the diagnostics accept what it produced.
  expect_silent(
    compute_smd_balance(fixture$data, trio$matched, X_vars = fixture$X_vars,
                        groups = fixture$groups)
  )
})

test_that("extract_matched_data still refuses an empty match", {
  # This one is not a report: there is no subject-level cohort to build, and
  # returning a zero-row frame would let an empty analysis run on silently.
  fixture <- empty_fixture()
  empty_result <- fixture$match
  empty_result$matched <- fixture$empty

  expect_error(
    extract_matched_data(fixture$data, fixture$search, empty_result,
                         treatment_var = "treatment"),
    "no matched sets"
  )
})
