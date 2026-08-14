################################################################################
### Shared result assembly and public API alignment (#8)
################################################################################

test_that("transform_ps is the one logit transform", {
  values <- matrix(c(0.5, 0.25, 0.5, 0.75), 2)

  expect_identical(transform_ps(values, "raw"), values)
  expect_equal(transform_ps(values, "logit"), stats::qlogis(values))

  # Boundary probabilities are clipped rather than mapped to infinities.
  boundary <- matrix(c(0, 1), 1)
  expect_true(all(is.finite(transform_ps(boundary, "logit"))))

  expect_error(transform_ps(values, "nope"))
})

test_that("an empty matched frame keeps the full column schema", {
  groups <- c("B", "C")
  empty <- empty_matched_frame(groups)

  expect_equal(nrow(empty), 0L)
  expect_identical(
    names(empty),
    c("matched_set_id", "anchor", "B", "C", "dist_B", "dist_C", "loss")
  )

  with_extra <- empty_matched_frame(groups, "rassen_perimeter")
  expect_identical(names(with_extra)[7:8], c("loss", "rassen_perimeter"))
})

test_that("both engines return the same column schema", {
  four <- sam_fixture()
  three <- sam_fixture("sample_3group")
  trio <- match_3way(three$data, three$search, three$gps,
                     treatment_var = "treatment")

  expect_identical(names(four$match$matched),
                   matched_frame_columns(four$groups))

  # match_3way carries one extra column, appended after loss.
  expect_identical(names(trio$matched),
                   matched_frame_columns(three$groups, "rassen_perimeter"))
})

test_that("groups_from_matched recovers the comparator groups", {
  four <- sam_fixture()
  expect_identical(groups_from_matched(four$match$matched), four$groups)

  three <- sam_fixture("sample_3group")
  trio <- match_3way(three$data, three$search, three$gps,
                     treatment_var = "treatment")
  # rassen_perimeter must not be mistaken for a group.
  expect_identical(groups_from_matched(trio$matched), three$groups)
})

test_that("max_possible_rate is the ceiling the group sizes impose", {
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")

  group_sizes <- vapply(fixture$groups, function(g) sum(labels == g), integer(1))
  n_anchor <- length(fixture$search$anchor_rows)
  expected <- min(group_sizes) / n_anchor

  expect_equal(fixture$match$max_possible_rate, expected)
  expect_lte(fixture$match$matching_rate,
             fixture$match$max_possible_rate + 1e-12)

  # On the bundled four-group data SAM saturates the ceiling exactly, which is
  # what makes a 13% rate a complete result rather than a poor one.
  expect_equal(fixture$match$matching_rate, fixture$match$max_possible_rate)
  expect_equal(min(group_sizes), 59L)
  expect_equal(n_anchor, 448L)
})

test_that("match_3way reports max_possible_rate too", {
  fixture <- sam_fixture("sample_3group")
  trio <- match_3way(fixture$data, fixture$search, fixture$gps,
                     treatment_var = "treatment")

  expect_true("max_possible_rate" %in% names(trio))
  expect_lte(trio$matching_rate, trio$max_possible_rate + 1e-12)
})

test_that("match_3way takes gps_space and deprecates ps_space", {
  fixture <- sam_fixture("sample_3group")

  # Previously: gps_space was an "unused argument" error, and the README used
  # both names within one section.
  by_new_name <- match_3way(fixture$data, fixture$search, fixture$gps,
                            treatment_var = "treatment", gps_space = "logit")

  expect_warning(
    by_old_name <- match_3way(fixture$data, fixture$search, fixture$gps,
                              treatment_var = "treatment", ps_space = "logit"),
    "deprecated"
  )
  expect_equal(by_new_name$matched, by_old_name$matched)
})

test_that("match_3way warns instead of discarding X_vars", {
  fixture <- sam_fixture("sample_3group")

  # Previously: silently identical, even for a column that does not exist.
  expect_warning(
    with_covariates <- match_3way(fixture$data, fixture$search, fixture$gps,
                                  X_vars = "nonexistent_column",
                                  treatment_var = "treatment"),
    "ignored by match_3way"
  )
  without <- match_3way(fixture$data, fixture$search, fixture$gps,
                        treatment_var = "treatment")
  expect_equal(with_covariates$matched, without$matched)
})

test_that("match_3way rejects a search that is not three-group", {
  fixture <- sam_fixture()

  expect_error(
    match_3way(fixture$data, fixture$search, fixture$gps,
               treatment_var = "treatment"),
    "exactly three treatment groups"
  )
})

test_that("an empty anchor group returns the schema, not a bare NaN", {
  fixture <- sam_fixture()
  empty_search <- fixture$search
  empty_search$anchor_rows <- integer(0)
  empty_search$candidates <- list()

  result <- sam_match(fixture$data, empty_search, X_vars = fixture$X_vars,
                      treatment_var = "treatment")

  expect_equal(nrow(result$matched), 0L)
  expect_identical(names(result$matched),
                   matched_frame_columns(fixture$groups))
  expect_true(is.nan(result$matching_rate))
  expect_true(is.nan(result$max_possible_rate))
  expect_length(result$unmatched_anchor_rows, 0L)
})

test_that("the diagnostics infer what their siblings default", {
  fixture <- sam_fixture()

  explicit <- compute_smd_balance(fixture$data, fixture$match$matched,
                                  X_vars = fixture$X_vars,
                                  groups = fixture$groups)
  inferred <- compute_smd_balance(fixture$data, fixture$match$matched,
                                  X_vars = fixture$X_vars)
  expect_equal(explicit, inferred)

  explicit_auc <- compute_pairwise_treatment_auc(
    fixture$gps, fixture$match$matched,
    groups = fixture$groups, anchor_level = fixture$anchor
  )
  inferred_auc <- compute_pairwise_treatment_auc(fixture$gps,
                                                 fixture$match$matched)
  expect_equal(explicit_auc, inferred_auc)
})

test_that("get_pooled_covariance has usable defaults", {
  data <- sam_small_data(n = 90L, n_groups = 3L)
  names(data)[1:3] <- paste0("X", 1:3)
  data <- data[, c(paste0("X", 1:3), "T")]

  # X_vars defaults to paste0("X", 1:10), which is not present here, so the
  # point is only that treatment_var defaults to "T" like everywhere else.
  pooled <- get_pooled_covariance(data, X_vars = paste0("X", 1:3))
  expect_equal(dim(pooled$S), c(3L, 3L))
})

test_that("calc_caliper_3way refuses input it cannot summarise", {
  ps <- matrix(stats::runif(60), 30, 2)

  expect_error(calc_caliper_3way(ps, rep(c("A", "B"), each = 15)),
               "exactly three treatment groups")

  # stats::var() of a single observation is NA, which would otherwise reach the
  # caliper silently.
  expect_error(
    calc_caliper_3way(ps, c(rep("A", 14), rep("B", 15), "C")),
    "at least two subjects"
  )

  non_finite <- ps
  non_finite[1, 1] <- NA_real_
  expect_error(calc_caliper_3way(non_finite, rep(c("A", "B", "C"), each = 10)),
               "only finite values")
})

test_that("match_3way validates the caliper and the reference level", {
  fixture <- sam_fixture("sample_3group")

  match_with <- function(...) {
    match_3way(fixture$data, fixture$search, fixture$gps,
               treatment_var = "treatment", ...)
  }

  expect_error(match_with(caliper = 0), "finite number greater than zero")
  expect_error(match_with(caliper = -1), "finite number greater than zero")
  expect_error(match_with(caliper = Inf), "finite number greater than zero")
  expect_error(match_with(caliper = c(1, 2)), "finite number greater than zero")

  # A reference level given as a number must name a level, not a position.
  expect_error(match_with(reference_level = 1))
  expect_silent(match_with(reference_level = fixture$anchor))
})
