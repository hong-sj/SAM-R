################################################################################
### Candidate-count parameter validation (#5)
################################################################################

test_that("require_positive_int accepts only positive whole numbers", {
  expect_identical(require_positive_int(10, "top_m"), 10L)
  expect_identical(require_positive_int(1L, "top_m"), 1L)
  expect_identical(require_positive_int(3, "top_m"), 3L)

  # `TRUE` coerces to the integer 1 in R, but is not a candidate count.
  expect_error(require_positive_int(TRUE, "top_m"), "positive integer")

  expect_error(require_positive_int(0, "top_m"), "positive integer")
  expect_error(require_positive_int(-1, "top_m"), "positive integer")
  expect_error(require_positive_int(2.7, "top_m"), "positive integer")
  expect_error(require_positive_int(c(1, 2), "top_m"), "positive integer")
  expect_error(require_positive_int(NA_integer_, "top_m"), "positive integer")
  expect_error(require_positive_int(Inf, "top_m"), "positive integer")
  expect_error(require_positive_int("5", "top_m"), "positive integer")
})

test_that("the error message names the offending parameter", {
  expect_error(require_positive_int(0, "top_n"), "top_n", fixed = TRUE)
})

test_that("gps_candidate_search rejects a nonsensical top_m", {
  fixture <- sam_fixture()

  search_with <- function(top_m) {
    gps_candidate_search(
      fixture$data, fixture$gps,
      treatment_var = "treatment",
      anchor_level = fixture$anchor,
      top_m = top_m,
      gps_space = "logit"
    )
  }

  # Previously: a silent 0% match, reported as a legitimate result.
  expect_error(search_with(0), "top_m")
  # Previously: an error from deep inside seq_len().
  expect_error(search_with(-1), "top_m")
  # Previously: silently truncated to 2.
  expect_error(search_with(2.7), "top_m")
  expect_error(search_with(TRUE), "top_m")

  expect_silent(search_with(3L))
})

test_that("match_3way rejects a nonsensical top_n", {
  fixture <- sam_fixture("sample_3group")

  match_with <- function(top_n) {
    match_3way(
      fixture$data, fixture$search, fixture$gps,
      treatment_var = "treatment",
      top_n = top_n
    )
  }

  expect_error(match_with(0), "top_n")
  expect_error(match_with(-1), "top_n")
  expect_error(match_with(2.5), "top_n")
  expect_error(match_with(TRUE), "top_n")
})

test_that("a validated top_m still keeps at most top_m candidates", {
  fixture <- sam_fixture()

  for (top_m in c(1L, 4L, 10L)) {
    search <- gps_candidate_search(
      fixture$data, fixture$gps,
      treatment_var = "treatment",
      anchor_level = fixture$anchor,
      top_m = top_m,
      gps_space = "logit"
    )
    sizes <- unlist(lapply(search$candidates, function(x) {
      vapply(x, length, integer(1))
    }))
    expect_true(all(sizes <= top_m))
    expect_true(max(sizes) == top_m)
  }
})
