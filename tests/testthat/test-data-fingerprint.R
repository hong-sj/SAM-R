################################################################################
### The frame every later stage must be handed (#3)
################################################################################

test_that("gps_candidate_search records a fingerprint", {
  fixture <- sam_fixture()

  expect_false(is.null(fixture$search$fingerprint))
  expect_equal(fixture$search$fingerprint$n_rows, nrow(fixture$data))
})

test_that("re-sorting data between stages is rejected", {
  fixture <- sam_fixture()
  reordered <- fixture$data[order(fixture$data[[fixture$X_vars[1]]]), ]

  # Previously: max_abs_smd moved from 0.104/0.393/0.224 to 0.449/0.486/0.523
  # and described a cohort that was never matched, with no warning.
  expect_error(
    sam_evaluate(reordered, fixture$search, fixture$match, fixture$gps,
                 X_vars = fixture$X_vars, treatment_var = "treatment"),
    "not the frame"
  )
  expect_error(
    sam_match(reordered, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment"),
    "not the frame"
  )
  expect_error(
    extract_matched_data(reordered, fixture$search, fixture$match,
                         treatment_var = "treatment"),
    "not the frame"
  )
})

test_that("a reorder followed by resetting the row names is still caught", {
  # Row names alone would not catch this: they are back to 1:n afterwards.
  fixture <- sam_fixture()
  reordered <- fixture$data[order(fixture$data[[fixture$X_vars[1]]]), ]
  rownames(reordered) <- NULL

  expect_error(
    sam_match(reordered, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment"),
    "treatment column changed"
  )
})

test_that("a reorder within one treatment group is still caught", {
  # The treatment column alone would not catch this: reversing the anchor rows
  # among themselves leaves every label exactly where it was.
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")
  anchor_rows <- which(labels == fixture$anchor)

  permutation <- seq_len(nrow(fixture$data))
  permutation[anchor_rows] <- rev(anchor_rows)
  shuffled <- fixture$data[permutation, ]

  expect_identical(treatment_labels(shuffled, "treatment"), labels)
  expect_error(
    sam_match(shuffled, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment"),
    "row names changed"
  )
})

test_that("rewriting row contents in place is a documented blind spot", {
  # The fingerprint records row identity and the treatment column, not the
  # covariates, so that attaching an outcome column stays legal. Overwriting
  # the *contents* of existing rows therefore passes: it changes neither the
  # row names nor the labels. This is a deliberate trade-off, pinned here so
  # that the limit is visible rather than assumed away.
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")
  anchor_rows <- which(labels == fixture$anchor)

  rewritten <- fixture$data
  rewritten[anchor_rows, ] <- rewritten[rev(anchor_rows), ]

  expect_silent(
    sam_match(rewritten, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment")
  )
})

test_that("filtering rows is rejected and the counts are named", {
  fixture <- sam_fixture()
  filtered <- fixture$data[-1, ]

  expect_error(
    sam_match(filtered, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment"),
    "row count changed"
  )
})

test_that("attaching a column between stages stays legal", {
  fixture <- sam_fixture()
  with_outcome <- fixture$data
  with_outcome$new_outcome <- rep(0L, nrow(with_outcome))

  expect_silent(
    sam_match(with_outcome, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment")
  )
  expect_equal(
    sam_match(with_outcome, fixture$search, X_vars = fixture$X_vars,
              treatment_var = "treatment")$matched,
    fixture$match$matched
  )
})

test_that("a search without a fingerprint is still accepted", {
  # Objects created by earlier versions of the package keep working.
  fixture <- sam_fixture()
  legacy <- fixture$search
  legacy$fingerprint <- NULL

  expect_silent(
    result <- sam_match(fixture$data, legacy, X_vars = fixture$X_vars,
                        treatment_var = "treatment")
  )
  expect_equal(result$matched, fixture$match$matched)
})

test_that("match_3way verifies the fingerprint too", {
  fixture <- sam_fixture("sample_3group")
  reordered <- fixture$data[rev(seq_len(nrow(fixture$data))), ]

  expect_error(
    match_3way(reordered, fixture$search, fixture$gps,
               treatment_var = "treatment"),
    "not the frame"
  )
})
