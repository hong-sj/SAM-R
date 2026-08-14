################################################################################
### The greedy engine, against a direct restatement of the matching rule (#7)
################################################################################

#' Reference implementation of Shared Anchor Matching
#'
#' A deliberately naive restatement of the rule the documentation gives:
#' repeatedly take the globally cheapest matched set still available, where a
#' matched set is one anchor plus its nearest still-available candidate from
#' each comparator group. No pointers, no reverse index, no incremental
#' updates -- everything is recomputed from scratch on every iteration.
reference_sam_match <- function(data, search, X_vars, treatment_var) {
  labels <- treatment_labels(data, treatment_var)
  anchor_rows <- as.integer(search$anchor_rows)
  groups <- as.character(search$groups)
  S_inv <- get_pooled_covariance(data, X_vars, treatment_var)$S_inv
  X <- covariate_matrix(data, X_vars)

  available <- stats::setNames(
    lapply(groups, function(g) which(labels == g)), groups
  )
  remaining <- seq_along(anchor_rows)
  results <- list()

  repeat {
    best <- NULL

    for (i in remaining) {
      pick <- list()
      total <- 0

      for (g in groups) {
        candidates <- intersect(search$candidates[[i]][[g]], available[[g]])
        if (length(candidates) == 0L) {
          pick <- NULL
          break
        }
        distances <- as.numeric(mahalanobis_distance_matrix(
          X[anchor_rows[i], , drop = FALSE],
          X[candidates, , drop = FALSE], S_inv
        ))
        # Ties resolve by row order, as a stable sort would.
        winner <- which.min(distances)
        pick[[g]] <- list(row = candidates[winner], dist = distances[winner])
        total <- total + distances[winner]
      }

      if (!is.null(pick) && (is.null(best) || total < best$total)) {
        best <- list(anchor = i, pick = pick, total = total)
      }
    }

    if (is.null(best)) {
      break
    }

    results[[length(results) + 1L]] <- best
    remaining <- setdiff(remaining, best$anchor)
    for (g in groups) {
      available[[g]] <- setdiff(available[[g]], best$pick[[g]]$row)
    }
    if (length(remaining) == 0L) {
      break
    }
  }

  frame <- data.frame(
    matched_set_id = seq_along(results),
    anchor = vapply(results, function(r) anchor_rows[r$anchor], integer(1))
  )
  for (g in groups) {
    frame[[g]] <- vapply(results, function(r) as.integer(r$pick[[g]]$row),
                         integer(1))
  }
  for (g in groups) {
    frame[[paste0("dist_", g)]] <- vapply(results,
                                          function(r) r$pick[[g]]$dist,
                                          numeric(1))
  }
  frame$loss <- vapply(results, function(r) r$total, numeric(1))
  frame
}

test_that("the engine reproduces the rule on small random data", {
  for (seed in 1:4) {
    data <- sam_small_data(n = 90L, n_groups = 3L, seed = seed)
    X_vars <- c("x1", "x2", "x3")

    fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                                 anchor_level = "A")
    search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                   anchor_level = "A", top_m = 4)

    engine <- sam_match(data, search, X_vars = X_vars, treatment_var = "T")
    reference <- reference_sam_match(data, search, X_vars, "T")

    expect_equal(engine$matched, reference, tolerance = 1e-10,
                 info = paste("seed", seed))
  }
})

test_that("the engine reproduces the rule on the bundled data", {
  fixture <- sam_fixture()

  reference <- reference_sam_match(fixture$data, fixture$search,
                                   fixture$X_vars, "treatment")
  expect_equal(fixture$match$matched, reference, tolerance = 1e-10)
})

test_that("globally cheapest wins, not anchor-by-anchor", {
  # Two anchors compete for one comparator subject. Taking them in row order
  # would give it to the anchor at x = 0; the rule gives it to the anchor at
  # x = 3, whose set is cheaper, and leaves the first with the distant one.
  data <- data.frame(
    x = c(0, 3, 2, 8),
    T = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  search <- list(
    anchor_rows = c(1L, 2L),
    groups = "B",
    candidates = list(list(B = c(3L, 4L)), list(B = c(3L, 4L)))
  )

  result <- sam_match(data, search, X_vars = "x", treatment_var = "T")

  expect_equal(nrow(result$matched), 2L)
  # The cheaper set (anchor at x = 3, distance 1) is taken first.
  expect_equal(result$matched$anchor, c(2L, 1L))
  expect_equal(result$matched$B, c(3L, 4L))
  expect_lt(result$matched$loss[1], result$matched$loss[2])
})

test_that("no subject is used twice and every matched set is complete", {
  fixture <- sam_fixture()
  matched <- fixture$match$matched

  expect_equal(anyDuplicated(matched$anchor), 0L)
  for (g in fixture$groups) {
    expect_equal(anyDuplicated(matched[[g]]), 0L)
    expect_false(anyNA(matched[[g]]))
  }
  expect_equal(matched$loss,
               rowSums(as.matrix(matched[, paste0("dist_", fixture$groups)])),
               tolerance = 1e-10)
})

test_that("an anchor whose candidates are all consumed is left unmatched", {
  # Both anchors can only reach one comparator subject, so exactly one matches.
  data <- data.frame(
    x = c(0, 0.5, 1),
    T = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  search <- list(
    anchor_rows = c(1L, 2L),
    groups = "B",
    candidates = list(list(B = 3L), list(B = 3L))
  )

  result <- sam_match(data, search, X_vars = "x", treatment_var = "T")

  expect_equal(nrow(result$matched), 1L)
  expect_equal(result$matched$anchor, 2L)
  expect_equal(result$unmatched_anchor_rows, 1L)
  expect_equal(result$matching_rate, 0.5)
  expect_equal(result$max_possible_rate, 0.5)
})

test_that("the engine is deterministic", {
  fixture <- sam_fixture()

  repeated <- sam_match(fixture$data, fixture$search,
                        X_vars = fixture$X_vars, treatment_var = "treatment")
  expect_identical(repeated$matched, fixture$match$matched)
})

test_that("the shortlist reproduces the rule when it is rebuilt constantly", {
  # sam_match() selects the cheapest anchor from a shortlist of the 1024
  # smallest losses rather than scanning every anchor. Test data never reaches
  # that many anchors, so the rebuild path would otherwise never run. Shrinking
  # the shortlist to three forces a rebuild every few matches.
  withr::local_options(SAMatch.shortlist_size = 3L)

  for (seed in 1:4) {
    data <- sam_small_data(n = 90L, n_groups = 3L, seed = seed)
    X_vars <- c("x1", "x2", "x3")

    fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                                 anchor_level = "A")
    search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                   anchor_level = "A", top_m = 4)

    engine <- sam_match(data, search, X_vars = X_vars, treatment_var = "T")
    reference <- reference_sam_match(data, search, X_vars, "T")

    expect_equal(engine$matched, reference, tolerance = 1e-10,
                 info = paste("seed", seed))
  }
})

test_that("the shortlist size does not change the result", {
  # Whatever the shortlist size, the selection rule is the same, so every
  # setting must produce the same matched sets.
  fixture <- sam_fixture()
  expected <- fixture$match$matched

  for (size in c(1L, 2L, 7L, 50L, 10000L)) {
    withr::local_options(SAMatch.shortlist_size = size)
    result <- sam_match(fixture$data, fixture$search, X_vars = fixture$X_vars,
                        treatment_var = "treatment")
    expect_equal(result$matched, expected, info = paste("size", size))
  }
})

test_that("a shortlist rebuild preserves the tie-break on tied data", {
  # Exact ties are where the shortlist's fence matters: an anchor sitting just
  # outside it at the same loss must not be overlooked.
  withr::local_options(SAMatch.shortlist_size = 2L)

  set.seed(31L)
  n <- 120L
  data <- as.data.frame(matrix(round(stats::rnorm(n * 3)), n, 3,
                               dimnames = list(NULL, c("x1", "x2", "x3"))))
  data$T <- rep(c("A", "B", "C"), length.out = n)

  fit <- estimate_gps_multinom(data, X_vars = c("x1", "x2", "x3"),
                               treatment_var = "T", anchor_level = "A")
  search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                 anchor_level = "A", top_m = 4)

  engine <- sam_match(data, search, X_vars = c("x1", "x2", "x3"),
                      treatment_var = "T")
  reference <- reference_sam_match(data, search, c("x1", "x2", "x3"), "T")

  expect_equal(engine$matched, reference, tolerance = 1e-10)
})
