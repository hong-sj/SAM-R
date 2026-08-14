################################################################################
### Candidate search against an explicit full-matrix reference (#6)
################################################################################

#' The implementation the k-d tree query replaced
#'
#' Builds the whole anchor-by-comparator distance matrix and fully sorts every
#' row. `order()` is a stable sort, so equidistant candidates come out in row
#' order -- which is the behaviour the tree query has to reproduce.
reference_candidates <- function(data, gps, treatment_var, anchor_level,
                                 top_m, gps_space) {
  labels <- treatment_labels(data, treatment_var)
  anchor_level <- as.character(anchor_level)

  gps_used <- if (gps_space == "logit") {
    eps <- 1e-6
    stats::qlogis(pmin(pmax(gps, eps), 1 - eps))
  } else {
    gps
  }

  groups <- setdiff(colnames(gps), anchor_level)
  anchor_rows <- which(labels == anchor_level)
  X_anchor <- gps_used[anchor_rows, , drop = FALSE]

  by_group <- lapply(groups, function(g) {
    group_rows <- which(labels == g)
    X_group <- gps_used[group_rows, , drop = FALSE]

    x_sq <- rowSums(X_anchor^2)
    y_sq <- rowSums(X_group^2)
    cross <- X_anchor %*% t(X_group)
    d2 <- outer(x_sq, y_sq, "+") - 2 * cross
    d2[d2 < 0] <- 0
    dist_mat <- sqrt(d2)

    m <- min(top_m, length(group_rows))
    lapply(seq_len(nrow(dist_mat)), function(i) {
      group_rows[order(dist_mat[i, ])[seq_len(m)]]
    })
  })
  names(by_group) <- groups

  lapply(seq_along(anchor_rows), function(i) {
    stats::setNames(lapply(groups, function(g) by_group[[g]][[i]]), groups)
  })
}

expect_matches_reference <- function(data, gps, top_m, gps_space,
                                     treatment_var = "T", anchor_level = "A") {
  search <- gps_candidate_search(data, gps, treatment_var = treatment_var,
                                 anchor_level = anchor_level, top_m = top_m,
                                 gps_space = gps_space)
  reference <- reference_candidates(data, gps, treatment_var, anchor_level,
                                    top_m, gps_space)

  expect_equal(search$candidates, reference,
               info = paste("top_m =", top_m, "space =", gps_space))
}

test_that("the tree query reproduces the full-matrix result", {
  data <- sam_small_data(n = 180L, n_groups = 4L)
  gps <- estimate_gps_multinom(data, X_vars = c("x1", "x2", "x3"),
                               treatment_var = "T", anchor_level = "A")$gps

  for (gps_space in c("raw", "logit")) {
    for (top_m in c(1L, 2L, 5L, 10L, 25L)) {
      expect_matches_reference(data, gps, top_m, gps_space)
    }
  }
})

test_that("it reproduces the full-matrix result on the bundled data", {
  data <- get(utils::data("sample_4group", package = "SAMatch",
                          envir = environment()))
  X_vars <- sam_covariates(data)
  anchor <- sam_anchor(data)
  gps <- estimate_gps_multinom(data, X_vars = X_vars,
                               treatment_var = "treatment",
                               anchor_level = anchor)$gps

  for (gps_space in c("raw", "logit")) {
    for (top_m in c(1L, 3L, 10L, 20L)) {
      expect_matches_reference(data, gps, top_m, gps_space,
                               treatment_var = "treatment",
                               anchor_level = anchor)
    }
  }
})

#' Data whose comparator subjects are exact duplicates of each other
#'
#' Every covariate is replicated, so the GPS rows are identical too and a large
#' number of candidates sit at exactly the same distance from each anchor. This
#' is the case where a k-d tree returns an arbitrary subset at the cutoff.
tied_candidate_data <- function(n_anchor = 30L, n_distinct = 6L,
                                copies = 6L, seed = 4242L) {
  set.seed(seed)
  covariates <- c("x1", "x2", "x3")

  anchor_block <- as.data.frame(
    matrix(stats::rnorm(n_anchor * 3), n_anchor, 3,
           dimnames = list(NULL, covariates))
  )
  anchor_block$T <- "A"

  distinct <- as.data.frame(
    matrix(stats::rnorm(n_distinct * 3), n_distinct, 3,
           dimnames = list(NULL, covariates))
  )
  # Replicate whole rows, so every copy is identical in all three covariates.
  comparator <- distinct[rep(seq_len(n_distinct), each = copies), ,
                         drop = FALSE]
  # Each distinct point contributes `copies / 2` duplicates to each group.
  comparator$T <- rep(rep(c("B", "C"), each = copies / 2), times = n_distinct)

  data <- rbind(anchor_block, comparator)
  rownames(data) <- NULL
  data
}

test_that("the tied fixture really does produce duplicate GPS rows", {
  # Guards the test below: if the fixture stopped producing exact ties, the
  # contested branch would never run and the comparison would prove nothing.
  data <- tied_candidate_data()
  gps <- estimate_gps_multinom(data, X_vars = c("x1", "x2", "x3"),
                               treatment_var = "T", anchor_level = "A")$gps

  comparator_gps <- gps[data$T != "A", , drop = FALSE]
  expect_lt(nrow(unique(comparator_gps)), nrow(comparator_gps) / 2)
})

test_that("exact ties are broken by row order, not by tree order", {
  data <- tied_candidate_data()
  gps <- estimate_gps_multinom(data, X_vars = c("x1", "x2", "x3"),
                               treatment_var = "T", anchor_level = "A")$gps

  # The tie-resolution branch must actually be reached, otherwise this test
  # only re-checks the ordinary path.
  labels <- treatment_labels(data, "T")
  X_anchor <- gps[labels == "A", , drop = FALSE]
  X_group <- gps[labels == "B", , drop = FALSE]
  n_contested <- vapply(c(1L, 2L, 3L, 7L, 12L), function(top_m) {
    m <- min(top_m, nrow(X_group))
    k_query <- min(m + 1L, nrow(X_group))
    if (k_query == m) return(0L)
    neighbours <- RANN::nn2(data = X_group, query = X_anchor, k = k_query)
    sum(neighbours$nn.dists[, k_query] <=
          neighbours$nn.dists[, m] * (1 + .SAM_TIE_TOL))
  }, integer(1))
  expect_true(any(n_contested > 0L))

  for (gps_space in c("raw", "logit")) {
    for (top_m in c(1L, 2L, 3L, 7L, 12L)) {
      expect_matches_reference(data, gps, top_m, gps_space)
    }
  }
})

test_that("matching on tied candidate pools still uses each subject once", {
  data <- tied_candidate_data()
  X_vars <- c("x1", "x2", "x3")
  gps <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "A")$gps
  search <- gps_candidate_search(data, gps, treatment_var = "T",
                                 anchor_level = "A", top_m = 5)
  result <- sam_match(data, search, X_vars = X_vars, treatment_var = "T")

  for (g in search$groups) {
    expect_equal(anyDuplicated(result$matched[[g]]), 0L)
  }
})

test_that("identical anchors receive identical candidate pools", {
  set.seed(99L)
  data <- sam_small_data(n = 120L, n_groups = 3L)
  # Make the first two anchors exact duplicates of each other.
  anchor_rows <- which(data$T == "A")
  data[anchor_rows[2], c("x1", "x2", "x3")] <-
    data[anchor_rows[1], c("x1", "x2", "x3")]

  gps <- estimate_gps_multinom(data, X_vars = c("x1", "x2", "x3"),
                               treatment_var = "T", anchor_level = "A")$gps
  search <- gps_candidate_search(data, gps, treatment_var = "T",
                                 anchor_level = "A", top_m = 6)

  expect_equal(search$candidates[[1]], search$candidates[[2]])
})

test_that("candidates always belong to their comparator group", {
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")

  for (g in fixture$search$groups) {
    rows <- unique(unlist(lapply(fixture$search$candidates, `[[`, g)))
    expect_true(all(labels[rows] == g))
  }
})

test_that("a group smaller than top_m returns the whole group", {
  set.seed(7L)
  data <- data.frame(x1 = stats::rnorm(40), x2 = stats::rnorm(40))
  data$T <- c(rep("A", 36), rep("B", 4))

  gps <- estimate_gps_multinom(data, X_vars = c("x1", "x2"),
                               treatment_var = "T", anchor_level = "A")$gps
  search <- gps_candidate_search(data, gps, treatment_var = "T",
                                 anchor_level = "A", top_m = 10)

  group_rows <- which(data$T == "B")
  for (candidates in search$candidates) {
    expect_setequal(candidates$B, group_rows)
  }
})
