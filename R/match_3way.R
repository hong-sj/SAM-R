################################################################################
### Three-Way Propensity Score Matching
##
## Author: Sungjun Hong
################################################################################

# Pure-R implementation of the three-way nearest-neighbor propensity score
# matching method described by Rassen et al. (2013) and used by
# Yoshida et al. (2017).

#' Calculate a caliper for three-way matching
#'
#' @description
#' Computes a data-driven caliper from the variance of the two propensity
#' score dimensions used for three-way matching.
#'
#' @param ps_used Numeric matrix with two propensity score columns.
#' @param treatment_var_values Treatment-group labels for each subject.
#'
#' @return A numeric caliper value.
#'
#' @export
calc_caliper_3way <- function(ps_used, treatment_var_values) {
  
  stopifnot(ncol(ps_used) == 2, nrow(ps_used) == length(treatment_var_values))
  var_by_group <- sapply(seq_len(2), function(j) {
    tapply(ps_used[, j], treatment_var_values, stats::var)
  })
  
  0.6 * sqrt(sum(rowMeans(var_by_group)) / 3)
}

# Build a 2D KD-tree
kdtree_build <- function(coords, idx = seq_len(nrow(coords)), depth = 0L) {
  n <- length(idx)
  if (n == 1L) {
    return(list(leaf = TRUE, point = idx[1]))
  }
  dim <- (depth %% 2L) + 1L
  ord <- idx[order(coords[idx, dim])]
  mid <- n %/% 2L
  left_idx  <- ord[seq_len(mid)]
  right_idx <- ord[(mid + 1L):n]
  split_val <- coords[right_idx[1], dim]
  list(
    leaf = FALSE, split_dim = dim, split_val = split_val,
    left  = kdtree_build(coords, left_idx,  depth + 1L),
    right = kdtree_build(coords, right_idx, depth + 1L)
  )
}

# Find the nearest active point in a KD-tree
kdtree_nearest <- function(node, coords, query, active) {
  if (node$leaf) {
    if (!active[node$point]) return(NULL)
    return(node$point)
  }
  gap <- query[node$split_dim] - node$split_val
  if (gap <= 0) { primary <- node$left; other <- node$right
  } else        { primary <- node$right; other <- node$left }

  best <- kdtree_nearest(primary, coords, query, active)
  best_d2 <- if (!is.null(best)) sum((coords[best, ] - query)^2) else Inf

  # Check the opposite branch when it could contain a closer point
  if (gap * gap < best_d2) {
    cand <- kdtree_nearest(other, coords, query, active)
    if (!is.null(cand)) {
      cand_d2 <- sum((coords[cand, ] - query)^2)
      if (cand_d2 < best_d2) { best <- cand; best_d2 <- cand_d2 }
    }
  }
  best
}

# Find active points within a given squared radius
kdtree_range <- function(node, coords, query, radius2, active) {
  if (node$leaf) {
    if (!active[node$point]) return(integer(0))
    d2 <- sum((coords[node$point, ] - query)^2)

    if (d2 < radius2) return(node$point) else return(integer(0))
  }
  gap <- query[node$split_dim] - node$split_val
  if (gap <= 0) { primary <- node$left; other <- node$right
  } else        { primary <- node$right; other <- node$left }

  result <- kdtree_range(primary, coords, query, radius2, active)
  if (gap * gap < radius2) {
    result <- c(result, kdtree_range(other, coords, query, radius2, active))
  }
  result
}

#' Three-way propensity score matching
#'
#' @description
#' Performs 1:1:1 nearest-neighbor propensity score matching for three
#' treatment groups using the Rassen three-way matching approach.
#'
#' The smallest treatment group is used as the matching base, and candidate
#' trios are selected using distances in a two-dimensional propensity score
#' space subject to a caliper.
#'
#' @param data A `data.frame` containing the treatment variable.
#' @param search Output from [gps_candidate_search()], used to identify the
#'   anchor and comparator treatment groups.
#' @param gps Numeric matrix of generalized propensity scores.
#' @param X_vars Character vector of covariate names. Retained for
#'   compatibility with [sam_match()]. Default `paste0("X", 1:10)`.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#' @param caliper Numeric caliper or `"auto"` to estimate it from the data.
#'   Default `"auto"`.
#' @param ps_space Propensity score scale used for matching, either `"raw"`
#'   or `"logit"`. Default `"raw"`.
#' @param top_n Maximum number of candidate trios retained per subject.
#'   Must be a positive integer. Default `10`.
#' @param reference_level Treatment group used as the propensity score
#'   reference. If `NULL`, the last treatment level is used.
#'
#' @return A list containing:
#'   \describe{
#'     \item{matched}{A `data.frame` containing the matched trios,
#'       propensity score distances, `loss`, and `rassen_perimeter`.}
#'     \item{unmatched_anchor_rows}{Row indices of unmatched anchor subjects.}
#'     \item{matching_rate}{Proportion of anchor subjects included in
#'       matched trios.}
#'     \item{caliper}{Caliper value used for matching.}
#'     \item{red_level}{Smallest treatment group used as the matching base.}
#'     \item{reference_level}{Treatment group used as the propensity score
#'       reference.}
#'   }
#'
#' @details
#' The `loss` returned by `match_3way()` is based on propensity score
#' distances and is not directly comparable with the Mahalanobis-distance
#' loss returned by [sam_match()]. Covariate balance measures and matching
#' rates can be compared across methods.
#'
#' @references
#' Rassen JA, Shelat AA, Franklin JM, Solomon DH, Bray M, Brookhart MA.
#' Matching by propensity score in cohort studies with three treatment
#' groups. Epidemiology. 2013.
#'
#' Yoshida K, et al. Matching weights to simultaneously compare three
#' treatment groups. Epidemiology. 2017.
#'
#' @examples
#' data(sample_3group)
#'
#' covariates <- setdiff(
#'   names(sample_3group),
#'   c("synthetic_id", "treatment", "mortality_28d")
#' )
#' anchor <- unique(as.character(sample_3group$treatment))[1]
#'
#' fit <- estimate_gps_multinom(
#'   sample_3group,
#'   X_vars = covariates,
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#'
#' search <- gps_candidate_search(
#'   sample_3group,
#'   fit$gps,
#'   treatment_var = "treatment",
#'   anchor_level = anchor,
#'   top_m = 10,
#'   gps_space = "logit"
#' )
#'
#' matched <- match_3way(
#'   sample_3group,
#'   search,
#'   fit$gps,
#'   treatment_var = "treatment"
#' )
#'
#' matched$matching_rate
#' head(matched$matched)
#'
#' @export
match_3way <- function(data, search, gps, X_vars = paste0("X", 1:10),
                        treatment_var = "T", caliper = "auto",
                        ps_space = c("raw", "logit"), top_n = 10L,
                        reference_level = NULL) {
  ps_space <- match.arg(ps_space)
  top_n <- require_positive_int(top_n, "top_n")
  check_fingerprint(search, data, treatment_var)
  labels <- treatment_labels(data, treatment_var)
  groups <- as.character(search$groups)
  if (length(groups) != 2L) {
    stop("match_3way() requires exactly three treatment groups ",
         "(1 anchor + 2 comparator groups).")
  }
  anchor_rows <- search$anchor_rows

  anchor_level <- unique(labels[anchor_rows])
  stopifnot(length(anchor_level) == 1)
  all_levels <- c(anchor_level, groups)
  stopifnot(all(all_levels %in% colnames(gps)))
  stopifnot(nrow(gps) == nrow(data))

  # --- Propensity score space
  # Define the two-dimensional propensity score space
  if (is.null(reference_level)) {
    reference_level <- all_levels[length(all_levels)]
  }
  stopifnot(reference_level %in% all_levels)
  ps_levels <- setdiff(all_levels, reference_level)
  stopifnot(length(ps_levels) == 2L)

  ps_raw <- gps[, ps_levels, drop = FALSE]
  if (ps_space == "logit") {
    eps <- 1e-6
    ps_used <- stats::qlogis(pmin(pmax(ps_raw, eps), 1 - eps))
  } else {
    ps_used <- ps_raw
  }
  colnames(ps_used) <- ps_levels

  # --- Caliper
  if (identical(caliper, "auto")) {
    caliper <- calc_caliper_3way(ps_used, labels)
  }
  stopifnot(is.numeric(caliper), length(caliper) == 1, caliper > 0)

  # --- Treatment groups
  # Use the smallest treatment group as the matching base
  rows_by_level <- stats::setNames(
    lapply(all_levels, function(lv) {
      require_rows(which(labels == lv), lv, treatment_var, available = labels)
    }),
    all_levels
  )
  sizes <- vapply(rows_by_level, length, integer(1))
  red_level    <- names(sizes)[which.min(sizes)]
  other_levels <- setdiff(all_levels, red_level)
  o1 <- other_levels[1]; o2 <- other_levels[2]
  n_red <- sizes[[red_level]]

  coords_red <- ps_used[rows_by_level[[red_level]], , drop = FALSE]
  coords_o1  <- ps_used[rows_by_level[[o1]],        , drop = FALSE]
  coords_o2  <- ps_used[rows_by_level[[o2]],        , drop = FALSE]

  # --- KD-trees and candidate tracking
  # Use the smallest treatment group as the matching base
  tree_o1 <- kdtree_build(coords_o1)
  tree_o2 <- kdtree_build(coords_o2)
  active_o1 <- rep(TRUE, nrow(coords_o1))
  active_o2 <- rep(TRUE, nrow(coords_o2))

  matched_flag <- rep(FALSE, n_red)
  exhausted    <- rep(FALSE, n_red)
  counters     <- integer(n_red)

  # Shared pool of candidate trios
  pool_red <- integer(0); pool_o1 <- integer(0); pool_o2 <- integer(0)
  pool_perim <- numeric(0); pool_d_o1 <- numeric(0); pool_d_o2 <- numeric(0)

  # --- Candidate generation
  # Build KD-trees for the two remaining treatment groups
  push_candidates <- function(i) {
    pr <- coords_red[i, ]
    nb <- kdtree_nearest(tree_o1, coords_o1, pr, active_o1)
    if (is.null(nb)) { exhausted[i] <<- TRUE; counters[i] <<- 0L; return(invisible(NULL)) }
    nbg <- kdtree_nearest(tree_o2, coords_o2, coords_o1[nb, ], active_o2)
    if (is.null(nbg)) { exhausted[i] <<- TRUE; counters[i] <<- 0L; return(invisible(NULL)) }

    # Initial candidate triangle
    d_pr_nb  <- sqrt(sum((pr - coords_o1[nb, ])^2))
    d_pr_nbg <- sqrt(sum((pr - coords_o2[nbg, ])^2))
    d_nb_nbg <- sqrt(sum((coords_o1[nb, ] - coords_o2[nbg, ])^2))
    small <- d_pr_nb + d_pr_nbg + d_nb_nbg

    cand_o1 <- integer(0); cand_o2 <- integer(0); cand_perim <- numeric(0)
    cand_d1 <- numeric(0); cand_d2 <- numeric(0)

    # Include the initial candidate when it satisfies the caliper
    if (small <= caliper) {
      cand_o1 <- nb; cand_o2 <- nbg; cand_perim <- small
      cand_d1 <- d_pr_nb; cand_d2 <- d_pr_nbg
    }

    # Search for additional candidates near the base subject
    radius2 <- (small / 2)^2
    nbs <- kdtree_range(tree_o1, coords_o1, pr, radius2, active_o1)
    ngs <- kdtree_range(tree_o2, coords_o2, pr, radius2, active_o2)

    if (length(nbs) > 0 && length(ngs) > 0) {
      A <- coords_o1[nbs, , drop = FALSE]; B <- coords_o2[ngs, , drop = FALSE]
      d_pr_A <- sqrt(rowSums((A - matrix(pr, nrow(A), 2, byrow = TRUE))^2))
      d_pr_B <- sqrt(rowSums((B - matrix(pr, nrow(B), 2, byrow = TRUE))^2))
      a_sq <- rowSums(A * A); b_sq <- rowSums(B * B)
      d2m <- outer(a_sq, b_sq, "+") - 2 * (A %*% t(B))
      
      # Guard against small negative values from floating-point error
      d2m[d2m < 0] <- 0
      d_AB <- sqrt(d2m)
      perim <- outer(d_pr_A, d_pr_B, "+") + d_AB
      keep <- which(perim < small & perim <= caliper, arr.ind = TRUE)
      if (nrow(keep) > 0) {
        cand_o1    <- c(cand_o1,    nbs[keep[, 1]])
        cand_o2    <- c(cand_o2,    ngs[keep[, 2]])
        cand_perim <- c(cand_perim, perim[keep])
        cand_d1    <- c(cand_d1,    d_pr_A[keep[, 1]])
        cand_d2    <- c(cand_d2,    d_pr_B[keep[, 2]])
      }
    }

    if (length(cand_perim) == 0) {
      exhausted[i] <<- TRUE; counters[i] <<- 0L
      return(invisible(NULL))
    }

    ord <- order(cand_perim)[seq_len(min(top_n, length(cand_perim)))]
    pool_red   <<- c(pool_red,   rep(i, length(ord)))
    pool_o1    <<- c(pool_o1,    cand_o1[ord])
    pool_o2    <<- c(pool_o2,    cand_o2[ord])
    pool_perim <<- c(pool_perim, cand_perim[ord])
    pool_d_o1  <<- c(pool_d_o1,  cand_d1[ord])
    pool_d_o2  <<- c(pool_d_o2,  cand_d2[ord])
    counters[i] <<- length(ord)
    invisible(NULL)
  }

  # Initialize candidate trios
  for (i in seq_len(n_red)) push_candidates(i)

  # --- Global greedy matching
  matched_rows <- list()

  while (length(pool_perim) > 0) {
    pop <- which.min(pool_perim)
    i <- pool_red[pop]; b <- pool_o1[pop]; g <- pool_o2[pop]
    perim_popped <- pool_perim[pop]

    pool_red   <- pool_red[-pop];   pool_o1   <- pool_o1[-pop];   pool_o2 <- pool_o2[-pop]
    pool_perim <- pool_perim[-pop]; pool_d_o1 <- pool_d_o1[-pop]; pool_d_o2 <- pool_d_o2[-pop]

    # Skip candidates belonging to an already matched base subject
    if (matched_flag[i]) next

    if (active_o1[b] && active_o2[g]) {
      # Accept the candidate trio
      matched_flag[i] <- TRUE
      active_o1[b] <- FALSE
      active_o2[g] <- FALSE

      row_red_global <- rows_by_level[[red_level]][i]
      row_o1_global  <- rows_by_level[[o1]][b]
      row_o2_global  <- rows_by_level[[o2]][g]

      # Relabel subjects according to the anchor and comparator groups
      rows_this_trio <- stats::setNames(
        list(row_red_global, row_o1_global, row_o2_global),
        c(red_level, o1, o2)
      )
      row_anchor <- rows_this_trio[[anchor_level]]

      row_df <- data.frame(matched_set_id = length(matched_rows) + 1L, anchor = row_anchor)
      
      for (gg in groups) row_df[[gg]] <- rows_this_trio[[gg]]

      # Compute anchor-to-comparator propensity score distances
      for (gg in groups) {
        row_df[[paste0("dist_", gg)]] <- sqrt(sum((ps_used[row_anchor, ] - ps_used[row_df[[gg]], ])^2))
      }
      
      # Compute anchor-to-comparator propensity score distances
      row_df$loss <- sum(vapply(groups, function(gg) row_df[[paste0("dist_", gg)]], numeric(1)))
      
      # Full propensity-score triangle perimeter used by this method
      row_df$rassen_perimeter <- perim_popped
      matched_rows[[length(matched_rows) + 1]] <- row_df
    } else {
      # Refresh candidates when all current options become unavailable
      counters[i] <- counters[i] - 1L
      if (counters[i] <= 0L && !exhausted[i]) {
        push_candidates(i)
      }
    }
  }

  # --- Output
  matched <- if (length(matched_rows) > 0) {
    out <- do.call(rbind, matched_rows)
    rownames(out) <- NULL
    out
  } else {
    empty <- data.frame(matched_set_id = integer(0), anchor = integer(0))
    for (g in groups) empty[[g]] <- integer(0)
    for (g in groups) empty[[paste0("dist_", g)]] <- numeric(0)
    empty$loss <- numeric(0)
    empty$rassen_perimeter <- numeric(0)
    empty
  }

  matched_anchor_rows <- if (nrow(matched) > 0) matched$anchor else integer(0)
  unmatched_anchor_rows <- setdiff(anchor_rows, matched_anchor_rows)

  list(
    matched = matched,
    unmatched_anchor_rows = unmatched_anchor_rows,
    matching_rate = nrow(matched) / length(anchor_rows),
    caliper = caliper,
    red_level = red_level,
    reference_level = reference_level
  )
}
