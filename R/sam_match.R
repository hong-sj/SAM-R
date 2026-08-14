################################################################################
### Shared Anchor Matching
##
## Author: Sungjun Hong
################################################################################

#' Match multiple treatment groups to a shared anchor
#'
#' @description
#' Performs global greedy matching of multiple treatment groups to a shared
#' anchor using GPS-screened candidate pools and Mahalanobis distance.
#'
#' For each anchor, the nearest available candidate from each comparator
#' group is combined into a matched set. At each iteration, the anchor with
#' the smallest total Mahalanobis distance is selected, and all subjects in
#' that matched set are removed from further matching.
#'
#' @param data A `data.frame` containing covariate and treatment variables.
#' @param search Output from [gps_candidate_search()].
#' @param X_vars Character vector of covariate column names used for
#'   Mahalanobis distance. Default `paste0("X", 1:10)`.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{matched}{A `data.frame` containing successfully matched sets,
#'       subject row indices, group-specific Mahalanobis distances, and
#'       total matched-set loss.}
#'     \item{unmatched_anchor_rows}{Row indices of anchor subjects that
#'       could not be matched.}
#'     \item{matching_rate}{Proportion of anchor subjects successfully
#'       matched.}
#'     \item{max_possible_rate}{The highest rate the group sizes allow. Each
#'       matched set consumes one subject from every comparator group, so the
#'       rate is capped by the smallest comparator group over the anchor
#'       count. Reported alongside `matching_rate` so that a rate which looks
#'       low can be recognised as saturated.}
#'   }
#'
#' @examples
#' data(sample_4group)
#'
#' covariates <- setdiff(
#'   names(sample_4group),
#'   c("synthetic_id", "treatment", "mortality_28d")
#' )
#' anchor <- unique(as.character(sample_4group$treatment))[1]
#'
#' fit <- estimate_gps_multinom(
#'   sample_4group,
#'   X_vars = covariates,
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#'
#' search <- gps_candidate_search(
#'   sample_4group,
#'   fit$gps,
#'   treatment_var = "treatment",
#'   anchor_level = anchor,
#'   top_m = 10,
#'   gps_space = "logit"
#' )
#'
#' matched <- sam_match(
#'   sample_4group,
#'   search,
#'   X_vars = covariates,
#'   treatment_var = "treatment"
#' )
#'
#' matched$matching_rate
#' head(matched$matched)
#'
#' @export
sam_match <- function(data, search, X_vars = paste0("X", 1:10),
                      treatment_var = "T") {
  
  # Validate the covariates first. A caller who both rewrites a value and makes
  # it unusable -- inserting an NA, say -- should be told about the NA, not
  # about the fingerprint mismatch the same edit produces.
  X <- covariate_matrix(data, X_vars)

  check_fingerprint(search, data, treatment_var)
  labels <- treatment_labels(data, treatment_var)
  anchor_rows <- as.integer(search$anchor_rows)
  groups <- as.character(search$groups)
  candidates <- search$candidates
  n_A <- length(anchor_rows)
  
  # Empty anchor group
  if (n_A == 0L) {
    return(c(
      list(matched = empty_matched_frame(groups)),
      summarize_matching(empty_matched_frame(groups), anchor_rows, list())
    ))
  }

  # Mahalanobis distance setup
  pooled <- get_pooled_covariance(data, X_vars, treatment_var)
  S_inv <- pooled$S_inv
  
  group_rows <- stats::setNames(
    lapply(groups, function(g) {
      require_rows(which(labels == g), g, treatment_var, available = labels)
    }),
    groups
  )
  
  # The covariates were materialised once above and are sliced here, rather
  # than rebuilding a frame per group.
  X_anchor <- X[anchor_rows, , drop = FALSE]

  X_group <- stats::setNames(
    lapply(groups, function(g) X[group_rows[[g]], , drop = FALSE]),
    groups
  )
  
  # --- Candidate distances
  # Compute Mahalanobis distances for GPS-screened candidates and sort
  # candidates by distance for each anchor and comparator group
  candidate_sorted <- stats::setNames(
    vector("list", length(groups)),
    groups
  )
  
  distance_sorted <- stats::setNames(
    vector("list", length(groups)),
    groups
  )
  
  for (g in groups) {
    candidate_sorted[[g]] <- vector("list", n_A)
    distance_sorted[[g]] <- vector("list", n_A)
    
    for (i in seq_len(n_A)) {
      # `match()` returns NA for a candidate that is not in this group; drop
      # those directly rather than through na.omit(), which builds an
      # attribute recording every removed position and dominated this loop.
      idx <- match(candidates[[i]][[g]], group_rows[[g]])
      idx <- as.integer(idx[!is.na(idx)])

      if (length(idx) == 0L) {
        candidate_sorted[[g]][[i]] <- integer(0)
        distance_sorted[[g]][[i]] <- numeric(0)
        next
      }
      
      d <- as.numeric(
        mahalanobis_distance_matrix(
          X_anchor[i, , drop = FALSE],
          X_group[[g]][idx, , drop = FALSE],
          S_inv
        )
      )
      
      # Preserve candidate order when distances are tied
      ord <- order(d, seq_along(d))
      
      candidate_sorted[[g]][[i]] <- idx[ord]
      distance_sorted[[g]][[i]] <- d[ord]
    }
  }
  
  # Candidate availability
  active <- stats::setNames(
    lapply(groups, function(g) rep(TRUE, length(group_rows[[g]]))),
    groups
  )
  
  pointer <- stats::setNames(
    lapply(groups, function(g) rep(1L, n_A)),
    groups
  )
  
  # `alive` replaces a shrinking `remaining_anchor` vector, which was
  # reallocated on every match. A dead anchor carries best_loss = Inf, so the
  # global which.min below never selects one.
  alive <- rep(TRUE, n_A)
  needs_update <- rep(TRUE, n_A)
  best_loss <- rep(Inf, n_A)


  best_choice <- matrix(
    NA_integer_,
    nrow = n_A,
    ncol = length(groups),
    dimnames = list(NULL, groups)
  )
  
  best_dist <- matrix(
    NA_real_,
    nrow = n_A,
    ncol = length(groups),
    dimnames = list(NULL, groups)
  )

  # --- Reverse index
  # Which anchors are currently claiming which comparator subject, as one
  # doubly-linked list per comparator subject. Consuming a subject then visits
  # only its own claimants, instead of rescanning every remaining anchor after
  # every match. The lists must be doubly linked: an anchor moves between
  # buckets as its pick is consumed, and a singly-linked list cannot unlink a
  # node from the middle -- an anchor that re-picks the same subject would link
  # to itself and the walk would never terminate.
  claim <- matrix(0L, nrow = n_A, ncol = length(groups))
  link_next <- matrix(0L, nrow = n_A, ncol = length(groups))
  link_prev <- matrix(0L, nrow = n_A, ncol = length(groups))
  bucket_head <- lapply(groups, function(g) integer(length(group_rows[[g]])))

  unlink_claim <- function(i, gi) {
    held <- claim[i, gi]
    if (held == 0L) {
      return(invisible(NULL))
    }

    previous <- link_prev[i, gi]
    following <- link_next[i, gi]

    if (previous == 0L) {
      bucket_head[[gi]][held] <<- following
    } else {
      link_next[previous, gi] <<- following
    }
    if (following != 0L) {
      link_prev[following, gi] <<- previous
    }

    claim[i, gi] <<- 0L
    link_next[i, gi] <<- 0L
    link_prev[i, gi] <<- 0L
    invisible(NULL)
  }

  set_claim <- function(i, gi, choice) {
    if (claim[i, gi] == choice) {
      return(invisible(NULL))
    }
    unlink_claim(i, gi)

    head <- bucket_head[[gi]][choice]
    link_next[i, gi] <<- head
    link_prev[i, gi] <<- 0L
    if (head != 0L) {
      link_prev[head, gi] <<- i
    }
    bucket_head[[gi]][choice] <<- i
    claim[i, gi] <<- choice
    invisible(NULL)
  }

  # Best available matched set for one anchor
  recompute <- function(i) {
    total <- 0

    for (gi in seq_along(groups)) {
      g <- groups[gi]
      idx <- candidate_sorted[[g]][[i]]
      d <- distance_sorted[[g]][[i]]
      p <- pointer[[g]][i]

      while (p <= length(idx) && !active[[g]][idx[p]]) {
        p <- p + 1L
      }

      pointer[[g]][i] <<- p

      if (p > length(idx)) {
        best_loss[i] <<- Inf
        best_choice[i, ] <<- NA_integer_
        best_dist[i, ] <<- NA_real_
        for (gj in seq_along(groups)) {
          unlink_claim(i, gj)
        }
        return(invisible(NULL))
      }

      best_choice[i, gi] <<- idx[p]
      best_dist[i, gi] <<- d[p]
      set_claim(i, gi, idx[p])
      total <- total + d[p]
    }

    best_loss[i] <<- total
    invisible(NULL)
  }
  
  # Global greedy matching.
  # Results accumulate into preallocated vectors and are assembled once at the
  # end; growing a list of one-row data frames and rbind-ing it is quadratic in
  # the number of matched sets.
  n_matched <- 0L
  matched_anchor <- integer(n_A)
  matched_loss <- numeric(n_A)
  matched_choice <- matrix(NA_integer_, n_A, length(groups),
                           dimnames = list(NULL, groups))
  matched_dist <- matrix(NA_real_, n_A, length(groups),
                         dimnames = list(NULL, groups))

  stale <- which(alive)

  repeat {
    if (length(stale) > 0L) {
      for (i in stale) {
        recompute(i)
      }
      needs_update[stale] <- FALSE
    }

    # Select the globally smallest-loss matched set. Retired anchors, and
    # anchors with no available candidate in some group, both carry Inf, so a
    # non-finite minimum means nothing matchable is left. Scanning the full
    # vector beats maintaining a pure-R heap: which.min is one C-level pass
    # with no allocation, whereas a heap would pay R-level call overhead per
    # sift.
    i_star <- which.min(best_loss)
    if (!is.finite(best_loss[i_star])) {
      break
    }

    choice_star <- best_choice[i_star, , drop = TRUE]

    n_matched <- n_matched + 1L
    matched_anchor[n_matched] <- anchor_rows[i_star]
    matched_loss[n_matched] <- best_loss[i_star]
    matched_dist[n_matched, ] <- best_dist[i_star, ]

    for (gi in seq_along(groups)) {
      matched_choice[n_matched, gi] <-
        group_rows[[gi]][choice_star[[gi]]]
    }

    # Retire the selected anchor.
    alive[i_star] <- FALSE
    best_loss[i_star] <- Inf
    for (gi in seq_along(groups)) {
      unlink_claim(i_star, gi)
    }

    # Consume the selected comparator subjects and visit only the anchors that
    # were claiming them.
    affected <- integer(0)

    for (gi in seq_along(groups)) {
      chosen <- choice_star[[gi]]
      active[[gi]][chosen] <- FALSE

      claimant <- bucket_head[[gi]][chosen]
      while (claimant != 0L) {
        following <- link_next[claimant, gi]
        unlink_claim(claimant, gi)
        affected <- c(affected, claimant)
        claimant <- following
      }
    }

    if (length(affected) > 0L) {
      affected <- unique(affected)
      affected <- affected[alive[affected]]
      needs_update[affected] <- TRUE
    }
    stale <- affected
  }

  # Output
  keep <- seq_len(n_matched)
  named_groups <- stats::setNames(groups, groups)

  matched <- build_matched_frame(
    anchor = matched_anchor[keep],
    group_choices = lapply(named_groups, function(g) matched_choice[keep, g]),
    group_distances = lapply(named_groups, function(g) matched_dist[keep, g]),
    loss = matched_loss[keep],
    groups = groups
  )

  c(list(matched = matched),
    summarize_matching(matched, anchor_rows, group_rows))
}