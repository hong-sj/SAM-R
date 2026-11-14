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

  if (!all(is.finite(ps_used))) {
    stop("`ps_used` must contain only finite values.", call. = FALSE)
  }

  groups <- unique(treatment_var_values)
  if (length(groups) != 3L) {
    stop("`treatment_var_values` must contain exactly three treatment groups, ",
         "got ", length(groups), ".", call. = FALSE)
  }

  # stats::var() returns NA for a single observation, which would propagate
  # into the caliper as a silent NA rather than a refusal.
  sizes <- table(treatment_var_values)
  too_small <- names(sizes)[sizes < 2L]
  if (length(too_small) > 0L) {
    stop("The automatic three-way caliper needs at least two subjects in ",
         "every treatment group; too few in: ",
         paste(too_small, collapse = ", "), ".", call. = FALSE)
  }

  var_by_group <- sapply(seq_len(2), function(j) {
    tapply(ps_used[, j], treatment_var_values, stats::var)
  })

  0.6 * sqrt(sum(rowMeans(var_by_group)) / 3)
}

# --- 2-D k-d tree with removable points
#
# Points are removed from the search as they are matched. Flagging them
# inactive and letting the walk find out at the leaves means that once most of
# a group is consumed, every query descends and backtracks through large dead
# subtrees to reach the few points still available -- the walk grows with the
# number of removals rather than with log(n).
#
# Each node therefore carries a count of the active points beneath it, which a
# removal decrements along the path from the leaf to the root. A subtree whose
# count has reached zero is skipped outright.
#
# The tree is held in an environment of flat vectors, so a removal updates the
# counts in place; a nested list would copy the path on every write. Traversal
# order is unchanged, which matters because callers rely on the order
# kdtree_range() returns points in, and on nearest ties resolving in favour of
# the primary branch.

kdtree_build <- function(coords) {
  # The splits alternate between two dimensions, so anything else would index
  # a column that is not there or silently ignore one. An empty group would
  # leave the recursion below with nothing to split.
  if (!is.matrix(coords) || ncol(coords) != 2L) {
    stop("`coords` must be a two-column matrix, got ",
         if (is.matrix(coords)) paste0(ncol(coords), " columns") else
           class(coords)[1], ".", call. = FALSE)
  }
  if (nrow(coords) == 0L) {
    stop("`coords` must contain at least one point.", call. = FALSE)
  }
  if (!all(is.finite(coords))) {
    stop("`coords` must contain only finite values.", call. = FALSE)
  }

  n <- nrow(coords)
  capacity <- max(1L, 2L * n)

  tree <- new.env(parent = emptyenv())
  tree$leaf <- logical(capacity)
  tree$point <- integer(capacity)
  tree$split_dim <- integer(capacity)
  tree$split_val <- numeric(capacity)
  tree$left <- integer(capacity)
  tree$right <- integer(capacity)
  tree$parent <- integer(capacity)
  tree$active <- integer(capacity)
  tree$point_leaf <- integer(n)

  next_id <- 0L

  build <- function(idx, depth, parent) {
    next_id <<- next_id + 1L
    id <- next_id
    tree$parent[id] <- parent
    tree$active[id] <- length(idx)

    if (length(idx) == 1L) {
      tree$leaf[id] <- TRUE
      tree$point[id] <- idx[1L]
      tree$point_leaf[idx[1L]] <- id
      return(id)
    }

    dim <- (depth %% 2L) + 1L
    ord <- idx[order(coords[idx, dim])]
    mid <- length(idx) %/% 2L
    left_idx <- ord[seq_len(mid)]
    right_idx <- ord[(mid + 1L):length(idx)]

    tree$leaf[id] <- FALSE
    tree$split_dim[id] <- dim
    tree$split_val[id] <- coords[right_idx[1L], dim]
    tree$left[id] <- build(left_idx, depth + 1L, id)
    tree$right[id] <- build(right_idx, depth + 1L, id)
    id
  }

  tree$root <- build(seq_len(n), 0L, 0L)
  tree
}

# Is a point still available?
kdtree_active <- function(tree, point) {
  tree$active[tree$point_leaf[point]] > 0L
}

# Remove a point, decrementing every ancestor's count. The path is collected
# first and updated in one assignment: taking a second reference to
# `tree$active` and writing through it would copy the whole vector.
kdtree_remove <- function(tree, point) {
  parent <- tree$parent
  id <- tree$point_leaf[point]

  path <- integer(64L)
  depth <- 0L
  while (id != 0L) {
    depth <- depth + 1L
    path[depth] <- id
    id <- parent[id]
  }

  path <- path[seq_len(depth)]
  tree$active[path] <- tree$active[path] - 1L
  invisible(NULL)
}

# Nearest active point to `query`, or NULL when the tree is empty.
kdtree_nearest <- function(tree, coords, query) {
  leaf <- tree$leaf; point <- tree$point
  left <- tree$left; right <- tree$right
  split_dim <- tree$split_dim; split_val <- tree$split_val

  descend <- function(id) {
    if (tree$active[id] == 0L) return(NULL)
    if (leaf[id]) return(point[id])

    gap <- query[split_dim[id]] - split_val[id]
    if (gap <= 0) { primary <- left[id]; other <- right[id]
    } else        { primary <- right[id]; other <- left[id] }

    best <- descend(primary)
    best_d2 <- if (!is.null(best)) sum((coords[best, ] - query)^2) else Inf

    # Check the opposite branch when it could contain a closer point
    if (gap * gap < best_d2) {
      cand <- descend(other)
      if (!is.null(cand)) {
        cand_d2 <- sum((coords[cand, ] - query)^2)
        if (cand_d2 < best_d2) { best <- cand; best_d2 <- cand_d2 }
      }
    }
    best
  }

  descend(tree$root)
}

# Active points strictly within a squared radius, in traversal order.
kdtree_range <- function(tree, coords, query, radius2) {
  leaf <- tree$leaf; point <- tree$point
  left <- tree$left; right <- tree$right
  split_dim <- tree$split_dim; split_val <- tree$split_val

  descend <- function(id) {
    if (tree$active[id] == 0L) return(integer(0))

    if (leaf[id]) {
      p <- point[id]
      d2 <- sum((coords[p, ] - query)^2)
      # Boundary points are excluded by design.
      if (d2 < radius2) return(p) else return(integer(0))
    }

    gap <- query[split_dim[id]] - split_val[id]
    if (gap <= 0) { primary <- left[id]; other <- right[id]
    } else        { primary <- right[id]; other <- left[id] }

    result <- descend(primary)
    if (gap * gap < radius2) {
      result <- c(result, descend(other))
    }
    result
  }

  descend(tree$root)
}

# --- Available-subject index for three-way matching
#
# Comparator subjects still available, sorted by their first coordinate. A
# query scans the vertical strip that can hold an answer -- located by binary
# search -- and decides membership on exact squared distances,
# dx * dx + dy * dy, the same arithmetic the Python implementation performs.
# Each query is a handful of vectorised operations instead of a recursive
# R-level tree walk, which was where three-way matching spent its time.
#
# Consumed subjects are flagged, and the sorted index is rebuilt over the
# survivors whenever half of the subjects it holds have been consumed, so a
# strip never fills up with dead entries.
.SAM_INDEX_RTOL <- 1e-9

.active_group <- function(coords) {
  group <- new.env(parent = emptyenv())
  # Plain unnamed vectors: findInterval() runs as.double() on its table, which
  # copies the whole vector on every query if it carries names.
  group$x <- unname(coords[, 1L])
  group$y <- unname(coords[, 2L])
  group$active <- rep(TRUE, nrow(coords))
  group$n_active <- nrow(coords)
  .ag_index(group)
  group
}

.ag_index <- function(group) {
  members <- which(group$active)
  ordered <- members[order(group$x[members], members)]
  group$members <- ordered
  group$sorted_x <- group$x[ordered]
  group$n_indexed <- length(ordered)
  # Starting half-width for nearest-neighbour strips: about one typical
  # spacing between points. Quadrupled until the answer is certain.
  if (length(members) > 0L) {
    extent <- max(diff(range(group$x[members])), diff(range(group$y[members])))
    group$start_width <- max(extent / sqrt(length(members)), 1e-8)
  }
  invisible(NULL)
}

.ag_deactivate <- function(group, point) {
  group$active[point] <- FALSE
  group$n_active <- group$n_active - 1L
  if (group$n_active * 2L <= group$n_indexed) {
    .ag_index(group)
  }
  invisible(NULL)
}

# Indexed subjects whose first coordinate lies within `half_width` of `qx`.
.ag_strip <- function(group, qx, half_width) {
  bounds <- findInterval(c(qx - half_width, qx + half_width), group$sorted_x,
                         checkSorted = FALSE, checkNA = FALSE)
  lo <- bounds[1L] + 1L
  hi <- bounds[2L]
  if (hi < lo) {
    return(integer(0))
  }
  group$members[lo:hi]
}

# Available subjects strictly inside the squared radius, by ascending row.
.ag_within <- function(group, query, radius2) {
  if (group$n_active == 0L) {
    return(integer(0))
  }
  points <- .ag_strip(group, query[1L],
                      sqrt(radius2) * (1 + .SAM_INDEX_RTOL) + 1e-300)
  points <- points[group$active[points]]
  dx <- group$x[points] - query[1L]
  dy <- group$y[points] - query[2L]
  sort(points[dx * dx + dy * dy < radius2])
}

# Lowest-row available subject at the minimum distance, or NA.
.ag_nearest <- function(group, query) {
  if (group$n_active == 0L) {
    return(NA_integer_)
  }
  qx <- query[[1L]]
  qy <- query[[2L]]
  width <- group$start_width
  n_indexed <- group$n_indexed
  repeat {
    bounds <- findInterval(c(qx - width, qx + width), group$sorted_x,
                           checkSorted = FALSE, checkNA = FALSE)
    lo <- bounds[1L] + 1L
    hi <- bounds[2L]
    whole <- lo <= 1L && hi >= n_indexed
    if (hi >= lo) {
      points <- group$members[lo:hi]
      points <- points[group$active[points]]
      if (length(points) > 0L) {
        dx <- group$x[points] - qx
        dy <- group$y[points] - qy
        d2 <- dx * dx + dy * dy
        best <- min(d2)
        # Everything outside the strip is at least `width` away. Once that
        # clears the best distance by the tolerance, nothing unseen can tie
        # or win.
        if (whole || sqrt(best) * (1 + .SAM_INDEX_RTOL) + 1e-300 < width) {
          return(min(points[d2 == best]))
        }
      } else if (whole) {
        return(NA_integer_)
      }
    } else if (whole) {
      return(NA_integer_)
    }
    width <- width * 4
  }
}

# Nearest subject for many queries while nothing has been consumed: one
# batched query to RANN's compiled tree, with any query whose runner-up is not
# clearly farther resolved exactly by .ag_nearest().
.nearest_untouched <- function(group, coords, queries) {
  k <- min(2L, nrow(coords))
  found <- RANN::nn2(data = coords, query = queries, k = k)
  result <- found$nn.idx[, 1L]
  if (k == 2L) {
    distance <- found$nn.dists
    unclear <- which(!(distance[, 2L] >
                         distance[, 1L] * (1 + .SAM_INDEX_RTOL) + 1e-300))
    for (j in unclear) {
      result[j] <- .ag_nearest(group, queries[j, ])
    }
  }
  as.integer(result)
}

# Perimeter of each seed trio, summed in the canonical order.
.seed_perimeter <- function(red, nearest_o1, nearest_o2) {
  point_distance <- function(a, b) {
    dx <- a[, 1L] - b[, 1L]
    dy <- a[, 2L] - b[, 2L]
    sqrt(dx * dx + dy * dy)
  }
  (point_distance(nearest_o1, red) + point_distance(nearest_o2, red)) +
    point_distance(nearest_o1, nearest_o2)
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
#' @param X_vars Ignored, and warned about if supplied. Accepted only so that
#'   a call written for [sam_match()] fails loudly rather than silently
#'   matching on something else: this algorithm works in propensity score
#'   space and never looks at covariates.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#' @param caliper Numeric caliper or `"auto"` to estimate it from the data.
#'   Default `"auto"`.
#' @param gps_space Propensity score scale used for matching, either `"raw"`
#'   or `"logit"`. Default `"raw"`. Named to match [gps_candidate_search()].
#' @param ps_space Deprecated alias for `gps_space`.
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
#'     \item{max_possible_rate}{The highest rate the group sizes allow; see
#'       [sam_match()].}
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
match_3way <- function(data, search, gps, X_vars = NULL,
                        treatment_var = "T", caliper = "auto",
                        gps_space = c("raw", "logit"), top_n = 10L,
                        reference_level = NULL, ps_space = NULL) {
  gps_space <- match.arg(gps_space)

  # `ps_space` named the same raw/logit toggle that gps_candidate_search calls
  # `gps_space`, and the README used both within one section.
  if (!is.null(ps_space)) {
    warning("`ps_space` is deprecated; use `gps_space` instead.",
            call. = FALSE)
    gps_space <- match.arg(ps_space, c("raw", "logit"))
  }

  # Covariates play no part in this algorithm -- it matches in two-dimensional
  # propensity score space. The argument was previously accepted and silently
  # discarded, so a caller passing covariates and expecting them to matter got
  # no signal at all.
  if (!is.null(X_vars)) {
    warning("`X_vars` is ignored by match_3way(): matching uses the ",
            "propensity score space, not covariate distances.", call. = FALSE)
  }

  top_n <- require_positive_int(top_n, "top_n")
  check_fingerprint(search, data, treatment_var)
  gps <- validate_gps(data, gps, treatment_var, "match_3way")
  check_gps_fingerprint(search, gps)
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
  } else {
    reference_level <- treatment_level(reference_level)
  }
  stopifnot(reference_level %in% all_levels)
  ps_levels <- setdiff(all_levels, reference_level)
  stopifnot(length(ps_levels) == 2L)

  ps_used <- transform_ps(gps[, ps_levels, drop = FALSE], gps_space)
  colnames(ps_used) <- ps_levels

  # --- Caliper
  if (identical(caliper, "auto")) {
    caliper <- calc_caliper_3way(ps_used, labels)
  }
  if (!is.numeric(caliper) || length(caliper) != 1L ||
      !is.finite(caliper) || caliper <= 0) {
    stop("`caliper` must be a single finite number greater than zero.",
         call. = FALSE)
  }

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

  # --- Candidate search
  # Candidate generation follows one canonical rule, written so that its result
  # depends only on the data and never on how a spatial index happens to order
  # equidistant points:
  #
  # * the seed partner is the nearest available subject, lowest row first
  #   among equidistant ones;
  # * candidate trios are ranked by (perimeter key, o1 row, o2 row), the key
  #   being order_key() of the perimeter;
  # * every distance is computed as sqrt(dx * dx + dy * dy) and every
  #   perimeter as (d_red_o1 + d_red_o2) + d_o1_o2, exactly as the Python
  #   implementation does, so both languages see bit-identical values.
  #
  # Spatial indexes only propose points; membership is decided on those exact
  # values. Ties previously followed the k-d tree's traversal order, which
  # depended on how the tree was built -- and differed between R and Python.
  group_o1 <- .active_group(coords_o1)
  group_o2 <- .active_group(coords_o2)

  matched_flag <- rep(FALSE, n_red)
  exhausted    <- rep(FALSE, n_red)
  counters     <- integer(n_red)

  # Search radius used whenever the caliper is tighter than a subject's seed
  # trio. Any trio kept needs perimeter < seed AND perimeter <= caliper, and
  # each side of a triangle is at most half its perimeter, so nothing farther
  # than min(seed, caliper) / 2 from the base subject can contribute. When the
  # caliper binds, the bound is inclusive, so it is widened by a slack far
  # above rounding error and far below any meaningful distance; the few extra
  # points gathered fail the perimeter test.
  caliper_half <- caliper / 2 + 1e-6 * (1 + max(abs(ps_used)))
  search_half_width <- function(seed) {
    ifelse(seed <= caliper, seed / 2, caliper_half)
  }

  x1 <- unname(coords_o1[, 1L]); y1 <- unname(coords_o1[, 2L])
  x2 <- unname(coords_o2[, 1L]); y2 <- unname(coords_o2[, 2L])

  # --- Shared pool of candidate trios, as a binary min-heap
  # The key is (perimeter key, insertion order); the second component keeps
  # equal perimeters in the order they were ranked. The raw perimeter travels
  # with the entry for the output.
  heap_capacity <- 256L
  heap_key <- numeric(heap_capacity)
  heap_perim <- numeric(heap_capacity)
  heap_seq <- integer(heap_capacity)
  heap_red <- integer(heap_capacity)
  heap_o1 <- integer(heap_capacity)
  heap_o2 <- integer(heap_capacity)
  heap_size <- 0L
  heap_pushed <- 0L

  heap_before_key <- function(a, key, seq_no) {
    heap_key[a] < key || (heap_key[a] == key && heap_seq[a] < seq_no)
  }

  # Entries move by carrying a hole up or down and writing the travelling
  # record once at the end.
  heap_move <- function(to, from) {
    heap_key[to] <<- heap_key[from]
    heap_perim[to] <<- heap_perim[from]
    heap_seq[to] <<- heap_seq[from]
    heap_red[to] <<- heap_red[from]
    heap_o1[to] <<- heap_o1[from]
    heap_o2[to] <<- heap_o2[from]
    invisible(NULL)
  }

  heap_place <- function(at, key, perimeter, seq_no, red, o1_index, o2_index) {
    heap_key[at] <<- key
    heap_perim[at] <<- perimeter
    heap_seq[at] <<- seq_no
    heap_red[at] <<- red
    heap_o1[at] <<- o1_index
    heap_o2[at] <<- o2_index
    invisible(NULL)
  }

  heap_push <- function(key, perimeter, red, o1_index, o2_index) {
    if (heap_size == heap_capacity) {
      grown <- heap_capacity * 2L
      length(heap_key) <<- grown; length(heap_perim) <<- grown
      length(heap_seq) <<- grown; length(heap_red) <<- grown
      length(heap_o1) <<- grown; length(heap_o2) <<- grown
      heap_capacity <<- grown
    }
    heap_size <<- heap_size + 1L
    heap_pushed <<- heap_pushed + 1L
    seq_no <- heap_pushed
    hole <- heap_size
    while (hole > 1L) {
      parent <- hole %/% 2L
      if (heap_before_key(parent, key, seq_no)) break
      heap_move(hole, parent)
      hole <- parent
    }
    heap_place(hole, key, perimeter, seq_no, red, o1_index, o2_index)
    invisible(NULL)
  }

  # Removes and returns the smallest entry as (perimeter, red, o1, o2).
  heap_pop <- function() {
    top <- c(heap_perim[1L], heap_red[1L], heap_o1[1L], heap_o2[1L])
    last <- heap_size
    moving_key <- heap_key[last]
    moving_perim <- heap_perim[last]
    moving_seq <- heap_seq[last]
    moving_red <- heap_red[last]
    moving_o1 <- heap_o1[last]
    moving_o2 <- heap_o2[last]
    heap_size <<- heap_size - 1L
    hole <- 1L
    repeat {
      left <- hole * 2L
      if (left > heap_size) break
      best <- left
      right <- left + 1L
      if (right <= heap_size &&
          (heap_key[right] < heap_key[left] ||
             (heap_key[right] == heap_key[left] &&
                heap_seq[right] < heap_seq[left]))) {
        best <- right
      }
      if (!heap_before_key(best, moving_key, moving_seq)) break
      heap_move(hole, best)
      hole <- best
    }
    if (heap_size > 0L) {
      heap_place(hole, moving_key, moving_perim, moving_seq, moving_red,
                 moving_o1, moving_o2)
    }
    top
  }

  # --- Candidate generation
  # Ranks one base subject's candidate trios, given its seed partners and the
  # available subjects inside the search radius. Returns NULL, and marks the
  # subject exhausted, when it has none.
  rank_candidates <- function(i, nb, nbg, seed, nbs, ngs) {
    qx <- coords_red[[i, 1L]]
    qy <- coords_red[[i, 2L]]
    cand_o1 <- integer(0); cand_o2 <- integer(0); cand_perim <- numeric(0)

    # Include the seed trio when it satisfies the caliper
    if (seed <= caliper) {
      cand_o1 <- nb; cand_o2 <- nbg; cand_perim <- seed
    }

    n_a <- length(nbs)
    if (n_a > 0L && length(ngs) > 0L) {
      dx <- x1[nbs] - qx; dy <- y1[nbs] - qy
      d_pr_A <- sqrt(dx * dx + dy * dy)
      dx <- x2[ngs] - qx; dy <- y2[ngs] - qy
      d_pr_B <- sqrt(dx * dx + dy * dy)
      DX <- outer(x1[nbs], x2[ngs], "-")
      DY <- outer(y1[nbs], y2[ngs], "-")
      perim <- outer(d_pr_A, d_pr_B, "+") + sqrt(DX * DX + DY * DY)
      keep <- which(perim < seed & perim <= caliper)
      if (length(keep) > 0L) {
        cand_o1    <- c(cand_o1,    nbs[(keep - 1L) %% n_a + 1L])
        cand_o2    <- c(cand_o2,    ngs[(keep - 1L) %/% n_a + 1L])
        cand_perim <- c(cand_perim, perim[keep])
      }
    }

    if (length(cand_perim) == 0L) {
      exhausted[i] <<- TRUE; counters[i] <<- 0L
      return(invisible(NULL))
    }

    keys <- order_key(cand_perim)
    ord <- order(keys, cand_o1, cand_o2)[seq_len(min(top_n, length(cand_perim)))]
    counters[i] <<- length(ord)
    list(key = keys[ord], perimeter = cand_perim[ord],
         o1 = cand_o1[ord], o2 = cand_o2[ord])
  }

  # Searches one subject again against the subjects still available.
  regenerate <- function(i) {
    q <- coords_red[i, ]
    nb <- .ag_nearest(group_o1, q)
    if (is.na(nb)) {
      exhausted[i] <<- TRUE; counters[i] <<- 0L
      return(invisible(NULL))
    }
    nbg <- .ag_nearest(group_o2, coords_o1[nb, ])
    if (is.na(nbg)) {
      exhausted[i] <<- TRUE; counters[i] <<- 0L
      return(invisible(NULL))
    }
    seed <- .seed_perimeter(coords_red[i, , drop = FALSE],
                            coords_o1[nb, , drop = FALSE],
                            coords_o2[nbg, , drop = FALSE])
    half <- search_half_width(seed)
    radius2 <- half * half
    ranked <- rank_candidates(i, nb, nbg, seed,
                              .ag_within(group_o1, q, radius2),
                              .ag_within(group_o2, q, radius2))
    for (k in seq_along(ranked$key)) {
      heap_push(ranked$key[k], ranked$perimeter[k], i, ranked$o1[k],
                ranked$o2[k])
    }
    invisible(NULL)
  }

  # Initial pass. Nothing has been consumed yet, so every subject's nearest
  # neighbours come from one batched query to a compiled tree.
  #
  # The initial candidates -- most of all candidates ever generated -- are
  # sorted once into a static queue instead of being pushed through the heap
  # one at a time; only regenerated candidates use the heap. The next trio is
  # whichever queue head is smaller by (key, insertion order), and since every
  # initial entry was inserted before any regenerated one, the order trios
  # come out in is exactly that of a single heap.
  initial <- vector("list", n_red)
  if (n_red > 0L) {
    nb_all <- .nearest_untouched(group_o1, coords_o1, coords_red)
    nbg_all <- .nearest_untouched(group_o2, coords_o2,
                                  coords_o1[nb_all, , drop = FALSE])
    seed_all <- .seed_perimeter(coords_red,
                                coords_o1[nb_all, , drop = FALSE],
                                coords_o2[nbg_all, , drop = FALSE])
    half_all <- search_half_width(seed_all)
    radius2_all <- half_all * half_all
    for (i in seq_len(n_red)) {
      q <- coords_red[i, ]
      ranked <- rank_candidates(i, nb_all[i], nbg_all[i], seed_all[i],
                                .ag_within(group_o1, q, radius2_all[i]),
                                .ag_within(group_o2, q, radius2_all[i]))
      if (!is.null(ranked)) {
        ranked$red <- rep(i, length(ranked$key))
        initial[[i]] <- ranked
      }
    }
  }
  # as.numeric()/as.integer() turn the NULL that unlist() gives when no
  # subject has a candidate into an empty vector.
  gather <- function(field) unlist(lapply(initial, `[[`, field), use.names = FALSE)
  static_key <- as.numeric(gather("key"))
  static_perim <- as.numeric(gather("perimeter"))
  static_red <- as.integer(gather("red"))
  static_o1 <- as.integer(gather("o1"))
  static_o2 <- as.integer(gather("o2"))
  rm(initial)
  static_n <- length(static_key)
  # Insertion order is the position in the concatenation, so ordering by
  # (key, position) is the heap's (key, insertion order).
  static_order <- order(static_key, seq_len(static_n))
  static_key <- static_key[static_order]
  static_perim <- static_perim[static_order]
  static_red <- static_red[static_order]
  static_o1 <- static_o1[static_order]
  static_o2 <- static_o2[static_order]
  static_pos <- 1L
  # Regenerated entries must sort after every initial one at an equal key.
  heap_pushed <- static_n

  # --- Global greedy matching
  # Results accumulate into preallocated vectors and are assembled once, rather
  # than growing a list of one-row data frames.
  n_matched <- 0L
  matched_anchor <- integer(n_red)
  matched_loss <- numeric(n_red)
  matched_perimeter <- numeric(n_red)
  matched_choice <- matrix(NA_integer_, n_red, length(groups),
                           dimnames = list(NULL, groups))
  matched_dist <- matrix(NA_real_, n_red, length(groups),
                         dimnames = list(NULL, groups))

  repeat {
    static_left <- static_pos <= static_n
    if (heap_size > 0L &&
        (!static_left || heap_key[1L] < static_key[static_pos])) {
      popped <- heap_pop()
      perim_popped <- popped[1L]
      i <- as.integer(popped[2L])
      b <- as.integer(popped[3L])
      g <- as.integer(popped[4L])
    } else if (static_left) {
      perim_popped <- static_perim[static_pos]
      i <- static_red[static_pos]
      b <- static_o1[static_pos]
      g <- static_o2[static_pos]
      static_pos <- static_pos + 1L
    } else {
      break
    }

    # Skip candidates belonging to an already matched base subject
    if (matched_flag[i]) next

    if (group_o1$active[b] && group_o2$active[g]) {
      # Accept the candidate trio
      matched_flag[i] <- TRUE
      .ag_deactivate(group_o1, b)
      .ag_deactivate(group_o2, g)

      row_red_global <- rows_by_level[[red_level]][i]
      row_o1_global  <- rows_by_level[[o1]][b]
      row_o2_global  <- rows_by_level[[o2]][g]

      # Relabel subjects according to the anchor and comparator groups
      rows_this_trio <- stats::setNames(
        list(row_red_global, row_o1_global, row_o2_global),
        c(red_level, o1, o2)
      )
      row_anchor <- rows_this_trio[[anchor_level]]

      n_matched <- n_matched + 1L
      matched_anchor[n_matched] <- row_anchor

      # Anchor-to-comparator propensity score distances
      for (gg in groups) {
        row_group <- rows_this_trio[[gg]]
        matched_choice[n_matched, gg] <- row_group
        matched_dist[n_matched, gg] <-
          sqrt(sum((ps_used[row_anchor, ] - ps_used[row_group, ])^2))
      }

      matched_loss[n_matched] <- sum(matched_dist[n_matched, ])

      # Full propensity-score triangle perimeter used by this method
      matched_perimeter[n_matched] <- perim_popped
    } else {
      # Refresh candidates when all current options become unavailable
      counters[i] <- counters[i] - 1L
      if (counters[i] <= 0L && !exhausted[i]) {
        regenerate(i)
      }
    }
  }

  # --- Output
  keep <- seq_len(n_matched)
  named_groups <- stats::setNames(groups, groups)

  matched <- build_matched_frame(
    anchor = matched_anchor[keep],
    group_choices = lapply(named_groups, function(g) matched_choice[keep, g]),
    group_distances = lapply(named_groups, function(g) matched_dist[keep, g]),
    loss = matched_loss[keep],
    groups = groups,
    extra = list(rassen_perimeter = matched_perimeter[keep])
  )

  c(
    list(matched = matched),
    summarize_matching(matched, anchor_rows, rows_by_level[groups]),
    list(
      caliper = caliper,
      red_level = red_level,
      reference_level = reference_level
    )
  )
}
