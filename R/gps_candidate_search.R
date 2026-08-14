################################################################################
### GPS-Guided Candidate Search
##
## Author: Sungjun Hong
################################################################################

#' GPS-guided candidate search
#'
#' @description
#' For each anchor subject, identifies the `top_m` nearest subjects from
#' each comparator group using Euclidean distance in generalized propensity
#' score (GPS) space. The resulting candidate pools are used by [sam_match()]
#' for subsequent Mahalanobis-distance matching.
#'
#' @param data A `data.frame` containing the treatment variable.
#' @param gps Numeric matrix of generalized propensity scores, with one row
#'   per subject and one column per treatment group.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#' @param anchor_level Anchor treatment group. Default `"A"`.
#' @param top_m Number of candidates retained per anchor and comparator group.
#'   Must be a positive integer. Default `10`.
#' @param gps_space GPS scale used to calculate Euclidean distance.
#'   Either `"raw"` or `"logit"`. Default `"raw"`.
#'
#' @return A list with:
#'   \describe{
#'     \item{anchor_rows}{Row indices of anchor subjects.}
#'     \item{groups}{Comparator treatment groups.}
#'     \item{candidates}{Candidate row indices for each anchor and
#'       comparator group.}
#'     \item{fingerprint}{Identifies the exact frame, in the exact row order,
#'       that the search ran on. Every later stage verifies it, because
#'       matched sets are stored as positional row indices: re-sorting,
#'       filtering or rewriting `data` in between would silently repoint them
#'       at other subjects. Adding a column, such as an outcome, is still
#'       allowed.}
#'     \item{gps_fingerprint}{Identifies the GPS matrix the candidate pools
#'       were selected against, so that a later stage handed a different one
#'       fails rather than reporting diagnostics for scores that never
#'       produced the match.}
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
#' length(search$anchor_rows)
#' search$candidates[[1]]
#'
#' @export
gps_candidate_search <- function(data, gps, treatment_var = "T", anchor_level = "A",
                                  top_m = 10, gps_space = c("raw", "logit")) {
  gps_space <- match.arg(gps_space)
  top_m <- require_positive_int(top_m, "top_m")
  labels <- treatment_labels(data, treatment_var)
  anchor_level <- treatment_level(anchor_level)
  gps <- validate_gps(data, gps, treatment_var, "gps_candidate_search")

  if (!(anchor_level %in% colnames(gps))) {
    stop("`anchor_level` \"", anchor_level, "\" is not a column of `gps`. ",
         "GPS columns: ",
         paste0("\"", colnames(gps), "\"", collapse = ", "), ".",
         call. = FALSE)
  }

  gps_used <- transform_ps(gps, gps_space)

  groups <- setdiff(colnames(gps), anchor_level)
  anchor_rows <- require_rows(which(labels == anchor_level), anchor_level,
                              treatment_var, available = labels)
  X_anchor <- gps_used[anchor_rows, , drop = FALSE]

  # Find the nearest candidates in each comparator group.
  #
  # A k-d tree query replaces the full anchor-by-comparator distance matrix,
  # which allocated the product of the two group sizes -- 1.3 GB at n = 20,000
  # -- and then fully sorted every row to keep top_m of it.
  candidates_by_group <- lapply(groups, function(g) {
    group_rows <- require_rows(which(labels == g), g, treatment_var,
                               available = labels)
    X_group <- gps_used[group_rows, , drop = FALSE]
    n_group <- length(group_rows)
    n_anchor <- nrow(X_anchor)
    m <- min(top_m, n_group)

    # One neighbour beyond the cutoff, so a tie straddling it is detectable.
    k_query <- min(m + 1L, n_group)
    neighbours <- RANN::nn2(data = X_group, query = X_anchor, k = k_query)

    # The tree's own distances are not used for ordering. ANN accumulates a
    # distance along the path it took to reach a point, so two points at
    # identical coordinates can come back differing in the last bit, and
    # sorting on that would follow the traversal rather than the row order the
    # full-matrix implementation used. Distances are recomputed here with that
    # implementation's formula, for every returned pair in one pass.
    anchor_sq <- rowSums(X_anchor^2)
    group_sq <- rowSums(X_group^2)

    flat_anchor <- rep(seq_len(n_anchor), times = k_query)
    flat_candidate <- as.integer(neighbours$nn.idx)

    d2 <- anchor_sq[flat_anchor] + group_sq[flat_candidate] -
      2 * rowSums(X_anchor[flat_anchor, , drop = FALSE] *
                    X_group[flat_candidate, , drop = FALSE])
    d2[d2 < 0] <- 0
    exact <- matrix(sqrt(d2), nrow = n_anchor, ncol = k_query)

    lapply(seq_len(n_anchor), function(i) {
      idx <- neighbours$nn.idx[i, ]
      distance <- exact[i, ]
      ranked <- order(distance, idx)

      # More candidates may be tied at the cutoff than there are slots left, in
      # which case the tree returned an arbitrary subset of them. Those anchors
      # are re-resolved against every candidate in the group.
      contested <- k_query > m &&
        distance[ranked[k_query]] <= distance[ranked[m]] * (1 + .SAM_TIE_TOL)

      if (contested) {
        # Accumulated column by column rather than through `%*%`. A BLAS
        # matrix-vector product may reduce different rows in different orders,
        # depending on blocking and vector width, so two points at identical
        # coordinates can come out differing in the last bit -- which is
        # precisely what this branch exists to rule out. Summing explicitly
        # gives every row the same accumulation order.
        cross <- numeric(n_group)
        for (k in seq_len(ncol(X_group))) {
          cross <- cross + X_group[, k] * X_anchor[i, k]
        }

        all_d2 <- anchor_sq[i] + group_sq - 2 * cross
        all_d2[all_d2 < 0] <- 0
        group_rows[order(sqrt(all_d2), seq_len(n_group))[seq_len(m)]]
      } else {
        group_rows[idx[ranked[seq_len(m)]]]
      }
    })
  })
  names(candidates_by_group) <- groups

  # Organize candidate lists by anchor
  candidates <- lapply(seq_along(anchor_rows), function(i) {
    stats::setNames(lapply(groups, function(g) candidates_by_group[[g]][[i]]), groups)
  })

  list(
    anchor_rows = anchor_rows,
    groups = groups,
    candidates = candidates,
    # Matched sets are stored as positional row indices, so every later stage
    # needs this exact frame in this exact row order, and the GPS the pools
    # were selected against.
    fingerprint = data_fingerprint(data, treatment_var),
    gps_fingerprint = gps_fingerprint(gps)
  )
}
