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
#'       matched sets are stored as positional row indices: re-sorting or
#'       filtering `data` in between would silently repoint them at other
#'       subjects. Adding a column, such as an outcome, is still allowed.}
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
  stopifnot(nrow(gps) == nrow(data))

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
    m <- min(top_m, n_group)

    # One neighbour beyond the cutoff, so a tie straddling it is detectable.
    k_query <- min(m + 1L, n_group)
    neighbours <- RANN::nn2(data = X_group, query = X_anchor, k = k_query)

    neighbour_idx <- neighbours$nn.idx[, seq_len(m), drop = FALSE]
    neighbour_dist <- neighbours$nn.dists[, seq_len(m), drop = FALSE]

    # The tree resolves equidistant candidates arbitrarily, and when more of
    # them are tied at the cutoff than there are slots left it returns an
    # arbitrary subset. A stable sort by row order does neither, so those
    # anchors are re-resolved against every candidate below.
    contested <- if (k_query > m) {
      neighbours$nn.dists[, k_query] <=
        neighbour_dist[, m] * (1 + .SAM_TIE_TOL)
    } else {
      rep(FALSE, nrow(X_anchor))
    }

    group_sq <- rowSums(X_group^2)

    lapply(seq_len(nrow(X_anchor)), function(i) {
      if (contested[i]) {
        # Recomputed exactly as the full-matrix implementation did, so that
        # the tie-break is identical to the last bit.
        d2 <- sum(X_anchor[i, ]^2) + group_sq -
          2 * as.numeric(X_group %*% X_anchor[i, ])
        d2[d2 < 0] <- 0
        group_rows[order(sqrt(d2), seq_len(n_group))[seq_len(m)]]
      } else {
        # Ties within the returned set still follow row order, matching the
        # stable sort this replaces.
        j <- neighbour_idx[i, ]
        group_rows[j[order(neighbour_dist[i, ], j)]]
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
    # needs this exact frame in this exact row order.
    fingerprint = data_fingerprint(data, treatment_var)
  )
}
