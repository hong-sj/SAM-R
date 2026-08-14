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

  # Transform GPS values to the logit scale if requested
  if (gps_space == "logit") {
    eps <- 1e-6
    gps_used <- stats::qlogis(pmin(pmax(gps, eps), 1 - eps))
  } else {
    gps_used <- gps
  }

  groups <- setdiff(colnames(gps), anchor_level)
  anchor_rows <- require_rows(which(labels == anchor_level), anchor_level,
                              treatment_var, available = labels)
  X_anchor <- gps_used[anchor_rows, , drop = FALSE]

  # Compute Euclidean GPS distances for each comparator group
  candidates_by_group <- lapply(groups, function(g) {
    group_rows <- require_rows(which(labels == g), g, treatment_var,
                               available = labels)
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
