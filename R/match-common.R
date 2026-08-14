################################################################################
### Result Assembly Shared by the Matching Algorithms
##
## sam_match() and match_3way() run genuinely different algorithms -- a
## Mahalanobis greedy in covariate space and a KD-tree perimeter greedy in
## two-dimensional propensity score space -- but they return the same shape of
## result. Only that shared plumbing lives here; the engines stay separate.
################################################################################

# Probabilities are clipped away from the boundary before the logit, which
# would otherwise map 0 and 1 to infinities.
.SAM_LOGIT_EPS <- 1e-6

# Relative tolerance for deciding that two candidate distances are tied. A k-d
# tree accumulates a distance in a different order than the matrix formula
# does, so the two do not agree to the last bit, and an exact comparison would
# silently drop ties that the full-matrix implementation kept.
.SAM_TIE_TOL <- 1e-9


#' Propensity scores on the requested scale
#'
#' @param values Numeric matrix of propensity scores.
#' @param space Either `"raw"` or `"logit"`.
#'
#' @return `values`, transformed.
#'
#' @keywords internal
#' @noRd
transform_ps <- function(values, space) {
  space <- match.arg(space, c("raw", "logit"))

  if (space == "raw") {
    return(values)
  }

  stats::qlogis(pmin(pmax(values, .SAM_LOGIT_EPS), 1 - .SAM_LOGIT_EPS))
}


#' Column order of a matched-set frame
#'
#' @param groups Character vector of comparator treatment groups.
#' @param extra_columns Additional columns appended after `loss`.
#'
#' @return Character vector of column names.
#'
#' @keywords internal
#' @noRd
matched_frame_columns <- function(groups, extra_columns = character(0)) {
  c("matched_set_id", "anchor", groups, paste0("dist_", groups), "loss",
    extra_columns)
}


#' An empty matched-set frame with the correct schema
#'
#' The column set must survive a run that matched nothing, otherwise every
#' caller has to special-case the empty result.
#'
#' @inheritParams matched_frame_columns
#'
#' @return A zero-row `data.frame`.
#'
#' @keywords internal
#' @noRd
empty_matched_frame <- function(groups, extra_columns = character(0)) {
  columns <- matched_frame_columns(groups, extra_columns)

  frame <- data.frame(matched_set_id = integer(0), anchor = integer(0))
  for (g in groups) {
    frame[[g]] <- integer(0)
  }
  for (column in c(paste0("dist_", groups), "loss", extra_columns)) {
    frame[[column]] <- numeric(0)
  }

  frame[, columns, drop = FALSE]
}


#' Assemble a matched-set frame from column vectors
#'
#' Accumulating one-row data frames and rbind-ing them is quadratic in the
#' number of matched sets, so the engines collect vectors and assemble once.
#'
#' @param anchor Integer vector of anchor row indices, one per matched set.
#' @param group_choices Named list of integer vectors, one per comparator
#'   group, holding the matched row index for each set.
#' @param group_distances Named list of numeric vectors, in the same shape.
#' @param loss Numeric vector of matched-set losses.
#' @param groups Character vector fixing the comparator group order.
#' @param extra Named list of additional columns appended after `loss`.
#'
#' @return A `data.frame` with the standard matched-set schema.
#'
#' @keywords internal
#' @noRd
build_matched_frame <- function(anchor, group_choices, group_distances, loss,
                                groups, extra = list()) {
  if (length(anchor) == 0L) {
    return(empty_matched_frame(groups, names(extra)))
  }

  frame <- data.frame(
    matched_set_id = seq_along(anchor),
    anchor = as.integer(anchor)
  )
  for (g in groups) {
    frame[[g]] <- as.integer(group_choices[[g]])
  }
  for (g in groups) {
    frame[[paste0("dist_", g)]] <- as.numeric(group_distances[[g]])
  }
  frame$loss <- as.numeric(loss)
  for (column in names(extra)) {
    frame[[column]] <- extra[[column]]
  }

  rownames(frame) <- NULL
  frame
}


#' Recover comparator group names from a matched-set frame
#'
#' Lets the diagnostic functions be called with just the matched sets, instead
#' of requiring the caller to carry `search$groups` alongside.
#'
#' @param matched A matched-set `data.frame`.
#'
#' @return Character vector of comparator group names.
#'
#' @keywords internal
#' @noRd
groups_from_matched <- function(matched) {
  reserved <- c("matched_set_id", "anchor", "loss", "rassen_perimeter")
  columns <- names(matched)

  groups <- columns[!(columns %in% reserved) & !startsWith(columns, "dist_")]

  if (length(groups) == 0L) {
    stop("Could not infer the comparator groups from `matched`; pass ",
         "`groups` explicitly.", call. = FALSE)
  }

  groups
}


#' Highest matching rate the group sizes allow
#'
#' Each matched set consumes one subject from every comparator group, so the
#' rate is capped by the smallest comparator group over the anchor count. On
#' the bundled four-group data that cap is 59/448 = 0.132 and SAM reaches it
#' exactly -- reporting the rate alone makes a saturated result read as a poor
#' 13%.
#'
#' @param group_rows Named list of row indices per comparator group.
#' @param anchor_rows Row indices of every anchor subject.
#'
#' @return `min(comparator group sizes) / n_anchor`.
#'
#' @keywords internal
#' @noRd
max_possible_rate <- function(group_rows, anchor_rows) {
  if (length(anchor_rows) == 0L) {
    return(NaN)
  }

  sizes <- vapply(group_rows, length, integer(1))
  smallest <- if (length(sizes) == 0L) 0L else min(sizes)

  min(1, smallest / length(anchor_rows))
}


#' Unmatched anchors, the matching rate, and its ceiling
#'
#' @param matched A matched-set `data.frame`.
#' @param anchor_rows Row indices of every anchor subject.
#' @param group_rows Named list of row indices per comparator group.
#'
#' @return A list with `unmatched_anchor_rows`, `matching_rate` and
#'   `max_possible_rate`.
#'
#' @keywords internal
#' @noRd
summarize_matching <- function(matched, anchor_rows, group_rows) {
  anchor_rows <- as.integer(anchor_rows)

  matched_anchor_rows <- if (nrow(matched) > 0L) {
    as.integer(matched$anchor)
  } else {
    integer(0)
  }

  list(
    unmatched_anchor_rows = setdiff(anchor_rows, matched_anchor_rows),
    matching_rate = if (length(anchor_rows) > 0L) {
      nrow(matched) / length(anchor_rows)
    } else {
      NaN
    },
    max_possible_rate = max_possible_rate(group_rows, anchor_rows)
  )
}
