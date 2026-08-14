################################################################################
### Shared Input Validation and Canonicalization
##
## Every SAM entry point routes its inputs through this file so that the same
## input is interpreted the same way at every stage. Historically each function
## made its own assumptions, which allowed a mismatch to pass through silently
## and produce an empty or wrong result rather than an error.
################################################################################


#' Treatment column as canonical character labels
#'
#' `estimate_gps_multinom()` builds its GPS matrix from a `factor()`, so every
#' GPS column label and every group name flowing through the pipeline is a
#' character string. Comparing those labels against a raw numeric or factor
#' treatment column relies on R's coercion rules, which differ between `==`,
#' `setdiff()` and `relevel()`. Canonicalizing here keeps every stage on one
#' representation.
#'
#' @param data A `data.frame` containing the treatment variable.
#' @param treatment_var Name of the treatment variable.
#'
#' @return Character vector of treatment labels, one per row.
#'
#' @keywords internal
#' @noRd
treatment_labels <- function(data, treatment_var) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (length(treatment_var) != 1L || !is.character(treatment_var)) {
    stop("`treatment_var` must be a single column name.", call. = FALSE)
  }
  if (!(treatment_var %in% names(data))) {
    stop("Treatment column not found: ", treatment_var, call. = FALSE)
  }

  labels <- as.character(data[[treatment_var]])

  # A missing treatment label is not a level. Left alone it compares to NA
  # against every group, so `which()` drops the row from all of them and the
  # subject silently disappears from the analysis.
  if (anyNA(labels)) {
    stop("`", treatment_var, "` contains ", sum(is.na(labels)),
         " missing treatment label(s). Drop or impute these rows before ",
         "matching.", call. = FALSE)
  }

  labels
}


#' Treatment level in the canonical representation
#'
#' Guards in particular against a numeric level reaching `stats::relevel()`,
#' which interprets a numeric `ref` as a *position* in `levels()` rather than
#' as a level name.
#'
#' @param level A treatment level, or `NULL`.
#'
#' @return `NULL`, or the level as a single character string.
#'
#' @keywords internal
#' @noRd
treatment_level <- function(level) {
  if (is.null(level)) {
    return(NULL)
  }
  if (length(level) != 1L) {
    stop("A treatment level must be a single value, got length ",
         length(level), ".", call. = FALSE)
  }
  if (is.na(level)) {
    stop("A treatment level must not be NA.", call. = FALSE)
  }

  as.character(level)
}


#' Require that a treatment level selected at least one row
#'
#' An empty selection means the requested level is not present under the
#' canonical character representation, which otherwise propagates silently as
#' an empty matched set and a `NaN` matching rate.
#'
#' @param rows Integer row indices selected by the level.
#' @param level The treatment level that was requested.
#' @param treatment_var Name of the treatment variable.
#' @param available Optional character vector of levels present in the data,
#'   listed in the error message.
#'
#' @return `rows`, invisibly unchanged.
#'
#' @keywords internal
#' @noRd
require_rows <- function(rows, level, treatment_var, available = NULL) {
  if (length(rows) == 0L) {
    message_text <- paste0(
      "No rows found with ", treatment_var, " == \"", level, "\". ",
      "Treatment levels are compared as character strings; check that the ",
      "level matches the values in the treatment column."
    )
    if (!is.null(available)) {
      message_text <- paste0(
        message_text, " Levels present: ",
        paste0("\"", unique(available), "\"", collapse = ", "), "."
      )
    }
    stop(message_text, call. = FALSE)
  }

  rows
}


#' Require that every covariate column is present
#'
#' @param data A `data.frame`.
#' @param X_vars Character vector of covariate column names.
#'
#' @return `X_vars`, unchanged.
#'
#' @keywords internal
#' @noRd
require_covariates <- function(data, X_vars) {
  if (!is.character(X_vars) || length(X_vars) == 0L) {
    stop("`X_vars` must be a non-empty character vector of column names.",
         call. = FALSE)
  }

  missing_columns <- setdiff(X_vars, names(data))
  if (length(missing_columns) > 0L) {
    stop("Covariate column(s) not found in data: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  duplicated_columns <- unique(X_vars[duplicated(X_vars)])
  if (length(duplicated_columns) > 0L) {
    stop("Duplicated covariate column(s) in `X_vars`: ",
         paste(duplicated_columns, collapse = ", "), call. = FALSE)
  }

  X_vars
}


#' Covariates as a validated numeric matrix
#'
#' Non-numeric and non-finite covariates are rejected here, with a message
#' naming the offending columns. Reaching the linear algebra with `NA` present
#' is not detectable downstream: `solve()` returns an all-`NA` matrix without
#' raising, so a single missing value silently destroys every Mahalanobis
#' distance and returns an empty match.
#'
#' @param data A `data.frame` containing the covariates.
#' @param X_vars Character vector of covariate column names.
#' @param rows Optional integer row indices to restrict the returned matrix
#'   to. Validation always covers the whole column, since the pooled
#'   covariance uses every row.
#'
#' @return Numeric matrix with `length(X_vars)` columns.
#'
#' @keywords internal
#' @noRd
covariate_matrix <- function(data, X_vars, rows = NULL) {
  require_covariates(data, X_vars)

  columns <- data[, X_vars, drop = FALSE]

  # Logical columns are accepted as 0/1 indicators; factors and characters are
  # not, because silently coercing a factor to its integer codes would treat
  # arbitrary level ordering as a distance.
  non_numeric <- X_vars[!vapply(
    columns, function(x) is.numeric(x) || is.logical(x), logical(1)
  )]
  if (length(non_numeric) > 0L) {
    stop("Covariate column(s) are not numeric: ",
         paste(non_numeric, collapse = ", "),
         ". Categorical covariates must be encoded (for example as indicator ",
         "columns via stats::model.matrix()) before being passed as X_vars.",
         call. = FALSE)
  }

  values <- as.matrix(data.frame(lapply(columns, as.numeric),
                                 check.names = FALSE))
  colnames(values) <- X_vars

  finite <- is.finite(values)
  if (!all(finite)) {
    offending <- vapply(
      which(!apply(finite, 2, all)),
      function(j) paste0(X_vars[j], " (", sum(!finite[, j]), " rows)"),
      character(1)
    )
    stop("Covariate column(s) contain missing or non-finite values: ",
         paste(offending, collapse = ", "),
         ". SAM does not impute; drop or impute these rows before matching.",
         call. = FALSE)
  }

  if (is.null(rows)) values else values[rows, , drop = FALSE]
}


#' Require a positive integer
#'
#' @param value The value to check.
#' @param name Parameter name, used in the error message.
#'
#' @return `value` as an integer.
#'
#' @keywords internal
#' @noRd
require_positive_int <- function(value, name) {
  # `TRUE` coerces to the integer 1 in R but is not a candidate count.
  if (is.logical(value) || !is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value != as.integer(value) || value < 1L) {
    stop("`", name, "` must be a single positive integer, got ",
         paste(deparse(value, nlines = 1L), collapse = " "), ".",
         call. = FALSE)
  }

  as.integer(value)
}


#' Fingerprint identifying the exact frame matching ran on
#'
#' Matched sets are stored as *positional* row indices, so every stage after
#' `gps_candidate_search()` must be handed the same frame in the same order.
#' Re-sorting or filtering in between silently repoints those indices at
#' different subjects.
#'
#' Both the row identity and the treatment column are recorded. Both parts are
#' needed: recording the row names alone misses a reorder followed by
#' `rownames(data) <- NULL`, and recording the treatment column alone misses a
#' reorder *within* one treatment group.
#'
#' Every column present when the fingerprint is taken is digested as well, so
#' that overwriting the *contents* of existing rows is caught too. Only the
#' recorded columns are checked later, which keeps attaching an outcome column
#' between stages legal.
#'
#' Row names are stored only when they are not R's automatic `1:n` sequence
#' (`.row_names_info()` reports a negative count for the compact internal
#' form), which keeps the fingerprint small for the common case.
#'
#' @param data A `data.frame`.
#' @param treatment_var Name of the treatment variable.
#' @param columns Columns to digest. Defaults to every column in `data`.
#'
#' @return A list with `n_rows`, `row_names`, `treatment`, `columns` and
#'   `digest`.
#'
#' @keywords internal
#' @noRd
data_fingerprint <- function(data, treatment_var, columns = NULL) {
  if (is.null(columns)) {
    columns <- names(data)
  }

  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop("Column(s) used to identify the original data are missing: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  automatic_row_names <- .row_names_info(data) < 0L

  list(
    n_rows = nrow(data),
    row_names = if (automatic_row_names) NULL else rownames(data),
    treatment = factor(treatment_labels(data, treatment_var)),
    columns = columns,
    digest = vapply(data[columns], column_digest, numeric(3L))
  )
}


#' Order-sensitive digest of one column
#'
#' R has no hash function in base, so each column is reduced to a
#' position-weighted triple instead. Two rows can only swap without changing
#' the weighted sum when their values are equal, in which case the swap changed
#' nothing -- so the digest is exactly as sensitive as it needs to be, at the
#' cost of one vectorised pass and 24 bytes per column.
#'
#' This is a guard against accidental reordering and rewriting, not against a
#' deliberately constructed collision.
#'
#' @param x A column.
#'
#' @return A numeric vector of length three.
#'
#' @keywords internal
#' @noRd
column_digest <- function(x) {
  values <- if (is.numeric(x) || is.logical(x)) {
    as.numeric(x)
  } else {
    # Factor codes of the sorted unique values: stable for a given column,
    # and independent of how the column happens to be stored.
    as.numeric(match(as.character(x), sort(unique(as.character(x)))))
  }

  finite <- is.finite(values)
  usable <- ifelse(finite, values, 0)

  c(
    sum(usable),
    sum(usable * seq_along(usable)),
    sum(!finite)
  )
}


#' Validate a GPS matrix against the data it describes
#'
#' Every stage that takes `gps` separately from `data` needs to know they line
#' up: the two are matched by row position, so a matrix with the wrong number of
#' rows, non-probability entries, or rows that do not sum to one produces a
#' wrong answer rather than an error.
#'
#' @param data A `data.frame`.
#' @param gps Numeric matrix of generalized propensity scores.
#' @param treatment_var Name of the treatment variable.
#' @param context Name of the calling function, used in error messages.
#'
#' @return `gps` as a numeric matrix.
#'
#' @keywords internal
#' @noRd
validate_gps <- function(data, gps, treatment_var, context) {
  if (is.data.frame(gps)) {
    gps <- as.matrix(gps)
  }
  if (!is.matrix(gps)) {
    stop("`gps` passed to ", context, "() must be a matrix or data.frame.",
         call. = FALSE)
  }
  if (nrow(gps) != nrow(data)) {
    stop("`gps` and `data` passed to ", context, "() must have the same ",
         "number of rows (", nrow(gps), " vs ", nrow(data), "). GPS values ",
         "are matched to subjects by row position.", call. = FALSE)
  }
  if (ncol(gps) == 0L) {
    stop("`gps` must have at least one treatment column.", call. = FALSE)
  }
  if (is.null(colnames(gps)) || anyDuplicated(colnames(gps)) > 0L) {
    stop("`gps` must have unique, named treatment columns.", call. = FALSE)
  }
  if (!is.numeric(gps)) {
    stop("`gps` values must be numeric.", call. = FALSE)
  }
  if (!all(is.finite(gps))) {
    stop("`gps` values must all be finite.", call. = FALSE)
  }
  if (any(gps < 0) || any(gps > 1)) {
    stop("`gps` values must lie in [0, 1].", call. = FALSE)
  }

  worst_deviation <- max(abs(rowSums(gps) - 1))
  if (worst_deviation > 1e-6) {
    stop("Each row of `gps` must sum to 1; the largest deviation is ",
         format(worst_deviation, digits = 3), ".", call. = FALSE)
  }

  # `gps` and `data` are matched by row position, never by name, so row names
  # that disagree mean one of them was reordered. Only checked when `gps`
  # carries names at all, which avoids materialising `data`'s automatic
  # sequence for callers that pass an unnamed matrix.
  if (!is.null(rownames(gps)) &&
      !identical(rownames(gps), rownames(data))) {
    stop("`gps` and `data` passed to ", context, "() must be in the same row ",
         "order; their row names differ. GPS values are matched to subjects ",
         "by row position, not by name.", call. = FALSE)
  }

  # Validate the treatment column here too, before callers derive group masks
  # from its canonical representation.
  labels <- treatment_labels(data, treatment_var)

  # A treatment level with no GPS column cannot be matched or weighted. Left
  # unchecked it simply drops out of `groups`, taking its subjects with it.
  missing_levels <- setdiff(unique(labels), colnames(gps))
  if (length(missing_levels) > 0L) {
    stop("Treatment group(s) not found in `gps`: ",
         paste(missing_levels, collapse = ", "),
         ". Every level present in `", treatment_var, "` needs a GPS column.",
         call. = FALSE)
  }

  gps
}


#' Fingerprint identifying a GPS matrix
#'
#' @param gps Numeric matrix of generalized propensity scores.
#'
#' @return A list identifying `gps`.
#'
#' @keywords internal
#' @noRd
gps_fingerprint <- function(gps) {
  gps <- as.matrix(gps)

  list(
    n_rows = nrow(gps),
    columns = colnames(gps),
    digest = vapply(seq_len(ncol(gps)),
                    function(j) column_digest(gps[, j]), numeric(3L))
  )
}


#' Verify that a later stage received the GPS the search was built from
#'
#' @param search Output from [gps_candidate_search()].
#' @param gps Numeric matrix of generalized propensity scores.
#'
#' @return `TRUE`, invisibly.
#'
#' @keywords internal
#' @noRd
check_gps_fingerprint <- function(search, gps) {
  recorded <- search$gps_fingerprint
  if (is.null(recorded)) {
    return(invisible(TRUE))
  }

  if (!identical(recorded, gps_fingerprint(gps))) {
    stop("`gps` is not the matrix `search` was built from. Use the same GPS ",
         "matrix, with the same rows, columns and values, at every stage: ",
         "candidate pools and matched sets were selected against the one the ",
         "search saw.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Verify that `data` is the frame a candidate search ran on
#'
#' A `search` without a fingerprint is accepted, so objects created by earlier
#' versions of the package keep working.
#'
#' @param search Output from [gps_candidate_search()].
#' @param data A `data.frame`.
#' @param treatment_var Name of the treatment variable.
#'
#' @return `TRUE`, invisibly.
#'
#' @keywords internal
#' @noRd
check_fingerprint <- function(search, data, treatment_var) {
  recorded <- search$fingerprint
  if (is.null(recorded)) {
    return(invisible(TRUE))
  }

  # Objects from earlier versions recorded neither the column list nor the
  # digests; comparing only what they hold keeps them usable.
  recorded_columns <- recorded$columns
  current <- data_fingerprint(data, treatment_var, columns = recorded_columns)

  mismatch <- if (!identical(recorded$n_rows, current$n_rows)) {
    paste0("row count changed (", recorded$n_rows, " -> ",
           current$n_rows, ")")
  } else if (!identical(recorded$row_names, current$row_names)) {
    "row names changed"
  } else if (!identical(as.character(recorded$treatment),
                        as.character(current$treatment))) {
    "the treatment column changed"
  } else if (!is.null(recorded$digest) &&
             !isTRUE(all.equal(recorded$digest, current$digest))) {
    changed <- recorded_columns[
      !vapply(seq_along(recorded_columns), function(j) {
        isTRUE(all.equal(recorded$digest[, j], current$digest[, j]))
      }, logical(1))
    ]
    paste0("column value(s) changed: ", paste(changed, collapse = ", "))
  } else {
    NULL
  }

  if (!is.null(mismatch)) {
    stop("`data` is not the frame `search` was built from: ", mismatch, ". ",
         "Matched sets are stored as positional row indices, so `data` must ",
         "be passed to every stage in the same row order that ",
         "gps_candidate_search() saw. Adding a column is fine; re-sorting or ",
         "filtering is not.", call. = FALSE)
  }

  invisible(TRUE)
}


#' The stack a golden fixture was generated on
#'
#' Package versions alone do not identify it. The same versions on a different
#' platform or BLAS move an iterative solver's last digits, which is exactly
#' what the strict golden tier compares, so the platform and the linear algebra
#' libraries are recorded too.
#'
#' Shared by `tests/testthat/golden/make-golden.R` and `test-golden.R` so the
#' record and the check cannot drift apart.
#'
#' @return A one-row `data.frame`.
#'
#' @keywords internal
#' @noRd
golden_environment <- function() {
  data.frame(
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    blas = basename(extSoftVersion()[["BLAS"]]),
    lapack = basename(La_library()),
    nnet = as.character(utils::packageVersion("nnet")),
    MASS = as.character(utils::packageVersion("MASS")),
    RANN = as.character(utils::packageVersion("RANN")),
    stringsAsFactors = FALSE
  )
}
