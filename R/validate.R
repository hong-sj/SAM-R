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

  as.character(data[[treatment_var]])
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
#' reorder *within* one treatment group. No other column is recorded, so
#' attaching an outcome column between stages stays legal.
#'
#' Row names are stored only when they are not R's automatic `1:n` sequence
#' (`.row_names_info()` reports a negative count for the compact internal
#' form), which keeps the fingerprint small for the common case.
#'
#' @param data A `data.frame`.
#' @param treatment_var Name of the treatment variable.
#'
#' @return A list with `n_rows`, `row_names` and `treatment`.
#'
#' @keywords internal
#' @noRd
data_fingerprint <- function(data, treatment_var) {
  automatic_row_names <- .row_names_info(data) < 0L

  list(
    n_rows = nrow(data),
    row_names = if (automatic_row_names) NULL else rownames(data),
    treatment = factor(treatment_labels(data, treatment_var))
  )
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

  current <- data_fingerprint(data, treatment_var)

  mismatch <- if (!identical(recorded$n_rows, current$n_rows)) {
    paste0("row count changed (", recorded$n_rows, " -> ",
           current$n_rows, ")")
  } else if (!identical(recorded$row_names, current$row_names)) {
    "row names changed"
  } else if (!identical(as.character(recorded$treatment),
                        as.character(current$treatment))) {
    "the treatment column changed"
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
