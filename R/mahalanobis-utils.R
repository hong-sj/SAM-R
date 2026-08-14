################################################################################
### Mahalanobis Distance Utilities
##
## Author: Sungjun Hong
################################################################################

#' Pooled within-group covariance matrix
#'
#' @description
#' Computes the pooled within-group covariance matrix across treatment groups
#' and its inverse for Mahalanobis distance calculations.
#'
#' @param data A `data.frame` containing covariate and treatment variables.
#' @param X_vars Character vector of covariate column names.
#' @param treatment_var Name of the treatment variable.
#'
#' @return A list containing the pooled covariance matrix (`S`) and its
#'   inverse (`S_inv`).
#'
#' @examples
#' data(sample_4group)
#'
#' covariates <- setdiff(
#'   names(sample_4group),
#'   c("synthetic_id", "treatment", "mortality_28d")
#' )
#'
#' pooled <- get_pooled_covariance(
#'   sample_4group,
#'   X_vars = covariates,
#'   treatment_var = "treatment"
#' )
#'
#' dim(pooled$S)
#'
#' @export
get_pooled_covariance <- function(data, X_vars, treatment_var) {
  labels <- treatment_labels(data, treatment_var)
  X <- covariate_matrix(data, X_vars)
  groups <- unique(labels)
  p <- length(X_vars)
  S_within <- matrix(0, p, p, dimnames = list(X_vars, X_vars))

  for (g in groups) {
    Xg <- X[labels == g, , drop = FALSE]
    Xg_centered <- sweep(Xg, MARGIN = 2, STATS = colMeans(Xg), FUN = "-")
    S_within <- S_within + crossprod(Xg_centered)
  }

  df <- nrow(X) - length(groups)
  if (df <= 0L) {
    stop("Residual degrees of freedom for the pooled covariance are ", df,
         " (", nrow(X), " rows, ", length(groups), " treatment groups). ",
         "At least one more row than treatment groups is required.",
         call. = FALSE)
  }
  S <- S_within / df

  # solve() does not raise on a non-finite matrix -- it returns an all-NA one --
  # so the result is checked explicitly. Without that check the fallback below
  # is dead code for missing data, and MASS::ginv() covers genuine, finite rank
  # deficiency only.
  S_inv <- tryCatch(solve(S), error = function(e) NULL)

  if (is.null(S_inv) || !all(is.finite(S_inv))) {
    warning("Pooled covariance matrix is numerically singular; using MASS::ginv().")
    S_inv <- MASS::ginv(S)
  }

  list(S = S, S_inv = S_inv)
}

#' Pairwise Mahalanobis distance matrix
#'
#' @description
#' Computes pairwise Mahalanobis distances between rows of `X_query`
#' and `X_reference` using the supplied inverse covariance matrix.
#'
#' @param X_query Numeric matrix with observations in rows.
#' @param X_reference Numeric matrix with observations in rows.
#' @param S_inv Inverse covariance matrix.
#'
#' @return A numeric matrix of pairwise Mahalanobis distances.
#'
#' @examples
#' data(sample_4group)
#'
#' covariates <- setdiff(
#'   names(sample_4group),
#'   c("synthetic_id", "treatment", "mortality_28d")
#' )
#' groups <- unique(as.character(sample_4group$treatment))
#'
#' pooled <- get_pooled_covariance(
#'   sample_4group,
#'   covariates,
#'   "treatment"
#' )
#'
#' X1 <- as.matrix(
#'   sample_4group[sample_4group$treatment == groups[1],
#'                 covariates, drop = FALSE]
#' )
#' X2 <- as.matrix(
#'   sample_4group[sample_4group$treatment == groups[2],
#'                 covariates, drop = FALSE]
#' )
#'
#' mahalanobis_distance_matrix(
#'   X1[1:5, , drop = FALSE],
#'   X2[1:5, , drop = FALSE],
#'   pooled$S_inv
#' )
#' @export
mahalanobis_distance_matrix <- function(X_query, X_reference, S_inv) {

  # Remove dimension names from distance calculations
  X_query <- unname(X_query)
  X_reference <- unname(X_reference)

  A <- X_query %*% S_inv
  B <- X_reference %*% S_inv
  x_sq <- rowSums(A * X_query)
  y_sq <- rowSums(B * X_reference)
  cross <- A %*% t(X_reference)

  d2 <- outer(x_sq, y_sq, "+") - 2 * cross
  
  # Guard against small negative values from floating-point error
  d2[d2 < 0] <- 0
  unname(sqrt(d2))
}

#' Build Mahalanobis distance matrices by comparator group
#'
#' @description
#' Computes Mahalanobis distance matrices between anchor subjects and
#' subjects in each comparator treatment group.
#'
#' @param data A `data.frame` containing covariate and treatment variables.
#' @param X_vars Character vector of covariate column names.
#' @param treatment_var Name of the treatment variable.
#' @param anchor_rows Row indices of anchor subjects.
#' @param groups Character vector of comparator treatment groups.
#'
#' @return A list containing the inverse covariance matrix (`S_inv`),
#'   comparator-group row indices (`group_rows`), and Mahalanobis distance
#'   matrices (`D`).
#'
#' @examples
#' data(sample_4group)
#'
#' covariates <- setdiff(
#'   names(sample_4group),
#'   c("synthetic_id", "treatment", "mortality_28d")
#' )
#'
#' treatment_levels <- unique(
#'   as.character(sample_4group$treatment)
#' )
#' anchor <- treatment_levels[1]
#' groups <- setdiff(treatment_levels, anchor)
#'
#' anchor_rows <- which(
#'   sample_4group$treatment == anchor
#' )
#'
#' built <- build_group_distance_matrices(
#'   sample_4group,
#'   X_vars = covariates,
#'   treatment_var = "treatment",
#'   anchor_rows = anchor_rows,
#'   groups = groups
#' )
#'
#' dim(built$D[[groups[1]]])
#'
#' @export
build_group_distance_matrices <- function(data, X_vars, treatment_var, anchor_rows, groups) {
  pooled <- get_pooled_covariance(data, X_vars, treatment_var)
  labels <- treatment_labels(data, treatment_var)
  groups <- as.character(groups)
  group_rows <- stats::setNames(
    lapply(groups, function(g) {
      require_rows(which(labels == g), g, treatment_var, available = labels)
    }), groups
  )
  X <- covariate_matrix(data, X_vars)
  X_anchor <- X[anchor_rows, , drop = FALSE]
  D <- stats::setNames(lapply(groups, function(g) {
    mahalanobis_distance_matrix(X_anchor, X[group_rows[[g]], , drop = FALSE],
                                pooled$S_inv)
  }), groups)
  list(S_inv = pooled$S_inv, group_rows = group_rows, D = D)
}
