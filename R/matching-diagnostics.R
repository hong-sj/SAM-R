################################################################################
### Matching Diagnostics
##
## Author: Sungjun Hong
################################################################################


#' Standardized mean difference balance
#'
#' @description
#' Computes standardized mean differences (SMDs) between the matched anchor
#' group and each comparator group for all specified covariates.
#'
#' @param data A `data.frame` containing the covariates.
#' @param matched Matched-set data returned by [sam_match()] or
#'   [match_3way()].
#' @param X_vars Character vector of covariate column names.
#' @param groups Character vector of comparator treatment groups.
#'
#' @return A list containing:
#'   \describe{
#'     \item{by_covariate}{SMD for each covariate and comparator group, with an
#'       `smd_defined` flag that is `FALSE` when a covariate has no variance in
#'       either arm and its SMD is therefore undefined.}
#'     \item{summary}{Mean and maximum absolute SMD for each comparator group,
#'       computed over the assessable covariates only, plus `n_undefined`.}
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
#' balance <- compute_smd_balance(
#'   sample_4group,
#'   matched$matched,
#'   X_vars = covariates,
#'   groups = search$groups
#' )
#'
#' balance$summary
#'
#' @export
compute_smd_balance <- function(data, matched, X_vars, groups) {
  require_covariates(data, X_vars)
  groups <- as.character(groups)

  x_anchor <- data[matched$anchor, X_vars, drop = FALSE]

  rows <- list()
  for (g in groups) {
    x_g <- data[matched[[g]], X_vars, drop = FALSE]
    for (v in X_vars) {
      mean_a <- mean(x_anchor[[v]]); mean_g <- mean(x_g[[v]])
      var_a <- stats::var(x_anchor[[v]]); var_g <- stats::var(x_g[[v]])

      # A covariate with no variance in either arm has no defined SMD.
      # Reporting 0 keeps it out of the summary statistics without emitting a
      # NaN, which propagates through mean() and max() and would otherwise
      # take the whole summary down with it.
      pooled_sd <- sqrt((var_a + var_g) / 2)
      defined <- isTRUE(is.finite(pooled_sd) && pooled_sd > 0)
      smd <- if (defined) (mean_a - mean_g) / pooled_sd else 0

      rows[[length(rows) + 1]] <- data.frame(
        group = g, covariate = v, smd = smd, smd_defined = defined
      )
    }
  }
  by_covariate <- do.call(rbind, rows)
  rownames(by_covariate) <- NULL

  summary_df <- summarize_smd(by_covariate, groups)

  list(by_covariate = by_covariate, summary = summary_df)
}


#' Summarize a per-covariate SMD table
#'
#' Shared by [compute_smd_balance()] and [compute_weighted_balance()] so that
#' both report the same columns under the same rules. Covariates that could not
#' be assessed are excluded from the mean and the maximum and counted in
#' `n_undefined`, so a degenerate covariate is visible rather than either
#' hidden or fatal.
#'
#' @param by_covariate A data frame with `group`, `smd` and `smd_defined`.
#' @param groups Character vector of comparator groups, in output order.
#'
#' @return A data frame with one row per group.
#'
#' @keywords internal
#' @noRd
summarize_smd <- function(by_covariate, groups) {
  summary_rows <- lapply(groups, function(g) {
    in_group <- by_covariate$group == g
    values <- abs(by_covariate$smd[in_group & by_covariate$smd_defined])

    data.frame(
      group = g,
      mean_abs_smd = if (length(values) > 0L) mean(values) else NA_real_,
      max_abs_smd = if (length(values) > 0L) max(values) else NA_real_,
      n_undefined = sum(in_group & !by_covariate$smd_defined)
    )
  })
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_df
}

#' Pairwise treatment-discrimination AUC
#'
#' @description
#' Computes pairwise treatment-discrimination AUCs among matched treatment
#' groups using contrasts of generalized propensity scores. Values close to
#' 0.5 indicate limited discrimination between treatment groups after matching.
#'
#' @param gps Numeric matrix of generalized propensity scores.
#' @param matched Matched-set data returned by [sam_match()] or
#'   [match_3way()].
#' @param groups Character vector of comparator treatment groups.
#' @param anchor_level Anchor treatment group.
#'
#' @return A list containing:
#'   \describe{
#'     \item{pairwise}{Pairwise AUC for each treatment-group comparison.}
#'     \item{mean_auc}{Mean pairwise AUC.}
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
#' auc <- compute_pairwise_treatment_auc(
#'   fit$gps,
#'   matched$matched,
#'   groups = search$groups,
#'   anchor_level = anchor
#' )
#'
#' auc$pairwise
#' auc$mean_auc
#'
#' @export
compute_pairwise_treatment_auc <- function(gps, matched, groups, anchor_level) {
  all_levels <- c(anchor_level, groups)
  rows_by_level <- stats::setNames(
    lapply(all_levels, function(g) if (g == anchor_level) matched$anchor else matched[[g]]),
    all_levels
  )

  pairs <- utils::combn(all_levels, 2, simplify = FALSE)
  pairwise_rows <- lapply(pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    score <- log(gps[, g1]) - log(gps[, g2])
    rows_g1 <- rows_by_level[[g1]]; rows_g2 <- rows_by_level[[g2]]
    combined_rows <- c(rows_g1, rows_g2)
    label <- c(rep(1L, length(rows_g1)), rep(0L, length(rows_g2)))
    auc <- auc_mannwhitney(score[combined_rows], label)
    data.frame(group_1 = g1, group_2 = g2, auc = auc)
  })
  pairwise <- do.call(rbind, pairwise_rows)
  rownames(pairwise) <- NULL

  list(pairwise = pairwise, mean_auc = mean(pairwise$auc))
}
