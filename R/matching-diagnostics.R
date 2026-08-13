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
#'     \item{by_covariate}{SMD for each covariate and comparator group.}
#'     \item{summary}{Mean and maximum absolute SMD for each comparator group.}
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
  rows <- list()
  for (g in groups) {
    x_anchor <- data[matched$anchor, X_vars, drop = FALSE]
    x_g <- data[matched[[g]], X_vars, drop = FALSE]
    for (v in X_vars) {
      mean_a <- mean(x_anchor[[v]]); mean_g <- mean(x_g[[v]])
      var_a <- stats::var(x_anchor[[v]]); var_g <- stats::var(x_g[[v]])
      smd <- (mean_a - mean_g) / sqrt((var_a + var_g) / 2)
      rows[[length(rows) + 1]] <- data.frame(group = g, covariate = v, smd = smd)
    }
  }
  by_covariate <- do.call(rbind, rows)
  rownames(by_covariate) <- NULL

  summary_rows <- lapply(groups, function(g) {
    vals <- abs(by_covariate$smd[by_covariate$group == g])
    data.frame(group = g, mean_abs_smd = mean(vals), max_abs_smd = max(vals))
  })
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  list(by_covariate = by_covariate, summary = summary_df)
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
