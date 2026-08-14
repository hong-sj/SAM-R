################################################################################
### Comparator methods: Stabilized IPTW, Overlap Weighting and Matching Weights
## 
## Author: Sungjun Hong
################################################################################

#' Generalized balancing weights for multiple treatment groups
#'
#' @description
#' Computes IPTW, overlap weights, or matching weights from generalized
#' propensity scores (GPS). If `gps` is not supplied, it is estimated with
#' [estimate_gps_multinom()].
#'
#' @param data A `data.frame` containing `treatment_var` (and `X_vars`, if
#'   `gps` is not supplied).
#' @param method One of `"iptw"`, `"overlap"`, `"matching"`.
#' @param gps Optional pre-computed GPS matrix (e.g.
#'   `estimate_gps_multinom(data)$gps`); estimated internally if `NULL`.
#' @param X_vars,treatment_var,anchor_level See [estimate_gps_multinom()].
#'   Only used when `gps` is `NULL` (`anchor_level` is passed through to
#'   the GPS estimation call for its reference-category convention).
#' @param stabilize Logical, whether to multiply IPTW weights by the
#'   empirical marginal treatment prevalence (Robins et al., 2000).
#'   Ignored for `method %in% c("overlap", "matching")`, which are already
#'   bounded and do not use this stabilization (see file header). Default
#'   `TRUE`.
#'
#' @return A list with elements `weights` (numeric vector, length
#'   `nrow(data)`), `h` (the tilting-function value per subject), `method`,
#'   and `gps` (the GPS matrix used).
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
#' w_ov <- compute_balancing_weights(
#'   sample_4group,
#'   method = "overlap",
#'   X_vars = covariates,
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#' summary(w_ov$weights)
#'
#' @export
compute_balancing_weights <- function(data, method = c("iptw", "overlap", "matching"),
                                       gps = NULL, X_vars = paste0("X", 1:10),
                                       treatment_var = "T", anchor_level = "A",
                                       stabilize = TRUE) {
  method <- match.arg(method)
  treat_chr <- treatment_labels(data, treatment_var)
  anchor_level <- treatment_level(anchor_level)
  if (is.null(gps)) {
    gps <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = treatment_var,
                                  anchor_level = anchor_level)$gps
  }
  stopifnot(all(treat_chr %in% colnames(gps)))

  own_gps <- gps[cbind(seq_len(nrow(gps)), match(treat_chr, colnames(gps)))]

  h <- switch(method,
    iptw = rep(1, nrow(gps)),
    overlap = 1 / rowSums(1 / gps),
    matching = apply(gps, 1, min)
  )

  weights <- h / own_gps

  if (method == "iptw" && isTRUE(stabilize)) {
    marg_p <- prop.table(table(treat_chr))
    weights <- weights * as.numeric(marg_p[treat_chr])
  }

  list(weights = as.numeric(weights), h = as.numeric(h), method = method, gps = gps)
}

#' Weighted covariate balance
#'
#' @description
#' Computes weighted standardized mean differences (SMDs) between the anchor
#' group and each comparator group for all covariates in `X_vars`.
#'
#' SMDs are reported as anchor minus comparator, the same orientation as
#' [compute_smd_balance()], so the two tables can be read side by side.
#'
#' @param data,X_vars,treatment_var,anchor_level See
#'   [compute_balancing_weights()].
#' @param weights Numeric vector, length `nrow(data)`, e.g.
#'   `compute_balancing_weights(...)$weights`.
#'
#' @return A list with `by_covariate` and `summary` data frames, in the
#'   same shape as [compute_smd_balance()]'s return value. `by_covariate`
#'   additionally carries `abs_smd`, and both carry the `smd_defined` flag and
#'   `n_undefined` count described there.
#'
#' @export
compute_weighted_balance <- function(data, weights, X_vars = paste0("X", 1:10),
                                      treatment_var = "T", anchor_level = "A") {
  stopifnot(length(weights) == nrow(data))
  treat_chr <- treatment_labels(data, treatment_var)
  anchor_level <- treatment_level(anchor_level)
  require_rows(which(treat_chr == anchor_level), anchor_level, treatment_var,
               available = treat_chr)
  groups <- setdiff(unique(treat_chr), anchor_level)

  wtd_mean <- function(x, w) sum(x * w) / sum(w)
  wtd_var <- function(x, w) {
    m <- wtd_mean(x, w)
    sum(w * (x - m)^2) / sum(w)
  }

  is_anchor <- treat_chr == anchor_level
  w_anchor <- weights[is_anchor]
  x_anchor <- data[is_anchor, X_vars, drop = FALSE]

  rows <- lapply(groups, function(g) {
    is_g <- treat_chr == g
    w_g <- weights[is_g]
    x_g <- data[is_g, X_vars, drop = FALSE]
    do.call(rbind, lapply(X_vars, function(v) {
      m_a <- wtd_mean(x_anchor[[v]], w_anchor); m_g <- wtd_mean(x_g[[v]], w_g)
      var_a <- wtd_var(x_anchor[[v]], w_anchor); var_g <- wtd_var(x_g[[v]], w_g)

      # Anchor minus group, matching compute_smd_balance(). Same guard as
      # there: see summarize_smd().
      pooled_sd <- sqrt((var_a + var_g) / 2)
      defined <- isTRUE(is.finite(pooled_sd) && pooled_sd > 0)
      smd <- if (defined) (m_a - m_g) / pooled_sd else 0

      data.frame(group = g, covariate = v, smd = smd, abs_smd = abs(smd),
                 smd_defined = defined)
    }))
  })
  by_covariate <- do.call(rbind, rows)
  rownames(by_covariate) <- NULL

  summary_df <- summarize_smd(by_covariate, groups)

  list(by_covariate = by_covariate, summary = summary_df)
}

#' Effective sample size and weight-distribution diagnostics, per group
#'
#' @description
#' Computes Kish's effective sample size (ESS) and simple weight-distribution
#' diagnostics separately for each treatment group.
#'
#' @param data,treatment_var See [compute_balancing_weights()].
#' @param weights Numeric vector, length `nrow(data)`.
#'
#' @return A `data.frame`, one row per treatment level, with columns
#'   `group`, `n`, `ess`, `mean_weight`, `max_weight`, `max_over_mean`.
#'
#' @export
compute_effective_sample_size <- function(data, weights, treatment_var = "T") {
  stopifnot(length(weights) == nrow(data))
  treat_chr <- treatment_labels(data, treatment_var)
  groups <- unique(treat_chr)
  rows <- lapply(groups, function(g) {
    w <- weights[treat_chr == g]
    data.frame(
      group = g, n = length(w), ess = sum(w)^2 / sum(w^2),
      mean_weight = mean(w), max_weight = max(w), max_over_mean = max(w) / mean(w)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Estimate + evaluate one weighting comparator method in a single call
#'
#' @description
#' Convenience wrapper that computes balancing weights, weighted covariate
#' balance, and effective sample size in one call.
#'
#' @inheritParams compute_balancing_weights
#'
#' @return A list with elements `weights` (from
#'   [compute_balancing_weights()]), `balance` (from
#'   [compute_weighted_balance()]), and `ess` (from
#'   [compute_effective_sample_size()]).
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
#' report <- evaluate_comparator_weighting(
#'   sample_4group,
#'   method = "matching",
#'   X_vars = covariates,
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#' report$balance$summary
#' report$ess
#'
#' @export
evaluate_comparator_weighting <- function(data, method = c("iptw", "overlap", "matching"),
                                           gps = NULL, X_vars = paste0("X", 1:10),
                                           treatment_var = "T", anchor_level = "A",
                                           stabilize = TRUE) {
  method <- match.arg(method)
  w <- compute_balancing_weights(data, method = method, gps = gps, X_vars = X_vars,
                                  treatment_var = treatment_var, anchor_level = anchor_level,
                                  stabilize = stabilize)
  balance <- compute_weighted_balance(data, w$weights, X_vars = X_vars,
                                       treatment_var = treatment_var, anchor_level = anchor_level)
  ess <- compute_effective_sample_size(data, w$weights, treatment_var = treatment_var)
  list(weights = w, balance = balance, ess = ess)
}
