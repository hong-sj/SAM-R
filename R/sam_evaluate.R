################################################################################
### SAM Matching Evaluation and Outcome Analysis
##
## Author: Sungjun Hong
################################################################################


#' Evaluate a completed SAM match
#'
#' @description
#' Evaluates a completed SAM match using matched-set distance distributions,
#' covariate balance, treatment-discrimination AUC, and the matching rate.
#'
#' @param data Original `data.frame` used for matching.
#' @param search Output from [gps_candidate_search()].
#' @param match_result Output from [sam_match()].
#' @param gps Numeric matrix of generalized propensity scores.
#' @param X_vars Character vector of covariate column names.
#'   Default `paste0("X", 1:10)`.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{loss_distribution}{Summary of matched-set losses.}
#'     \item{dispersion_distribution}{Summary of Mahalanobis dispersion.}
#'     \item{smd_balance}{Covariate balance based on standardized mean
#'       differences.}
#'     \item{treatment_discrimination_auc}{Pairwise treatment-discrimination
#'       AUCs.}
#'     \item{matching_rate}{Proportion of anchor subjects successfully
#'       matched.}
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
#' match_result <- sam_match(
#'   sample_4group,
#'   search,
#'   X_vars = covariates,
#'   treatment_var = "treatment"
#' )
#'
#' report <- sam_evaluate(
#'   sample_4group,
#'   search,
#'   match_result,
#'   fit$gps,
#'   X_vars = covariates,
#'   treatment_var = "treatment"
#' )
#'
#' report$loss_distribution
#' report$smd_balance$summary
#' report$treatment_discrimination_auc$mean_auc
#'
#' @export
sam_evaluate <- function(data, search, match_result, gps,
                          X_vars = paste0("X", 1:10), treatment_var = "T") {
  groups <- as.character(search$groups)
  matched <- match_result$matched

  check_fingerprint(search, data, treatment_var)
  gps <- validate_gps(data, gps, treatment_var, "sam_evaluate")
  check_gps_fingerprint(search, gps)
  labels <- treatment_labels(data, treatment_var)
  anchor_level <- unique(labels[search$anchor_rows])
  stopifnot(length(anchor_level) == 1)

  # Summarize a matched-set distance metric
  loss_dist_of <- function(x) {
    if (length(x) == 0) {
      return(data.frame(mean = NA_real_, median = NA_real_, sd = NA_real_,
                         p95 = NA_real_, max = NA_real_))
    }
    data.frame(
      mean = mean(x), median = stats::median(x), sd = stats::sd(x),
      p95 = stats::quantile(x, 0.95, names = FALSE), max = max(x)
    )
  }

  # Matched-set loss distribution
  loss_distribution <- loss_dist_of(matched$loss)

  # Sum of squared anchor-to-comparator Mahalanobis distances
  dispersion_vals <- if (nrow(matched) > 0) {
    rowSums(as.matrix(matched[, paste0("dist_", groups), drop = FALSE])^2)
  } else {
    numeric(0)
  }
  dispersion_distribution <- loss_dist_of(dispersion_vals)

  # Covariate balance
  smd_balance <- compute_smd_balance(data, matched, X_vars, groups)

  # Pairwise treatment-discrimination AUC
  treatment_discrimination_auc <- compute_pairwise_treatment_auc(gps, matched, groups, anchor_level)

  matching_rate <- match_result$matching_rate

  list(
    loss_distribution = loss_distribution,
    dispersion_distribution = dispersion_distribution,
    smd_balance = smd_balance,
    treatment_discrimination_auc = treatment_discrimination_auc,
    matching_rate = matching_rate
  )
}

#' Extract matched subject-level data
#'
#' @description
#' Reconstructs the subject-level matched cohort from a completed SAM match.
#' Each retained subject is assigned its matched-set identifier, treatment
#' role, and original row index.
#'
#' @param data Original `data.frame` used for matching.
#' @param search Output from [gps_candidate_search()].
#' @param match_result Output from [sam_match()].
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#' @param anchor_level Optional anchor treatment group. If `NULL`, it is
#'   inferred from `search`.
#'
#' @return A subject-level `data.frame` containing the matched cohort with
#'   `matched_set_id`, `matched_role`, and `original_row` appended.
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
#' match_result <- sam_match(
#'   sample_4group,
#'   search,
#'   X_vars = covariates,
#'   treatment_var = "treatment"
#' )
#'
#' matched_data <- extract_matched_data(
#'   sample_4group,
#'   search,
#'   match_result,
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#'
#' head(matched_data)
#'
#' @export
extract_matched_data <- function(data, search, match_result,
                                     treatment_var = "T",
                                     anchor_level = NULL) {
  matched <- match_result$matched
  groups <- as.character(search$groups)

  check_fingerprint(search, data, treatment_var)
  labels <- treatment_labels(data, treatment_var)
  anchor_level <- treatment_level(anchor_level)

  if (is.null(matched) || nrow(matched) == 0L) {
    stop("`match_result` contains no matched sets.")
  }

  required_match_cols <- c("anchor", groups)
  missing_cols <- setdiff(required_match_cols, names(matched))
  if (length(missing_cols) > 0L) {
    stop("Missing matched-set column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  if (is.null(anchor_level)) {
    anchor_level <- unique(labels[search$anchor_rows])
    if (length(anchor_level) != 1L) {
      stop("Could not uniquely determine `anchor_level` from ",
           "`search$anchor_rows`.")
    }
  }

  # Support legacy matched-set identifiers
  if ("matched_set_id" %in% names(matched)) {
    set_ids <- matched$matched_set_id
  } else if ("quartet_id" %in% names(matched)) {
    warning("Legacy `quartet_id` detected; treating it as `matched_set_id`.")
    set_ids <- matched$quartet_id
  } else {
    set_ids <- seq_len(nrow(matched))
  }

  # Expand each K-way matched set into K subject-level records
  matched_index <- do.call(
    rbind,
    lapply(seq_len(nrow(matched)), function(i) {
      subject_rows <- as.integer(
        unlist(
          matched[i, c("anchor", groups), drop = FALSE],
          use.names = FALSE
        )
      )
      data.frame(
        matched_set_id = set_ids[i],
        matched_role = c(anchor_level, groups),
        original_row = subject_rows,
        stringsAsFactors = FALSE
      )
    })
  )

  if (anyNA(matched_index$original_row) ||
      any(matched_index$original_row < 1L |
          matched_index$original_row > nrow(data))) {
    stop("Invalid row index detected in `match_result$matched`.")
  }

  matched_data <- data[matched_index$original_row, , drop = FALSE]
  matched_data$matched_set_id <- matched_index$matched_set_id
  matched_data$matched_role <- matched_index$matched_role
  matched_data$original_row <- matched_index$original_row

  treatment_order <- c(anchor_level, groups)
  matched_data <- matched_data[
    order(
      matched_data$matched_set_id,
      match(labels[matched_index$original_row], treatment_order)
    ),
    , drop = FALSE
  ]
  rownames(matched_data) <- NULL

  # Integrity check:
  # exactly one subject from every treatment group must appear in every matched set
  K <- length(treatment_order)
  set_sizes <- table(matched_data$matched_set_id)
  if (any(set_sizes != K)) {
    stop("At least one matched set does not contain exactly ", K,
         " subjects.")
  }

  group_counts_by_set <- table(
    matched_data$matched_set_id,
    as.character(matched_data[[treatment_var]])
  )
  missing_groups <- setdiff(treatment_order, colnames(group_counts_by_set))
  if (length(missing_groups) > 0L) {
    stop("Matched data are missing treatment group(s): ",
         paste(missing_groups, collapse = ", "))
  }
  group_counts_by_set <- group_counts_by_set[
    , treatment_order, drop = FALSE
  ]
  if (any(group_counts_by_set != 1L)) {
    stop("Each matched set must contain exactly one subject from every ",
         "treatment group.")
  }

  attr(matched_data, "anchor_level") <- anchor_level
  attr(matched_data, "groups") <- groups
  attr(matched_data, "K") <- K

  matched_data
}


# Matched-set cluster-robust covariance matrix with finite-sample correction.
.sam_cluster_vcov <- function(fit, cluster) {
  if (!inherits(fit, "glm")) {
    stop("`fit` must be a fitted glm object.")
  }
  if (length(cluster) != stats::nobs(fit)) {
    stop("`cluster` length must equal the number of observations used ",
         "in the model.")
  }

  X <- stats::model.matrix(fit)
  mf <- stats::model.frame(fit)
  y <- stats::model.response(mf)
  mu <- stats::fitted(fit)

  W <- mu * (1 - mu)
  XtWX <- crossprod(X, X * W)
  bread <- tryCatch(solve(XtWX), error = function(e) NULL)
  if (is.null(bread)) {
    stop("Could not invert the logistic-regression information matrix.")
  }

  # Subject-level score contributions, summed within matched sets.
  score_i <- X * (y - mu)
  score_cluster <- rowsum(score_i, factor(cluster), reorder = FALSE)
  meat <- crossprod(score_cluster)
  vcov_cr <- bread %*% meat %*% bread

  G <- nrow(score_cluster)
  N <- nrow(X)
  P <- ncol(X)
  if (G > 1L && N > P) {
    correction <- (G / (G - 1)) * ((N - 1) / (N - P))
    vcov_cr <- correction * vcov_cr
  }

  dimnames(vcov_cr) <- list(colnames(X), colnames(X))
  vcov_cr
}


# Numerical gradient for delta-method inference
.sam_numeric_gradient <- function(fun, beta, rel_step = 1e-6) {
  p <- length(beta)
  grad <- numeric(p)

  for (j in seq_len(p)) {
    h <- rel_step * max(1, abs(beta[j]))
    beta_hi <- beta
    beta_lo <- beta
    beta_hi[j] <- beta_hi[j] + h
    beta_lo[j] <- beta_lo[j] - h
    grad[j] <- (fun(beta_hi) - fun(beta_lo)) / (2 * h)
  }

  names(grad) <- names(beta)
  grad
}


#' Estimate treatment effects after SAM matching
#'
#' @description
#' Estimates treatment effects in the matched cohort using an unadjusted
#' logistic regression of the binary outcome on treatment. Confidence
#' intervals use a matched-set cluster-robust covariance matrix.
#'
#' Reports observed treatment-group risks and comparator-versus-anchor
#' odds ratios (OR), risk ratios (RR), and risk differences (RD).
#'
#' @param matched_data Subject-level matched data returned by
#'   [extract_matched_data()].
#' @param outcome_var Name of the binary outcome variable.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#' @param set_id_var Name of the matched-set identifier.
#'   Default `"matched_set_id"`.
#' @param anchor_level Optional anchor treatment group. If `NULL`, it is
#'   obtained from `matched_data`.
#' @param conf_level Confidence level. Default `0.95`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{analysis_summary}{Matched-cohort summary.}
#'     \item{group_risk}{Observed risk estimates by treatment group.}
#'     \item{contrasts}{Comparator-versus-anchor OR, RR, and RD estimates
#'       with confidence intervals, and a `separation` flag marking rows whose
#'       odds ratio is unreliable because one of the two arms has no events or
#'       no non-events.}
#'     \item{model}{Fitted logistic regression model.}
#'     \item{vcov_cluster}{Matched-set cluster-robust covariance matrix.}
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
#' match_result <- sam_match(
#'   sample_4group,
#'   search,
#'   X_vars = covariates,
#'   treatment_var = "treatment"
#' )
#'
#' matched_data <- extract_matched_data(
#'   sample_4group,
#'   search,
#'   match_result,
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#'
#' effects <- sam_estimate_effects(
#'   matched_data,
#'   outcome_var = "mortality_28d",
#'   treatment_var = "treatment",
#'   anchor_level = anchor
#' )
#'
#' effects$group_risk
#' effects$contrasts
#'
#' @export
sam_estimate_effects <- function(matched_data, outcome_var,
                                 treatment_var = "T",
                                 set_id_var = "matched_set_id",
                                 anchor_level = NULL,
                                 conf_level = 0.95) {
  if (!is.data.frame(matched_data)) {
    stop("`matched_data` must be a data.frame.")
  }

  required_cols <- c(outcome_var, treatment_var, set_id_var)
  missing_cols <- setdiff(required_cols, names(matched_data))
  if (length(missing_cols) > 0L) {
    stop("Missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number strictly between 0 and 1.")
  }

  y <- matched_data[[outcome_var]]
  if (is.logical(y)) y <- as.integer(y)
  if (!all(stats::na.omit(y) %in% c(0, 1))) {
    stop("`outcome_var` must be binary and coded as 0/1.")
  }

  # Restrict the analysis to complete matched *sets*. Dropping individual rows
  # would leave a set with fewer than K members in the analysis, which breaks
  # the 1:1:...:1 structure that the matched-set cluster-robust variance and
  # the group comparison both assume.
  row_complete <- stats::complete.cases(
    matched_data[, required_cols, drop = FALSE]
  )
  incomplete_sets <- unique(matched_data[[set_id_var]][!row_complete])
  analysis_data <- matched_data[
    !(matched_data[[set_id_var]] %in% incomplete_sets), , drop = FALSE
  ]

  if (nrow(analysis_data) == 0L) {
    stop("No complete matched sets are available for outcome analysis.")
  }

  expected_size <- attr(matched_data, "K")
  if (is.null(expected_size)) {
    expected_size <- length(unique(as.character(matched_data[[treatment_var]])))
  }
  retained_sizes <- table(analysis_data[[set_id_var]])
  if (any(retained_sizes != expected_size)) {
    stop("Outcome analysis contains incomplete matched sets: ",
         sum(retained_sizes != expected_size), " of ", length(retained_sizes),
         " do not have ", expected_size, " subjects.", call. = FALSE)
  }

  analysis_data[[outcome_var]] <- as.integer(analysis_data[[outcome_var]])

  if (is.null(anchor_level)) {
    anchor_level <- attr(matched_data, "anchor_level")
  }
  if (is.null(anchor_level) || length(anchor_level) != 1L) {
    stop("`anchor_level` must be supplied or available as an attribute ",
         "from `extract_matched_data()`.")
  }
  anchor_level <- treatment_level(anchor_level)

  treatment_levels <- unique(
    as.character(analysis_data[[treatment_var]])
  )
  if (!(anchor_level %in% treatment_levels)) {
    stop("`anchor_level` is not present in the matched dataset.")
  }

  # Preserve the comparator order recorded by extract_matched_data()
  # whenever possible; otherwise retain observed order.
  stored_groups <- attr(matched_data, "groups")
  if (!is.null(stored_groups)) {
    stored_groups <- as.character(stored_groups)
    comparator_levels <- stored_groups[stored_groups %in% treatment_levels]
  } else {
    comparator_levels <- setdiff(treatment_levels, anchor_level)
  }

  analysis_data[[treatment_var]] <- factor(
    analysis_data[[treatment_var]],
    levels = c(anchor_level, comparator_levels)
  )

  # A treatment group with no events, or no non-events, separates the
  # treatment-only model completely. glm() reports converged = TRUE for such a
  # fit -- it stops on its deviance tolerance having pushed the coefficient
  # toward infinity -- so the convergence check below never fires. The affected
  # rows are flagged as well as warned about, so callers filtering results
  # programmatically have something to test.
  separated_levels <- character(0)

  for (level in c(anchor_level, comparator_levels)) {
    outcomes <- analysis_data[[outcome_var]][
      as.character(analysis_data[[treatment_var]]) == level
    ]
    if (length(outcomes) > 0L &&
        (sum(outcomes) == 0L || sum(outcomes) == length(outcomes))) {
      separated_levels <- c(separated_levels, level)
      warning("Treatment group '", level, "' has ", sum(outcomes), "/",
              length(outcomes), " events. The treatment-only logistic model ",
              "may exhibit complete separation and OR inference may be ",
              "unstable.", call. = FALSE)
    }
  }

  # Fit the unadjusted treatment-outcome model
  effect_formula <- stats::reformulate(
    treatment_var,
    response = outcome_var
  )
  fit <- stats::glm(
    effect_formula,
    data = analysis_data,
    family = stats::binomial()
  )

  if (!isTRUE(fit$converged) || isTRUE(fit$boundary)) {
    warning("The logistic outcome model did not converge cleanly or ",
            "reached the parameter boundary.")
  }

  vcov_cluster <- .sam_cluster_vcov(
    fit,
    cluster = analysis_data[[set_id_var]]
  )
  beta <- stats::coef(fit)
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Marginal risk by treatment group
  all_levels <- c(anchor_level, comparator_levels)
  group_risk_rows <- lapply(all_levels, function(g) {
    x_g <- stats::setNames(rep(0, length(beta)), names(beta))
    x_g["(Intercept)"] <- 1

    if (g != anchor_level) {
      coef_name <- paste0(treatment_var, g)
      if (!(coef_name %in% names(beta))) {
        stop("Could not find treatment coefficient: ", coef_name)
      }
      x_g[coef_name] <- 1
    }

    eta_g <- sum(x_g * beta)
    se_eta_g <- sqrt(as.numeric(t(x_g) %*% vcov_cluster %*% x_g))
    risk_g <- stats::plogis(eta_g)

    group_data <- analysis_data[
      analysis_data[[treatment_var]] == g, , drop = FALSE
    ]

    data.frame(
      treatment = g,
      n = nrow(group_data),
      events = sum(group_data[[outcome_var]]),
      risk = risk_g,
      risk_ci_low = stats::plogis(eta_g - z_value * se_eta_g),
      risk_ci_high = stats::plogis(eta_g + z_value * se_eta_g),
      stringsAsFactors = FALSE
    )
  })
  group_risk <- do.call(rbind, group_risk_rows)
  rownames(group_risk) <- NULL

  # Anchor-referenced OR, RR, and RD for every comparator
  contrast_rows <- lapply(comparator_levels, function(g) {
    coef_name <- paste0(treatment_var, g)

    get_risks <- function(b) {
      p_A <- stats::plogis(b["(Intercept)"])
      p_g <- stats::plogis(b["(Intercept)"] + b[coef_name])
      c(p_A = unname(p_A), p_g = unname(p_g))
    }

    log_or_fun <- function(b) unname(b[coef_name])
    log_rr_fun <- function(b) {
      p <- get_risks(b)
      log(p["p_g"] / p["p_A"])
    }
    rd_fun <- function(b) {
      p <- get_risks(b)
      unname(p["p_g"] - p["p_A"])
    }

    # Odds ratio
    log_or <- log_or_fun(beta)
    grad_log_or <- .sam_numeric_gradient(log_or_fun, beta)
    se_log_or <- sqrt(as.numeric(
      t(grad_log_or) %*% vcov_cluster %*% grad_log_or
    ))
    OR <- exp(log_or)
    OR_ci_low <- exp(log_or - z_value * se_log_or)
    OR_ci_high <- exp(log_or + z_value * se_log_or)

    # Risk ratio
    log_rr <- log_rr_fun(beta)
    grad_log_rr <- .sam_numeric_gradient(log_rr_fun, beta)
    se_log_rr <- sqrt(as.numeric(
      t(grad_log_rr) %*% vcov_cluster %*% grad_log_rr
    ))
    RR <- exp(log_rr)
    RR_ci_low <- exp(log_rr - z_value * se_log_rr)
    RR_ci_high <- exp(log_rr + z_value * se_log_rr)

    # Risk difference
    RD <- rd_fun(beta)
    grad_rd <- .sam_numeric_gradient(rd_fun, beta)
    se_RD <- sqrt(as.numeric(
      t(grad_rd) %*% vcov_cluster %*% grad_rd
    ))
    RD_ci_low <- RD - z_value * se_RD
    RD_ci_high <- RD + z_value * se_RD

    risks <- get_risks(beta)

    anchor_idx <- analysis_data[[treatment_var]] == anchor_level
    comparator_idx <- analysis_data[[treatment_var]] == g

    data.frame(
      anchor = anchor_level,
      comparator = g,
      n_anchor = sum(anchor_idx),
      n_comparator = sum(comparator_idx),
      events_anchor = sum(analysis_data[[outcome_var]][anchor_idx]),
      events_comparator = sum(analysis_data[[outcome_var]][comparator_idx]),
      risk_anchor = unname(risks["p_A"]),
      risk_comparator = unname(risks["p_g"]),
      log_or = log_or,
      se_log_or = se_log_or,
      OR = OR,
      OR_ci_low = OR_ci_low,
      OR_ci_high = OR_ci_high,
      log_rr = log_rr,
      se_log_rr = se_log_rr,
      RR = RR,
      RR_ci_low = RR_ci_low,
      RR_ci_high = RR_ci_high,
      RD = RD,
      se_RD = se_RD,
      RD_ci_low = RD_ci_low,
      RD_ci_high = RD_ci_high,
      separation = anchor_level %in% separated_levels ||
        g %in% separated_levels,
      stringsAsFactors = FALSE
    )
  })

  contrasts <- do.call(rbind, contrast_rows)
  rownames(contrasts) <- NULL

  analysis_summary <- data.frame(
    anchor_level = anchor_level,
    n_matched_sets = length(unique(analysis_data[[set_id_var]])),
    K = length(levels(analysis_data[[treatment_var]])),
    n_subjects = nrow(analysis_data),
    confidence_level = conf_level,
    stringsAsFactors = FALSE
  )

  list(
    analysis_summary = analysis_summary,
    group_risk = group_risk,
    contrasts = contrasts,
    model = fit,
    vcov_cluster = vcov_cluster
  )
}

