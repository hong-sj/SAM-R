################################################################################
### Multinomial Logistic Regression: Estimate the Generalized Propensity Score
## 
## Author: Sungjun Hong
################################################################################

#' Estimate generalized propensity scores
#'
#' @description
#' Estimates generalized propensity scores using multinomial logistic
#' regression with `anchor_level` as the reference treatment group.
#'
#' @param data A `data.frame` containing treatment and covariate variables.
#' @param X_vars Character vector of covariate column names.
#'   Default `paste0("X", 1:10)`.
#' @param treatment_var Name of the treatment variable. Default `"T"`.
#' @param anchor_level Reference treatment group. Default `"A"`.
#'
#' @return A list containing the fitted multinomial model (`model`) and
#'   predicted treatment probabilities (`gps`).
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
#' head(fit$gps)
#' rowSums(fit$gps)[1:5]
#'
#' @export
estimate_gps_multinom <- function(data, X_vars = paste0("X", 1:10),
                                   treatment_var = "T", anchor_level = "A") {
  labels <- treatment_labels(data, treatment_var)
  anchor_level <- treatment_level(anchor_level)
  require_rows(which(labels == anchor_level), anchor_level, treatment_var,
               available = labels)

  # Anchor as the reference category. The level order of an existing factor is
  # preserved, but `ref` must be passed as a character label: relevel() reads a
  # numeric `ref` as a position in levels(), so a numeric anchor_level would
  # otherwise select a different reference category than the caller named.
  treatment_factor <- stats::relevel(factor(data[[treatment_var]]),
                                     ref = anchor_level)

  # Standardize covariates for numerical stability
  X_raw <- covariate_matrix(data, X_vars)
  X_mean <- colMeans(X_raw)
  X_sd <- apply(X_raw, 2, stats::sd)
  X_sd[X_sd == 0] <- 1
  X_std <- scale(X_raw, center = X_mean, scale = X_sd)

  # data.frame() applies make.names() unless told otherwise, so a covariate
  # named "age yrs" would land in the model frame as "age.yrs" while the
  # formula below still asked for "age yrs". Keeping the names verbatim means
  # the formula has to quote them.
  if ("treatment_factor" %in% X_vars) {
    stop("A covariate may not be named \"treatment_factor\": it collides with ",
         "the internal response column.", call. = FALSE)
  }

  model_data <- data.frame(treatment_factor = treatment_factor, X_std,
                           check.names = FALSE)
  formula <- stats::as.formula(paste(
    "treatment_factor ~",
    paste(sprintf("`%s`", X_vars), collapse = " + ")
  ))

  model <- nnet::multinom(formula, data = model_data, trace = FALSE,
                           maxit = 5000, reltol = 1e-10)

  # Predict on the fitted scale, before the coefficients below are rewritten,
  # so the GPS is exactly what the fit produced.
  gps <- stats::predict(model, newdata = model_data, type = "probs")

  model <- unstandardize_multinom(model, X_vars, X_mean, X_sd)

  # Ensure a consistent GPS matrix for two-group treatments
  if (is.null(dim(gps))) {
    other_level <- setdiff(levels(treatment_factor), anchor_level)
    gps <- cbind(1 - gps, gps)
    colnames(gps) <- c(anchor_level, other_level)
  }

  # Place the anchor group first
  gps <- gps[, c(anchor_level, setdiff(colnames(gps), anchor_level)), drop = FALSE]

  # The rewrite is exact algebra, so a disagreement here means the weight
  # layout assumed above no longer holds -- which must not pass silently, since
  # the whole point is that `model` can be used to score new subjects.
  check <- stats::predict(model, newdata = data[, X_vars, drop = FALSE],
                          type = "probs")
  if (is.null(dim(check))) {
    check <- cbind(1 - check, check)
    colnames(check) <- c(anchor_level, setdiff(levels(treatment_factor),
                                               anchor_level))
  }
  drift <- max(abs(check[, colnames(gps), drop = FALSE] - gps))
  if (!is.finite(drift) || drift > 1e-8) {
    stop("Internal error: rescaling the fitted coefficients changed the ",
         "predictions by ", format(drift, digits = 3), ". Please report this ",
         "with your nnet version (", utils::packageVersion("nnet"), ").",
         call. = FALSE)
  }

  list(model = model, gps = gps)
}


#' Rewrite a multinom fit so it applies to unstandardized covariates
#'
#' The covariates are standardized to condition the fit, which leaves the
#' returned model expecting standardized input. Nothing about the object says
#' so, and the obvious next step -- scoring subjects with it -- is then wrong:
#' on the bundled four-group data, predicting from the covariates as they
#' appear in `data` differed from the returned GPS by up to 0.91.
#'
#' For a linear predictor `b0 + sum_j b_j (x_j - m_j) / s_j`, folding the
#' shift and scale back in gives `b_j / s_j` and an intercept of
#' `b0 - sum_j b_j m_j / s_j`.
#'
#' `nnet::multinom` stores its weights as one block per outcome level, each
#' holding a bias, the model-matrix intercept and then one entry per covariate
#' in column order. The caller verifies the result rather than trusting that
#' layout.
#'
#' @param model A fitted `nnet::multinom` object.
#' @param X_vars Covariate names, in design-matrix column order.
#' @param X_mean,X_sd The centring and scaling that were applied.
#'
#' @return `model`, with its weights on the original covariate scale.
#'
#' @keywords internal
#' @noRd
unstandardize_multinom <- function(model, X_vars, X_mean, X_sd) {
  p <- length(X_vars)
  block <- p + 2L                      # bias, intercept, one per covariate
  weights <- model$wts

  if (length(weights) %% block != 0L) {
    stop("Unexpected nnet::multinom weight layout: ", length(weights),
         " weights for ", p, " covariates.", call. = FALSE)
  }

  for (start in seq(0L, length(weights) - block, by = block)) {
    intercept <- start + 2L
    slopes <- start + 2L + seq_len(p)

    weights[slopes] <- weights[slopes] / X_sd
    weights[intercept] <- weights[intercept] -
      sum(weights[slopes] * X_mean)
  }

  model$wts <- weights
  model
}
