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
  stopifnot(treatment_var %in% names(data), all(X_vars %in% names(data)))

  # Anchor A as the reference category
  treatment_factor <- stats::relevel(factor(data[[treatment_var]]), ref = anchor_level)

  # Standardize covariates for numerical stability
  X_raw <- as.matrix(data[, X_vars, drop = FALSE])
  X_mean <- colMeans(X_raw)
  X_sd <- apply(X_raw, 2, stats::sd)
  X_sd[X_sd == 0] <- 1
  X_std <- scale(X_raw, center = X_mean, scale = X_sd)
  model_data <- data.frame(treatment_factor = treatment_factor, X_std)
  formula <- stats::as.formula(paste(
    "treatment_factor ~", paste(X_vars, collapse = " + ")
  ))

  model <- nnet::multinom(formula, data = model_data, trace = FALSE,
                           maxit = 5000, reltol = 1e-10)
  gps <- stats::predict(model, newdata = model_data, type = "probs")

  # Ensure a consistent GPS matrix for two-group treatments
  if (is.null(dim(gps))) {
    other_level <- setdiff(levels(treatment_factor), anchor_level)
    gps <- cbind(1 - gps, gps)
    colnames(gps) <- c(anchor_level, other_level)
  }

  # Place the anchor group first
  gps <- gps[, c(anchor_level, setdiff(colnames(gps), anchor_level)), drop = FALSE]

  list(model = model, gps = gps)
}
