################################################################################
### Mathematical Utilities
##
## Author: Sungjun Hong
################################################################################


#' Logistic function
#'
#' @description
#' Computes the inverse-logit transformation.
#'
#' @param x Numeric vector on the log-odds scale.
#'
#' @return Numeric vector of probabilities between 0 and 1.
#'
#' @examples
#' expit(0)
#' expit(c(-2, 0, 2))
#'
#' @export
expit <- function(x) {
  stats::plogis(x)
}

#' Logit function
#'
#' @description
#' Computes the logit transformation of probabilities.
#'
#' @param p Numeric vector of probabilities between 0 and 1.
#'
#' @return Numeric vector on the log-odds scale.
#'
#' @examples
#' logit(0.5)
#' logit(c(0.1, 0.5, 0.9))
#'
#' @export
logit <- function(p) {
  stats::qlogis(p)
}

#' Area under the ROC curve using the Mann-Whitney statistic
#'
#' @description
#' Computes the area under the ROC curve (AUC) for a continuous score and
#' binary outcome using the rank-based Mann-Whitney U statistic.
#'
#' @param score Numeric vector of continuous scores.
#' @param label Binary vector coded as 0 and 1.
#'
#' @return Numeric AUC value between 0 and 1, or `NA` if one outcome class
#'   is absent.
#'
#' @examples
#' score <- c(0.1, 0.3, 0.6, 0.9)
#' label <- c(0, 0, 1, 1)
#'
#' auc_mannwhitney(score, label)
#'
#' @export
auc_mannwhitney <- function(score, label) {

  # Rank scores using average ranks for ties
  r <- rank(score)
  n1 <- sum(label == 1L)
  n0 <- sum(label == 0L)
  if (n1 == 0L || n0 == 0L) {
    return(NA_real_)
  }

  # Convert the Mann-Whitney U statistic to AUC
  sum_ranks_pos <- sum(r[label == 1L])
  u <- sum_ranks_pos - n1 * (n1 + 1) / 2
  u / (n1 * n0)
}
