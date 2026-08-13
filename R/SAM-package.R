#' SAM (Shared Anchor Matching): A Scalable Matching Framework for Multiple Treatment Groups
#'
#' @description
#' SAM (Shared Anchor Matching) provides a scalable framework for
#' simultaneously matching multiple treatment groups to a shared anchor group.
#'
#' SAM combines generalized propensity score estimation, GPS-guided
#' candidate screening, and Mahalanobis-distance-based matching. The package
#' also provides covariate-balance diagnostics, matched-cohort outcome
#' analysis, three-way propensity score matching, and multi-arm weighting
#' methods including IPTW, overlap weighting, and matching weights.
#'
#' The primary SAM workflow is:
#'
#' [estimate_gps_multinom()] ->
#' [gps_candidate_search()] ->
#' [sam_match()] ->
#' [sam_evaluate()]
#'
#' Matched subject-level data can be reconstructed with
#' [extract_matched_data()] and analyzed using [sam_estimate_effects()].
#'
#' @keywords internal
"_PACKAGE"