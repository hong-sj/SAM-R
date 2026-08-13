################################################################################
### Tests for SAM Matching
################################################################################

data("sample_4group", package = "SAMatch")

covariates <- setdiff(
  names(sample_4group),
  c("synthetic_id", "treatment", "mortality_28d")
)
treatment_levels <- unique(as.character(sample_4group$treatment))
anchor <- treatment_levels[1]


# Generalized propensity scores ------------------------------------------------

test_that("estimate_gps_multinom returns a valid GPS matrix", {
  fit <- estimate_gps_multinom(
    sample_4group,
    X_vars = covariates,
    treatment_var = "treatment",
    anchor_level = anchor
  )
  
  expect_equal(nrow(fit$gps), nrow(sample_4group))
  expect_setequal(colnames(fit$gps), treatment_levels)
  expect_equal(colnames(fit$gps)[1], anchor)
  expect_true(all(fit$gps >= 0 & fit$gps <= 1))
  expect_equal(
    unname(rowSums(fit$gps)),
    rep(1, nrow(sample_4group)),
    tolerance = 1e-6
  )
})


# GPS candidate search ---------------------------------------------------------

test_that("gps_candidate_search returns valid candidate pools", {
  fit <- estimate_gps_multinom(
    sample_4group,
    X_vars = covariates,
    treatment_var = "treatment",
    anchor_level = anchor
  )
  
  search <- gps_candidate_search(
    sample_4group,
    fit$gps,
    treatment_var = "treatment",
    anchor_level = anchor,
    top_m = 5,
    gps_space = "logit"
  )
  
  expect_equal(length(search$candidates), length(search$anchor_rows))
  expect_setequal(search$groups, setdiff(treatment_levels, anchor))
  
  candidate_lengths <- vapply(
    search$candidates,
    function(x) all(vapply(x, length, integer(1)) <= 5),
    logical(1)
  )
  expect_true(all(candidate_lengths))
  
  # Every candidate must belong to the corresponding comparator group
  for (g in search$groups) {
    candidate_rows <- unique(unlist(lapply(search$candidates, `[[`, g)))
    expect_true(
      all(as.character(sample_4group$treatment[candidate_rows]) == g)
    )
  }
})


# SAM matching -----------------------------------------------------------------

test_that("sam_match returns consistent matched sets without duplicate use", {
  fit <- estimate_gps_multinom(
    sample_4group,
    X_vars = covariates,
    treatment_var = "treatment",
    anchor_level = anchor
  )
  
  search <- gps_candidate_search(
    sample_4group,
    fit$gps,
    treatment_var = "treatment",
    anchor_level = anchor,
    top_m = 10,
    gps_space = "logit"
  )
  
  result <- sam_match(
    sample_4group,
    search,
    X_vars = covariates,
    treatment_var = "treatment"
  )
  
  expect_true(result$matching_rate >= 0 && result$matching_rate <= 1)
  expect_equal(
    result$matching_rate,
    nrow(result$matched) / length(search$anchor_rows)
  )
  
  if (nrow(result$matched) > 0L) {
    dist_cols <- paste0("dist_", search$groups)
    
    # Total loss must equal the sum of group-specific distances
    expect_equal(
      result$matched$loss,
      rowSums(as.matrix(result$matched[, dist_cols, drop = FALSE])),
      tolerance = 1e-8
    )
    
    # Each subject can be used at most once
    expect_false(anyDuplicated(result$matched$anchor) > 0)
    for (g in search$groups) {
      expect_false(anyDuplicated(result$matched[[g]]) > 0)
    }
    
    # Matched subjects must belong to the corresponding treatment groups
    expect_true(
      all(as.character(sample_4group$treatment[result$matched$anchor]) == anchor)
    )
    for (g in search$groups) {
      expect_true(
        all(as.character(sample_4group$treatment[result$matched[[g]]]) == g)
      )
    }
  }
  
  # Matched and unmatched anchors must partition all anchor subjects
  matched_anchor_rows <- if (nrow(result$matched) > 0L) {
    result$matched$anchor
  } else {
    integer(0)
  }
  
  expect_setequal(
    c(matched_anchor_rows, result$unmatched_anchor_rows),
    search$anchor_rows
  )
  expect_length(
    intersect(matched_anchor_rows, result$unmatched_anchor_rows),
    0
  )
})


# SAM evaluation ---------------------------------------------------------------

test_that("sam_evaluate returns valid matching diagnostics", {
  fit <- estimate_gps_multinom(
    sample_4group,
    X_vars = covariates,
    treatment_var = "treatment",
    anchor_level = anchor
  )
  
  search <- gps_candidate_search(
    sample_4group,
    fit$gps,
    treatment_var = "treatment",
    anchor_level = anchor,
    top_m = 10,
    gps_space = "logit"
  )
  
  result <- sam_match(
    sample_4group,
    search,
    X_vars = covariates,
    treatment_var = "treatment"
  )
  
  report <- sam_evaluate(
    sample_4group,
    search,
    result,
    fit$gps,
    X_vars = covariates,
    treatment_var = "treatment"
  )
  
  expected_names <- c(
    "loss_distribution",
    "dispersion_distribution",
    "smd_balance",
    "treatment_discrimination_auc",
    "matching_rate"
  )
  expect_true(all(expected_names %in% names(report)))
  
  # AUC values must lie between 0 and 1
  auc_vals <- report$treatment_discrimination_auc$pairwise$auc
  expect_true(all(auc_vals >= 0 & auc_vals <= 1, na.rm = TRUE))
  
  mean_auc <- report$treatment_discrimination_auc$mean_auc
  if (!is.na(mean_auc)) {
    expect_true(mean_auc >= 0 && mean_auc <= 1)
  }
  
  expect_equal(nrow(report$loss_distribution), 1)
  expect_equal(nrow(report$dispersion_distribution), 1)
  expect_true(
    all(c("mean", "median", "sd", "p95", "max") %in%
          names(report$loss_distribution))
  )
})


# Factor treatment handling ----------------------------------------------------

test_that("SAM handles a factor treatment variable", {
  data_factor <- sample_4group
  data_factor$treatment <- factor(data_factor$treatment)
  
  fit <- estimate_gps_multinom(
    data_factor,
    X_vars = covariates,
    treatment_var = "treatment",
    anchor_level = anchor
  )
  
  search <- gps_candidate_search(
    data_factor,
    fit$gps,
    treatment_var = "treatment",
    anchor_level = anchor,
    top_m = 10,
    gps_space = "logit"
  )
  
  result <- sam_match(
    data_factor,
    search,
    X_vars = covariates,
    treatment_var = "treatment"
  )
  
  report <- sam_evaluate(
    data_factor,
    search,
    result,
    fit$gps,
    X_vars = covariates,
    treatment_var = "treatment"
  )
  
  expect_true(is.list(report$treatment_discrimination_auc))
})


# Mahalanobis distance ----------------------------------------------------------

test_that("mahalanobis_distance_matrix returns clean numeric output", {
  pooled <- get_pooled_covariance(
    sample_4group,
    X_vars = covariates,
    treatment_var = "treatment"
  )
  
  anchor_rows <- which(as.character(sample_4group$treatment) == anchor)
  comparator <- setdiff(treatment_levels, anchor)[1]
  comparator_rows <- which(
    as.character(sample_4group$treatment) == comparator
  )
  
  X_anchor <- as.matrix(
    sample_4group[anchor_rows, covariates, drop = FALSE]
  )
  X_comparator <- as.matrix(
    sample_4group[comparator_rows, covariates, drop = FALSE]
  )
  
  D <- mahalanobis_distance_matrix(
    X_anchor,
    X_comparator,
    pooled$S_inv
  )
  
  expect_true(is.matrix(D))
  expect_true(is.numeric(D))
  expect_equal(dim(D), c(length(anchor_rows), length(comparator_rows)))
  expect_null(dimnames(D))
  
  scalar_value <- D[1, 1] + 0
  expect_null(names(scalar_value))
})