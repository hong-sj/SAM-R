################################################################################
### Matched-cohort extraction and outcome analysis (#12)
################################################################################

# A completely separated fit drives one fitted risk to the boundary, so the
# delta-method gradient for log(RR) legitimately evaluates log(0) and R emits
# "NaNs produced". That is a consequence of the separation the test is about,
# not a second defect, so it is muffled explicitly rather than globally.
muffle_boundary_nans <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("NaN", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

matched_cohort <- function() {
  fixture <- sam_fixture()
  extract_matched_data(fixture$data, fixture$search, fixture$match,
                       treatment_var = "treatment",
                       anchor_level = fixture$anchor)
}

test_that("extract_matched_data expands each set into K subject records", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  K <- length(fixture$groups) + 1L
  expect_equal(nrow(cohort), nrow(fixture$match$matched) * K)
  expect_equal(attr(cohort, "K"), K)
  expect_identical(attr(cohort, "anchor_level"), fixture$anchor)
  expect_identical(attr(cohort, "groups"), fixture$groups)

  # Exactly one subject per treatment group in every set.
  counts <- table(cohort$matched_set_id,
                  as.character(cohort$treatment))
  expect_true(all(counts == 1L))

  # original_row points back at the row it came from.
  expect_equal(cohort$treatment,
               fixture$data$treatment[cohort$original_row])
})

test_that("the cohort is ordered by set, then anchor first", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  expect_equal(cohort$matched_set_id, sort(cohort$matched_set_id))

  order_within_set <- match(as.character(cohort$treatment),
                            c(fixture$anchor, fixture$groups))
  by_set <- split(order_within_set, cohort$matched_set_id)
  expect_true(all(vapply(by_set, function(x) identical(x, sort(x)),
                         logical(1))))
})

test_that("extract_matched_data refuses inconsistent input", {
  fixture <- sam_fixture()

  empty <- fixture$match
  empty$matched <- empty_matched_frame(fixture$groups)
  expect_error(
    extract_matched_data(fixture$data, fixture$search, empty,
                         treatment_var = "treatment"),
    "no matched sets"
  )

  truncated <- fixture$match
  truncated$matched <- truncated$matched[, setdiff(names(truncated$matched),
                                                   fixture$groups[1])]
  expect_error(
    extract_matched_data(fixture$data, fixture$search, truncated,
                         treatment_var = "treatment"),
    "Missing matched-set column"
  )

  out_of_range <- fixture$match
  out_of_range$matched[[fixture$groups[1]]][1] <- nrow(fixture$data) + 1L
  expect_error(
    extract_matched_data(fixture$data, fixture$search, out_of_range,
                         treatment_var = "treatment"),
    "Invalid row index"
  )
})

test_that("the unadjusted model reproduces the observed risks exactly", {
  # An unadjusted logistic regression on a factor treatment is saturated, so
  # its fitted risk for each group must equal that group's event proportion.
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)

  for (i in seq_len(nrow(effects$group_risk))) {
    group <- effects$group_risk$treatment[i]
    outcomes <- cohort$mortality_28d[
      as.character(cohort$treatment) == group
    ]
    expect_equal(effects$group_risk$n[i], length(outcomes))
    expect_equal(effects$group_risk$events[i], sum(outcomes))
    expect_equal(effects$group_risk$risk[i], mean(outcomes),
                 tolerance = 1e-8, info = group)
  }
})

test_that("OR, RR and RD agree with the two-by-two table", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)

  anchor_outcomes <- cohort$mortality_28d[
    as.character(cohort$treatment) == fixture$anchor
  ]
  risk_anchor <- mean(anchor_outcomes)

  for (i in seq_len(nrow(effects$contrasts))) {
    comparator <- effects$contrasts$comparator[i]
    comparator_outcomes <- cohort$mortality_28d[
      as.character(cohort$treatment) == comparator
    ]
    risk_comparator <- mean(comparator_outcomes)

    odds <- function(p) p / (1 - p)
    expect_equal(effects$contrasts$OR[i],
                 odds(risk_comparator) / odds(risk_anchor),
                 tolerance = 1e-6, info = comparator)
    expect_equal(effects$contrasts$RR[i], risk_comparator / risk_anchor,
                 tolerance = 1e-6, info = comparator)
    expect_equal(effects$contrasts$RD[i], risk_comparator - risk_anchor,
                 tolerance = 1e-6, info = comparator)
  }
})

test_that("confidence intervals bracket the point estimates", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)
  contrasts <- effects$contrasts

  expect_true(all(contrasts$OR_ci_low < contrasts$OR))
  expect_true(all(contrasts$OR < contrasts$OR_ci_high))
  expect_true(all(contrasts$RR_ci_low < contrasts$RR))
  expect_true(all(contrasts$RR < contrasts$RR_ci_high))
  expect_true(all(contrasts$RD_ci_low < contrasts$RD))
  expect_true(all(contrasts$RD < contrasts$RD_ci_high))

  # A wider confidence level gives a wider interval.
  wide <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                               treatment_var = "treatment",
                               anchor_level = fixture$anchor,
                               conf_level = 0.99)$contrasts
  expect_true(all(wide$OR_ci_low < contrasts$OR_ci_low))
  expect_true(all(wide$OR_ci_high > contrasts$OR_ci_high))
})

test_that("the cluster-robust covariance is symmetric and positive definite", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)
  vcov_cluster <- effects$vcov_cluster

  expect_equal(vcov_cluster, t(vcov_cluster))
  expect_true(all(eigen(vcov_cluster, only.values = TRUE)$values > 0))
  expect_identical(dim(vcov_cluster),
                   rep(length(stats::coef(effects$model)), 2L))
})

test_that("sam_estimate_effects validates its input", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  expect_error(sam_estimate_effects(as.matrix(1:4), outcome_var = "y"),
               "must be a data.frame")
  expect_error(
    sam_estimate_effects(cohort, outcome_var = "nope",
                         treatment_var = "treatment"),
    "Missing required column"
  )
  expect_error(
    sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                         treatment_var = "treatment", conf_level = 1),
    "strictly between 0 and 1"
  )

  non_binary <- cohort
  non_binary$mortality_28d[1] <- 7
  expect_error(
    sam_estimate_effects(non_binary, outcome_var = "mortality_28d",
                         treatment_var = "treatment"),
    "binary and coded as 0/1"
  )

  no_anchor <- cohort
  attr(no_anchor, "anchor_level") <- NULL
  expect_error(
    sam_estimate_effects(no_anchor, outcome_var = "mortality_28d",
                         treatment_var = "treatment"),
    "must be supplied or available as an attribute"
  )
})

test_that("the analysis summary describes the cohort it analysed", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)
  summary_row <- effects$analysis_summary

  expect_equal(summary_row$n_subjects, nrow(cohort))
  expect_equal(summary_row$n_matched_sets, nrow(fixture$match$matched))
  expect_equal(summary_row$K, length(fixture$groups) + 1L)
  expect_identical(summary_row$anchor_level, fixture$anchor)
  expect_equal(summary_row$confidence_level, 0.95)
})

test_that("a missing value drops the whole matched set, not just the row", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()
  K <- length(fixture$groups) + 1L

  complete <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                   treatment_var = "treatment",
                                   anchor_level = fixture$anchor)

  with_gap <- cohort
  with_gap$mortality_28d[1] <- NA

  # Previously: the row was dropped and its set stayed in the analysis with
  # K - 1 members, silently breaking the structure the cluster-robust variance
  # assumes.
  reduced <- sam_estimate_effects(with_gap, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)

  expect_equal(reduced$analysis_summary$n_matched_sets,
               complete$analysis_summary$n_matched_sets - 1L)
  expect_equal(reduced$analysis_summary$n_subjects,
               complete$analysis_summary$n_subjects - K)
})

test_that("every retained matched set still has K subjects", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()
  K <- length(fixture$groups) + 1L

  # A cohort whose sets are already broken must be refused rather than
  # analysed, which is the check extract_matched_data() already performs.
  broken <- cohort[-1, , drop = FALSE]
  attr(broken, "K") <- K

  expect_error(
    sam_estimate_effects(broken, outcome_var = "mortality_28d",
                         treatment_var = "treatment",
                         anchor_level = fixture$anchor),
    "incomplete matched sets"
  )
})

test_that("no complete matched set left is an error, not an empty model", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()
  cohort$mortality_28d <- NA

  expect_error(
    sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                         treatment_var = "treatment",
                         anchor_level = fixture$anchor),
    "No complete matched sets"
  )
})

test_that("complete outcome separation is warned about and flagged", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()
  separated_group <- fixture$groups[1]
  cohort$mortality_28d[
    as.character(cohort$treatment) == separated_group
  ] <- 0

  # Previously: glm() reports converged = TRUE for a separated fit, so the
  # existing convergence check never fired and an OR of 1.4e-08 was returned
  # as though it were an estimate.
  muffle_boundary_nans(
    expect_warning(
      effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                      treatment_var = "treatment",
                                      anchor_level = fixture$anchor),
      "complete separation"
    )
  )

  contrasts <- effects$contrasts
  expect_true(contrasts$separation[contrasts$comparator == separated_group])
  expect_false(any(contrasts$separation[
    contrasts$comparator != separated_group
  ]))
})

test_that("a separated anchor flags every contrast", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()
  cohort$mortality_28d[
    as.character(cohort$treatment) == fixture$anchor
  ] <- 1

  muffle_boundary_nans(
    expect_warning(
      effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                      treatment_var = "treatment",
                                      anchor_level = fixture$anchor),
      "complete separation"
    )
  )
  expect_true(all(effects$contrasts$separation))
})

test_that("a healthy cohort reports no separation", {
  fixture <- sam_fixture()
  cohort <- matched_cohort()

  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = fixture$anchor)
  expect_false(any(effects$contrasts$separation))
})
