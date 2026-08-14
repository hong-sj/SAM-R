################################################################################
### Golden fixtures for the bundled pipelines (#12)
##
## The fixtures pin the output of two iterative solvers -- BFGS for the GPS,
## IRLS for the outcome model -- which stop on a convergence tolerance rather
## than at machine precision. Their last digits therefore differ legitimately
## across versions of nnet, MASS and R itself.
##
## The comparison is tiered accordingly:
##
##   estimates   full precision on the environment recorded in ENVIRONMENT.dcf,
##               loosely elsewhere -- with an absolute floor as well as a
##               relative tolerance, because a standardized mean difference
##               sits near zero and a purely relative tolerance is meaningless
##               there
##   decisions   matched sets are integer row indices; there is no meaningful
##               loose comparison of one, so they are compared exactly on the
##               reference environment and skipped elsewhere
##   shape       checked everywhere
################################################################################

golden_path <- function(name) {
  testthat::test_path("golden", paste0(name, ".rds"))
}

golden <- function(name) {
  path <- golden_path(name)
  if (!file.exists(path)) {
    testthat::skip(paste("golden fixture not found:", name))
  }
  readRDS(path)
}

recorded_environment <- function() {
  path <- testthat::test_path("golden", "ENVIRONMENT.dcf")
  if (!file.exists(path)) {
    return(NULL)
  }
  as.data.frame(read.dcf(path), stringsAsFactors = FALSE)
}

current_environment <- function() {
  data.frame(
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    nnet = as.character(utils::packageVersion("nnet")),
    MASS = as.character(utils::packageVersion("MASS")),
    RANN = as.character(utils::packageVersion("RANN")),
    stringsAsFactors = FALSE
  )
}

#' Differences between the recorded and current stacks, as a readable string
environment_drift <- function() {
  recorded <- recorded_environment()
  if (is.null(recorded)) {
    return("no ENVIRONMENT.dcf recorded")
  }

  current <- current_environment()
  drift <- vapply(names(recorded), function(component) {
    if (identical(recorded[[component]], current[[component]])) {
      NA_character_
    } else {
      paste0(component, " ", recorded[[component]], " != ", current[[component]])
    }
  }, character(1))
  drift <- drift[!is.na(drift)]

  if (length(drift) == 0L) "" else paste(drift, collapse = ", ")
}

is_reference_environment <- function() {
  identical(environment_drift(), "")
}

skip_unless_reference <- function() {
  if (!is_reference_environment()) {
    testthat::skip(paste0(
      "not the environment the fixtures were generated on: ",
      environment_drift()
    ))
  }
}

#' Compare an estimate against its fixture at the appropriate tier
expect_golden_estimate <- function(current, name) {
  expected <- golden(name)

  if (is_reference_environment()) {
    expect_equal(current, expected, tolerance = 1e-12, info = name)
    return(invisible(NULL))
  }

  # Loose tier: catch a gross regression -- a sign flip, a shifted formula --
  # not solver drift. The absolute floor matters: an SMD of -0.0165 moving by
  # 5e-05 is 0.3% relatively but 2,000 times below the 0.1 threshold anyone
  # reads an SMD against.
  expect_equal(current, expected, tolerance = 1e-5, info = name)

  numeric_current <- unlist(Filter(is.numeric, as.list(current)))
  numeric_expected <- unlist(Filter(is.numeric, as.list(expected)))
  if (length(numeric_current) > 0L) {
    expect_lt(max(abs(numeric_current - numeric_expected), na.rm = TRUE), 1e-3)
  }
}

#' Compare a matching decision against its fixture
expect_golden_decision <- function(current, name) {
  skip_unless_reference()
  expect_equal(current, golden(name), tolerance = 0, info = name)
}


test_that("the recorded environment is legible", {
  recorded <- recorded_environment()
  skip_if(is.null(recorded), "no ENVIRONMENT.dcf recorded")

  expect_setequal(names(recorded), names(current_environment()))
  expect_true(all(nzchar(unlist(recorded))))
})

for (dataset in c("four", "three")) {
  local({
    tag <- dataset
    name <- if (tag == "four") "sample_4group" else "sample_3group"

    test_that(paste(tag, "- GPS matches the fixture"), {
      fixture <- sam_fixture(name)
      expect_golden_estimate(fixture$gps, paste0(tag, "_gps"))
    })

    test_that(paste(tag, "- candidate pools match the fixture"), {
      fixture <- sam_fixture(name)
      # Candidate pools are row indices: decisions, not estimates.
      expect_golden_decision(fixture$search$candidates,
                             paste0(tag, "_candidates"))
    })

    test_that(paste(tag, "- matched sets match the fixture"), {
      fixture <- sam_fixture(name)
      matched <- fixture$match$matched

      # Shape is checked on every environment.
      expect_identical(names(matched), matched_frame_columns(fixture$groups))
      expect_equal(nrow(matched), nrow(golden(paste0(tag, "_matched"))))

      index_columns <- c("matched_set_id", "anchor", fixture$groups)
      expect_golden_decision(
        matched[, index_columns],
        paste0(tag, "_matched_indices")
      )
    })

    test_that(paste(tag, "- distances and rates match the fixture"), {
      fixture <- sam_fixture(name)
      expected <- golden(paste0(tag, "_matched"))

      distance_columns <- c(paste0("dist_", fixture$groups), "loss")
      expect_golden_estimate(fixture$match$matched[, distance_columns],
                             paste0(tag, "_distances"))

      rates <- c(matching_rate = fixture$match$matching_rate,
                 max_possible_rate = fixture$match$max_possible_rate)
      expect_golden_estimate(rates, paste0(tag, "_rates"))
    })

    test_that(paste(tag, "- balance and AUC match the fixture"), {
      fixture <- sam_fixture(name)
      report <- sam_evaluate(fixture$data, fixture$search, fixture$match,
                             fixture$gps, X_vars = fixture$X_vars,
                             treatment_var = "treatment")

      expect_golden_estimate(report$smd_balance, paste0(tag, "_smd"))
      expect_golden_estimate(report$treatment_discrimination_auc,
                             paste0(tag, "_auc"))
      expect_golden_estimate(report$loss_distribution, paste0(tag, "_loss"))
    })

    test_that(paste(tag, "- outcome contrasts match the fixture"), {
      fixture <- sam_fixture(name)
      cohort <- extract_matched_data(fixture$data, fixture$search,
                                     fixture$match,
                                     treatment_var = "treatment",
                                     anchor_level = fixture$anchor)
      effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                      treatment_var = "treatment",
                                      anchor_level = fixture$anchor)

      expect_golden_estimate(effects$contrasts, paste0(tag, "_contrasts"))
      expect_golden_estimate(effects$group_risk, paste0(tag, "_group_risk"))
    })

    test_that(paste(tag, "- weighting output matches the fixture"), {
      fixture <- sam_fixture(name)

      for (method in c("iptw", "overlap", "matching")) {
        weighting <- evaluate_comparator_weighting(
          fixture$data, method = method, X_vars = fixture$X_vars,
          treatment_var = "treatment", anchor_level = fixture$anchor
        )
        expect_golden_estimate(weighting$balance,
                               paste0(tag, "_weighted_balance_", method))
        expect_golden_estimate(weighting$ess, paste0(tag, "_ess_", method))
      }
    })
  })
}

test_that("three-way matching matches the fixture", {
  fixture <- sam_fixture("sample_3group")
  trio <- match_3way(fixture$data, fixture$search, fixture$gps,
                     treatment_var = "treatment", gps_space = "logit")

  expect_identical(names(trio$matched),
                   matched_frame_columns(fixture$groups, "rassen_perimeter"))
  expect_golden_decision(
    trio$matched[, c("matched_set_id", "anchor", fixture$groups)],
    "three_way_indices"
  )
  expect_golden_estimate(
    c(matching_rate = trio$matching_rate,
      max_possible_rate = trio$max_possible_rate,
      caliper = trio$caliper),
    "three_way_meta"
  )
})

test_that("the pipeline is deterministic, with no fixture needed", {
  # This is the one golden-style check that holds on every environment: the
  # same input must give the same output within a single session.
  data <- sam_small_data(n = 120L, n_groups = 3L)
  X_vars <- c("x1", "x2", "x3")

  run_once <- function() {
    fit <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                                 anchor_level = "A")
    search <- gps_candidate_search(data, fit$gps, treatment_var = "T",
                                   anchor_level = "A", top_m = 6,
                                   gps_space = "logit")
    sam_match(data, search, X_vars = X_vars, treatment_var = "T")$matched
  }

  expect_identical(run_once(), run_once())
})
