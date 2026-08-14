################################################################################
### Regenerate the golden fixtures
##
## Run from the package root:
##
##   Rscript tests/testthat/golden/make-golden.R
##
## Fixtures pin the output of two iterative solvers -- BFGS for the GPS,
## IRLS for the outcome model -- which stop on a convergence tolerance rather
## than at machine precision. The environment that produced them is therefore
## recorded alongside, and test-golden.R compares strictly only there.
################################################################################

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# golden_environment() is defined in R/validate.R so that the generator and
# the comparison in test-golden.R cannot drift apart.

golden_dir <- "tests/testthat/golden"
dir.create(golden_dir, showWarnings = FALSE, recursive = TRUE)

save_fixture <- function(name, value) {
  saveRDS(value, file.path(golden_dir, paste0(name, ".rds")), version = 3)
  invisible(NULL)
}

for (dataset in c("four", "three")) {
  data <- get(utils::data(
    list = if (dataset == "four") "sample_4group" else "sample_3group",
    package = "SAMatch", envir = environment()
  ))
  X_vars <- setdiff(names(data),
                    c("synthetic_id", "treatment", "mortality_28d"))
  anchor <- unique(as.character(data$treatment))[1]

  fit <- estimate_gps_multinom(data, X_vars = X_vars,
                               treatment_var = "treatment",
                               anchor_level = anchor)
  save_fixture(paste0(dataset, "_gps"), fit$gps)

  search <- gps_candidate_search(data, fit$gps, treatment_var = "treatment",
                                 anchor_level = anchor, top_m = 10,
                                 gps_space = "logit")
  save_fixture(paste0(dataset, "_candidates"), search$candidates)

  match_result <- sam_match(data, search, X_vars = X_vars,
                            treatment_var = "treatment")
  matched <- match_result$matched
  groups <- search$groups

  # Split the matched frame along the tiering boundary: row indices are
  # decisions and compared exactly, distances are estimates and compared
  # loosely off the reference environment.
  save_fixture(paste0(dataset, "_matched"), matched)
  save_fixture(paste0(dataset, "_matched_indices"),
               matched[, c("matched_set_id", "anchor", groups)])
  save_fixture(paste0(dataset, "_distances"),
               matched[, c(paste0("dist_", groups), "loss")])
  save_fixture(paste0(dataset, "_rates"), c(
    matching_rate = match_result$matching_rate,
    max_possible_rate = match_result$max_possible_rate
  ))

  report <- sam_evaluate(data, search, match_result, fit$gps,
                         X_vars = X_vars, treatment_var = "treatment")
  save_fixture(paste0(dataset, "_smd"), report$smd_balance)
  save_fixture(paste0(dataset, "_auc"), report$treatment_discrimination_auc)
  save_fixture(paste0(dataset, "_loss"), report$loss_distribution)

  cohort <- extract_matched_data(data, search, match_result,
                                 treatment_var = "treatment",
                                 anchor_level = anchor)
  effects <- sam_estimate_effects(cohort, outcome_var = "mortality_28d",
                                  treatment_var = "treatment",
                                  anchor_level = anchor)
  save_fixture(paste0(dataset, "_contrasts"), effects$contrasts)
  save_fixture(paste0(dataset, "_group_risk"), effects$group_risk)

  for (method in c("iptw", "overlap", "matching")) {
    weighting <- evaluate_comparator_weighting(
      data, method = method, X_vars = X_vars, treatment_var = "treatment",
      anchor_level = anchor
    )
    save_fixture(paste0(dataset, "_weighted_balance_", method),
                 weighting$balance)
    save_fixture(paste0(dataset, "_ess_", method), weighting$ess)
  }

  if (dataset == "three") {
    trio <- match_3way(data, search, fit$gps, treatment_var = "treatment",
                       gps_space = "logit")
    save_fixture("three_way_matched", trio$matched)
    save_fixture("three_way_indices",
                 trio$matched[, c("matched_set_id", "anchor", groups)])
    save_fixture("three_way_meta", c(
      matching_rate = trio$matching_rate,
      max_possible_rate = trio$max_possible_rate,
      caliper = trio$caliper
    ))
  }
}

# The stack that produced the fixtures, so test-golden.R can tell whether it is
# comparing against the environment they were generated on. The platform and
# the linear algebra libraries are part of that: the same package versions on a
# different BLAS move an IRLS standard error in its tenth digit, which is well
# inside what the strict tier compares.
write.dcf(golden_environment(), file.path(golden_dir, "ENVIRONMENT.dcf"))

cat("Wrote", length(list.files(golden_dir, pattern = "[.]rds$")),
    "fixtures and ENVIRONMENT.dcf\n")
