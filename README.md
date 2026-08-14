# SAM-R

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21926956.svg)](https://doi.org/10.5281/zenodo.21926956)

## Shared Anchor Matching: A Scalable Matching Framework for Multiple Treatment Groups

**SAM** is an R package for matching observational data with multiple treatment groups using a shared anchor group.

The framework combines generalized propensity score (GPS) estimation, GPS-guided candidate screening, and Mahalanobis-distance-based matching. Each subject in the anchor group is matched simultaneously to one subject from each comparator group, producing matched sets across multiple treatment groups.

The package also provides tools for covariate-balance assessment, treatment-discrimination diagnostics, matched-cohort outcome analysis, three-way propensity score matching, and multi-arm weighting.

---

## Installation

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hong-sj/SAM-R")
```

Then load the package:

```r
library(SAMatch)
```

---

## Quick Start

The package includes `sample_4group`, an example dataset with four treatment groups.

```r
data("sample_4group")

head(sample_4group)
table(sample_4group$treatment)
```

Define the covariates and anchor treatment:

```r
covariates <- setdiff(
  names(sample_4group),
  c("synthetic_id", "treatment", "mortality_28d")
)

anchor <- unique(as.character(sample_4group$treatment))[1]
```

### 1. Estimate generalized propensity scores

Estimate the GPS using multinomial logistic regression:

```r
fit <- estimate_gps_multinom(
  sample_4group,
  X_vars = covariates,
  treatment_var = "treatment",
  anchor_level = anchor
)

head(fit$gps)
```

`fit$gps` contains the estimated treatment-assignment probabilities for each subject.

### 2. Search for candidate matches

For each anchor subject, identify candidate subjects from each comparator group in GPS space:

```r
search <- gps_candidate_search(
  sample_4group,
  fit$gps,
  treatment_var = "treatment",
  anchor_level = anchor,
  top_m = 10,
  gps_space = "logit"
)
```

SAM retains the 10 nearest candidate subjects from each comparator group for subsequent matching.

### 3. Perform Shared Anchor Matching

Match each anchor subject simultaneously to one subject from every comparator group:

```r
matched <- sam_match(
  sample_4group,
  search,
  X_vars = covariates,
  treatment_var = "treatment"
)

matched$matching_rate
matched$max_possible_rate
head(matched$matched)
```

The resulting object contains the matched sets, group-specific Mahalanobis distances, total matching loss, unmatched anchor subjects, and overall matching rate.

Read `matching_rate` against `max_possible_rate`. Every matched set consumes one subject from each comparator group, so the rate cannot exceed the smallest comparator group divided by the anchor count. On `sample_4group` that ceiling is 59/448 = 0.132, and SAM reaches it exactly — the rate looks like 13% but the matching is saturated, not poor.

### 4. Evaluate matching quality

Evaluate covariate balance and matching diagnostics:

```r
report <- sam_evaluate(
  sample_4group,
  search,
  matched,
  fit$gps,
  X_vars = covariates,
  treatment_var = "treatment"
)

report$matching_rate
report$loss_distribution
report$smd_balance$summary
report$treatment_discrimination_auc
```

The diagnostic output includes matching loss, covariate balance based on standardized mean differences (SMDs), and pairwise treatment-discrimination AUCs.

`sam_evaluate()` bundles the diagnostics, but each can also be called directly. All three infer what they can from the matched sets:

```r
compute_smd_balance(sample_4group, matched$matched, X_vars = covariates)

compute_pairwise_treatment_auc(fit$gps, matched$matched)

get_pooled_covariance(sample_4group, X_vars = covariates,
                      treatment_var = "treatment")
```

`compute_smd_balance()` reports one row per covariate and comparator group, oriented as anchor minus comparator. A covariate with no variance in either arm has no defined SMD; those rows carry `smd_defined = FALSE` and are counted in the summary's `n_undefined` rather than being folded into the mean and maximum.

**Note on row order.** Matched sets are stored as positional row indices, so `data` must be passed to every stage in the row order `gps_candidate_search()` saw. Adding a column, such as an outcome, is fine; re-sorting or filtering between stages is not, and now raises rather than silently reporting balance for a cohort that was never matched.

---

## Matched Cohort and Outcome Analysis

The matched observations can be extracted into a subject-level dataset:

```r
matched_data <- extract_matched_data(
  sample_4group,
  search,
  matched,
  treatment_var = "treatment",
  anchor_level = anchor
)

head(matched_data)
```

Treatment effects can then be estimated in the matched cohort using:

```r
effects <- sam_estimate_effects(
  matched_data,
  outcome_var = "mortality_28d",
  treatment_var = "treatment",
  anchor_level = anchor
)

effects
```

---

## Three-Group Matching

The package also provides `match_3way()` for three-treatment-group propensity score matching.

The included `sample_3group` dataset can be used as an example:

```r
data("sample_3group")

covariates3 <- setdiff(
  names(sample_3group),
  c("synthetic_id", "treatment", "mortality_28d")
)

anchor3 <- unique(as.character(sample_3group$treatment))[1]

fit3 <- estimate_gps_multinom(
  sample_3group,
  X_vars = covariates3,
  treatment_var = "treatment",
  anchor_level = anchor3
)

search3 <- gps_candidate_search(
  sample_3group,
  fit3$gps,
  treatment_var = "treatment",
  anchor_level = anchor3,
  top_m = 10,
  gps_space = "logit"
)

matched3 <- match_3way(
  sample_3group,
  search3,
  fit3$gps,
  treatment_var = "treatment",
  gps_space = "logit"
)

matched3$matching_rate
head(matched3$matched)
```

`match_3way()` matches in two-dimensional propensity score space and takes the same `gps_space` argument as `gps_candidate_search()`. It never uses covariates, so it does not accept `X_vars`; passing one warns rather than discarding it silently.

---

## Multi-Arm Weighting

SAM also includes functions for multi-arm comparator weighting:

- `compute_balancing_weights()`
- `compute_weighted_balance()`
- `compute_effective_sample_size()`
- `evaluate_comparator_weighting()`

These functions provide a weighting-based approach for balancing multiple treatment groups and evaluating weighted covariate balance and effective sample size.

```r
report <- evaluate_comparator_weighting(
  sample_4group,
  method = "overlap",
  X_vars = covariates,
  treatment_var = "treatment",
  anchor_level = anchor
)

report$balance$summary
report$ess
report$n_trimmed
```

Because the GPS model is fitted without regularization, near-separation can drive a propensity score toward zero and produce an unbounded IPTW weight. `compute_balancing_weights()` therefore bounds the scores at `trim` (default `1e-3`) before dividing by them, and reports `n_trimmed` — the number of subjects with at least one score below the floor. A nonzero `n_trimmed` is a positivity warning about the data, not about the weights. Set `trim = 0` to disable it.

---

## Main Functions

| Function | Description |
|---|---|
| `estimate_gps_multinom()` | Estimate generalized propensity scores using multinomial logistic regression |
| `gps_candidate_search()` | Identify candidate matches in GPS space |
| `sam_match()` | Perform Shared Anchor Matching |
| `sam_evaluate()` | Evaluate matching quality and covariate balance |
| `extract_matched_data()` | Extract the matched subject-level cohort |
| `sam_estimate_effects()` | Estimate treatment effects in the matched cohort |
| `match_3way()` | Perform three-group propensity score matching |
| `compute_smd_balance()` | Standardized mean differences in the matched cohort |
| `compute_pairwise_treatment_auc()` | Pairwise treatment-discrimination AUCs |
| `get_pooled_covariance()` | Pooled within-group covariance and its inverse |
| `compute_balancing_weights()` | Compute multi-arm balancing weights |
| `compute_weighted_balance()` | Assess weighted covariate balance |
| `compute_effective_sample_size()` | Calculate effective sample size after weighting |

---

## Example Data

Two example datasets are included:

- `sample_4group`: four-treatment-group example data for Shared Anchor Matching
- `sample_3group`: three-treatment-group example data for three-way matching

Dataset documentation can be accessed in R:

```r
?sample_4group
?sample_3group
```

---

## Citation

If you use SAM-R in your research, please cite:

> Hong S, Hong S, Lee KH, Cha N. **SAM-R: Shared Anchor Matching**. Version 0.2.0. Zenodo. https://doi.org/10.5281/zenodo.21926956

A formal citation for the associated methodological paper will be added upon publication.

---

## License

This package is distributed under the MIT License.
