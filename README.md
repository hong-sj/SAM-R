# SAM

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
head(matched$matched)
```

The resulting object contains the matched sets, group-specific Mahalanobis distances, total matching loss, unmatched anchor subjects, and overall matching rate.

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
  treatment_var = "treatment"
)

matched3$matching_rate
head(matched3$matched)
```

---

## Multi-Arm Weighting

SAM also includes functions for multi-arm comparator weighting:

- `compute_balancing_weights()`
- `compute_weighted_balance()`
- `compute_effective_sample_size()`
- `evaluate_comparator_weighting()`

These functions provide a weighting-based approach for balancing multiple treatment groups and evaluating weighted covariate balance and effective sample size.

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

If you use SAM in your research, please cite the associated methodological work.

A formal citation will be added upon publication.

---

## License

This package is distributed under the MIT License.