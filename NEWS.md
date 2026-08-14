# SAMatch 0.2.0

A correctness, performance and validation pass over the whole package. Results
on the bundled `sample_4group` and `sample_3group` datasets are unchanged
throughout: every matching decision, distance, balance value, weight and
outcome contrast is identical to 0.1.0, apart from the new columns and the one
sign change listed under "Breaking changes".

## Breaking changes

* Input that was previously accepted and silently mishandled now raises. This
  is the point of the release, but it does mean a script that was quietly
  producing a wrong answer will now stop:
  - a covariate containing a missing or non-finite value (#2),
  - a `data` frame re-sorted or filtered between pipeline stages (#3),
  - `top_m` or `top_n` that is not a positive integer (#5),
  - an `anchor_level` that matches no rows (#1),
  - a matched cohort with incomplete matched sets (#13).
* `compute_weighted_balance()` now reports SMDs as anchor minus comparator, the
  same orientation as `compute_smd_balance()`. Only the sign of
  `by_covariate$smd` changes; `abs_smd` and the summary are unaffected (#10).
* `compute_balancing_weights()` gains `trim = 1e-3`, which bounds the
  propensity scores before dividing by them. Pass `trim = 0` for the previous
  behaviour. The bundled data trims nothing (#9).
* `match_3way()` renames `ps_space` to `gps_space`, matching
  `gps_candidate_search()`. `ps_space` still works and warns (#8).
* `match_3way(X_vars=)` now warns instead of silently discarding the argument
  (#8).
* New columns: `smd_defined` on both balance tables, `abs_smd` on
  `compute_smd_balance()`, `n_undefined` on both summaries (#4, #10),
  `n_trimmed` on the weighting output (#9), `max_possible_rate` on both
  matching engines (#8), and `separation` on the outcome contrasts (#14).
* `RANN` is a new dependency (#6).
* `gps` is validated at every stage that accepts it, and `gps_candidate_search()`
  records a `gps_fingerprint` that later stages verify.

## Bug fixes

* A single missing covariate value no longer silently empties the match.
  `solve()` returns an all-`NA` matrix rather than raising, so the
  singular-matrix fallback was dead code for missing data and the matching rate
  fell from 0.132 to 0 without an error or a warning (#2).
* Re-sorting `data` between stages is rejected. Matched sets store positional
  row indices, so a reorder silently repointed them at other subjects and the
  balance report described a cohort that was never matched (#3).
* A zero-variance covariate no longer turns the whole balance summary into
  `NaN`. It is reported as unassessable and excluded from the mean and maximum
  (#4).
* A numeric `anchor_level` names a level rather than a position.
  `stats::relevel()` reads a numeric `ref` as an index into `levels()`, so
  `anchor_level = 0` raised for a level that was plainly present (#1).
* Covariate names that are not syntactic R names survive GPS fitting (#11).
* The outcome analysis drops whole matched sets, not individual rows, when a
  value is missing (#13).
* Complete outcome separation is detected, warned about, and flagged (#14).
* An unbounded IPTW weight from near-separation is trimmed; the largest weight
  in a separated example falls from 9,478 to 334 (#9).
* A GPS matrix that does not describe `data` is rejected rather than silently
  producing a different match: values outside `[0, 1]`, non-finite values, rows
  that do not sum to one, duplicated or unnamed columns, a row order that
  disagrees with `data`, and a treatment level with no GPS column. That last
  one previously dropped the level out of `groups`, so 59 subjects on the
  bundled data were never matched and nothing said so.
* Handing a *different* GPS to `sam_evaluate()` or `match_3way()` than the one
  the candidate pools were selected against is rejected.
* A missing treatment label is rejected. It previously compared as `NA` against
  every group, so the subject silently vanished from all of them.
* The data fingerprint now digests every column recorded at search time, so
  rewriting the contents of existing rows is caught as well as reordering them.
  Adding a column between stages remains legal.
* `compute_weighted_balance()` rejects weights that are non-finite, negative or
  the wrong length, and a treatment group whose weights sum to zero — which
  previously produced a table of `NaN`.
* `calc_caliper_3way()` requires exactly three groups, finite scores and at
  least two subjects per group; `stats::var()` of a single observation is `NA`,
  which reached the caliper silently. `match_3way()` requires a finite caliper
  and reads a numeric `reference_level` as a level name.

## Performance

* `gps_candidate_search()` queries a k-d tree instead of building the full
  anchor-by-comparator distance matrix. At n = 20,000 this is 5.32 s / 1224 MB
  to 0.20 s / 33 MB; n = 100,000, previously not feasible, takes 1.12 s /
  112 MB. Tie handling is preserved exactly (#6).
* `sam_match()` maintains a reverse index of which anchors claim which
  comparator subject, instead of rescanning every remaining anchor after each
  match. At n = 20,000 this is 2.84 s to 1.63 s (#7).
* `match_3way()` tombstones popped candidates rather than deleting them from
  six parallel vectors (#7).
* `sam_match()` builds one reverse index per comparator group instead of
  `match()`-ing every anchor's candidates against the whole group, which
  rebuilt a hash table of the group on every call and accounted for 74% of the
  function's runtime.
* `sam_match()` computes the Mahalanobis distances for all screened pairs in
  one pass rather than one call per anchor. Distances move by at most 2.3e-14
  relative, from BLAS summation order; no matching decision changes.
* `sam_match()` selects the cheapest anchor from a shortlist of the smallest
  losses rather than scanning every remaining anchor. Losses only ever rise, so
  an anchor outside the shortlist cannot drop below it; this removes the last
  quadratic term.
* `extract_matched_data()` assembles the subject-level cohort in one pass
  instead of building and row-binding one data frame per matched set.

Together, on synthetic data with four balanced groups:

    n = 50,000    sam_match             12.05 s -> 1.44 s
                  extract_matched_data   1.32 s -> 0.03 s
    n = 200,000   sam_match             13.34 s -> 5.76 s

Time per anchor in `sam_match()` is now flat at roughly 115 us across
n = 50,000 to 200,000, where it previously grew with the sample.

## Documentation and infrastructure

* `compute_smd_balance()`, `compute_pairwise_treatment_auc()` and
  `get_pooled_covariance()` gain defaults and infer their groups, and are
  documented in the README, which previously listed none of them (#8).
* The test suite grows from 6 blocks to per-issue regression tests, coverage
  for all 20 exported functions, reference implementations for the candidate
  search and the matching engine, hand-computed values, and golden fixtures
  with a tiered comparison (#12).
* `R CMD check` now runs in CI across R devel, release and oldrel-1 on Ubuntu,
  plus macOS and Windows (#12).
