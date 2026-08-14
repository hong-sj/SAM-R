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
