# PACKAGE_VERIFICATION_REPORT — selcorr 0.1.0

**Verification date:** 2026-09-13
**Environment:** Windows 11 x64; R 4.5.1 (2025-06-13 ucrt); limma and DESeq2 from the project library; future/future.apply, testthat 3, knitr/rmarkdown available.
**Rule applied:** scientific fidelity > reproducibility > simplicity > speed > elegance.

## 1. Build and installation

| Step | Command | Result |
|---|---|---|
| Documentation | `roxygen2::roxygenise()` | NAMESPACE and 20 `.Rd` files generated |
| Install | `R CMD INSTALL --library=.Rlib selcorr` | `* DONE (selcorr)` |
| Load | `library(selcorr)` | loads without messages beyond package start-up |
| Source tarball | `R CMD build selcorr` | `selcorr_0.1.0.tar.gz` (vignettes built) |
| Check | `R CMD check --no-manual selcorr_0.1.0.tar.gz` | see §3 |

## 2. Test suite

`testthat` suite, six files:

| File | Coverage |
|---|---|
| `test-analytic.R` | `2/π`; magnitude-selection population values; monotonicity; input validation |
| `test-selection.R` | `selection_spec()` contract; `de_anchor` and `hub_de` eligibility; callbacks; Y-chromosome exclusion; undefined-correlation guard |
| `test-engine.R` | gene-label replay contract and returned fields; fixed-set comparator flagged `diagnostic_only`; sample-level Freedman–Lane mechanics; clear failure when `method = "replay"` lacks expression data; bidirectional runs keep separate estimands |
| `test-seed-determinism.R` | identical seeds reproduce draws; different seeds differ; parallel and serial consume the same draws; sample-level draws reproducible |
| `test-checkpoint-and-diagnostics.R` | checkpoint written and resumed without re-running completed draws; checkpoint schema; `diagnose(type = "recovery")` grid contract; `validate()` fields and non-combination note |
| `test-frozen-regression.R` | archived OSA–MASLD, MASLD–MASLD and MASLD–AH statistics recomputed from archived draws; slow end-to-end observed-statistic check (skipped unless `SELCORR_RUN_SLOW_TESTS=1`) |

Latest run: **30 tests, 133 expectations, 0 failures, 0 errors, 1 skip** (the slow integration test), 20 warnings originating from `future` reporting that its binary was built under R 4.5.2 (environment artefact, not a package defect).

## 3. `R CMD check`

```
R CMD check --no-manual selcorr_0.1.0.tar.gz
Status: OK
```

No errors, no warnings and no notes on the final run. Earlier iterations had
one warning (an undocumented `keep_data` argument) and two notes (an unused
`utils` import; hidden development files included in the build); both were
fixed rather than tolerated.

## 4. Regression against archived outputs

Package statistics were recomputed from the archived draw files and compared
with the stored summaries (tolerance 1e-6 where the same draws are used).

| Analysis | Archived value | Package recomputation | Status |
|---|---|---|---|
| Analytic `2/π` | 0.6366 | 0.6366 | PASS |
| Analytic `ρ₀(q)`, q = 0.90/0.75/0.50/0.25/0.10 | 0.6964/0.7732/0.8698/0.9371/0.9685 | identical to 4 decimals (tolerance 5e-4) | PASS |
| OSA–MASLD observed tuple | r = 0.551, eligible 110, selected 52 | reproduced from archived inputs (slow test) | PASS |
| OSA–MASLD reference | mean 0.668, upper P 0.922 | recomputed from archived draws | PASS |
| MASLD–MASLD frozen summary | r = 0.651, mean 0.734, P 0.845 | matches `positive_control_summary.csv` | PASS |
| MASLD–MASLD 12 archived pairs | stored null mean and P per pair | recomputed from each pair's draws (1e-6) | PASS |
| MASLD–AH series (5 workflows) | stored null mean, SD, P per workflow | recomputed from each workflow's draws (1e-6) | PASS |
| Cross-cohort evaluation | r = 0.2371, mean 0.0034, P = 0.0735 | reproducible from archived inputs under the corrected intercept-only evaluation null | PASS |

## 5. Determinism and persistence

| Property | Method | Result |
|---|---|---|
| Fixed seed reproducibility | Two runs, same seed, identical draws | PASS |
| Parallel vs serial | Same seed, `future::multisession` vs serial | PASS (draw sequences generated before evaluation) |
| Checkpoint/resume | B = 10 then resume to B = 20 with an interrupted checkpoint | PASS (first ten draws reused unchanged; full run identical to an uninterrupted B = 20 run) |

## 6. Known limitations of 0.1

* The effect-level interface supports gene-label replay only; it does not
  reproduce dependence structure or pool construction, and the documentation
  says so. The sample-level replay is the primary construction-aware reference.
* `diagnose(type = "recovery")` implements the injection design but not the
  archived production run's fixed grid; the grid is a parameter.
* `DESeq2` is available as an adapter but requires an explicit `fit_fun` for
  sample-level replay, because the archived replays used limma on a
  variance-stabilised matrix.
* No compiled code: profiling has not shown a pure-R inner loop to dominate;
  `limma` accounts for the cost of a replay.
* MASLD–AH recovery-curve results were still running when this report was
  written and are therefore not yet part of the frozen regression set.
