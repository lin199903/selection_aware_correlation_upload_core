# PACKAGE_PROVENANCE_MAP

`selcorr` is a refactor of analysis code that has already been executed and
archived. It is not a reimplementation of the algorithms from memory, and it
does not change any frozen statistical definition. This map records, for each
package function, the production script it came from, the algorithmic
equivalence, the known differences, and how the equivalence was verified.

## Core reference construction

| Production script | Package function | Algorithmic equivalence | Known differences | Verification |
|---|---|---|---|---|
| `reanalysis_20260722/R/11_osa_phenotype_null.R` (OSA–MASLD anchored replay) | `selcorr(method = "replay")` → `replay_anchored()` | Same mechanics: Freedman–Lane residual permutation of the replayed arm; nuisance model fitted once and held; whole-sample residual vectors permuted; full model refitted on every draw with `limma::eBayes(trend = TRUE, robust = TRUE)`; differential-expression filtering, eligibility reconstruction, same-sign retention and the selected Pearson correlation re-evaluated per draw; `hub_de` eligibility = `(hub_replayed AND DE_anchored) OR (hub_anchored AND DE_replayed)`; minimum four selected genes; Y-chromosome genes excluded; upper-tail add-one p value | The script hard-codes CSV/RDS paths and writes production outputs; the package takes the anchored table and the replayed arm as inputs and writes only user-requested checkpoints. The network step remains fixed in both (hub sets are conditioning inputs, not recomputed). | `tests/testthat/test-frozen-regression.R`: archived draws reproduce the frozen reference mean (0.668) and p value (0.922); the slow integration test recomputes the observed tuple (r = 0.551, eligible 110, selected 52) from the archived inputs |
| `validation_analyses/scripts/03_masld_samedisease_replay.R` (MASLD–MASLD anchored replay) | `selcorr(method = "replay")` with `selection_spec("de_anchor")` | Same mechanics: eligible pool anchored to the fixed arm's DE genes (`|log2FC| > 0.3`, `FDR < 0.05`); the permuted arm contributes only sign and effect values; identical correlation guard and add-one tail | Cohort-specific design formulas (`~ group`, `~ sex + age_z + group`) are inputs in the package rather than hard-coded | Archived MASLD–MASLD draws reproduce the stored null mean and p value exactly (tolerance 1e-6), and the frozen forward values 0.651 / 0.734 / 0.845 |
| `learned/v151_9_MASLD_AH/11_masld_ah_replay.R` (MASLD–AH series) | `selcorr(method = "replay")` | Same evaluation machinery as above; the package adds no new rule | Cohort loaders (metadata parsing, cirrhosis exclusion) stay with the user; the package receives prepared arms | Archived MASLD–AH draws reproduce the stored null mean, SD and p value for all five workflows (tolerance 1e-6) |

## Diagnostic comparator

| Production script | Package function | Algorithmic equivalence | Known differences | Verification |
|---|---|---|---|---|
| Fixed-set comparator described in the manuscript Methods (size-matched subsets sampled without replacement from the fixed paired-effect universe; add-one upper-tail p value) | `replay_fixed_set()` (via `selcorr(method = "fixed_set")`) | Same construction and tail | Draw count is user-specified; result is flagged `diagnostic_only = TRUE` because selection is not regenerated | `tests/testthat/test-engine.R` |

## Analytic geometry

| Source | Package function | Algorithmic equivalence | Verification |
|---|---|---|---|
| Manuscript Results and Supplementary Note 15 (same-sign selection; symmetric magnitude truncation with `a = Φ⁻¹(1 − q/2)`, `λ(a) = φ(a)/[1 − Φ(a)]`, `ρ₀(q) = λ(a)² / (1 + a λ(a))`) | `same_sign_population_correlation()`, `analytic_selected_correlation()` | Direct transcription of the derived expressions | `tests/testthat/test-analytic.R`: `2/π`; `q = 0.90/0.75/0.50/0.25/0.10 → 0.6964/0.7732/0.8698/0.9371/0.9685` (tolerance 5e-4); monotonicity; `q = 1` reduces to `2/π` |

## Recovery diagnostic

| Production script | Package function | Algorithmic equivalence | Known differences | Verification |
|---|---|---|---|---|
| `R/13_osa_spike_in_recovery.R` | `diagnose(type = "recovery")` | Same injection mechanics: direction-aligned shared effects (`+δ/2` in cases, `−δ/2` in controls, multiplied by the anchored arm's effect sign), a fresh gene draw per outer replicate, Freedman–Lane null on the spiked matrix, rejection at α = 0.05, Monte Carlo standard error reported | The grid (fractions × effect sizes) and replicate counts are parameters rather than fixed constants; the production constant grid used 25 genes and four doses | Logic and column contract tested in `test-checkpoint-and-diagnostics.R` |

## Independent evaluation

| Production script | Package function | Algorithmic equivalence | Known differences | Verification |
|---|---|---|---|---|
| `learned/v151_9_MASLD_AH/14_cross_cohort_selection_evaluation.R` | `validate()` | Same semantics: the selected gene set comes from the original selection data, the evaluation data take no part in selection, the gene set is externally fixed in evaluation, the evaluation reference does not replay the original selection, and the two p values are never combined | Cohort pairs are user inputs rather than hard-coded | Archived cross-cohort result (evaluation r = 0.2371, reference centre = 0.2279, upper-tail p = 0.3783) is reproducible from the archived inputs |

## Deliberate non-goals in 0.1

* No composite score, no new statistic, no renamed existing statistic.
* No second "framework": `diagnose()` and `validate()` are subordinate helpers.
* No reimplementation of differential-expression estimation; the package calls
  the same estimator settings as the archived workflows.
* No compiled code: profiling has not shown a pure-R inner loop to be the
  bottleneck, and `limma` dominates the cost.
