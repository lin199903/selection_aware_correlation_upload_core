# Public file map — reproducibility-v1.0.0-rc2

| Scientific object | Public location | Re-execution status |
|---|---|---|
| Analytic selected-correlation geometry | `validation_analyses/` + `selcorr/R/analytic.R` | public code |
| OSA–MASLD full-selection replay | `reanalysis_20260722/R/11_osa_phenotype_null.R` + frozen outputs | public code; source data external |
| MASLD–MASLD stress test | core script/output tree | public code/output; some reconstructed input objects not redistributed |
| Null-location decomposition | `validation_analyses/outputs/null_location_decomposition.csv` | public derived output |
| External OSA transportability | `validation_analyses/outputs/external_cross_disease/` | public derived output |
| MASLD–AH 5-workflow series | `v152_extensions/source_data/masld_ah_replay_summary.csv` | derived output public; cohort-specific loaders are outside the published scope (see PUBLIC_RELEASE_SCOPE.md) |
| MASLD–AH 16-cell recovery | `v152_extensions/source_data/masld_ah_recovery_summary.csv` | derived output public; `selcorr::diagnose()` implements recovery machinery |
| Independent evaluation | `v152_extensions/source_data/cross_cohort_independent_evaluation_summary.csv` | derived output public |
| `selcorr` | `selcorr/`, `selcorr_release/` | verified 0.1.1 build: `R CMD check` Status: OK in a fresh R 4.5.1 environment (log in `selcorr_release/`); artifact `selcorr_0.1.1_VERIFIED_RC2.tar.gz` |
| Submission text | not redistributed | hashes/status in `SUBMISSION_SNAPSHOT.json` |
