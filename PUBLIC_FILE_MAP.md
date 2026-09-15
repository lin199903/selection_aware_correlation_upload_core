# Public file map — reproducibility-v1.0.0-rc3

| Scientific object | Public location | Re-execution status |
|---|---|---|
| Analytic selected-correlation geometry | `validation_analyses/` + `selcorr/R/analytic.R` | public code |
| OSA–MASLD full-selection replay | `reanalysis_20260722/R/11_osa_phenotype_null.R` + frozen outputs | public code; source data external |
| MASLD–MASLD stress test | core script/output tree | public code/output; some reconstructed input objects not redistributed |
| Null-location decomposition | `validation_analyses/outputs/null_location_decomposition.csv` | public derived output |
| External OSA transportability | `validation_analyses/outputs/external_cross_disease/` | public derived output |
| MASLD–AH 5-workflow series | `paper_extensions/source_data/masld_ah_replay_summary.csv` | derived output public; cohort-specific loaders are outside the published scope |
| MASLD–AH 16-cell recovery | `paper_extensions/source_data/masld_ah_recovery_summary.csv` | derived output public; `selcorr::diagnose()` implements recovery machinery |
| Independent evaluation | `paper_extensions/source_data/cross_cohort_independent_evaluation_summary.csv` | derived output public |
| Figure 4 source layer | `paper_extensions/source_data/figure4_reference_position_source_data.csv` | public derived output |
| `selcorr` | `selcorr/`, `selcorr_release/` | verified 0.1.1 build; `R CMD check` Status: OK in fresh R 4.5.1 environment |
| Submission text | not redistributed | hashes/status in `SUBMISSION_SNAPSHOT.json` |
