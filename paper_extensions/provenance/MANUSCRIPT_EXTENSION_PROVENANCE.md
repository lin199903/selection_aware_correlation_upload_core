# Manuscript-extension provenance

This public provenance record covers the liver-only MASLD–AH replay series, empirical-covariance recovery grid, independent cross-cohort evaluation and the synchronized `selcorr` 0.1.1 implementation. These analyses were closed without adding a new cohort, metric, sensitivity grid or package feature after result inspection.

- MASLD–AH liver-only replay: 5 workflows, 2,000 reference draws each.
- Empirical-covariance recovery: fixed 4 × 4 grid (effect size 0.2/0.4/0.6/1.0 SD × fraction 5/10/25/50%); 30 outer replicates per cell; 3,000 reference draws per cell in the archived full run; rejection rate 0 in all 16 cells.
- Independent cross-cohort evaluation: 2,000 evaluation-side draws; r = 0.2371, reference mean = 0.2279, percentile = 62.2, upper-tail P = 0.3783.
- `selcorr`: version 0.1.1; 32 tests/143 expectations in the verified closure; R CMD check status OK.

Internal manuscript draft numbers are intentionally omitted from the public path and title of this provenance record. Historical draft mapping is retained only in `SUBMISSION_SNAPSHOT.json` and `RELEASE_METADATA.json`.
