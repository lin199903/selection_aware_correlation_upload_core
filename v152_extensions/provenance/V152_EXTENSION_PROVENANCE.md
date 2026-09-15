# v152 extension provenance

The v152 extension was closed on 2026-09-14 without adding a new cohort, metric, sensitivity grid or package feature after result inspection.

- MASLD–AH liver-only replay: 5 workflows, 2,000 reference draws each.
- Empirical-covariance recovery: fixed 4 × 4 grid (effect size 0.2/0.4/0.6/1.0 SD × fraction 5/10/25/50%); 30 outer replicates per cell; 3,000 reference draws per cell in the archived full run; rejection rate 0 in all 16 cells.
- Independent cross-cohort evaluation: 2,000 evaluation-side draws; r = 0.2371, reference mean = 0.2279, percentile = 62.2, upper-tail P = 0.3783.
- `selcorr` release state recorded by the original closure: version 0.1.1, 32 tests/143 expectations, 0 failures/errors, one slow-test skip, `R CMD check` status OK.

The original verified `selcorr_0.1.1.tar.gz` was not mounted in the environment used to assemble this GitHub RC. See `selcorr_release/RECONSTRUCTION_STATUS.md`.
