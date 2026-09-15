# GitHub release audit — reproducibility-v1.0.0-rc1

## Overall status

**RC PASS / FINAL HOLD.** The repository is suitable for author/reviewer inspection and GitHub staging, but should not yet be tagged as the final public `reproducibility-v1.0.0` release.

## Passed

- Verified v151.8 base re-tested before staging: PASS.
- No new scientific analysis introduced during packaging.
- v152 MASLD–AH five-workflow summary present.
- v152 16-cell recovery grid present; all 16 cells match the frozen reported values.
- Independent cross-cohort evaluation present.
- Figure 4 machine-readable source layer present.
- `selcorr` source tree included; core R/man/NAMESPACE/vignettes byte-identical to verified 0.1.0: True.
- No raw third-party biological data added by the v152 extension.
- Internal manuscript/reviewer/decision ledgers are excluded from the public root.

## Hard stops before final tag

1. **Final v152 Supplementary Information synchronization.** Current `SUBMISSION_SNAPSHOT.json` deliberately marks the available SI as v151.8-era and pending v152 sync.
2. **Exact `selcorr` 0.1.1 artifact.** The original project record confirms a verified tarball (63,597 bytes; R CMD check OK), but that artifact is not mounted here. Replace `selcorr_release/selcorr_0.1.1_RECONSTRUCTED_RC.tar.gz` with the original verified tarball, or re-build/re-check the reconstructed source in R before final tag.

## Non-blocking production items

- Update release URL/tag after pushing to GitHub.
- Final AH cohort bibliography/source-publication cross-check in the manuscript/SI.
- Optionally refresh figure binaries if journal production changes labels only; source data are already locked.

## No-analysis rule

Closing these items requires no new cohort, endpoint, threshold, sensitivity analysis or statistic.

## Post-build packaging checks

- Public-surface scan: 370 non-archive files scanned; 0 credential/secret-pattern hits; 0 internal-risk filenames outside `archive/`.
- Fresh-extraction verification: PASS (`python verify_release.py`).
- The release ZIP is reproducible as a staging artifact and was verified after fresh extraction.
- Final RC ZIP SHA-256 is distributed alongside the archive in `reproducibility-v1.0.0_GITHUB_RELEASE_CANDIDATE_RC1.sha256`.
