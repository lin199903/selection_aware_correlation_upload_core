# GitHub release audit — reproducibility-v1.0.0-rc3

## Purpose

Keep the public repository/release identity independent of manuscript draft numbering while synchronizing the corrected independent-evaluation null and the current v153 provenance snapshot. The paper's central methodological claim and the `selcorr` API are unchanged.

## Version-decoupling changes

- `v152_extensions/` was renamed to version-neutral `paper_extensions/`.
- Public README, file map, scope, reproducibility notes, analysis-stage index and verifier use version-neutral paths.
- `RELEASE_METADATA.json` and `SUBMISSION_SNAPSHOT.json` retain internal manuscript versions only as provenance metadata.
- `CITATION.cff` points to the current repository and the `reproducibility-v1.0.0-rc3` release identity.

## Scientific correction synchronized in rc3

1. **Independent-evaluation null corrected.** The earlier evaluation construction retained the observed group contrast in the fitted component and produced a shifted, artificially narrow reference (mean 0.2279, SD 0.0285; upper-tail P = 0.3783). The corrected intercept-only evaluation null is centred near zero (mean 0.0034, SD 0.1621); the observed r = 0.2371 lies at the 92.7th percentile (upper-tail P = 0.0735) and does not reach the prespecified one-sided significance threshold. Public source data and package provenance use the corrected construction.
2. **GSE135251 estimator sensitivity classified correctly.** The DESeq2 sensitivity had been executed previously; the issue was missing archived execution outputs, not a missing analysis. Re-execution reproduced 79 versus 78 selected genes, Jaccard 0.8690, and selected r = 0.6507 versus 0.6296. The public repository carries the summary; full execution outputs remain in the full author/release archive tier.
3. **Manuscript snapshot advanced to v153.** The repository records the author-designated v153 manuscript and SI only as provenance hashes. Public repository identity remains `reproducibility-v1.0.0-rc3` and is independent of manuscript draft numbering.

## Scientific invariants

- MASLD–AH replay summary remains 5 workflows.
- Recovery grid remains 16 cells with rejection rate 0 in every archived cell.
- Strongest recovery cell remains median observed r = 0.7430 and reference mean = 0.6687.
- `selcorr` remains version 0.1.1; no statistical API was changed by this repository synchronization.

## Release boundary

This is an RC, not the final `reproducibility-v1.0.0` release. Final release URL/DOI and a fresh-checkout final archive remain production steps.
