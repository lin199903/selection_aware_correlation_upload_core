# GitHub release audit — reproducibility-v1.0.0-rc3

## Purpose

Decouple the public repository/release identity from internal manuscript draft numbering without changing any scientific result, statistical definition or `selcorr` API.

## Changes audited

- `v152_extensions/` renamed to version-neutral `paper_extensions/`.
- Public README, file map, scope, reproducibility notes, analysis-stage index and verifier updated to version-neutral paths.
- `RELEASE_METADATA.json` and `SUBMISSION_SNAPSHOT.json` retain `v153` only as provenance metadata.
- `CITATION.cff` now points to the current repository and uses `reproducibility-v1.0.0-rc3`.
- Source-data blobs were moved without content changes.
- Extension provenance was renamed and rewritten to remove public draft-number identity.

## Scientific invariants

- MASLD–AH replay summary remains 5 workflows.
- Recovery grid remains 16 cells with rejection rate 0 in every archived cell.
- Strongest recovery cell remains median observed r 0.7430 and reference centre 0.6687.
- `selcorr` remains version 0.1.1; no statistical code or API was changed in this release-identity update.

## Release boundary

This is an RC, not the final `reproducibility-v1.0.0` release. Final release URL/DOI and a fresh-checkout final archive remain production steps.
