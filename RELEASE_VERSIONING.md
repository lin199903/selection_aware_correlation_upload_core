# Release versioning policy

This repository uses **public reproducibility-release versions**, independent of internal manuscript draft numbers.

## Public tags

- Release candidate: `reproducibility-v1.0.0-rc1`
- First stable public release: `reproducibility-v1.0.0`
- Later repository-only corrections that do not alter reported scientific results: patch releases, e.g. `reproducibility-v1.0.1`
- Additions that materially expand the reproducibility archive while preserving the paper's scientific identity: minor releases, e.g. `reproducibility-v1.1.0`

Internal manuscript identifiers such as `v152`, `v153`, etc. are provenance metadata only and are **not** used as repository or release identities.

## Software package version

The R package retains its own independent semantic version (`selcorr` 0.1.1). Repository-release versions and package versions must not be conflated.

## Manuscript snapshot mapping

A release may record which internal manuscript snapshot it accompanied, but later wording-only manuscript revisions do not require renaming the repository or invalidating the release.
