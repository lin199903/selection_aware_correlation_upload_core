# Reproducibility — reproducibility-v1.0.0-rc3

## Level 1: integrity and claim audit
Run `python verify_release.py`. It checks the paper-extension tables, package-version markers, release-boundary metadata and current snapshot contract.

## Level 2: verified historical base
The verified historical candidate and its release metadata remain preserved under `archive/v151_8_release_metadata/`.

## Level 3: R package
`selcorr` 0.1.1 was built and checked in a fresh R 4.5.1 environment during RC2 assembly. The verified source artifact and R CMD check log are retained under `selcorr_release/`. No statistical definition is changed by the repository-versioning update.

## Level 4: third-party source-data re-execution
Retrieve external source data under their own terms. Core scripts remain public. Late MASLD–AH cohort-specific loader/parsing scripts are outside this public release boundary; the derived outputs and package-level replay machinery are public under `paper_extensions/` and `selcorr/`.
