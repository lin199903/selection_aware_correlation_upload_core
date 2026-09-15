# Reproducibility — reproducibility-v1.0.0-rc1

## Level 1: integrity and claim audit
Run `python verify_release.py`. It checks the v152 extension tables, package-version markers, absence of final-release claims while hard stops remain, and SHA-256 manifest integrity.

## Level 2: verified historical base
Before extension, the v151.8 public candidate passed its own `verify_manifest.py` byte-level/provenance checks. Its original release metadata and verifier are preserved under `archive/v151_8_release_metadata/`.

## Level 3: R package
The original project closure records a verified `selcorr` 0.1.1 build with 32 tests/143 expectations and `R CMD check` status OK. The exact verified 0.1.1 tarball was not available in this RC build environment. The included source is a reconstruction from the verified 0.1.0 source plus documented 0.1.1 delta; no core `R/`, `man/`, `NAMESPACE` or vignette file changed. See `selcorr_release/RECONSTRUCTION_STATUS.md`.

## Level 4: third-party source-data re-execution
Retrieve external source data under their own terms. The v151.x core scripts remain public. Late MASLD–AH cohort-specific loader/parsing scripts are not present in this RC; v152 derived outputs and package-level replay machinery are public.
