# Final upload checklist

- [x] Verify `selcorr` 0.1.1 build/check artifact (`selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz`; R CMD check Status: OK).
- [x] Remove the archive-dependent hard-coded local path from package tests.
- [x] Decouple public repository paths and release identity from internal manuscript versions.
- [x] Record the current v153 manuscript/SI hashes as provenance-only snapshot metadata.
- [ ] Change release version from `reproducibility-v1.0.0-rc3` to `reproducibility-v1.0.0` at final release.
- [ ] Add final GitHub release URL / archival DOI to citation metadata when available.
- [ ] Run `python verify_release.py` from a fresh checkout.
- [ ] Create final SHA-256 manifest and release archive from a fresh git checkout.
