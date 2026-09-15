# Final upload checklist

- [x] Replace/recheck the reconstructed `selcorr` 0.1.1 RC with a freshly built and checked artifact
      (`selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz`; `R CMD check` Status: OK; log in `selcorr_release/`).
- [x] Remove the hard-coded archive path from the archive-dependent tests (portability defect found in review).
- [x] Synchronize final v152 SI and record its file name, SHA-256 and status in `SUBMISSION_SNAPSHOT.json`.
- [ ] Change release version from `reproducibility-v1.0.0-rc2` to `reproducibility-v1.0.0` in README/CITATION/snapshot.
- [ ] Add final GitHub release URL to CITATION/README.
- [ ] Run `python verify_release.py` from a fresh checkout.
- [ ] Create final SHA256 manifest and release ZIP from a fresh git checkout.
