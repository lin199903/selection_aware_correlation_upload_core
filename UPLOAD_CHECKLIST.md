# Final upload checklist

- [x] Verify `selcorr` 0.1.1 build/check artifact (`selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz`; R CMD check Status: OK).
- [x] Remove the archive-dependent hard-coded local path from package tests.
- [x] Decouple public repository paths and release identity from internal manuscript versions.
- [x] Correct the independent-evaluation null to the intercept-only construction and synchronize the public summary/provenance files.
- [x] Classify the GSE135251 DESeq2 sensitivity correctly as executed; retain the public summary and full execution outputs in the full archive tier.
- [x] Record the current `v153` manuscript and SI SHA-256 values as provenance-only snapshot metadata, with explicit statement that the manuscript bytes are not redistributed in this repository.
- [ ] Change release version from `reproducibility-v1.0.0-rc3` to `reproducibility-v1.0.0` at final release.
- [ ] Add final GitHub release URL / archival DOI to `CITATION.cff` when available.
- [ ] Run `python verify_release.py` from a fresh checkout immediately before final tagging.
- [ ] Create the final SHA-256 manifest and release archive from that fresh checkout.
