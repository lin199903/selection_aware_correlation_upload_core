# selcorr 0.1.1 — verified release summary

The original project closure records the following completed state on 2026-09-14:

- version 0.1.1;
- 7 test files, 32 tests, 143 expectations, 0 failures, 0 errors, 1 slow-test skip;
- serial/parallel determinism re-verified;
- checkpoint/resume re-verified;
- `R CMD check --no-manual selcorr_0.1.1.tar.gz`: 0 errors, 0 warnings, 0 notes;
- source tarball size recorded as 63,597 bytes;
- recovery fixture `inst/extdata/masld_ah_recovery_summary.csv`;
- no API or core statistical-definition change from 0.1.0.

The original verified 0.1.1 tarball is not present in this RC build environment. The included `selcorr/` tree was reconstructed from the verified 0.1.0 source plus the documented 0.1.1 delta. It must either be replaced by the original verified tarball/source snapshot or re-run through `R CMD build` and `R CMD check` before the final GitHub tag.
