# selcorr 0.1.1 build status (rc2)

**Status: closed.** The 0.1.1 source tree in this release was built and checked in a fresh
R environment during rc2 assembly.

| Item | Value |
|---|---|
| Environment | R 4.5.1 (2025-06-13 ucrt), Windows 11 x64 |
| `R CMD build selcorr` | `selcorr_0.1.1.tar.gz` produced |
| `R CMD check --no-manual` | **Status: OK** (0 errors, 0 warnings, 0 notes) |
| Log | `selcorr_release/R_CHECK_LOG_v0.1.1_RC2.txt` |
| Artifact | `selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz` |
| Test suite | 32 tests / 143 expectations, 0 failures, 0 errors, 1 skip (slow integration test) |

## Delta from the 0.1.0 baseline

1. Version 0.1.0 → 0.1.1.
2. Added `tests/testthat/test-recovery-regression.R` and the fixture
   `inst/extdata/masld_ah_recovery_summary.csv` (16-cell recovery grid).
3. `NEWS.md` records both.
4. **Portability fix (rc2 only):** `tests/testthat/test-frozen-regression.R` no longer
   carries a hard-coded default archive path. The archive-dependent tests read
   `SELCORR_ARCHIVE_ROOT` from the environment and skip when it is unset, so the
   package no longer embeds any local path.

No file under `R/`, `man/`, `NAMESPACE` or `vignettes/` differs from the verified
0.1.0 source (`CORE_SOURCE_IDENTITY.json`), and no statistical definition changed.
