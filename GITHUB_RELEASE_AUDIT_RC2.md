# GitHub release audit — reproducibility-v1.0.0-rc2

## Overall status

**RC PASS / FINAL TAG HOLD (one step).** The two rc1 hard stops are closed; only the
release-identifier bump, manifest regeneration from a fresh checkout and the final release
URL remain.

## Closed since rc1

| rc1 hard stop | How it was closed |
|---|---|
| Final v152 SI synchronization | The v152.6 SI (`SUPPLEMENTARY_INFORMATION_v152_6_REVIEWER_STRENGTHENED_FINAL.md`) is recorded in `SUBMISSION_SNAPSHOT.json` with its SHA-256; the correspondence audit found three content gaps (liver-only series, recovery grid, independent evaluation) and closed them by merging into Supplementary Notes 11 and 14 without adding numbering |
| Verified `selcorr` 0.1.1 artifact | The 0.1.1 source was built and checked in a fresh R 4.5.1 environment during rc2 assembly: `R CMD check --no-manual` **Status: OK**; artifact and log in `selcorr_release/` |

## Portability defect fixed

`selcorr/tests/testthat/test-frozen-regression.R` carried a hard-coded default archive path
(a local absolute path). It now reads `SELCORR_ARCHIVE_ROOT` from the environment and skips the
archive-dependent tests when unset. No other file under `R/`, `man/`, `NAMESPACE` or
`vignettes/` changed; no statistical definition changed.

## Still open (before the final tag)

1. Bump the release identifier to `reproducibility-v1.0.0` in README, CITATION.cff and the snapshot.
2. Add the final GitHub release URL to README/CITATION.
3. Regenerate `SHA256_MANIFEST.txt` and the release ZIP from a fresh git checkout.

## Post-rc2 external audit fixes (2026-09-14, pre-push)

Independent pre-push audit found and fixed five packaging defects; no code, data or
release-decision changes:

1. `.gitignore` `*.gz` would have silently dropped `selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz`
   on `git add .` - a fresh clone would then fail `verify_release.py`. Added an explicit
   whitelist exception for the verified tarball.
2. `selcorr/build/vignette.rds` was silently dropped by `*.rds`; `selcorr/build/` is now
   excluded explicitly as a build artifact.
3. README referenced `GITHUB_RELEASE_AUDIT_reproducibility-v1.0.0_RC1.md` (no such file);
   now references the actual RC1/RC2 audit files.
4. README repository map listed a top-level `figures/` directory that does not exist;
   the row now points to the real figure locations and the Figure 4 source-data file.
5. `RELEASE_METADATA.json` / `CITATION.cff` still carried rc1 identifiers; aligned to
   `reproducibility-v1.0.0-rc2`. The final-tag bump, final release URL and
   fresh-checkout manifest/zip remain open author steps (see UPLOAD_CHECKLIST.md).
