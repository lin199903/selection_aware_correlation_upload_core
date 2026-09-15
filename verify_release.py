from pathlib import Path
import csv, json, re, sys, tarfile

root = Path(__file__).resolve().parent
errors = []
TAG = 'reproducibility-v1.0.0-rc3'
MAIN_SHA = 'c906033a010449d415ae531cf7376e6dfa55694f575487e4fe753fb599238332'
SI_SHA = '4cf28fffc3319dad961a703d99a01dbe8d4fada5ad23a7088dbc75b6d4634e6b'

# Paper-extension invariants
ext = root / 'paper_extensions' / 'source_data'
rec = list(csv.DictReader(open(ext / 'masld_ah_recovery_summary.csv', encoding='utf-8')))
if len(rec) != 16:
    errors.append('recovery grid must have 16 rows')
if any(float(r['rejection_rate']) != 0 for r in rec):
    errors.append('recovery rejection must be zero in all archived cells')
strong = [r for r in rec if float(r['effect_size_sd']) == 1.0 and float(r['fraction']) == 0.5]
if len(strong) != 1 or abs(float(strong[0]['median_observed_r']) - 0.7430) > 1e-9 or abs(float(strong[0]['reference_mean']) - 0.6687) > 1e-9:
    errors.append('strongest recovery cell mismatch')
rep = list(csv.DictReader(open(ext / 'masld_ah_replay_summary.csv', encoding='utf-8')))
if len(rep) != 5:
    errors.append('MASLD-AH replay summary must have 5 workflows')

# F0: corrected independent-evaluation null
rows = list(csv.DictReader(open(ext / 'cross_cohort_independent_evaluation_summary.csv', encoding='utf-8')))
if len(rows) != 1:
    errors.append('independent-evaluation summary must have exactly one row')
else:
    r = rows[0]
    expected = {
        'evaluation_r': 0.2371,
        'reference_mean': 0.0034,
        'reference_sd': 0.1621,
        'reference_q95': 0.2544,
        'percentile': 92.7,
        'p_upper_add_one': 0.0735,
        'z': 1.442,
        'valid_draws': 2000.0,
    }
    for key, value in expected.items():
        try:
            observed = float(r[key])
        except Exception:
            errors.append('independent-evaluation field missing/non-numeric: ' + key)
            continue
        if abs(observed - value) > 1e-6:
            errors.append(f'independent-evaluation {key} mismatch: {observed} != {value}')

# P0-2: executed estimator sensitivity, not a missing-analysis placeholder
rows = list(csv.DictReader(open(ext / 'gse135251_estimator_sensitivity_reported_summary.csv', encoding='utf-8')))
est = {r.get('estimator'): r for r in rows}
if set(est) != {'limma', 'DESeq2'}:
    errors.append('estimator-sensitivity summary must contain limma and DESeq2 rows')
else:
    if int(est['limma']['selected_genes']) != 79 or abs(float(est['limma']['selected_r']) - 0.6507) > 1e-6:
        errors.append('limma estimator-sensitivity summary mismatch')
    if int(est['DESeq2']['selected_genes']) != 78 or abs(float(est['DESeq2']['selected_r']) - 0.6296) > 1e-6:
        errors.append('DESeq2 estimator-sensitivity summary mismatch')
    if '0.8690' not in est['DESeq2'].get('comparison_note', ''):
        errors.append('DESeq2 Jaccard declaration missing')

# Package state
_tp = root / 'selcorr_release' / 'selcorr_0.1.1_VERIFIED_RC2.tar.gz'
if not _tp.exists():
    errors.append('verified selcorr tarball missing')
else:
    with tarfile.open(_tp, 'r:gz') as tf:
        names = tf.getnames()
        desc = next((m for m in names if m.endswith('/DESCRIPTION')), None)
        if desc is None:
            errors.append('tarball has no DESCRIPTION')
        else:
            txt = tf.extractfile(desc).read().decode('utf-8', errors='ignore')
            if 'Version: 0.1.1' not in txt:
                errors.append('tarball DESCRIPTION not 0.1.1')
        if not any('inst/extdata/masld_ah_recovery_summary.csv' in m for m in names):
            errors.append('tarball missing recovery fixture')
        if not any('test-recovery-regression.R' in m for m in names):
            errors.append('tarball missing recovery regression test')

# Release-boundary and manuscript provenance
snap = json.loads((root / 'SUBMISSION_SNAPSHOT.json').read_text(encoding='utf-8'))
meta = json.loads((root / 'RELEASE_METADATA.json').read_text(encoding='utf-8'))
if snap.get('release_status') != 'RC_NOT_FINAL':
    errors.append('RC must not claim final status')
if snap.get('public_release_tag') != TAG or meta.get('release_tag_candidate') != TAG:
    errors.append('release tag mismatch')
if snap.get('internal_manuscript_version') != 'v153' or meta.get('internal_manuscript_version') != 'v153':
    errors.append('current internal manuscript version must be v153')
if snap.get('internal_manuscript_version_role') != 'provenance_only':
    errors.append('manuscript version must be provenance only')
if snap.get('main_manuscript_sha256') != MAIN_SHA:
    errors.append('v153 manuscript hash declaration mismatch')
if snap.get('supplementary_information_sha256') != SI_SHA:
    errors.append('v153 SI hash declaration mismatch')
for key in ('main_manuscript_sha256', 'supplementary_information_sha256'):
    if not re.fullmatch(r'[0-9a-f]{64}', str(snap.get(key, ''))):
        errors.append(key + ' must be lowercase SHA-256')
if (root / 'v152_extensions').exists():
    errors.append('version-coupled v152_extensions path must not exist')
if not (root / 'paper_extensions' / 'provenance' / 'MANUSCRIPT_EXTENSION_PROVENANCE.md').exists():
    errors.append('version-neutral extension provenance missing')

if errors:
    print('FAIL')
    for e in errors:
        print('-', e)
    sys.exit(1)
print('PASS: reproducibility-v1.0.0-rc3 F0/P0-2, package, version-decoupling and v153 provenance checks verified')
