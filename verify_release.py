from pathlib import Path
import csv, hashlib, json, sys
root=Path(__file__).resolve().parent
errors=[]
# v152 tables
rec=list(csv.DictReader(open(root/'v152_extensions/source_data/masld_ah_recovery_summary.csv',encoding='utf-8')))
if len(rec)!=16: errors.append('recovery grid must have 16 rows')
if any(float(r['rejection_rate'])!=0 for r in rec): errors.append('frozen recovery rejection must be zero in all cells')
strong=[r for r in rec if float(r['effect_size_sd'])==1.0 and float(r['fraction'])==0.5]
if len(strong)!=1 or abs(float(strong[0]['median_observed_r'])-0.7430)>1e-9 or abs(float(strong[0]['reference_mean'])-0.6687)>1e-9: errors.append('strongest recovery cell mismatch')
rep=list(csv.DictReader(open(root/'v152_extensions/source_data/masld_ah_replay_summary.csv',encoding='utf-8')))
if len(rep)!=5: errors.append('MASLD-AH replay summary must have 5 workflows')
# Package state
# Core set ships the verified source tarball instead of the unpacked package tree.
import tarfile as _tar
_tp = root/'selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz'
if not _tp.exists(): errors.append('verified selcorr tarball missing')
else:
    with _tar.open(_tp, 'r:gz') as _tf:
        _names = _tf.getnames()
        _desc = next((m for m in _names if m.endswith('/DESCRIPTION')), None)
        if _desc is None: errors.append('tarball has no DESCRIPTION')
        else:
            _txt = _tf.extractfile(_desc).read().decode('utf-8', errors='ignore')
            if 'Version: 0.1.1' not in _txt: errors.append('tarball DESCRIPTION not 0.1.1')
        if not any('inst/extdata/masld_ah_recovery_summary.csv' in m for m in _names):
            errors.append('tarball missing recovery fixture')
        if not any('test-recovery-regression.R' in m for m in _names):
            errors.append('tarball missing recovery regression test')
# RC honesty
snap=json.loads((root/'SUBMISSION_SNAPSHOT.json').read_text())
if snap.get('release_status')!='RC_NOT_FINAL': errors.append('RC must not claim final status')
if 'SYNCED_WITH_V152' not in snap.get('supplementary_information_status',''): errors.append('SI synchronization status missing')
if not str(snap.get('supplementary_information_filename','')).startswith('SUPPLEMENTARY_INFORMATION_v152_6'):
    errors.append('snapshot must record the v152.6 SI file name')
if not (root/'selcorr_release/RECONSTRUCTION_STATUS.md').exists(): errors.append('package reconstruction status missing')
import glob as _glob
if not _glob.glob(str(root/'selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz')): errors.append('verified selcorr artifact missing')
if not _glob.glob(str(root/'selcorr_release/R_CHECK_LOG_v0.1.1_RC2.txt')): errors.append('R CMD check log missing')
# Manifest
manifest=root/'SHA256_MANIFEST.txt'
if manifest.exists():
    for line in manifest.read_text().splitlines():
        if not line.strip(): continue
        h, rel=line.split('  ',1)
        p=root/rel
        if not p.exists(): errors.append('manifest missing file: '+rel); continue
        hh=hashlib.sha256(p.read_bytes()).hexdigest()
        if hh!=h: errors.append('hash mismatch: '+rel)
if errors:
    print('FAIL')
    for e in errors: print('-',e)
    sys.exit(1)
print('PASS: reproducibility-v1.0.0-rc2 integrity and release-boundary checks verified')
