from pathlib import Path
import csv, json, sys, tarfile
root=Path(__file__).resolve().parent
errors=[]

# Paper-extension tables
ext=root/'paper_extensions/source_data'
rec=list(csv.DictReader(open(ext/'masld_ah_recovery_summary.csv',encoding='utf-8')))
if len(rec)!=16: errors.append('recovery grid must have 16 rows')
if any(float(r['rejection_rate'])!=0 for r in rec): errors.append('recovery rejection must be zero in all archived cells')
strong=[r for r in rec if float(r['effect_size_sd'])==1.0 and float(r['fraction'])==0.5]
if len(strong)!=1 or abs(float(strong[0]['median_observed_r'])-0.7430)>1e-9 or abs(float(strong[0]['reference_mean'])-0.6687)>1e-9: errors.append('strongest recovery cell mismatch')
rep=list(csv.DictReader(open(ext/'masld_ah_replay_summary.csv',encoding='utf-8')))
if len(rep)!=5: errors.append('MASLD-AH replay summary must have 5 workflows')

# Package state
_tp=root/'selcorr_release/selcorr_0.1.1_VERIFIED_RC2.tar.gz'
if not _tp.exists(): errors.append('verified selcorr tarball missing')
else:
    with tarfile.open(_tp,'r:gz') as _tf:
        names=_tf.getnames()
        desc=next((m for m in names if m.endswith('/DESCRIPTION')),None)
        if desc is None: errors.append('tarball has no DESCRIPTION')
        else:
            txt=_tf.extractfile(desc).read().decode('utf-8',errors='ignore')
            if 'Version: 0.1.1' not in txt: errors.append('tarball DESCRIPTION not 0.1.1')
        if not any('inst/extdata/masld_ah_recovery_summary.csv' in m for m in names): errors.append('tarball missing recovery fixture')
        if not any('test-recovery-regression.R' in m for m in names): errors.append('tarball missing recovery regression test')

# Release-boundary metadata
snap=json.loads((root/'SUBMISSION_SNAPSHOT.json').read_text())
meta=json.loads((root/'RELEASE_METADATA.json').read_text())
if snap.get('release_status')!='RC_NOT_FINAL': errors.append('RC must not claim final status')
if snap.get('public_release_tag')!='reproducibility-v1.0.0-rc3': errors.append('snapshot public tag mismatch')
if meta.get('release_tag_candidate')!='reproducibility-v1.0.0-rc3': errors.append('metadata public tag mismatch')
if snap.get('internal_manuscript_version_role')!='provenance_only': errors.append('manuscript version must be provenance only')
if meta.get('internal_manuscript_version_role')!='provenance only; not a public release identifier': errors.append('metadata manuscript-role contract mismatch')
if (root/'v152_extensions').exists(): errors.append('version-coupled v152_extensions path must not exist')
if not (root/'paper_extensions/source_data').exists(): errors.append('paper_extensions/source_data missing')
if not (root/'paper_extensions/provenance/MANUSCRIPT_EXTENSION_PROVENANCE.md').exists(): errors.append('version-neutral extension provenance missing')

if errors:
    print('FAIL')
    for e in errors: print('-',e)
    sys.exit(1)
print('PASS: reproducibility-v1.0.0-rc3 integrity and version-decoupling checks verified')
