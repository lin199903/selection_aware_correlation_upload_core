import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(SCRIPT_DIR))  # repository root
# Current public-release analysis tree.
ARCH = os.path.join(BASE, 'reanalysis_20260722')
CANON_OUT = os.path.join(BASE, 'sensitivity_analyses', 'outputs')
RUNTIME_ROOT = os.environ.get('MASLD_OSA_RUNTIME_OUT', os.path.join(BASE, 'runtime_regenerated'))
OUT = os.path.join(RUNTIME_ROOT, 'sensitivity_analyses', 'outputs')
PRIMARY_OUT = os.path.join(RUNTIME_ROOT, 'reanalysis_20260722', 'outputs', 'primary')
os.makedirs(OUT, exist_ok=True)
os.makedirs(PRIMARY_OUT, exist_ok=True)

it = pd.read_csv(os.path.join(ARCH, 'outputs', 'phenotype_null', 'phenotype_null_iterations.csv'))
prev = pd.read_csv(os.path.join(ARCH, 'outputs', 'phenotype_null', 'phenotype_null_summary.csv'))

obs_common = 12958
obs_elig = 110
obs_sig_osa = 0
obs_same = 52
obs_r = 0.550719414222321

valid = it[it['Converged'] == True]
assert len(valid) == 2000, len(valid)
rnull = valid['SelectedPearson'].values

n_eff = len(rnull)
upr = int(np.sum(rnull >= obs_r))
p = (upr + 1) / (n_eff + 1)
mcse = float(prev.iloc[0]['MonteCarloSE'])

row = {
    'Analysis': 'OSA_Freedman_Lane_full_selection_replay',
    'PlannedIterations': 2000,
    'ValidIterations': n_eff,
    'FailedIterations': 0,
    'ObservedCommonGenes': obs_common,
    'ObservedSignificantOSAGenes': obs_sig_osa,
    'ObservedEligiblePoolGenes': obs_elig,
    'ObservedSameSignSelectedGenes': obs_same,
    'ObservedSelectedPearson': obs_r,
    'NullMeanPearson': rnull.mean(),
    'NullSDPearson': rnull.std(ddof=1),
    'NullQ025': np.quantile(rnull, 0.025),
    'NullMedian': np.quantile(rnull, 0.5),
    'NullQ975': np.quantile(rnull, 0.975),
    'UpperTailCount': upr,
    'PUpperAddOne': p,
    'MonteCarloSE': mcse,
}

out = pd.DataFrame([row])
out.to_csv(os.path.join(PRIMARY_OUT, 'rule_replay_summary.csv'), index=False)

prev_row = prev.iloc[0]
checks = []
for k in row:
    same = np.isclose(row[k], prev_row[k]) if isinstance(row[k], float) else str(row[k]) == str(prev_row[k])
    checks.append((k, row[k], prev_row[k], same))
ok = all(c[3] for c in checks)

diag = []
diag.append('rule_replay_summary.csv regenerated from iterations: %d valid' % n_eff)
diag.append('MATCHES phenotype_null_summary.csv: %s' % ok)
for c in checks:
    diag.append('  %-26s new=%-22s prev=%-22s same=%s' % (c[0], c[1], c[2], c[3]))
diag.append('')

# ---- Part 1: null-standardized z and null summary (terminology corrected) ----
null_mean = rnull.mean()
null_sd = rnull.std(ddof=1)
null_med = np.quantile(rnull, 0.5)
null_q95 = np.quantile(rnull, 0.95)
null_q025 = np.quantile(rnull, 0.025)
null_q975 = np.quantile(rnull, 0.975)
z = (obs_r - null_mean) / null_sd
pct_below = 100.0 * np.mean(rnull <= obs_r)

diag.append('NULL-STANDARDIZED z: z = (obs - null_mean) / null_sd = %.4f' % z)
diag.append('OBSERVED r empirical percentile within null: %.1f%% (below null median %.4f)' % (pct_below, null_med))
diag.append('OBSERVED r one-tailed upper-tail P (add-one, 2,000 iterations): %.4f (MCSE %.4f)' % (p, mcse))
diag.append('one-tailed 95%% critical value (q95): %.4f' % null_q95)
diag.append('two-sided 95%% null interval: [%.4f, %.4f]' % (null_q025, null_q975))
diag.append('NULL mean sd: %.4f (%.4f)' % (null_mean, null_sd))
diag.append('')

# ---- Part 2: module block bootstrap CI of the observed correlation ----
pool = pd.read_csv(os.path.join(ARCH, 'outputs', 'primary', 'masld_phenotype_hub_sensitivity_eligible_pool.csv'))
sel = pool[pool['SameSignObserved'].astype(bool)].copy()
assert len(sel) == 52, len(sel)
r_obs_check = np.corrcoef(sel['MASLD_log2FC'], sel['OSA_log2FC'])[0, 1]
assert np.isclose(r_obs_check, obs_r, atol=1e-12), r_obs_check

osa_mod = pd.read_csv(os.path.join(ARCH, 'outputs', 'osa_network', 'OSA_module_assignments.csv'))
masld_mod = pd.read_csv(os.path.join(ARCH, 'outputs', 'network', 'MASLD_module_assignments.csv'))


def block_bootstrap_ci(sel_df, mod_df, b=10000, seed=20260726):
    m = mod_df.set_index('Gene')['Module'].reindex(sel_df['Gene']).fillna('_unassigned')
    df = pd.DataFrame({
        'x': sel_df['MASLD_log2FC'].values,
        'y': sel_df['OSA_log2FC'].values,
        'mod': m.values,
    })
    mods = df['mod'].unique()
    clusters = [df.index[df['mod'] == mm].values for mm in mods]
    rng = np.random.default_rng(seed)
    k = len(clusters)
    rs = np.empty(b)
    sizes = np.empty(b, dtype=int)
    for j in range(b):
        picks = rng.integers(0, k, size=k)
        idx = np.concatenate([clusters[pi] for pi in picks])
        x = df['x'].iloc[idx].values
        y = df['y'].iloc[idx].values
        rs[j] = np.corrcoef(x, y)[0, 1]
        sizes[j] = len(idx)
    ci = np.quantile(rs, [0.025, 0.975])
    return ci, rs, sizes


ci_osa, rs_osa, sizes_osa = block_bootstrap_ci(sel, osa_mod)
ci_masld, rs_masld, sizes_masld = block_bootstrap_ci(sel, masld_mod)

cov_osa = sel['Gene'].isin(osa_mod['Gene']).mean()
cov_masld = sel['Gene'].isin(masld_mod['Gene']).mean()

z_f = np.arctanh(obs_r)
se_z = 1.0 / np.sqrt(52 - 3)
fisher_ci = (np.tanh(z_f - 1.96 * se_z), np.tanh(z_f + 1.96 * se_z))

diag.append('OBSERVED r 95%% CI (module block bootstrap, OSA modules, B=10,000, cluster resampling): [%.4f, %.4f]' % (ci_osa[0], ci_osa[1]))
diag.append('  OSA module coverage of the 52 selected genes: %.1f%%; clusters=%d; bootstrap resample n: median %d, range %d-%d' % (
    cov_osa * 100, len(osa_mod.set_index('Gene')['Module'].reindex(sel['Gene']).dropna().unique()), np.median(sizes_osa), sizes_osa.min(), sizes_osa.max()))
diag.append('OBSERVED r 95%% CI (module block bootstrap, MASLD modules, sensitivity): [%.4f, %.4f]' % (ci_masld[0], ci_masld[1]))
diag.append('  MASLD module coverage of the 52 selected genes: %.1f%% (unassigned genes form singleton blocks)' % (cov_masld * 100))
diag.append('OBSERVED r 95%% CI (Fisher, n=52, independence assumption): [%.4f, %.4f]' % fisher_ci)
diag.append('')

# ---- Part 3: metric sensitivity (observed side) and fixed-set permutation null ----
g126 = pd.read_csv(os.path.join(ARCH, 'outputs', 'de', 'GSE126848_DESeq2_primary.csv'))
meta_m = pd.read_csv(os.path.join(ARCH, 'outputs', 'de', 'MASLD_fixed_effect_meta_secondary.csv'))
g135 = pd.read_csv(os.path.join(ARCH, 'outputs', 'osa_de', 'GSE135917_limma_primary.csv'))
g387 = pd.read_csv(os.path.join(ARCH, 'outputs', 'external_osa', 'GSE38792_case_control_limma.csv'))
g750 = pd.read_csv(os.path.join(ARCH, 'outputs', 'external_osa', 'GSE75097_case_control_limma.csv'))

pool2 = pool.merge(g126[['Gene', 'lfcSE']], on='Gene', how='left')
pool2 = pool2.merge(meta_m[['Gene', 'meta_SE']], on='Gene', how='left')
pool2['se_m_meta'] = pool2['meta_SE'].fillna(pool2['meta_SE'].median())
pool2['se_m_same'] = pool2['lfcSE'].fillna(pool2['lfcSE'].median())

g135_t = g135.set_index('Gene')['t']
pool2['se_o_same'] = (pool2['OSA_log2FC'] / pool2['Gene'].map(g135_t)).abs().replace([np.inf], np.nan)
pool2['se_o_same'] = pool2['se_o_same'].fillna(pool2['se_o_same'].median())


def per_gene_se(g, fc_col):
    out = {}
    for _, r in g.iterrows():
        gene = r['Gene']
        if fc_col in r and r['t'] != 0 and not np.isnan(r['t']):
            out[gene] = abs(r[fc_col] / r['t'])
    return out


se135 = per_gene_se(g135, 'logFC')
se387 = per_gene_se(g387, 'logFC')
se750 = per_gene_se(g750, 'logFC')


def meta_se(se_dicts, gene):
    s = [d[gene] for d in se_dicts if gene in d and d[gene] > 0]
    if not s:
        return np.nan
    return np.sqrt(1.0 / np.sum(1.0 / np.asarray(s) ** 2))


pool2['se_o_3'] = [meta_se([se135, se387, se750], g) for g in pool2['Gene']]
pool2['se_o_3'] = pool2['se_o_3'].fillna(pool2['se_o_3'].median())

sel2 = pool2[pool2['SameSignObserved'].astype(bool)]
x = sel2['MASLD_log2FC'].values
y = sel2['OSA_log2FC'].values
xs = sel2['MASLD_log2FC'].values / sel2['se_m_same'].values
ys = sel2['OSA_log2FC'].values / sel2['se_o_same'].values
r_pear = np.corrcoef(x, y)[0, 1]


def spearman(a, b):
    return np.corrcoef(pd.Series(a).rank().values, pd.Series(b).rank().values)[0, 1]


r_spear = spearman(x, y)
r_std = np.corrcoef(xs, ys)[0, 1]

rng = np.random.default_rng(20260727)
B = 10000
cnt_p = cnt_s = cnt_t = 0
for j in range(B):
    yp = rng.permutation(y)
    cnt_p += (np.corrcoef(x, yp)[0, 1] >= r_pear)
    cnt_s += (spearman(x, yp) >= r_spear)
    cnt_t += (np.corrcoef(xs, rng.permutation(ys))[0, 1] >= r_std)

diag.append('METRIC SENSITIVITY (observed selected set, n=52):')
diag.append('  Pearson (log2FC):        %.4f' % r_pear)
diag.append('  Spearman (rank):         %.4f' % r_spear)
diag.append('  standardized-effect:     %.4f  (log2FC / same-source SE)' % r_std)
diag.append('FIXED-SET PAIRWISE-PERMUTATION NULL (B=10,000, pairing broken):')
diag.append('  Pearson P(add-one): %.4f; Spearman P(add-one): %.4f; standardized P(add-one): %.4f' % (
    (cnt_p + 1) / (B + 1), (cnt_s + 1) / (B + 1), (cnt_t + 1) / (B + 1)))
diag.append('  (naive fixed-set significance is robust to the metric; the selection-adjusted null above is decisive)')
diag.append('')

# ---- Part 4: replayed null under alternative metrics (v113_00, if available) ----
ms_path = os.path.join(CANON_OUT, 'null_metric_sensitivity_iterations.csv')
if os.path.exists(ms_path):
    ms = pd.read_csv(ms_path)
    ms_v = ms[ms['Converged'] == True]
    arch_p = valid['SelectedPearson'].values
    arch_p_500 = arch_p[:len(ms_v)]
    n_match = int(np.sum(np.abs(ms_v['PearsonR'].dropna().values - arch_p_500[:len(ms_v['PearsonR'].dropna())]) < 1e-12))
    diag.append('NULL METRIC SENSITIVITY (500-iteration Freedman-Lane replay, seed prefix 50260723):')
    diag.append('  per-iteration Pearson agreement with archived first-500: %d/500' % n_match)
    diag.append('  (not bit-exact: archived run session state/RNGkind/BLAS not archived; validity at distribution level)')
    diag.append('  distribution comparison, Pearson r under null:')
    diag.append('    replay : mean=%.4f sd=%.4f median=%.4f q95=%.4f' % (
        ms_v['PearsonR'].mean(), ms_v['PearsonR'].std(ddof=1), ms_v['PearsonR'].median(), np.quantile(ms_v['PearsonR'], 0.95)))
    diag.append('    archived(all %d): mean=%.4f sd=%.4f median=%.4f q95=%.4f' % (
        len(arch_p), arch_p.mean(), arch_p.std(ddof=1), np.median(arch_p), np.quantile(arch_p, 0.95)))
    obs_map = {'PearsonR': r_pear, 'SpearmanR': r_spear, 'StdR': r_std}
    for col in ['PearsonR', 'SpearmanR', 'StdR']:
        v = ms_v[col].dropna().values
        if len(v) == 0:
            diag.append('  %-12s no finite values' % col)
            continue
        p_up = (int(np.sum(v >= obs_map[col])) + 1) / (len(v) + 1)
        diag.append('  %-12s mean=%.4f sd=%.4f median=%.4f q95=%.4f P_upper(add-one)=%.4f' % (
            col, v.mean(), v.std(ddof=1), np.median(v), np.quantile(v, 0.95), p_up))
    diag.append('')
else:
    diag.append('NULL METRIC SENSITIVITY: not available yet (v113_00 output missing)')
    diag.append('')

# ---- Part 5: selected-set size, pool statistics (P1-5 wording support) ----
diag.append('SELECTED-SET SIZE under null: mean %.2f, sd %.2f, range %d-%d (observed 52)' % (
    valid['SameSignSelectedGenes'].mean(), valid['SameSignSelectedGenes'].std(ddof=1),
    valid['SameSignSelectedGenes'].min(), valid['SameSignSelectedGenes'].max()))
diag.append('SELECTED PEARSON vs SELECTED SIZE correlation under null: %.3f' % np.corrcoef(
    valid['SelectedPearson'].values, valid['SameSignSelectedGenes'].values)[0, 1])

diag.append('')
diag.append('ELIGIBLE POOL (n=%d, MASLD-direction selection): MASLD log2FC sd=%.4f; OSA log2FC sd=%.4f; ratio=%.2f' % (
    len(pool2), pool2['MASLD_log2FC'].std(ddof=1), pool2['OSA_log2FC'].std(ddof=1),
    pool2['OSA_log2FC'].std(ddof=1) / pool2['MASLD_log2FC'].std(ddof=1)))
diag.append('ELIGIBLE POOL per-gene noise: MASLD meta SE median=%.4f; OSA 3-cohort meta SE median=%.4f' % (
    np.median(pool2['se_m_meta']), np.median(pool2['se_o_3'])))
diag.append('ELIGIBLE POOL same-sign fraction: %.3f (%d/110)' % (
    pool2['SameSignObserved'].mean(), int(pool2['SameSignObserved'].sum())))

with open(os.path.join(OUT, 'replay_null_diagnostics.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(diag))
print('\n'.join(diag))
