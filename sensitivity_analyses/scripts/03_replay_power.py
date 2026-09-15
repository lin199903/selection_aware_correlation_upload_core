import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(SCRIPT_DIR))  # repository root
# Current public-release analysis tree.
ARCH = os.path.join(BASE, 'reanalysis_20260722')
RUNTIME_ROOT = os.environ.get('MASLD_OSA_RUNTIME_OUT', os.path.join(BASE, 'runtime_regenerated'))
OUT = os.path.join(RUNTIME_ROOT, 'sensitivity_analyses', 'outputs')
FIG = os.path.join(RUNTIME_ROOT, 'sensitivity_analyses', 'figures')
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

rcParams['font.family'] = 'Times New Roman'
rcParams['font.size'] = 9
rcParams['axes.linewidth'] = 0.8
rcParams['pdf.fonttype'] = 42

pool = pd.read_csv(os.path.join(ARCH, 'outputs', 'primary', 'masld_phenotype_hub_sensitivity_eligible_pool.csv'))
meta_m = pd.read_csv(os.path.join(ARCH, 'outputs', 'de', 'MASLD_fixed_effect_meta_secondary.csv'))
g135 = pd.read_csv(os.path.join(ARCH, 'outputs', 'osa_de', 'GSE135917_limma_primary.csv'))
g387 = pd.read_csv(os.path.join(ARCH, 'outputs', 'external_osa', 'GSE38792_case_control_limma.csv'))
g750 = pd.read_csv(os.path.join(ARCH, 'outputs', 'external_osa', 'GSE75097_case_control_limma.csv'))

df = pool[['Gene', 'MASLD_log2FC', 'OSA_log2FC', 'SameSignObserved']].merge(
    meta_m[['Gene', 'meta_log2FC', 'meta_SE']], on='Gene', how='left')
df['se_m'] = df['meta_SE'].fillna(df['meta_SE'].median())

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

df['se_o'] = [meta_se([se135, se387, se750], g) for g in df['Gene']]
df['se_o'] = df['se_o'].fillna(df['se_o'].median())

t_m = df['MASLD_log2FC'].values
t_o = df['OSA_log2FC'].values
se_m = df['se_m'].values
se_o = df['se_o'].values
s = np.sign(t_m)
s[s == 0] = 1

it = pd.read_csv(os.path.join(ARCH, 'outputs', 'phenotype_null', 'phenotype_null_iterations.csv'))
rnull = it[it['Converged'] == True]['SelectedPearson'].values
null_mean = float(np.mean(rnull))
null_q95 = float(np.quantile(rnull, 0.95))
null_sd = float(np.std(rnull, ddof=1))

def replay_r(tm, to):
    sel = np.sign(tm) == np.sign(to)
    n = int(sel.sum())
    if n < 3:
        return np.nan, n
    r = np.corrcoef(tm[sel], to[sel])[0, 1]
    return r, n

deltas = np.arange(0.0, 0.85, 0.05)
det_rows = []
for d in deltas:
    tm = t_m + d * s
    to = t_o + d * s
    r, n = replay_r(tm, to)
    det_rows.append((d, r, n))
det = pd.DataFrame(det_rows, columns=['delta', 'det_r', 'det_n'])

rng = np.random.default_rng(20260724)
N = 2000
sim_rows = []
for d in deltas:
    rs = []
    ns = []
    for _ in range(N):
        tm = t_m + d * s + rng.normal(0, se_m)
        to = t_o + d * s + rng.normal(0, se_o)
        r, n = replay_r(tm, to)
        rs.append(r)
        ns.append(n)
    rs = np.asarray(rs)
    ns = np.asarray(ns)
    valid = ~np.isnan(rs)
    rs = rs[valid]
    ns = ns[valid]
    power = float(np.mean(rs > null_q95))
    se_p = np.sqrt(power * (1 - power) / len(rs))
    sim_rows.append({
        'delta': d,
        'sim_mean_r': rs.mean(),
        'sim_sd_r': rs.std(ddof=1),
        'sim_median_r': np.median(rs),
        'sim_q025': np.quantile(rs, 0.025),
        'sim_q975': np.quantile(rs, 0.975),
        'power_vs_null_q95': power,
        'power_se': se_p,
        'sim_mean_n': ns.mean(),
        'det_r': det.loc[det['delta'] == d, 'det_r'].iloc[0],
        'det_n': det.loc[det['delta'] == d, 'det_n'].iloc[0],
    })
sim = pd.DataFrame(sim_rows)
sim.to_csv(os.path.join(OUT, 'replay_power_sensitivity.csv'), index=False)

pvals = sim['power_vs_null_q95'].values
above = int(np.argmax(pvals >= 0.80))
d1 = sim['delta'].iloc[above]
p1 = pvals[above]
d0 = sim['delta'].iloc[above - 1]
p0 = pvals[above - 1]
d80_interp = d0 + (0.80 - p0) / (p1 - p0) * (d1 - d0)

d0row = sim.loc[sim['delta'] == 0.0].iloc[0]

lines = ['REPLAY POWER SENSITIVITY: observed-effect-conditioned (per-gene shared effect, log2FC units)']
lines.append('this simulation adds a shared effect and independent estimation noise to the OBSERVED MASLD/OSA effect vectors;')
lines.append('it is conditional on the observed effect architecture and is NOT the phenotype-replay null')
lines.append('null reference (phenotype-replay null): mean=%.4f; q95=%.4f (from 2,000 real replay iterations)' % (null_mean, null_q95))
lines.append('per-gene noise calibration: MASLD meta SE median=%.4f; OSA 3-cohort meta SE median=%.4f' % (
    np.median(se_m), np.median(se_o)))
lines.append('')
lines.append('%-6s %-8s %-8s %-10s %-10s %-8s %-6s' % ('delta', 'det_r', 'sim_md', 'sim_q025', 'sim_q975', 'power', 'se'))
for _, r in sim.iterrows():
    lines.append('%-6.2f %-8.3f %-8.3f %-10.3f %-10.3f %-8.3f %-6.4f' % (
        r['delta'], r['det_r'], r['sim_median_r'], r['sim_q025'], r['sim_q975'], r['power_vs_null_q95'], r['power_se']))

lines.append('')
lines.append('delta=0 status: the observed-effect-conditioned simulation does NOT reproduce the phenotype-replay null;')
lines.append('  median r = %.3f and %.3f of iterations exceed the replay-null q95 = %.3f' % (
    d0row['sim_median_r'], d0row['power_vs_null_q95'], null_q95))
lines.append('  (under a true null, the exceedance probability at the null critical value would be ~0.05)')
lines.append('  the position of the observed r is therefore evaluated against the real replay null, not this simulation')
lines.append('')
lines.append('replay delta80 (per-gene): grid = %.2f (first grid point with power >= 0.80, conservative);' % d1)
lines.append('  linear interpolation between power %.3f at delta=%.2f and power %.3f at delta=%.2f gives delta80 ~ %.3f' % (
    p0, d0, p1, d1, d80_interp))
lines.append('  formula: delta80 = %.2f + (0.80 - %.3f)/(%.3f - %.3f) * %.2f' % (d0, p0, p1, p0, d1 - d0))
lines.append('')
lines.append('observed r = 0.551; replay-null mean = %.3f; replay-null q95 = %.3f; observed r below null mean by z = %.2f' % (
    null_mean, null_q95, (0.5507194 - null_mean) / null_sd))
lines.append('scale note: the pathway-level gate threshold (0.455 SD) is on a different scale and is not compared here')

with open(os.path.join(OUT, 'replay_power_sensitivity.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n'.join(lines))

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
ax = axes[0]
ax.plot(sim['delta'], sim['power_vs_null_q95'], color='#c00000', marker='o', ms=3, lw=1.4)
ax.axhline(0.80, color='#404040', linestyle=':', lw=1.0)
ax.axhline(0.05, color='#7f7f7f', linestyle='--', lw=1.0)
ax.annotate('power = 0.80', xy=(0.42, 0.83), fontsize=7, color='#404040')
ax.annotate('alpha = 0.05', xy=(0.42, 0.07), fontsize=7, color='#7f7f7f')
ax.set_xlabel('Per-gene shared effect \u03b4 (log2FC)')
ax.set_ylabel('Replay power (P < 0.05)')
ax.set_ylim(-0.03, 1.03)
ax.set_title('A  Replay power, observed-effect-conditioned', fontsize=9)

ax = axes[1]
ax.plot(sim['delta'], sim['sim_median_r'], color='#1f4e79', marker='o', ms=3, lw=1.4,
        label='median r (observed-effect-conditioned)')
ax.fill_between(sim['delta'], sim['sim_q025'], sim['sim_q975'], color='#1f4e79', alpha=0.15)
ax.axhline(null_q95, color='#404040', linestyle=':', lw=1.0)
ax.axhline(null_mean, color='#404040', linestyle=':', lw=1.0)
ax.axhline(0.551, color='#c00000', linestyle='--', lw=1.0)
ax.text(0.02, null_q95 + 0.012, 'replay-null q95 = 0.791', fontsize=7, color='#404040')
ax.text(0.02, null_mean - 0.028, 'replay-null mean = 0.668', fontsize=7, color='#404040')
ax.text(0.02, 0.551 + 0.012, 'observed r = 0.551', fontsize=7, color='#c00000')
ax.text(0.48, 0.36, 'delta = 0: median 0.585;\n0.001 of iterations exceed\nreplay-null q95 (not the null)', fontsize=6.5, color='#1f4e79')
ax.set_xlabel('Per-gene shared effect \u03b4 (log2FC)')
ax.set_ylabel('Cross-arm Pearson r')
ax.set_title('B  Expected r vs \u03b4 (observed-effect-conditioned)', fontsize=9)
ax.legend(loc='lower right', fontsize=6.5, frameon=False)

fig.tight_layout(w_pad=2.0)
fig.savefig(os.path.join(FIG, 'Supplementary_Figure_5_replay_power.pdf'))
fig.savefig(os.path.join(FIG, 'Supplementary_Figure_5_replay_power.png'), dpi=300)
fig.savefig(os.path.join(FIG, 'Supplementary_Figure_5_replay_power.tiff'), dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
print('figures saved')
