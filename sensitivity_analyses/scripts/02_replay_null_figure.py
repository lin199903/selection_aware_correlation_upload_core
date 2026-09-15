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
OUT = os.path.join(RUNTIME_ROOT, 'sensitivity_analyses', 'figures')
os.makedirs(OUT, exist_ok=True)

rcParams['font.family'] = 'Times New Roman'
rcParams['font.size'] = 9
rcParams['axes.linewidth'] = 0.8
rcParams['pdf.fonttype'] = 42

it = pd.read_csv(os.path.join(ARCH, 'outputs', 'phenotype_null', 'phenotype_null_iterations.csv'))
valid = it[it['Converged'] == True]
rnull = valid['SelectedPearson'].values
snull = valid['SameSignSelectedGenes'].values
obs_r = 0.550719414222321
null_mean = rnull.mean()
null_sd = rnull.std(ddof=1)
q95 = np.quantile(rnull, 0.95)
d = (obs_r - null_mean) / null_sd
ci = (0.3014, 0.7679)

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

ax = axes[0]
ax.hist(rnull, bins=40, color='#8ea4c2', edgecolor='white', linewidth=0.4, density=True)
ax.axvline(null_mean, color='#1f4e79', linestyle='--', linewidth=1.2)
ax.axvline(obs_r, color='#c00000', linewidth=1.6)
ax.axvline(q95, color='#7f7f7f', linestyle=':', linewidth=1.2)
ax.axvspan(ci[0], ci[1], color='#c00000', alpha=0.12)
ax.set_xlabel('Cross-arm Pearson r of selected genes')
ax.set_ylabel('Density')
ax.set_xlim(0.30, 1.00)
ax.text(null_mean + 0.006, ax.get_ylim()[1] * 0.96, 'null mean\n0.668', fontsize=7, color='#1f4e79', va='top')
ax.text(obs_r - 0.006, ax.get_ylim()[1] * 0.60, 'observed\n0.551', fontsize=7, color='#c00000', ha='right', va='top')
ax.text(q95 - 0.006, ax.get_ylim()[1] * 0.30, 'null 95th (one-tailed\ncritical value) 0.791', fontsize=7, color='#404040', ha='right')
ax.set_title('A  Replay null distribution (2,000 iterations)', fontsize=9)

ax = axes[1]
ax.hist(snull, bins=16, color='#8ea4c2', edgecolor='white', linewidth=0.4, density=True)
ax.axvline(52, color='#c00000', linewidth=1.6)
ax.axvline(snull.mean(), color='#1f4e79', linestyle='--', linewidth=1.2)
ax.set_xlabel('Same-sign selected genes')
ax.set_ylabel('Density')
ax.text(52 + 1.5, ax.get_ylim()[1] * 0.90, 'observed\n52', fontsize=7, color='#c00000', va='top')
ax.text(snull.mean() + 1.5, ax.get_ylim()[1] * 0.50, 'null mean 55.0', fontsize=7, color='#1f4e79', va='top')
ax.set_title('B  Selected-set size under null', fontsize=9)

fig.text(0.01, 0.01, 'Shaded band: observed r 95% CI by module block bootstrap (B=10,000, cluster resampling): [0.301, 0.768]; '
                     'null 95th = 0.791 (one-tailed critical value).',
         fontsize=7, color='#404040')

fig.tight_layout(w_pad=2.0)
fig.savefig(os.path.join(OUT, 'Supplementary_Figure_4_replay_null_diagnostics.pdf'))
fig.savefig(os.path.join(OUT, 'Supplementary_Figure_4_replay_null_diagnostics.png'), dpi=300)
fig.savefig(os.path.join(OUT, 'Supplementary_Figure_4_replay_null_diagnostics.tiff'), dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
print('figure saved; d=%.3f' % d)
