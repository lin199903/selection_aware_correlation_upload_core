import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(SCRIPT_DIR))  # repository root
RUNTIME_ROOT = os.environ.get('MASLD_OSA_RUNTIME_OUT', os.path.join(BASE, 'runtime_regenerated'))
OUT = os.path.join(RUNTIME_ROOT, 'validation_analyses', 'outputs')
FIG = os.path.join(RUNTIME_ROOT, 'validation_analyses', 'figures')
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

rcParams['font.family'] = 'Times New Roman'
rcParams['font.size'] = 9
rcParams['axes.linewidth'] = 0.8
rcParams['pdf.fonttype'] = 42

G = 120
qs = np.array([0.10, 0.25, 0.50, 0.75, 0.90])
a = stats.norm.ppf(1 - qs / 2)
lam = stats.norm.pdf(a) / (1 - stats.norm.cdf(a))
rho_theory = lam**2 / (1 + a * lam)
e_m_theory = G * qs**2 / 2

rng = np.random.default_rng(20260801)
N = 4_000_000
x = rng.standard_normal(N)
y = rng.standard_normal(N)
r_same_sign = np.corrcoef(x[x * y > 0], y[x * y > 0])[0, 1]
rho_sim = np.zeros(len(qs))
for i, q in enumerate(qs):
    m = (np.abs(x) > a[i]) & (np.abs(y) > a[i]) & (x * y > 0)
    rho_sim[i] = np.corrcoef(x[m], y[m])[0, 1]

sim = pd.read_csv(os.path.join(BASE, 'reanalysis_20260722', 'outputs', 'simulation',
                               'selection_fraction_sensitivity_summary.csv'))
sub = sim[(sim['rho'] == 0.1) & (sim['confounding'] == 'none') &
          (sim['method'] == 'naive_fixed_set_gene_label')].sort_values('selection_fraction')

out_rows = []
for i, q in enumerate(qs):
    row = sub[sub['selection_fraction'] == q].sort_values('n_group')
    meds = row['median_statistic'].tolist()
    ns = row['mean_selected_n'].tolist()
    out_rows.append({
        'selection_fraction': q,
        'a': a[i],
        'lambda': lam[i],
        'rho_theory': rho_theory[i],
        'E_M_theory': e_m_theory[i],
        'rho_sim_independent_gaussian': rho_sim[i],
        'sim_median_corr_n8': meds[0],
        'sim_median_corr_n16': meds[1],
        'sim_median_corr_n32': meds[2],
        'sim_median_corr_mean': float(np.mean(meds)),
        'sim_mean_n_n8': ns[0],
        'sim_mean_n_n16': ns[1],
        'sim_mean_n_n32': ns[2],
    })
tab = pd.DataFrame(out_rows)
tab.to_csv(os.path.join(OUT, 'analytic_null_predictions.csv'), index=False)

r_same_sign_str = f'{r_same_sign:.4f}'
summary = []
summary.append('ANALYTIC NULL FOR DIRECTIONAL SELECTION (independent Gaussian null)')
summary.append('')
summary.append('1. Same-sign selection only: E[r | XY>0] = 2/pi = 0.6366')
summary.append(f'   independent-Gaussian simulation (N=4e6): {r_same_sign_str}')
summary.append('')
summary.append('2. Magnitude-threshold + same-sign selection (a = Phi^-1(1-q/2)):')
summary.append('   rho_sel(q) = lambda(a)^2 / (1 + a*lambda(a)),  lambda(a) = phi(a)/(1-Phi(a))')
summary.append('   E[M] = G*q^2/2')
summary.append('')
summary.append('q     a         rho_theory  rho_sim(independent)  sim_median_corr (grid, rho=0.1; n8/n16/n32; mean)')
for i, q in enumerate(qs):
    summary.append(f'{q:.2f}  {a[i]:7.4f}  {rho_theory[i]:9.4f}  {rho_sim[i]:19.4f}  '
                   f'{tab["sim_median_corr_n8"].iloc[i]:8.3f}/{tab["sim_median_corr_n16"].iloc[i]:.3f}/'
                   f'{tab["sim_median_corr_n32"].iloc[i]:.3f}  ({tab["sim_median_corr_mean"].iloc[i]:.3f})')
summary.append('')
summary.append('3. Real-data replay null mean 0.668 (archived 2,000-iteration phenotype replay)')
summary.append('   vs ideal same-sign benchmark 2/pi = 0.6366: difference +0.031')
summary.append('   (consistent with gene-gene dependence in real transcriptomes)')
summary.append('')
summary.append('4. Simulation grid caveat: median correlation is set to 0 when the selected set')
summary.append('   contains <4 genes (R implementation rule); at q=0.10 the expected set size is')
summary.append('   0.6 genes (degenerate boundary), at q=0.25 it is 3.75 (partially degenerate).')
with open(os.path.join(OUT, 'analytic_null_sensitivity.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(summary))

qgrid = np.linspace(0.20, 1.00, 400)
ag = stats.norm.ppf(1 - qgrid / 2)
lamg = stats.norm.pdf(ag) / (1 - stats.norm.cdf(ag))
rho_grid = lamg**2 / (1 + ag * lamg)

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

ax = axes[0]
ax.plot(qgrid, rho_grid, color='#1f4e79', linewidth=1.6, label='Analytic null $\\rho_{sel}(q)$')
mk = tab['selection_fraction'] >= 0.5
for col in ['sim_median_corr_n8', 'sim_median_corr_n16', 'sim_median_corr_n32']:
    ax.scatter(tab['selection_fraction'][mk], tab[col][mk], color='#c00000', s=12, zorder=3, alpha=0.55)
ax.scatter(tab['selection_fraction'][mk], tab['sim_median_corr_mean'][mk], color='#c00000', s=32,
           zorder=4, marker='D', label='Simulation median (mean over n=8/16/32)')
ax.axhline(2 / np.pi, color='#7f7f7f', linestyle=':', linewidth=1.1)
ax.text(0.985, 2 / np.pi + 0.010, '2/$\\pi$ = 0.637 (same-sign only)', fontsize=7, color='#404040', ha='right')
ax.axvspan(0.20, 0.50, color='#c00000', alpha=0.06)
ax.text(0.35, 0.97, 'partially\\ndegenerate', fontsize=6.5, color='#c00000', ha='center', va='top')
ax.set_xlabel('Selection fraction q')
ax.set_ylabel('Selected-set correlation')
ax.set_xlim(0.20, 1.00)
ax.set_ylim(0.55, 1.02)
ax.legend(fontsize=6.5, loc='lower left')
ax.set_title('A  Directional selection manufactures correlation', fontsize=9)

ax = axes[1]
ax.plot(qgrid, 120 * qgrid**2 / 2, color='#1f4e79', linewidth=1.6, label='Analytic $Gq^2/2$ (G=120)')
for col in ['sim_mean_n_n8', 'sim_mean_n_n16', 'sim_mean_n_n32']:
    ax.scatter(tab['selection_fraction'], tab[col], color='#c00000', s=12, zorder=3, alpha=0.55)
ax.scatter(tab['selection_fraction'], tab[['sim_mean_n_n8', 'sim_mean_n_n16', 'sim_mean_n_n32']].mean(axis=1),
           color='#c00000', s=32, zorder=4, marker='D', label='Simulation mean')
ax.axvline(0.258, color='#c00000', linestyle='--', linewidth=0.9)
ax.text(0.258, 46, 'set size 4\n(estimability floor)', fontsize=6.5, color='#c00000', ha='center')
ax.set_xlabel('Selection fraction q')
ax.set_ylabel('Expected selected-set size (genes)')
ax.set_xlim(0.20, 1.00)
ax.set_ylim(0, 64)
ax.legend(fontsize=6.5, loc='upper left')
ax.set_title('B  Information retained by the selected set', fontsize=9)

fig.text(0.01, 0.01,
         'Analytic predictions follow classical moments of truncated bivariate normal distributions '
         '(Gupta & Tracy, 1980). Simulation points: median/mean across replicates in the '
         'near-independent null grid (within-block correlation 0.1, no confounding). '
         'Correlations for selected sets with fewer than 4 genes are undefined and were set to 0.',
         fontsize=6.5, color='#404040')

fig.tight_layout(w_pad=2.0)
fig.savefig(os.path.join(FIG, 'Supplementary_Figure_6_analytic_null.pdf'))
fig.savefig(os.path.join(FIG, 'Supplementary_Figure_6_analytic_null.png'), dpi=300)
fig.savefig(os.path.join(FIG, 'Supplementary_Figure_6_analytic_null.tiff'), dpi=300)
print('\n'.join(summary))
print('figure saved')
