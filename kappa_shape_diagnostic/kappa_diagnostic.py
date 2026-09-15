# Kappa shape diagnostic (post-freeze, reviewer-triggered)
# Question: is the replay reference inflated because residual permutation makes the
# test-arm effect distribution near-Gaussian while observed effects are heavier-tailed?
# Under pure same-sign selection with independent arms: E[r] = kappa_X * kappa_Y,
# kappa_Z = E|Z| / sqrt(E|Z|^2) = 1/sqrt(1+CV^2). Gaussian kappa = sqrt(2/pi) = 0.7979.
import pandas as pd, numpy as np, os, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "reanalysis_20260722" / "outputs"
RUNTIME_ROOT = Path(os.environ.get('MASLD_OSA_RUNTIME_OUT', ROOT / 'runtime_regenerated'))
OUT = RUNTIME_ROOT / 'kappa_shape_diagnostic'
os.makedirs(OUT, exist_ok=True)

def kappa(mean_abs, sd_abs):
    return mean_abs / np.sqrt(sd_abs**2 + mean_abs**2)

results = {}

# ---------- 1. Primary OSA-MASLD replay (per-draw summaries archived) ----------
pd_draws = pd.read_csv(os.path.join(BASE, "phenotype_null_geometry", "replay_geometry_per_draw.csv"))
pd_obs   = pd.read_csv(os.path.join(BASE, "phenotype_null_geometry", "marginal_geometry_comparison.csv")).set_index("statistic")["observed"]

ka_d = kappa(pd_draws["absMASLD_mean"], pd_draws["absMASLD_sd"])
kt_d = kappa(pd_draws["absOSA_mean"],  pd_draws["absOSA_sd"])
q_d  = pd_draws["r"] / (ka_d * kt_d)

ka_o = kappa(pd_obs["absMASLD_mean"], pd_obs["absMASLD_sd"])
kt_o = kappa(pd_obs["absOSA_mean"],  pd_obs["absOSA_sd"])
r_o  = pd_obs["r"]
q_o  = r_o / (ka_o * kt_o)

B = len(pd_draws)
p_r   = (1 + (pd_draws["r"] >= r_o).sum()) / (1 + B)
p_q   = (1 + (q_d >= q_o).sum()) / (1 + B)
results["primary_OSA_MASLD"] = dict(
    B=B, r_obs=r_o, replay_r_mean=pd_draws["r"].mean(), P_r=p_r,
    kappa_anchor_obs=ka_o, kappa_test_obs=kt_o,
    kappa_test_replay_mean=kt_d.mean(), kappa_anchor_replay_mean=ka_d.mean(),
    corrected_reference=ka_o*kt_o, q_obs=q_o,
    replay_q_mean=q_d.mean(), replay_q_median=q_d.median(),
    replay_q_sd=q_d.std(), P_shape_corrected=p_q,
    q_percentile=(q_d < q_o).mean())

# ---------- 2. Forward stress test MASLD-MASLD (per-gene per-draw archived) ----------
def stress_kappas(gene_files, iter_files):
    frames = [pd.read_csv(f) for f in gene_files]
    g = pd.concat(frames, ignore_index=True)
    def k(x): return np.mean(np.abs(x)) / np.sqrt(np.mean(np.asarray(x)**2))
    per = g.groupby("Iteration").agg(ka=("A_log2FC", k), kt=("B_log2FC", k), n=("Gene","size"))
    it  = pd.concat([pd.read_csv(f) for f in iter_files], ignore_index=True).set_index("Iteration")
    per = per.join(it[["SelectedPearson"]])
    per["q"] = per["SelectedPearson"] / (per["ka"] * per["kt"])
    return per

st_dir = os.path.join(BASE, "positive_control_masld_masld")
st = stress_kappas(
    [os.path.join(st_dir, f) for f in ["replay_genes_0001_0010.csv","replay_genes_0011_1000.csv","replay_genes_1001_2000.csv"]],
    [os.path.join(st_dir, f) for f in ["replay_iterations_0001_0010.csv","replay_iterations_0011_1000.csv","replay_iterations_1001_2000.csv"]])
st_obs_genes = pd.read_csv(os.path.join(st_dir, "observed_selected_set.csv"))
def kser(x): x=np.asarray(x); return np.mean(np.abs(x))/np.sqrt(np.mean(x**2))
st_ka_o = kser(st_obs_genes["A_log2FC"]); st_kt_o = kser(st_obs_genes["B_log2FC"])
st_r_o  = float(pd.read_csv(os.path.join(st_dir,"observed_summary.csv"))["ObservedSelectedPearson"].iloc[0])
st_q_o  = st_r_o/(st_ka_o*st_kt_o)
B2 = len(st)
results["stress_MASLD_MASLD"] = dict(
    B=B2, r_obs=st_r_o, replay_r_mean=st["SelectedPearson"].mean(),
    P_r=(1+(st["SelectedPearson"]>=st_r_o).sum())/(1+B2),
    kappa_anchor_obs=st_ka_o, kappa_test_obs=st_kt_o,
    kappa_test_replay_mean=st["kt"].mean(), kappa_anchor_replay_mean=st["ka"].mean(),
    corrected_reference=st_ka_o*st_kt_o, q_obs=st_q_o,
    replay_q_mean=st["q"].mean(), replay_q_median=st["q"].median(), replay_q_sd=st["q"].std(),
    P_shape_corrected=(1+(st["q"]>=st_q_o).sum())/(1+B2),
    q_percentile=(st["q"]<st_q_o).mean())

# ---------- 3. Reverse stress test (no per-gene archive; observed-side only) ----------
rv_dir = os.path.join(BASE, "positive_control_masld_masld_reverse")
rv_obs_genes = pd.read_csv(os.path.join(rv_dir, "observed_selected_set.csv"))
cols = list(rv_obs_genes.columns)
rv_ka_o = kser(rv_obs_genes[cols[1]]); rv_kt_o = kser(rv_obs_genes[cols[2]])
rv_sum = pd.read_csv(os.path.join(rv_dir, "reverse_summary.csv"))
rv_r_o = float(rv_sum["ObservedSelectedPearson"].iloc[0]) if "ObservedSelectedPearson" in rv_sum.columns else None
rv_it = pd.concat([pd.read_csv(os.path.join(rv_dir,f)) for f in ["replay_iterations_0001_1000.csv","replay_iterations_1001_2000.csv"]])
results["stress_reverse"] = dict(
    B=len(rv_it), r_obs=rv_r_o, replay_r_mean=rv_it["SelectedPearson"].mean(),
    kappa_anchor_obs=rv_ka_o, kappa_test_obs=rv_kt_o,
    corrected_reference=rv_ka_o*rv_kt_o,
    q_obs=(rv_r_o/(rv_ka_o*rv_kt_o) if rv_r_o else None),
    note="no per-draw gene effects archived; per-draw q distribution unavailable")

# ---------- save ----------
with open(os.path.join(OUT, "kappa_diagnostic_summary.json"), "w") as f:
    json.dump(results, f, indent=2, default=float)
pd.DataFrame({"ka":ka_d,"kt":kt_d,"r":pd_draws["r"],"q":q_d}).to_csv(os.path.join(OUT,"primary_per_draw_kappa.csv"), index=False)
st.to_csv(os.path.join(OUT,"stress_per_draw_kappa.csv"))

for k,v in results.items():
    print("="*70); print(k)
    for kk,vv in v.items(): print(f"  {kk}: {vv}")
