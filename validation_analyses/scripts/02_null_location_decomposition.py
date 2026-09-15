#!/usr/bin/env python3
"""Diagnostic decomposition of the OSA–MASLD replay-null location.

This is a fixed-margin, gene-label sign-selection diagnostic. It asks how much
selected correlation is expected after restricting the gene universe while
preserving the observed marginal effect distributions. It is not a causal or
additive decomposition of the full Freedman–Lane workflow replay.
"""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNTIME_ROOT = Path(os.environ.get("MASLD_OSA_RUNTIME_OUT", ROOT / "runtime_regenerated"))
OUT = RUNTIME_ROOT / "validation_analyses" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260913
Y_GENES = {
    "UTY", "USP9Y", "DDX3Y", "KDM5D", "EIF1AY", "ZFY", "SRY",
    "NLGN4Y", "RPS4Y1", "RPS4Y2", "TSPY1", "RBMY1A1", "DAZ1",
    "PRKY", "AMELY", "TBL1Y", "PCDH11Y", "TMSB4Y", "VCY",
    "CDY1", "CDY2A", "HSFY1", "TXLNGY", "BPY2", "PRY",
}

masld = pd.read_csv(ROOT / "reanalysis_20260722/outputs/de/GSE126848_DESeq2_primary.csv")
osa = pd.read_csv(ROOT / "reanalysis_20260722/outputs/osa_de/GSE135917_limma_primary.csv")
masld_hubs = set(pd.read_csv(ROOT / "reanalysis_20260722/outputs/network/MASLD_topology_hubs_primary.csv")["Gene"].astype(str))
osa_hubs = set(pd.read_csv(ROOT / "reanalysis_20260722/outputs/osa_network/OSA_topology_hubs_primary.csv")["Gene"].astype(str))
full_replay = pd.read_csv(ROOT / "reanalysis_20260722/outputs/phenotype_null/phenotype_null_summary.csv").iloc[0]

m = masld[["Gene", "log2FoldChange", "padj"]].copy()
m.columns = ["Gene", "MASLD", "MASLD_FDR"]
o = osa[["Gene", "logFC", "adj.P.Val"]].copy()
o.columns = ["Gene", "OSA", "OSA_FDR"]
common = m.merge(o, on="Gene", how="inner")
common = common[np.isfinite(common.MASLD) & np.isfinite(common.OSA)]
common = common[~common.Gene.isin(Y_GENES)].drop_duplicates("Gene").copy()
common["MASLD_DE"] = (common.MASLD_FDR < 0.05) & (common.MASLD.abs() > 0.3)
common["OSA_DE"] = (common.OSA_FDR < 0.05) & (common.OSA.abs() > 0.3)
common["OSA_HUB"] = common.Gene.isin(osa_hubs)
common["MASLD_HUB"] = common.Gene.isin(masld_hubs)
common["ELIGIBLE"] = (common.OSA_HUB & common.MASLD_DE) | (common.MASLD_HUB & common.OSA_DE)
common["SAME_SIGN"] = common.MASLD * common.OSA > 0


def corr(x, y):
    if len(x) < 3 or np.std(x, ddof=1) == 0 or np.std(y, ddof=1) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def perm_sign_reference(df, n_perm, seed):
    x = df.MASLD.to_numpy(float)
    y = df.OSA.to_numpy(float)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm, float)
    ns = np.empty(n_perm, int)
    for i in range(n_perm):
        yp = rng.permutation(y)
        keep = x * yp > 0
        ns[i] = int(keep.sum())
        vals[i] = corr(x[keep], yp[keep]) if ns[i] >= 3 else np.nan
    v = vals[np.isfinite(vals)]
    n = ns[np.isfinite(vals)]
    return {
        "permutations": int(n_perm),
        "valid_permutations": int(len(v)),
        "sign_selected_reference_mean": float(np.mean(v)),
        "sign_selected_reference_sd": float(np.std(v, ddof=1)),
        "sign_selected_reference_median": float(np.median(v)),
        "sign_selected_reference_q025": float(np.quantile(v, 0.025)),
        "sign_selected_reference_q975": float(np.quantile(v, 0.975)),
        "mean_selected_n": float(np.mean(n)),
    }

sets = [
    ("all_common", common, 3000),
    ("OSA_hubs_only", common[common.OSA_HUB], 50000),
    ("MASLD_DE_only", common[common.MASLD_DE], 5000),
    ("joint_eligible_pool", common[common.ELIGIBLE], 100000),
]
rows = []
for j, (name, df, nperm) in enumerate(sets):
    same = df[df.MASLD * df.OSA > 0]
    ref = perm_sign_reference(df, nperm, SEED + j * 1009)
    rows.append({
        "stage": name,
        "genes": int(len(df)),
        "observed_paired_correlation": corr(df.MASLD.to_numpy(), df.OSA.to_numpy()),
        "observed_same_sign_n": int(len(same)),
        "observed_same_sign_correlation": corr(same.MASLD.to_numpy(), same.OSA.to_numpy()),
        **ref,
    })

eligible = common[common.ELIGIBLE]
unselected = eligible[~eligible.SAME_SIGN]
rows.append({
    "stage": "eligible_pool_unselected_complement_observed_only",
    "genes": int(len(unselected)),
    "observed_paired_correlation": corr(unselected.MASLD.to_numpy(), unselected.OSA.to_numpy()),
    "observed_same_sign_n": 0,
    "observed_same_sign_correlation": np.nan,
    "permutations": 0,
    "valid_permutations": 0,
    "sign_selected_reference_mean": np.nan,
    "sign_selected_reference_sd": np.nan,
    "sign_selected_reference_median": np.nan,
    "sign_selected_reference_q025": np.nan,
    "sign_selected_reference_q975": np.nan,
    "mean_selected_n": np.nan,
})

dfout = pd.DataFrame(rows)
dfout.to_csv(OUT / "null_location_decomposition.csv", index=False)

joint_mean = float(dfout.loc[dfout.stage == "joint_eligible_pool", "sign_selected_reference_mean"].iloc[0])
summary = {
    "analysis": "OSA_MASLD_null_location_diagnostic_decomposition",
    "seed": SEED,
    "theoretical_independent_gaussian_sign_only": float(2 / np.pi),
    "full_freedman_lane_replay_mean": float(full_replay["NullMeanPearson"]),
    "joint_eligible_fixed_margin_sign_only_mean": joint_mean,
    "full_minus_joint_fixed_margin": float(full_replay["NullMeanPearson"] - joint_mean),
    "interpretation": (
        "Diagnostic, not causal/additive: most of the elevation of the observed-workflow replay mean "
        "above 2/pi is already present after empirical candidate-pool construction and fixed-margin "
        "same-sign selection; the remaining difference to the full sample-level replay is small."
    ),
}
(OUT / "null_location_decomposition_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(dfout.to_string(index=False))
print(json.dumps(summary, indent=2))
