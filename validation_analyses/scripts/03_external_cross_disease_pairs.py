#!/usr/bin/env python3
"""Independent-effect cross-disease stress tests using fixed topology priors.

MASLD effects come from GSE130970 (independent of primary GSE126848), OSA effects
from external GSE38792 or GSE75097 (independent of primary GSE135917). The
candidate pool applies the frozen primary topology priors and cohort-specific DE
rules. Because the external OSA cohorts contribute no genes at the frozen DE
threshold, the pool is anchored by primary OSA topology hubs intersected with
GSE130970 MASLD DE genes. Gene-label permutation of external OSA effects within
that fixed pool regenerates same-sign selection while conditioning on the pool.
This is an external-pair stress test, not a full sample-level workflow replay.
"""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNTIME_ROOT = Path(os.environ.get("MASLD_OSA_RUNTIME_OUT", ROOT / "runtime_regenerated"))
OUT = RUNTIME_ROOT / "validation_analyses" / "outputs" / "external_cross_disease"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260914
B = 100000
Y_GENES = {"UTY","USP9Y","DDX3Y","KDM5D","EIF1AY","ZFY","SRY","NLGN4Y","RPS4Y1","RPS4Y2","TSPY1","RBMY1A1","DAZ1","PRKY","AMELY","TBL1Y","PCDH11Y","TMSB4Y","VCY","CDY1","CDY2A","HSFY1","TXLNGY","BPY2","PRY"}

masld = pd.read_csv(ROOT / "reanalysis_20260722/outputs/de/GSE130970_DESeq2_replication.csv")
osa_hubs = set(pd.read_csv(ROOT / "reanalysis_20260722/outputs/osa_network/OSA_topology_hubs_primary.csv")["Gene"].astype(str))
masld_hubs = set(pd.read_csv(ROOT / "reanalysis_20260722/outputs/network/MASLD_topology_hubs_primary.csv")["Gene"].astype(str))

m = masld[["Gene","log2FoldChange","padj"]].copy()
m.columns = ["Gene","MASLD","MASLD_FDR"]
m = m.drop_duplicates("Gene")
m["MASLD_DE"] = (m.MASLD_FDR < 0.05) & (m.MASLD.abs() > 0.3)


def corr(x,y):
    if len(x)<3 or np.std(x,ddof=1)==0 or np.std(y,ddof=1)==0: return np.nan
    return float(np.corrcoef(x,y)[0,1])

def evaluate(path, accession, seed):
    o = pd.read_csv(path)[["Gene","logFC","adj.P.Val"]].copy()
    o.columns = ["Gene","OSA","OSA_FDR"]
    o = o.drop_duplicates("Gene")
    o["OSA_DE"] = (o.OSA_FDR < .05) & (o.OSA.abs() > .3)
    d = m.merge(o,on="Gene",how="inner")
    d = d[np.isfinite(d.MASLD)&np.isfinite(d.OSA)&~d.Gene.isin(Y_GENES)].copy()
    d["OSA_HUB_MASLD_DE"] = d.Gene.isin(osa_hubs) & d.MASLD_DE
    d["MASLD_HUB_OSA_DE"] = d.Gene.isin(masld_hubs) & d.OSA_DE
    pool = d[d.OSA_HUB_MASLD_DE | d.MASLD_HUB_OSA_DE].copy()
    keep = pool.MASLD.to_numpy()*pool.OSA.to_numpy()>0
    obs = corr(pool.MASLD.to_numpy()[keep], pool.OSA.to_numpy()[keep])
    rng=np.random.default_rng(seed)
    x=pool.MASLD.to_numpy(float); y=pool.OSA.to_numpy(float)
    null=np.empty(B,float); ns=np.empty(B,int)
    for i in range(B):
        yp=rng.permutation(y); k=x*yp>0; ns[i]=k.sum(); null[i]=corr(x[k],yp[k]) if ns[i]>=3 else np.nan
    v=null[np.isfinite(null)]
    p=(1+np.sum(v>=obs))/(1+len(v)) if np.isfinite(obs) else np.nan
    out={
      "MASLD_accession":"GSE130970", "OSA_accession":accession,
      "common_genes":int(len(d)), "MASLD_DE_genes":int(d.MASLD_DE.sum()), "OSA_DE_genes":int(d.OSA_DE.sum()),
      "eligible_pool":int(len(pool)), "same_sign_selected":int(keep.sum()), "observed_selected_r":obs,
      "reference_mean":float(np.mean(v)), "reference_sd":float(np.std(v,ddof=1)),
      "reference_q025":float(np.quantile(v,.025)), "reference_q975":float(np.quantile(v,.975)),
      "upper_tail_add_one_p":float(p), "permutations":B,
      "mean_selected_n":float(np.mean(ns[np.isfinite(null)])),
      "pool_route_OSA_hub_x_MASLD_DE":int(pool.OSA_HUB_MASLD_DE.sum()),
      "pool_route_MASLD_hub_x_OSA_DE":int(pool.MASLD_HUB_OSA_DE.sum()),
      "interpretation_role":"external independent-effect cross-disease fixed-pool sign-replay stress test; not full sample-level workflow replay"
    }
    pool.to_csv(OUT/f"{accession}_eligible_pool.csv",index=False)
    return out

rows=[]
for j,acc in enumerate(["GSE38792","GSE75097"]):
    path=ROOT/f"reanalysis_20260722/outputs/external_osa/{acc}_case_control_limma.csv"
    rows.append(evaluate(path,acc,SEED+j*1009))
pd.DataFrame(rows).to_csv(OUT/"external_cross_disease_summary.csv",index=False)
(OUT/"external_cross_disease_summary.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
print(pd.DataFrame(rows).to_string(index=False))
