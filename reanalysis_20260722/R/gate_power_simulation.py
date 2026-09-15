#!/usr/bin/env python3
"""
Gate Power Simulation: six-condition shared-programme gate power analysis
for the OSA–MASLD worked application in the current manuscript.

RECALIBRATED v2 — 2026-07-24
Uses observed cohort-level SEs from cohort_program_effects.csv and empirical
tau2 from the meta-analysis rather than naive sqrt(1/n1 + 1/n2).

Key calibration targets:
  - MASLD four-cohort meta SE (mTORC1) = 0.0490
  - OSA three-cohort meta SE (mTORC1) = 0.1226
  - MASLD tau2 = 0.0013, OSA tau2 = 0.0000

The simulation asserts that meta-analytic SEs match these targets within +/-10%
before proceeding to the power curve.

Author: Junyi Lin
Date: 2026-07-24
"""

import numpy as np
from scipy import stats
import sys

np.random.seed(20260724)

# =============================================================================
# 1. OBSERVED COHORT-LEVEL SEs (from cohort_program_effects.csv, mTORC1)
#    AND EMPIRICAL TAU2 FROM META-ANALYSIS
# =============================================================================

MASLD_COHORT_SES = np.array([0.0686, 0.1408, 0.0934, 0.0927])   # GSE126848, GSE130970, GSE89632, GSE135251
OSA_COHORT_SES = np.array([0.1965, 0.1997, 0.2532])              # GSE135917, GSE38792, GSE75097
ORTHOGONAL_SES = np.array([0.0949, 0.0647, 0.7276])              # AHI, CPAP meta, IH (from CSV)

TAU2_MASLD = 0.0013     # from program_meta_analysis.csv, MASLD_all mTORC1
TAU2_OSA   = 0.0000     # from program_meta_analysis.csv, OSA_all_disease mTORC1

TARGET_MASLD_SE = 0.0490   # from program_meta_analysis.csv, MASLD_all mTORC1
TARGET_OSA_SE   = 0.1226   # from program_meta_analysis.csv, OSA_all_disease mTORC1

N_PATHWAYS = 12
SE_TOLERANCE = 0.10       # +/-10% calibration tolerance


# =============================================================================
# 2. RANDOM-EFFECTS META-ANALYSIS
# =============================================================================

def random_effects_meta(effects, ses, fixed_tau2=None):
    """Random-effects meta-analysis using DerSimonian-Laird estimator.
    
    Parameters
    ----------
    effects : array-like
        Cohort-level effect estimates
    ses : array-like
        Cohort-level standard errors
    fixed_tau2 : float or None
        If provided, use this tau2 instead of estimating from data.
        This ensures calibration: the simulation uses empirical tau2
        (TAU2_MASLD=0.0013 from program_meta_analysis.csv,
         TAU2_OSA=0.0000) rather than re-estimating from each replicate.
    """
    k = len(effects)
    if k == 0:
        return {"effect": np.nan, "se": np.nan, "p_value": np.nan, "tau2": np.nan,
                "I2": np.nan, "weights": None, "q_stat": np.nan}

    if fixed_tau2 is not None:
        tau2 = fixed_tau2
    else:
        w_fixed = 1.0 / ses**2
        mu_fixed = np.sum(w_fixed * effects) / np.sum(w_fixed)
        if k > 1:
            q_stat = np.sum(w_fixed * (effects - mu_fixed)**2)
            c = np.sum(w_fixed) - np.sum(w_fixed**2) / np.sum(w_fixed)
            tau2 = max(0.0, (q_stat - (k - 1)) / c) if c > 0 else 0.0
        else:
            tau2 = 0.0

    w_random = 1.0 / (ses**2 + tau2)
    mu_random = np.sum(w_random * effects) / np.sum(w_random)
    se_random = np.sqrt(1.0 / np.sum(w_random))

    p_value = 2.0 * stats.norm.sf(abs(mu_random / se_random))

    return {
        "effect": mu_random,
        "se": se_random,
        "p_value": p_value,
        "tau2": tau2,
        "weights": w_random / np.sum(w_random),
    }


# =============================================================================
# 3. BH FDR
# =============================================================================

def bh_fdr(p_values):
    n = len(p_values)
    order = np.argsort(p_values)
    sorted_p = p_values[order]
    adjusted = np.minimum(1.0, sorted_p * n / (np.arange(n) + 1))
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    result = np.zeros(n)
    result[order] = adjusted
    rejected = result < 0.05
    return result, rejected


# =============================================================================
# 4. GATE CONDITION CHECKS (unchanged logic)
# =============================================================================

def check_gate(masld_cohort_effects, masld_cohort_ses,
               osa_cohort_effects, osa_cohort_ses,
               masld_meta, osa_meta,
               masld_pvals_all, osa_pvals_all,
               ortho_effects, ortho_ses,
               target_idx=0):
    c = {}

    masld_fdr, masld_rej = bh_fdr(masld_pvals_all)
    c["c1_masld_fdr"] = masld_rej[target_idx]

    osa_fdr, osa_rej = bh_fdr(osa_pvals_all)
    same_dir = np.sign(masld_meta["effect"]) == np.sign(osa_meta["effect"])
    c["c2_osa_fdr_same_dir"] = osa_rej[target_idx] and same_dir

    masld_core_signs = np.sign(masld_cohort_effects[:2])
    osa_signs = np.sign(osa_cohort_effects)
    shared_sign = np.sign(masld_meta["effect"])
    c["c3_sign_concordance"] = (
        np.sum(masld_core_signs == shared_sign) >= 2 and
        np.sum(osa_signs == shared_sign) >= 2
    )

    ortho_support = False
    for i in range(len(ortho_effects)):
        z = ortho_effects[i] / ortho_ses[i]
        p_one = 1.0 - stats.norm.cdf(z)
        if ortho_effects[i] * shared_sign > 0 and p_one < 0.05:
            ortho_support = True
            break
    c["c4_orthogonal_support"] = ortho_support

    masld_max_w = np.max(masld_meta["weights"])
    osa_max_w = np.max(osa_meta["weights"])
    c["c5_weight_balance"] = not (masld_max_w > 0.80 and osa_max_w > 0.80)

    masld_loo_ok = True
    for i in range(len(masld_cohort_effects)):
        idx = np.ones(len(masld_cohort_effects), dtype=bool)
        idx[i] = False
        loo = random_effects_meta(masld_cohort_effects[idx], masld_cohort_ses[idx])
        if np.sign(loo["effect"]) != shared_sign:
            masld_loo_ok = False
            break
    osa_loo_ok = True
    for i in range(len(osa_cohort_effects)):
        idx = np.ones(len(osa_cohort_effects), dtype=bool)
        idx[i] = False
        loo = random_effects_meta(osa_cohort_effects[idx], osa_cohort_ses[idx])
        if np.sign(loo["effect"]) != shared_sign:
            osa_loo_ok = False
            break
    c["c6_leave_one_out"] = masld_loo_ok and osa_loo_ok

    c["all_pass"] = all([
        c["c1_masld_fdr"], c["c2_osa_fdr_same_dir"],
        c["c3_sign_concordance"], c["c4_orthogonal_support"],
        c["c5_weight_balance"], c["c6_leave_one_out"]
    ])
    return c


# =============================================================================
# 5. SINGLE REPLICATE
# =============================================================================

def simulate_one_replicate(true_shared_effect, tau2_masld=TAU2_MASLD, tau2_osa=TAU2_OSA):
    masld_ses = MASLD_COHORT_SES
    osa_ses   = OSA_COHORT_SES
    ortho_ses = ORTHOGONAL_SES

    n_masld = len(masld_ses)
    n_osa   = len(osa_ses)
    n_ortho = len(ortho_ses)

    masld_cohort_effects = np.zeros((N_PATHWAYS, n_masld))
    osa_cohort_effects   = np.zeros((N_PATHWAYS, n_osa))

    for p in range(N_PATHWAYS):
        delta = true_shared_effect if p == 0 else 0.0
        masld_cohort_effects[p, :] = delta + np.random.normal(0, np.sqrt(tau2_masld), n_masld)
        osa_cohort_effects[p, :]   = delta + np.random.normal(0, np.sqrt(tau2_osa), n_osa)

    masld_observed = masld_cohort_effects + np.random.normal(0, masld_ses, (N_PATHWAYS, n_masld))
    osa_observed   = osa_cohort_effects   + np.random.normal(0, osa_ses, (N_PATHWAYS, n_osa))

    masld_pvals = np.zeros(N_PATHWAYS)
    osa_pvals   = np.zeros(N_PATHWAYS)
    masld_metas = []
    osa_metas   = []

    for p in range(N_PATHWAYS):
        m_meta = random_effects_meta(masld_observed[p, :], masld_ses, fixed_tau2=tau2_masld)
        o_meta = random_effects_meta(osa_observed[p, :], osa_ses, fixed_tau2=tau2_osa)
        masld_pvals[p] = m_meta["p_value"]
        osa_pvals[p]   = o_meta["p_value"]
        masld_metas.append(m_meta)
        osa_metas.append(o_meta)

    ortho_effects = np.random.normal(true_shared_effect, ortho_ses)

    gate = check_gate(
        masld_observed[0, :], masld_ses,
        osa_observed[0, :], osa_ses,
        masld_metas[0], osa_metas[0],
        masld_pvals, osa_pvals,
        ortho_effects, ortho_ses,
        target_idx=0
    )

    null_gate_pass = False
    for p in range(1, N_PATHWAYS):
        null_gate = check_gate(
            masld_observed[p, :], masld_ses,
            osa_observed[p, :], osa_ses,
            masld_metas[p], osa_metas[p],
            masld_pvals, osa_pvals,
            np.random.normal(0, ortho_ses), ortho_ses,
            target_idx=p
        )
        if null_gate["all_pass"]:
            null_gate_pass = True
            break

    return {
        "gate_pass":       gate["all_pass"],
        "c1": gate["c1_masld_fdr"],
        "c2": gate["c2_osa_fdr_same_dir"],
        "c3": gate["c3_sign_concordance"],
        "c4": gate["c4_orthogonal_support"],
        "c5": gate["c5_weight_balance"],
        "c6": gate["c6_leave_one_out"],
        "masld_effect": masld_metas[0]["effect"],
        "masld_se":     masld_metas[0]["se"],
        "masld_p":      masld_pvals[0],
        "osa_effect":   osa_metas[0]["effect"],
        "osa_se":       osa_metas[0]["se"],
        "osa_p":        osa_pvals[0],
        "null_gate_pass": null_gate_pass,
    }


# =============================================================================
# 6. CALIBRATION ASSERTION
# =============================================================================

def assert_calibration(n_cal=5000):
    """
    Verify that null-simulation meta-analytic SEs reproduce the observed
    meta SEs (0.0490 MASLD, 0.1226 OSA) within +/-10%.
    """
    masld_ses_sim = []
    osa_ses_sim = []
    for _ in range(n_cal):
        rep = simulate_one_replicate(0.0)
        masld_ses_sim.append(rep["masld_se"])
        osa_ses_sim.append(rep["osa_se"])

    masld_se_mean = np.mean(masld_ses_sim)
    osa_se_mean   = np.mean(osa_ses_sim)

    masld_rel_err = abs(masld_se_mean - TARGET_MASLD_SE) / TARGET_MASLD_SE
    osa_rel_err   = abs(osa_se_mean - TARGET_OSA_SE) / TARGET_OSA_SE

    ok = True
    if masld_rel_err > SE_TOLERANCE:
        print(f"  FAIL: simulated MASLD meta SE = {masld_se_mean:.4f} vs target {TARGET_MASLD_SE:.4f} "
              f"(relative error {masld_rel_err:.3f} > {SE_TOLERANCE:.2f})")
        ok = False
    else:
        print(f"  PASS: simulated MASLD meta SE = {masld_se_mean:.4f} (target {TARGET_MASLD_SE:.4f}, "
              f"rel err {masld_rel_err:.3f})")

    if osa_rel_err > SE_TOLERANCE:
        print(f"  FAIL: simulated OSA meta SE = {osa_se_mean:.4f} vs target {TARGET_OSA_SE:.4f} "
              f"(relative error {osa_rel_err:.3f} > {SE_TOLERANCE:.2f})")
        ok = False
    else:
        print(f"  PASS: simulated OSA meta SE = {osa_se_mean:.4f} (target {TARGET_OSA_SE:.4f}, "
              f"rel err {osa_rel_err:.3f})")

    if not ok:
        print("\nERROR: Calibration assertion failed. SE mismatch > 10%.")
        print("This indicates a discrepancy between simulated and observed precision.")
        sys.exit(1)

    print(f"\n  Calibration assertion passed: simulated meta SEs replicate observed values "
          f"within {SE_TOLERANCE*100:.0f}% tolerance.\n")


# =============================================================================
# 7. POWER CURVE
# =============================================================================

def compute_power_curve(effect_grid, n_reps=2000, verbose=True):
    results = []
    for delta in effect_grid:
        passes = c1_p = c2_p = c3_p = c4_p = c5_p = c6_p = null_p = 0

        for _ in range(n_reps):
            rep = simulate_one_replicate(delta)
            if rep["gate_pass"]: passes += 1
            if rep["c1"]: c1_p += 1
            if rep["c2"]: c2_p += 1
            if rep["c3"]: c3_p += 1
            if rep["c4"]: c4_p += 1
            if rep["c5"]: c5_p += 1
            if rep["c6"]: c6_p += 1
            if rep["null_gate_pass"]: null_p += 1

        power = passes / n_reps
        se = np.sqrt(power * (1 - power) / n_reps)
        ci_l = max(0, power - 1.96 * se)
        ci_u = min(1, power + 1.96 * se)

        row = {
            "true_effect": delta, "power": power, "se": se,
            "ci_lower": ci_l, "ci_upper": ci_u, "replicates": n_reps,
            "c1_masld_fdr": c1_p / n_reps, "c2_osa_fdr": c2_p / n_reps,
            "c3_sign_concordance": c3_p / n_reps, "c4_orthogonal": c4_p / n_reps,
            "c5_weight": c5_p / n_reps, "c6_loo": c6_p / n_reps,
            "global_fpr": null_p / n_reps,
        }
        results.append(row)

        if verbose:
            print(f"  delta={delta:.2f}: power={power:.3f} [{ci_l:.3f}, {ci_u:.3f}], "
                  f"C1={row['c1_masld_fdr']:.3f}, C2={row['c2_osa_fdr']:.3f}, "
                  f"C3={row['c3_sign_concordance']:.3f}, C4={row['c4_orthogonal']:.3f}")
    return results


def find_threshold(results, target_power=0.80):
    for r in results:
        if r["power"] >= target_power:
            return r["true_effect"], r["power"], r["ci_lower"], r["ci_upper"]
    return None, None, None, None



def self_test():
    """Fast runtime/portability check using the same frozen core functions.

    This does not replace the configured 5,000/2,000-replicate analysis; it verifies
    that the public script imports, calibrates its fixed-SE meta-analysis machinery,
    simulates replicates, evaluates the six-condition gate and serializes a small
    result in the current environment.
    """
    print("GATE POWER SCRIPT SELF-TEST")
    assert_calibration(n_cal=2)
    rep = simulate_one_replicate(0.25)
    required = ("gate_pass", "masld_se", "osa_se", "c1", "c2", "c3", "c4", "c5", "c6")
    missing = [k for k in required if k not in rep]
    if missing:
        raise RuntimeError(f"self-test missing keys: {missing}")
    if not (np.isfinite(rep["masld_se"]) and np.isfinite(rep["osa_se"])):
        raise RuntimeError("self-test produced non-finite meta SE")
    tiny = compute_power_curve([0.0, 0.5], n_reps=3, verbose=False)
    if len(tiny) != 2 or any("power" not in r for r in tiny):
        raise RuntimeError("self-test power-curve failure")
    print("PASS: gate_power_simulation.py self-test")

# =============================================================================
# 8. MAIN
# =============================================================================

def main():
    print("=" * 72)
    print("GATE POWER SIMULATION (RECALIBRATED v2)")
    print("Six-condition shared-programme gate, observed cohort SEs")
    print("OSA-MASLD worked application: six-condition shared-program gate")
    print("=" * 72)

    print("\nObserved cohort-level SEs (from cohort_program_effects.csv, mTORC1):")
    print(f"  MASLD: {MASLD_COHORT_SES}")
    print(f"  OSA:   {OSA_COHORT_SES}")
    print(f"  Orthogonal layers: {ORTHOGONAL_SES}")
    print(f"\n  tau2 MASLD = {TAU2_MASLD}, tau2 OSA = {TAU2_OSA}")
    print(f"  Target meta SE: MASLD = {TARGET_MASLD_SE}, OSA = {TARGET_OSA_SE}")

    # --- Calibration check ---
    print("\n" + "-" * 72)
    print("CALIBRATION ASSERTION: simulated meta SEs must replicate observed values")
    print("-" * 72)
    assert_calibration(n_cal=5000)

    # --- Null calibration ---
    print("-" * 72)
    print("NULL CALIBRATION (delta=0)")
    print("-" * 72)
    n_null = 5000
    null = compute_power_curve([0.0], n_reps=n_null, verbose=False)[0]
    print(f"  delta=0: gate FPR = {null['power']:.4f} [{null['ci_lower']:.4f}, {null['ci_upper']:.4f}]")
    print(f"  Condition-level: C1={null['c1_masld_fdr']:.4f}, C2={null['c2_osa_fdr']:.4f}, "
          f"C3={null['c3_sign_concordance']:.4f}, C4={null['c4_orthogonal']:.4f}")
    print(f"  Family-wise FPR: {null['global_fpr']:.4f}")

    # --- Main power curve ---
    print("\n" + "=" * 72)
    print("MAIN POWER CURVE")
    print("=" * 72)
    effect_grid = np.arange(0.0, 1.6, 0.05)
    n_reps = 2000
    print(f"Delta grid: {effect_grid[::4]}, {n_reps} reps/point\n")

    results = compute_power_curve(effect_grid, n_reps=n_reps)

    # --- Thresholds ---
    print("\n" + "-" * 72)
    print("POWER THRESHOLDS")
    print("-" * 72)
    for target in [0.50, 0.80, 0.90, 0.95]:
        d, p, lo, hi = find_threshold(results, target)
        if d is not None:
            print(f"  {target*100:.0f}% power: delta >= {d:.2f} SD (power={p:.3f} [{lo:.3f}, {hi:.3f}])")

    # --- What this excludes ---
    print("\n" + "-" * 72)
    print("WHAT THE 0/12 RESULT EXCLUDES")
    print("-" * 72)
    d80, p80, lo80, hi80 = find_threshold(results, 0.80)
    if d80 is not None:
        print(f"  The gate has >=80% power to detect shared pathway effects of delta >= {d80:.2f} SD.")
        print(f"  The 0/12 result is therefore compatible with true shared effects < {d80:.2f} SD.")
        print(f"  At delta = {d80:.2f}, power = {p80:.3f} [{lo80:.3f}, {hi80:.3f}].")
        print(f"  The observed OSA mTORC1 estimate (0.236, SE 0.123, 95% CI [-0.004, 0.476])")
        print(f"  lies within the unresolved range, consistent with this power analysis.")

    # --- Condensed table ---
    print("\n" + "-" * 72)
    print("CONDENSED POWER TABLE (for manuscript)")
    print("-" * 72)
    print(f"  {'delta':>6s}  {'Power':>7s}  {'95% CI':>16s}  {'C1_masld':>9s}  {'C2_osa':>8s}")
    print(f"  {'-'*6}  {'-'*7}  {'-'*16}  {'-'*9}  {'-'*8}")
    for r in results:
        if r["true_effect"] in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5]:
            print(f"  {r['true_effect']:6.2f}  {r['power']:7.3f}  [{r['ci_lower']:6.3f}, {r['ci_upper']:6.3f}]  "
                  f"{r['c1_masld_fdr']:9.3f}  {r['c2_osa_fdr']:8.3f}")

    # --- Save results ---
    import json
    from pathlib import Path
    
    null_p = int(round(null["power"] * n_null))
    ci_lower = 0.0
    ci_upper = 1 - 0.05**(1/n_null)
    if hasattr(stats, 'binomtest'):
        bt = stats.binomtest(null_p, n_null)
        ci = bt.proportion_ci(confidence_level=0.95)
        ci_lower = float(ci.low)
        ci_upper = float(ci.high)
    
    null_calibration = {
        "replicates": n_null,
        "gate_passages": null_p,
        "fpr": null_p / n_null,
        "exact_ci_lower": ci_lower,
        "exact_ci_upper": ci_upper,
    }
    
    output = {
        "simulation": "gate_power_v2_recalibrated",
        "calibration": {
            "target_masld_se": TARGET_MASLD_SE,
            "target_osa_se": TARGET_OSA_SE,
            "tolerance": SE_TOLERANCE,
            "cohort_ses_masld": MASLD_COHORT_SES.tolist(),
            "cohort_ses_osa": OSA_COHORT_SES.tolist(),
            "tau2_masld": TAU2_MASLD,
            "tau2_osa": TAU2_OSA,
        },
        "power_curve": results,
        "null_calibration": null_calibration,
    }
    out_path = Path(__file__).resolve().parents[1] / "outputs" / "programs" / "gate_power_simulation_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")
    print("Done.")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        self_test()
    else:
        main()
