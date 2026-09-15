#!/usr/bin/env python3
"""
Analytical verification of gate power thresholds.
Confirms that with OSA meta SE = 0.1226 and BH-12 rank-1 two-sided
threshold z = 2.865, the 50% and 80% power points are ~0.35 and ~0.46 SD.
"""

import numpy as np
from scipy import stats

OSA_META_SE = 0.1226

# BH(12) rank-1 two-sided threshold
alpha = 0.05
n_tests = 12
p_crit_rank1 = alpha / n_tests  # 0.05/12 = 0.004167
z_crit = stats.norm.ppf(1 - p_crit_rank1 / 2)  # two-sided
print(f"BH(12) rank-1 two-sided threshold: p_crit = {p_crit_rank1:.6f}, z = {z_crit:.4f}")

# Power: P(|Z| > z_crit) where Z ~ N(delta/SE, 1)
# For delta > 0: P(Z > z_crit) = 1 - Phi(z_crit - delta/SE)
# P(Z < -z_crit) is negligible
# So power = 1 - Phi(z_crit - delta/SE)

def power_at_delta(delta, se):
    return 1 - stats.norm.cdf(z_crit - delta / se)

# Find delta for target power
for target_power in [0.50, 0.80, 0.90, 0.95]:
    # Analytical: z_crit - delta/SE = Phi^{-1}(1 - target_power)
    # delta/SE = z_crit - Phi^{-1}(1 - target_power)
    z_target = stats.norm.ppf(1 - target_power)
    delta_se = z_crit - z_target
    delta = delta_se * OSA_META_SE

    # Verify
    power = power_at_delta(delta, OSA_META_SE)

    print(f"\n{target_power*100:.0f}% power:")
    print(f"  z_target = Phi^-1(1-{target_power}) = {z_target:.4f}")
    print(f"  delta/SE = z_crit - z_target = {z_crit:.4f} - ({z_target:.4f}) = {delta_se:.4f}")
    print(f"  delta = {delta_se:.4f} * {OSA_META_SE:.4f} = {delta:.4f} SD")
    print(f"  verification: power = {power:.4f}")

# Also compute power for key delta values
print("\n" + "-" * 60)
print("Power at key effect sizes (analytical, C2 only):")
for delta in [0.1, 0.2, 0.234, 0.3, 0.35, 0.4, 0.45, 0.46, 0.5, 0.6, 0.8, 1.0]:
    p = power_at_delta(delta, OSA_META_SE)
    print(f"  delta = {delta:.3f} SD: power(C2 only) = {p:.4f}")

# Sanity check against user's numbers
print("\n" + "-" * 60)
print("Cross-check with user's analytical derivation:")
print(f"  50% power delta = (z_crit - 0) * SE_osa = {z_crit:.3f} * {OSA_META_SE:.4f} = {z_crit * OSA_META_SE:.4f} SD")
z_80 = stats.norm.ppf(0.20)  # Phi(z_80) = 0.20, so z_80 = -0.8416
delta_80_analytical = (z_crit - z_80) * OSA_META_SE
print(f"  80% power: z_crit - Phi^-1(0.20) = {z_crit:.3f} - ({z_80:.4f}) = {z_crit - z_80:.4f}")
print(f"  80% power delta = {z_crit - z_80:.4f} * {OSA_META_SE:.4f} = {delta_80_analytical:.4f} SD")
print(f"\n  Ratio 80%/50% = {delta_80_analytical / (z_crit * OSA_META_SE):.2f}")
print(f"  (User reports simulation ratio = 1.33, analytical = 1.30 — consistent)")

print("\nDone. Analytical verification complete.")
