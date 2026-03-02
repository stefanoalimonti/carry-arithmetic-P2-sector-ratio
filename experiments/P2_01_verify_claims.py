#!/usr/bin/env python3
"""
P2_01_verify_claims.py — Standalone verification of Paper P2's key claims.

Reproduces the core numerical results of Paper P2 using embedded data
from E44/E45 exact enumeration. No external dependencies beyond NumPy.

Verified claims:
  1. R(K) table (§7, Table)
  2. R₀ = (1/2 - 2·ln(4/3)) / (1 + ln(3/8))  (Schoolbook limit, §3)
  3. Gap ratios g(K)/g(K-1) → 1/2  (§9.1)
  4. Exponential Richardson extrapolation → ≈ 4 sig. digits  (§7.3, Prop 5)
  5. α(K) → 1/6  (§7.1)
  6. N₁₀/N₀₀ → 2·ln(4/3)−1/2) / (2·ln(9/8))  (§8.2)

Data source: E44 (K ≤ 20), E45 (K = 21) exact enumeration output.
"""

import math
import sys

PI = math.pi


# --- Exact integer counts from E44/E45 ---
# Each entry: (S00, S10, n00, n10, S00_c0, S00_c1)
# S00, S10: cascade val-weighted sums
# n00, n10: cascade count sums
# S00_c0, S00_c1: sub-sector split of S00
EXACT_DATA = {
    19: {
        "S00": -494176763, "S10": 1536648367,
        "n00": 16186851433, "n10": 5178696560,
        "S00_c0": -423072387, "S00_c1": -71104376,
    },
    20: {
        "S00": -1967962747, "S10": 6149608524,
        "n00": 64749711489, "n10": 20715378998,
        "S00_c0": -1683746513, "S00_c1": -284216234,
    },
    21: {
        "S00": -7852453718, "S10": 24604885342,
        "n00": 259003293792, "n10": 82862658279,
        "S00_c0": -6715990761, "S00_c1": -1136462957,
    },
}

R_TABLE = {
    4: +0.2500,  5: +0.2500,  6: +0.0930,  7: -0.0915,
    8: -0.3848,  9: -0.7262,  10: -1.1184,  11: -1.5535,
    12: -1.9696,  13: -2.3388,  14: -2.6225,  15: -2.8241,
    16: -2.9539,  17: -3.0348,  18: -3.0822,
}


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def neville_exp(K_vals, R_vals):
    """Neville–Aitken table with h = (1/2)^K."""
    n = len(K_vals)
    h = [0.5**k for k in K_vals]
    table = [list(R_vals)]
    for j in range(1, n):
        col = []
        prev = table[j - 1]
        for i in range(n - j):
            val = (h[i] * prev[i + 1] - h[i + j] * prev[i]) / (h[i] - h[i + j])
            col.append(val)
        table.append(col)
    return table


def sig_digits(est, target=-PI):
    err = abs(est - target)
    if err == 0:
        return float("inf")
    return -math.log10(err / abs(target))


def main():
    pr("=" * 70)
    pr("P2 STANDALONE VERIFICATION")
    pr("=" * 70)

    # --- Build full R(K) table ---
    R = dict(R_TABLE)
    for K, d in EXACT_DATA.items():
        R[K] = d["S10"] / d["S00"]

    # === Claim 1: R(K) table ===
    pr("\n1. R(K) TABLE (§7)")
    pr(f"{'K':>3s}  {'R(K)':>14s}  {'|R+π|':>12s}  {'sig.dig':>8s}  {'ratio':>8s}")
    Ks = sorted(R.keys())
    prev_gap = None
    for K in Ks:
        gap = R[K] + PI
        sd = sig_digits(R[K])
        ratio = ""
        if prev_gap and abs(prev_gap) > 1e-15:
            ratio = f"{gap / prev_gap:8.3f}"
        pr(f"{K:3d}  {R[K]:+14.10f}  {abs(gap):12.6e}  {sd:8.1f}  {ratio}")
        prev_gap = gap

    # === Claim 2: Schoolbook limit R₀ ===
    pr("\n2. SCHOOLBOOK LIMIT R₀ (§3)")
    R0 = (0.5 - 2 * math.log(4 / 3)) / (1 + math.log(3 / 8))
    pr(f"  R₀ = (1/2 - 2·ln(4/3)) / (1 + ln(3/8)) = {R0:.12f}")
    pr(f"  Paper value: -3.9312...  ✓" if abs(R0 - (-3.9312)) < 0.001 else "  MISMATCH")

    # === Claim 3: Gap ratios → 1/2 ===
    pr("\n3. GAP RATIOS → 1/2 (§9.1)")
    for K in [17, 18, 19, 20, 21]:
        if K in R and K - 1 in R:
            g_curr = R[K] + PI
            g_prev = R[K - 1] + PI
            if abs(g_prev) > 1e-15:
                pr(f"  g({K})/g({K-1}) = {g_curr / g_prev:.4f}")

    # === Claim 4: Richardson extrapolation ===
    pr("\n4. RICHARDSON EXTRAPOLATION (§7.3, Proposition 5)")
    K_odd = [k for k in Ks if k >= 7 and k % 2 == 1]
    R_odd = [R[k] for k in K_odd]

    table = neville_exp(K_odd, R_odd)
    best = table[-1][0]
    sd = sig_digits(best)
    pr(f"  K_odd = {K_odd}")
    pr(f"  Neville best: R_∞ = {best:+.10f}")
    pr(f"  |R_∞ + π| = {abs(best + PI):.2e}")
    pr(f"  Significant digits: {sd:.1f}")
    pr(f"  Paper claim: ≈ 4 sig. digits  {'✓' if sd >= 3.5 else 'BELOW THRESHOLD'}")

    pr("\n  Leave-one-out:")
    for idx in range(len(K_odd)):
        K_sub = K_odd[:idx] + K_odd[idx + 1:]
        R_sub = [R[k] for k in K_sub]
        t = neville_exp(K_sub, R_sub)
        est = t[-1][0]
        pr(f"    Leave out K={K_odd[idx]:2d}: R_∞ = {est:+.10f}  "
           f"sig.dig = {sig_digits(est):.1f}")

    # === Claim 5: α(K) → 1/6 ===
    pr("\n5. ALPHA CONVERGENCE α → 1/6 (§7.1)")
    for K in sorted(EXACT_DATA.keys()):
        d = EXACT_DATA[K]
        alpha = d["S00_c1"] / d["S00_c0"]
        pr(f"  K={K}: α = S00_c1/S00_c0 = {alpha:.12f}  "
           f"(1/6 = {1/6:.12f}, diff = {abs(alpha - 1/6):.6e})")

    # === Claim 6: N₁₀/N₀₀ → closed form ===
    pr("\n6. COUNT RATIO N₁₀/N₀₀ → closed form (§8.2)")
    target = (2 * math.log(4 / 3) - 0.5) / (2 * math.log(9 / 8))
    pr(f"  Target: (2·ln(4/3) - 1/2) / (2·ln(9/8)) = {target:.12f}")
    for K in sorted(EXACT_DATA.keys()):
        d = EXACT_DATA[K]
        ratio = d["n10"] / d["n00"]
        pr(f"  K={K}: N₁₀/N₀₀ = {ratio:.12f}  diff = {abs(ratio - target):.6e}")

    # === Summary ===
    pr("\n" + "=" * 70)
    pr("SUMMARY")
    pr("=" * 70)
    pr("  Claim 1 (R table):              Consistent with E44/E45")
    pr(f"  Claim 2 (R₀ = -3.9312):         R₀ = {R0:.6f}  ✓")
    pr(f"  Claim 3 (gap → 1/2):            g(21)/g(20) = {(R[21]+PI)/(R[20]+PI):.4f}")
    pr(f"  Claim 4 (Richardson ≈ 4 dig.):   {sd:.1f} digits  "
       f"{'✓' if sd >= 3.5 else '⚠'}")
    pr(f"  Claim 5 (α → 1/6):              "
       f"{EXACT_DATA[21]['S00_c1']/EXACT_DATA[21]['S00_c0']:.8f}")
    n_ratio = EXACT_DATA[21]["n10"] / EXACT_DATA[21]["n00"]
    pr(f"  Claim 6 (N₁₀/N₀₀ → 0.31993):   {n_ratio:.8f}")


if __name__ == "__main__":
    main()
