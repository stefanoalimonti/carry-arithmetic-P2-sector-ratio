#!/usr/bin/env python3
"""
P2_02: Three sector-ratio objects and their numerical values.

Distinguishes and computes three quantities that are all called
"sector ratios" in related work, but are distinct objects:

  (a) R = S10 / S00  (total-sum ratio, Paper P2's conjecture -> -pi)
      Uses unweighted cascade valuations aggregated over the sector.
      One-step Richardson: 2*R(21) - R(20) gives 3.5 sig. digits.

  (b) <val>_10 / <val>_00  (per-sector mean ratio -> -pi / omega ~ -9.82)
      Equal to R / omega where omega = N10/N00 -> 0.31993.

  (c) mu_00^chi4(s) / mu_10^chi4(s)  (chi4-weighted, s-dependent, ~ -0.13)
      The chi4-weighted carry amplitude ratio on the critical line;
      developed in Paper L §7.4.

Relation: R(K) = omega(K) * (<val>_10 / <val>_00).

References: Paper P2 §7.3 (remark), §10 open problems 4 and 7.
            Paper L §7.4 for object (c).
"""

import sys
import os
import math
import numpy as np

PI = math.pi

# ── Embedded exact enumeration data (E44/E45) ─────────────────────────────
# R(K) = S10(K) / S00(K) from exhaustive enumeration.
R_EXACT = {
    19: -3.10951, 20: -3.12486, 21: -3.13340,
}
# omega(K) = N10(K) / N00(K) from enumeration; limit 0.31993 (Theorem 7).
OMEGA_EXACT = {
    19: 0.31990, 20: 0.31991, 21: 0.31992,
}
OMEGA_LIM = (2 * math.log(4/3) - 0.5) / (2 * math.log(9/8))


def main():
    print("=" * 70)
    print("P2_02: THREE SECTOR-RATIO OBJECTS")
    print("=" * 70)

    # ── Object (a): R(K) = S10/S00 -> -pi ─────────────────────────────────
    print("\n--- Object (a): R(K) = S10/S00  (Paper P2 conjecture) ---")
    print(f"\n  {'K':>4}  {'R(K)':>12}  {'gap |R+pi|':>14}  {'sig. digits':>12}")
    for K in sorted(R_EXACT):
        g = abs(R_EXACT[K] + PI)
        sd = -math.log10(g / PI) if g > 0 else float('inf')
        print(f"  {K:>4d}  {R_EXACT[K]:>12.6f}  {g:>14.6e}  {sd:>12.1f}")

    R20, R21 = R_EXACT[20], R_EXACT[21]
    R_rich = 2 * R21 - R20
    gap_rich = abs(R_rich + PI)
    sd_rich  = -math.log10(gap_rich / PI)
    print(f"\n  One-step Richardson  2*R(21)-R(20) = {R_rich:.8f}")
    print(f"  -pi                               = {-PI:.8f}")
    print(f"  Gap  = {gap_rich:.2e}  ({sd_rich:.1f} significant digits)")

    # ── Object (b): per-sector mean ratio ─────────────────────────────────
    print("\n--- Object (b): <val>_10 / <val>_00 = R / omega ---")
    print(f"\n  Limiting omega = N10/N00 = {OMEGA_LIM:.8f}")
    print(f"  -pi / omega               = {-PI / OMEGA_LIM:.6f}")
    print(f"\n  {'K':>4}  {'R(K)':>12}  {'omega(K)':>10}  {'R/omega':>12}  {'pi/omega':>12}")
    for K in sorted(R_EXACT):
        om = OMEGA_EXACT[K]
        per = R_EXACT[K] / om
        print(f"  {K:>4d}  {R_EXACT[K]:>12.6f}  {om:>10.6f}  {per:>+12.4f}  {-PI/om:>+12.4f}")

    per_rich = R_rich / OMEGA_LIM
    print(f"\n  Richardson per-sector mean: R_rich/omega = {per_rich:.4f}")
    print(f"  Limit -pi/omega            = {-PI/OMEGA_LIM:.4f}")
    print(f"  (Object (b) carries no information beyond R; it is determined by R and omega.)")

    # ── Object (c): chi4-weighted on the critical line ────────────────────
    print("\n--- Object (c): mu_00^chi4(s) / mu_10^chi4(s)  (s-dependent) ---")

    paper_L_exp = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "carry-arithmetic-L-dirichlet-bridge", "experiments"
    )
    sys.path.insert(0, os.path.abspath(paper_L_exp))

    try:
        from _shared import build_highk_bank, primitive_characters, L_hurwitz, n_of_tau

        def chi4(n):
            r = n % 4
            if r % 2 == 0: return 0
            return 1 if r == 1 else -1

        def mu_sector_chi(u_sector, s, tau0=3):
            import cmath
            total = 0j
            for tau, u in u_sector.items():
                if tau < tau0: continue
                n = n_of_tau(tau, tau0=tau0)
                c = chi4(n)
                if c == 0: continue
                total += u * c * cmath.exp(-s * math.log(n))
            return total

        bank = build_highk_bank()
        rec  = bank[999]
        u00, u10 = rec["u00"], rec["u10"]

        print(f"\n  Computing mu_00^chi4(s) / mu_10^chi4(s) on critical line ...")
        print(f"  {'sigma':>6}  {'t':>6}  {'R_chi(s)':>14}  {'|R_chi|':>10}")

        sample_pts = [(0.5, 2.0), (0.5, 5.0), (0.5, 14.1), (0.5, 25.0),
                      (1.5, 0.0), (2.0, 0.0), (5.0, 0.0)]
        R_chi_vals = []
        for sigma, t in sample_pts:
            s = complex(sigma, t)
            m00 = mu_sector_chi(u00, s)
            m10 = mu_sector_chi(u10, s)
            if abs(m10) > 1e-20:
                R_s = m00 / m10
                R_chi_vals.append(R_s.real)
                print(f"  {sigma:>6.1f}  {t:>6.1f}  {R_s.real:>+14.6f}  {abs(R_s):>10.6f}")

        ts_crit = np.arange(1.0, 55.0, 0.5)
        R_crit = []
        for t in ts_crit:
            m00 = mu_sector_chi(u00, complex(0.5, t))
            m10 = mu_sector_chi(u10, complex(0.5, t))
            if abs(m10) > 1e-20:
                R_crit.append((m00 / m10).real)

        if R_crit:
            print(f"\n  Statistics on critical line (t in [1, 55)):")
            print(f"    Mean Re(R_chi) = {np.mean(R_crit):.4f}")
            print(f"    Std  Re(R_chi) = {np.std(R_crit):.4f}")
            print(f"  (Not -pi; this object is s-dependent and developed in Paper L §7.4.)")

    except ImportError as e:
        print(f"  (skipping object (c): cannot import carry bank — {e})")

    print("\n--- Summary ---")
    print(f"  (a) R = S10/S00: Richardson gives {R_rich:.5f}  ({sd_rich:.1f} sig. digits from -pi)")
    print(f"  (b) <val>_10/<val>_00 = R/omega -> {-PI/OMEGA_LIM:.4f}  (= -pi/omega, not -pi)")
    print(f"  (c) mu_00^chi4(s)/mu_10^chi4(s) ~ -0.13 on critical line  (s-dependent)")


if __name__ == "__main__":
    main()
