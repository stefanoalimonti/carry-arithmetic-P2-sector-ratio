# Experiments — Paper P2

## Local Verification Script

| Script | Description | Claims verified |
|--------|-------------|-----------------|
| `P2_01_verify_claims.py` | Standalone verification of P2's key numerical claims | §2, §7, §7.1, §7.3, §9.1, §9.2 |

Run:

```bash
python3 P2_01_verify_claims.py
```

This script uses embedded data from E44/E45 exact enumeration (K = 19–21 integer
counts) plus the published R(K) table (K = 4–18), and verifies:

1. **R(K) table** (§7) — consistent with E44/E45
2. **Schoolbook limit** R₀ = (1/2 − 2·ln(4/3)) / (1 + ln(3/8)) = −3.9312…  (§2)
3. **Gap ratios** g(K)/g(K−1) → 1/2  (§7.3, Proposition 4)
4. **Richardson extrapolation** → ≈ 4 sig. digits  (§7.3, Proposition 5)
5. **α(K) convergence**  (§7.2)
6. **N₁₀/N₀₀ → closed form**  (§9.2, Theorem 7)

## Full Reproduction Pipeline

P2's R(K) data originates from exact enumeration in the companion repository
`carry-arithmetic-E-trace-anomaly`. To reproduce from scratch:

1. **E44/E45** (C programs in `carry-arithmetic-E-trace-anomaly/experiments/`):
   exhaustive sector enumeration for K = 4–21.
2. **P1_07** (in `carry-arithmetic-P1-pi-spectral/experiments/`):
   Richardson extrapolation of R(K).

## Requirements

- Python ≥ 3.8, NumPy (only for full pipeline; `P2_01_verify_claims.py` uses
  only the standard library).
