# carry-arithmetic-P2-sector-ratio

**The Sector Ratio in Binary Multiplication: From Markov Failure to Transcendence**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

Computational companion paper providing enumeration data and structural decomposition backing Paper P1. Establishes that the Markov model fails to predict the sector ratio at $K \geq 7$, and traces the failure to bit-sharing correlations at the D-parity boundary. The cascade valuation $R^{\text{cas}}(K)$ isolates the stopping-time contribution whose limit is conjectured to be $-\pi$ (Conjecture 1 of [P1]).

## Status

Computational companion paper. Numerical claims rely on E44/E45 exact enumeration data.

## Repository Structure

```
paper/sector_ratio.md                The paper
experiments/
  P2_01_verify_claims.py             Standalone verification of key claims (§3, §7, §8, §9)
  P2_02_sector_ratio_types.py        Per-sector mean ratio and weight separation (§10, OP 7)
  README.md                          Experiment index and reproduction pipeline
```

## Reproduction

```bash
python experiments/P2_01_verify_claims.py    # Verifies R(K) table, R₀, Richardson, α→1/6
```

`P2_01_verify_claims.py` uses embedded E44/E45 data and requires only the Python standard library.

Full reproduction of R(K) for K = 4–21 requires the C enumerator in `carry-arithmetic-E-trace-anomaly` (see `experiments/README.md`).

## Dependencies

- Python ≥ 3.8 (standard library only for P2_01)
- NumPy (for full pipeline only)

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [E] | The Trace Anomaly of Binary Multiplication (**experiments**) | [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly) |
| [P1] | Pi from Pure Arithmetic | [`carry-arithmetic-P1-pi-spectral`](https://github.com/stefanoalimonti/carry-arithmetic-P1-pi-spectral) |
| [G] | The Angular Uniqueness of Base 2 | [`carry-arithmetic-G-angular-uniqueness`](https://github.com/stefanoalimonti/carry-arithmetic-G-angular-uniqueness) |
| [L] | The Carry–Dirichlet Bridge | [`carry-arithmetic-L-dirichlet-bridge`](https://github.com/stefanoalimonti/carry-arithmetic-L-dirichlet-bridge) |

### Citation

```bibtex
@article{alimonti2026sector_ratio,
  author  = {Alimonti, Stefano},
  title   = {The Sector Ratio in Binary Multiplication: From Markov Failure to Transcendence},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-P2-sector-ratio}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
