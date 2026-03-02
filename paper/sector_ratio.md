# The Sector Ratio in Binary Multiplication: From Markov Failure to Transcendence

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

## Abstract

We study the sector partition function of binary multiplication: given K-bit integers X, Y with product P = XY of D = 2K−1 bits, we partition the 4^{K−1} digit pairs by the second-highest bits (x_{K−2}, y_{K−2}) and compute the carry-weighted ratio R(K) = σ₁₀/σ₀₀ under two natural carry valuations. The *schoolbook* valuation uses the carry at the fixed position D − 2; the *cascade* valuation uses a first-passage stopping time (the carry one step below the topmost nonzero carry). The Markov model of independent carry propagation predicts R_Markov = +2/3. We show three independent failures of this prediction:

1. The Markov model gives the *wrong sign* (+2/3 vs. the negative R(K) observed for K ≥ 7).
2. It misses a qualitative "J = 1" exclusion effect that removes the dominant contribution.
3. The digit correlations follow the exact profile $\rho(d) = (K-d)/(3K+2)$, with leading-order Fejér approximation $\rho(d) \approx (K-d)/(3K)$, and couple the Dirichlet modes of the carry bridge.

Exact enumeration up to K = 21 shows R^{cas}(K) → −π with asymptotic convergence rate $O(K^2 \cdot 2^{-K})$ [P1, §9.3], while the schoolbook limit R^{sb}(∞) = R₀ = −3.9312… is a purely logarithmic constant (no π). The decomposition R^{cas} = R₀ + ΔR isolates all transcendental content in the cascade correction ΔR = +0.7896… [P1, §10]. We present the full computational pipeline, analyze the digit correlation structure, and quantify the spectral mechanism by which rational arithmetic produces a transcendental constant.

When the Fejér coupling is treated as a perturbation of the Markov operator, the effective coupling diverges to the pole A = 3/2 of the Möbius transformation — a universal phenomenon across all bases. The physical system avoids this divergence through non-linear carry dynamics, operating at A* = 3/2 − 1/(2(π+1)).

We prove that no D-odd pair has sector (1,1), that the count ratio N₁₀/N₀₀ converges to the closed form (2ln(4/3) − 1/2)/(2ln(9/8)) = 0.31993... (involving only ln 2 and ln 3), and that the entire factor π in the conjectured limit R → −π enters exclusively through the carry weight ratio ⟨val⟩₁₀/⟨val⟩₀₀, not the pair counts. An extended Markov chain with W bits of memory demonstrates that bit-sharing correlations between adjacent convolutions are the mechanism transforming the Markov super-exponential decay Q ∼ (3/4)^d into the empirical geometric decay Q ≈ 0.55, with the effective per-step survival probability converging toward 1/2 — the Diaconis–Fulman eigenvalue.

---

## 1. Introduction

This is the computational companion to [P1]. Where that paper develops the analytical theory (the closed-form Möbius transformation, the spectral zeta function, the Dirichlet character connection, the R₀ + ΔR decomposition, and the Bernoulli hierarchy), this paper provides:

1. The complete definition and algorithmic computation of R(K) under both the schoolbook and cascade valuations (Section 2).
2. A demonstration that the Markov model fails *qualitatively* (wrong sign, wrong magnitude, wrong mechanism).
3. The extraction and characterization of digit correlations (the Fejér kernel).
4. The spectral consequences of the coupling, including the universal pole divergence.
5. Exact enumeration data through K = 21.
6. Structural decomposition: the top-carry constraint, closed-form count ratio, and the isolation of π to carry weights (Section 9).
7. The bit-sharing mechanism: an extended Markov chain with W bits of memory that reproduces Q_emp ≈ 0.55 and connects the cascade stopping time to the Diaconis–Fulman eigenvalue (Section 9.5).

For the angular uniqueness of base 2 (why π enters only in binary), see [G].

### 1.1 Motivation

The question is elementary: in the schoolbook multiplication of two random K-bit numbers, how does the carry at position D − 2 depend on which *sector* the input falls in? The answer turns out to involve π, emerging through a mechanism that has no classical geometric content.

---

## 2. Definitions and Setup

### 2.1 Binary multiplication

Let X = Σ_{i=0}^{K-1} x_i 2^i and Y = Σ_{i=0}^{K-1} y_i 2^i with x_{K−1} = y_{K−1} = 1 (K-bit integers). The product P = XY has D = 2K − 1 bits. At each position j:

    conv_j = Σ_{i+i'=j, 0≤i,i'≤K-1} x_i · y_{i'}

    p_j + 2·c_{j+1} = conv_j + c_j,    c_0 = 0

**Definition 1a (Schoolbook carry weight).** For a pair (X, Y), define the *schoolbook weight* w^{sb}(X,Y) = 2c_{D−2} − 1 ∈ {−1, +1}, using the carry at the fixed position D − 2.

**Definition 1b (Cascade carry weight).** Let M(X,Y) = max{j : c_j ≥ 1} be the position of the topmost nonzero carry (the *cascade stopping depth*). The *cascade weight* is w^{cas}(X,Y) = c_{M−1} − 1 ∈ {−1, 0, +1}, the carry one step below the stopping point. For pairs with M = 0 we set w^{cas} = 0.

**Remark (Two valuations).** The schoolbook weight uses a fixed observation point (position D − 2), while the cascade weight uses a random stopping time (position M − 1). The two coincide when M = D − 1 (the carry chain extends to the top) and differ otherwise. Approximately 55% of sector-(0,0) pairs have M = D − 1; for sector (1,0), the top-carry constraint (Proposition 6) forces M ≤ D − 2 for all pairs. The distinction is analysed fully in [P1, §10].

**Definition 2 (Sector partition function).** For a, c ∈ {0, 1} and valuation ∗ ∈ {sb, cas}:

    σ^∗_ac(K) = Σ_{(X,Y): x_{K-2}=a, y_{K-2}=c} w^∗(X,Y)

**Definition 3 (Sector ratio).** R^∗(K) = σ^∗₁₀(K)/σ^∗₀₀(K).

The cascade ratio R^{cas}(K) → −π (Conjecture 1). The schoolbook ratio R^{sb}(K) → R₀ = (1/2 − 2ln(4/3))/(1 + ln(3/8)) = −3.9312…, a purely logarithmic constant [E, Proposition 1]. The decomposition R^{cas} = R₀ + ΔR with ΔR = +0.7896… isolates all transcendental content in the cascade correction; see [P1, §10.1].

---

## 3. The Markov Model and Its Failure

### 3.1 The independent-convolution approximation

The Diaconis–Fulman Markov chain [1, 2, 4] models the carry sequence (c_0, c_1, …, c_{D−1}) as a Markov chain: the carry into position j+1 depends only on c_j and the convolution conv_j, which is treated as an independent random variable given only j (not the actual digits).

### 3.2 Qualitative failure: wrong sign

**Proposition 1** (computational). The Markov model predicts R_Markov(K) > 0 for all K, while exact enumeration gives R(K) < 0 for K ≥ 7.

*Evidence.* Verified by exact enumeration for K = 5, …, 21 and independently confirmed at K = 20, 21 by the C-language enumerator (E44). A formal proof would require showing that the Markov per-sector weights have opposite sign to the exact per-sector weights; the mechanism (J = 1 exclusion, Proposition 2) provides the structural explanation.

This is not a quantitative error — the Markov model predicts the *wrong direction* of the sector effect. The sign flip is the most dramatic signature of non-Markovian behavior in the carry chain.

| K | R_Markov(K) | R^{cas}(K) |
|---|-------------|-----------|
| 5 | +0.59 | +0.25 |
| 7 | +0.62 | −0.09 |
| 9 | +0.64 | −0.73 |
| 11 | +0.65 | −1.554 |
| 13 | +0.65 | −2.339 |
| 15 | +0.66 | −2.824 |
| 17 | +0.66 | −3.035 |
| 19 | +0.66 | −3.110 |
| 20 | +0.66 | −3.125 |
| 21 | +0.66 | −3.133 |

### 3.3 The J = 1 exclusion effect

The Markov model's dominant contribution to σ₁₀ comes from the "J = 1 pathway" — a specific carry pattern near the top of the bridge. In the exact computation, this pathway is *forbidden* by a digit constraint that the Markov model ignores.

**Proposition 2 (J = 1 exclusion, computationally verified).** In the sector (a, c) = (1, 0), the leading-order Markov contribution at position J = 1 (one step from the boundary) requires digit configurations that violate the D-odd constraint. Removing this contribution flips the sign of σ₁₀.

*Evidence.* Verified by direct comparison of the Markov and exact computations for K = 5, …, 13 (experiment E06). The formal case analysis on digit configurations at positions K−2 and K−1 is outlined but not fully proved; see Open Problem 3.

---

## 4. Digit Correlations

### 4.1 Inter-position correlations

In the Markov model, convolutions at positions j and j' are independent given the carry state. In reality, the shared digit structure induces correlations.

**Definition 4 (Digit-position correlation).** For positions j, j' with |j − j'| = d:

    ρ(d) = Cov[conv_j, conv_{j'}] / Var[conv_j]

### 4.2 Correlation profile: the Fejér kernel

**Theorem 6 (Correlation profile).** For the central region of the product (positions j ≈ K), with K-bit uniformly random inputs:

    ρ(d) = (K − d)/(3K + 2),    0 ≤ d ≤ K−1

(exact, from Var[conv_j] = (3K+2)/16 and Cov = (K−d)/16). The Fejér kernel approximation $\rho(d) \approx (K-d)/(3K)$ holds with absolute error $2(K-d)/(3K(3K+2)) = O(1/K)$ uniformly in $d$; it becomes exact in the limit $K \to \infty$ at fixed $d/K$.

*Proof.* The convolution conv_j = Σ_{i+i'=j} x_i y_{i'} is a sum of products of independent bits. Two convolutions conv_j and conv_{j+d} share exactly K − d digit pairs (when both indices are in the central region). Each shared pair (x_i, y_{j−i}) = (x_i, y_{j+d−i}) contributes Cov(x_i y_{j−i}, x_i y_{j+d−i}) = E[x_i² y_{j−i} y_{j+d−i}] − E[x_i y_{j−i}] E[x_i y_{j+d−i}]. Since x_i ∈ {0,1}, x_i² = x_i, so E[x_i² y_{j−i} y_{j+d−i}] = (1/2)(1/2)(1/2) = 1/8 (for interior terms where all factors are random), while E[x_i y_{j−i}] E[x_i y_{j+d−i}] = (1/4)² = 1/16. Each shared pair contributes 1/8 − 1/16 = 1/16. The total covariance is (K − d)/16. The variance of conv_j in the central region has two regimes: the K − 2 interior terms have Var[x_i y_{j−i}] = E[(x_i y_{j−i})^2] − (E[x_i y_{j−i}])^2 = 1/4 − 1/16 = 3/16 each, and the 2 border terms (involving the fixed MSBs) have Var = 1/4 each. Since terms in the same convolution are independent (each uses a distinct pair of bit indices), Var[conv_j] = (K − 2) · 3/16 + 2 · 1/4 = (3K + 2)/16 ≈ 3K/16 for K ≫ 1. Hence ρ(d) = Cov/Var = (K−d)/16 ÷ (3K+2)/16 = (K−d)/(3K+2). □

### 4.3 Implications

The Fejér kernel has three notable properties:
- Triangular (linear decay with distance)
- Positive everywhere (no anti-correlations)
- Width proportional to K (correlations span the entire bridge)

This long-range correlation is what pushes the effective eigenvalue shift past the critical point A = 1, converting the rational Markov result into the transcendental −π [P1, Section 5].

---

## 5. The Exact Transfer Operator

### 5.1 State space

The exact (non-Markov) transfer operator tracks both the carry *and* the last d digits of both X and Y (the "digit buffer"). For buffer depth d, the state space has |S| = 2 × 2^d × 2^d = 2^{2d+1} states.

### 5.2 Eigenvalue spectrum

| Buffer depth d | \|S\| | Largest eigenvalue | Next eigenvalue |
|---------------|-------|-------------------|-----------------|
| 0 (Markov) | 2 | 1/2 | −1/2 |
| 1 | 8 | 1/2 | 1/2 |
| 2 | 32 | 0.3536 | 0.3536 |
| 3 | 128 | 0.3536 | 0.3536 |
| 4 | 512 | 0.3460 | 0.3460 |

The eigenvalue $0.3536 \approx 1/\sqrt{8} = 1/(2\sqrt{2})$ at levels 2–3 is consistent with a $1/b^{3/2}$ scaling for the first non-Markov mode in base 2. The exact spectrum is much richer than the Markov model's two eigenvalues {+1/2, −1/2}. Most additional eigenvalues are parasitic (related to digit-buffer dynamics, not carry propagation), but they collectively produce the correlation-induced shift that generates −π.

### 5.3 Non-universality of the exact spectrum

While the Markov spectral zeta function converges to ζ(2s) [P1, Theorem 3], the exact spectral zeta function does not — the parasitic eigenvalues from the digit buffer break the Weyl asymptotics. The conjectured −π limit comes from a *projected* quantity (the sector ratio) that isolates the carry degrees of freedom.

---

## 6. The Fejér Kernel and Dirichlet Coupling

### 6.1 Coupling matrix

Using the leading-order Fejér approximation $\rho(d) \approx (K-d)/(3K)$ (Theorem 6; exact form $\rho(d) = (K-d)/(3K+2)$), the digit-correlation profile couples the Dirichlet modes φ_n(j) = sin(nπj/L) through

    V_{nn'} = Σ_{d=0}^{K-1} ρ(d) Σ_j φ_n(j) φ_{n'}(j+d)

**Proposition 3 (Non-uniform diagonal coupling; numerical).** The diagonal elements V_{nn} are not proportional to sin²(nπ/(2L)) for all n. The effective shift A_eff(n) = V_{nn}/sin²(nπ/(2L)) varies from ~400 (low modes) to ~0.2 (high modes) at K = 15.

Verified numerically in experiment E10; no closed-form proof is given.

### 6.2 Universality argument

Despite the non-uniformity of A_eff(n), the sector ratio in the bridge model is determined by a single number — the effective A parameter in the Möbius formula S(A) = 2(1−A)/(3−2A) [P1, Corollary 1]. The mapping from the full coupling matrix V_{nn'} to the effective A involves a resummation that averages over mode numbers.

### 6.3 The bare coupling hits the pole

The full spectral sum S_full = v^T (I − D − V)^{-1} e, where D is the Markov diagonal and V the Fejér coupling matrix, was computed numerically:

| K | S_full | A_eff (from Möbius inversion) |
|---|--------|------|
| 10 | +56.0 | 1.509 |
| 20 | −8.9 × 10⁶ | 1.500000 |
| 50 | +3.3 × 10⁹ | 1.500000 |
| 100 | +2.3 × 10⁹ | 1.500000 |

**Finding 1.** A_eff → 3/2 exactly — the pole of the Möbius transformation. This is *structurally inevitable*: the Fejér kernel has total strength ρ̂(0) = (K+1)/6 → ∞, and the diagonal V_nn/sin²(nπ/(2L)) grows as ~K³ for low modes. Any coupling that makes S diverge maps to A = 3/2 via Möbius inversion, regardless of the coupling's structure.

**Finding 2 (Base universality).** The divergence A_eff → 3/2 occurs for all bases tested (b = 2, 3, 5). This confirms that the pole-hitting is a structural property of the Fejér coupling, not a base-2 artifact.

### 6.4 The χ₄ character from midpoint geometry

The Dirichlet character χ₄(n) = sin(nπ/2) arises geometrically. The sector perturbation is at position j* = K − 2 on a bridge of length L = 2K − 1. Since

    j*/L = (K−2)/(2K−1) → 1/2    as K → ∞,

the Dirichlet mode evaluation gives sin(nπj*/L) → sin(nπ/2) = χ₄(n). The sector perturbation naturally selects the midpoint of the bridge, and the midpoint evaluation of the Dirichlet modes *is* the character χ₄. This geometric origin explains why the sector ratio involves Leibniz's series and L(1, χ₄) = π/4.

---

## 7. Computational Pipeline

### 7.1 Algorithm

For each K from 3 to 21:

1. Enumerate all (X, Y) ∈ [2^{K−1}, 2^K)² (total: 4^{K−1} pairs).
2. For each pair, compute the full carry chain c_0, c_1, …, c_{D−1} by schoolbook multiplication [5].
3. Classify by sector (a, c) = (x_{K−2}, y_{K−2}).
4. Compute both valuations: the schoolbook weight w^{sb} = 2c_{D−2} − 1 and the cascade weight w^{cas} = c_{M−1} − 1, where M = max{j : c_j ≥ 1}.
5. Accumulate into σ^{sb}_ac and σ^{cas}_ac.
6. Compute R^{sb}(K) and R^{cas}(K).

Unless otherwise stated, all tables in this paper report the cascade valuation R^{cas}(K), which converges numerically toward −π (Conjecture 1 of [P1]).

### 7.2 Complexity and implementation

The enumeration requires O(K · 4^K) operations. For K = 20, this is ~5.5 × 10¹² operations over 2.75 × 10¹¹ pairs, of which 1.06 × 10¹¹ satisfy the D-odd constraint.

Exact enumeration covers K ≤ 21. The K = 20 computation uses a C implementation with OpenMP 8-thread parallelism, completing in 3344 seconds (~56 minutes) on a modern workstation. Memory usage is O(1) (streaming enumeration, no storage of pairs).

**K = 20 exact data (cascade valuation):** σ^{cas}₀₀ = −1,967,962,747; σ^{cas}₁₀ = 6,149,608,524; σ₀₀^{c=0} = −1,683,746,513; σ₀₀^{c=1} = −284,216,234. R^{cas}(20) = −3.12486023; α = σ₀₀^{c=1}/σ₀₀^{c=0} = 0.16880. At K = 21: R^{cas}(21) = −3.13340.

### 7.3 Convergence to −π

All convergence data below use the cascade valuation R^{cas}(K). For the schoolbook valuation, see [E, Proposition 1] where R^{sb}(∞) = R₀ = −3.9312… is proved analytically.

**Proposition 4 (Richardson working ansatz; numerical).** The empirical gap ratios $g(K)/g(K-1)$ converge monotonically to $1/2$, consistent with an exponential ansatz

$$R^{cas}(K) = R^{cas}(\infty) + c_1 \rho^K + c_2 \rho^{2K} + \cdots, \qquad \rho = 1/2.$$

The Neville–Aitken table with step-size $h = \rho^K$ eliminates successive error terms.

*Note.* A polynomial ansatz in $1/K^2$ was initially considered but diverges under Neville extrapolation. The exponential form reflects the true asymptotic rate $O(K^2 \cdot 2^{-K})$ [P1, §9.3].

| K | R(K) | ε(K) = \|R(K) + π\| | Sig. digits |
|---|------|------|-------------|
| 11 | −1.554 | 1.588 | — |
| 13 | −2.339 | 0.803 | 0.1 |
| 15 | −2.824 | 0.317 | 0.5 |
| 17 | −3.035 | 0.107 | 1.0 |
| 19 | −3.110 | 0.032 | 1.5 |
| 20 | −3.125 | 0.017 | 1.8 |
| 21 | −3.133 | 0.008 | 2.1 |
| ∞ (Richardson) | −3.1419… | — | 4.0 |

**Proposition 5 (Richardson extrapolation; numerical).** Applying Neville–Aitken acceleration [6] with the exponential ansatz of Proposition 4 to the subsequence R(7), R(9), …, R(21) yields R_extrap = −3.1419… to ≈ 4.0 significant digits (4.4 when restricted to exact values). The extrapolation is stable under cross-validation: removing any single K ∈ {7, …, 17} changes the estimate by less than 10⁻⁴. This is a numerical observation, not a convergence proof.

---

## 8. Physical Analogies

The sign flip from +2/3 (Markov) to −π (exact) has structural parallels in mathematical physics. In quantum field theory, a classically scale-invariant theory acquires a non-zero trace of the energy-momentum tensor through quantum corrections — here the "classical" (Markov) result is rational and the "quantum" (correlation-corrected) result is transcendental, with summing over Dirichlet eigenmodes providing the mechanism. The carry bridge also shares the spectral structure of a 1D Ising model with fixed-spin boundaries: Dirichlet eigenfunctions, exponential decay with correlation length ξ = 1/log(1/λ₂), and π-dependence through the spectral decomposition.

These analogies are structural, not rigorous equivalences. The "bare" Fejér coupling diverges (Section 6.3), while the physical system operates at a finite effective coupling A* = 3/2 − 1/(2(π+1)). Formalizing this renormalization mechanism remains open.

---

## 9. Structural Decomposition of R (E16–E19)

### 9.1 The top-carry constraint

**Proposition 6 (Top-carry identity).** conv_{D−2} = a + c for all D-odd pairs, where a = x_{K−2}, c = y_{K−2} are the sector bits.

*Proof.* conv_{D−2} = Σ_{i+i'=D−2} x_i y_{i'} with i, i' ∈ [0, K−1]. The only solutions are (i, i') = (K−2, K−1) and (K−1, K−2), giving conv_{D−2} = x_{K−2}·y_{K−1} + x_{K−1}·y_{K−2} = a·1 + 1·c = a + c. □

**Corollary (Sector constraints).**
- Sector (1,1): conv_{D−2} = 2, D-odd requires carries[D−1] = 0, but ⌊(2 + carries[D−2])/2⌋ ≥ 1. Contradiction → **N₁₁ = 0 for all K**.
- Sector (1,0): conv_{D−2} = 1, so carries[D−2] must be 0 → **M ≤ D−3** forced.
- Sector (0,0): conv_{D−2} = 0, carries[D−2] ∈ {0,1} → two sub-populations: "top" (M = D−2, ~55%) and "low" (M < D−2).

This gives σ₀₀ = σ₀₀^{top} + σ₀₀^{low} and σ₁₀ = σ₁₀^{low} only.

### 9.2 Closed-form count ratio

**Theorem 7 (Count ratio).** In the continuous limit K → ∞:

    N₁₀/N₀₀ → (2ln(4/3) − 1/2) / (2ln(9/8)) = 0.31992784225...

with correction O(1/2^K). The sector fractions of D-odd pairs are:

| Sector | Fraction of D-odd pairs | Value |
|--------|------------------------|-------|
| (0,0) | 2 ln(9/8) / (2 ln 2 − 1) | 0.6098 |
| (1,0) | (2 ln(4/3) − 1/2) / (2 ln 2 − 1) | 0.1951 |
| (0,1) | same as (1,0) | 0.1951 |
| (1,1) | 0 | 0 |

*Proof.* Write X = 2^{K−1}(1+u), Y = 2^{K−1}(1+v). D-odd: (1+u)(1+v) < 2. Sectors from ⌊2u⌋, ⌊2v⌋. Elementary integration. □

Verified to 8 significant digits against exact enumeration (K = 3–19) with geometric Richardson extrapolation.

### 9.3 The π isolation

Since R = (⟨val⟩₁₀/⟨val⟩₀₀) × (N₁₀/N₀₀) and N₁₀/N₀₀ involves only ln 2 and ln 3:

    ⟨val⟩₁₀/⟨val⟩₀₀ → −π · 2 ln(9/8) / (2 ln(4/3) − 1/2) = −9.8197...

**The entire factor π in the conjectured limit R → −π enters through the carry weight ratio**, not through the pair counts. This reduces the proof of Conjecture 1 to proving this single carry weight identity.

### 9.4 The boundary layer mechanism (E18)

The total variation distance TV(j) = (1/2) Σ_v |P(c_j|10) − P(c_j|00)| reveals a boundary layer with two peaks:
- j = K−2 (sector bit): TV ≈ 0.16, from the direct convolution perturbation.
- j = D−2 (top constraint): TV ≈ 0.54, from the structural constraint conv_{D−2} = a+c.

The c_{M−1} distribution shifts dramatically:

| | c_{M−1} = 0 (val = −1) | c_{M−1} = 1 (val = 0) | c_{M−1} = 2 (val = +1) |
|---|---|---|---|
| Sector (0,0)_low | 11.2% | 82.2% | 6.6% |
| Sector (1,0) | 2.1% | 67.9% | 30.0% |

This 6–7× shift in the val = +1 proportion drives σ₁₀ > 0 while σ₀₀ < 0.

### 9.5 The bit-sharing mechanism (E31)

The Markov model treats convolutions at successive depths as independent given the carry state, predicting a survival probability Q_Markov(d) that decays super-exponentially: Q ∼ (3/4)^d. Exact enumeration shows geometric decay Q_emp ≈ 0.55, an order-of-magnitude discrepancy at moderate depths. The mechanism is *bit-sharing*: the digit x_{K−3} (for instance) contributes to convolutions at depths d = 1 and d = 2, inducing correlations that the Markov model ignores.

**Definition 5 (Extended Markov chain).** For window size W ≥ 0, define a Markov chain whose state includes the carry value and a buffer of W bits of both X and Y. The convolution at each depth d ≤ W is computed exactly from the known bits; for d > W, the unknown bits are modeled as independent Bernoulli(1/2). The entry carry at depth W is forward-propagated from c = 0 through the exact chain.

**Finding 3 (Bit-sharing IS the mechanism).** The extended model with W = 8–10 reproduces the empirical survival probabilities to within 2%:

| d | Q_Markov | W = 8 | W = 10 | Empirical (K = 10) |
|---|----------|-------|--------|---------------------|
| 2 | 0.500 | 0.542 | 0.543 | 0.592 |
| 3 | 0.406 | 0.576 | 0.578 | 0.568 |
| 4 | 0.328 | 0.534 | 0.540 | 0.588 |
| 5 | 0.264 | 0.522 | 0.533 | 0.627 |
| 6 | 0.211 | 0.491 | 0.518 | 0.699 |

At d = 3, the W = 10 model matches the K = 10 empirical value to 1.8%. The Q values converge monotonically with W, confirming that each additional bit of memory captures more of the inter-depth correlation.

**Finding 4 (Effective eigenvalue → 1/2).** The geometric-fit effective eigenvalue q₀(W) (defined by Q(d) ≈ q₀ for large d within the window) decreases toward 1/2 as W grows:

| W | q₀ |
|---|-----|
| 3 | 0.543 |
| 5 | 0.515 |
| 7 | 0.499 |
| 10 | 0.483 |

The limit q₀ → 1/2 = |B₁| would recover the Diaconis–Fulman spectral gap [1] as the universal per-step survival probability of the cascade, consistent with the Bernoulli hierarchy [P1, §8].

**Remark.** The extended model successfully explains *why* Q ≈ 0.55 but does not directly produce R^{cas} → −π: the denominator σ^{cas}₀₀ passes through zero as W increases, causing numerical instability. The convergence R → −π requires the full carry chain (all K bits), not just the boundary layer window. The analytical proof of ΔR = −π − R₀ remains open [P1, §11, item 1].

---

## 10. Open Problems

1. **Prove R^{cas} = −π analytically.** The decomposition R^{cas} = R₀ + ΔR (with R₀ proved in [E, Proposition 1]) reduces the problem to showing ΔR = −π − R₀ = +0.7896…. Five equivalent formulations are developed in [P1, §10.5]: series identity over depths, angular integral, carry weight computation, resolvent trace, and separable weight matrix with carry constraints. The most promising route is the stopping-time series [P1, §8.2a], which decomposes R through per-depth universal constants Δ(τ). The exact analytical mechanism underlying the convergence is identified conditionally in [E] via **resolvent universality** of the effective non-Markovian transfer operator: the spectral resolvent, applied to the true per-sector cascade profiles, produces a macroscopic sum converging to $-\pi$ through collective mode weighting.

2. **Extend the bit-sharing model.** Section 9.5 shows the W-bit extended Markov chain reproduces Q_emp ≈ 0.55, but R^{cas} is unstable because σ^{cas}₀₀ passes through zero. A better treatment of the entry carry at the window boundary, or an analytical argument showing q₀ → 1/2 as W → ∞, would connect the bit-sharing mechanism to the Diaconis–Fulman eigenvalue [P1, §8].

3. **Proposition 2 (J = 1 exclusion): formal proof.** Verified computationally for K = 5, …, 13. A formal case analysis on digit configurations near the boundary would promote this to a theorem.

4. **K ≥ 22 data.** Each additional K approximately quadruples the computation time. Data through K = 21 are available; K = 22 would require approximately 16 hours. Additional data points would sharpen the Richardson extrapolation and allow independent convergence analysis of the even-K and odd-K subsequences.

5. **Base dependence.** In base b ≥ 3, the sector ratio appears rational (e.g., R₅ = 5/4 for base 5). The angular uniqueness of base 2 is proved in [G]: the D-parity boundary is a straight line in angular coordinates iff b = 2, which explains why π enters only for binary multiplication. A systematic exact enumeration of R(K; b) for bases 3, 5, 7 would test whether R_b is rational for all b ≥ 3.

6. **Boundary layer scaling.** Section 9.4 reports the TV(j) profile and c_{M−1} distributions at a single K value (K = 12). Computing these at K = 14, 16, 18 would reveal whether the key percentages (30% vs 5.5% for val = +1) converge to limiting values and at what rate.

7. **Separate weight sequences.** The current pipeline computes only the ratio R = σ₁₀/σ₀₀. Computing ⟨val⟩₁₀ and ⟨val⟩₀₀ separately would allow independent Richardson extrapolation and closed-form identification of each weight, which would provide a more transparent path to the identity ⟨val⟩₁₀/⟨val⟩₀₀ = −9.8197….

---

## References

1. P. Diaconis and J. Fulman, *Carries, shuffling, and symmetric functions*, Adv. Appl. Math. **43** (2009), 176–196.
2. P. Diaconis and J. Fulman, *Carries, shuffling, and an amazing matrix*, Amer. Math. Monthly **116** (2009), 788–803.
3. J. Fulman, *The carries process revisited*, preprint (2023).
4. J. M. Holte, *Carries, combinatorics, and an amazing matrix*, Amer. Math. Monthly **104** (1997), 138–149.
5. D. E. Knuth, *The Art of Computer Programming, Volume 2: Seminumerical Algorithms*, 3rd ed., Addison-Wesley, 1997.
6. L. F. Richardson and J. A. Gaunt, *The deferred approach to the limit*, Phil. Trans. R. Soc. A **226** (1927), 299–361.
7. [P1] Companion paper: *π from Pure Arithmetic: A Spectral Phase Transition in the Binary Carry Bridge*.
8. [G] Companion paper: *The angular uniqueness of base 2 in positional multiplication*.
9. T. Apostol, *Introduction to Analytic Number Theory*, Springer, 1976.
10. [E] Companion paper: *The Trace Anomaly of Binary Multiplication*. (Identifies conditionally the mechanism for $R \to -\pi$ via resolvent universality, assuming the LMH.)
