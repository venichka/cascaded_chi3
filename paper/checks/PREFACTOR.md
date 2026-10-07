# The prefactor of the cascaded χ⁽³⁾: re-derivation and checks

*2026-10-06. Supports `paper/PLAN.md` §3 (F1, F2). Numerical checks: `prefactor_check.py` in this folder.*

## Summary

The boxed results of the main note (eq. `chi3casc`, §6 and §14) and of the SH-only note (eq. `chi3-1R-boxed`) are **4× too small**. The correct, normalisation-invariant results (Boyd amplitudes) are:

$$
\chi^{(3)}_\text{casc}
 = -\frac{2\varepsilon_0\omega_1}{3}\,
   \frac{|M|^2}{N_2\,(\Delta\omega+i\gamma_2/2)\,\int|\mathbf f_1|^4 d^3r},
\qquad
\chi^{(3)}_{\text{casc},1R}
 = -\frac{2\varepsilon_0\omega_1}{3}\,
   \frac{|\widetilde M|^2}{N_2\,(\Delta\omega+i\gamma_2/2)\,V_\text{cell}},
$$

with
- $M=\int\chi^{(2)}_{ijk}f^*_{2,i}f_{1,j}f_{1,k}\,d^3r$
- $\widetilde M=\int\chi^{(2)}_{ijk}f^*_{2,i}\hat e_j\hat e_k\,d^3r$
- $N_n=\int\varepsilon_0\tilde\varepsilon_n|\mathbf f_n|^2 d^3r$

In explicit form,

$$
\chi^{(3)}_\text{casc}\simeq-\frac{2\,\omega_1}{3\,n_2^2\,\Delta\omega}\,|\chi^{(2)}|^2\,\eta_\text{geom},
\qquad
\eta_\text{geom}\equiv\frac{\varepsilon_0 n_2^2|M|^2}{|\chi^{(2)}|^2 N_2\int|\mathbf f_1|^4},
$$

where $\eta_\text{geom}=1$ is exactly the bulk bound-wave limit.

**The cause is two slips.** In the projection of the wave equation, the main note writes β = ε₀k₂²c²M/(4ω₂) where its own preceding line gives /(2ω₂) (line 521). It then "adopts" ε₀ω₁/4 instead of the resulting ε₀ω₁/2 (line 535). The bare Kerr coefficient α₃ = 3ε₀ω₁/8 (line 544) carries the same overall deficit, justified by an "additional factor of 1/2" (line 552) that does not exist in Boyd's convention. Since χ³_casc ∝ |β|²/α₃, the net error is a factor of 4.

## 1. Conventions

**Boyd's convention:** Ẽ(r,t) = E(r)e^{−iωt} + c.c., with no ½. With it:
- P(2ω₁) = ε₀χ²E₁E₁ (SHG)
- P(ω₁) = 2ε₀χ²E₂E₁* (DFG)
- P(ω₁) = 3ε₀χ³|E₁|²E₁ (SPM)

These are the relations the notes use, and the tabulated χ² = 2d₃₃ and χ³ values belong to this convention.

The main note *states* the other convention, ℰ = ½(E e^{−iωt} + c.c.) (lines 96–97). In that convention the same relations would read P(2ω) = ½ε₀χ²E², P(ω) = ε₀χ²E₂E₁* and P(ω) = ¾ε₀χ³|E|²E. **The stated convention and the formulas used do not match (F2).**

- **χ³_casc is unaffected.** It is a ratio of material constants, and either convention gives the same value when used consistently.
- **The energy reading of |a|² is affected.** In Boyd's convention ⟨Ẽ²⟩ = 2|E|², so the stored energy is U = ∫(ε₀ε̃|E|² + μ₀|H|²) = 2N|a|². The note's normalisation ½∫[ε₀ε̃|f|² + μ₀|h|²] = 1 (line 156) means N = 1 and therefore U = 2|a|², not the stated U = |a|² (line 660). Every metasurface quantity quoted in joules or watts inherits this factor of 2.

## 2. Projection of the wave equation

**Equations.** The positive-frequency part of ∇×∇×Ẽ + c⁻²∂ₜ²(εẼ) = −μ₀∂ₜ²P̃, with Ẽ = a(t)f(r)e^{−iΩt} + c.c., ∇×∇×f = (ω_n²/c²)εf, slowly varying a, and Ṗ neglected, is

$$\varepsilon_0\varepsilon\,[(\omega_n^2-\Omega^2)a-2i\Omega\dot a]\,\mathbf f=\Omega^2\mathbf P .$$

**Projection.** Project onto f* and use (ω_n² − Ω²)/(2Ω) ≃ ω_n − Ω:

$$\boxed{\;\dot a=i(\Omega-\omega_n)\,a+\frac{i\Omega}{2N}\int\mathbf f^*\!\cdot\mathbf P\,d^3r\;}$$

**Cross-check.** This is the standard first-order perturbation result, δω = −(ω/2)∫E*·δP/∫ε₀ε|E|² (Johnson *et al.*, PRE 65, 066611 (2002)), applied to a nonlinear δP.

**Where the note departs.** In the note's notation the LHS is −2iω₂ȧ₂/ε₀ and the RHS is k₂²Ma₁² (lines 505–515). That gives β = ε₀c²k₂²M/(2ω₂) = ε₀ω₂M/2 = ε₀ω₁M. Line 521 writes /(4ω₂).

## 3. Coupled-mode coefficients

| Process | Source P | Coefficient |
|---|---|---|
| SHG into a₂ (Ω = 2ω₁) | ε₀χ²a₁²f₁f₁ | **β = ε₀ω₁M/N₂** |
| DFG into a₁ (Ω = ω₁) | 2ε₀χ²a₂a₁*f₂f₁* | iβ₁*a₁*a₂ with **β₁ = ε₀ω₁M/N₁** |
| SPM into a₁ | 3ε₀χ³\|a₁\|²a₁\|f₁\|²f₁ | **α₃ = (3ε₀ω₁/2)·∫χ³\|f₁\|⁴/N₁** |

**Manley–Rowe.** With the table's coefficients, N₁ d|a₁|²/dt = 2Im(ε₀ω₁M a₁²a₂*) = −N₂ d|a₂|²/dt. So U₁ + U₂ = 2(N₁|a₁|² + N₂|a₂|²) is conserved, and the photon numbers obey dn₂/dt = −½dn₁/dt. The DFG factor 2 is what makes β₁N₁ = β₂N₂. This is the note's point, and it survives the correction.

**Conclusion.** With N = 1 the note has β/4 and α₃/4: the same overall deficit in both couplings. That is why the TCMT dynamics in the code (dimensionless) never noticed it. The absolute scale enters only when χ³_casc, R or a physical-unit observable is formed.

## 4. Adiabatic elimination and χ³_casc

**Steady state.** In the steady state, a₂ = −βa₁²/(Δω + iγ₂/2), so the FM equation gains

$$i\beta_1^*a_1^*a_2=-i\,\frac{\varepsilon_0^2\omega_1^2|M|^2}{N_1N_2(\Delta\omega+i\gamma_2/2)}\,|a_1|^2a_1 .$$

**Matching.** Match this to α₃ evaluated with a uniform χ³ (the note's definition):

$$\chi^{(3)}_\text{casc}=-\frac{2\varepsilon_0\omega_1}{3}\,\frac{|M|^2}{N_2(\Delta\omega+i\gamma_2/2)\int|\mathbf f_1|^4}.$$

N₁ cancels, and N₂ cancels the f₂-scaling of |M|², so the result is normalisation-invariant. At N₂ = 1 it is 4× the note's −ε₀ω₁|M|²/(6Δω∫|f₁|⁴).

**Design-point form.** At Δω = 2γ₂ = 4ω₁/Q₂ (with Δω ≫ γ₂/2):

$$|\chi^{(3)}_\text{casc}|=\frac{Q_2}{6n_2^2}|\chi^{(2)}|^2\eta_\text{geom}\quad(\text{note: }Q_2/24n_2^2).$$

## 5. Uniform-medium limit: the bulk bound wave

**Exact bulk.** A plane wave E₁e^{ik₁z} in a uniform medium drives E₂e^{2ik₁z}. The Helmholtz equation gives exactly

$$E_2=\frac{\chi^{(2)}E_1^2}{n_1^2-n_2^2},\qquad
P(\omega_1)=2\varepsilon_0\chi^{(2)}E_2E_1^*
\;\Rightarrow\;
\chi^{(3)}_\text{casc,bulk}=-\frac{2|\chi^{(2)}|^2}{3(n_2^2-n_1^2)}.$$

This is the note's eq. `chi3casc-bulk` with F = 1, since cn₂Δk ≃ ω₁(n₂² − n₁²).

**Resonator formula in the same geometry.** Take a periodic box of volume V with single plane-wave modes (f_n = 1/√(ε₀n_n²V), N = 1, ω₂ = 2ω₁n₁/n₂, so n₂² − n₁² ≃ n₂²Δω/ω₁). §4 gives −2ω₁|χ²|²/(3n₂²Δω), the same to first order in n₂ − n₁, with η_geom = 1. The note's formula gives ¼ of it.

## 6. Singly resonant (SH-only)

**Derivation.** The FM is a prescribed local field E_loc ê.
- SH amplitude: a₂ = −ε₀ω₁E_loc²M̃/[N₂(Δω + iγ₂/2)].
- Cell-integrated DFG polarisation: ∫ê·P = 2ε₀a₂E_loc*M̃* = −2ε₀²ω₁|M̃|²|E_loc|²E_loc/[N₂(Δω + iγ₂/2)].
- Matching to 3ε₀χ³|E_loc|²E_loc V_cell gives the χ³_casc,1R of the Summary.

**Comparison with the note.** The SH-only note uses β̂ = (ε₀ω₁/4)E_loc²M̃ (SH-only note, lines 165–174), which is 4× too small. There is no α₃ in this chain, so the factor 4 enters through β̂ alone and the 1R result is also 4× too small.

## 7. Photon units and the microresonator literature

**Mapping.** In Boyd's convention with N = 1, U = 2|a|² = ħω n, so a = √(ħω/2)â. The SH equation becomes ḃ = i(β/2)√(ħω₁)â². Comparing with H/ħ = g[(a†)²b + a²b†], which gives ḃ = −iga², yields:

$$|g|=\tfrac12\beta\sqrt{\hbar\omega_1}=\frac{\varepsilon_0\omega_1}{2}\sqrt{\hbar\omega_1}\,|M|
=\sqrt{\frac{\hbar\omega_1^2\,(2\omega_1)}{8\varepsilon_0}}\;
\frac{|\int\chi^{(2)}_{ijk}E_{1i}E_{1j}E^*_{2k}|}{\int\varepsilon_1|E_1|^2\sqrt{\int\varepsilon_2|E_2|^2}} .$$

**Literature check.** This is exactly the standard microresonator coupling rate, g = √(ħω_a²ω_b/8ε₀)·∫χ²E_a*E_a*E_b/(∫ε_a|E_a|²√∫ε_b|E_b|²) (Zhao *et al.*, InGaP, arXiv:2105.12705, Eq. (1)), with χ = 2d as in Boyd's convention. The note's chain (β/4 together with the |a|² = U mapping g = β√(ħω₁/2)) gives a g that is **2√2 too small**.

**Consequence for the microring note.** Its scalar first-principles estimate of g for Lu *et al.*'s device becomes 0.56 × 2.83 = **1.6 MHz** at Θ = 1, against 1.78 MHz from Lu's full calculation and 1.2 MHz measured, i.e. Θ ≈ 0.75. This is consistent. It replaces the "factor 2.1 below" paragraph and caveat 8.

**Native Kerr check.** α₃ = (3ε₀ω₁/2)χ³∫|f₁|⁴ gives a per-photon shift U = α₃ħω₁/2 = 3ħω₁²χ³ζ/(4ε₀n⁴V). This is exactly the microring note's U_nat = ħω₁²c n₂^Kerr/(n_g²V_eff), which was derived from the round-trip phase with n₂^Kerr = 3χ³/(4ε₀cn²). The note's α₃ = 3ε₀ω₁/8, with its |a|² = U mapping, gives half of that.

## 8. Numerical checks (`prefactor_check.py`)

The setup is a 1-D Bloch-periodic two-layer cell (n = 2.2 / 1.45, with χ² only in the high-index layer). Every check solves the wave equation exactly, by finite differences, with no modal truncation.

Each test reports ratio = exact/formula. "Corrected" is the formula of this document; "note" is the note's formula.

| Test | What is exact | Compared with | Δω/2ω₁ = 2×10⁻² | 5×10⁻³ | 10⁻³ |
|---|---|---|---|---|---|
| A: β | Driven SH field from the full Helmholtz solve, projected on the resonant mode | −βa₁²/Δω | 1.0101 / 4.04 | 1.0025 / 4.01 | **1.0005 / 4.00** |
| B: full 2R chain | Exact SH field → exact DFG polarisation → χ³_eff by the note's matching definition | χ³_casc (§4) | 1.0077 / 4.03 | 1.0019 / 4.01 | **1.0004 / 4.00** |
| C: 1R chain | Prescribed uniform FM field (Bloch k = 0) → exact SH field → cell-integrated DFG polarisation | χ³_casc,1R (§6) | 1.109 / 4.43 | 1.026 / 4.10 | **1.005 / 4.02** |

Entries are corrected / note.

| Test | What is exact | Compared with | Result |
|---|---|---|---|
| D: α₃ | Self-consistent nonlinear eigenfrequency with a native local Kerr (iterated eigenproblems, no projection formula) | −α₃\|a\|² | \|a\| = 0.1: **0.9994** (note: 3.997); \|a\| = 0.2: **0.99999** (note: 4.0000) |
| E: uniform limit | Bulk bound wave | χ³_casc (§4) | **1.011** (note: 4.04); the residual is first order in n₂ − n₁ |

How to read the table:
- The residual in A–C is the expected O(Δω/ω) linearisation error plus the contribution of non-resonant SH modes. It vanishes as Δω → 0, and the ratio to the note's formula converges to **4**.
- C converges more slowly because a uniform drive at k = 0 also excites a strongly overlapping lower band.
- The first version of test D used |a| = 10⁻³. There the shift (~10⁻¹⁰) was below the eigensolver's round-off, and the run gave meaningless ratios (1.08 and 0.68). At |a| ≥ 0.1 the shift is resolved.
- Full output: `prefactor_check_output.txt`.

## 9. Independent re-derivation

A separate agent derived (1)–(6) from scratch. It was given only the conventions and the task, not this document, the notes or the expected answer. Its results:

| Quantity | Blind derivation | This document |
|---|---|---|
| Projection prefactor | iΩ/(2N_n) | same |
| β₂, β₁ | ε₀ω₁M/N₂, ε₀ω₁M*/N₁ (N₁β₁ = N₂β₂*) | same |
| α₃ | (3ε₀ω₁/2N₁)∫χ³\|f₁\|⁴ | same |
| χ³_casc | −(2/3)ε₀ω₁\|M\|²/(N₂Δω∫\|f₁\|⁴) | same |
| χ³_casc,1R | −(2/3)ε₀ω₁\|M̃\|²/(N₂ΔωV_cell) | same |
| Stored energy | U = 2N\|a\|² (Boyd) | same |
| g | ½ε₀M√(ħω₁³), i.e. K = 8 in √(ħω_a²ω_b/(Kε₀ε_a²ε_bV)) | same; matches the literature formula of §7 |
| Uniform limit | Reproduces the exact bound wave (sympy: difference 0 with the unlinearised detuning) | same to first order |

**Additional numerical checks from the blind derivation:**
1. A 1-D inhomogeneous cavity with all-mode driven SH: ratio 1.020, 1.00076 and 1.000076 at Δω/ω₁ = 3×10⁻², 10⁻³ and 10⁻⁴.
2. **A time-domain integration of the two-mode Maxwell equations** with the full ∂ₜ²P term (no slowly-varying approximation), using deliberately non-unit N₁ = 1.7 and N₂ = 0.55. The measured FM shift agrees with TCMT to 0.9995–1.0016 over six cases, with the correct sign: blue shift for Δω > 0.
3. A singly-resonant 1-D case: 1.041, 1.0041 and 1.0004 at Δω/ω₁ = 10⁻², 10⁻³ and 10⁻⁴.
4. Conservation of N₁|a₁|² + N₂|a₂|² to 12 digits.

**One physics remark for the paper.** For a general tensor the bulk bound wave also has a longitudinal piece, set by ∇·D = 0. It is −(2/3)d_z²/n₂² with d_m = χ²_mpq e_p e_q. It is non-resonant, it is dropped by single-mode TCMT, and it vanishes for transverse d. Mention it in one sentence alongside the Bosshard rectification channel as a non-resonant background.

## 10. Verdict

**The factor of 4 is confirmed.** The evidence:

- two independent analytic derivations;
- the exact uniform-medium limit;
- five frequency-domain numerical tests (2R chain, 1R chain, α₃ by an exact nonlinear eigen-solve);
- a time-domain Maxwell integration without the slowly-varying approximation;
- the standard microresonator coupling formula from the literature;
- the microring note's independently derived native-Kerr U_nat.

The note's boxed χ³_casc, its SH-only counterpart, and every R derived from them should be multiplied by 4. Quantities in joules or watts additionally need the U = 2N|a|² correction (F2).
