# Questions for the metasurface designer — SH-only samples 1–4

We've run the singly-resonant cascaded-χ³ model on the four samples (`reports/SH_only_sample_analysis.md`).
The **rankings, lineshapes, and Q/η scalings are already robust**; the items below pin the **absolute
magnitudes** and refine the resonance/field model. Ordered by impact.

## A. Sets the absolute scale (most important)

1. **One absolute SHG efficiency would settle it.** Could you give a single absolute SHG conversion
   efficiency — η_SHG = P(2ω)/P(ω) at a stated pump intensity (W/cm²) and detuning — from a nonlinear
   simulation of any one sample? This calibrates our overlap-normalisation "bridge" directly and removes
   the largest absolute-scale unknown in one number.
   - *If a nonlinear sim isn't handy:* what exactly are the field amplitudes in the denominator of your
     `overlap` definition — field **maxima, RMS, or energy-normalised eigenmode** amplitudes? And is the
     `∬ … dA` a **2-D cross-section at one z-plane** or **integrated through the film thickness**? (That
     fixes the bridge analytically.)

2. **What is the `symmetry` field?** Every `.mat` has a complex `symmetry` array (≈0.02–0.04·area near
   resonance) that isn't in your field list (freq/S21/overlap/area). Is it a second (DFG/parity) overlap,
   the symmetry-forbidden component, or a diagnostic? Should we use it, or is `overlap` the only nonlinear
   coupling that enters the cascade?

## B. Refines the resonance model

3. **Is Q(2ω) loaded or intrinsic?** Our FWHM read of your SH transmission dips gives ≈0.6× the tabulated
   Q (consistently, all four samples). The model needs the **total/loaded** energy-decay rate γ₂. Is the
   quoted Q the loaded value or the intrinsic/unloaded one?

4. **Radiative vs absorptive split.** What is γ_rad vs γ_nonrad (or the coupling regime) per SH resonance?
   The dips differ a lot — S1/S2 reach ~10⁻⁴ (near-critical), S4 only ~0.16 (under-coupled) — which sets
   the SHG out-coupling and the saturable-absorber loss.

## C. Geometry / field (refines observable prefactors)

5. **FM near-field enhancement.** Can you give |E_z(ω)|_inside / |E_z(ω)|_incident inside the LiNbO₃ at
   the operating wavelength (or the local-field factor)? It sets the single-pass-phase prefactor; we
   currently assume T_geom = 1.

6. **qBIC mechanism and overlap headroom (our top design lever).** Confirm the SH mode is a symmetry-broken
   qBIC, and which asymmetry sets the radiative coupling (the H1/H2 base/pillar split, or a lateral one).
   **Can the overlap be pushed toward |overlap| ≳ 0.3 at Δω = ±2γ₂** — via stronger asymmetry and/or
   sub-domain (unit-cell) poling of χ⁽²⁾ — without collapsing Q? χ³_casc scales linearly with |overlap|²,
   so this is the single most valuable improvement.

## D. Confirmations / practical

7. **Table typo.** Please confirm sample-2 `overlap(−2γ) = 0.9000` is a typo for ≈ **0.0900** (your raw
   `.mat` gives 0.086).

8. **FM band featureless?** We see flat FM transmission |t|² ≈ 0.73 (no FM resonance). Confirm there is no
   FM-band feature near ω₁ for any sample — the SH-only model assumes the fundamental is non-resonant.

9. **Film & material.** Confirm total LiNbO₃ thickness = H_t = H₁(base) + H₂(pillars) = 180–220 nm. Are the
   films **MgO-doped or undoped congruent**? (The SH lands in the green, 488–532 nm — photorefractive risk.)

10. **Pump source.** What pump wavelength / source, and CW vs pulsed? The designs cluster at 975–1064 nm
    (S4 ≈ 1064 nm matches Nd:YAG). Can the SH resonance be **trimmed post-fabrication** to set Δω?
