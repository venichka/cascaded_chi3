# Simulation suite for §11 and §12 of the cascaded-χ³ note

This directory contains 8 figures (PNG) and 9 Python scripts that
generate them. The figures are organised in two groups corresponding
to §11 (architecture comparison) and §12 (application landscape) of
the LaTeX note.

## Group 1: §11 Architecture comparison (3 figures)

### G1_A1_phi_comparison.png
**Where it goes**: §11.2 "Per-photon nonlinear phase: the Q-scaling".
After the boxed definition of R (eq. 37).
**Shows**: SPM phase shift vs input power for four architectures:
A (native only), B (cascade only with Δω<0), A+B constructive (Δω<0,
cascade reinforces native), A+B destructive (Δω>0, cascade cancels
native). Panel (b) shows the ratio (A+B)/A approaching 1+R = 11 in
the constructive case and 1−R = −9 in the destructive case at low
power, with saturation dropping both above |s₊|² ≈ 0.1.
**Caption suggestion**: "Architecture A (native-only) vs A+B
(combined), showing the sign-tunability of the cascaded contribution.
Open markers: TCMT measurement; solid lines: small-signal predictions.
The constructive (Δω<0) and destructive (Δω>0) configurations give
ratios exactly $1+R$ and $1-R$ in the low-power limit."

### G1_A2_phase_bandwidth.png
**Where it goes**: §11.3 "The phase-bandwidth locus". After the table
of the three families.
**Shows**: Panel (a) — LiNbO₃ phase-bandwidth locus with three families
overlaid. Native (blue) has slope −2 in log-log; cascade-with-fixed-Q₁
(green) has the sum-rule slope −1; cascade-with-fixed-Q₂ (red) parallel
to native at slope −2 but shifted up by the structural factor.
Panel (b) — material survey of R vs Q₂ for four cascade-capable
materials, with crossover Q₂ marked.
**Caption suggestion**: "Phase-bandwidth tradeoffs. (a) For LiNbO₃,
the cascade family (green/red) lies above the native family (blue) at
all bandwidths of practical interest. (b) Material survey of R vs Q₂
at Δω/γ₂=2 and η_geom=0.3."

### G1_A3_Q_scaling.png
**Where it goes**: §11.2, right after eq. (38)-(39) showing the
Q-scaling for both architectures, OR as a stand-alone "design knob"
illustration.
**Shows**: At fixed FM cavity Q₁ and fixed material parameters, sweep
the SH cavity Q₂. Native phase is independent of Q₂ (flat horizontal
line); cascade phase rises linearly. Panel (b) shows the ratio
(A+B)/A growing as 1 + R(Q₂) = 1 + α·Q₂.
**Caption suggestion**: "Cascade has a free extra Q-knob. Native
architecture is independent of Q₂; cascade rises linearly with Q₂.
Combined (A+B) follows the cascade once R ≳ 1."

## Group 2: §12 Applications (5 figures)

### G2_B2_bistability_comparison.png
**Where it goes**: §12.1.3 "Bistability threshold and all-optical
switching energy". After the boxed (A+B)/A = 1/(1+R) result.
**Shows**: Panel (a) hysteresis loops for four values of β (i.e. four
values of R from 0 to 0.8). Panel (b) ratio of measured threshold to
the R=0 value, falling on the predicted 1/(1+R) curve. Panel (c)
absolute threshold vs R with the same prediction.
**Caption suggestion**: "Bistability threshold lowered by the
cascade as 1/(1+R), in constructive configuration. Note that the
upper branch shows parametric instability at the highest β, where
|βa₁| approaches |Δ_SH| — the edge of TCMT validity."

### G2_B3_fwm_comparison.png
**Where it goes**: §12.1.4 "Four-wave-mixing wavelength conversion".
After the boxed (1+R)² result.
**Shows**: Panel (a) FWM conversion efficiency vs pump power, ∝ |s₊|⁴
in both architectures, with cascade shifted up. Panel (b) ratio
plateau at (1+R)² = 9 for R=2, dropping above |s₊|² ≈ 0.05 as the
cascade saturates.
**Caption suggestion**: "FWM is the most differentiating observable:
(1+R)² in cascade leverage. Plateau at low-to-moderate pump matches
prediction exactly; high-power roll-off is cascade saturation."

### G2_B4_saturable_absorber.png
**Where it goes**: §12.2.5 "Sign-tunable saturable absorber".
**Shows**: Panel (a) the dispersive (Kerr) and absorptive (TPA-like)
cascade coefficients vs Δ_SH, illustrating that the imaginary part
peaks at Δ_SH=0. Panel (b) intensity-dependent FM transmission for
three Δ_SH (0, γ₂, 3γ₂) showing how the SA character can be turned
off by detuning. Panel (c) the loss/Kerr ratio crossing 1 at
|Δ_SH| = γ₂/2.
**Caption suggestion**: "Cascaded saturable absorber. (a) The
TPA-like loss peaks Lorentzian at Δ_SH=0 (peak 2|β|²/γ₂) and is
suppressed by detuning. (c) Loss/Kerr ratio at |Δ_SH|=γ₂/2 marks
the operating boundary."

### G2_B5_AM_modulator.png
**Where it goes**: §12.2.6 "All-optical amplitude modulator via Δω
tuning".
**Shows**: Panel (a) clean time-domain co-modulation of Δ_SH(t) and
FM transmission at Ω_mod=γ₂/2 in linear-response regime. Panel (b)
transfer function: low-pass roll-off at γ₂/2 with **AC peak at
Ω = Δ_SH,₀**. Panel (c) modulation depth vs working detuning, falling
as 1/Δ_SH².
**Caption suggestion**: "Δω(t) AM modulator. (b) Transfer function
has both a low-pass roll-off at Ω=γ₂/2 and an AC band-pass peak at
Ω = Δ_SH (the resonant cascade response from the sum-rule analysis).
Operating at Ω_mod ≈ Δ_SH gives 3× enhancement over DC."

### G2_B6_sensor_sensitivity.png
**Where it goes**: §12.2.7 "Δω-readout transducer for sensing".
After the boxed Q₂² scaling.
**Shows**: Panel (a) the Lorentzian derivative d α_3_casc / d Δω,
with TCMT points exactly on the analytic curve. Zero slope at
|Δω|=γ₂/2, peak magnitude at Δω=0. Panel (b) **peak sensitivity vs Q₂
on log-log: TCMT points on Q₂² line (black), with linear-cavity Q₂
scaling for comparison (gray dashed)**. Panel (c) phase sensitivity
dφ_NL/dΔω vs working detuning showing the same Lorentzian-derivative
shape.
**Caption suggestion**: "Sensor sensitivity scales as Q₂², not Q₂.
Peak sensitivity at Δω=0 (panel a); the Q₂² scaling (panel b) reflects
the Lorentzian-derivative origin of the response: an extra factor of
1/γ₂ beyond the coefficient itself."

## Code structure

- `core.py`: shared infrastructure — material parameter dataclass,
  R-ratio helpers, plotting style. Imports the TCMT integrator from
  the earlier sim suite at `/home/claude/sim/cascaded_chi3_simulation.py`.
- `G1_A*.py`, `G2_B*.py`: one self-contained script per figure. Each
  has a `main` guard and can be run independently.
- All figures land in `figures/` with names matching the script names.

## Two notes from doing the simulation

1. **Note's worked example has an arithmetic inconsistency.** Row 655
   of §10 claims Δω = 2γ₂ = 2ω₁/Q₂ = 6×10¹¹ rad/s. The middle equality
   uses γ₂ = ω₁/Q₂ instead of γ₂ = ω₂/Q₂ = 2ω₁/Q₂. With the correct γ₂,
   either Δω = 6×10¹¹ rad/s ≡ γ₂ (not 2γ₂), giving R≈25 at the edge of
   adiabatic validity, OR the operating point should be Δω = 2γ₂ =
   1.22×10¹² rad/s giving R≈12. The simulations use the latter
   (genuinely safe adiabatic, R≈12).

2. **AM modulator has a band-pass character, not just low-pass.** The
   transfer function in B5(b) has both a low-pass roll-off at γ₂/2
   AND a resonant peak at Ω_mod = Δ_SH (the "AC peak" from the
   sum-rule analysis: |α₃_eff(Ω)| peaks when the modulation frequency
   matches the SH detuning in the rotating frame). This is new
   physics that wasn't in the original §12.6 derivation — operating
   at Ω_mod ≈ Δ_SH gives a ~3× sensitivity boost. Worth mentioning
   in the AM modulator subsection.

3. **Sensor sensitivity peak is at Δω=0, not Δω=γ₂/(2√3).** My early
   message in this conversation said the peak was at Δω/(2√3); the
   note correctly fixes this to Δω=0 (eq. 51 region). The simulation
   confirms this.

4. **B4 panel (b) "linear cavity" reference line** at T=1 is what the
   FM transmission would be far from resonance. On resonance with
   critical coupling, T=0 in the linear case; with γ_rad = γ₁ (the
   normalisation used in this code), the linear-cavity reflection
   coefficient at resonance is unity. The figure label is slightly
   misleading and should say "uncoupled cavity / far off resonance"
   rather than "linear cavity".
