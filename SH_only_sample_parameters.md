# SH-only sample parameter extraction → model inputs

**Scope.** How the metasurface designer's deliverables (`data/sample1-4/` + slide-1 of
`presentation/Presentation3.pptx`) are turned into inputs for the singly-resonant
(SH-only) cascaded-χ³ model of `Cascaded_chi3_SH_only_note.tex`. Every number below is
produced and cross-checked by `cascaded_chi3/sample_params.py` (run it to reproduce).

The four files are **not** four different designs — they are the realised geometries for
the four samples 1–4 of the earlier "singly resonant" table (verified by f₀ and Q below).

---

## 0. Final per-sample model parameters (the table the simulations consume)

| Sample | file | λ₁ (pump) | λ₂ (SH) | Q₂ | Δf = γ₂/2π | n₁ = nₑ(λ₁) | n₂ = nₑ(λ₂) | L_meta = H_t | A_cell | fill | η_geom(+2γ) | η_geom(−2γ) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | design4_alpha9 | 996.7 nm | 498.3 nm | 3442 | 174.8 GHz | 2.159 | 2.250 | 190 nm | (400 nm)² | 0.298 | 0.0144 | 0.0088 |
| 2 | design4_H1_70 | 974.9 nm | 487.5 nm | 3486 | 176.4 GHz | 2.161 | 2.256 | 180 nm | (400 nm)² | 0.273 | 0.0111 | 0.0074 |
| 3 | design4_H1_120_H2_100 | 1057.6 nm | 528.8 nm | 4641 | 122.1 GHz | 2.156 | 2.236 | 220 nm | (400 nm)² | 0.389 | 0.0374 | 0.0303 |
| 4 | design4_H1_140_H2_80 | 1064.5 nm | 532.2 nm | 5633 | 100.0 GHz | 2.156 | 2.234 | 220 nm | (400 nm)² | 0.421 | 0.0435 | 0.0442 |

Material constants (supplied by us, **not** the designer; same for all samples):
`|χ⁽²⁾| = 2 d₃₃ = 5×10⁻¹¹ m/V`, `χ⁽³⁾_native = 2×10⁻²¹ m²/V²`, bridge factor `C_bridge = 1` (baseline).

> **These designs operate at λ₁ ≈ 1 µm**, not the 1.55 µm the note/`params.py` were anchored at.
> The material constants barely move (see §5), but every per-photon-phase / damage number must be
> re-evaluated at ~1 µm.

---

## 1. Inputs and provenance

| Source | What it gives |
|---|---|
| `data/sample1-4/design4_*.mat` | per-sample COMSOL output (below) |
| `data/sample1-4/ref_idx_LN.m` | LiNbO₃ Sellmeier table nₑ, nₒ vs λ [µm], 0.4–5 µm |
| slide 1 (`image4.png`) | f₀(2ω), Q₂, γ, overlap(±2γ), area, H1/H2/H_t table |
| slide 1 (`image6.png`) | **definition** of `overlap` (see §6) |
| slide 1 (`image7.png`) | axis convention: `E_z` = out-of-plane = `d₃₃`, `nₑ` |

**`.mat` fields** (each a column over a frequency sweep):

| field | meaning | use |
|---|---|---|
| `freq` | SH frequency [THz] | resonance axis; `2ω₁ = 2π·freq` |
| `freq_fun` | FM frequency [THz] = `freq/2` | sanity check |
| `S21_2`, `S21_phase` | SH transmission \|t\|², phase | Q cross-check, coupling regime |
| `S21_2_fun`, `S21_phase_fun` | FM transmission \|t\|², phase | confirms FM is featureless |
| `area` | LiNbO₃ cross-section [m²] | overlap normalisation, filling factor |
| `overlap` | **complex**, freq-resolved | → η_geom (the nonlinear coupling) |
| `symmetry` | **complex**, freq-resolved | **undocumented — open item (§11)** |

---

## 2. Sample ↔ file ↔ geometry

`H_t = H₁ (base slab) + H₂ (pillars)` = total LiNbO₃ thickness. The two-pillar asymmetry per
unit cell is what makes the SH a radiating qBIC. Confirmed against the table:

| Sample | file | f₀(2ω) [THz] | Q₂ | H₁ (base) | H₂ (pillars) | H_t |
|---|---|---|---|---|---|---|
| 1 | design4_alpha9 | 601.5751 | 3442 | 80 | 110 | 190 nm |
| 2 | design4_H1_70 | 615.0085 | 3486 | 70 | 110 | 180 nm |
| 3 | design4_H1_120_H2_100 | 566.9170 | 4641 | 120 | 100 | 220 nm |
| 4 | design4_H1_140_H2_80 | 563.2812 | 5633 | 140 | 80 | 220 nm |

---

## 3. Conventions and unit conversions

- **Frequencies.** Table f₀(2ω) is the **SH** centre. Pump ω₁ = π·f₀; ω₂ = 2π·f₀ = 2ω₁.
  Wavelengths λ[µm] = 299.792 / f[THz]: λ₁ from f₀/2, λ₂ from f₀.
- **Linewidth.** Designer's `γ[THz] = f₀/Q₂` is the **FWHM in ordinary frequency**. The model's
  `γ₂` is the SH **energy decay rate in rad/s**: `γ₂ = ω₂/Q₂ = 2π·γ[THz]`. Hence the **usable Kerr
  bandwidth** `Δf = γ₂/2π = γ[THz]` directly (note §`sec:sum-1R`).
- **Boyd / degeneracies.** SHG deg. 1, DFG deg. 2, SPM deg. 3 (note §2) — already baked into the
  boxed χ³_casc; nothing to set per sample.
- **Tensor axis.** `E_z` (out-of-plane) ⇒ the active component is `χ⁽²⁾_zzz = 2 d₃₃` and the relevant
  index is the **extraordinary** `nₑ`. So `n₁ = nₑ(λ₁)`, `n₂ = nₑ(λ₂)`.

---

## 4. Direct parameters + consistency checks

Computed in `sample_params.py` §(A). `γ = f₀/Q₂` reproduces the tabulated γ to **0.0 %** for all four,
and `freq_fun = freq/2` exactly — the table is internally consistent.

| Sample | ω₁ [rad/s] | ω₂ [rad/s] | γ₂ [rad/s] | λ₁ | λ₂ | n₁ | n₂ |
|---|---|---|---|---|---|---|---|
| 1 | 2π·300.79 THz | 2π·601.58 THz | 2π·0.1748 THz | 996.7 nm | 498.3 nm | 2.1593 | 2.2499 |
| 2 | 2π·307.50 THz | 2π·615.01 THz | 2π·0.1764 THz | 974.9 nm | 487.5 nm | 2.1608 | 2.2560 |
| 3 | 2π·283.46 THz | 2π·566.92 THz | 2π·0.1221 THz | 1057.6 nm | 528.8 nm | 2.1559 | 2.2355 |
| 4 | 2π·281.64 THz | 2π·563.28 THz | 2π·0.1000 THz | 1064.5 nm | 532.2 nm | 2.1555 | 2.2338 |

---

## 5. Material parameters (supplied by us)

- **Indices** `nₑ` interpolated from the designer's `ref_idx_LN.m` at λ₁, λ₂ (§4 table). Note that
  `n₂ ≈ 2.234–2.256` is essentially the note's old value (2.26) — the wavelength change barely moves it.
- **`|χ⁽²⁾| = 2 d₃₃ = 5×10⁻¹¹ m/V`.** d₃₃(LiNbO₃) ≈ 25 pm/V (telecom/1064-nm SHG). At SH ≈ 488–532 nm
  Miller's rule lifts this ~10–20 %; baseline keeps the note value, treated as ±20 %.
- **`χ⁽³⁾_native = 2×10⁻²¹ m²/V²`** (note baseline at 1.55 µm). At ~1 µm the 2-photon energy at ~500 nm
  (≈2.48 eV) approaches the ~3.9 eV gap, so this may be somewhat **larger** — which would *reduce* R.
  Flagged as the main material uncertainty (R ∝ 1/χ⁽³⁾_native).

These two material numbers and `C_bridge` (§6) set the **absolute** scale of every result; the
per-sample *ranking* and the Q/η scaling are far more robust than the absolute magnitudes.

---

## 6. The overlap → η_geom mapping (the crux)

**Designer's definition** (slide `image6.png`):

```
overlap = (1/Area) · ∬ E_z(ω) E_z(ω) E_z*(2ω) dA  /  ( |E_z(ω)| |E_z(ω)| |E_z(2ω)| )
```

with the field amplitudes in the denominator acting as scalar mode-norms (so the stored `.mat`
`overlap` has area units ~10⁻¹⁴ m², and `|overlap|/area` is the dimensionless number she tabulates).

**Mapping to the model.** This is exactly the SH-only SHG overlap `M̃ = ∫ χ⁽²⁾_zzz f*_{2,z} ê ê dV`
of note `eq:Mtilde-def`, normalised. The pump field `E_z(ω)²` plays the role of `ê ê` (and is actually
*better* than the note's uniform-pump approximation, since it keeps the real near-field profile), and
`E_z*(2ω)/|E_z(2ω)|` is the SH mode direction `f*_{2,z}`. Because χ³_casc ∝ |M̃|² (note `eq:chi3-1R-explicit`):

```
η_geom^(1R)(Δω) = ( |overlap_complex(f)| / area )² · C_bridge ,
                  evaluated at the SH frequency f = f₀ + Δω/(2π).
```

Two powers of the overlap (it enters χ³ twice: SHG up-conversion **and** DFG back-conversion, equal by
Kleinman). `C_bridge = 1` baseline absorbs the O(1) mismatch between her **2-D, mode-norm-normalised**
overlap and the model's **3-D Brillouin-normalised** `M̃` — see the caveat below.

**Do not double-count the resonance.** Her normalisation divides out `|E_z(2ω)|`, i.e. the resonant SH
*amplitude*. So `overlap` is the **spatial-coherence factor only**; the resonant build-up (the Lorentzian
`1/(Δω + iγ₂/2)` and the Q₂ enhancement) is supplied separately by the model. This is why `overlap` is
*small at exact resonance* (Δω = 0: 0.02–0.04, a Fano normalisation dip) and grows with detuning — that
frequency dependence is spatial, not the resonant lineshape.

**Why ±2γ.** The two tabulated values are `overlap` evaluated at the conservative operating points
Δω = ±2γ₂ (f = f₀ ± 2γ[THz]); the two signs are the self-focusing (Δω<0) and self-defocusing (Δω>0)
branches. The full `.mat` array lets us evaluate η at any Δω for a sweep.

**Verification (`sample_params.py` §C).** Reproducing `|overlap_complex|/area` at f₀±2γ from the raw
arrays matches the slide table to within ±6 % on every point **except** sample 2's −2γ:

| Sample | calc(−2γ) | table(−2γ) | calc(+2γ) | table(+2γ) |
|---|---|---|---|---|
| 1 | 0.094 | 0.100 | 0.120 | 0.120 |
| 2 | 0.086 | **0.900** ⚠ | 0.106 | 0.100 |
| 3 | 0.174 | 0.180 | 0.193 | 0.190 |
| 4 | 0.210 | 0.200 | 0.209 | 0.220 |

⚠ **Sample 2 `overlap(−2γ) = 0.9000` is a decimal typo for ≈ 0.0900** (raw-data value 0.086). Used the
raw value, not the table value.

> **Normalisation bridge (dominant absolute-scale uncertainty).** `C_bridge` is currently 1. To pin it
> we need either (i) the designer's exact field-norm convention in the `overlap` denominator (max vs RMS
> vs energy-normalised; 2-D slice vs z-integrated), or (ii) **one absolute SHG conversion-efficiency
> number** from her, which calibrates `M̃` directly and removes the bridge entirely. Until then, absolute
> χ³_casc carries an O(1)–O(few) uncertainty (on top of the |χ⁽²⁾|, χ⁽³⁾_native uncertainties).

---

## 7. Operating point (Δω)

Δω ≡ 2ω₁ − ω₂ is **our knob** (set by the pump wavelength), not a design parameter; tuning the pump
shifts λ₁ by only ~γ₂/ω₁ ~ 10⁻⁴, negligible. Default = conservative **Δω = +2γ₂** (note default:
adiabatic clean, TPA suppressed ×1/16). Sweeps run Δω/γ₂ ∈ [−5, 5] with η_geom(Δω) from the overlap
array, trustworthy for |Δω| ≳ γ₂ (single-pole validity, note §`sec:validity` item 1).

---

## 8. Geometry for observable quantities

- **L_meta = H_t** = 180–220 nm (total LiNbO₃ thickness) → single-pass phase length for φ_NL, ΔT/T, FWM.
- **A_cell = (400 nm)² = 1.6×10⁻¹³ m²** (slide 1: period 400 nm for all).
- **`area`** (LiNbO₃ cross-section) = 4.4–6.7×10⁻¹⁴ m² → **filling factor 0.27–0.42**.
- **T_geom** (local-field factor, note `eq:Eloc-magnitude`): baseline **1**. The FM transmission
  `|t|² ≈ 0.73` (≈15 % reflection) is consistent with T_geom ~ O(1); a sharper value needs the actual
  inside/incident field ratio (open item §11, optional).
- **V_cell = A_cell · L_meta** (only needed if we use the boxed `|M̃|²/V_cell` form instead of the
  data-derived η_geom).

---

## 9. Verification summary

| Check (`sample_params.py`) | Result |
|---|---|
| γ = f₀/Q₂ vs tabulated γ | ✅ 0.0 % all four |
| freq_fun = freq/2 | ✅ exact |
| `.mat` `area` vs table `Area` | ✅ exact all four |
| FM transmission flat (SH-only premise) | ✅ \|t\|² ≈ 0.73, ripple <1.3 % |
| overlap(±2γ) reproduces table | ✅ ±6 %, except S2 −2γ typo (0.9000→0.0900) |
| SH dip depth (coupling) | S1,S2 → ~10⁻⁴ (near-critical); S3 → 0.034; S4 → 0.16 (under-coupled) |
| loaded Q from FWHM vs quoted Q | ⚠ Q_fit ≈ 0.6× quoted (Fano/loaded — see §11) |

---

## 10. Preliminary derived outputs (baseline C_bridge = 1, Δω = +2γ₂)

From `chi3_casc(Δω) = −ω₁|χ⁽²⁾|² η_geom / [6 n₂² (Δω + iγ₂/2)]` (`sample_params.py` §E):

| Sample | η_geom(+2γ) | \|χ³_casc\| [m²/V²] | R = \|χ³_casc\|/χ³_native | Δf [GHz] | \|χ³\|·Δf |
|---|---|---|---|---|---|
| 1 | 0.0144 | 9.9×10⁻²² | **0.49** | 174.8 | 1.7×10⁻¹⁰ |
| 2 | 0.0111 | 7.7×10⁻²² | **0.39** | 176.4 | 1.4×10⁻¹⁰ |
| 3 | 0.0374 | 3.5×10⁻²¹ | **1.76** | 122.1 | 4.3×10⁻¹⁰ |
| 4 | 0.0435 | 5.0×10⁻²¹ | **2.48** | 100.0 | 5.0×10⁻¹⁰ |

Sign flip across Δω = 0 confirmed: Re χ³(+2γ) < 0 (defocusing), Re χ³(−2γ) > 0 (focusing); Im χ³ > 0
(two-photon-like loss). **Ranking S4 > S3 > S1 > S2.**

> **Headline caveat.** The data-derived η_geom (0.011–0.044) is **2–10× below** the note's *assumed*
> η_geom^(1R) ≈ 0.1, so R^(1R) ≈ 0.4–2.5 rather than the note's ~4. The cascade is a **modest** material
> enhancement for these designs (only S3, S4 beat native Kerr), and — as the note already warns — the
> single-pass per-photon Kerr phase will be µrad-level. The most *observable* experiments are therefore
> SHG conversion efficiency, the Re χ³ sign-flip vs detuning, and pulsed saturable absorption, not CW Kerr
> phase. Quantifying all of these per sample is the next step.

---

## 11. Open items / assumptions

| # | Item | Impact | Action |
|---|---|---|---|
| 1 | `symmetry` complex field (undocumented; \|·\|/area ≈ 0.02–0.04 at resonance) | unknown | ask designer what it is |
| 2 | `C_bridge` (overlap norm convention) | O(1)–O(few) on absolute χ³_casc | ask norm convention **or** one absolute SHG-efficiency point |
| 3 | Quoted Q loaded vs intrinsic (FWHM check gives ~0.6×) | γ₂, and χ³_casc ∝ Q₂ | confirm Q is total/loaded |
| 4 | `χ⁽³⁾_native`, `d₃₃` at ~1 µm / 500 nm SH | R ∝ \|χ⁽²⁾\|²/χ³_native | refine from literature (Miller's rule) |
| 5 | `T_geom` | φ_NL, ΔT/T prefactor (not χ³_casc itself) | optional: inside/incident field ratio |

None of these block the simulations; items 2–4 set the absolute scale and are flagged in every output.

---

## 12. Effects to compute next (SH-only observable suite, per sample)

Survive in 1R (single-pass / material observables) — all to be computed vs Δω and pump power:
χ³_casc(Δω) Re/Im, R^(1R), two-photon-like loss / saturable absorption ΔT/T, self-focus↔defocus sign
flip, sum-rule \|χ³\|·Δf, per-photon phase φ_NL (+ damage-threshold check), broadband cascaded FWM
efficiency, and the SHG conversion efficiency (the cleanest direct measurement + bridge calibrator).
Die in 1R (need an FM cavity): bistability, all-optical switching, parametric limiter, AM-modulator
readout (note §`sec:applications-1R`).

*Reproduce everything here with* `conda run -n ml_base python cascaded_chi3/sample_params.py`.
