# Host-material constants: cascaded vs native Kerr, pump 1.55 µm, SH 0.775 µm

> **How the paper uses this (2026-10-07):**
> - The LiNbO₃ range quoted in the paper is R = 11–30. It uses n₂ᴵ = 1.8–3.2e-19 and excludes the cascade-containing 0.91e-19 Z-scan.
> - GaP is presented as an upper estimate, since d₁₄ is taken at 1.064 µm.
> - The table adds an n₁ column (pump index; called n0 below).

> **Update 2026-11-09 (literature check; no adopted value changed):**
> - **LiNbO₃ n2_I now has primary sources.** He et al., Optica 6, 1138 (2019) (`he2019`) inferred 1.8e-19 m²/W from the OPO threshold of a z-cut ring, quasi-TE (ordinary), simulated A_eff, no error bar. This is the earliest traceable origin of 1.8e-19. Shams-Ansari et al., APL Photon. 7, 081301 (2022) (`shamsansari2022`, App. A–B) measured 1.74e-19 (TM, extraordinary) and 1.61e-19 (TE, ordinary) at 1550 nm by pump–probe XPM in a z-cut TFLN ring. They estimate a cascaded SFG term of opposite sign "below 10%" and do not subtract it. The adopted 1.8e-19 stays (R 20 → 20.7 with 1.74e-19). Supplement §5 cites both.
> - **GaP d14: still unverified.** Search-engine extracts of Shoji 1997 Table 11 give 36.8 / 70.6 / 159 pm/V at 1.313 / 1.064 / 0.852 µm, with Δ36 = 0.53 / 0.85 / 1.31. The extracted columns were jumbled, so the wavelength alignment is unconfirmed. If 36.8 pm/V at 1.313 µm is right, Miller scaling gives 33.7 pm/V at 1.55 µm and F(GaP) falls 2.8× (R 0.91 → 0.33). GaP papers at 1.55 µm all assume a value: 50 pm/V (Pantzas 2022, Anthur 2020), or χ2 = 82 pm/V (Wilson 2018 arXiv, unsourced). **Needs a library copy of Shoji 1997 Table 11.** Until then the paper's "upper estimate" wording stands, and the lower value would only widen LiNbO₃'s lead.
> - At 1.55 µm the bulk LiNbO₃ cascade is −3.3 to −4.8e-19 m²/W (e-pol, d33 = 18.8–22.7 pm/V), comparable to the native term. A bulk Z-scan at 1.55 µm would therefore be dominated by cascading; only waveguide and ring measurements give a near-native value.

Compiled 2026-10-07. Every derived number below is reproduced by `compute.py`
(`conda run -n ml_base python compute.py`; full output in `compute_output.txt`).
Bibliographic data for every key are in `new_refs.bib` and were checked against
Crossref. `bache2012` is the exception: it is an arXiv preprint, checked with the arXiv API.

## Conventions

- Boyd convention, chi2_ijk = 2 d_il. |chi2| = max |chi2_ijk u_i* v_j v_k| over
  complex unit vectors. A numerical maximisation in `compute.py` gives:
  - **2 d33** for LiNbO3 (3m) and AlN (6mm). The result holds for either sign of
    d31 and for |d22| ≤ 3 pm/V.
  - **(2/√3)·2 d14 = 1.1547 × 2 d14** for zincblende, with the optimum along [111].
- n2 = smallest principal index of the host at 0.775 µm. This is n_e for LiNbO3
  (negative uniaxial) and **n_o for AlN** (positive uniaxial).
- n0 = index at 1.55 µm for the pump polarisation, which is extraordinary for the
  d33 hosts.
- chi_nat = 4 eps0 c n0² n2_I / 3, F = |chi2|²/(n2² chi_nat),
  R = Q2 η F/6 = 200 F (Q2 = 4000, η = 0.3), Q2(R=1) = 6/(η F) = 20/F.

## Summary (adopted values)

| Host | d at 1.55 µm (pm/V) | \|chi2\| (pm/V) | n2 (SH) | n0 (pump) | n2_I (m²/W) | chi_nat (m²/V²) | F | R | Q2(R=1) |
|---|---|---|---|---|---|---|---|---|---|
| LiNbO3 (cong., e) | d33 = 18.8 | 37.5 | 2.178 (n_e) | 2.138 (n_e) | 1.8e-19 | 2.91e-21 | 0.102 | 20 | 197 |
| GaP | d14 = 56.1 | 129.7 | 3.195 | 3.054 | 1.1e-17 | 3.63e-19 | 4.5e-3 | 0.91 | 4.4e3 |
| AlN (e, d33) | d33 = 4.04 | 8.1 | 2.140 (n_o) | 2.160 (n_e) | 3.5e-19 | 5.78e-21 | 2.5e-3 | 0.49 | 8.1e3 |
| Al0.2Ga0.8As | d14 = 92.6 | 213.9 | 3.533 | 3.267 | 2.6e-17 | 9.82e-19 | 3.7e-3 | 0.75 | 5.4e3 |

Range across sources (combining the low and high ends of d and n2_I; details in the
per-host tables):

| Host | d range (pm/V) | n2_I range (m²/W) | F | R (Q2=4000) | Q2(R=1) |
|---|---|---|---|---|---|
| LiNbO3 | 18.8–22.7 | 0.91–3.2e-19 | 0.057–0.29 | 11–59 | 68–350 |
| GaP | 50–56 | 0.7–1.4e-17 | 2.8–7.1e-3 | 0.57–1.4 | 2.8–7.1e3 |
| AlN | 3.8–4.3 | 2.3–3.5e-19 | 2.1–4.3e-3 | 0.43–0.86 | 4.7–9.4e3 |
| Al0.2Ga0.8As | 79–112 | 1.5–2.6e-17 | 2.7–9.4e-3 | 0.54–1.9 | 2.1–7.4e3 |

LiNbO3 sensitivity, R at Q2 = 4000, with Q2(R=1) in brackets:

| d33 \ n2_I | 0.91e-19 | 1.8e-19 | 3.2e-19 |
|---|---|---|---|
| 18.75 (Miller from 1.313 µm) | 40 [99] | 20 [197] | 11 [349] |
| 22.67 (Miller from 1.064 µm) | 59 [68] | 30 [134] | 17 [239] |
| 25.2 (1.064-µm value, unscaled) | 73 [55] | 37 [109] | 21 [193] |

---

## LiNbO3 (congruent; extraordinary pump and SH; d33)

| Quantity | Value | λ | Source | Notes |
|---|---|---|---|---|
| d33 | 19.5 pm/V | SHG of 1.313 µm | shoji1997 | Congruent LN. From the abstract (verified). |
| d33 | 25.2 pm/V | SHG of 1.064 µm | shoji1997 | 25.7 at 0.852 µm. Shoji et al. state that Miller's rule does not reproduce their wavelength dependence. |
| d33 (adopted) | **18.75 pm/V** | 1.55 µm | Miller from 19.5 @ 1.313 | Uses chi1 = n_e² − 1 from zelmon1997, and the measurement closest to 1.55 µm. |
| d33 (alt.) | 22.67 pm/V | 1.55 µm | Miller from 25.2 @ 1.064 | Shoji's 1.064→1.313 drop is steeper than Miller predicts, so 18.75 is not a lower bound. The 25–27 pm/V often quoted at 1.55 µm in the TFLN literature appears to be the 1.064-µm value carried over without dispersion. |
| n2_I (adopted) | **1.8e-19 m²/W** | 1.55 µm | (none found) | The value used across TFLN Kerr-comb papers. I could not trace it to a measurement: kovach2020 attributes it to wang2019, but wang2019 uses 0.91e-19. The microring note gives no source for it either. |
| n2_I | 0.93e-19 m²/W (9.1 ± 1.3e-20) | 1.064 µm, e-pol | desalvo1996 | Z-scan with 30-ps pulses. The measured value still contains a negative cascaded χ(2):χ(2) contribution. wang2019 uses this number at 1.55 µm. |
| n2_I (non-cascaded) | 3.86 ± 0.32e-19 m²/W | 1.064 µm, e-pol | bache2012 (also used in phillips2011) | DeSalvo's value minus a computed n2,casc = −2.94e-19 (d33 = 25.2). Electronic plus DC Raman. Unpublished preprint. |
| n2_I (non-cascaded) | ≈2.6e-19 m²/W | 2.05 µm | jankowski2020 | Two-band scaling of the phillips2011 value. |
| n2_I (estimate) | 3.2–3.6e-19 m²/W | 1.55 µm | this work: sheikbahae1991 G2 scaling of 3.86e-19 | Uses Eg = 3.8–4.0 eV. The lower end scales the whole value; the upper end keeps a 50 % Raman part fixed. The same scaling gives 3.0e-19 at 2.05 µm, against 2.6e-19 quoted by jankowski2020, so the method is good to about 15 %. |
| n_e, n_o | n_e 2.1376 / 2.1784; n_o 2.2111 / 2.2587 | 1.55 / 0.775 µm | zelmon1997 | 21 °C, congruent. Check: n_e(1.064) = 2.1555, against 2.1558 used by bache2012. |
| Absorption | none intrinsic | — | — | The band gap is about 3.8–4 eV. Neither 2ħω1 = 1.6 eV, 2ħω2 = 3.2 eV, nor the pump–SH sum ħω1 + ħω2 = 2.4 eV reaches it. The practical issue at 775 nm is photorefraction in congruent material (MgO doping reduces it), plus the photorefractive and thermal drifts typical of TFLN. |

**Flags.**
- The current Supplement table uses n2 = 2.26 for LiNbO3. That is n_o(0.775 µm), but
  the paper defines n2 as the smallest principal index, which is n_e = 2.178.
- Its chi_nat = 2e-21 corresponds to n2_I = 1.24e-19 with n0 = n_e. That is
  inconsistent with the microring note's 3.1e-21, which comes from 1.8e-19 with
  n1 = 2.21 = n_o.
- The microring note calls 1.8e-19 "the conservative, large end" of a 0.9–1.8e-19
  range. That is not right. The 0.9e-19 value is a raw 1.064-µm Z-scan that
  includes cascading. The cascading-corrected values are larger: 3.9e-19 at
  1.064 µm and about 3.2e-19 at 1.55 µm.

## GaP (zincblende; d14)

| Quantity | Value | λ | Source | Notes |
|---|---|---|---|---|
| d14 (= d36) | 70.6 pm/V | SHG of 1.064 µm | shoji1997 | Corroborated only by citing papers: the OP-GaP OPO literature, and arXiv:2102.02406, which quotes d36 ≈ 70 pm/V near 1 µm. Shoji's semiconductor tables (Tables 11 and 13) are subscriber-only, and I could not read them. |
| d14 (adopted) | **56.1 pm/V** | 1.55 µm | Miller from 70.6 @ 1.064 | Uses Bond 1965 indices. The 532-nm SH lies above the indirect gap (2.26 eV), so the 1.064-µm value may be resonantly raised. |
| d14 (cross-check) | 50 pm/V | 1.55 µm | pantzas2022 | Value assumed for CW OP-GaP SHG at 1.55 µm, cited to the Shoji–Kondo–Ito 2002 compilation, which I did not read. |
| d14 (unverified) | 36.8 pm/V | SHG of 1.313 µm | attributed to shoji1997 by a search summary only | Not verified. If correct, Miller scaling gives 33.7 pm/V at 1.55 µm, and F falls by a factor of 2.8. |
| \|chi2\| | 129.7 pm/V | 1.55 µm | (4/√3) d14 | |
| n2_I (adopted) | **1.1(3)e-17 m²/W** | 1.56 µm | wilson2020 | GaP-on-insulator ring, CW modulation transfer normalised to the photothermal shift. A comb-threshold estimate gives 1.2(5)e-17, which is the arXiv-v1 abstract value. Wilson notes that pulsed measurements give lower values. |
| n2_I | 7e-18 m²/W | 1.04 µm | liu2010 (as tabulated in vermeulen2023) | Bulk (100) GaP, femtosecond Z-scan. This paper also measured 3PA at 1.04 µm. |
| n | 3.0536 / 3.1953 | 1.55 / 0.775 µm | bond1965 | Room temperature, tabulated data with a cubic spline. The temperature-dependent Sellmeier of Wei et al. 2018 (0.7–12.5 µm, Opt. Mater. Express 8, 485) was behind a captcha and not used. |
| Absorption | see notes | — | vurgaftman2001, vermeulen2023 | Indirect gap ≈ 2.26 eV (549 nm), direct gap ≈ 2.78 eV at room temperature. **Pump:** no 2PA (1.6 eV); phonon-assisted 3PA is possible (2.4 eV > 2.26 eV; liu2010). **SH:** linearly transparent, but degenerate 2PA at 775 nm is allowed (3.2 eV > 2.78 eV), with β ≈ 4–8e-11 m/W near 800 nm (vermeulen2023 Table 2B: Liu 2010 4.2e-11, Grinblat 2019 8e-11). Pump–SH non-degenerate 2PA (2.4 eV) lies above only the indirect gap, so it is phonon-assisted and weak. |

Kerr anisotropy caveat: zincblende n2 depends on field direction. The [111]
polarisation that maximises |chi2| is not the in-plane TE geometry of `wilson2020`.
The correction could be tens of percent and has not been applied.

## AlN (wurtzite; d33; extraordinary pump, TM-like)

| Quantity | Value | λ | Source | Notes |
|---|---|---|---|---|
| d33 | 4.3 ± 0.3 pm/V; \|d31\| = d33/(45 ± 5) | SHG of 1.030 µm | majkic2017 | Bulk PVT single crystals, Maker fringes referenced to a LiNbO3 plate. Taken from the abstract via two independent search snippets; full text not read (paywalled). |
| d33 (adopted) | **4.04 pm/V** | 1.55 µm | Miller from 4.3 @ 1.030 | Uses Pastrňák n_e. |
| d33, thin films | χ(2) = 6.2 ± 0.4 pm/V; 8.6 pm/V for MOCVD films at 1.064 µm (as cited there) | 1.55 µm | bruch2018 | From SHG in a crystalline AlN ring. Whether "χ(2)" means d33 or 2 d33 there is not stated clearly. Film values run higher than bulk; the bulk value is adopted. |
| \|chi2\| | 8.1 pm/V | 1.55 µm | 2 d33 | |
| n2_I (adopted) | **3.5e-19 m²/W** | ~1.56 µm | liu2018 | Crystalline AlN-on-sapphire. liu2020 uses this value for a TM00 mode at 1559 nm and gets a 44 mW threshold, against about 50 mW measured. Polarisation (TM, extraordinary) matches the d33 pump. How liu2018 obtained the value is unknown (paywalled). |
| n2_I | (2.3 ± 1.5)e-19 m²/W | ~1.55 µm | jung2013 | TE-like mode in the (001) plane, i.e. ordinary polarisation. Sputtered c-axis AlN on SiO2, comb/OPO threshold. |
| n_e, n_o | n_e 2.1602 / 2.1830; n_o 2.1199 / 2.1399 | 1.55 / 0.775 µm | pastrnak1966 | Sellmeier for 0.22–5 µm (refractiveindex.info record of this DOI). n2 = n_o(0.775) = 2.140; n0 = n_e(1.55) = 2.160. |
| Absorption | none intrinsic | — | vurgaftman2003 | The gap is about 6.0–6.2 eV, so there is no 2PA or 3PA at either wavelength (3ħω2 = 4.8 eV). Losses come from defects, impurities and the film. |

## Al0.2Ga0.8As (zincblende; d14)

| Quantity | Value | λ | Source | Notes |
|---|---|---|---|---|
| Composition | x = 0.20 | — | chang2020 | AlGaAs-on-insulator Kerr combs use x = 0.20 (chang2020). pu2016 used x = 0.17 with Eg = 1.63 eV, as summarised in kovach2020. |
| Band gap | 1.673 eV (741 nm) | 300 K | adachi1985 (1.424 + 1.247x) | The Gehrsitz model's gap parameter gives 1.715 eV (723 nm). The SH photon (1.600 eV) is therefore **74 meV (~33 nm) below the gap** (115 meV with the Gehrsitz value). The half-gap is 0.837 eV (1.482 µm), so the pump is 37 meV below it. |
| d14(GaAs) anchor | d14/d33(LN) = 5.01 ± 0.3 | SHG of 4.1 µm | skauli2002 | Combined with d33(LN) = 17.1 pm/V at 4.135 µm (Miller from Shoji's 25.2 pm/V), this gives d14(GaAs) = 85.6 pm/V. Skauli et al. quote 94 ± 10 pm/V using d33(LN) = 27 pm/V at 1.064 µm; that number is in Pinguet's Stanford thesis, not the Opt. Lett. abstract. |
| d14 (adopted) | **92.6 pm/V** | 1.55 µm | Miller from the GaAs anchor | Indices from skauli2003 (GaAs, 4.1 and 2.07 µm) and gehrsitz2000 (x = 0.2, 1.55 and 0.775 µm). Scaling across composition with Miller's rule follows ohashi1993. Starting from 94 pm/V gives 101.7. gili2016 assumes about 100 pm/V for Al0.18Ga0.82As at 1.55 µm. The near-gap SH could add resonant enhancement that Miller's rule misses, so the uncertainty is ±15–20 %. |
| \|chi2\| | 213.9 pm/V | 1.55 µm | (4/√3) d14 | |
| n2_I (adopted) | **2.6e-17 m²/W** | 1.55 µm | pu2016 (also quoted by chang2020 for x = 0.2) | I could not read pu2016, so whether the value was measured there or cited is unknown. |
| n2_I | (1.3–1.5)e-17 m²/W | 1.55 µm | aitchison1997, via vermeulen2023 | Al0.18Ga0.82As waveguides near the half-gap. vermeulen2023 calls this range the "most cited"; waveguide reports span 0.21–3.3e-17. For x = 0.2 rather than 0.17, the Eg⁻⁴ factor alone lowers n2 by about 8 %. |
| n | 3.2672 / 3.5329 | 1.55 / 0.775 µm | gehrsitz2000 | 20 °C. Cross-checks: 3.27 at 1.55 µm (Xie et al. 2020, arXiv:2004.14537; not in the bib); at x = 0, Gehrsitz agrees with skauli2003 to 0.2 % at 2–4 µm. The 0.775-µm value is close to the upper limit of the model's below-gap range. |
| Absorption | see notes | — | vermeulen2023, dolgaleva2011, vurgaftman2001 | **Pump:** degenerate 2PA is nominally forbidden (1.600 < 1.673 eV), but band-tail 2PA of 2.6–6e-13 m/W is reported below the half-gap. 3PA is 8e-26 m³/W² at 1550 nm (dolgaleva2011, Al0.2 core). **SH:** Urbach-tail linear absorption only 74 meV below the gap; the gap red-shifts by about 0.45 meV/K (GaAs Varshni parameters), so heating eats into the margin. Degenerate SH 2PA (3.2 eV) is allowed. **Pump–SH non-degenerate 2PA** (2.4 eV > Eg) is allowed: it is a loss on the pump proportional to the SH intensity, competing with Im chi_casc. |

---

## Uncertainty

- **LiNbO3, n2_I.** Reports and estimates span a factor of 3.5 (0.9–3.2e-19). The
  low end includes cascading; the high end is a cascading-corrected, two-band-scaled
  estimate.
- **LiNbO3, d33.** The value at 1.55 µm is uncertain by about 20 % (18.8–22.7 pm/V),
  depending on which Shoji wavelength is scaled. Together, F spans 0.06–0.29 and
  Q2(R=1) spans 70–350.
- **GaP.** Uncertainty is ±30 % in n2_I (wilson2020 quotes ±27 %) and about 10 % in
  d14 (50–56 pm/V). F spans 2.8–7.1e-3. The unverified 1.313-µm Shoji value would
  lower F by a further factor of 2.8.
- **AlN.** n2_I differs between polarisations and platforms, 2.3 ± 1.5 versus
  3.5e-19. The bulk d33 is ±7 %; thin films report up to about 2× larger values,
  which would raise F by up to 4×. F spans 2.1–4.3e-3 with the bulk d33.
- **Al0.2Ga0.8As.** n2_I is uncertain by a factor of 1.7 (1.5–2.6e-17), and the
  Miller-scaled d14 by about ±20 %. F spans 2.7–9.4e-3. The SH lies 74–115 meV below
  the gap, so near-resonant dispersion of both n(0.775 µm) and d14 is the main
  systematic.
- **Overall.** The three semiconductors (GaP, AlN, AlGaAs) are indistinguishable
  within these spreads, with F ≈ (2–9)e-3 and Q2(R=1) ≈ 2–9e3. At the adopted
  values LiNbO3 is 20–40× better; across the full spreads the factor is 6–140.

## Changes relative to the current Supplement table (Table S-hosts)

| Host | Current \|chi2\|, chi_nat, n2 | Current F / Q2(R=1) | New F / Q2(R=1) | Main cause |
|---|---|---|---|---|
| LiNbO3 | 50, 2e-21, 2.26 | 0.245 / 82 | 0.10 / 197 | d33 at 1.55 µm (18.8, not 25); chi_nat from 1.8e-19 with n_e; n2 = n_e(0.775) = 2.178 |
| GaP | 100, 1e-19, 3.1 | 0.0104 / 1.9e3 | 0.0045 / 4.4e3 | The old chi_nat implies n2_I = 3.0e-18, 3.6× below wilson2020; the (4/√3) norm factor partly compensates |
| AlN | 10, 2e-21, 2.1 | 0.0113 / 1.8e3 | 0.0025 / 8.1e3 | d33 = 4.0 (not 5); the old chi_nat implies n2_I = 1.2e-19, below both measurements |
| AlGaAs | 100, 1e-18, 3.3 | 0.00092 / 2.2e4 | 0.0037 / 5.4e3 | \|chi2\| = (4/√3) × 92.6 = 214 (not 100); n2 = 3.53 |

The main-text sentence (main.tex around l. 286) gives F = 0.24 / 0.01 / 0.01 / 9e-4
and "Q2 ≈ 80, 2e3, 2e4". With the adopted values it would read
F ≈ 0.10 (LiNbO3) and (2.5–4.5)e-3 for GaP, AlN and AlGaAs, with Q2(R=1) ≈ 200 for
LiNbO3 and ≈ (4–8)e3 for the semiconductors. AlGaAs is no longer an order of
magnitude worse than GaP and AlN. The project memory's "R ≈ 49" for LiNbO3 at
Q2 = 4000 becomes R ≈ 20, with a range of 11–59.

## Not verified / open items

1. Shoji 1997 Tables 11 and 13 (GaP and GaAs) are subscriber-only. GaP 70.6 pm/V at
   1.064 µm rests on secondary citations. The "36.8 pm/V at 1.313 µm" claim is
   unverified. Shoji's 1997 abstract also lists SHG at 1.548 and 1.533 µm for some
   material, possibly a direct GaP or GaAs value near 1.55 µm. Worth checking in the
   library copy.
2. LiNbO3 n2_I = 1.8e-19 has no primary source that I could find.
3. For liu2018 (AlN n2 = 3.5e-19), the attribution is confirmed by liu2020 and by
   Brunetta et al. 2026 (arXiv:2511.05666; not in the bib), but the method was not read.
4. For pu2016 (AlGaAs n2 = 2.6e-17), the full text was not read. aitchison1997's
   value comes via vermeulen2023.
5. For majkic2017, only the abstract-level numbers were seen. It is unknown which
   LiNbO3 d33 served as the reference.
6. pastrnak1966 page range (K5–K8) is not in Crossref.
7. bache2012 is an unpublished preprint. phillips2011 (published) uses its 38.6e-20
   value.
8. bruch2018's χ(2) convention is ambiguous.
9. The Gehrsitz code is transcribed from refractiveindex.info-scripts, not from the
   paper. It was cross-checked at x = 0 against Skauli and at x = 0.2, 1.55 µm
   against a literature value.
