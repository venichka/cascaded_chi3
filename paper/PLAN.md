# Paper plan: resonant cascaded χ⁽²⁾:χ⁽²⁾ Kerr nonlinearity in metasurfaces (and microrings)

*Created 2026-10-06. Owners: Nikita Nefedkin (theory), Sedigheh Esfahani (metasurface design, completion). Status (2026-10-21): the draft is an outline plus a filled theory part, edited for logic and language on 2026-10-21. §0 and §10 are current. §§1–9 and 11–12 are the long-form plan of 2026-10-06, kept as background; where their numbers differ from §0, §7 and the manuscript, the manuscript is right.*

**Decisions already taken (2026-10-06):**

- Format: long-form theory paper (~10–15 pp main text plus Supplement).
- Metasurface scope: the general theory plus the singly-resonant (SH-only) implementation, shown on **one** design. Sedigheh will make a new design or pick one of the existing ones.
- Microring: one section and one figure, framed as "the same theory, different applications".
- Writing is split: Nikita drafts the theory sections before 2026-10-26, and Sedigheh writes the metasurface-design sections and finishes the paper.

**Further decisions (2026-10-06):**
- **Authors:** Nefedkin, Esfahani, Khurgin, Alù; plus the designers if their design is used.
- **Venue shortlist:** Phys. Rev. Applied, ACS Photonics, or Optica.
- **Microring:** show both the demonstrated device (Lu *et al.*) and the target ring.
- **`paper/` is tracked in git** (not ignored).
- **F1 has been re-derived and double-checked:** see `checks/PREFACTOR.md`.

---

## 0. Optica first (decided 2026-10-06)

The paper goes to **Optica as a Research Article** first. The long-form structure of §5–§6 below is the fallback for Phys. Rev. Applied or ACS Photonics, and its outline is kept in `main_prapplied.tex`. The Optica outline is `main.tex` plus `supplement.tex`. Both build with the official Optica templates (`optica-article.cls`, `osa-supplemental-document.cls`, `opticajnl.bst`), which are copied into `paper/`.

### 0.1 Optica requirements (checked 2026-10-06 against opg.optica.org and the official template)

| Item | Requirement |
|---|---|
| Length | Research Article 6–8 journal pages (Letters 4, Mini-reviews 8–12). Estimate: universal-template pages × 0.6, so ~11 template pages is about 7 journal pages. |
| Abstract | About 100 words; no reference numbers; "first/new/novel" must be qualified ("to our knowledge"). The outline's abstract is 103 words. |
| Criteria | Significance, potential impact, originality; technical quality and rigour; readability and interest to the broader optics community. Incremental work is declined. |
| Cover letter | 250 words: why Optica, what is novel, why it matters to the broader optics community. |
| Authors | Exactly one corresponding author, with an e-mail. |
| Figures and tables | Captions under figures; tables have titles. Don't rely on colour alone (use line styles and markers). |
| References | Numbered in order of first citation, with article titles, up to three authors then et al. (`opticajnl.bst` does this). |
| Back matter | Funding (filled in from Prism), Acknowledgment (optional), Disclosures ("The authors declare no conflicts of interest."), Data availability (required), Supplemental document ("See Supplement 1 for supporting content."). |
| Supplement | A separate PDF (Supplement 1), peer-reviewed, not copy-edited. Items are numbered S1, S2, … The article must be coherent without it. |
| Submission | Line numbers on (`lineno`); files at the archive root. |

### 0.2 What changes relative to the long-form plan

- **Length.** About 11 template pages including figures and references, versus ~13–15 two-column pages. That means ~4,500 words, five figures, one table and ~35 references in the main text. Page budget per section:
  - Introduction 0.8;
  - Theory 1.2;
  - Limits and hosts 1.0;
  - Metasurface 2.2;
  - Microring 0.8;
  - Discussion 0.6;
  - references ~1.
- **Sections.** §§2–4 of the long form merge into two sections ("Cascaded susceptibility at a second-harmonic resonance"; "Limits and host materials"). The metasurface section keeps two subsections: doubly resonant, and singly resonant + proposed measurement. The microring is one paragraph-length section. Discussion and Conclusions merge.
- **Moved to Supplement 1:**
  - every derivation (projection, QNM kernel, uniform limit);
  - the theory-check figure (Fig. S1);
  - the sum-rule and threshold derivations;
  - the host table (Table S1);
  - simulation methods;
  - the singly-resonant overlap extraction;
  - the microring model (photon units, Bogoliubov threshold, thermal/photorefractive model, EO tuning);
  - all parameters.
- **Figures (5):**
  1. concept and lineshape;
  2. sum rule + host comparison;
  3. doubly resonant metasurface (cell and modes, R(Q₂, η), bistability, FWM);
  4. singly resonant design + simulated Z-scan for ±Δω + pulse-length dependence;
  5. microring (R vs Q₂, tilt reversal, thresholds).

  Table 1 is the platform comparison.
- **Framing for a broad readership.** Lead with one formula that covers bulk crystals, metasurfaces and microrings (the bulk limit makes it recognisable), followed by two design rules (sum rule, host rule) and two concrete platforms.
  - The metasurface is the original part.
  - The microring is presented as the CW-device limit and credits Cui *et al.* 2022/2025 explicitly.
  - No "first/new/novel".

### 0.3 Optica-specific risks

| Risk | Mitigation |
|---|---|
| Theory-only paper judged not significant enough | Make the metasurface part concrete: a full-wave design with M̃ computed from its fields, and **a simulated Z-scan signal with realistic picosecond pulses, showing it is detectable** (Fig. 4c). Without that panel the proposal is weak. |
| "Incremental over Cui *et al.*" | Keep the ring short; the novelty rests on the metasurface, the general formula with its bulk limit, the sum rule and the host comparison. |
| Length overrun | Cut the microring to one figure with three panels; keep only one table in the main text. |
| Desk rejection | The Optica Publishing Group template is shared by Optics Express and Photonics Research, so a transfer is straightforward. Otherwise PR Applied / ACS Photonics with `main_prapplied.tex`. |

### 0.4 Cover letter (250 words), points to make

1. A resonance at the second harmonic replaces the phase mismatch of bulk cascading, which gives a sign-tunable Kerr response in χ⁽²⁾ hosts, with one closed formula from bulk to resonators (overlap η ≤ 1, exact bulk limit). This extends the resonant idea of Khurgin and Ding (1994–95) from surface-emitting structures to arbitrary resonators.
2. Two rules follow for any resonator: the strength–bandwidth product does not depend on Q₂, and the gain over the native Kerr response depends on the host only through F. For LiNbO₃ the cascade exceeds the native response by ~20 in a doubly resonant metasurface (Q₂ = 4000, η = 0.3) and by 400–3000 in poled microrings, at bandwidths of 0.1 THz and below 1 GHz. In a microring the cascaded shift follows from the measured SHG efficiency with no free parameter (Eq. 7), which makes a direct test possible on existing poled rings.
3. A concrete metasurface design and a picosecond Z-scan protocol that isolates the cascade (the part odd in Δω) from the native response. *This point needs the design and its detectability budget.*
4. Relevance to metasurface nonlinear optics, integrated LiNbO₃ photonics and all-optical switching.

### 0.5 Decisions this adds

- **Material constants.** **Decided 2026-10-07: sourced 1.55-µm values** (see §7 and `checks/constants/`). The LiNbO₃ headline is R ≈ 20; the abstract says "about 20".
- **Corresponding author.** **Decided 2026-10-06: Andrea Alù** (aalu@gc.cuny.edu; set in `main.tex` and `main_prapplied.tex`). Optica requires exactly one, and Nikita leaves CUNY on 10-26.
- **Spelling.** Optica Publishing Group's style guide does not fix British versus American spelling. The outline keeps British spelling consistently, as in the notes.
- **Data availability.** The default is "available on reasonable request". Depositing the scripts (e.g. on figshare/Zenodo) and citing them would strengthen it.
- **Review of 2026-10-08 (`checks/REVIEW_2026-10-08.md`), decided by Nikita:**
  1. *Prior art is flagged.* Khurgin–Ding 1994–95 (`khurgin1994`, `lee1995`) proposed a SH resonance in place of Δk; cavity effective Kerr (Drummond, Majumdar, Fryett, Gustin) and LiNbO₃ cavities/waveguides at 10³–10⁴× (Musgrave, Nie, Li 2026) are cited. The abstract and introduction now present the work as the general resonator formulation (vector overlap, η ≤ 1, bulk limit, F, strength–bandwidth relation) and its metasurface implementation. A `% FLAG (positioning)` comment marks the spot in the intro.
  2. *Venue, the concrete design, absolute figures of merit and the Z-scan SNR budget are Sedigheh's* (review items A2 and A5; see §10).
  3. *Fig. 3(c) moved to γ₁ = γ₂* (was 10γ₂). The SH then stores 0.5–0.8× the pump energy at the folds, and the native cross terms are ≤ 0.16|Δω|. Validity condition (main §2, supplement §4.B): cascaded shift ≪ |Δω|(R/2)^{1/2}.
  4. *Microring R is quoted as a range:* 400–3000 in the abstract, 6×10²–2.8×10³ (A) and 4×10²–1.8×10³ (B) in §5, from g taken from η_SHG versus the fitted lineshape.
- **Status:** the paper is a draft with an outline and a filled theory part; the design, Fig. 1, the introduction prose and Section 6 are open.

### 0.6 Story check for Optica (2026-10-21)

**The story in four steps.** The introduction (item 4) and the discussion follow it, and the cover letter should too.

1. In a bulk crystal the cascaded Kerr response is set by the material phase mismatch, cΔk/n₂ ≈ 0.04ω₁ in LiNbO₃. A detuned SH resonance replaces it with a designed detuning, Δω ≈ 10⁻³ω₁ at Q₂ = 4000 and Δω = 2γ₂.
2. One closed formula (Eq. 3) covers bulk crystals, metasurfaces and rings. Its overlap η ≤ 1 equals 1 in the bulk limit, so it is the fraction of the bulk bound-wave coupling that a resonator keeps.
3. Two rules follow for any resonator. The product of strength and bandwidth is fixed (Eq. 5), and the gain over the native Kerr response is R ≈ Q₂ηF/6, where the host enters only through F (0.1 for LiNbO₃, about 4×10⁻³ for the III–V hosts).
4. The platforms are different points of the same trade-off (Table 1):
   - metasurface: R ≈ 20 at 0.1 THz;
   - microring: R ≈ 400–3000 below 1 GHz;
   - bulk: R ≈ 2 above 10 THz.

   The metasurface is the free-space, sub-micrometre, picosecond option, and a Z-scan of the part odd in Δω isolates its cascade.

**Verdict.** The logic from §2 to §6 holds: formula, then the bound and F, then the two platforms, then the comparison. The weak point for Optica is the gap between the title and the evidence.
- **The metasurface results are the least concrete.** The title and abstract lead with metasurfaces, but the doubly resonant result rests on an assumed Q₂ = 4000 and η = 0.3, and the singly resonant case is a placeholder. The ring has the larger numbers and the only CW device, and its physics is prior art.
- **The metasurface does not win on size.** |χ_casc| ≈ 6×10⁻²⁰ m²/V² is below the native χ⁽³⁾ of Si or GaP, and bulk cascading has the larger |χ|Δf. What it offers is a Kerr sign set by Δω (tunable in situ), in a sub-µm free-space layer at picosecond speed, plus a parity measurement that isolates the cascade. The text now argues this (intro item 4, §6), but the argument needs a design and a detectability number to convince a referee.
- **Optica therefore depends on two deliverables.**
  1. One real design with M̃ computed from its fields. It fills Fig. 4(a), §4.2 and row 4 of Table 1.
  2. The Z-scan signal of that design in absolute units: phase per GW/cm², substrate background and SNR.

  With both, Optica is a reasonable attempt. Without them this is a theory paper with idealised metasurface numbers, better sent to Optics Express or PR Applied; the text carries over unchanged.

**Recommendations, not applied (decide before submission):**
- **Title (D3).** Keep the metasurface title only if the design arrives. Otherwise use a general one, e.g. *Cascaded Kerr nonlinearity at a second-harmonic resonance: from metasurfaces to microrings*.
- **Positioning.** Confirm the wording with J. Khurgin, since his 1994–95 papers are the direct precedent for replacing Δk by a resonance (intro item 2, `% FLAG`).
- ~~**li2026.** Check its mechanism.~~ Done 2026-11-09, see §0.7.
- **Measurement.** Weigh the alternative measurement (the pump-resonance shift versus Δω in a doubly resonant sample), which the 10-08 review rated more sensitive than the singly resonant Z-scan.

### 0.7 Literature check of the open items (2026-11-09)

Five parallel searches, with full texts where accessible. All new `refs.bib` entries were checked against Crossref or the arXiv API. Main, supplement and the fallback build clean (main 11 pp, supplement 12 pp).

**Closed and applied:**
- **li2026.** arXiv:2512.10462 full text: the mechanism is cascaded SHG–DFG, analysed as phase matched, in a straight 1-cm PPLN TFLN waveguide with CW pumps. n₂,eff is fitted from the FWM idler efficiency with an assumed seed, Pᵢ(0) = 10⁻⁸ W, and no sign or Δk dependence is reported. It is not an R comparable to ours.
  - Table 1 row replaced by the measured AlN microring of cui2025: SPM tuned −2.7 to +4.7× intrinsic, so R = 3.7; Q₂ = 1.8×10⁵; pump linewidth 0.33 GHz. li2026 is kept in intro item 3 with the scoped wording.
  - Fallback `main_prapplied.tex` fixed the same way.
- **musgrave2026 / nie2026** are free-space bulk-PPLN cavities, not microresonators. eLight is doubly resonant, with the Lorentzian SH-detuning form; Nat. Commun. is singly resonant, with the sinc-Δk form. Their ">10³" and ">100×" are model ratios. Intro item 3 now says "estimated at 10² to above 10³ times … drives dissipative quadratic solitons".
- **cui2022 / cui2025** used unpoled AlN rings, with g calculated from an assumed χ⁽²⁾ ≈ 1 pm/V. Cui 2022 measured the FWM enhancement, not a shift; Cui 2025 measured the SPM shift versus Δ, with no SHG efficiency reported. **No experiment has yet measured both η_SHG and the cascaded shift on one device**, which is the opening for the ring section (see the microring discussion).
- **Positioning against Khurgin–Ding.**
  - What they have (1997 JOSA B 14, 1977 full text; 1994–97 abstracts): the SH-only resonance taken exactly on resonance, with χ⁽²⁺²⁾ ∝ Q_i, scalar unbounded overlaps, and "the enhancement … is equal to the Q factor of the intermediate wave". There is no detuning dependence or sign change, η ≤ 1 bound, closed bulk limit, sum rule or F.
  - Applied: intro item 2 now adds "with an enhancement … of the order of the quality factor of that resonance [khurgin1997]", and the FLAG comment lists what is and is not in their papers.
  - Applied: `khurgin2023` (Nat. Photon. 17, 545, "Nonlinear optics from the viewpoint of interaction time": resonant enhancement limited by bandwidth) is cited next to Eq. (5).
  - Still open: confirm the wording with Khurgin; the 1994/95 formulas are paywalled.
- **Scoop check (2024 to Nov 2026): no direct scoop** of the metasurface cascaded Kerr, the sum rule or the odd-in-Δω Z-scan. Newly cited: `liu2019ol` (OL 44, 1456: cascaded effective FWM in an LNOI microdisk, the earliest LN microcavity case) and `musgrave2026opo` (arXiv:2604.03474, sign-tunable effective Kerr above 10³ in a doubly resonant OPO).
  - Candidates not added:
    - Son & Chekhova, arXiv:2601.23137 (in a bare thin χ⁽²⁾ layer, direct χ⁽³⁾ SFWM beats the cascade; supports our gap);
    - Panov, OL 50, 5849 (2025) (GaP metasurface effective n₂ about 700× bulk, native);
    - Skryabin, JOSA B 37, 2604 (2020) (CMT for QPM χ⁽²⁾ rings);
    - Talenti 2025, JEOS-RP;
    - Gutierrez/Marandi, arXiv:2610.30079.
  - Coverage of 2026 Optica and Nature-family journals was weak (bot walls). Rerun a manual search before submission.
- **LiNbO₃ n₂ᴵ now has primary sources.**
  - He et al., Optica 6, 1138 (2019) inferred 1.8e-19 from an OPO threshold, ordinary polarisation.
  - Shams-Ansari et al., APL Photon. 7, 081301 (2022) measured 1.74e-19 (extraordinary) and 1.61e-19 (ordinary) by XPM at 1550 nm, with the cascade estimated below 10% and not subtracted.
  - Supplement §5 cites both. The adopted 1.8e-19 and every number stay.
  - Details are in `checks/constants/constants.md`.
- **Q₂ limits of the SH BIC (main §4.1, literature text).**
  - Q_dis ∝ σ⁻² (kuhne2021, maslova2021) and Q_fin ∝ N² (bulgakov2017).
  - Measured Q: ~10⁴ for Si₃N₄ near 740 nm (hsu2013), 800 for the GaN SH mode at 775 nm (wang2020), >8×10³ for LNOI metasurfaces at 1.55 µm (difrancescantonio2025). Q₂ = 4000 at 775 nm has not been reported in LN.
  - **New design constraint:** a Γ-point BIC above the diffraction threshold is not symmetry protected (cerjan2021). On silica the SH needs p < 535 nm, where a first-order Γ-point pump mode would need n_eff > 2.9, out of reach in LN. The doubly resonant LN design therefore needs a period-doubled pump mode (weak SH leakage ∝ δ²) or a membrane (p < 775 nm). **The singly resonant design is not affected**, which supports putting the measurement there. Stated in §4.1, with a `FOR THE DESIGN` note for Sedigheh.

**Microring section reframed (2026-11-09, user decision): the ring is now the place where the theory can be tested quantitatively.**
- **Main §5 rewritten (about 380 words).** It leads with the new **Eq. (7)**: δω_casc/κ₁ = η_SHG P_in κ₂Δω/(Δω² + κ₂²) at critical coupling and locked pump, to first order. Other couplings multiply it by (κ₁/2κ₁ᵉ)(κ₂/2κ₂ᵉ).
  - Eq. (7) has no overlap, material constant or mode volume. It is ρ written with the measured SH power.
  - Lu 2020 device: 0.09 half-widths = 17 MHz at 5 µW and Δω = 2γ₂, against 0.03 MHz native. That is below its parametric and photorefractive limits. The fitted g gives 4.7× more.
  - Claim: "to our knowledge, no device has been characterised by both η_SHG and the cascaded shift". The odd-in-Δω EO-dither extraction is tied to the §4.2 Z-scan.
  - Kept: R ranges and the parametric threshold 8 µW (vs 9 mW native, below 37 µW bistability).
  - Moved out of the main text (still in supplement §8): the response-time 4%, 8 fJ switching and the thermal numbers.
- **Supplement §8, new subsection "Cascaded shift from measured quantities".**
  - Exact relation δω_casc = (Δω/κ₂)·P_SH^tot/(2W₁), valid for any resonator within the elimination.
  - Eq. (S-shift-eta) with the general coupling factors. Lu device: coupling ratios 0.70 / 0.22, factor 1.59.
  - The locked-pump nonlinear-loss correction is 0.96 at 2γ₂ and 0.76 at |Δω| = κ₂.
  - All numbers checked numerically, including the energy form, which reproduces −U n exactly.
- **Fig. 5 now has two panels** (`figures/fig5_ring.py`):
  - (a) R vs Q₂ with a band between the two g values for design B, unpoled C, the Lu device at both g, and the AlN point of cui2025 (3.7 at |Δω| = γ₂/2, scaled by 8/17 to 1.7 at 2γ₂);
  - (b) predicted shift vs Δω for the Lu device at 5 µW, from Eq. (7) and from the locked-pump steady state.
  - The old threshold-vs-Q₂ panel was dropped; its numbers stay in Table S2.
- **Table 1.** The last column is now "CW threshold": bistability / bistability / none / – / "8 µW, parametric". The AlN R is 1.7 at 2γ₂, with a footnote giving the measured range. The footnote also records White's pump-resonant cavity and the 37 µW bistability.
- **Intro item 4** adds: "In a microring the cascaded shift follows from the measured second-harmonic efficiency, so the formula can be tested with no free parameter."
- Main is now 12 pp, supplement 12 pp; both build clean.

**Final check before Overleaf and the handoff (2026-11-09, later).** Four independent reviews (theory analytics and numerics; microring; claims, citations and consistency; language) plus mechanical checks: every `refs.bib` DOI verified against Crossref (title, volume, first page, year, first author), all `\cite` keys present, writing-rule grep clean. Everything below is applied; main 12 pp, supplement 12 pp, fallback 9 pp, all clean.
- **Numerics corrected:** pulsed Z-scan window "within 7%" → "within 8%" (the reduction factor peaks at 1.073 at γ₂τ = 2.25); "overestimates R by 6%" → "by one part in 16"; ring modulation bandwidth at 20 µW "8%" → "7%"; Table 1 ring row R = 4×10²–1.8×10³ (design B; 2.8×10³ belongs to design A); Table S2 footnote and §8.E cross-references moved to 8.E/8.F after the new subsection.
- **Analytics sharpened:** Eq. (7) holds "when the nonlinear loss is neglected" (exact in the shift for a locked pump); its Lu-device numbers include the coupling factor 1.6, now stated; Eq. (7) described as the SH power measuring the loss and ρ converting it; "a further fraction equal to the cascaded shift divided by |Δω|"; the quasi-static bound sentence no longer calls Δf ≤ |Δω+iγ₂/2|/2π "the quasi-static range"; the bulk strength–bandwidth comparison says Eq. (5) does not bound the bulk response (no single pole; the η = 1 bound would be 3.8e-8 against the quoted >5e-8).
- **Claims scoped:** "thermal and photorefractive shifts comparable with or larger than" (design A at 5 µW: 0.6κ₁ vs 0.09κ₁); "we have found no report of Q₂ = 4000 at 775 nm in LiNbO₃"; N² scaling marked as one-dimensional; Cui 2025 extremes "at the extrema of the tuning curve, |Δω| = γ₂/2"; Fig. 4(b) "an illustrative value R = 3"; gorbounova1996 "appear in" (bandwidth trade-off only implicit in its abstract). Verified and kept: White 1996 is monolithic; Lu 2023 reports neither η_SHG nor a measured shift, so "no device characterised by both" stands; khurgin1997 does state enhancement ≈ Q.
- **Consistency:** Table 1 bandwidth footnote names which resonance (SH for the metasurface, pump for the rings); Discussion says the SH linewidth sets the response time; "intrinsic" → "native" for the Kerr response throughout; Table 1 caption in `\bf` like the supplement tables; Fig. 3(c) result sentence moved from the caption into §4.1; Fig. 3(d) "for Δω < 0"; Fig. S1 drive given as |s|² = 10⁻⁴γ₂³/|β|²; `fig3_double.py` docstring fixed.
- **Language:** abstract 122 → 113 words; ~60 line edits (long sentences split, Oxford commas removed, duplicate intro phrases cut, "Here we" dropped, measurements put in the conditional, supplement symbols Q₂/δ₁/R/D₂ defined at first use, supplement equation references as "Eq. (S\,n)", the standard g formula displayed, subsection "Electro-optic and thermal tuning of Δω").
- **refs.bib:** pastrnak1966 → Phys. Status Solidi B; Nat. Photonics and APL Photonics spelled as in the other entries. fathi2026 has no volume yet (Crossref, 2026-09-02).
- **Left as found (the user's own edits since 10-21):** the corresponding author is now Nefedkin (nnefedkin@gc.cuny.edu) with a comment to discuss the order with Alù and Khurgin; Data availability reads "Write here something."; Funding is empty; "See Supplement for supporting content." These must be settled before Sedigheh sees the file.

**Still open after the check:**
- **GaP d₁₄.** Search-engine extracts of Shoji 1997 Table 11 suggest 36.8 pm/V at 1.313 µm, which gives 34 pm/V at 1.55 µm and F 2.8× lower. Unverified; needs a library copy (Sedigheh has CUNY access). The paper's "upper estimate" wording is safe either way.
- **Positioning wording:** confirm with Khurgin.
- **Funding statement** (grant numbers from Alù), **title** (D3), **venue** and **data deposit**: decisions, not research.

## 1. The paper in one paragraph

Cascading two χ⁽²⁾ processes (SHG, then DFG) produces an effective Kerr response on the pump. In bulk crystals the strength of this cascade is set by the phase mismatch, χ³_casc ∝ 1/Δk. **We replace the phase mismatch with a resonance:** a dark (bound-state-in-the-continuum) second-harmonic mode detuned by Δω = 2ω₁ − ω₂ gives χ³_casc ∝ Q₂/Δω. We derive a closed, normalisation-invariant expression for χ³_casc in an arbitrary resonator from two independent routes. One is the χ⁽²⁾·G·χ⁽²⁾ Green's-function route; the other is temporal coupled-mode theory with exact elimination of the SH. In the uniform-medium limit the expression reduces exactly to the textbook bulk result. The resonant cascade is sign-tunable through Δω and comes with a fixed loss partner, ρ = γ₂/(2|Δω|). It obeys a bandwidth–enhancement sum rule, so a resonance redistributes cascaded nonlinear action rather than creating it, which places it in the narrowband/CW niche. In LiNbO₃, whose native Kerr is weak, a doubly-resonant BIC metasurface beats the native χ⁽³⁾ by R ≈ 20 (11–30 across the reported constants) at a conservative operating point *(updated 2026-10-21; was 50 with the pre-10-07 constants)*. A singly-resonant design gives a smaller R, to be computed for the chosen design, and its fingerprint is that the nonlinear refraction reverses sign as Δω crosses zero. For a periodically poled thin-film LN microring, the same theory reduces to the cavity expression already used to demonstrate sign-tunable cascaded self-phase modulation in AlN rings (Cui, Zhang & Fan 2022, 2025). Applied there quantitatively, it gives R ≈ 400–3000 and microwatt-level CW bistability, switching and parametric thresholds. Δω can be set electro-optically, and a measurement antisymmetric in Δω separates the cascade from thermal and photorefractive shifts.

**Working titles** (decide later, D3):

1. *Resonant cascaded Kerr nonlinearity in χ⁽²⁾ metasurfaces*
2. *Trading bandwidth for Kerr: sign-tunable cascaded χ⁽³⁾ from dark second-harmonic resonances*
3. *Cascaded χ⁽²⁾:χ⁽²⁾ nonlinearity at bound states in the continuum: theory, limits, and implementations*

---

## 2. Key messages (what a reader must leave with)

| # | Message | Main evidence |
|---|---|---|
| M1 | **Resonant lever.** χ³_casc = −(2ε₀ω₁/3)·\|M\|²/[N₂(Δω + iγ₂/2)∫\|f₁\|⁴] is exact within single-mode theory, has two independent derivations, and reduces to the bulk bound-wave cascade in the uniform limit. | Theory §2; Fig. 2 |
| M2 | **Built-in limits.** The response is sign-tunable (odd in Δω) and its loss partner is even in Δω. The sum rule says \|χ³_casc\|·Δf is independent of Q₂. A parametric-instability ceiling sets in at \|βa₁\| ~ \|Δω\|. The niche is narrowband and CW, not ultrafast. | §3; Fig. 2c |
| M3 | **Host selection.** R ∝ Q₂·η·\|χ²\|²/(n₂²χ³_nat), so χ⁽²⁾ hosts with weak native Kerr win (LiNbO₃, F = 0.1). III–V hosts, with their band-edge Kerr, reach R = 1 only at Q₂ ≈ 4–8×10³ (F ≈ 4×10⁻³). | §4; Table 1; Fig. 3 |
| M4 | **Metasurfaces (centre of mass).** A dark SH BIC is the natural resonant lever. The doubly-resonant device gives device-level leverage, (1+R) on phase and (1+R)² on FWM. The singly-resonant metasurface is the fabricable first step, and its sign flip can be measured by a picosecond Z-scan versus pump detuning, the resonant analogue of DeSalvo 1992. The signature must be defined as the part that is **odd in Δω = 2ω₁ − ω₂ at fixed FW detuning**: native χ⁽³⁾ alone can invert the effective n₂ near a Mie resonance (Panov 2020). | §5; Figs. 4–5 |
| M5 | **Microrings: the same theory, quantified.** Sign-tunable cascaded Kerr in χ⁽²⁾ rings is established experimentally (Cui *et al.* 2022, 2025, AlN; Lu *et al.* 2023, PPLN). Our contribution is quantitative and protocol-level, for poled TFLN: R ≈ 400–3000 against the same ring's native Kerr (≈13 unpoled); thresholds evaluated at Δ_eff with multimode Bogoliubov stability; electro-optic Δω tuning; and an antisymmetric-in-Δω lock-in, valid only at δ₁ = 0, that rejects thermal and photorefractive backgrounds. | §6; Fig. 6 |

---

## 3. What the analysis found before any writing (read this first)

### F1 (critical): the boxed χ³_casc in both notes is 4× too small

**Where it comes from.** In note §5 the projection of the wave equation onto f₂* reads −2iω₂ȧ₂/ε₀ = k₂²Ma₁², which gives β = ε₀ω₂M/2 = ε₀ω₁M. The note writes ε₀ω₂/4 at that step, an arithmetic slip, and then "adopts" ε₀ω₁/4 as a convention, which halves it again. α₃ = 3ε₀ω₁/8 carries the same factor-4 deficit. The result scales as χ³_casc ∝ \|β\|²/α₃, so it comes out 4× too small. The SH-only note uses the same β̂ = ε₀ω₁/4 (lines 165–174), so χ³_casc^(1R) is also 4× too small. There the factor comes through β̂ alone, because the 1R back-projection has no α₃.

**Correct, normalisation-invariant form** (Boyd amplitudes, N_n ≡ ∫ε₀ε̃_n\|f_n\|²d³r):

- Couplings: β = ε₀ω₁M/N₂ and α₃ = (3ε₀ω₁/2)·∫χ³f₁*f₁f₁f₁*/N₁.
- Doubly resonant:
  ```
  χ³_casc = −(2ε₀ω₁/3) · |M|² / [ N₂ (Δω + iγ₂/2) ∫|f₁|⁴ ]
          ≈ −(2ω₁ / 3n₂²Δω) · |χ²|² · η_geom,
  η_geom ≡ ε₀ n₂² |M|² / (|χ²|² N₂ ∫|f₁|⁴)
  ```
  - η_geom equals 1 for phase-matched uniform plane waves, so it reads as *the fraction of the bulk bound-wave cascade the resonator recovers*.
  - At the conservative point Δω = 2γ₂: \|χ³_casc\| = Q₂\|χ²\|²η/(6n₂²). The old value was Q₂/(24n₂²).
- Singly resonant: χ³_casc^(1R) = −(2ε₀ω₁/3)·\|M̃\|²/[N₂(Δω + iγ₂/2)V_cell].

**How it was verified (2026-10-06; full write-up in `paper/checks/PREFACTOR.md`):**

1. **Two independent analytic derivations.** One is mine. The other is blind: a separate agent was given only the conventions, and it reproduced every coefficient.
2. **The uniform-medium limit.** The corrected formula gives χ³ = −2χ²²/[3(n₂² − n₁²)], the exact bulk bound wave and the note's own eq. `chi3casc-bulk`. The note's boxed form gives 1/4 of it.
3. **Five 1-D tests in a non-uniform cell** (`paper/checks/prefactor_check.py`, exact finite-difference solves). At Δω/2ω₁ = 10⁻³, exact/corrected and exact/note are:

   | Test | exact/corrected | exact/note |
   |---|---|---|
   | β | 1.0005 | 4.00 |
   | Full doubly-resonant χ³ chain | 1.0004 | 4.00 |
   | SH-only χ³ chain | 1.005 | 4.02 |
   | α₃, from an exact nonlinear eigen-solve | 1.0002 | 4.00 |

4. **A time-domain two-mode Maxwell integration** without the slowly-varying approximation, from the blind derivation: 0.9995–1.0016.
5. **External checks.** The corrected chain reproduces the standard microresonator coupling g = √(ħω_a²ω_b/8ε₀)·(overlap) (Zhao *et al.*, arXiv:2105.12705, Eq. 1). It also reproduces the microring note's independently derived native U_nat. The note's chain misses both, by 2√2 and by 2.

**What changes:**

| Quantity | Old | Corrected |
|---|---|---|
| R, LiNbO₃ design point (Q₂ = 4000, η = 0.3, Δω = 2γ₂) | 12 | **≈ 50** |
| R at Δω = γ₂ (edge of adiabatic) / γ₂/2 | 25 / 49 | ≈ 100 / ≈ 200 |
| Q₂^cross (LN) | ≈ 310 | ≈ 80 |
| Host survey R at Q₂ = 4000: GaP / AlN / GaAs / AlGaAs | 0.54 / 0.57 / 0.36 / 0.05 | 2.2 / 2.3 / 1.4 / 0.2 |
| SH-only R⁽¹ᴿ⁾ at C_bridge = 1: S1 / S2 / S3 / S4 / D6 | 0.49 / 0.39 / 1.76 / 2.49 / 1.81 | 2.0 / 1.6 / 7.0 / 10.0 / 7.2 |
| Bulk PPLN advantage in \|χ³\|·Δf | ≈ 40× | ≈ 10× |
| Microring η_geom inferred from measured g (poled) | 0.93 ("at the physical ceiling") | **0.23** (the anomaly disappears) |

**What does not change:**

- Microring R, thresholds, n₂_eff and switching numbers. They are anchored on the measured g, not on β.
- Every Δω-dependence and lineshape, the sign flip, ρ = γ₂/(2\|Δω\|), the sum-rule *shape*, the Q-scaling exponents and the observable leverage factors expressed in R (1+R, (1+R)²).
- Dynamic range in decades (onset and saturation both scale as 1/\|β\|²).

**What must be recomputed:** every absolute power, energy or phase in physical units on the metasurface side (thresholds in W, φ_NL in rad, saturation energy). F2 shows the note's |a|²-to-joule mapping is itself off by a factor of 2.

**Propagation list** (code hard-codes the old prefactor):

- `cascaded_chi3/sample_params.py:204`
- `design6_effects.py:103,123`
- `fig_cascaded_vs_native.py:43,48`
- `fig_bulk_vs_metasurface.py:52`
- `fig_overlap_eta.py:41`
- `fig_n2_scaling.py:39`
- the note's §5 (β, α₃), §6 (boxed), §10 (numbers; *delete* the "Prefactor uncertainty" paragraph), §11 table, §14
- SH-only note §3, §5, §10, Appendix
- simulation report; sample report; `data/design6_parameters.md`
- microring note §4 (g-estimate paragraph) and caveat 8

The MURI deck needs no change: it is delivered, and its numbers carried an explicit "O(1) prefactor" disclosure.

**Why this matters for the paper:** the main weakness of the earlier notes (an admitted factor-4 uncertainty) becomes a strength. The paper can state an exact prefactor and validate it against the bulk limit, and the metasurface case gets 4× stronger. **Status: confirmed and propagated (2026-10-06).** Updated:
- the main note, the SH-only note and the compact derivation;
- the simulation and samples reports and the microring note;
- the Markdown records;
- the code and the regenerated figures.

All documents recompile with no errors or undefined references.

Side fixes made along the way:
- the report's operating-point table now uses the full Lorentzian (its Δω = γ₂/2 row had overstated R by 2×);
- a dimensional slip in the SH-only sensor slope;
- the γ₃ factor-2 mismatch;
- a pre-existing crash in `ring_effects.summary()`.

The 2026-08 MURI deck and the SH-only beamer decks were left as delivered. If recompiled they would pick up the new figures while their text still quotes the old R values.

### F2: the field convention is stated inconsistently

The note says real fields are ½(E e^{−iωt} + c.c.). Yet it uses Boyd's polarisation relations (SHG with no prefactor, DFG ×2, SPM ×3), which belong to Boyd's convention, Ẽ = E e^{−iωt} + c.c. (no ½). χ³_casc is unaffected, since it is a ratio of material constants and is convention-free when each convention is used consistently. The *energy* interpretation of \|a\|² is affected, and so is every quantity quoted in joules or watts.

**Fix (done 2026-10-06):** the notes now state Boyd's convention, normalise with ∫[ε₀ε̃|f|² + μ₀|h|²] = 1 (so \|a\|² = U and N_n = ½), and give the principal results with N₂ explicit. The SH-only local-field relation now uses Boyd's I = 2nε₀c\|E\|². That cancels the 4× in χ³ for the single-pass phase, loss and FWM, which are therefore unchanged, while η_SHG doubles.

### F3: the disorder model does not fit a collective BIC

Note §13 item 7 models disorder as inhomogeneous broadening of independent elements, γ_array = √(γ₂² + σ²). A BIC or qBIC is a collective lattice mode, so disorder mainly opens radiation channels (Q_dis ∝ 1/σ²), and finite array size limits Q as well. Replace the model with one based on the BIC-disorder and finite-size literature (to be sourced), or with full-wave supercell statistics from Sedigheh's design. This is the main fabrication argument for the metasurface, so it needs to be right.

### F4: smaller inconsistencies to clean before reusing text (all fixed 2026-10-06)

**Main note:**

- §6(ii) defines γ₃^eff = \|β\|²γ₂/(2Δω²), while the §11 caveat uses \|β\|²γ₂/Δω². This is an amplitude-versus-energy-rate factor of 2; pick one.
- The §14 closing sentence says "can outperform bulk by 1–2 orders of magnitude" without naming the metric. It holds for CW threshold power, not for \|χ³\|·Δf.
- §1's literature-audit sentence is replaced by the prior-art section below.

**Microring note** (from the 2026-10-06 read):

- Summary item 3 (`:1357`) still claims the measured g "favour[s] the larger branch". After F1, rewrite §4's scalar estimate with the consistent chain (g = ½ε₀M√(ħω₁³), the standard formula). It gives g/2π ≈ 1.6 MHz at Θ = 1 for Lu's device, against 1.78 MHz from Lu's own calculation and 1.2 MHz measured, i.e. Θ ≈ 0.75. That is consistent. The "factor 2.1 below" paragraph and caveat 8 go, and the inferred η_geom falls from 0.93 to 0.23.
- §1 (`:118`) still says Q₂ is reached "routinely".
- Step 5 (`:1174`) and the `fig:thresholds` caption (`:833`) quote 3.7γ₂ as the optimum. The cold-cavity optimum is 3.3γ₂ (the figure annotates 3.4γ₂), and the caption's "divergence" applies only to the frozen model.
- The `ring_bistability` legend says "3.7γ₂ (optimum)" (`ring_effects.py:590`).
- Step 3 quotes dΔω/dT ≈ 0.9 GHz/K, while §7.1 says −0.34 to +1.6 GHz/K.
- Caveat 5 says "comparable to or larger"; §8 says "within ×2 either way".

### What is already solid and can be reused as is

| Result | Status |
|---|---|
| Closed form agrees with full two-mode TCMT integration | Validated (simulation report §2.1). Both share β, so this is prefactor-blind, which is why F1 went unnoticed. |
| Manley–Rowe (equal SHG/DFG coupling; the DFG factor 2) | Derived and checked numerically |
| Sum rule (area and peak × FWHM conserved) | Derived and checked numerically |
| Bistability, saturation, parametric onset, pulse-bandwidth limit | TCMT-validated |
| SH-only pulsed response (quasi-static for τ_p ≳ 2/γ₂; sign flip survives) | Time-domain TCMT, cross-checked |
| Microring: exact CW cubic, 4×4 stability, multimode Bogoliubov MI threshold, antisymmetric extraction under δ₁ = 0 | Validated; reviewed 2026-09-10 |

---

## 4. Prior art and positioning

*Search run 2026-10-06 across arXiv, OpenAlex, Crossref and publisher pages; every entry below carries a verified DOI or arXiv ID. The two Cui et al. papers were re-checked independently. Unverified items are listed at the end of 4.1.*

**Bottom line:**

- **Metasurface side: open.** No paper derives or measures an effective χ⁽³⁾, n₂, nonlinear phase, Z-scan, bistability or sign flip from cascaded χ⁽²⁾ in a resonant metasurface, Mie resonator or nanoresonator.
- **Microring side: core physics published.** Sign-tunable cascaded SPM, with g anchored on the measured SHG efficiency, is in Cui, Zhang & Fan (Nat. Phys. 2022; PRL 2025). Related LN-cavity results come from the Huang group (2025–26).
- **Consequence:** this confirms your choice to make the metasurface the centre of mass, and it means §6 is framed as quantification plus a measurement protocol, not as discovery.

**Background that was known before this search** (verified in earlier passes):

- Ostrovskii 1967 (prediction).
- DeSalvo *et al.* 1992 (Z-scan sign reversal through Δk = 0 in KTP).
- Stegeman, Hagan & Torner 1996 (review).
- Khurgin, Obeidat, Lee & Ding 1997 (microscopic picture).
- Bache *et al.* 2007 (scaling laws).
- Gennaro *et al.*, ACS Photonics 9, 1026 (2022). Metasurface cascading for *harmonic generation* (THG via SHG + SFG), not a Kerr response on the pump.
- Lu *et al.*, Optica 2019 / 2020 (PPLN rings; resonant anharmonicity).
- Tang *et al.*, Opt. Lett. 49, 2449 (2024) (cascaded quadratic combs in thin-film LN).
- Doubly-resonant or high-Q SHG metasurfaces: Carletti 2018, Koshelev 2020, Anthur 2020.

### 4.1 Verdict per claim and how to position it

| Claim | Verdict | Decisive prior work | Positioning in the paper |
|---|---|---|---|
| C1. Closed-form χ³_casc (TCMT + χ²·G·χ² routes) | **Partially anticipated** (cavity form) | Cui *et al.* 2022/2025: g_c = −i\|g₂\|²/(iΔ + γ/2), real part odd, imaginary part an even loss. Gustin *et al.* 2025: the two-photon loss from SH elimination. Hagmann *et al.* 2008 and Kolleck 2004: effective χ⁽³⁾ in structured and local-field settings. Tu *et al.* 2024: back-conversion in a qBIC metasurface, but only as depletion. | Claim the **first metasurface/QNM formulation**: tensorial vector overlap, normalisation-invariant, two routes, exact bulk limit (F1 check). Show that it reduces to Cui's cavity expression. Justify neglecting the optical-rectification + EO cascade channel (Bosshard *et al.* 1995): it is not resonantly enhanced, so it is a small constant background. |
| C2. Bandwidth–enhancement sum rule | **Concept anticipated; the specific rule appears new** | Khurgin's general strength–speed trade-off (SPIE 2024 names BICs; Khurgin & Sun 2013). Gorbounova, Ding, Khurgin *et al.* 1996: an early resonant-cascade power–bandwidth trade. Shcherbakov *et al.* 2019. Bache *et al.* 2007 (bulk side). | Present it as a **quantitative, cascade-specific instance** of Khurgin's trade-off. Claim the explicit sum rule and the result that resonance cannot beat bulk phase-mismatched cascading in \|χ³\|·Δf, not the principle itself. Khurgin is a co-author, so this framing is natural. |
| C3. Cascade-versus-native host survey | **Partially anticipated; the cross-host metasurface comparison appears new** | Bulk comparisons: Bache, Moses & Wise 2007; Bache & Wise 2010 (LN). Device benchmarks: Cui (AlN, ~6×), Nie *et al.* 2026 (~100×, *check by hand*), Musgrave *et al.* 2026 (>10³×), Li *et al.* 2026 (PPLN waveguide, 1.6×10⁴×). Native-χ⁽³⁾ metasurface enhancement: Panov 2020/2023/2025. | Cite the benchmarks. Add engineered χ⁽²⁾ hosts (Ramesh *et al.* 2026, coupled QWs up to 2750 pm/V; Fathi *et al.* Nat. Nanotechnol. 2026, MQW metasurface) as a caveat to "LN is favoured", since R ∝ \|χ²\|². This is a natural outlook for the team's own CQW work. |
| C4. Sign flip, bistability, ceiling as signatures | **Anticipated in bulk and cavities; new for metasurfaces** | DeSalvo 1992 (bulk Z-scan sign flip); Cui 2025 (ring sign flip versus SH detuning); White, Mlynek & Schiller 1996 (cascaded-χ⁽²⁾ cavity bistability); Drummond, McNeil & Walls 1980 and Hashemi *et al.* 2009 (instabilities); Stegeman 1996 review. **Confound: Panov 2020**, where native χ⁽³⁾ alone inverts the sign of n₂,eff near a Mie resonance. | Define the observable as the part **odd in Δω at fixed FW detuning** (tune ω₂, or fit antisymmetrically). In the 1R design the FW is non-resonant: transmission is flat (0.98 for D6) across the ~1 nm scan, so the native contribution is a constant offset. Say this explicitly in §5.3. |
| C5. Singly-resonant (SH-only) metasurface | **Not found: new as presented** | Nearest: cavities resonant at the FW only (White 1996; Nie 2026), the opposite limit. SH-BIC designs used only for SHG (Minkov, Gerace & Fan 2019; Wang *et al.* 2020). | Present it as the SH-only limit of C1 and cite those works. |
| C6. PPLN microring | **Substantially anticipated in core physics; specifics new** | Cui 2022/2025 (sign-tunable SPM −2.7 to +4.7×, Kerr suppression, Fano-shaped odd shift with even broadening, g from SHG efficiency). Lu *et al.* 2023 (χ⁽²⁾-dominated tilt in a PPLN TFLN ring). Musgrave *et al.* 2026 (>10³×, sign via detuning). Li *et al.* 2026. µW OPO thresholds (Lu 2021; McKenna 2022). | **Drop any "first sign-tunable/cascaded Kerr" wording; cite Cui prominently.** Contribution: quantified R for poled versus unpoled TFLN against the same ring's native Kerr; thresholds at Δ_eff; multimode Bogoliubov MI; EO Δω tuning; the antisymmetric lock-in with its δ₁ = 0 condition. Keep §6 short. |

### 4.2 References to cite (beyond the known background above)

- **Metasurfaces and nanoresonators:**
  - Gennaro, Sarma & Brener, Nanotechnology 33, 402001 (2022) (review).
  - Celebrano *et al.*, Nano Lett. 19, 7013 (2019) (cascaded THG).
  - McLaughlin *et al.*, JOSA B 39, 1853 (2022) (GaP microdisk, SH-detuning dependence).
  - **Tu *et al.*, PRA 109, 063522 (2024)** (nearest metasurface paper).
  - Panov, Ann. Phys. 532, 1900574 (2020) (confound) and OL 50, 5849 (2025).
  - Yang *et al.*, Nano Lett. 24, 1679 (2024) (native ultrafast qBIC switching benchmark).
  - Son & Chekhova, arXiv:2601.23137 (thin-film baseline).
- **Doubly-resonant SHG:**
  - Minkov, Gerace & Fan, Optica 6, 1039 (2019) and **Wang *et al.*, Optica 7, 1126 (2020)**: SH BIC in a photonic-crystal slab, the direct precedent for §5.1's "dark SH resonance".
  - Zanotti *et al.* 2021.
  - Fedotova *et al.*, Nano Lett. 20, 8608 (2020) (LN metasurface, d₃₃).
  - Kalinic *et al.* 2025.
  - Hou *et al.*, arXiv:2602.00830 (TiO₂ on TFLN, which is D6-like).
- **Cavity cascade:**
  - **Cui *et al.*, Nat. Phys. 18, 497 (2022); PRL 134, 103803 (2025).**
  - Lu *et al.*, Nat. Commun. 14, 2798 (2023).
  - Musgrave, Nie & Huang, eLight 6, 15 (2026) and arXiv:2604.03474.
  - Nie *et al.*, Nat. Commun. 17, 502 (2026).
  - White, Mlynek & Schiller, EPL 35, 425 (1996).
  - Gustin *et al.*, QST 10, 025035 (2025).
  - Liu *et al.*, OL 42, 3626 (2017) and OL 44, 1456 (2019).
  - Li *et al.*, Appl. Phys. Rev. 13, 011428 (2026).
  - Gao *et al.* (Vahala), Nat. Commun. 13, 3323 (2022) (separating thermal from Kerr).
- **Theory:**
  - Hagmann *et al.*, PRA 77, 023809 (2008); Kolleck, PRA 69, 053812 (2004); Bosshard *et al.*, PRL 74, 2816 (1995).
  - QNM and TCMT: Lalanne *et al.* LPR 2018; Gigli *et al.*, ACS Photonics 7, 1197 (2020) (QNM SHG overlap); Wu & Lalanne, OE 2024; Christopoulos *et al.*, JAP 2024; Rodriguez *et al.*, OE 15, 7303 (2007); Lin *et al.*, Optica 3, 233 (2016).
- **Limits:**
  - Khurgin, SPIE 2024 (DOI 10.1117/12.3009925); Khurgin & Sun 2013; Khurgin, Clerici & Kinsey, LPR 2021.
  - Gorbounova *et al.*, OL 21, 558 (1996).
  - Shcherbakov *et al.*, PRA 100, 063847 (2019).
  - Hutchings *et al.* 1992 (nonlinear Kramers–Kronig).
- **Hosts:** Ramesh *et al.*, arXiv:2602.23246; Fathi *et al.*, Nat. Nanotechnol. (2026).

**Do not cite or check first:**

- "Khurgin, *Performance of nonlinear photonic crystal devices at high bit rates*" was not found in any index.
- The "~100× intrinsic" figure for Nie *et al.* came from an automated full-text read.
- arXiv:2505.12571 ↔ eLight 6, 15 is identified by abstract only.
- The material of Gennaro 2022 (GaAs) is inferred from the group's review.
- Coverage of Semantic Scholar was rate-limited, so a conference abstract could have been missed.

**Read in full before drafting:** Cui 2022 and 2025, Tu 2024, Panov 2020, and Minkov/Wang. These five decide the framing.

---

## 5. Paper structure (long-form fallback; for Optica see §0)

Page estimates assume a two-column Phys. Rev.-style layout. "Source" points to existing text that can be condensed. "Owner": **N** = Nikita, **S** = Sedigheh.

| § | Title | Content | Source | Owner | ~pp |
|---|---|---|---|---|---|
| 1 | Introduction | Why a larger Kerr in χ⁽²⁾ hosts; the cascading lineage (Ostrovskii → DeSalvo → Stegeman/Khurgin); bulk lever 1/Δk versus resonant lever Q₂/Δω; **cavity cascaded Kerr is established in rings** (Cui 2022/2025, Lu 2023, Musgrave 2026); metasurfaces have reached only cascaded harmonic generation (Gennaro 2022) and back-conversion as depletion (Tu 2024); the gap: a resonant cascaded *Kerr* in metasurfaces; contributions list. | MURI slides 2–4; note §1; §4 above | N (draft), S (polish) | 1.2 |
| 2 | Theory of the resonant cascade | (a) conventions (Boyd, Kleinman; F2 fixed); (b) modes and normalisation-invariant notation; (c) Green's-function route χ²·G·χ² with QNM G (Lalanne 2018; Gigli 2020 for SHG overlaps); (d) TCMT with exact steady-state elimination; (e) **boxed χ³_casc** (Kerr + cascaded loss), explicit form and η_geom; (f) **uniform-medium limit gives the bulk result** (new; F1 check); (g) **reduction to the cavity expression of Cui *et al.*** (single-mode, photon units); (h) non-resonant backgrounds dropped by single-mode theory: optical rectification + EO (Bosshard 1995) and the longitudinal bound-wave term −(2/3)d_z²/n₂² from ∇·D = 0, both small constant offsets; (i) validity: adiabaticity, saturation. | `notes/Cascaded_chi3_derivation.tex` §§1–5 (cleanest; has explicit η_geom); note §§2–6 | N | 2.5 |
| 3 | Limits | Sign tunability and ρ = γ₂/2\|Δω\|; sum rule (peak × FWHM, area), framed as the cascade-specific instance of Khurgin's strength–speed trade-off (Khurgin SPIE 2024; Gorbounova 1996; Shcherbakov 2019); saturation and parametric ceiling (Drummond 1980); speed; \|χ³\|·Δf against bulk PPLN (≈10× after F1). | note §§6–7, §10 end, §13 items 1–3; derivation §6 | N | 1.2 |
| 4 | When does cascading beat native χ⁽³⁾? | R and Q₂^cross; host survey (Table 1, corrected); per-photon phase scaling Q₁²Q₂ versus Q₁². | note §11; sim report §7 | N | 1.0 |
| 5 | **Metasurface implementation** (centre of mass) | See 5.1–5.3. | | | **3.5** |
| 5.1 | Engineering the dark SH resonance | Symmetry-protected SH BIC (absorption-limited Q₂) plus FM qBIC (critical coupling), with SH-BIC precedents for SHG (Minkov 2019; Wang 2020, photonic-crystal slab); tensor selection rules for non-zero M (d₃₃ orientation, unit-cell symmetry); Γ-point momentum matching at normal incidence; tuning Δω by geometry, temperature or EO; disorder and finite size (F3); design targets as R(Q₂, η) maps. | note §§8–9, §13 items 5–7; MURI "platform" and "fabrication constraints" slides; design study (2R Q₂ ≈ 5000 / ~10 % overlap; 1R Q₂ ≈ 8000 / ~20 %) | S (with N on selection rules) | 1.0 |
| 5.2 | Doubly-resonant metasurface: observables | Leverage table: SPM/XPM (1+R), bistability threshold 1/(1+R), FWM (1+R)²; cascade-only: Δω-tunable saturable absorption (Im χ³); short mention of the AM modulator, limiter and sensor (Tier 2; D5). | note §12; sim report §§3–4 | N | 1.0 |
| 5.3 | Singly-resonant metasurface: one design | What survives without the FM resonance (form intact; Q₁² leverage lost); **the chosen design**: geometry, SH mode, M̃ and η⁽¹ᴿ⁾ computed in the paper's normalisation, χ³⁽¹ᴿ⁾(Δω) with sign flip and Im part, R⁽¹ᴿ⁾; pulse-bandwidth and τ_p window; Z-scan sign-flip protocol and SHG calibration. **One paragraph on the native-χ⁽³⁾ confound** (Panov 2020): the observable is the part odd in Δω at fixed FW detuning; the FW is non-resonant here, so the native term is a constant offset over the ~1 nm scan. | SH-only note §§1–9; samples report; `data/design6_parameters.md`; MURI Z-scan slide | S (design numbers, Fig. 5); N (theory paragraph) | 1.5 |
| 6 | Microring: the CW-device limit | **Open by citing Cui 2022/2025 (AlN), Lu 2023 (PPLN tilt) and Musgrave 2026, then quantify for poled TFLN.** Photon-unit map (same expression as Cui's); R from a measured η_SHG; poled ≈ 10³ versus unpoled ≈ 13 against the same ring's native Kerr; no bandwidth penalty (γ₂/γ₁ = 2Q₁/Q₂); µW thresholds at Δ_eff with multimode Bogoliubov MI; EO Δω tuning; antisymmetric lock-in with the δ₁ = 0 condition; Kerr nulling (cf. Cui's suppression) and sign-tunable MI for combs; not a blockade route. **Applications angle:** CW all-optical switching and memory, comb MI engineering, Kerr-free high-Q rings. | `notes/Cascaded_chi3_microring_note.tex` (after the F4 fixes) | N | 1.2 |
| 7 | Discussion | Platform table (bulk PPLN / 2R metasurface / 1R metasurface / PPLN ring); what each is for; outlook (III–V below half-gap, MQW intersubband with quantum suppression of direct χ⁽³⁾, BTO-on-Si as EO knob). | microring §16 table; MURI outlook | N + S | 0.7 |
| 8 | Conclusions | | | S | 0.3 |
| SM | Supplement | A: full QNM/TCMT derivation with normalisation; B: uniform limit and 1-D numerical check (F1); C: TCMT validation figures; D: SH-only derivation and overlap definition; E: microring details (Bogoliubov MI, thermal/PR model, EO tuning); F: parameter tables and code pointer. | | N (A–C, E), S (D, F) | — |

---

## 6. Figures and tables (long-form fallback; for Optica see §0.2)

| Fig. | Content (panels) | Source | Owner |
|---|---|---|---|
| 1 | **Concept.** (a) Cascade energy diagram with the virtual SH (χ²·G·χ²). (b) Doubly-resonant LN metasurface: qBIC FM plus dark SH BIC, pump in, no SH out. (c) Poled microring. (d) Re/Im χ³_casc versus Δω: sign flip and even loss. | New; reuse the TikZ from the MURI "idea" slide; `figures/sign_flip.pdf` | N |
| 2 | **Theory checks.** (a) Closed form versus full TCMT. (b) Uniform-limit / 1-D Bloch check of the prefactor (new). (c) Sum rule: peak ∝ Q₂, FWHM ∝ 1/Q₂, constant product. | `validation.pdf`, `paper/checks/prefactor_check.py`, `sum_rule.pdf` | N |
| 3 | **When cascading wins.** (a) R versus Q₂ for hosts (corrected). (b) LN R(Q₂, η) contour with marks for the idealised 2R design point, the chosen 1R design and the microring η. | `cascaded_vs_native.pdf` (regenerate), `overlap_eta.pdf` (rework) | N |
| 4 | **Doubly-resonant observables.** (a) Bistability. (b) Saturation and parametric onset. (c) FWM (1+R)². | `bistability.pdf`, `parametric_onset.pdf`, `fwm_efficiency.pdf` (regenerate with corrected R) | N |
| 5 | **Singly-resonant design** (one design). (a) Unit cell and SH mode (COMSOL). (b) SH transmission and overlap versus frequency. (c) Predicted Re/Im χ³⁽¹ᴿ⁾ versus Δω. (d) Pulsed r_φ versus τ_p. | Generalise the `design6_params.py` / `design6_effects.py` pipeline to "the chosen design" | S (data, a–b), N or S (c–d via pipeline) |
| 6 | **Microring.** (a) R versus Q₂, poled/unpoled, with metasurface reference lines. (b) U_casc(Δω) against native and thermal. (c) Bistable tilt reversal under Δω → −Δω. (d) Thresholds against native. | `figures/ring/ring_R_landscape`, `ring_sign_flip`, `ring_bistability`, `ring_thresholds` (merge into one figure; fix the 3.7γ₂ label) | N |

**Tables:**

- **Table 1:** host survey (corrected).
- **Table 2:** observable leverage across 2R / 1R / ring.
- **Table 3:** platform comparison (bulk PPLN, 2R, 1R, ring): Q₂, η mechanism, R, bandwidth bottleneck, CW device?, sign-flip readout, principal risk.

---

## 7. Numbers the paper will quote

**Updated 2026-10-07: host constants are now sourced 1.55-µm values** (user decision). Provenance: `checks/constants/constants.md`; the numbers are reproduced by `checks/constants/compute.py`. **(V)** = verified; **(A)** = assumption to state.

| Quantity | Value | Note |
|---|---|---|
| LiNbO₃ constants | 2d₃₃ = 37.5 pm/V (Shoji 1997, Miller-scaled from 1.313 µm), n₂ = n_e(775) = 2.178, n₂ᴵ = 1.8×10⁻¹⁹ m²/W → χ_nat = 2.9×10⁻²¹ m²/V², F = 0.10 | n₂ᴵ = 1.8×10⁻¹⁹ is the thin-film value quoted by Kovach 2020, but **no primary measurement was found**. The estimates used span 1.8–3.2×10⁻¹⁹. The 0.91×10⁻¹⁹ Z-scan is excluded because it includes the bulk cascade. With d₃₃ between 18.8 and 22.7 pm/V, R spans 11–30, as stated in the paper. |
| LN design point | Q₂ = 4000 (A), η = 0.3 (A), Δω = 2γ₂ → **R ≈ 20** from Q₂ηF/6 (19 with the full Lorentzian), Q₂ at R = 1 ≈ 200; Δf = γ₂/2π ≈ 0.1 THz | Was R ≈ 49 with the 1.064-µm d₃₃, χ_nat = 2×10⁻²¹ and n_o. The notes, reports, code and decks still use the old constants. |
| Host survey | F: GaP 4.5×10⁻³, AlN 2.5×10⁻³, Al₀.₂Ga₀.₈As 3.7×10⁻³ → Q₂ at R = 1 between 4×10³ and 8×10³; the three cannot be ranked within their spreads | |χ²| is the tensor norm: (2/√3)·2d₁₄ for zincblende. GaAs is dropped, since its SH at 775 nm lies above the gap. **GaP is an upper estimate:** d₁₄ = 70.6 pm/V at 1.064 µm, with its SH above the indirect gap. Shoji's 1.313-µm value, reportedly 36.8 pm/V, would give F = 1.6×10⁻³ and Q₂(R=1) ≈ 1.2×10⁴. **Check Shoji 1997 Table 11 in a library copy.** |
| Bulk comparison | Bulk LN \|χ_casc\| ≈ 5×10⁻²¹ m²/V², bandwidth > 10 THz → product > 5×10⁻⁸; metasurface bound (Eq. 5) 1×10⁻⁸ m²V⁻²Hz | Bulk R ≈ 2 |
| 1R designs | **Removed from the text 2026-10-08** (was "0.6–4", the samples report × 0.39, undocumented designs). The main text and Table 1 now carry a design placeholder. | To be computed for the chosen design with η⁽¹⁾ of Supplement §7 (projection onto e₁, ∫_NL\|e₁\|⁴; with a substrate, e₁^rec replaces e₁* in one factor) |
| Microring, both designs (D9) | **A, demonstrated** (Lu *et al.* 2020): R ≈ 2.8×10³ with the fitted g = 1.2 MHz, ≈600 with g = 0.56 MHz from η_SHG; P_bist 73 µW, P_MI 19 µW. **B, target ring**: R ≈ 1.8×10³ (≈400 with g from η_SHG), P_bist 32–37 µW, P_MI 8.0 µW (spread < 2% for D₂/2π = 1–50 MHz; 11.6 µW for D₂/2π = −10 MHz). Shift by κ₁ costs 8.4 fJ intracavity (native 15 pJ). At 7 µW (below P_MI): bandwidth −4% vs instantaneous Kerr, thermal shift 1.1κ₁ vs cascaded 0.7κ₁. Native, same ring: 14 mW / 9.4 mW. **Quoted range: R = 400–3000.** | Uses n₂ᴵ = 1.8×10⁻¹⁹, the same value as above. Operating points are now below the parametric threshold (Fig. 5b at 7 µW); the old 20–30 µW points were above it. |
| Bistability penalty | Threshold falls by (1+R)/4.7 (R = 10), (1+R)/5.2 (R = 20), (1+R)/5.8 (R → ∞) at Δω = −2γ₂, ρ = 1/4 | `checks/bistability_penalty.py`; supplement §4.C |
| Bulk check | Bulk cascade at 1.064 µm (d₃₃ = 25.2 pm/V, Zelmon n_e) = −3.0×10⁻¹⁹ m²/W; + native 3.9×10⁻¹⁹ → 0.89×10⁻¹⁹ vs DeSalvo's measured 0.91×10⁻¹⁹ | `checks/bulk_desalvo_check.py`; supplement §5 |

## 8. What Sedigheh needs to deliver for §5.3 (design checklist)

**Selection rule for the design.** Maximise R⁽¹ᴿ⁾ ∝ Q₂·η⁽¹ᴿ⁾ subject to three constraints:

- **Pulse window:** the SH lifetime 1/γ₂ must stay well below the available τ_p. The laser gives 1–3 ps, so Q₂ ≲ ~2000 at 1.76 µm. D6 (Q₂ = 924, 1/γ₂ = 0.43 ps) is robust; high-Q S4 loses up to 2.5× at 1 ps.
- **Disorder tolerance:** σ_ω/ω ≲ 1/Q₂.
- **Reachable pump wavelength:** 1.76 µm is confirmed.

D6 is the default if nothing better comes out; a new design should beat Q₂·η ≈ 170.

**Data to deliver:**

1. Geometry table and unit-cell volume and LN volume. Name the tensor orientation (which crystal axis is in-plane) explicitly.
2. SH eigenmode at f₀ with **loaded and intrinsic Q₂**, separating radiative from absorptive decay. Use a supercell if the cell is a 2×1 supercell.
3. The linear FM near-field at ω₁ for unit incident amplitude, at the operating pump wavelength(s).
4. **M̃ computed directly** with the paper's definition, M̃ = ∫_LN χ²_ijk f₂ᵢ* e₁ⱼ e₁ₖ d³r, where e₁ is the FM field per unit incident amplitude. Report it together with N₂ = ∫ε₀ε̃|f₂|², so that η⁽¹ᴿ⁾ = ε₀n₂²\|M̃\|²/(\|χ²\|²N₂V_cell) needs no bridge constant. *Nikita provides a reference script before handoff* (extend `design6_params.py`).
5. Effective indices n₁, n₂ (or the dispersion model) and any non-LN material that carries field.
6. Optional: a statistical disorder study (Q₂ versus σ of the asymmetry parameter) for F3.

**Authorship note:** if the design comes from Dee & Yajun (D6) or another collaborator, they are co-authors or need to consent to publication of unpublished data (D1).

---

## 9. Work plan and timeline

The constraints come from the vault plan: Boulder trip ~10-10 to 10-15; CUNY last day to be confirmed (~10-16 to 10-23); Atom start Mon 2026-10-26.

| When | What | Who |
|---|---|---|
| **Tue 10-06 to Fri 10-09** | Approve this plan. ~~Independently confirm F1~~ (done 2026-10-06, `checks/PREFACTOR.md`). Read the five framing papers listed at the end of §4.2. Send this plan plus the §8 checklist to Sedigheh; she chooses or starts a design. ~~Propagate F1/F2 into the notes and code, regenerate figures, fix the microring-note items in F4~~ (done 2026-10-06). | N (+ Claude) |
| Sat 10-10 to Thu 10-15 | Boulder trip: no writing assumed | — |
| **Fri 10-16 to Fri 10-23** | In `main.tex` numbering: draft Sections 2, 3, 4.1 (observables) and 5, and Supplement S1–S6 and S8. Draft Section 1. Make Figs. 1, 2, 3 and 5, Table 1 and Table S1. Write the M̃ reference script. Assemble the handoff package (§10). | N |
| From 10-26 | Section 4.1 (design paragraph) and 4.2, Fig. 4, Supplement S7 and S9; then integration, Section 6, the cover letter and intro polish. | S |
| Then | Internal review (Khurgin; Alù), revision, submission | S (lead), N (theory responses) |

**Manuscript location:** `paper/` in this repo, tracked in git: `main.tex` + `supplement.tex` (Optica templates); `main_prapplied.tex` (REVTeX) is the fallback. Note that the repository is public on GitHub, so anything committed under `paper/` is public.

---

## 10. Handoff package for Sedigheh (by Fri 10-23)

### Status of the manuscript (2026-10-08)

**The manuscript is a draft: an outline plus a filled theory part.** The six-agent review of 2026-10-08 (`checks/REVIEW_2026-10-08.md`) has been applied except for the items assigned to Sedigheh below.

**Written (theory, checked by independent referee agents):**
- Main text: Sections 2–3.
- Supplement:
  - §1–4 (derivations, the uniform limit, the layered-cell test, response time, loss, bistability and parametric gain);
  - §5 (sourced host constants);
  - §7 (singly resonant theory, plus the recipe for extracting M̃ and η⁽¹⁾ from full-wave fields);
  - §8 (microring: photon-number form, g from η_SHG, designs A/B/C with all thresholds, bistability at the pump-frame Δω, the multimode parametric threshold, response time, thermal and photorefractive shifts, antisymmetric extraction, tuning of Δω).
- Fig. S1 is plotted by `figures/fig_s1_checks.py`. Every microring number reproduces from `cascaded_chi3/ring_effects.summary()`.

**Draft figures (2026-10-08):**
- Figs. 2–5 are generated by `figures/fig{2_limits,3_double,4_single,5_ring}.py` with the shared style `figures/paperstyle.py` (Times/STIX to match the main template, TrueType fonts, 300 dpi, 5.25 in wide; Fig. S1 stays Palatino). Run each with `conda run -n ml_base python <script>`; Fig. 5 takes a few minutes for the parametric thresholds.
- Panels 3(a) and 4(a) are placeholders for the chosen design.
- Fig. 4 uses a draft R = 3 at |Δω| = 2γ₂; replace it with the design's value.
- **Physics noted while drawing Fig. 3(c):** with cascaded loss at ρ = 1/4 the bistability threshold falls by (1+R)/5 to (1+R)/6, not 1+R. Main §4.1 and supplement §4.C are corrected. The notes and the derivation still say 1/(1+R).
- Fig. 3(c) now uses γ₁ = γ₂ (cold mismatch 2.4γ₂; folds at 7.2/7.8 vs native 22.6/27.3 in γ₂³/|β|²). Fig. 5(b) is at 7 µW, below the parametric threshold.

**Editing pass of 2026-10-21 (main text and supplement; both build clean, main 11 pp, supplement 11 pp):**
- **Introduction.** It stays an enumerated outline. Item 1 adds that the native sign is fixed by the material. Item 3 now runs from cavities to metasurfaces and ends on the gap ("to our knowledge"). Item 4 states the contribution relative to Khurgin–Ding and lists the two rules (strength–bandwidth, F), the LiNbO₃ numbers with their bandwidths, and the Z-scan.
- **Abstract.** It reads "a second-harmonic quality factor of 4000" instead of the undefined "Q = 4000".
- **§3 reordered.** The order is now response time and Eq. (5), bulk comparison, R and hosts, then loss and stability, so Fig. 2 is discussed in panel order and the stability paragraph sits next to §4.1, which uses it. The cascaded shift S is defined in the main text.
- **§4.1.** R ≈ 20 comes first, then the observables; added that for Δω > 0 the net Kerr response reverses sign when R > 1.
- **§4.2.** It now says why the singly resonant structure suits the measurement (pump field and native response flat across the Δω scan). The parity sentence no longer says "fixed pump detuning", since there is no pump resonance.
- **§5.** The thresholds sentence now names its ring (design B). The 8 fJ figure now carries the 15 pJ native value, and the comb sentence is separate.
- **§6 and Table 1.** The discussion closes on Eq. (5) ordering the platforms by bandwidth. The "CW device" column became "Bistability" (in a cavity / – / predicted / no / 37 µW); the waveguide's "no" was unsupported.
- **§2.** The validity condition gives both bounds (|Δω|R/4 and |Δω|(R/2)^{1/2}).
- **Supplement.** S is defined at first use (§4.B), "by a factor" was added (§4.C), the DeSalvo paragraph was reordered (§5), and the two extraction conditions are labelled first and second (§8.G); a few repeated "same" were removed.

**Still Nikita's (before 10-26):**
- Section 1 (introduction prose from the outline items; the items were rewritten 2026-10-21).
- Fig. 1 (concept); polish of the draft Figs. 2, 3 and 5.
- Supplement §6 (simulation methods) and the TCMT part of §9.
- Optional: a reference script that takes a COMSOL field export to M̃, η⁽¹⁾ and R⁽¹⁾, following the recipe in supplement §7.

**Sedigheh's (the design part, marked `% FOR THE DESIGN` in the sources):**
- Main §4.1, design paragraph: unit cell, symmetry, selection rules for d₃₃, and how Δω is tuned.
- Main §4.2: the chosen singly resonant design. Give its geometry, Q₂, M̃, η⁽¹⁾, and R⁽¹⁾ at Δω = ±2γ₂, computed with the supplement §7 recipe. That replaces the bracketed placeholder in §4.2 and the "design" cells of Table 1.
- **From the 2026-10-08 review (items A2 and A5), Sedigheh's to decide or add:**
  - venue: the simulated referee would reject at Optica without a concrete 3D design, absolute figures of merit and an SNR-budgeted proposal; fallbacks are PR Applied (`main_prapplied.tex`), Optics Express, ACS Photonics, Nanophotonics, or a condensed Optics Letters;
  - absolute figures of merit next to R: |χ_casc| ≈ 6×10⁻²⁰ m²/V² is below the native χ⁽³⁾ of Si, GaP and AlGaAs, and bulk wins on |χ|Δf;
  - Z-scan detectability budget: phase per GW/cm² (~10⁻⁴ rad at 10 GW/cm²), substrate background (10–50× the film's native phase), beam waist versus mode nonlocality and SH diffraction orders, ps source, damage;
  - SH diffraction: at 1.55 µm on silica the SH lies above the first diffraction order, which limits Q₂;
  - reciprocity with a substrate (supplement §7 now states the condition e₁^rec ∝ e₁*);
  - the alternative measurement: pump-resonance shift and tilt of a doubly resonant sample versus Δω (metasurface analogue of Cui 2025); relations δω/Γ_SH = Δω/γ₂ and φ′_max = η_SH,gen(0)/4;
  - the OR+EO background (≈ −5×10⁻²¹ m²/V², flat in Δω, shifts the native baseline) for the chosen crystal cut.
- Fig. 4: unit cell and SH mode; χ_casc(Δω); simulated Z-scan for ±Δω; pulse-length dependence.
- Supplement §7 design paragraph and Z-scan model.
- Supplement §9 metasurface parameters.
- Then integration, the cover letter (points in §0.4), and submission with Andrea as corresponding author.

**Handoff items:**
- This plan.
- `checks/` (prefactor and constants provenance).
- `figures/` scripts.
- `cascaded_chi3/README.md`, for which script produces which number.

## 11. Open decisions (for Nikita)

| ID | Decision | Default / recommendation |
|---|---|---|
| D1 | Authors and order | **Decided:** Nefedkin, Esfahani, Khurgin, Alù; plus the design's authors if their design is used |
| D2 | Venue | **Decided 2026-10-06:** Optica first (Research Article, §0); fallback Phys. Rev. Applied or ACS Photonics (`main_prapplied.tex`). |
| D3 | Title | §1 options |
| D4 | 2R design point: idealised only (η as a parameter) or anchored to the designer's 2R candidate (Q₂ = 5235)? | Idealised plus an R(Q₂, η) map. The 2R candidate's overlap normalisation is not cross-calibrated. |
| D5 | Tier-2 applications (saturable absorber, AM modulator, limiter, sensor): main text or SM? | One table row each in the main text; derivations in SM |
| D6 | Affiliation and timing (CUNY last day) | Submit with CUNY affiliation if done before the last day; otherwise dual |
| D7 | Repo visibility for `paper/` | **Decided:** tracked, not ignored |
| D8 | Include the 1-D prefactor check as an SM figure? | Yes: it is the cleanest demonstration that the prefactor is exact |
| D9 | Microring anchor | **Decided:** mention both (A demonstrated, B target) |

---

## 12. Risks

| Risk | Mitigation |
|---|---|
| Prior art (checked 2026-10-06): the metasurface side is open; the microring side is largely published | §6 framed as quantification plus protocol; cite Cui prominently; no "first" claims on rings |
| Being scooped on the metasurface version: active cavity-cascade groups (L. Fan; Y.-P. Huang, 2025–26) and SH-BIC groups (S. Fan / Gerace) are one step away | Prioritise speed. An arXiv posting as soon as the theory sections and one design are done is worth more than polish. |
| Reviewer: "the sign flip could be native χ⁽³⁾ near a Mie resonance" (Panov 2020) | The odd-in-Δω-at-fixed-FW-detuning definition, plus the flat-FW argument for 1R (§5.3) |
| F1 is contested by a co-author | The three independent checks in `paper/checks/`; Nikita's own re-derivation |
| No new design in time | Fall back to D6 (already analysed; needs only the direct M̃ computation) |
| Microring section outshines the metasurface | Keep §6 to 1.5 pp and one figure; frame it as "the CW-device limit of the same theory" |
| Handoff loses context after 10-26 | §10 package; the vault progress note stays canonical |
