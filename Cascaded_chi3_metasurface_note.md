# Cascaded $\chi^{(2)}{:}\chi^{(2)} \to \chi^{(3)}_\text{eff}$ in a Doubly-Resonant Dielectric Metasurface

## 1. Setting and goal

The objective is to engineer a metasurface whose effective third-order nonlinear susceptibility, *driven only at the fundamental frequency* $\omega_1$, exceeds the bulk $\chi^{(3)}$ of the constituent material by orders of magnitude. The mechanism is the well-known *cascade* of two $\chi^{(2)}$ processes: $\omega_1+\omega_1\to 2\omega_1$ followed by $2\omega_1-\omega_1\to\omega_1$. When the second-harmonic (SH) channel is *off-resonant* (or phase-mismatched), the SH amplitude follows the fundamental adiabatically and its back-action produces a Kerr-like response on the fundamental [Ostrovskii 1967; DeSalvo *et al.* 1992; Stegeman *et al.* 1996].

The metasurface platform we have in mind consists of an array of high-index dielectric resonators (LiNbO$_3$, GaAs, GaP, Si) supporting two co-localised resonances:

- a fundamental resonance at $\omega_1$ — a *quasi*-bound state in the continuum (qBIC) with controlled radiative coupling $\gamma_1^{(\text{rad})}$ to free space, designed for in-coupling at $\omega_1$;
- a second resonance at $\omega_2 \approx 2\omega_1$ — a symmetry-protected BIC (or strongly under-coupled qBIC) whose total decay $\gamma_2$ is dominated by intrinsic absorption, giving the largest accessible quality factor $Q_2 = \omega_2/\gamma_2$.

The relevant tensor of $\chi^{(2)}$ is patterned in space (via the geometry of the meta-atoms, or by periodic poling, or by symmetry-breaking perturbations) so that the SHG overlap integral does not vanish (Section 7). High-Q metasurface SHG of exactly this kind has been demonstrated in AlGaAs and LiNbO$_3$ [Carletti *et al.* 2018; Koshelev *et al.* 2020; Anthur *et al.* 2020]; the *cascaded* $\chi^{(3)}$ regime has, to my knowledge, not yet been explicitly demonstrated in a metasurface and is the experimental target.

---

## 2. Conventions

Several factors of 2 and 3 in the final answer come from book-keeping in the field/polarisation conventions; we fix them once and for all.

**Real fields.** Physical real fields are written as a sum of complex envelopes:

$$
\mathcal{E}_i(\mathbf r,t)=\tfrac12\sum_n\!\bigl[E_{n,i}(\mathbf r,t)\,e^{-i\omega_n t}+\text{c.c.}\bigr],\qquad
\mathcal{P}_i(\mathbf r,t)=\tfrac12\sum_n\!\bigl[P_{n,i}(\mathbf r,t)\,e^{-i\omega_n t}+\text{c.c.}\bigr],
$$

with $E_{n,i},P_{n,i}$ slowly varying compared with $1/\omega_n$.

**Susceptibility tensors (Boyd convention).** With this field convention,

$$
P_i^{(2)}(\omega_\sigma)=\varepsilon_0\!\!\sum_{(p,q)}\!\chi^{(2)}_{ijk}(\omega_\sigma;\omega_p,\omega_q)\,E_j(\omega_p)\,E_k(\omega_q),
$$

where the sum runs over distinct ordered pairs of frequencies summing to $\omega_\sigma$. Indices are summed by Einstein convention.

For the two processes at hand:

$$
\boxed{\;P_i^{(2)}(2\omega_1)=\varepsilon_0\,\chi^{(2)}_{ijk}(2\omega_1;\omega_1,\omega_1)\,E_j(\omega_1)\,E_k(\omega_1)\;}\qquad\text{(SHG)}
$$

$$
\boxed{\;P_i^{(2)}(\omega_1)=2\varepsilon_0\,\chi^{(2)}_{ijk}(\omega_1;2\omega_1,-\omega_1)\,E_j(2\omega_1)\,E_k^*(\omega_1)\;}\qquad\text{(DFG)}
$$

The factor $2$ in DFG is the number of orderings $(2\omega_1,-\omega_1)$ and $(-\omega_1,2\omega_1)$. **This is the factor missing in the original note**, and dropping it violates Manley–Rowe (Section 4).

The native Kerr response on the fundamental, for comparison, is

$$
P_i^{(3)}(\omega_1)=3\varepsilon_0\,\chi^{(3)}_{ijkl}(\omega_1;\omega_1,\omega_1,-\omega_1)\,E_j(\omega_1)E_k(\omega_1)E_l^*(\omega_1),
$$

where the factor $3$ counts the distinct orderings of $(\omega_1,\omega_1,-\omega_1)$.

We assume Kleinman symmetry (full permutation symmetry of all index pairs $(i,\omega_\sigma),(j,\omega_p),(k,\omega_q)$), which is good for both LiNbO$_3$ and GaAs at telecom-band fundamentals well below the electronic gap. In particular,

$$
\chi^{(2)}_{ijk}(2\omega_1;\omega_1,\omega_1)=\chi^{(2)}_{kij}(\omega_1;2\omega_1,-\omega_1)\equiv\chi^{(2)}_{ijk}(\mathbf r),
$$

i.e. we use a single symmetric tensor $\chi^{(2)}_{ijk}(\mathbf r)$ in what follows.

---

## 3. Modes and energy normalisation

Let $\mathbf f_n(\mathbf r),\mathbf h_n(\mathbf r)$ be the eigenmodes of the linear, lossless, dispersive Maxwell problem,

$$
\nabla\times\nabla\times\mathbf f_n(\mathbf r)=\frac{\omega_n^2}{c^2}\,\varepsilon_r(\mathbf r,\omega_n)\,\mathbf f_n(\mathbf r),
$$

with $\mathbf f_n$ a complex vector field encoding the polarisation pattern of mode $n$ (and the magnetic mode $\mathbf h_n=\frac{i}{\omega_n\mu_0}\nabla\times\mathbf f_n$). Note that both the *vectorial* nature of $\mathbf f_n$ and the dispersion of $\varepsilon_r$ are essential — the original note's scalar, dispersionless $f(\mathbf r)$ is a simplification that drops information needed for the tensorial overlap below.

The mode is energy-normalised such that

$$
\boxed{\;\frac12\!\int\!\!\bigl[\,\varepsilon_0\,\widetilde\varepsilon_n(\mathbf r)\,|\mathbf f_n|^2+\mu_0|\mathbf h_n|^2\,\bigr]\,d^3r=1,\qquad
\widetilde\varepsilon_n(\mathbf r)\equiv\frac{\partial(\omega\varepsilon_r(\mathbf r,\omega))}{\partial\omega}\biggr|_{\omega=\omega_n}\;}
$$

(Brillouin's expression for the energy density of a time-harmonic field in a dispersive medium [Landau & Lifshitz §80]). Defining the modal amplitude $a_n(t)$ by

$$
E_{n,i}(\mathbf r,t)=a_n(t)\,f_{n,i}(\mathbf r),
$$

the time-averaged energy stored in mode $n$ is exactly $|a_n|^2$ (units $\sqrt{\text{J}}$).

For non-dispersive, lossless dielectrics one may replace $\widetilde\varepsilon_n\to\varepsilon_r(\mathbf r,\omega_n)\equiv n^2(\mathbf r,\omega_n)$. **It is important not to use a single $n$ for both modes**: in LiNbO$_3$ at $\lambda_1=1550$ nm, $n_o(\omega_1)\simeq 2.21$ and $n_o(2\omega_1)\simeq 2.26$ (a $\sim$2 % difference), and the magnitude of $\chi^{(2)}_{ijk}$ also differs at the two frequencies via Miller's rule.

---

## 4. Coupled-mode equations

Substituting the modal expansion into Maxwell's equations and projecting onto $\mathbf f_n^*$ (slowly-varying-envelope approximation; e.g. [Suh, Wang & Fan 2004; Rodriguez *et al.* 2007]) gives, after adding phenomenological loss rates and the input port for the FM:

$$
\boxed{\;\frac{da_2}{dt}=\bigl(i\Delta\omega-\tfrac{\gamma_2}{2}\bigr)a_2+i\beta\,a_1^2\;}
$$

$$
\boxed{\;\frac{da_1}{dt}=-\tfrac{\gamma_1}{2}a_1+i\beta^*\,a_1^*\,a_2+i\alpha_3\,|a_1|^2 a_1+\sqrt{\gamma_1^{(\text{rad})}}\,s_+\;}
$$

where $\Delta\omega\equiv 2\omega_1-\omega_2$ and $s_+$ is the input wave amplitude normalised so that $|s_+|^2$ is incident power; the input–output relation is $s_-=-s_++\sqrt{\gamma_1^{(\text{rad})}}\,a_1$. The total FM decay is $\gamma_1=\gamma_1^{(\text{rad})}+\gamma_1^{(\text{nr})}$; for the SH mode we have, by design, $\gamma_2\simeq\gamma_2^{(\text{nr})}$.

**The SHG and DFG couplings carry the same coupling constant $\beta$.** This is the Manley–Rowe / energy-conservation requirement and can be checked directly: in the lossless limit ($\gamma_n,s_+\to 0$, no Kerr), the rate of change of total mode energy is

$$
\frac{d}{dt}(|a_1|^2+|a_2|^2)=2\,\mathrm{Im}\!\bigl[\beta\,a_1^{*2}a_2-\beta^{(1)}\,a_1^2a_2^*\bigr],
$$

where $\beta^{(1)}$ denotes the (a priori unknown) coupling appearing in the FM equation. This vanishes for arbitrary $a_1,a_2$ if and only if $\beta^{(1)}=\beta^*$. The same conclusion follows from the Hamiltonian formulation in [Rodriguez *et al.* 2007, Eqs. 8–9], where energy conservation is built in from the start. **It is the factor 2 in the DFG polarisation (Section 2) that produces this equality** — without it the SHG and DFG couplings would differ by a factor of 2 and energy would not be conserved. The original note's "$\kappa$ vs. $\kappa/2$" form violates Manley–Rowe.

The explicit, tensorial expression for $\beta$ is

$$
\boxed{\;\beta=\frac{\omega_1}{4}\!\int\!\chi^{(2)}_{ijk}(\mathbf r)\,f_{2,i}^*(\mathbf r)\,f_{1,j}(\mathbf r)\,f_{1,k}(\mathbf r)\,d^3r\;}\qquad[\text{units:}\ \text{J}^{-1/2}\,\text{s}^{-1}]
$$

and analogously for the bare Kerr coefficient,

$$
\boxed{\;\alpha_3=\frac{3\omega_1}{8}\!\int\!\chi^{(3)}_{ijkl}(\mathbf r)\,f_{1,i}^*(\mathbf r)f_{1,j}(\mathbf r)f_{1,k}(\mathbf r)f_{1,l}^*(\mathbf r)\,d^3r\;}\qquad[\text{units:}\ \text{J}^{-1}\,\text{s}^{-1}].
$$

The numerical prefactors $\omega_1/4$ and $3\omega_1/8$ depend on the chosen normalisation and field convention; the structure of the integrals — full contraction of the tensor with the vectorial mode profiles — is convention-independent and is the practical input for full-wave numerical evaluation. The factor 3 in $\alpha_3$ is the SPM degeneracy from Section 2 and was absent in the original note.

---

## 5. Adiabatic elimination and the effective $\chi^{(3)}$

Assume $|\Delta\omega|\gg\gamma_2/2$ and $|\dot a_2|\ll|\Delta\omega\,a_2|$ (slow modulation of the FM compared with the SH response time). Setting $da_2/dt\to 0$ in the SH equation,

$$
a_2 = \frac{i\beta\,a_1^2}{-i\Delta\omega+\gamma_2/2}=\frac{-\beta\,a_1^2}{\Delta\omega+i\gamma_2/2}\simeq\frac{-\beta\,a_1^2}{\Delta\omega}\!\left(1-\frac{i\gamma_2}{2\Delta\omega}\right).
$$

Substituting into the FM equation,

$$
\frac{da_1}{dt}=\Bigl[-\tfrac{\gamma_1}{2}-\underbrace{\tfrac{|\beta|^2\gamma_2}{2\Delta\omega^2}}_{\text{cascaded loss}}|a_1|^2\Bigr]a_1+i\Bigl[\underbrace{\alpha_3-\tfrac{|\beta|^2}{\Delta\omega}}_{\text{effective Kerr}}\Bigr]|a_1|^2 a_1+\sqrt{\gamma_1^{(\text{rad})}}\,s_+.
$$

Two physical effects emerge automatically:

(i) **Effective Kerr coefficient:**

$$
\boxed{\;\alpha_3^{\text{eff}}=\alpha_3-\frac{|\beta|^2}{\Delta\omega}\;}
$$

The cascaded contribution can have *either* sign relative to the bare Kerr response, depending on whether $\Delta\omega\equiv 2\omega_1-\omega_2$ is positive or negative. This sign-tunability is one of the hallmarks of cascaded nonlinearity: the same metasurface can be tuned through zero, used for self-defocusing, or used for self-focusing simply by walking the SH resonance through twice the FM frequency. This was used dramatically in soliton compression experiments [DeSalvo 1992; Liu *et al.* 1999].

(ii) **Two-photon-like loss:** $\gamma_\text{eff}^{(2)}=|\beta|^2\gamma_2/\Delta\omega^2$, an intensity-dependent loss term that comes from the residual width of the SH resonance. It is the analogue of two-photon absorption (here phenomenological, governed by $\gamma_2$ rather than by the imaginary part of $\chi^{(3)}$) and is unavoidable as the price of the resonant enhancement.

To compare with bulk material parameters, invert the definition of $\alpha_3$: a bulk Kerr coefficient $\chi^{(3)}_\text{eff}$ uniformly filling the mode produces $\alpha_3=(3\omega_1/8)\chi^{(3)}_\text{eff}\!\int\!|\mathbf f_1|^4 d^3r$, so the cascaded contribution is equivalent to

$$
\boxed{\;\chi^{(3)}_{\text{casc}}=-\frac{|\beta|^2/\Delta\omega}{(3\omega_1/8)\!\int\!|\mathbf f_1|^4\,d^3r}=-\frac{8}{3\omega_1\Delta\omega}\,\frac{|\beta|^2}{\int|\mathbf f_1|^4 d^3r}.\;}
$$

Inserting the explicit $\beta$,

$$
\chi^{(3)}_{\text{casc}}\;=\;-\frac{\omega_1}{6\Delta\omega}\;\frac{\Bigl|\!\int\chi^{(2)}_{ijk}(\mathbf r)\,f_{2,i}^*\,f_{1,j}\,f_{1,k}\,d^3r\,\Bigr|^2}{\int|\mathbf f_1|^4\,d^3r}\,.
$$

The structure is

$$
\chi^{(3)}_\text{casc}\;\sim\;\underbrace{\frac{\omega_1}{|\Delta\omega|}}_{\text{resonance enhancement}}\times\;\underbrace{|\chi^{(2)}|^2}_{\text{material}}\;\times\;\underbrace{\eta}_{\text{mode-overlap factor}},
$$

where $\eta$ is a dimensionless overlap whose maximum is achieved when the SH mode profile is the symmetric square of the FM mode profile, properly contracted with $\chi^{(2)}_{ijk}$.

---

## 6. The bandwidth–enhancement sum rule

The peak enhancement is bounded above by the SH linewidth: as $|\Delta\omega|\to\gamma_2/2$ the adiabatic approximation breaks down and $a_2$ saturates at $\sim\beta a_1^2/\gamma_2$. The maximum useful enhancement is therefore

$$
\bigl|\chi^{(3)}_\text{casc}\bigr|_\text{max}\;\sim\;\frac{\omega_1}{\gamma_2}\,|\chi^{(2)}|^2\,\eta\;=\;Q_2\;\times\;\frac{|\chi^{(2)}|^2\,\eta}{\omega_1}\,.
$$

But the bandwidth over which $\chi^{(3)}_\text{casc}$ is "available" (i.e. the FM modulation frequency over which $a_2$ tracks $a_1^2$ adiabatically) is itself $\sim\gamma_2$. Hence

$$
\boxed{\;\bigl|\chi^{(3)}_\text{casc}\bigr|\,\times\,\Delta f_\text{usable}\;\lesssim\;|\chi^{(2)}|^2\,\eta\;\;\text{(independent of }Q_2\text{)}.\;}
$$

This is a sum rule: $Q_2$ does not enter the product. **Resonant enhancement redistributes nonlinear action into a tall, narrow peak; it does not increase the area under the curve.** This is the same rule that limits resonant linear absorption (oscillator-strength sum rule), purely classical, and it is what disqualifies high-Q double-resonance schemes from delivering "fast" $\chi^{(3)}$ in the ultrafast sense (cf. our earlier discussion). For the nominal numbers below ($Q_2\sim 4\times 10^3$ at $\lambda_1=1.55\,\mu$m), $\Delta f_\text{usable}\sim 50$ GHz — narrowband by Kerr-electronics standards (THz), broadband by photonic-cavity standards.

For non-resonant (phase-mismatched bulk) cascading the same sum rule applies with $\gamma_2$ replaced by the inverse coherence time $|\Delta k| v_g$, which gives modest enhancement over a much larger bandwidth — this is the regime exploited in soliton compression and bulk all-optical switching.

---

## 7. Tensor structure and phase matching

The overlap that defines $\beta$,

$$
\beta\propto\!\int\chi^{(2)}_{ijk}(\mathbf r)\,f_{2,i}^*\,f_{1,j}\,f_{1,k}\,d^3r,
$$

vanishes identically for several reasons that must each be addressed in the metasurface design:

1. **Bulk inversion symmetry.** If a meta-atom has inversion symmetry, $\chi^{(2)}_{ijk}\equiv 0$. Hence one needs either a non-centrosymmetric crystal substrate (LiNbO$_3$ — point group $3m$; GaAs, GaP — point group $\bar{4}3m$ ($T_d$); AlGaAs, etc.) *and* a meta-atom geometry that does not impose inversion symmetry on the mode profile.

2. **Mode-symmetry mismatch.** Even with $\chi^{(2)}_{ijk}\neq 0$ in the material, the integrand has a definite parity under any spatial symmetry shared by both modes. If the FM mode is, say, magnetic-dipole-like (odd under a mirror) and the SH is electric-dipole-like (even), the integrand may be odd under that mirror and the integral vanishes by the symmetry of the unit cell. Quasi-BIC engineering is precisely the controlled breaking of one such symmetry to lift this constraint while keeping $Q_2$ large.

3. **Quasi-phase matching across the metasurface period.** For an extended metasurface, the integrand carries Bloch phase factors $e^{i(2\mathbf k_1-\mathbf k_2)\cdot\mathbf R}$. Periodic poling or spatial modulation of $\chi^{(2)}_{ijk}(\mathbf r)$ across unit cells provides a reciprocal-lattice vector $\mathbf G$ that compensates the residual momentum $2\mathbf k_1-\mathbf k_2-\mathbf G=0$. In a single-unit-cell-thick metasurface with both modes at $\mathbf k_\parallel=0$ (high-symmetry $\Gamma$ point), this momentum-matching is automatic and the only requirement is on the *internal* symmetry of the mode profiles.

For LiNbO$_3$ (point group $3m$), the non-zero independent components of $\chi^{(2)}_{ijk}$ in contracted notation are $d_{15}=d_{31}, d_{22}=-d_{21}=-d_{16}$, and $d_{33}$, with $|d_{33}|\simeq 25$ pm/V at 1064 nm being the largest. With the $\hat z$ axis parallel to the optic axis, the largest overlap is achieved by mode profiles with strong $|f_{1,z}|^2 f_{2,z}^*$ content over the meta-atom volume. Z-cut LiNbO$_3$ with field components engineered along the optic axis is therefore preferred. The same considerations apply, with the appropriate tensor structure, for $\bar{4}3m$ semiconductors (largest component $d_{14}=d_{25}=d_{36}\simeq 100$ pm/V for GaAs).

For LiNbO$_3$, $\chi^{(2)}\equiv 2d_{33}\simeq 50$ pm/V $=5\times 10^{-11}$ m/V; for GaAs, $\chi^{(2)}\simeq 2\times 10^{-10}$ m/V.

---

## 8. Metasurface implementation: BIC and quasi-BIC platforms

The high-$Q_2$ requirement is the key technical hurdle, and bound states in the continuum (BIC) provide a particularly elegant route [Hsu *et al.* 2016]. A symmetry-protected BIC is a guided-wave-like solution that lies inside the radiation continuum but is *symmetry-forbidden* from coupling to plane waves at $\mathbf k_\parallel=0$. Realistic structures are quasi-BICs (qBICs): finite arrays, slight disorder, or deliberate symmetry breaking convert the ideal BIC into a high-$Q$ leaky resonance, whose $Q$ scales as $1/\alpha^2$ with the asymmetry parameter $\alpha$ [Koshelev *et al.* 2018]. Demonstrated $Q$ factors in dielectric metasurfaces are routinely $10^3$–$10^4$ and have reached $\sim 10^5$ in carefully optimised devices.

Two practical architectures:

(a) **Two-mode qBIC pair.** Both the FM and SH are qBICs of the same array, with the SH qBIC very close to the symmetry-protected limit (large $Q_2$) and the FM qBIC moderately broadened ($Q_1\sim 10^2$–$10^3$) to allow efficient in-coupling at $\omega_1$. The SH leakage to the radiation continuum is irrelevant for the cascaded process — we *don't want* SH out-coupling here.

(b) **FM Mie–qBIC + SH symmetry-protected BIC.** The FM is a low-$Q$ Mie resonance (dipolar, easy to pump); the SH lives in a true BIC mode of the array, with $\gamma_2$ set by absorption only. Best for maximum enhancement, hardest to fabricate.

In both cases the requirement on the SH mode is purely *internal*: high $Q_2$, large field overlap with the FM-squared, no need for any external connection.

A subtle point: the cascaded $\chi^{(3)}$ effect manifests at the FM, so the experimental signature is **intensity-dependent reflection/transmission at $\omega_1$** (self-phase modulation, intensity-dependent line shape, optical bistability for sufficient pump), not SH light leaving the sample. In fact, a *good* implementation produces *no measurable* SH radiation while the back-action on the fundamental is maximum. This is also the experimental fingerprint that distinguishes the cascaded process from the bare Kerr response.

---

## 9. A worked numerical estimate (LiNbO$_3$ metasurface)

Take

| quantity | value |
|---|---|
| $\lambda_1$ | $1.55\,\mu$m |
| $\omega_1$ | $1.22\times 10^{15}$ rad/s |
| $n_1, n_2$ | $2.21,\,2.26$ |
| $|\chi^{(2)}|$ (=$2d_{33}$) | $5\times 10^{-11}$ m/V |
| $|\chi^{(3)}|_\text{native}$ | $\sim 2\times 10^{-21}$ m$^2$/V$^2$ |
| $Q_2$ | $4\times 10^3$ |
| $\Delta\omega$ | $2\gamma_2 = 2\omega_1/Q_2 = 6\times 10^{11}$ rad/s |
| $\eta\equiv \bigl|\!\int\!\chi^{(2)}f_{2}^*f_1^2\bigr|^2/(\chi^{(2)2}\!\int|\mathbf f_1|^4)$ (normalised overlap) | $\sim 0.1$ |

With these numbers,

$$
\bigl|\chi^{(3)}_\text{casc}\bigr|\;\sim\;\frac{\omega_1}{6\,\Delta\omega}\,|\chi^{(2)}|^2\,\eta\;\sim\;\frac{1}{6}\times\frac{Q_2}{4}\times 0.1\times|\chi^{(2)}|^2\;\sim\;15\times|\chi^{(2)}|^2.
$$

That is

$$
\bigl|\chi^{(3)}_\text{casc}\bigr|\sim 15\times(5\times 10^{-11})^2\text{ m}^2/\text{V}^2\;\simeq\;4\times 10^{-20}\text{ m}^2/\text{V}^2,
$$

i.e. roughly **20× the native bulk** $\chi^{(3)}$ of LiNbO$_3$, available over a usable bandwidth $\sim\gamma_2/(2\pi)\simeq 50$ GHz centred at the FM resonance. With $Q_2=10^4$ (achievable in optimised qBIC LiNbO$_3$ metasurfaces) and $\eta=0.3$, the same estimate gives a $\sim 10^{2}$ enhancement at the cost of a $\sim 20$ GHz bandwidth.

For comparison, a phase-mismatched bulk PPLN slab at the same $|\chi^{(2)}|$ gives non-resonant cascaded $\chi^{(3)}_\text{casc}\sim 2\,|\chi^{(2)}|^2/(\Delta k\,L)\sim$ few $\times 10^{-21}$ m$^2$/V$^2$ over THz bandwidth — comparable to native, broadband. The metasurface trades 1–2 orders of magnitude in $|\chi^{(3)}_\text{eff}|$ for 4 orders of magnitude in bandwidth. Both regimes have their applications; the high-$Q$ regime is suited to CW-pumped bistability, parametric oscillation thresholds, and frequency-comb stabilisation rather than to ultrafast switching.

---

## 10. Validity, caveats, and what to verify

1. **Adiabatic-elimination breakdown.** The expression for $\chi^{(3)}_\text{casc}$ diverges as $\Delta\omega\to 0$. The correct, non-divergent expression is obtained by *not* eliminating $a_2$: the system reduces to coupled SHG and only on resonance ($\Delta\omega=0$) the steady-state SH amplitude is bounded by $|a_2|=|\beta||a_1|^2/(\gamma_2/2)$. The cascaded-Kerr description is valid for $|\Delta\omega|\gtrsim 2\gamma_2$, i.e. several SH linewidths off resonance.

2. **Pump depletion / saturation.** The derivation is small-signal: it assumes $|a_2|\ll|a_1|$, equivalently $|\beta||a_1|^2/|\Delta\omega|\ll|a_1|$, i.e. $|\beta||a_1|\ll|\Delta\omega|$. Above this threshold the SH is no longer slaved to the FM and the response saturates; in fact, exactly at the saturation point one finds an analogue of the parametric oscillation threshold. The maximum useful FM intensity is therefore set by $|\Delta\omega|/|\beta|$, which scales as $1/\sqrt{Q_2}$ at fixed $\Delta\omega/\gamma_2$ — high-$Q$ devices saturate at lower powers.

3. **Two-photon-like loss.** The imaginary part of $\alpha_3^\text{eff}$ — the cascaded two-photon-absorption analogue — is unavoidable and bounded below by the same overlap. It does *not* dump energy outside the system in the limit of vanishing $\gamma_2^\text{rad}$, but rather absorbs into the SH mode and from there into the material. This is identical to two-photon absorption from an experimental standpoint and limits the achievable Kerr-driven phase shift before bleaching.

4. **Thermo-optic and free-carrier nonlinearities** in dielectric metasurfaces are notoriously easy to confuse with "Kerr". At the field intensities required to see a measurable cascaded $\chi^{(3)}$ phase shift (several MW/cm$^2$ in the resonator), thermal effects in absorbing materials (Si, GaAs near the band edge) often dominate. LiNbO$_3$ is comparatively benign in this respect.

5. **Tensor-overlap evaluation.** The dimensionless overlap $\eta$ is assumption-driven in this note. It must be evaluated by full-wave simulation (e.g. COMSOL or MEEP eigenmode solver, then post-processing the tensor contraction in the volume of the meta-atoms). The 10% figure is plausible for a properly engineered qBIC pair but is by no means generic — random mode pairings typically give $\eta$ of order $10^{-3}$ or less.

6. **Avoid the "single $n$" simplification.** Use $\widetilde\varepsilon(\mathbf r,\omega_1)$ and $\widetilde\varepsilon(\mathbf r,2\omega_1)$ separately in the energy normalisation. The dispersion of $\chi^{(2)}$ between $\omega_1$ and $2\omega_1$ is a few percent in LiNbO$_3$ but should be carried explicitly when comparing with experiment.

---

## 11. Summary of the principal result and where it improves on the prior derivation

The main result of the analysis is

$$
\boxed{\;\chi^{(3)}_\text{casc}(\omega_1)=-\frac{\omega_1}{6\,\Delta\omega}\;\frac{\bigl|\!\int\!\chi^{(2)}_{ijk}(\mathbf r)\,f_{2,i}^*\,f_{1,j}\,f_{1,k}\,d^3r\,\bigr|^2}{\int|\mathbf f_1(\mathbf r)|^4\,d^3r}\,,\;}
$$

valid for $\gamma_2/2\ll|\Delta\omega|\ll\omega_1$ and $|\beta\,a_1|\ll|\Delta\omega|$, with $\beta$, $\alpha_3$ as defined in Section 4. Compared with the original note, the corrected derivation:

- uses the full vector mode profile with proper energy normalisation (Brillouin) accounting for dispersion, rather than a scalar $f(\mathbf r)$ with single $n$;
- carries the tensorial $\chi^{(2)}_{ijk}$ (and $\chi^{(3)}_{ijkl}$) and gives an unambiguous overlap integral;
- includes the factor $2$ in the DFG polarisation and the factor $3$ in SPM, ensuring Manley–Rowe is satisfied and the SHG/DFG couplings carry the *same* coupling constant $\beta$;
- gives the correct power of $n$ (the original was off by one factor of $n$ in the bare-Kerr coupling), although the dependence on $n$ disappears in the final formula above where the answer is expressed directly through the overlap and the bulk $\chi^{(2)}$;
- includes both the cascaded Kerr *and* the cascaded two-photon-loss term, which is essential for understanding the saturation behaviour;
- exposes the bandwidth–enhancement sum rule, which is the correct framing for comparison with bulk (or any non-resonant) Kerr-equivalent processes;
- adds the symmetry/phase-matching constraints that govern $\eta$ in a metasurface platform.

Whether the resulting "fast Kerr from slow SHG" claim survives depends on the application: for narrowband, CW-pumped nonlinear photonics (bistability, frequency-comb interactions, low-threshold parametric oscillation) the metasurface route can outperform bulk by 1–2 orders of magnitude; for ultrafast applications the sum rule precludes any improvement over bulk.

---

## References

1. L. A. **Ostrovskii**, *Self-action of light in crystals*, JETP Lett. **5**, 272 (1967). — Original cascading idea.
2. R. **DeSalvo**, D. J. Hagan, M. Sheik-Bahae, G. Stegeman, E. W. Van Stryland, *Self-focusing and self-defocusing by cascaded second-order effects in KTP*, Opt. Lett. **17**, 28 (1992). — First measurement.
3. G. I. **Stegeman**, D. J. Hagan, L. Torner, *$\chi^{(2)}$ cascading phenomena and their applications to all-optical signal processing, mode-locking, pulse compression and solitons*, Opt. Quantum Electron. **28**, 1691 (1996). — Comprehensive review.
4. M. **Bache**, J. Moses, F. W. Wise, *Scaling laws for soliton pulse compression by cascaded quadratic nonlinearities*, J. Opt. Soc. Am. B **24**, 2752 (2007). — Practical scaling.
5. R. W. **Boyd**, *Nonlinear Optics*, 4th ed. (Academic Press 2020). — Conventions and tensor formalism in Sections 1–2.
6. L. D. **Landau** & E. M. Lifshitz, *Electrodynamics of Continuous Media*, 2nd ed. §80. — Brillouin energy density for dispersive media.
7. W. **Suh**, Z. Wang, S. Fan, *Temporal coupled-mode theory and the presence of non-orthogonal modes in lossless multimode cavities*, IEEE J. Quantum Electron. **40**, 1511 (2004). — Linear TCMT framework.
8. A. **Rodriguez**, M. Soljačić, J. D. Joannopoulos, S. G. Johnson, *$\chi^{(2)}$ and $\chi^{(3)}$ harmonic generation at a critical power in inhomogeneous doubly resonant cavities*, Opt. Express **15**, 7303 (2007). — Nonlinear TCMT with explicit overlap integrals; closest analogue to the present derivation.
9. P. T. **Kristensen**, K. Herrmann, F. Intravaia, K. Busch, *Modeling electromagnetic resonators using quasinormal modes*, Adv. Opt. Photon. **12**, 612 (2020). — Mode normalisation in lossy cavities (essential for finite-Q).
10. C. W. **Hsu**, B. Zhen, A. D. Stone, J. D. Joannopoulos, M. Soljačić, *Bound states in the continuum*, Nat. Rev. Materials **1**, 16048 (2016). — BIC review.
11. K. **Koshelev**, S. Lepeshov, M. Liu, A. Bogdanov, Y. Kivshar, *Asymmetric metasurfaces with high-$Q$ resonances governed by bound states in the continuum*, Phys. Rev. Lett. **121**, 193903 (2018). — qBIC scaling $Q\sim 1/\alpha^2$.
12. L. **Carletti**, K. Koshelev, C. De Angelis, Y. Kivshar, *Giant nonlinear response at the nanoscale driven by bound states in the continuum*, Phys. Rev. Lett. **121**, 033903 (2018). — qBIC-enhanced SHG in AlGaAs.
13. K. **Koshelev**, S. Kruk, E. Melik-Gaykazyan, J.-H. Choi, A. Bogdanov, H.-G. Park, Y. Kivshar, *Subwavelength dielectric resonators for nonlinear nanophotonics*, Science **367**, 288 (2020). — Demonstration with $Q\sim 10^3$.
14. A. P. **Anthur**, H. Zhang, R. Paniagua-Dominguez, D. A. Kalashnikov, S. T. Ha, T. W. W. Maß, A. I. Kuznetsov, L. Krivitsky, *Continuous wave second harmonic generation enabled by quasi-bound-states in the continuum on gallium phosphide metasurfaces*, Nano Lett. **20**, 8745 (2020).
15. K. **Koshelev** & Y. Kivshar, *Dielectric resonant metaphotonics*, ACS Photonics **8**, 102 (2021). — Comprehensive review of the platform.
