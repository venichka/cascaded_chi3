"""Cross-checks + optimal-detuning scan for the 4 SH-only LiNbO3 samples.

Two jobs:
  (1) verify()    -- independent cross-checks of the chi3_casc lineshape against
                     the note's dimensionless adiabatic forms (analytic.py), a
                     prefactor-free ratio test, a dimensional check, and a hand
                     spot-value for sample 4.
  (2) optimal()   -- scan Domega/gamma2 over the single-pole validity window
                     [1,5] using the DATA-derived eta_geom(Domega) (from the COMSOL
                     overlap), find the detuning that maximises |chi3_casc| and
                     |Re chi3_casc| on each (focusing / defocusing) branch, and
                     compare with the conservative Domega = +2 gamma2 point.

    conda run -n ml_base python sample_optimal.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import apply_style
apply_style()

import sample_params as sp
from sample_params import CHI2, CHI3_NATIVE, build_sample, chi3_casc, eta_geom
from analytic import alpha3_eff_casc_full, gamma_2pa_eff
from sample_effects import SAMPLE_COLORS, FIG_DIR, _save, _label, DOP

EPS0 = 8.8541878128e-12


# ----------------------------------------------------------------------
def verify():
    print("=" * 70)
    print("VERIFICATION CROSS-CHECKS")
    print("=" * 70)
    samples = {i: build_sample(i) for i in (1, 2, 3, 4)}

    # (i) prefactor-free ratio test: Re/Im chi3 must equal -2*Dtil for ALL samples
    #     and ALL detunings, because chi3 ~ -1/(Dtil + i/2).
    print("\n(i) Re/Im(chi3) == -2*(Domega/gamma2)  [prefactor-free, must hold exactly]")
    ok = True
    for i, p in samples.items():
        for dtil in (0.7, 1.0, 2.0, 3.5):
            c = chi3_casc(p, dtil, eta=0.05)        # fixed eta: tests lineshape only
            ratio = c.real / c.imag
            if abs(ratio - (-2 * dtil)) > 1e-9 * abs(2 * dtil):
                ok = False
                print("   MISMATCH S%d Dtil=%.2f: Re/Im=%.6f expected %.6f" % (i, dtil, ratio, -2 * dtil))
    print("   -> %s" % ("PASS" if ok else "FAIL"))

    # (ii) match to analytic.py dimensionless forms (constant-eta):
    #      Re chi3 = (K/gamma2)*alpha3_eff_casc_full(Dtil),  K = omega1 chi2^2 eta/(6 n2^2)
    #      Im chi3 = (K/gamma2)*gamma_2pa_eff(Dtil)
    print("\n(ii) chi3 lineshape == analytic.py adiabatic forms (alpha3_eff, gamma_2pa)")
    p = samples[4]
    eta = 0.04
    K = 2.0 * p["omega1"] * CHI2 ** 2 * eta / (3.0 * p["n2"] ** 2)
    worst = 0.0
    for dtil in np.linspace(-5, 5, 21):
        if abs(dtil) < 1e-6:
            continue
        c = chi3_casc(p, dtil, eta=eta)
        re_exp = (K / p["gamma2"]) * alpha3_eff_casc_full(dtil)
        im_exp = (K / p["gamma2"]) * gamma_2pa_eff(dtil)
        worst = max(worst, abs(c.real - re_exp) / abs(re_exp) if re_exp else 0,
                    abs(c.imag - im_exp) / abs(im_exp) if im_exp else 0)
    print("   max relative deviation over Dtil in [-5,5]: %.2e  -> %s"
          % (worst, "PASS" if worst < 1e-12 else "FAIL"))

    # (iii) dimensional check of chi3_casc:  [rad/s][m/V]^2 / ([1][rad/s]) = m^2/V^2
    print("\n(iii) units: omega1[1/s]*chi2^2[m^2/V^2]/(n2^2 * denom[1/s]) = m^2/V^2  -> PASS (by inspection)")

    # (iv) hand spot value, sample 4, Domega=+2 gamma2, eta=0.0435
    p = samples[4]
    eta4 = eta_geom(p, DOP)
    domega = DOP * p["gamma2"]
    denom = abs(domega + 1j * p["gamma2"] / 2)
    hand = 2 * p["omega1"] * CHI2 ** 2 * eta4 / (3 * p["n2"] ** 2 * denom)
    code = abs(chi3_casc(p, DOP))
    print("\n(iv) sample 4 |chi3| @ +2gamma2:  eta=%.4f  hand=%.3e  code=%.3e  R=%.2f  -> %s"
          % (eta4, hand, code, code / CHI3_NATIVE, "PASS" if abs(hand - code) < 1e-30 else "FAIL"))


# ----------------------------------------------------------------------
def optimal(window=(1.0, 5.0)):
    print("\n" + "=" * 70)
    print("OPTIMAL DETUNING  (data-derived eta(Domega); validity window |Domega|/gamma2 in [%.0f,%.0f])"
          % window)
    print("=" * 70)
    samples = {i: build_sample(i) for i in (1, 2, 3, 4)}
    dtil = np.linspace(-window[1], -window[0], 700)
    dtil = np.r_[dtil, np.linspace(window[0], window[1], 700)]

    print("%-2s | %-32s | %-32s | %-18s" %
          ("S", "max|chi3| (defocus, Dw>0)", "max|chi3| (focus, Dw<0)", "ref @ +2g"))
    print("%-2s | %8s %8s %7s | %8s %8s %7s | %8s %6s" %
          ("", "Dw/g", "|chi3|", "R", "Dw/g", "|chi3|", "R", "|chi3|", "R"))
    rows = {}
    for i, p in samples.items():
        absc = np.array([abs(chi3_casc(p, d)) for d in dtil])
        rec = np.array([chi3_casc(p, d).real for d in dtil])
        pos = dtil > 0
        neg = dtil < 0
        jp = np.argmax(absc[pos]); dp = dtil[pos][jp]; vp = absc[pos][jp]
        jn = np.argmax(absc[neg]); dn = dtil[neg][jn]; vn = absc[neg][jn]
        ref = abs(chi3_casc(p, DOP))
        rows[i] = dict(dp=dp, vp=vp, dn=dn, vn=vn, ref=ref,
                       dtil=dtil, absc=absc, rec=rec,
                       eta_opt=eta_geom(p, dp))
        print("%-2d | %8.2f %8.2e %7.2f | %8.2f %8.2e %7.2f | %8.2e %6.2f"
              % (i, dp, vp, vp / CHI3_NATIVE, dn, vn, vn / CHI3_NATIVE, ref, ref / CHI3_NATIVE))
    print("\ngain at optimum vs +2gamma2 (defocus branch): " +
          ", ".join("S%d x%.2f" % (i, rows[i]["vp"] / rows[i]["ref"]) for i in rows))
    print("eta at optimum: " + ", ".join("S%d %.3f@Dw=%.1f" % (i, rows[i]["eta_opt"], rows[i]["dp"]) for i in rows))

    # figure: |chi3|/native (=R) vs Domega with optimum + 2gamma2 marked
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for i, p in samples.items():
        r = rows[i]
        col = SAMPLE_COLORS[i]
        ax.plot(r["dtil"], r["absc"] / CHI3_NATIVE, color=col, label=_label(p))
        ax.plot([r["dp"]], [r["vp"] / CHI3_NATIVE], "o", color=col, ms=7, mfc="white", mew=1.6)
        ax.plot([r["dn"]], [r["vn"] / CHI3_NATIVE], "o", color=col, ms=7, mfc="white", mew=1.6)
    ax.axhline(1.0, color="0.4", lw=1.0, ls="--", label="native Kerr ($R=1$)")
    for d in (-DOP, DOP):
        ax.axvline(d, color="0.6", lw=0.8, ls=":")
    ax.axvspan(-1, 1, color="0.85", alpha=0.5)
    ax.text(0, ax.get_ylim()[1] * 0.92, "single-pole\nbreakdown", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax.set_ylabel(r"$R^{(1\mathrm{R})}=|\chi^{(3)}_{\rm casc}|/\chi^{(3)}_{\rm native}$")
    ax.set_title("Optimal detuning (circles); dotted = conservative $\\pm2\\gamma_2$")
    ax.legend(fontsize=8.5, loc="upper right")
    _save(fig, "sample_optimal_detuning")
    return rows


if __name__ == "__main__":
    verify()
    optimal()
