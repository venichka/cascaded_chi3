"""Extract SH-only model parameters for the designer's 4 LiNbO3 samples.

Reads the designer's deliverables in ``data/sample1-4/``:
  - design4_*.mat   : per-sample COMSOL output (SH/FM transmission, overlap, area)
  - ref_idx_LN.m    : LiNbO3 Sellmeier table n_e, n_o vs wavelength [um]
  - slide-1 table   : f0(2w), Q2, gamma, H1/H2/H_t, area, overlap(+-2g)

and maps them onto the singly-resonant (SH-only) cascaded-chi3 model of
``notes/Cascaded_chi3_SH_only_note.tex``:

    chi3_casc_1R(Domega) = - omega1 |chi2|^2 eta_geom
                            / [ 6 n2^2 (Domega + i gamma2/2) ]          (eq:chi3-1R-explicit, full complex)

    eta_geom_1R(Domega)  = ( |overlap_complex(f)| / area )^2  *  C_bridge
                           evaluated at SH frequency f = f0 + Domega/(2*pi)   (overlap def, slide image6)

Conventions (note section 2 / 3):
  * Boyd field convention; SHG deg. 1, DFG deg. 2, SPM deg. 3.
  * gamma2 is the SH ENERGY decay rate in rad/s = omega2/Q2 = 2*pi*gamma[THz]
    (gamma[THz] = f0/Q2 is the FWHM linewidth in ordinary frequency).
  * E_z = out-of-plane = extraordinary axis -> d33 / n_e  (slide image7).

Run directly to print the full extraction + verification report.
"""

import os
import re
import numpy as np
import scipy.io as sio
from scipy.interpolate import PchipInterpolator

# ----------------------------------------------------------------------
# locations
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data", "sample1-4"))

# ----------------------------------------------------------------------
# physical constants
C0 = 299_792_458.0          # m/s
TWO_PI = 2.0 * np.pi

# ----------------------------------------------------------------------
# material constants (SUPPLIED BY US, not by the designer).
#   chi2  = 2*d33,  d33(LiNbO3) ~ 25 pm/V at 1064/532 SHG  ->  5e-11 m/V
#           (note baseline; Miller-rule dispersion to ~500 nm SH is a ~+10-20% effect)
#   chi3_native ~ 2e-21 m^2/V^2 (note baseline at 1.55 um; at ~1 um it may be
#           somewhat larger as 2-photon @ ~500 nm = 2.48 eV approaches the ~3.9 eV gap)
CHI2 = 5.0e-11              # m/V
CHI3_NATIVE = 2.0e-21       # m^2/V^2
C_BRIDGE = 1.0             # 2D-normalised-overlap -> 3D Brillouin-M~ bridge, O(1), baseline 1

# unit cell
PERIOD = 400e-9            # m  (slide 1: "Lattice period for all samples is 400 nm")
A_CELL = PERIOD ** 2       # m^2

# operating point default
DOMEGA_OVER_GAMMA = 2.0    # conservative Domega = 2 gamma2 (note default)

# ----------------------------------------------------------------------
# sample <-> file <-> slide-1 table  (table values transcribed from image4.png)
SAMPLES = {
    1: dict(file="design4_alpha9.mat",        f0=601.5751, Q2=3441.6, g_tab=0.1748,
            ov_m2g=0.1000, ov_p2g=0.1200, H1=80,  H2=110, Ht=190, area_tab=4.7675e-14),
    2: dict(file="design4_H1_70.mat",         f0=615.0085, Q2=3486.1, g_tab=0.1764,
            ov_m2g=0.9000, ov_p2g=0.1000, H1=70,  H2=110, Ht=180, area_tab=4.3675e-14),
    3: dict(file="design4_H1_120_H2_100.mat", f0=566.9170, Q2=4641.4, g_tab=0.1221,
            ov_m2g=0.1800, ov_p2g=0.1900, H1=120, H2=100, Ht=220, area_tab=6.2250e-14),
    4: dict(file="design4_H1_140_H2_80.mat",  f0=563.2812, Q2=5632.8, g_tab=0.1000,
            ov_m2g=0.2000, ov_p2g=0.2200, H1=140, H2=80,  Ht=220, area_tab=6.7400e-14),
}


# ----------------------------------------------------------------------
def load_LN_dispersion():
    """Parse data_e / data_o tables (wavelength [um], n) from ref_idx_LN.m."""
    txt = open(os.path.join(DATA, "ref_idx_LN.m")).read()

    def block(name):
        m = re.search(name + r"\s*=\s*\[(.*?)\];", txt, re.S)
        rows = re.findall(r"([0-9.]+)\s*,\s*([0-9.]+)", m.group(1))
        arr = np.array([[float(a), float(b)] for a, b in rows])
        return arr[np.argsort(arr[:, 0])]

    return block("data_e"), block("data_o")


_DATA_E, _DATA_O = load_LN_dispersion()


def n_LN(wavelength_um, axis="e"):
    """Linear interpolation of LiNbO3 index at wavelength [um]."""
    tab = _DATA_E if axis == "e" else _DATA_O
    return float(np.interp(wavelength_um, tab[:, 0], tab[:, 1]))


# ----------------------------------------------------------------------
def _interp_complex(x, xp, fp):
    return np.interp(x, xp, fp.real) + 1j * np.interp(x, xp, fp.imag)


def loaded_Q_from_lineshape(freq, S21_2, f0):
    """Rough loaded Q from the FWHM of the SH transmission resonance (dip)."""
    bg = np.median(np.r_[S21_2[:5], S21_2[-5:]])      # off-resonance baseline
    i0 = np.argmin(S21_2)
    depth = bg - S21_2[i0]
    if depth <= 0:
        return np.nan, bg, S21_2[i0]
    half = bg - depth / 2.0
    # walk left/right from the dip to the half-depth crossings
    left = i0
    while left > 0 and S21_2[left] < half:
        left -= 1
    right = i0
    while right < len(S21_2) - 1 and S21_2[right] < half:
        right += 1
    fwhm = abs(freq[right] - freq[left])
    return (f0 / fwhm if fwhm > 0 else np.nan), bg, S21_2[i0]


def build_sample(idx):
    """Return a dict of fully-derived model parameters for sample ``idx``."""
    s = SAMPLES[idx]
    m = sio.loadmat(os.path.join(DATA, s["file"]))
    freq = np.squeeze(m["freq"]).astype(float)        # SH frequency [THz]
    ffun = np.squeeze(m["freq_fun"]).astype(float)    # FM frequency [THz]
    S21_2 = np.squeeze(m["S21_2"]).astype(float)      # |t_SH|^2
    S21_2_fun = np.squeeze(m["S21_2_fun"]).astype(float)
    overlap = np.squeeze(m["overlap"]).astype(complex)
    symmetry = np.squeeze(m["symmetry"]).astype(complex)
    area = float(np.squeeze(m["area"]))
    order = np.argsort(freq)
    freq, ffun, S21_2, S21_2_fun, overlap, symmetry = (
        freq[order], ffun[order], S21_2[order], S21_2_fun[order],
        overlap[order], symmetry[order])

    f0 = s["f0"]
    Q2 = s["Q2"]
    g_THz = f0 / Q2                                   # FWHM linewidth [THz]
    omega1 = TWO_PI * (f0 / 2.0) * 1e12               # rad/s
    omega2 = TWO_PI * f0 * 1e12                       # rad/s
    gamma2 = omega2 / Q2                              # rad/s (energy decay)

    lam1 = C0 / ((f0 / 2.0) * 1e12) * 1e6             # um
    lam2 = C0 / (f0 * 1e12) * 1e6                     # um
    n1 = n_LN(lam1, "e")
    n2 = n_LN(lam2, "e")

    # normalised overlap |overlap|/area as a function of SH frequency.
    # Shape-preserving (PCHIP) interpolant: smooth C^1 curve through every data
    # point, no overshoot -- removes the piecewise-linear kinks of np.interp when
    # the coarse COMSOL frequency grid is resampled for Domega sweeps.
    norm_ov = np.abs(overlap) / area
    ov_interp = PchipInterpolator(freq, norm_ov, extrapolate=False)

    def eta_at(domega_over_g):
        """eta_geom_1R at detuning Domega = k*gamma2: evaluate overlap at f0 + k*g_THz."""
        f_eval = f0 + domega_over_g * g_THz
        ov = abs(_interp_complex(f_eval, freq, overlap)) / area
        return (ov ** 2) * C_BRIDGE, ov

    eta_p2, ov_p2 = eta_at(+2.0)
    eta_m2, ov_m2 = eta_at(-2.0)
    eta_0, ov_0 = eta_at(0.0)

    Qfit, bg, Tmin = loaded_Q_from_lineshape(freq, S21_2, f0)

    return dict(
        idx=idx, file=s["file"], f0=f0, Q2=Q2, g_THz=g_THz, g_tab=s["g_tab"],
        freq_arr=freq, overlap_arr=overlap, S21_2=S21_2, ov_interp=ov_interp,
        omega1=omega1, omega2=omega2, gamma2=gamma2,
        lam1_um=lam1, lam2_um=lam2, n1=n1, n2=n2, n2_o=n_LN(lam2, "o"),
        area=area, area_tab=s["area_tab"], fill=area / A_CELL,
        H1=s["H1"], H2=s["H2"], Ht=s["Ht"], L_meta=s["Ht"] * 1e-9,
        ov_p2=ov_p2, ov_m2=ov_m2, ov_0=ov_0,
        ov_p2_tab=s["ov_p2g"], ov_m2_tab=s["ov_m2g"],
        eta_p2=eta_p2, eta_m2=eta_m2, eta_0=eta_0,
        norm_ov_min=norm_ov.min(), norm_ov_max=norm_ov.max(),
        Tfun=float(np.median(S21_2_fun)), Tfun_ripple=float(S21_2_fun.max() - S21_2_fun.min()),
        Qfit=Qfit, T_bg=bg, T_min=Tmin,
        Df_GHz=g_THz * 1e3,
    )


def eta_geom(p, domega_over_g):
    """eta_geom_1R(Domega) = (|overlap(f)|/area)^2 * C_bridge, f = f0 + Domega/2pi.

    Domega = domega_over_g * gamma2, so f = f0 + domega_over_g * g_THz.
    """
    f_eval = p["f0"] + domega_over_g * p["g_THz"]
    f_eval = min(max(f_eval, p["freq_arr"].min()), p["freq_arr"].max())  # no extrapolation
    ov = float(p["ov_interp"](f_eval))
    return (ov ** 2) * C_BRIDGE


def chi3_casc(p, domega_over_g, eta=None):
    """Complex chi3_casc_1R [m^2/V^2] at Domega = domega_over_g * gamma2.

    If eta is None it is taken from the overlap data at this detuning.
    """
    if eta is None:
        eta = eta_geom(p, domega_over_g)
    domega = domega_over_g * p["gamma2"]
    denom = domega + 1j * p["gamma2"] / 2.0
    return -2.0 * p["omega1"] * CHI2 ** 2 * eta / (3.0 * p["n2"] ** 2 * denom)  # prefactor 2/3 (2026-10-06 fix; was 1/6)


# ----------------------------------------------------------------------
def report():
    print("#" * 72)
    print("# SH-only sample parameter extraction + verification")
    print("# data dir:", DATA)
    print("# material baseline: |chi2|=%.2e m/V, chi3_native=%.2e m^2/V^2, C_bridge=%.1f"
          % (CHI2, CHI3_NATIVE, C_BRIDGE))
    print("# unit cell: period=%.0f nm, A_cell=%.3e m^2" % (PERIOD * 1e9, A_CELL))
    print("#" * 72)

    ps = {i: build_sample(i) for i in SAMPLES}

    print("\n=== (A) DIRECT PARAMETERS + CONSISTENCY CHECKS ===")
    hdr = ("S", "file", "f0[THz]", "Q2", "g=f0/Q2", "g_tab", "dG%", "lam1um", "lam2um", "n1", "n2")
    print("%-2s %-26s %9s %7s %8s %7s %5s %7s %7s %6s %6s" % hdr)
    for i, p in ps.items():
        dG = 100 * (p["g_THz"] - p["g_tab"]) / p["g_tab"]
        print("%-2d %-26s %9.4f %7.0f %8.5f %7.4f %5.1f %7.4f %7.4f %6.4f %6.4f"
              % (i, p["file"], p["f0"], p["Q2"], p["g_THz"], p["g_tab"], dG,
                 p["lam1_um"], p["lam2_um"], p["n1"], p["n2"]))

    print("\n=== (B) GEOMETRY ===")
    print("%-2s %5s %5s %6s %10s %10s %8s %8s" %
          ("S", "H1", "H2", "Ht=L", "area[m2]", "area_tab", "fill", "Tfun"))
    for i, p in ps.items():
        print("%-2d %5d %5d %6d %10.3e %10.3e %8.3f %8.3f"
              % (i, p["H1"], p["H2"], p["Ht"], p["area"], p["area_tab"], p["fill"], p["Tfun"]))
    print("  (Ht = H1 base + H2 pillars = total LiNbO3 thickness; fill = area/A_cell;")
    print("   Tfun = median FM transmission |t|^2, ripple = %s)"
          % {i: round(ps[i]["Tfun_ripple"], 4) for i in ps})

    print("\n=== (C) OVERLAP DECODE: |overlap_complex|/area  vs slide-1 table ===")
    print("%-2s %10s %10s %6s | %10s %10s %6s" %
          ("S", "calc(-2g)", "tab(-2g)", "d%", "calc(+2g)", "tab(+2g)", "d%"))
    for i, p in ps.items():
        dm = 100 * (p["ov_m2"] - p["ov_m2_tab"]) / p["ov_m2_tab"]
        dp = 100 * (p["ov_p2"] - p["ov_p2_tab"]) / p["ov_p2_tab"]
        flag = "  <-- table likely typo (0.9000 -> 0.0900)" if i == 2 else ""
        print("%-2d %10.4f %10.4f %6.0f | %10.4f %10.4f %6.0f%s"
              % (i, p["ov_m2"], p["ov_m2_tab"], dm, p["ov_p2"], p["ov_p2_tab"], dp, flag))
    print("  overlap @ resonance (Delta=0): %s  (dips by Fano normalisation)"
          % {i: round(ps[i]["ov_0"], 4) for i in ps})

    print("\n=== (D) LINESHAPE Q CROSS-CHECK (FWHM of SH transmission dip) ===")
    print("%-2s %8s %8s %10s %10s" % ("S", "Q_quoted", "Q_fit", "T_bg", "T_min"))
    for i, p in ps.items():
        print("%-2d %8.0f %8.0f %10.4f %10.2e" % (i, p["Q2"], p["Qfit"], p["T_bg"], p["T_min"]))
    print("  (FWHM-Q is a rough Fano-affected estimate; quoted Q is authoritative.)")

    print("\n=== (E) DERIVED MODEL OUTPUTS at Domega = +2 gamma2 (baseline, C_bridge=1) ===")
    print("%-2s %8s %10s %12s %8s %10s" %
          ("S", "eta(+2g)", "|chi3|", "R=|c3|/nat", "Df[GHz]", "|c3|*Df"))
    for i, p in ps.items():
        c3 = chi3_casc(p, +2.0, p["eta_p2"])
        R = abs(c3) / CHI3_NATIVE
        print("%-2d %8.4f %10.3e %12.2f %8.1f %10.3e"
              % (i, p["eta_p2"], abs(c3), R, p["Df_GHz"], abs(c3) * p["Df_GHz"] * 1e9))
    print("  (also reporting -2g branch for the self-focusing/defocusing sign flip:)")
    print("%-2s %8s %12s %12s %12s" % ("S", "eta(-2g)", "Re c3(+2g)", "Re c3(-2g)", "Im c3(+2g)"))
    for i, p in ps.items():
        cp = chi3_casc(p, +2.0, p["eta_p2"])
        cm = chi3_casc(p, -2.0, p["eta_m2"])
        print("%-2d %8.4f %12.3e %12.3e %12.3e" % (i, p["eta_m2"], cp.real, cm.real, cp.imag))

    print("\nRanking by R(+2g): " +
          ", ".join("S%d=%.2f" % (i, abs(chi3_casc(ps[i], 2.0, ps[i]["eta_p2"])) / CHI3_NATIVE)
                    for i in sorted(ps, key=lambda k: -abs(chi3_casc(ps[k], 2.0, ps[k]["eta_p2"])))))
    return ps


if __name__ == "__main__":
    report()
