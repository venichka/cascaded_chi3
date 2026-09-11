"""Extract SH-only model parameters for Design 6 (TiO2/LiNbO3 hybrid, single design).

Design 6 (by "Dee & Yajun", 2026-07-06) is a NEW single SH-only design, NOT one of
samples 1-4.  It differs in four load-bearing ways -- see data/design6_parameters.md:

  1. TiO2 pillars (two elliptical holes, etch-angle asymmetry alpha=0 vs 90 deg) on a
     SOLID LiNbO3 slab.  Only the LN is nonlinear.
  2. The COMSOL overlap is now VOLUME-normalised (m^3), not area-normalised (m^2):
        overlap_dimensionless(f) = |overlap_complex(f)| / V_LN
        eta_geom^(1R)(Domega)    = ( |overlap|/V_LN )^2 * C_bridge      (2 powers, as before)
     The ``.mat`` carries ``volume`` instead of ``area``.  This 3-D normalisation is the
     correct match to the model's 3-D M~ (note eq:Mtilde-def), so C_bridge=1 is better
     justified than for samples 1-4 (it no longer bridges a 2-D->3-D dimensional gap).
  3. Optical axis is in-plane along x: E_x drive, n_x = n_e, only d33/E_x for the overlap.
  4. lam1 ~ 1.757 um (SH lam2 ~ 878 nm), not the ~1 um of samples 1-4.

Reuses sample_params for the LiNbO3 Sellmeier table (n_LN) and for the p-keyed
chi3_casc / eta_geom (identical functional form; they read p["ov_interp"], p["n2"], ...).

Run directly to print the extraction + verification report:
    conda run -n ml_base python design6_params.py
"""

import os
import numpy as np
import scipy.io as sio
from scipy.interpolate import PchipInterpolator

import sample_params as sp
from sample_params import (CHI2, CHI3_NATIVE, C0, TWO_PI, C_BRIDGE,
                           n_LN, chi3_casc, eta_geom)

# ----------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
MAT = os.path.join(DATA, "design6.mat")

# ----------------------------------------------------------------------
# pptx-quoted scalars (design6 slides -- verified in data/design6_parameters.md)
Q2      = 924.0            # SH quality factor (image7; our Fano fit = 920, 0.4%)
H_LN    = 150e-9          # m  LiNbO3 slab thickness            (table image124)
H_TIO2  = 246e-9         # m  TiO2 layer thickness
A_LONG  = 176.4e-9       # m  ellipse long axis (full)
A_SHORT = 58.1e-9        # m  ellipse short axis (full)
PERIOD  = 219e-9         # m  single-pillar period p
# supercell holds the two alpha=0/90 pillars -> area = 2 p^2 (verified vs .mat volume)
A_CELL  = 2.0 * PERIOD ** 2
# only the LN is nonlinear -> cascade interaction length = LN thickness
L_META  = H_LN

# Representative TiO2 thin-film index (amorphous/poly, deposition-dependent ~+-0.1).
# Used only for the effective-medium hybrid index estimate (its effect is <1%; see below).
N_TIO2_SH = 2.35          # @ ~878 nm
N_TIO2_FM = 2.28          # @ ~1757 nm

DOMEGA_OVER_GAMMA = 2.0   # conservative default Domega = 2 gamma2


# ----------------------------------------------------------------------
def _hole_area():
    """Area of one elliptical air hole (full-axis convention: a_long,a_short are diameters)."""
    return np.pi * (A_LONG / 2.0) * (A_SHORT / 2.0)


def effective_index(n_LN_val, n_TiO2_val, V_LN):
    """Volume-fill-weighted (effective-medium, zeroth order) index of the TiO2/LN/air stack.

    n_eff = sqrt( <eps> ),  <eps> = (V_LN eps_LN + V_TiO2 eps_TiO2 + V_air) / V_stack.
    The LN is a solid slab (V_LN given); the TiO2 layer has two elliptical air holes.
    Returns (n_eff, dict of volume fractions).
    """
    A_hole2 = 2.0 * _hole_area()
    V_air = A_hole2 * H_TIO2
    V_TiO2 = (A_CELL - A_hole2) * H_TIO2
    V_stack = A_CELL * (H_LN + H_TIO2)
    eps = (V_LN * n_LN_val ** 2 + V_TiO2 * n_TiO2_val ** 2 + V_air * 1.0) / V_stack
    fr = dict(f_LN=V_LN / V_stack, f_TiO2=V_TiO2 / V_stack, f_air=V_air / V_stack)
    return float(np.sqrt(eps)), fr


# ----------------------------------------------------------------------
def _refine_argmin(freq, y):
    """Parabolic-refined location of the minimum of y(freq)."""
    i0 = int(np.argmin(y))
    if 0 < i0 < len(y) - 1:
        y0, y1, y2 = y[i0 - 1], y[i0], y[i0 + 1]
        den = (y0 - 2 * y1 + y2)
        dx = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        return freq[i0] + dx * (freq[i0 + 1] - freq[i0]), i0
    return freq[i0], i0


def loaded_Q_from_lineshape(freq, S21_2, f0):
    """Rough loaded Q from FWHM of the SH transmission dip (grid-limited; quoted Q authoritative)."""
    bg = np.median(np.r_[S21_2[:5], S21_2[-5:]])
    i0 = int(np.argmin(S21_2))
    depth = bg - S21_2[i0]
    if depth <= 0:
        return np.nan, bg, S21_2[i0]
    half = bg - depth / 2.0
    L = i0
    while L > 0 and S21_2[L] < half:
        L -= 1
    R = i0
    while R < len(S21_2) - 1 and S21_2[R] < half:
        R += 1
    fwhm = abs(freq[R] - freq[L])
    return (f0 / fwhm if fwhm > 0 else np.nan), bg, S21_2[i0]


# ----------------------------------------------------------------------
def build_design6():
    """Return a p-dict of fully-derived Design-6 model parameters.

    Keys match build_sample()'s output so the sample_effects / sample_pulsed physics
    (all p-keyed) can be reused directly, plus p["A_cell"], p["volume"], hybrid-index
    diagnostics, and the raw arrays.
    """
    m = sio.loadmat(MAT)
    freq  = np.squeeze(m["freq"]).astype(float)          # SH frequency [THz]
    ffun  = np.squeeze(m["freq_fun"]).astype(float)      # FM frequency [THz]
    S21_2 = np.squeeze(m["S21_2"]).astype(float)
    S21_2_fun = np.squeeze(m["S21_2_fun"]).astype(float)
    overlap = np.squeeze(m["overlap"]).astype(complex)
    symmetry = np.squeeze(m["symmetry"]).astype(complex)
    vol = float(np.squeeze(m["volume"]))                 # LiNbO3 volume [m^3]
    order = np.argsort(freq)
    freq, ffun, S21_2, S21_2_fun, overlap, symmetry = (
        freq[order], ffun[order], S21_2[order], S21_2_fun[order],
        overlap[order], symmetry[order])

    f0, i0 = _refine_argmin(freq, S21_2)                 # SH resonance centre [THz]
    g_THz = f0 / Q2                                       # FWHM linewidth [THz] = usable Kerr BW
    omega1 = TWO_PI * (f0 / 2.0) * 1e12                   # rad/s
    omega2 = TWO_PI * f0 * 1e12
    gamma2 = omega2 / Q2                                  # rad/s (SH energy decay)

    lam1 = C0 / ((f0 / 2.0) * 1e12) * 1e6                 # um
    lam2 = C0 / (f0 * 1e12) * 1e6

    # bare LiNbO3 extraordinary index at the two wavelengths
    n1_LN = n_LN(lam1, "e")
    n2_LN = n_LN(lam2, "e")
    # realistic effective-medium hybrid indices (TiO2/LN/air) -- see effective_index()
    n1_eff, fr1 = effective_index(n1_LN, N_TIO2_FM, vol)
    n2_eff, fr2 = effective_index(n2_LN, N_TIO2_SH, vol)

    # volume-normalised dimensionless overlap -> PCHIP interpolant (shape-preserving)
    norm_ov = np.abs(overlap) / vol
    norm_sym = np.abs(symmetry) / vol
    ov_interp = PchipInterpolator(freq, norm_ov, extrapolate=False)

    Qfit, bg, Tmin = loaded_Q_from_lineshape(freq, S21_2, f0)

    return dict(
        idx=6, name="Design 6", file="design6.mat",
        f0=f0, Q2=Q2, g_THz=g_THz,
        freq_arr=freq, ffun_arr=ffun, overlap_arr=overlap, symmetry_arr=symmetry,
        S21_2=S21_2, S21_2_fun=S21_2_fun, ov_interp=ov_interp,
        norm_ov=norm_ov, norm_sym=norm_sym,
        omega1=omega1, omega2=omega2, gamma2=gamma2,
        lam1_um=lam1, lam2_um=lam2,
        n1=n1_eff, n2=n2_eff, n1_LN=n1_LN, n2_LN=n2_LN,
        n2_o=n_LN(lam2, "o"), fills=fr2, n_TiO2_SH=N_TIO2_SH, n_TiO2_FM=N_TIO2_FM,
        volume=vol, A_cell=A_CELL, fill=vol / (A_CELL * H_LN),
        H_LN=H_LN * 1e9, H_TiO2=H_TIO2 * 1e9, L_meta=L_META,
        period_nm=PERIOD * 1e9, a_long_nm=A_LONG * 1e9, a_short_nm=A_SHORT * 1e9,
        Df_GHz=g_THz * 1e3,
        Tfun=float(np.median(S21_2_fun)),
        Tfun_ripple=float(S21_2_fun.max() - S21_2_fun.min()),
        Qfit=Qfit, T_bg=bg, T_min=Tmin,
    )


# ----------------------------------------------------------------------
def write_spectrum_csv(p=None, path=None):
    """(Re)generate data/design6_spectrum.csv (full-precision 14-column frequency array).

    Single source of truth for the array referenced by data/design6_parameters.md.
    """
    if p is None:
        p = build_design6()
    if path is None:
        path = os.path.join(DATA, "design6_spectrum.csv")
    f = p["freq_arr"]
    ov = p["overlap_arr"]
    sym = p["symmetry_arr"]
    vol = p["volume"]
    dw_over_g = (f - p["f0"]) / p["g_THz"]
    rows = np.column_stack([
        f, p["ffun_arr"], dw_over_g, p["S21_2"],
        np.zeros_like(f),                                  # placeholder if phase not carried
        p["S21_2_fun"], np.zeros_like(f),
        ov.real, ov.imag, p["norm_ov"],
        sym.real, sym.imag, p["norm_sym"], p["norm_ov"] ** 2 * C_BRIDGE])
    # re-read the raw phases straight from the .mat to keep all 14 columns exact
    m = sio.loadmat(MAT)
    order = np.argsort(np.squeeze(m["freq"]).astype(float))
    rows[:, 4] = np.squeeze(m["S21_phase"]).astype(float)[order]
    rows[:, 6] = np.squeeze(m["S21_phase_fun"]).astype(float)[order]
    hdr = ("freq_SH_THz,freq_FM_THz,domega_over_gamma2,S21_2_SH,S21_phase_SH_rad,"
           "S21_2_FM,S21_phase_FM_rad,overlap_re_m3,overlap_im_m3,overlap_abs_over_V,"
           "symmetry_re_m3,symmetry_im_m3,symmetry_abs_over_V,eta_geom_1R")
    np.savetxt(path, rows, delimiter=",", header=hdr, comments="", fmt="%.8g")
    print("wrote %s  (%d rows x %d cols)" % (path, rows.shape[0], rows.shape[1]))


def report():
    print("#" * 72)
    print("# Design 6  (TiO2/LiNbO3 hybrid)  parameter extraction + verification")
    print("# data:", MAT)
    print("# material baseline: |chi2|=%.2e m/V, chi3_native=%.2e m^2/V^2, C_bridge=%.1f"
          % (CHI2, CHI3_NATIVE, C_BRIDGE))
    print("#" * 72)
    p = build_design6()

    print("\n=== (A) DIRECT PARAMETERS ===")
    print("  f0(SH)        = %.4f THz   (argmin-refined; Fano fit 341.29, Q 920 ~ quoted 924)" % p["f0"])
    print("  Q2 (quoted)   = %.0f       Df = gamma_THz = f0/Q2 = %.4f THz (%.1f GHz)"
          % (p["Q2"], p["g_THz"], p["Df_GHz"]))
    print("  lam1(pump)    = %.1f nm    lam2(SH) = %.1f nm" % (p["lam1_um"] * 1e3, p["lam2_um"] * 1e3))
    print("  omega1=%.4e  omega2=%.4e  gamma2=%.4e rad/s" % (p["omega1"], p["omega2"], p["gamma2"]))

    print("\n=== (B) INDICES  (realistic hybrid effective-medium estimate) ===")
    print("  bare LiNbO3 n_e:  n1_LN(lam1)=%.4f   n2_LN(lam2)=%.4f  (n_o(lam2)=%.4f)"
          % (p["n1_LN"], p["n2_LN"], p["n2_o"]))
    print("  TiO2 (assumed) :  n_TiO2(FM)=%.2f    n_TiO2(SH)=%.2f   (amorphous/poly, +-0.1)"
          % (p["n_TiO2_FM"], p["n_TiO2_SH"]))
    print("  fill fractions :  LN=%.3f  TiO2=%.3f  air-holes=%.3f"
          % (p["fills"]["f_LN"], p["fills"]["f_TiO2"], p["fills"]["f_air"]))
    print("  -> hybrid n1_eff=%.4f (%.1f%% vs bare)  n2_eff=%.4f (%.1f%% vs bare)"
          % (p["n1"], 100 * (p["n1"] / p["n1_LN"] - 1), p["n2"], 100 * (p["n2"] / p["n2_LN"] - 1)))
    print("     (TiO2's higher index and the air holes' lower index nearly offset -> <1%% shift;")
    print("      n2_eff is used in chi3_casc; using bare LN n2 changes R by <1%%.)")

    print("\n=== (C) GEOMETRY / VOLUME CROSS-CHECK ===")
    vol_geom = A_CELL * H_LN
    print("  period p=%.0f nm  supercell 2p x p  A_cell=2p^2=%.4e m^2" % (p["period_nm"], p["A_cell"]))
    print("  H_LN=%.0f  H_TiO2=%.0f nm  ellipse a_long/a_short=%.1f/%.1f nm"
          % (p["H_LN"], p["H_TiO2"], p["a_long_nm"], p["a_short_nm"]))
    print("  V_LN(geom=2p^2 H_LN)=%.5e  V_LN(.mat)=%.5e  ratio=%.5f"
          % (vol_geom, p["volume"], vol_geom / p["volume"]))
    print("  L_meta (cascade length = LN thickness) = %.0f nm" % (p["L_meta"] * 1e9))

    print("\n=== (D) OVERLAP DECODE  (dimensionless = |overlap|/volume) ===")
    ov_p2 = float(p["ov_interp"](p["f0"] + 2 * p["g_THz"]))
    ov_m2 = float(p["ov_interp"](p["f0"] - 2 * p["g_THz"]))
    ov_0 = float(p["ov_interp"](p["f0"]))
    print("  |overlap|/V range = [%.4f, %.4f]" % (p["norm_ov"].min(), p["norm_ov"].max()))
    print("  at -2g: calc=%.4f  (pptx tab 0.41)   at +2g: calc=%.4f  (pptx tab 0.44)" % (ov_m2, ov_p2))
    print("  at   0: calc=%.4f  (Fano dip at resonance)" % ov_0)
    print("  |symmetry|/V range = [%.4f, %.4f]  (undocumented 'mode symmetry' diagnostic)"
          % (p["norm_sym"].min(), p["norm_sym"].max()))

    print("\n=== (E) SH-ONLY PREMISE + Q CROSS-CHECK ===")
    print("  FM transmission median=%.4f  ripple=%.4f  (flat -> SH-only holds)"
          % (p["Tfun"], p["Tfun_ripple"]))
    print("  SH dip: bg=%.4f  min=%.4f (deep, near-critical)   Q_fit(FWHM,rough)=%.0f vs quoted %.0f"
          % (p["T_bg"], p["T_min"], p["Qfit"], p["Q2"]))

    print("\n=== (F) DERIVED chi3_casc at Domega=+-2 gamma2  (baseline C_bridge=1) ===")
    for s, dop in (("+2g", +2.0), ("-2g", -2.0)):
        e = eta_geom(p, dop)
        c3 = chi3_casc(p, dop)
        print("  %s: eta=%.4f  |chi3|=%.3e  R=%.2f  Re=%+.3e  Im=%+.3e"
              % (s, e, abs(c3), abs(c3) / CHI3_NATIVE, c3.real, c3.imag))
    c3p = chi3_casc(p, +2.0)
    print("  |chi3|*Df=%.3e m^2 Hz/V^2   (sign flip: Re>0 at -2g, Re<0 at +2g)"
          % (abs(c3p) * p["Df_GHz"] * 1e9))
    return p


if __name__ == "__main__":
    p = report()
    print()
    write_spectrum_csv(p)
