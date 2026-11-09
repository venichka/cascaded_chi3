"""Material constants for the cascaded chi2:chi2 versus native Kerr comparison.

Reproduces every derived number in constants.md.  Run with
    conda run -n ml_base python compute.py

Conventions (paper, Supplement 1 Sec. 1 and Sec. "Host materials"):
  * Boyd convention, chi2_ijk = 2 d_il (contracted notation).
  * |chi2| = max |chi2_ijk u_i^* v_j v_k| over complex unit vectors u, v.
  * n2 = smallest principal refractive index at the SH (0.775 um).
  * chi_nat = 4 eps0 c n0^2 n2_I / 3, n0 = index at 1.55 um for the pump
    polarisation, n2_I = intensity Kerr index at 1.55 um.
  * F = |chi2|^2 / (n2^2 chi_nat),  R = Q2 eta F / 6,  Q2(R=1) = 6/(eta F).

Refractive-index models (sources in new_refs.bib):
  LiNbO3  Zelmon et al. 1997, congruent, 21 C   (zelmon1997)
  AlN     Pastrnak & Roskovcova 1966            (pastrnak1966)
  GaP     Bond 1965, tabulated, cubic spline    (bond1965)
  GaAs    Skauli et al. 2003 (only for the 4.1-um Miller anchor)  (skauli2003)
  AlGaAs  Gehrsitz et al. 2000 (x, T model; transcription of the
          refractiveindex.info-scripts implementation)            (gehrsitz2000)
Coefficients were taken from the refractiveindex.info database
(github.com/polyanskiy/refractiveindex.info-database), which cites the
same primary papers.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize

EPS0 = 8.8541878128e-12   # F/m
C0 = 299792458.0          # m/s
HC_EV_UM = 1.23984198     # eV um

LAM1, LAM2 = 1.55, 0.775  # um, pump and SH
Q2, ETA = 4000.0, 0.3


# ----------------------------------------------------------------------
# Refractive-index models
# ----------------------------------------------------------------------
def ri_formula1(lam, c):
    """refractiveindex.info 'formula 1': n^2 - 1 = C1 + sum Ci l^2/(l^2 - C_{i+1}^2)."""
    lam2 = lam ** 2
    n2m1 = c[0]
    for i in range(1, len(c), 2):
        n2m1 += c[i] * lam2 / (lam2 - c[i + 1] ** 2)
    return np.sqrt(1.0 + n2m1)


def ri_formula2(lam, c):
    """refractiveindex.info 'formula 2': n^2 - 1 = C1 + sum Ci l^2/(l^2 - C_{i+1})."""
    lam2 = lam ** 2
    n2m1 = c[0]
    for i in range(1, len(c), 2):
        n2m1 += c[i] * lam2 / (lam2 - c[i + 1])
    return np.sqrt(1.0 + n2m1)


# Zelmon, Small, Jundt, JOSA B 14, 3319 (1997); congruent LN, 21 C, 0.4-5 um
ZELMON_E = [0, 2.9804, 0.02047, 0.5981, 0.0666, 8.9543, 416.08]
ZELMON_O = [0, 2.6734, 0.01764, 1.2290, 0.05914, 12.614, 474.60]
# Pastrnak & Roskovcova, phys. stat. sol. 14, K5 (1966); AlN, 0.22-5 um
PASTRNAK_O = [2.1399, 1.3786, 0.1715, 3.861, 15.03]
PASTRNAK_E = [2.0729, 1.6173, 0.1746, 4.139, 15.03]
# Skauli et al., JAP 94, 6447 (2003); GaAs, 22 C, 0.97-17 um
SKAULI = [4.372514, 5.466742, 0.4431307, 0.02429960, 0.8746453,
          1.957522, 36.9166]
# Bond, JAP 36, 1674 (1965); GaP, room temperature (tabulated)
BOND_GAP = np.array([
    [0.5, 3.4595], [0.6, 3.3495], [0.7, 3.2442], [0.8, 3.183], [0.9, 3.143],
    [1.0, 3.1192], [1.1, 3.0981], [1.2, 3.0844], [1.4, 3.0646],
    [1.6, 3.0509], [1.8, 3.0439], [2.0, 3.0379], [2.2, 3.0331],
    [2.4, 3.0296], [2.6, 3.0271], [2.8, 3.0236], [3.0, 3.0215],
    [3.2, 3.0197], [3.4, 3.0181], [3.6, 3.0166], [3.8, 3.0159],
    [4.0, 3.0137]])
_bond_spline = CubicSpline(BOND_GAP[:, 0], BOND_GAP[:, 1])


def n_LN_e(lam):
    return ri_formula2(lam, ZELMON_E)


def n_LN_o(lam):
    return ri_formula2(lam, ZELMON_O)


def n_AlN_o(lam):
    return ri_formula1(lam, PASTRNAK_O)


def n_AlN_e(lam):
    return ri_formula1(lam, PASTRNAK_E)


def n_GaAs(lam):
    return ri_formula1(lam, SKAULI)


def n_GaP(lam):
    return float(_bond_spline(lam))


def n_AlGaAs(x, lam_um, T_C=20.0):
    """Gehrsitz et al., JAP 87, 7825 (2000), Eqs. 11-13 and Tables II, IV.
    Transcribed from refractiveindex.info-scripts
    ('Gehrsitz 2000 - Al_x_Ga_1-x_As.py').  Valid below the band gap."""
    T = T_C + 273.15
    um = 1.239856
    kB = 0.0861708e-3
    VT = 2 * kB * T
    E = 1.0 / lam_um  # photon energy in um^-1 (E[eV]/1.239856)
    coth = lambda z: np.cosh(z) / np.sinh(z)
    EG0, EDeb, ETO, S, STO = 1.5192, 15.9e-3, 33.6e-3, 1.8, 1.1
    EG_GaAs = (EG0 + S * EDeb * (1 - coth(EDeb / VT))
               + STO * ETO * (1 - coth(ETO / VT)))
    A0 = 5.9613 + 7.7178e-4 * T - 0.953e-6 * T ** 2
    A = A0 - 16.159 * x + 43.511 * x ** 2 - 71.317 * x ** 3 \
        + 57.535 * x ** 4 - 17.451 * x ** 5
    E1sq_GaAs = 4.7171 - 3.237e-4 * T - 1.358e-6 * T ** 2
    E2, C2 = 0.724e-3, 1.55e-3
    E3, C3 = 1.331e-3, 2.61e-3
    E0 = (EG_GaAs / um + 1.1308 * x + 0.1436 * x ** 2) ** 2
    C0c = 1.0 / (50.535 - 150.7 * x - 62.209 * x ** 2 + 797.16 * x ** 3
                 - 1125 * x ** 4 + 503.79 * x ** 5)
    E1 = E1sq_GaAs + 11.006 * x - 3.08 * x ** 2
    C1 = 21.5647 + 113.74 * x - 122.5 * x ** 2 + 108.401 * x ** 3 \
        - 47.318 * x ** 4
    R = (1 - x) * C2 / (E2 - E ** 2) + x * C3 / (E3 - E ** 2)
    n2 = A + C0c / (E0 - E ** 2) + C1 / (E1 - E ** 2) + R
    return float(np.sqrt(n2)), float(np.sqrt(E0) * um)  # (n, model E0 in eV)


# ----------------------------------------------------------------------
# Miller's rule  d(2w; w, w) proportional to chi1(2w) chi1(w)^2
# ----------------------------------------------------------------------
def miller(d0, nfun0, lam0, nfun1, lam1):
    """Scale d0 measured for SHG of lam0 (index model nfun0) to SHG of
    lam1 (index model nfun1); chi1 = n^2 - 1."""
    chi = lambda f, l: f(l) ** 2 - 1.0
    num = chi(nfun1, lam1 / 2) * chi(nfun1, lam1) ** 2
    den = chi(nfun0, lam0 / 2) * chi(nfun0, lam0) ** 2
    return d0 * num / den


# ----------------------------------------------------------------------
# Two-band (Sheik-Bahae 1991) dispersion of the bound-electronic n2
# ----------------------------------------------------------------------
def G2(x):
    t = 2 * (1 - 2 * x) ** 1.5 if x < 0.5 else 0.0
    return (-2 + 6 * x - 3 * x ** 2 - x ** 3 - 0.75 * x ** 4
            - 0.75 * x ** 5 + t) / (64 * x ** 6)


def two_band_ratio(lam_a, lam_b, Eg, nfun):
    """n2(lam_b)/n2(lam_a) for a purely bound-electronic Kerr index,
    n2 proportional to G2(hbar w/Eg)/n0^2 (Eg fixed)."""
    xa, xb = HC_EV_UM / lam_a / Eg, HC_EV_UM / lam_b / Eg
    return (G2(xb) / nfun(lam_b) ** 2) / (G2(xa) / nfun(lam_a) ** 2)


# ----------------------------------------------------------------------
# Tensor norm  max |chi_ijk u_i^* v_j v_k|  (Boyd chi = 2 d)
# ----------------------------------------------------------------------
VOIGT = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 2): 3, (2, 1): 3,
         (0, 2): 4, (2, 0): 4, (0, 1): 5, (1, 0): 5}


def chi_from_d(d):
    chi = np.zeros((3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                chi[i, j, k] = 2.0 * d[i, VOIGT[(j, k)]]
    return chi


def tensor_norm(chi, ntrial=200, seed=1):
    """max over complex unit v of |P(v)|, P_i = chi_ijk v_j v_k
    (the optimum u is P/|P|)."""
    rng = np.random.default_rng(seed)

    def negP(p):
        v = p[:3] + 1j * p[3:]
        v = v / np.linalg.norm(v)
        P = np.einsum('ijk,j,k->i', chi, v, v)
        return -np.linalg.norm(P)

    best = 0.0
    for _ in range(ntrial):
        res = minimize(negP, rng.normal(size=6), method='Nelder-Mead',
                       options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
        best = max(best, -res.fun)
    return best


def d_3m(d33, d31, d22):
    d15 = d31  # Kleinman
    return np.array([[0, 0, 0, 0, d15, -d22],
                     [-d22, d22, 0, d15, 0, 0],
                     [d31, d31, d33, 0, 0, 0]], float)


def d_6mm(d33, d31):
    d15 = d31  # Kleinman
    return np.array([[0, 0, 0, 0, d15, 0],
                     [0, 0, 0, d15, 0, 0],
                     [d31, d31, d33, 0, 0, 0]], float)


def d_43m(d14):
    return np.array([[0, 0, 0, d14, 0, 0],
                     [0, 0, 0, 0, d14, 0],
                     [0, 0, 0, 0, 0, d14]], float)


# ----------------------------------------------------------------------
# Figure of merit
# ----------------------------------------------------------------------
def chi_nat(n0, n2I):
    return 4 * EPS0 * C0 * n0 ** 2 * n2I / 3.0


def figure_of_merit(chi2_norm, nSH, n0, n2I):
    cn = chi_nat(n0, n2I)
    F = chi2_norm ** 2 / (nSH ** 2 * cn)
    return cn, F, Q2 * ETA * F / 6.0, 6.0 / (ETA * F)


def main():
    pm = 1e-12
    out = []
    p = out.append

    # ---------------- refractive indices ------------------------------
    nLN_e1, nLN_e2 = n_LN_e(LAM1), n_LN_e(LAM2)
    nLN_o1, nLN_o2 = n_LN_o(LAM1), n_LN_o(LAM2)
    nGaP1, nGaP2 = n_GaP(LAM1), n_GaP(LAM2)
    nAlN_o1, nAlN_o2 = n_AlN_o(LAM1), n_AlN_o(LAM2)
    nAlN_e1, nAlN_e2 = n_AlN_e(LAM1), n_AlN_e(LAM2)
    XAL = 0.20
    nAGA1, E0model = n_AlGaAs(XAL, LAM1)
    nAGA2, _ = n_AlGaAs(XAL, LAM2)
    p("== Refractive indices ==")
    p(f"LiNbO3 (Zelmon, cong.)  n_e(1.55)={nLN_e1:.4f}  n_e(0.775)={nLN_e2:.4f}"
      f"  n_o(1.55)={nLN_o1:.4f}  n_o(0.775)={nLN_o2:.4f}")
    p(f"GaP    (Bond 1965)      n(1.55)={nGaP1:.4f}  n(0.775)={nGaP2:.4f}")
    p(f"AlN    (Pastrnak 1966)  n_e(1.55)={nAlN_e1:.4f}  n_e(0.775)={nAlN_e2:.4f}"
      f"  n_o(1.55)={nAlN_o1:.4f}  n_o(0.775)={nAlN_o2:.4f}")
    p(f"Al0.2Ga0.8As (Gehrsitz, 20 C)  n(1.55)={nAGA1:.4f}  n(0.775)={nAGA2:.4f}"
      f"  [model E0 = {E0model:.4f} eV]")
    # band gap of Al0.2Ga0.8As (Adachi 1985, 300 K): 1.424 + 1.247 x
    Eg_AGA = 1.424 + 1.247 * XAL
    E_SH = HC_EV_UM / LAM2
    p(f"LiNbO3 n_e(1.064) = {n_LN_e(1.064):.4f} (Bache & Schiek used 2.1558)")
    p(f"Al0.2Ga0.8As Gehrsitz gap parameter E0 = {E0model:.3f} eV -> SH is "
      f"{(E0model-HC_EV_UM/LAM2)*1e3:.0f} meV below it")
    p(f"Al0.2Ga0.8As Eg(300 K, Adachi 1985) = {Eg_AGA:.3f} eV "
      f"({HC_EV_UM/Eg_AGA*1e3:.0f} nm); SH photon {E_SH:.3f} eV -> "
      f"{(Eg_AGA-E_SH)*1e3:.0f} meV below gap; half-gap "
      f"{Eg_AGA/2:.3f} eV ({HC_EV_UM/(Eg_AGA/2):.3f} um), pump "
      f"{HC_EV_UM/LAM1:.3f} eV -> {(Eg_AGA/2-HC_EV_UM/LAM1)*1e3:.0f} meV "
      f"below half-gap")

    # ---------------- second-order coefficients ------------------------
    p("\n== Second-order coefficients, Miller's rule to SHG of 1.55 um ==")
    # LiNbO3, Shoji 1997 (congruent): d33 = 19.5 @1.313, 25.2 @1.064
    d33_LN_a = miller(19.5, n_LN_e, 1.313, n_LN_e, LAM1)
    d33_LN_b = miller(25.2, n_LN_e, 1.064, n_LN_e, LAM1)
    d31_LN = miller(3.2, n_LN_o, 1.313, n_LN_o, LAM1)  # d31: o,o -> e approx
    p(f"LiNbO3 d33: 19.5 @1.313 -> {d33_LN_a:.2f} pm/V ; "
      f"25.2 @1.064 -> {d33_LN_b:.2f} pm/V  (Miller, Zelmon n_e)")
    # AlN, Majkic 2017: |d33| = 4.3 +- 0.3 @ 1.030 um, |d31| = |d33|/45
    d33_AlN = miller(4.3, n_AlN_e, 1.030, n_AlN_e, LAM1)
    p(f"AlN d33: 4.3 @1.030 -> {d33_AlN:.2f} pm/V (Miller, Pastrnak n_e)")
    # GaP, Shoji 1997: d36 = 70.6 @ 1.064
    d14_GaP = miller(70.6, n_GaP, 1.064, n_GaP, LAM1)
    p(f"GaP d14: 70.6 @1.064 -> {d14_GaP:.2f} pm/V (Miller, Bond n)")
    d14_GaP_alt = miller(36.8, n_GaP, 1.313, n_GaP, LAM1)
    p(f"   [unverified 36.8 @1.313 would give {d14_GaP_alt:.2f} pm/V]")
    # GaAs anchor, Skauli 2002: d14(GaAs)/d33(LN) = 5.01 +- 0.3 at 4.1 um
    lamM = 4.135
    d33_LN_41 = miller(25.2, n_LN_e, 1.064, n_LN_e, lamM)
    d33_LN_41_r = miller(27.0, n_LN_e, 1.064, n_LN_e, lamM)
    d14_GaAs_41 = 5.01 * d33_LN_41
    p(f"LN d33 @4.135 (Miller from 25.2@1.064) = {d33_LN_41:.2f};"
      f" from 27@1.064 (Skauli's reference) = {d33_LN_41_r:.2f}")
    p(f"GaAs d14 @4.135 = 5.01 x {d33_LN_41:.2f} = {d14_GaAs_41:.1f} pm/V"
      f" (Skauli/Pinguet quote 94 +- 10 with d33(LN)=27 @1.064)")
    nAGA = lambda l: n_AlGaAs(XAL, l)[0]
    d14_AGA = miller(d14_GaAs_41, n_GaAs, lamM, nAGA, LAM1)
    d14_AGA_94 = miller(94.0, n_GaAs, lamM, nAGA, LAM1)
    p(f"Al0.2Ga0.8As d14 @1.55: {d14_AGA:.1f} pm/V (from {d14_GaAs_41:.1f});"
      f" {d14_AGA_94:.1f} pm/V (from 94)  [Miller across x: Ohashi 1993]")
    # Gehrsitz at x=0 vs Skauli at the anchor wavelengths (model consistency)
    p(f"   GaAs index check @4.135/2.0675: Skauli {n_GaAs(lamM):.4f}/"
      f"{n_GaAs(lamM/2):.4f}; Gehrsitz x=0 {n_AlGaAs(0, lamM)[0]:.4f}/"
      f"{n_AlGaAs(0, lamM/2)[0]:.4f}")

    # ---------------- tensor norms -------------------------------------
    p("\n== Tensor norm checks: |chi2| / (2 d_ref) ==")
    for d31 in (+3.0, -3.0):
        for d22 in (0.0, 2.0, -2.0, 3.0):
            nrm = tensor_norm(chi_from_d(d_3m(18.9, d31, d22)), ntrial=40)
            p(f"LN 3m  d33=18.9 d31={d31:+.1f} d22={d22:+.1f}: "
              f"norm/(2 d33) = {nrm/(2*18.9):.4f}")
    for s in (+1, -1):
        nrm = tensor_norm(chi_from_d(d_6mm(4.0, s * 4.0 / 45)), ntrial=40)
        p(f"AlN 6mm d31={s*4.0/45:+.3f}: norm/(2 d33) = {nrm/(2*4.0):.4f}")
    nrm = tensor_norm(chi_from_d(d_43m(1.0)), ntrial=60)
    p(f"43m: norm/(2 d14) = {nrm/2:.5f}  (2/sqrt3 = {2/np.sqrt(3):.5f})")

    # ---------------- Kerr: LN two-band estimates ----------------------
    p("\n== LiNbO3 n2_I: two-band (Sheik-Bahae) scaling, e-pol ==")
    n2_tot_1064 = 38.6e-20  # Bache & Schiek: DeSalvo Z-scan minus cascading
    for Eg in (3.8, 4.0):
        r155 = two_band_ratio(1.064, 1.55, Eg, n_LN_e)
        r205 = two_band_ratio(1.064, 2.05, Eg, n_LN_e)
        p(f"Eg={Eg} eV: n2(1.55)/n2(1.064)={r155:.3f} -> "
          f"{n2_tot_1064*r155:.2e} m^2/W ; n2(2.05)/n2(1.064)={r205:.3f} -> "
          f"{n2_tot_1064*r205:.2e} (Jankowski 2020 quote 2.6e-19 at 2.05)")
        # Raman part held fixed (f_R = 0.5 at 1.064, cf. Bache & Schiek)
        fR = 0.5
        mix = n2_tot_1064 * ((1 - fR) * r155 + fR)
        p(f"      with f_R=0.5 Raman part held fixed: {mix:.2e} m^2/W")

    # ---------------- figure of merit ---------------------------------
    p("\n== Figure of merit (Q2 = 4000, eta = 0.3) ==")
    hosts = [
        # name, d_adopted [pm/V], norm factor, nSH, n0, n2I [m^2/W]
        ("LiNbO3", d33_LN_a, 2.0, nLN_e2, nLN_e1, 1.8e-19),
        ("GaP", d14_GaP, 4 / np.sqrt(3), nGaP2, nGaP1, 1.1e-17),
        ("AlN", d33_AlN, 2.0, nAlN_o2, nAlN_e1, 3.5e-19),
        ("Al0.2Ga0.8As", d14_AGA, 4 / np.sqrt(3), nAGA2, nAGA1, 2.6e-17),
    ]
    p(f"{'host':14s} {'d(pm/V)':>8s} {'|chi2|':>8s} {'nSH':>6s} {'n0':>6s}"
      f" {'n2_I':>9s} {'chi_nat':>10s} {'F':>9s} {'R':>8s} {'Q2(R=1)':>9s}")
    rows = {}
    for name, d, fac, nSH, n0, n2I in hosts:
        c2 = fac * d * pm
        cn, F, R, Qr1 = figure_of_merit(c2, nSH, n0, n2I)
        rows[name] = (d, c2, nSH, n0, n2I, cn, F, R, Qr1)
        p(f"{name:14s} {d:8.2f} {c2/pm:8.1f} {nSH:6.3f} {n0:6.3f}"
          f" {n2I:9.2e} {cn:10.3e} {F:9.3e} {R:8.3g} {Qr1:9.3g}")

    # ---------------- spread ------------------------------------------
    p("\n== Spread: F range from the source spread ==")
    spread = [
        # name, (d_lo, d_hi), (n2I_lo, n2I_hi), fac, nSH, n0
        ("LiNbO3", (d33_LN_a, d33_LN_b), (0.91e-19, 3.2e-19), 2.0,
         nLN_e2, nLN_e1),
        ("GaP", (50.0, d14_GaP), (0.7e-17, 1.4e-17), 4 / np.sqrt(3),
         nGaP2, nGaP1),
        ("AlN", (d33_AlN * 4.0 / 4.3, d33_AlN * 4.6 / 4.3), (2.3e-19, 3.5e-19),
         2.0, nAlN_o2, nAlN_e1),
        ("Al0.2Ga0.8As", (d14_AGA * 0.85, d14_AGA_94 * 1.1), (1.5e-17, 2.6e-17),
         4 / np.sqrt(3), nAGA2, nAGA1),
    ]
    for name, (dlo, dhi), (klo, khi), fac, nSH, n0 in spread:
        Flo = figure_of_merit(fac * dlo * pm, nSH, n0, khi)[1]
        Fhi = figure_of_merit(fac * dhi * pm, nSH, n0, klo)[1]
        p(f"{name:14s} d {dlo:6.1f}-{dhi:6.1f} pm/V, n2_I {klo:.2e}-{khi:.2e}:"
          f"  F {Flo:.3g}-{Fhi:.3g},  R {200*Flo:.3g}-{200*Fhi:.3g},"
          f"  Q2(R=1) {20/Fhi:.3g}-{20/Flo:.3g}")

    # ---------------- LiNbO3 sensitivity grid --------------------------
    p("\n== LiNbO3 sensitivity: R at Q2=4000, eta=0.3 (Q2(R=1) in brackets) ==")
    p("d33 \\ n2_I   " + "  ".join(f"{k:9.2e}" for k in (0.91e-19, 1.8e-19, 3.2e-19)))
    for d in (d33_LN_a, d33_LN_b, 25.2):
        cells = []
        for k in (0.91e-19, 1.8e-19, 3.2e-19):
            _, F, R, Qr1 = figure_of_merit(2 * d * pm, nLN_e2, nLN_e1, k)
            cells.append(f"{R:5.1f} [{Qr1:4.0f}]")
        p(f"{d:6.2f} pm/V  " + "  ".join(cells))

    # ---------------- previous paper table, for comparison -------------
    p("\n== Current Supplement table inputs, recomputed ==")
    for name, c2, cn, nSH in [("LiNbO3", 50, 2e-21, 2.26), ("GaP", 100, 1e-19, 3.1),
                              ("AlN", 10, 2e-21, 2.1), ("AlGaAs", 100, 1e-18, 3.3)]:
        F = (c2 * pm) ** 2 / (nSH ** 2 * cn)
        p(f"{name:8s} F={F:.3g}  R={200*F:.3g}  Q2(R=1)={20/F:.3g}")

    print("\n".join(out))


if __name__ == "__main__":
    main()
