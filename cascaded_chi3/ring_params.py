"""Cascaded chi2:chi2 -> chi3_eff in a lithium-niobate microring resonator.

Parameter builder and closed-form observables for the *travelling-wave
microring* platform (Lipson-type thin-film LiNbO3 rings), as opposed to the
metasurface platform of the main note.

Conventions
-----------
We work in the photon-number (Hamiltonian) convention standard in the chi2
microring literature, which is the one used by Lu et al., Optica 7, 1654
(2020) [arXiv:2007.07411]:

    H/hbar = w1 a^dag a + w2 b^dag b + g ( (a^dag)^2 b + a^2 b^dag )

so that, in a frame rotating at w_p for a and 2 w_p for b,

    da/dt = -(i d1 + k1) a - 2 i g conj(a) b + sqrt(2 k1x) s_in
    db/dt = -(i d2 + k2) b - i g a^2

with
    d1 = w1 - w_p,          d2 = w2 - 2 w_p,
    k_n = w_n / (2 Q_n)     (HWHM *amplitude* decay rate),
    gamma_n = 2 k_n = w_n/Q_n  (the energy decay rate used in the main note).

The main note's energy-normalised amplitudes are a_n^note = sqrt(hbar w_n) a_n,
which maps the note's modal coupling beta onto g via

    g = beta * sqrt(hbar w1 / 2).

Detuning sign convention (identical to the main note):

    Domega = 2 w1 - w2   ->   d2 = -Domega + 2 d1.

Kerr sign convention: a per-photon coefficient U > 0 means SELF-FOCUSING,
i.e. the FM resonance shifts DOWN by U per photon,

    w1_eff(n) = w1 - U n.

Native Kerr is always self-focusing in LiNbO3 (U_native > 0). The cascaded
contribution is

    U_casc(Domega) = -2 g^2 Domega / (Domega^2 + k2^2),

negative (self-DEfocusing) for Domega > 0 and positive for Domega < 0: this
is the sign-tunability that is the signature of the cascade.

Everything here is SI, angular frequencies in rad/s.
"""

from __future__ import annotations

from dataclasses import dataclass, replace, field
import numpy as np

# ---------------------------------------------------------------- constants
HBAR = 1.054_571_817e-34      # J s
C0 = 2.997_924_58e8           # m/s
EPS0 = 8.854_187_8128e-12     # F/m
TWOPI = 2.0 * np.pi


# ============================================================== ring params
@dataclass
class RingParams:
    """A doubly-resonant chi2 microring.

    Attributes
    ----------
    label        : human-readable design name
    lambda1      : FM (pump) vacuum wavelength [m]
    Q1_int, Q1_load : intrinsic / loaded quality factor of the FM mode
    Q2_int, Q2_load : intrinsic / loaded quality factor of the SH mode
    R_ring       : ring radius [m]
    A_eff        : effective modal area of the FM mode [m^2]
    n1, n2       : phase indices at FM / SH
    ng1, ng2     : group indices at FM / SH
    n2_kerr      : native nonlinear index n2 [m^2/W] of the host at the FM
    chi2         : |chi^(2)| = 2 d_eff [m/V] (annotation / eta_geom bridge)
    g            : chi2 photon-photon coupling rate [rad/s].  If None it is
                   computed from `eta_norm_ref` via `g_from_eta_norm`.
    eta_norm_ref : measured normalised SHG efficiency [1/W] used to pin g
                   (P_SH/P_FM^2 on double resonance).  Optional.
    """

    label: str = "PPLN ring (target)"
    lambda1: float = 1.55e-6
    Q1_int: float = 2.0e6
    Q1_load: float = 1.0e6
    Q2_int: float = 6.0e5
    Q2_load: float = 3.0e5
    R_ring: float = 50e-6
    A_eff: float = 1.0e-12
    n1: float = 2.21
    n2: float = 2.26
    ng1: float = 2.30
    ng2: float = 2.45
    n2_kerr: float = 1.8e-19
    chi2: float = 5.0e-11
    g: float | None = 1.42e6 * TWOPI
    eta_norm_ref: float | None = None
    notes: str = ""

    # ---------------------------------------------------------- geometry
    @property
    def omega1(self) -> float:
        return TWOPI * C0 / self.lambda1

    @property
    def omega2(self) -> float:
        return 2.0 * self.omega1

    @property
    def lambda2(self) -> float:
        return self.lambda1 / 2.0

    @property
    def L_ring(self) -> float:
        return TWOPI * self.R_ring

    @property
    def V_eff(self) -> float:
        """FM mode volume, A_eff * circumference."""
        return self.A_eff * self.L_ring

    @property
    def FSR(self) -> float:
        """Free spectral range [rad/s]."""
        return TWOPI * C0 / (self.ng1 * self.L_ring)

    # ------------------------------------------------------------- rates
    @property
    def kappa1(self) -> float:
        """FM total (loaded) HWHM amplitude decay rate [rad/s]."""
        return self.omega1 / (2.0 * self.Q1_load)

    @property
    def kappa1_0(self) -> float:
        return self.omega1 / (2.0 * self.Q1_int)

    @property
    def kappa1_ex(self) -> float:
        return self.kappa1 - self.kappa1_0

    @property
    def kappa2(self) -> float:
        return self.omega2 / (2.0 * self.Q2_load)

    @property
    def kappa2_0(self) -> float:
        return self.omega2 / (2.0 * self.Q2_int)

    @property
    def kappa2_ex(self) -> float:
        return self.kappa2 - self.kappa2_0

    @property
    def gamma1(self) -> float:
        """FM energy decay rate = FWHM linewidth [rad/s] (note convention)."""
        return 2.0 * self.kappa1

    @property
    def gamma2(self) -> float:
        return 2.0 * self.kappa2

    @property
    def coupling_eta1(self) -> float:
        """External coupling fraction kappa_ex/kappa at the FM (1/2 = critical)."""
        return self.kappa1_ex / self.kappa1

    # ------------------------------------------------------------ chi2 g
    @property
    def g_rate(self) -> float:
        """chi2 coupling rate [rad/s], from `g` or inferred from eta_norm_ref."""
        if self.g is not None:
            return self.g
        if self.eta_norm_ref is None:
            raise ValueError("need either g or eta_norm_ref")
        return g_from_eta_norm(self, self.eta_norm_ref)

    # ------------------------------------------------------ native Kerr
    @property
    def U_native(self) -> float:
        """Native Kerr per-photon resonance shift [rad/s per photon], > 0.

        w1_eff = w1 - U_native * n.  Derived from the round-trip nonlinear
        phase: gamma_NL = w1 n2 /(c A_eff), P_circ = U_cav c/(ng L), and
        dphi/dw = L ng / c, giving U = hbar w1^2 c n2 / (ng^2 V_eff).
        """
        return (HBAR * self.omega1**2 * C0 * self.n2_kerr
                / (self.ng1**2 * self.V_eff))

    @property
    def chi3_native(self) -> float:
        """Native chi3 [m^2/V^2] equivalent to n2_kerr (Boyd convention)."""
        return (4.0 / 3.0) * EPS0 * C0 * self.n1**2 * self.n2_kerr

    # --------------------------------------------------- adiabatic checks
    @property
    def adiabatic_ratio(self) -> float:
        """kappa2/kappa1 = 2 Q1/Q2.  Adiabatic elimination needs this >~ 1."""
        return self.kappa2 / self.kappa1

    def replace(self, **kw) -> "RingParams":
        return replace(self, **kw)


# ================================================== cascaded coefficients
def Domega_of(p: RingParams, x: float | np.ndarray) -> np.ndarray:
    """Convert a detuning given in units of gamma2 into rad/s."""
    return np.asarray(x, dtype=float) * p.gamma2


def U_casc(p: RingParams, Domega) -> np.ndarray:
    """Cascaded per-photon Kerr shift [rad/s/photon].

    U_casc = -2 g^2 Domega / (Domega^2 + kappa2^2).
    Negative (self-defocusing) for Domega > 0.  Exact (no |Domega|>>k2
    expansion).
    """
    D = np.asarray(Domega, dtype=float)
    g = p.g_rate
    return -2.0 * g**2 * D / (D**2 + p.kappa2**2)


def kappa_casc(p: RingParams, Domega) -> np.ndarray:
    """Cascaded two-photon-like loss, per photon [rad/s/photon], >= 0.

    kappa_eff = kappa1 + kappa_casc * n.
    """
    D = np.asarray(Domega, dtype=float)
    g = p.g_rate
    return 2.0 * g**2 * p.kappa2 / (D**2 + p.kappa2**2)


def loss_to_kerr(p: RingParams, Domega) -> np.ndarray:
    """rho = kappa_casc/|U_casc| = kappa2/|Domega| = gamma2/(2|Domega|).

    Diverges as Domega -> 0 (pure loss, no Kerr), which is the physically
    correct limit.
    """
    D = np.abs(np.asarray(Domega, dtype=float))
    with np.errstate(divide="ignore"):
        return np.where(D == 0, np.inf, p.kappa2 / np.where(D == 0, 1.0, D))


def R_ratio(p: RingParams, Domega) -> np.ndarray:
    """Cascaded enhancement ratio R = |U_casc| / U_native."""
    return np.abs(U_casc(p, Domega)) / p.U_native


def n2_casc_eff(p: RingParams, Domega) -> np.ndarray:
    """Equivalent nonlinear index [m^2/W] of the cascaded response (signed).

    Inverts U = hbar w1^2 c n2 /(ng^2 V).  Sign: positive = self-focusing,
    so it carries the opposite sign of the resonance shift, i.e. the same
    sign as U_casc.
    """
    U = U_casc(p, Domega)
    return U * p.ng1**2 * p.V_eff / (HBAR * p.omega1**2 * C0)


def chi3_casc_eff(p: RingParams, Domega) -> np.ndarray:
    """Equivalent chi3 [m^2/V^2] of the cascaded response (signed)."""
    return (4.0 / 3.0) * EPS0 * C0 * p.n1**2 * n2_casc_eff(p, Domega)


def eta_geom_implied(p: RingParams, Domega) -> np.ndarray:
    """Geometric overlap eta_geom implied by the main note's boxed formula.

        chi3_casc = - w1 |chi2|^2 eta_geom / (6 n2^2 Domega)

    Inverting gives eta_geom = 6 n2^2 |Domega| |chi3_casc| / (w1 |chi2|^2).
    A sanity bound: eta_geom must be <= 1.
    """
    D = np.abs(np.asarray(Domega, dtype=float))
    return (6.0 * p.n2**2 * D * np.abs(chi3_casc_eff(p, Domega))
            / (p.omega1 * p.chi2**2))


# ======================================================= SHG efficiency <-> g
def eta_norm_from_g(p: RingParams, g: float | None = None) -> float:
    """Normalised on-chip SHG efficiency P_SH/P_FM^2 [1/W] on double resonance.

    Derived from the same TCMT (undepleted pump, d1 = d2 = 0):

        eta = 16 g^2 k2x k1x^2 / (hbar w1 k2^2 k1^4).

    This reproduces Eq. (4) of Lu et al., Optica 7, 1654 (2020) exactly:
    their P_b,out = g^2 P^2/(hbar w_p^2) * w_b * [2 k_b1/(d_b^2+k_b^2)]
                    * [2 k_a1/(d_a^2+k_a^2)]^2, with w_b = 2 w_p.
    """
    g = p.g_rate if g is None else g
    return (16.0 * g**2 * p.kappa2_ex * p.kappa1_ex**2
            / (HBAR * p.omega1 * p.kappa2**2 * p.kappa1**4))


def g_from_eta_norm(p: RingParams, eta_norm: float) -> float:
    """Invert `eta_norm_from_g` for g [rad/s]."""
    return float(np.sqrt(eta_norm * HBAR * p.omega1 * p.kappa2**2 * p.kappa1**4
                         / (16.0 * p.kappa2_ex * p.kappa1_ex**2)))


def g_geometric(p: RingParams, overlap: float = 1.0, qpm: bool = True) -> float:
    """First-principles estimate of g from chi2 and the mode volume.

    From the main note's beta = (eps0 w1/4) * Int chi2 f2* f1 f1 with
    |f_n| ~ 1/sqrt(eps0 n_n^2 V) and g = beta sqrt(hbar w1/2):

        g ~ (w1 |chi2| / (4 sqrt(eps0) n2 n1^2 sqrt(V))) * sqrt(hbar w1/2)
            * overlap * (2/pi if quasi-phase-matched)

    NOTE this inherits the O(1) prefactor ambiguity flagged in note Sec. 10
    (eps0 w1/4 vs eps0 w1/2 -> factor 2 in g).  It is a cross-check, not the
    anchor; the anchor is `g_from_eta_norm` on a measured efficiency.
    """
    fac = (2.0 / np.pi) if qpm else 1.0
    return (p.omega1 * p.chi2 * np.sqrt(HBAR * p.omega1 / 2.0)
            / (4.0 * np.sqrt(EPS0) * p.n2 * p.n1**2 * np.sqrt(p.V_eff))
            * overlap * fac)


# ================================================== power / photon numbers
def n_from_power(p: RingParams, P_in: float, delta1: float = 0.0) -> float:
    """Intracavity FM photon number for input power P_in [W] (linear cavity)."""
    return (2.0 * p.kappa1_ex * P_in / (HBAR * p.omega1)
            / (delta1**2 + p.kappa1**2))


def power_from_n(p: RingParams, n: float, delta1: float = 0.0) -> float:
    """Inverse of `n_from_power` (linear cavity)."""
    return (n * HBAR * p.omega1 * (delta1**2 + p.kappa1**2)
            / (2.0 * p.kappa1_ex))


# ============================================ bistability / MI thresholds
def bistability_ok(U: float, kc: float) -> bool:
    """Does a Kerr + nonlinear-loss cavity admit bistability at all?

    The steady state is P(n) = n [ (d1 - U n)^2 + (k1 + kc n)^2 ].  Two
    turning points require 4 (k1 kc - d1 U)^2 > 3 (U^2 + kc^2)(d1^2 + k1^2)
    for some d1, which as |d1| -> inf reduces to 4 U^2 > 3 (U^2 + kc^2),
    i.e. rho = kc/|U| < 1/sqrt(3).
    """
    if U == 0.0:
        return False
    return abs(kc / U) < 1.0 / np.sqrt(3.0)


def bistability_delta_min(p: RingParams, U: float, kc: float) -> float:
    """Smallest |detuning| (in rad/s) that supports bistability.

    Solving 4 (rho - s)^2 = 3 (1 + rho^2)(s^2 + 1) for s = |d1|/k1 with the
    pump on the side the resonance moves toward, sign(d1) = +sign(U)
    (delta1_eff = delta1 - U n reaches zero only then):

        (4 - 3(1+rho^2)) s^2 - 8 rho s + (4 rho^2 - 3(1+rho^2)) = 0
    """
    if not bistability_ok(U, kc):
        return np.inf
    rho = abs(kc / U)
    A = 4.0 - 3.0 * (1.0 + rho**2)
    B = -8.0 * rho
    Cc = 4.0 * rho**2 - 3.0 * (1.0 + rho**2)
    disc = B**2 - 4.0 * A * Cc
    s = (-B + np.sqrt(disc)) / (2.0 * A)
    return s * p.kappa1


def bistability_threshold_power(p: RingParams, U: float, kc: float = 0.0,
                                delta1: float | None = None) -> float:
    """Input power [W] at the onset of bistability, FIXED-COEFFICIENT model.

    WARNING: this treats U and kc as constants.  In a real pump scan the
    pump-frame mismatch is Delta_eff = Domega - 2 delta1, so U and kc move
    with delta1; use `bistability_threshold_exact` for the physical
    threshold at a given cold-cavity Domega.  Kept because the closed-form
    cusp algebra it implements is what the analytic condition
    rho < 1/sqrt(3) refers to, and because it is exact in the limit
    |delta1| << |Domega|.

    With kc = 0 and delta1 = -sign(U) sqrt(3) k1 the closed form is
        P_th = 4 hbar w1 k1^3 / (3 sqrt(3) |U| k1x).
    With kc > 0 the turning points are found from the quadratic
        3(U^2+kc^2) n^2 + 4(k1 kc - d1 U) n + (d1^2 + k1^2) = 0
    at the minimum-|delta1| that still admits a double root.
    """
    if not bistability_ok(U, kc):
        return np.inf
    if delta1 is None:
        # The pump must sit on the side the resonance moves toward:
        # delta1_eff = delta1 - U n reaches zero only if sign(delta1) = sign(U).
        delta1 = np.sign(U) * bistability_delta_min(p, U, kc)
    a = 3.0 * (U**2 + kc**2)
    b = 4.0 * (p.kappa1 * kc - delta1 * U)
    c = delta1**2 + p.kappa1**2
    disc = max(b**2 - 4.0 * a * c, 0.0)   # == 0 exactly at threshold detuning
    n_turn = (-b - np.sqrt(disc)) / (2.0 * a)   # lower turning point
    if n_turn <= 0:
        n_turn = (-b + np.sqrt(disc)) / (2.0 * a)
    if n_turn <= 0:
        return np.inf
    P = (n_turn * HBAR * p.omega1 / (2.0 * p.kappa1_ex)
         * ((delta1 - U * n_turn)**2 + (p.kappa1 + kc * n_turn)**2))
    return P


def mi_threshold_power(p: RingParams, U: float, kc: float = 0.0) -> float:
    """Parametric (modulation-instability / comb) threshold input power [W].

    Gain-equals-loss estimate: |U| n = kappa1 + kc n, so
        n_th = kappa1 / (|U| - kc),
    finite only for rho = kc/|U| < 1.  Evaluated at the detuning that keeps
    the pump mode on effective resonance (delta1_eff = 0), where the input
    power is n_th hbar w1 kappa1^2/(2 kappa1x).
    """
    if abs(U) <= kc:
        return np.inf
    n_th = p.kappa1 / (abs(U) - kc)
    return n_th * HBAR * p.omega1 * p.kappa1**2 / (2.0 * p.kappa1_ex)


def switching_photons(p: RingParams, U: float) -> float:
    """Intracavity photons for a one-linewidth resonance shift (switching)."""
    return p.kappa1 / abs(U)


def saturation_photons(p: RingParams, Domega: float) -> float:
    """Analytic ceiling: |beta a1| = |Domega| maps to sqrt(2) g sqrt(n) = |Domega|.

    n_sat = Domega^2 / (2 g^2), at which |U_casc| n_sat -> |Domega| in the
    far-detuned limit.  Checked numerically in ring_effects.stability_ceiling.
    """
    return Domega**2 / (2.0 * p.g_rate**2)


# ================================================== quantum figure of merit
def anharmonicity_fom(p: RingParams, Domega) -> np.ndarray:
    """Dispersive (cascaded) single-photon anharmonicity |U_casc| / kappa1."""
    return np.abs(U_casc(p, Domega)) / p.kappa1


def resonant_fom(p: RingParams) -> float:
    """Resonant chi2 figure of merit g/kappa1 (Lu et al. definition)."""
    return p.g_rate / p.kappa1


# ============================================================ thermal model
@dataclass
class ThermalModel:
    """Simple lumped thermal model for the discrimination analysis."""
    dn_dT: float = 3.3e-5          # 1/K, LiNbO3 extraordinary thermo-optic
    R_th: float = 3.0e4            # K/W, microring-on-oxide thermal resistance
    absorbed_fraction: float = 0.3  # of the dropped power that becomes heat
    tau_th: float = 1.0e-6         # s, thermal response time

    def dw_dT(self, p: RingParams) -> float:
        """Thermal tuning of the FM resonance [rad/s per K] (negative)."""
        return -p.omega1 * self.dn_dT / p.ng1

    def shift(self, p: RingParams, P_in: float) -> float:
        """Thermal resonance shift magnitude [rad/s] at input power P_in."""
        P_drop = P_in * 4.0 * p.kappa1_ex * p.kappa1_0 / p.kappa1**2
        P_abs = P_drop * self.absorbed_fraction
        return abs(self.dw_dT(p)) * P_abs * self.R_th


# ============================================================ design points
def design_points() -> dict[str, RingParams]:
    """The three reference designs used throughout the note."""
    # (1) Demonstrated PPLN microring: Lu et al., Optica 7, 1654 (2020).
    lu = RingParams(
        label="PPLN, demonstrated (Lu 2020)",
        lambda1=1.560e-6,
        Q1_int=1.8e6, Q1_load=5.4e5,
        Q2_int=5.8e5, Q2_load=4.5e5,
        R_ring=70e-6, A_eff=1.0e-12,
        g=1.2e6 * TWOPI,
        eta_norm_ref=5.0e4,          # 5,000,000 %/W
        notes="z-cut, radially poled, TM0->TM0, Lambda = 2.95 um",
    )
    # (2) Same poling technology on a Lipson-type high-Q ring.  g is scaled
    #     from (1) as 1/sqrt(V) at fixed A_eff, i.e. sqrt(70/50).
    g_target = lu.g_rate * np.sqrt(lu.R_ring / 50e-6)
    ppln = RingParams(
        label="PPLN, target platform",
        lambda1=1.55e-6,
        Q1_int=2.0e6, Q1_load=1.0e6,   # critical coupling
        Q2_int=6.0e5, Q2_load=3.0e5,
        R_ring=50e-6, A_eff=1.0e-12,
        g=g_target,
        notes="Q1 loaded 1e6 (intrinsic 2e6, critically coupled); Q2 = 3e5",
    )
    # (3) No poling: modal / cyclic phase matching.  Penalty is d31-vs-d33
    #     plus the transverse mode-order mismatch; ~12x down in g, i.e. a
    #     normalised efficiency ~1e4 %/W rather than ~5e6 %/W.
    modal = ppln.replace(
        label="Modal / cyclic PM, no poling",
        g=g_target / 12.0,
        notes="d31 or mode-order-mismatched d33; g ~ 1/12 of the poled ring",
    )
    return {"lu": lu, "ppln": ppln, "modal": modal}


# ================================================================= verify
def _fmt(x, unit="", sig=3):
    return f"{x:.{sig}g} {unit}".strip()


def verify(p: RingParams | None = None) -> None:
    """Print a self-check report.  Run as `python ring_params.py`."""
    dps = design_points()
    lu, ppln, modal = dps["lu"], dps["ppln"], dps["modal"]
    p = ppln if p is None else p

    line = "=" * 78
    print(line)
    print("ring_params.py verification report")
    print(line)

    # ---- 1. eta_norm <-> g round trip -------------------------------
    print("\n[1] SHG efficiency <-> g inversion (round trip)")
    for q in (lu, ppln, modal):
        eta = eta_norm_from_g(q)
        g_back = g_from_eta_norm(q, eta)
        rel = abs(g_back - q.g_rate) / q.g_rate
        print(f"    {q.label:32s} eta = {eta:10.3e} /W  "
              f"({eta*100:9.3e} %/W)   round-trip rel.err = {rel:.2e}")
        assert rel < 1e-12, "eta<->g inversion broken"

    # ---- 2. cross-check against the measured Lu et al. device --------
    print("\n[2] Cross-check vs Lu et al., Optica 7, 1654 (2020)")
    print(f"    quoted g/2pi           = 1.2 MHz")
    print(f"    quoted eta_norm        = 5.0e6 %/W  (= 5.0e4 /W)")
    eta_pred = eta_norm_from_g(lu)
    print(f"    eta predicted from g   = {eta_pred*100:.3e} %/W")
    print(f"    ratio predicted/quoted = {eta_pred/5.0e4:.2f}")
    g_from_eta = g_from_eta_norm(lu, 5.0e4)
    print(f"    g inferred from eta    = {g_from_eta/TWOPI/1e6:.3f} MHz")
    print("    -> the two anchors differ by ~2x in g (4x in eta); the paper")
    print("       fits g from the full lineshape including a residual SH")
    print("       detuning d2 != 0, which lowers eta at fixed g.  We carry")
    print("       g/2pi = 1.2 MHz (their fitted value) and treat the spread")
    print("       as the O(1) uncertainty on the absolute scale.")
    print(f"    kappa_a/2pi (loaded)   = {lu.kappa1/TWOPI/1e6:.1f} MHz"
          f"   [paper: 184 MHz]")
    print(f"    kappa_a0/2pi           = {lu.kappa1_0/TWOPI/1e6:.1f} MHz"
          f"   [paper: 55.4 MHz]")
    print(f"    kappa_b0/2pi           = {lu.kappa2_0/TWOPI/1e6:.1f} MHz"
          f"   [paper: 343.8 MHz]")

    # ---- 3. geometric estimate of g ---------------------------------
    print("\n[3] First-principles g vs anchored g (O(1) prefactor check)")
    for q in (lu, ppln):
        g_perf = g_geometric(q, overlap=1.0, qpm=True)
        print(f"    {q.label:32s} g_geom(perfect overlap, QPM)/2pi = "
              f"{g_perf/TWOPI/1e6:.2f} MHz   vs  anchored "
              f"{q.g_rate/TWOPI/1e6:.2f} MHz  -> implied overlap^2 = "
              f"{(q.g_rate/g_perf)**2:.2f}")

    # ---- 4. native Kerr, two independent routes ---------------------
    print("\n[4] Native Kerr per-photon shift, two routes")
    U1 = p.U_native
    # route B: from chi3 through the note's alpha3 with zeta1 = 1 and the
    # standard n2 <-> chi3 conversion; must agree by construction of chi3.
    n2_back = (p.chi3_native * 3.0) / (4.0 * EPS0 * C0 * p.n1**2)
    print(f"    U_native/2pi           = {U1/TWOPI:.3f} Hz per photon")
    print(f"    chi3_native            = {p.chi3_native:.3e} m^2/V^2")
    print(f"    n2 round trip          = {n2_back:.4e} m^2/W "
          f"(input {p.n2_kerr:.4e})")
    assert abs(n2_back - p.n2_kerr) / p.n2_kerr < 1e-12

    # ---- 5. cascaded coefficients at the reference operating point ---
    print("\n[5] Cascaded response at Domega = 2 gamma2 (reference point)")
    for q in (lu, ppln, modal):
        D = 2.0 * q.gamma2
        U = U_casc(q, D)
        kc = kappa_casc(q, D)
        print(f"    {q.label:32s}")
        print(f"        gamma2/2pi = {q.gamma2/TWOPI/1e9:8.3f} GHz   "
              f"Domega/2pi = {D/TWOPI/1e9:8.3f} GHz")
        print(f"        U_casc/2pi = {U/TWOPI:10.1f} Hz/photon   "
              f"U_native/2pi = {q.U_native/TWOPI:6.3f} Hz/photon")
        print(f"        R          = {float(R_ratio(q, D)):10.1f}   "
              f"rho = kc/|U| = {float(loss_to_kerr(q, D)):.3f}")
        print(f"        n2_eff     = {float(n2_casc_eff(q, D)):.3e} m^2/W   "
              f"chi3_eff = {float(chi3_casc_eff(q, D)):.3e} m^2/V^2")
        print(f"        eta_geom implied by note Eq. (chi3casc_explicit) = "
              f"{float(eta_geom_implied(q, D)):.3f}   (must be <= 1)")
        print(f"        adiabatic ratio kappa2/kappa1 = "
              f"{q.adiabatic_ratio:.2f}  (= 2 Q1/Q2)")

    # ---- 6. analytic consistency of the exact vs far-detuned forms ---
    print("\n[6] Exact vs far-detuned cascaded Kerr")
    for x in (0.5, 1.0, 2.0, 5.0):
        D = x * p.gamma2
        exact = abs(float(U_casc(p, D)))
        approx = 2.0 * p.g_rate**2 / D
        print(f"    Domega = {x:4.1f} gamma2:  exact/approx = "
              f"{exact/approx:.4f}   (expect 1/(1+1/(4x^2)) = "
              f"{1.0/(1.0+1.0/(4*x*x)):.4f})")
        assert abs(exact / approx - 1.0 / (1.0 + 1.0 / (4 * x * x))) < 1e-12

    # ---- 7. bistability algebra -------------------------------------
    print("\n[7] Bistability condition")
    print("    rho = 0     -> |delta1|/kappa1 = "
          f"{bistability_delta_min(p, 1.0, 0.0)/p.kappa1:.4f}  (expect sqrt3 ="
          f" {np.sqrt(3.0):.4f})")
    for rho in (0.1, 0.25, 0.5, 0.577):
        U = -1.0
        kc = rho * abs(U)
        d = bistability_delta_min(p, U, kc)
        ok = bistability_ok(U, kc)
        print(f"    rho = {rho:5.3f}  bistable = {str(ok):5s}  "
              f"|delta1|/kappa1 = {d/p.kappa1 if np.isfinite(d) else np.inf:.3f}")
    print(f"    limiting rho = 1/sqrt(3) = {1/np.sqrt(3):.4f}"
          f"  -> |Domega| > sqrt(3)/2 * gamma2 = "
          f"{np.sqrt(3)/2:.3f} gamma2")

    # ---- 8. units -----------------------------------------------------
    print("\n[8] Dimensional spot checks")
    print(f"    V_eff (target ring)    = {p.V_eff:.3e} m^3 "
          f"= {p.V_eff*1e18:.1f} um^3")
    print(f"    FSR/2pi                = {p.FSR/TWOPI/1e9:.1f} GHz")
    print(f"    kappa1/2pi             = {p.kappa1/TWOPI/1e6:.2f} MHz")
    print(f"    kappa2/2pi             = {p.kappa2/TWOPI/1e6:.1f} MHz")
    print(f"    1 photon energy        = {HBAR*p.omega1:.3e} J")
    print(line)


if __name__ == "__main__":
    verify()


# ==================================================== tuning of Domega
def dDomega_dT(p: RingParams, dn_dT1: float = 3.3e-5,
               dn_dT2: float = 4.1e-5) -> float:
    """Thermal tuning of Domega = 2 w1 - w2 [rad/s per K].

    Each resonance tunes as dw_n/dT = -w_n (dn/dT)_n / ng_n, so
    d(Domega)/dT = 2 dw1/dT - dw2/dT.  The two terms very nearly cancel,
    which is why Domega is a slow, controllable knob even though each
    resonance moves fast.
    """
    dw1 = -p.omega1 * dn_dT1 / p.ng1
    dw2 = -p.omega2 * dn_dT2 / p.ng2
    return 2.0 * dw1 - dw2


def dDomega_dE(p: RingParams, r33_1: float = 30.8e-12,
               r33_2: float = 34.0e-12) -> float:
    """Electro-optic tuning of Domega [rad/s per (V/m)] via the Pockels effect.

    dn_n = -(1/2) n_n^3 r33(w_n) E, hence dw_n/dE = -w_n dn_n/(ng_n E).
    """
    dn1 = -0.5 * p.n1**3 * r33_1
    dn2 = -0.5 * p.n2**3 * r33_2
    dw1 = -p.omega1 * dn1 / p.ng1
    dw2 = -p.omega2 * dn2 / p.ng2
    return 2.0 * dw1 - dw2


def Q2_crossover(p: RingParams, x_gamma2: float = 2.0) -> float:
    """SH quality factor at which R = 1, at Domega = x_gamma2 * gamma2.

    R = 2 g^2 |Domega| / [(Domega^2 + kappa2^2) U_native] with
    Domega = x gamma2 = 2 x kappa2 and kappa2 = w2/(2 Q2):

        R = 4 x g^2 / [(4 x^2 + 1) kappa2 U_native]
          => Q2 = (4 x^2 + 1) w2 U_native / (8 x g^2).
    """
    return ((4.0 * x_gamma2**2 + 1.0) * p.omega2 * p.U_native
            / (8.0 * x_gamma2 * p.g_rate**2))


# ============================ exact (self-consistent) bistability threshold
def _fold_powers(p: RingParams, U: float, kc: float, delta1: float):
    """Input powers at the two folds of the CW cubic, at fixed delta1.

    P(n) = n[(k1 + kc n)^2 + (d1 - U n)^2] has turning points where
        a n^2 + b n + c = 0,   a = 3(U^2+kc^2), b = 4(k1 kc - d1 U),
                               c = d1^2 + k1^2.
    The smaller root n_- is a local MAXIMUM of P(n) and the larger root
    n_+ a local MINIMUM, so the bistable window is P in (P(n_+), P(n_-)).
    The ONSET of bistability as the power is raised is therefore P(n_+).
    Returns (P_onset, P_upper) or (inf, inf) if there is no fold.
    """
    a = 3.0 * (U**2 + kc**2)
    b = 4.0 * (p.kappa1 * kc - delta1 * U)
    c = delta1**2 + p.kappa1**2
    disc = b**2 - 4.0 * a * c
    if a <= 0 or disc <= 0 or b >= 0:
        return np.inf, np.inf
    rt = np.sqrt(disc)
    n_lo = (-b - rt) / (2.0 * a)          # local max of P(n)
    n_hi = (-b + rt) / (2.0 * a)          # local min of P(n)
    if n_lo <= 0:
        return np.inf, np.inf
    def P(n):
        return (n * HBAR * p.omega1 / (2.0 * p.kappa1_ex)
                * ((delta1 - U * n)**2 + (p.kappa1 + kc * n)**2))
    return P(n_hi), P(n_lo)


def bistability_threshold_exact(p: RingParams, Domega: float,
                                n_grid: int = 4001, span: float = 60.0):
    """Bistability onset power [W], minimised over pump detuning, with the
    cascaded coefficients evaluated at the PUMP-FRAME mismatch

        Delta_eff = Domega - 2 delta1,

    exactly as in the full CW solver.  This is the physically correct
    threshold for a pump scan at fixed cold-cavity mismatch Domega; the
    fixed-coefficient version `bistability_threshold_power` freezes U and
    kappa_casc at Domega and is only a local approximation.

    Returns (P_threshold, delta1_opt, Delta_eff_opt).
    """
    d1 = np.linspace(-span, span, n_grid) * p.kappa1
    best = (np.inf, np.nan, np.nan)
    for d in d1:
        Deff = Domega - 2.0 * d
        U = float(U_casc(p, Deff)) + p.U_native
        kc = float(kappa_casc(p, Deff))
        P_on, _ = _fold_powers(p, U, kc, d)
        if P_on < best[0]:
            best = (P_on, d, Deff)
    if not np.isfinite(best[0]):
        return best
    # golden-section style refinement around the grid minimum
    lo, hi = best[1] - 2.0 * span * p.kappa1 / n_grid, \
             best[1] + 2.0 * span * p.kappa1 / n_grid
    for _ in range(80):
        m1, m2 = lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0
        vals = []
        for d in (m1, m2):
            Deff = Domega - 2.0 * d
            U = float(U_casc(p, Deff)) + p.U_native
            kc = float(kappa_casc(p, Deff))
            vals.append(_fold_powers(p, U, kc, d)[0])
        if vals[0] < vals[1]:
            hi = m2
        else:
            lo = m1
    d = 0.5 * (lo + hi)
    Deff = Domega - 2.0 * d
    U = float(U_casc(p, Deff)) + p.U_native
    kc = float(kappa_casc(p, Deff))
    P_on, _ = _fold_powers(p, U, kc, d)
    return (P_on, d, Deff) if P_on < best[0] else best


def mi_threshold_selfconsistent(p: RingParams, Domega: float,
                                iters: int = 400):
    """Single-mode effective-Kerr MI estimate with self-consistent Delta_eff.

    Operating point: pump held on the SHIFTED resonance, delta1 = U n, so
    Delta_eff = Domega - 2 U n depends on n.  Iterates
        n = kappa1 / (|U(Delta_eff)| - kappa_casc(Delta_eff))
    under damping.  This is a HEURISTIC gain-equals-loss estimate for a
    single sideband pair treated as an ordinary Kerr medium; the multimode
    stability threshold is `ring_effects.mi_threshold_bogoliubov`.
    """
    n = 0.0
    U = float(U_casc(p, Domega)) + p.U_native
    kc = float(kappa_casc(p, Domega))
    for _ in range(iters):
        Deff = Domega - 2.0 * U * n
        U = float(U_casc(p, Deff)) + p.U_native
        kc = float(kappa_casc(p, Deff))
        if abs(U) <= kc:
            return np.inf, np.nan
        n_new = p.kappa1 / (abs(U) - kc)
        if abs(n_new - n) < 1e-12 * max(n_new, 1.0):
            n = n_new
            break
        n = 0.7 * n + 0.3 * n_new
    delta1 = U * n
    P = (n * HBAR * p.omega1 / (2.0 * p.kappa1_ex)
         * ((delta1 - U * n) ** 2 + (p.kappa1 + kc * n) ** 2))
    return P, n


def fsr_walkoff(p: RingParams) -> float:
    """Omega_2 - Omega_1 [rad/s], the FM/SH free-spectral-range mismatch."""
    return p.FSR * (p.ng1 / p.ng2 - 1.0)
