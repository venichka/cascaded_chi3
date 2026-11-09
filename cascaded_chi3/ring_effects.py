"""Observables of the cascaded chi2:chi2 Kerr response in a LiNbO3 microring.

Companion to `ring_params.py`.  Produces the figure set in
`figures/ring/` and prints a numerical summary table.

Run:  conda run -n ml_base python ring_effects.py
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from style import apply_style, COLORS
import ring_params as rp
from ring_params import (
    RingParams, HBAR, TWOPI, design_points,
    U_casc, kappa_casc, loss_to_kerr, R_ratio, n2_casc_eff, chi3_casc_eff,
    eta_norm_from_g, n_from_power, power_from_n,
    bistability_ok, bistability_delta_min, bistability_threshold_power,
    mi_threshold_power, switching_photons, saturation_photons,
    anharmonicity_fom, resonant_fom, ThermalModel,
    bistability_threshold_exact, mi_threshold_selfconsistent, fsr_walkoff,
)

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "figures", "ring")
os.makedirs(FIGDIR, exist_ok=True)
apply_style()


def _save(fig, name):
    path = os.path.join(FIGDIR, name)
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=160)
    plt.close(fig)
    print(f"    wrote {os.path.relpath(path)}.pdf/.png")


# =====================================================================
#  Exact CW steady state (no adiabatic approximation needed: eliminating
#  b in steady state is exact).
# =====================================================================
def cw_photon_numbers(p: RingParams, P_in, delta1, Domega, include_kerr=True,
                      cascade=True):
    """Real positive roots n = |a|^2 of the CW steady state.

    |p|^2 = n [ (kappa1 + kc n)^2 + (delta1 - U n)^2 ],  |p|^2 = 2 k1x P/(hbar w1)

    with the *pump-frame* SH detuning Delta_eff = Domega - 2 delta1, so that
    U = U_casc(Delta_eff) (+ U_native if include_kerr) and kc = kappa_casc.
    Returns a sorted array of roots.
    """
    Deff = Domega - 2.0 * delta1
    U = float(U_casc(p, Deff)) if cascade else 0.0
    kc = float(kappa_casc(p, Deff)) if cascade else 0.0
    if include_kerr:
        U += p.U_native
    rhs = 2.0 * p.kappa1_ex * P_in / (HBAR * p.omega1)
    # n[(k1+kc n)^2 + (d1 - U n)^2] - rhs = 0  ->  cubic in n
    c3 = kc**2 + U**2
    c2 = 2.0 * (p.kappa1 * kc - delta1 * U)
    c1 = p.kappa1**2 + delta1**2
    roots = np.roots([c3, c2, c1, -rhs]) if c3 != 0 else np.roots([c2, c1, -rhs])
    real = roots[np.abs(roots.imag) < 1e-9 * (1.0 + np.abs(roots))].real
    return np.sort(real[real > -1e-12].clip(0.0))


def transmission(p: RingParams, n, P_in, delta1, Domega, include_kerr=True,
                 cascade=True):
    """Bus transmission |s_out/s_in|^2 for a given steady-state n."""
    Deff = Domega - 2.0 * delta1
    U = float(U_casc(p, Deff)) if cascade else 0.0
    kc = float(kappa_casc(p, Deff)) if cascade else 0.0
    if include_kerr:
        U += p.U_native
    s_in = np.sqrt(P_in / (HBAR * p.omega1))
    denom = (p.kappa1 + kc * n) + 1j * (delta1 - U * n)
    a = np.sqrt(2.0 * p.kappa1_ex) * s_in / denom
    s_out = -s_in + np.sqrt(2.0 * p.kappa1_ex) * a
    return np.abs(s_out / s_in) ** 2


def n_selfconsistent(p: RingParams, P_in, delta1, Domega, cascade=True):
    """Steady-state photon number including the cascaded nonlinear loss.

    At delta1 = 0 the cubic has a single positive root; off resonance the
    lowest branch is returned.  Using the *linear* build-up formula instead
    badly overestimates n whenever kappa_casc * n is comparable to kappa1,
    which happens already at tens of microwatts at the reference point.
    """
    r = cw_photon_numbers(p, P_in, delta1, Domega, cascade=cascade)
    return float(r.min()) if r.size else 0.0


# =====================================================================
#  Time-domain TCMT in photon units (validation of the adiabatic result)
# =====================================================================
def _rhs(t, y, d1, d2, k2, G, P, Un):
    """Dimensionless TCMT: time in 1/kappa1, amplitudes in sqrt(n_scale).

    Un is the native Kerr coefficient U_native * n_scale / kappa1; it enters
    the FM equation as +i Un |a|^2 a, i.e. a red shift of the resonance.
    """
    x, yy, u, v = y
    nn = x * x + yy * yy
    Fx = -x + d1 * yy + 2 * G * (x * v - yy * u) - Un * nn * yy + P
    Fy = -d1 * x - yy - 2 * G * (x * u + yy * v) + Un * nn * x
    Fu = -k2 * u + d2 * v + 2 * G * x * yy
    Fv = -d2 * u - k2 * v - G * (x * x - yy * yy)
    return [Fx, Fy, Fu, Fv]


def tcmt_steady(p: RingParams, P_in, delta1, Domega, t_end_factor=400.0):
    """Integrate the full two-mode TCMT to steady state; return (n, nb, T).

    Non-dimensionalised (tau = kappa1 t, amplitudes in units of
    A = sqrt(n_lin) with n_lin the on-resonance linear photon number) so
    that every quantity fed to the integrator is O(1).
    """
    delta2 = -(Domega - 2.0 * delta1)
    s_in = np.sqrt(P_in / (HBAR * p.omega1))
    pdrive = np.sqrt(2.0 * p.kappa1_ex) * s_in
    A = np.sqrt(n_from_power(p, P_in))            # amplitude scale
    G = p.g_rate * A / p.kappa1
    P = pdrive / (p.kappa1 * A)
    sol = solve_ivp(_rhs, (0.0, t_end_factor), [0.0, 0.0, 0.0, 0.0],
                    args=(delta1 / p.kappa1, delta2 / p.kappa1,
                          p.kappa2 / p.kappa1, G, P,
                          p.U_native * A**2 / p.kappa1),
                    method="DOP853", rtol=1e-12, atol=1e-13)
    x, yy, u, v = sol.y[:, -1] * A
    a = x + 1j * yy
    n = abs(a) ** 2
    nb = u * u + v * v
    s_out = -s_in + np.sqrt(2.0 * p.kappa1_ex) * a
    return n, nb, abs(s_out / s_in) ** 2


# =====================================================================
#  Linear stability of the CW fixed point (independent implementation)
# =====================================================================
def jacobian(p: RingParams, n, delta1, Domega):
    """4x4 real Jacobian at the fixed point, gauge a = sqrt(n) (real > 0)."""
    g = p.g_rate
    delta2 = -(Domega - 2.0 * delta1)
    x, yv = np.sqrt(n), 0.0
    D = delta2**2 + p.kappa2**2
    # b = -i g a^2/(i d2 + k2) with a^2 = n real
    b = -1j * g * n / (1j * delta2 + p.kappa2)
    u, v = b.real, b.imag
    J = np.array([
        [-p.kappa1 + 2 * g * v, delta1 - 2 * g * u, -2 * g * yv, 2 * g * x],
        [-delta1 - 2 * g * u, -p.kappa1 - 2 * g * v, -2 * g * x, -2 * g * yv],
        [2 * g * yv, 2 * g * x, -p.kappa2, delta2],
        [-2 * g * x, 2 * g * yv, -delta2, -p.kappa2],
    ])
    return J


def stability_ceiling(p: RingParams, Domega, delta1=0.0,
                      n_lo=1e2, n_hi=1e12, tol=1e-4):
    """Smallest intracavity photon number at which the CW fixed point loses
    stability, found by bisection on max Re(eig J).  Returns (n_c, shift/k1).
    """
    def maxre(n):
        return np.max(np.linalg.eigvals(jacobian(p, n, delta1, Domega)).real)

    if maxre(n_hi) < 0:
        return np.inf, np.inf
    lo, hi = n_lo, n_hi
    if maxre(lo) > 0:
        return lo, 0.0
    while hi / lo > 1.0 + tol:
        mid = np.sqrt(lo * hi)
        if maxre(mid) > 0:
            hi = mid
        else:
            lo = mid
    n_c = np.sqrt(lo * hi)
    shift = abs(float(U_casc(p, Domega))) * n_c
    return n_c, shift / p.kappa1


# =====================================================================
#  Multimode (Bogoliubov) parametric-instability threshold
# =====================================================================
def bogoliubov_matrix(p: RingParams, n, delta1, Domega, q, D2):
    """4x4 complex stability matrix for the sideband pair +-q.

    Modes: pump a0 (m_p), signal a_+ (m_p+q), idler a_- (m_p-q), SH b0
    (2 m_p) driven by the pump, and SH sidebands b_+- (2 m_p +- q).
    With equal three-wave couplings g throughout,

        H/hbar = g(a0^dag^2 b0 + h.c.) + 2g(a_+^dag a_-^dag b0 + h.c.)
                 + 2g(a0^dag a_+^dag b_+ + h.c.)
                 + 2g(a0^dag a_-^dag b_- + h.c.)

    Linearising about (a0, b0) with a0 = sqrt(n) real (gauge) gives a closed
    system for v = (da_+, da_-^*, db_+, db_-^*).  The parametric gain enters
    through b0; the pump-induced cross-phase and cascaded loss on the
    sidebands enter through b_+-, whose detunings walk off with q because
    the FM and SH free spectral ranges differ.

    Assumes equal modal overlaps for every triplet (see the note) and
    neglects the native Kerr on the sidebands, which is ~10^-3 of the
    cascaded term here.
    """
    g = p.g_rate
    a0 = np.sqrt(n)
    d2 = -(Domega - 2.0 * delta1)                 # pump-frame SH mismatch
    b0 = -1j * g * n / (1j * d2 + p.kappa2)
    dW = rp.fsr_walkoff(p)                        # Omega_2 - Omega_1
    dint = 0.5 * D2 * q**2                        # integrated FM dispersion
    dp = delta1 + dint                            # signal and idler detuning
    dbp = d2 + q * dW                             # SH mode 2m_p + q
    dbm = d2 - q * dW                             # SH mode 2m_p - q
    k1, k2 = p.kappa1, p.kappa2
    M = np.array([
        [-(1j * dp + k1), -2j * g * b0,      -2j * g * a0,        0.0],
        [2j * g * np.conj(b0), (1j * dp - k1), 0.0,               2j * g * a0],
        [-2j * g * a0,        0.0,            -(1j * dbp + k2),   0.0],
        [0.0,                 2j * g * a0,    0.0,                (1j * dbm - k2)],
    ], dtype=complex)
    return M


def mi_gain(p: RingParams, n, delta1, Domega, q, D2):
    """Largest real part of the Bogoliubov eigenvalues for sideband index q."""
    return float(np.max(np.linalg.eigvals(
        bogoliubov_matrix(p, n, delta1, Domega, q, D2)).real))


def mi_threshold_bogoliubov(p: RingParams, Domega, D2,
                            q_max=40, P_lo=1e-9, P_hi=1e-1,
                            d1_grid=61, d1_span=12.0, tol=1e-3):
    """Multimode parametric threshold: lowest input power at which ANY
    sideband pair goes unstable, minimised over pump detuning.

    Returns (P_threshold, delta1_opt, q_opt).
    """
    qs = np.arange(1, q_max + 1)
    d1s = np.linspace(-d1_span, d1_span, d1_grid) * p.kappa1

    def unstable(P):
        best = (-np.inf, np.nan, 0)
        for d1 in d1s:
            n = n_selfconsistent(p, P, d1, Domega)
            if n <= 0:
                continue
            for q in qs:
                gmax = mi_gain(p, n, d1, Domega, q, D2)
                if gmax > best[0]:
                    best = (gmax, d1, q)
        return best

    if unstable(P_hi)[0] <= 0:
        return np.inf, np.nan, 0
    lo, hi = P_lo, P_hi
    if unstable(lo)[0] > 0:
        return lo, np.nan, 0
    while hi / lo > 1.0 + tol:
        mid = np.sqrt(lo * hi)
        if unstable(mid)[0] > 0:
            hi = mid
        else:
            lo = mid
    g, d1, q = unstable(hi)
    return np.sqrt(lo * hi), d1, q


# =====================================================================
#  Small-signal modulation transfer function (proper linear response)
# =====================================================================
def _jacobian_full(p: RingParams, x, y, u, v, delta1, delta2, Un):
    """Real 4x4 Jacobian of the driven two-mode system incl. native Kerr."""
    g = p.g_rate
    nn = x * x + y * y
    J = np.empty((4, 4))
    J[0] = [-p.kappa1 + 2 * g * v - Un * 2 * x * y,
            delta1 - 2 * g * u - Un * (nn + 2 * y * y),
            -2 * g * y, 2 * g * x]
    J[1] = [-delta1 - 2 * g * u + Un * (nn + 2 * x * x),
            -p.kappa1 - 2 * g * v + Un * 2 * x * y,
            -2 * g * x, -2 * g * y]
    J[2] = [2 * g * y, 2 * g * x, -p.kappa2, delta2]
    J[3] = [-2 * g * x, 2 * g * y, -delta2, -p.kappa2]
    return J


def modulation_transfer(p: RingParams, P_in, delta1, Domega, Omega,
                        cascade=True):
    """|dn(Omega)/dn(0)| for a small modulation of the input power.

    cascade=True  : full driven two-mode system, SH dynamics retained.
    cascade=False : MATCHED reference -- a single-mode cavity carrying the
                    same photon number, the same total per-photon Kerr
                    U_native + U_casc and the same nonlinear loss
                    kappa_casc, but responding instantaneously.  Any
                    difference between the two is therefore due solely to
                    the finite response time of the SH mode, which is the
                    quantity the bandwidth claim is about.

    The drive perturbation points along p (not along the real axis): with
    a gauged real, p = sqrt(n)[(k1 + kc n) + i(d1 - U n)].
    """
    Deff = Domega - 2.0 * delta1
    Uc = float(U_casc(p, Deff))
    kc = float(kappa_casc(p, Deff))
    Utot = p.U_native + Uc
    n = n_selfconsistent(p, P_in, delta1, Domega, cascade=True)
    if n <= 0:
        return np.zeros_like(np.atleast_1d(Omega))
    x = np.sqrt(n)
    delta2 = -Deff
    p_re, p_im = x * (p.kappa1 + kc * n), x * (delta1 - Utot * n)
    nrm = np.hypot(p_re, p_im)

    if cascade:
        b = -1j * p.g_rate * n / (1j * delta2 + p.kappa2)
        J = _jacobian_full(p, x, 0.0, b.real, b.imag, delta1, delta2,
                           p.U_native)
        drive = np.array([p_re, p_im, 0.0, 0.0]) / nrm
        dim = 4
    else:
        # instantaneous effective Kerr + nonlinear loss, same n, U, kc
        J = np.array([
            [-p.kappa1 - 3.0 * kc * n, delta1 - Utot * n],
            [-delta1 + 3.0 * Utot * n, -p.kappa1 - kc * n],
        ])
        drive = np.array([p_re, p_im]) / nrm
        dim = 2

    Om = np.atleast_1d(np.asarray(Omega, dtype=float))
    out = np.empty(Om.shape)
    I = np.eye(dim)
    for i, w in enumerate(Om):
        dy = np.linalg.solve(-1j * w * I - J, drive)
        out[i] = abs(2.0 * x * dy[0])
    dy0 = np.linalg.solve(-J, drive)
    return out / abs(2.0 * x * dy0[0])


def jacobian_poles(p: RingParams, P_in, delta1, Domega):
    """Eigenvalues of the full 4x4 Jacobian: two FM-like, two SH-like."""
    Deff = Domega - 2.0 * delta1
    n = n_selfconsistent(p, P_in, delta1, Domega)
    x = np.sqrt(n)
    delta2 = -Deff
    b = -1j * p.g_rate * n / (1j * delta2 + p.kappa2)
    J = _jacobian_full(p, x, 0.0, b.real, b.imag, delta1, delta2, p.U_native)
    return np.linalg.eigvals(J)


# =====================================================================
#  FIGURES
# =====================================================================
def fig_R_landscape(dps):
    ppln, modal, lu = dps["ppln"], dps["modal"], dps["lu"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # -- (a) R vs Q2 at the reference operating point Domega = 2 gamma2 --
    Q2 = np.logspace(4, 6.1, 400)
    for q, col, ls in ((ppln, COLORS["tcmt"], "-"),
                       (modal, COLORS["second"], "--")):
        Rv = []
        for Q in Q2:
            qq = q.replace(Q2_load=Q, Q2_int=2.0 * Q)
            Rv.append(float(R_ratio(qq, 2.0 * qq.gamma2)))
        ax[0].loglog(Q2, Rv, ls, color=col, label=q.label)
    # demonstrated device marker
    ax[0].plot([lu.Q2_load], [float(R_ratio(lu, 2.0 * lu.gamma2))], "*",
               ms=13, color=COLORS["highlight"], mec="k", mew=0.6, zorder=5,
               label="PPLN demonstrated (Lu 2020)")
    ax[0].axhline(1.0, color="k", lw=0.9, ls=":")
    # metasurface references after the 2026-10-06 prefactor correction (x4)
    for y, lab in ((49.0, r"metasurface note, $Q_2\!=\!4\times10^3$"),
                   (10.0, "best SH-only sample (S4)")):
        ax[0].axhline(y, color=COLORS["linear"], lw=0.9, ls="-.")
        ax[0].text(1.05e4, y * 1.25, lab, fontsize=8.5,
                   color="0.35")
    ax[0].axvspan(3e5, 6e5, color=COLORS["highlight"], alpha=0.10)
    ax[0].text(4.0e5, 4.0e3, "measured\n$Q_2$ range\nat 780 nm",
               fontsize=8.5, ha="center", color="0.35")
    ax[0].set_xlabel(r"SH quality factor $Q_2$")
    ax[0].set_ylabel(r"$R=|\alpha_3^{\rm casc}|/\alpha_3^{\rm native}$")
    ax[0].set_title(r"(a) Cascade enhancement, $\Delta\omega=2\gamma_2$")
    ax[0].legend(fontsize=8.6, loc="upper left")
    ax[0].set_ylim(1e0, 2e4)

    # -- (b) R and loss ratio vs operating point --
    x = np.linspace(0.3, 8.0, 500)
    D = x * ppln.gamma2
    ax[1].plot(x, R_ratio(ppln, D), color=COLORS["tcmt"], label=r"$R$")
    ax[1].set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax[1].set_ylabel(r"$R$", color=COLORS["tcmt"])
    ax[1].tick_params(axis="y", labelcolor=COLORS["tcmt"])
    ax2 = ax[1].twinx()
    ax2.plot(x, loss_to_kerr(ppln, D), color=COLORS["wrong"],
             label=r"$\rho=\kappa_{\rm casc}/|U_{\rm casc}|$")
    ax2.axhline(1 / np.sqrt(3), color=COLORS["wrong"], ls=":", lw=1.2)
    ax2.text(5.2, 1 / np.sqrt(3) + 0.03,
             r"$\rho=1/\sqrt{3}$: local bistability/MI limit",
             fontsize=8.5, color=COLORS["wrong"], ha="right")
    ax2.set_ylabel(r"loss-to-Kerr ratio $\rho$", color=COLORS["wrong"])
    ax2.tick_params(axis="y", labelcolor=COLORS["wrong"])
    ax2.set_ylim(0, 1.0)
    ax2.grid(False)
    ax[1].axvline(2.0, color=COLORS["highlight"], lw=1.2, ls="--")
    ax[1].text(2.12, R_ratio(ppln, 2 * ppln.gamma2) * 0.55,
               "reference\noperating point", fontsize=8.5,
               color=COLORS["highlight"])
    ax[1].set_title(r"(b) Enhancement vs cascaded loss")
    fig.tight_layout()
    _save(fig, "ring_R_landscape")


def fig_sign_flip(dps):
    p = dps["ppln"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # -- (a) per-photon shift, cascaded + native --
    x = np.linspace(-8, 8, 1200)
    D = x * p.gamma2
    Uc = U_casc(p, D) / TWOPI
    Un = p.U_native / TWOPI
    ax[0].plot(x, Uc, color=COLORS["tcmt"], label=r"cascaded $U_{\rm casc}$")
    ax[0].axhline(Un, color=COLORS["second"], ls="--",
                  label=r"native $U_{\rm native}$ ($\times 1$)")
    ax[0].plot(x, Uc + Un, color="k", lw=1.1, ls=":",
               label="total")
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].axvline(0, color="k", lw=0.8)
    ax[0].set_yscale("symlog", linthresh=1.0)
    ax[0].set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax[0].set_ylabel(r"per-photon shift $U/2\pi$  (Hz)")
    ax[0].set_title("(a) Sign-tunable Kerr (self-focusing $>0$)")
    ax[0].legend(fontsize=8.6, loc="upper right")
    ax[0].annotate("self-defocusing", xy=(3, -600), fontsize=9,
                   color="0.3")
    ax[0].annotate("self-focusing", xy=(-6.5, 600), fontsize=9, color="0.3")

    # -- (b) measurable: FM resonance shift vs pump power, both signs --
    P = np.logspace(-8, -3.4, 300)
    for xs, col, lab in ((+2.0, COLORS["tcmt"], r"$\Delta\omega=+2\gamma_2$"),
                         (-2.0, COLORS["highlight"],
                          r"$\Delta\omega=-2\gamma_2$")):
        D0 = xs * p.gamma2
        U = float(U_casc(p, D0))
        n = np.array([n_selfconsistent(p, PP, 0.0, D0) for PP in P])
        # resonance shift is  delta_omega_1 = -U n
        ax[1].plot(P * 1e6, -U * n / p.kappa1, color=col, label=lab)
    n_lin = np.array([n_from_power(p, PP) for PP in P])
    ax[1].plot(P * 1e6, -p.U_native * n_lin / p.kappa1, color=COLORS["second"],
               ls="--", label="native Kerr only")
    th = ThermalModel()
    ax[1].plot(P * 1e6,
               -np.array([thermal_shift_full(p, th, PP, 2.0 * p.gamma2)
                          for PP in P]) / p.kappa1,
               color=COLORS["wrong"], ls="-.", label="thermal (always red)")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xscale("log")
    ax[1].set_yscale("symlog", linthresh=1e-2)
    ax[1].set_xlabel(r"input power $P_{\rm in}$  ($\mu$W)")
    ax[1].set_ylabel(r"resonance shift  $\delta\omega_1/\kappa_1$")
    ax[1].set_title("(b) What the experiment sees")
    ax[1].legend(fontsize=8.4, loc="lower left")
    fig.tight_layout()
    _save(fig, "ring_sign_flip")


def threshold_vs_detuning(p: RingParams, xs, D2=10e6 * TWOPI, mi=False):
    """Exact bistability (and optionally multimode MI) threshold vs Domega."""
    Pb = np.array([bistability_threshold_exact(p, x * p.gamma2)[0] for x in xs])
    if not mi:
        return Pb, None
    Pm = np.array([mi_threshold_bogoliubov(p, x * p.gamma2, D2, q_max=40,
                                           d1_grid=41)[0] for x in xs])
    return Pb, Pm


def fig_thresholds(dps):
    ppln = dps["ppln"]
    D2 = 10e6 * TWOPI
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.1))
    Q2 = np.logspace(4, 6.1, 22)

    Pb_c, Pm_c, Psw_lin, Psw_sc = [], [], [], []
    for Q in Q2:
        q = ppln.replace(Q2_load=Q, Q2_int=2.0 * Q)
        D = 2.0 * q.gamma2
        U = float(U_casc(q, D)) + q.U_native
        Pb_c.append(bistability_threshold_exact(q, D)[0])
        Pm_c.append(mi_threshold_bogoliubov(q, D, D2, q_max=40, d1_grid=41)[0])
        n_sw = switching_photons(q, U)
        Psw_lin.append(power_from_n(q, n_sw))
        Psw_sc.append(_drive_for_photons(q, n_sw, D))
    Pb_n = bistability_threshold_power(ppln, ppln.U_native, 0.0)
    Pm_n = (ppln.kappa1 / ppln.U_native) * HBAR * ppln.omega1 * ppln.kappa1

    ax[0].loglog(Q2, np.array(Pb_c) * 1e6, "o-", ms=3.5, color=COLORS["tcmt"],
                 label="cascaded: bistability (exact)")
    ax[0].loglog(Q2, np.array(Pm_c) * 1e6, "s--", ms=3.5,
                 color=COLORS["second"], label="cascaded: parametric (multimode)")
    ax[0].axhline(Pb_n * 1e6, color=COLORS["tcmt"], ls=":", lw=1.4,
                  label="native: bistability")
    ax[0].axhline(Pm_n * 1e6, color=COLORS["second"], ls="-.", lw=1.4,
                  label="native: parametric")
    ax[0].axhspan(10.0, 1e4, color=COLORS["wrong"], alpha=0.07)
    ax[0].text(1.1e4, 1.6e3, "photorefractive onset\nreported above $\\sim\\!10\\,\\mu$W",
               fontsize=8.0, color=COLORS["wrong"])
    ax[0].set_xlabel(r"SH quality factor $Q_2$")
    ax[0].set_ylabel(r"threshold input power  ($\mu$W)")
    ax[0].set_title(r"(a) Threshold powers, $\Delta\omega=2\gamma_2$")
    ax[0].legend(fontsize=8.0, loc="lower left")

    # -- (b) operating point: exact vs frozen-coefficient --
    xs = np.linspace(0.05, 12.0, 220)
    Pb, _ = threshold_vs_detuning(ppln, xs)
    Pfroz = []
    for x in xs:
        D = x * ppln.gamma2
        Pfroz.append(bistability_threshold_power(
            ppln, float(U_casc(ppln, D)) + ppln.U_native,
            float(kappa_casc(ppln, D))))
    ax[1].semilogy(xs, Pb * 1e6, color=COLORS["tcmt"],
                   label=r"exact: $U,\kappa_{\rm casc}$ at $\Delta_{\rm eff}$")
    ax[1].semilogy(xs, np.array(Pfroz) * 1e6, color=COLORS["adiabatic"],
                   ls=":", label=r"frozen at cold $\Delta\omega$")
    ib = int(np.argmin(Pb))
    ax[1].plot([xs[ib]], [Pb[ib] * 1e6], "*", ms=13,
               color=COLORS["highlight"], mec="k", mew=0.6, zorder=5)
    ax[1].annotate(rf"optimum $\Delta\omega={xs[ib]:.1f}\gamma_2$"
                   "\n" rf"$P_{{\rm th}}={Pb[ib]*1e6:.0f}\,\mu$W",
                   xy=(xs[ib], Pb[ib] * 1e6), xytext=(xs[ib] + 1.1,
                                                      Pb[ib] * 1e6 * 2.6),
                   fontsize=8.6, color=COLORS["highlight"])
    ax[1].axvline(np.sqrt(3) / 2, color=COLORS["wrong"], ls="-.", lw=1.1)
    ax[1].text(np.sqrt(3) / 2 + 0.2, Pb[ib] * 1e6 * 22,
               r"$\frac{\sqrt{3}}{2}\gamma_2$: where the" "\n"
               r"frozen model gives up", fontsize=8.0, color=COLORS["wrong"])
    ax[1].set_xlabel(r"cold-cavity $\Delta\omega/\gamma_2$")
    ax[1].set_ylabel(r"bistability threshold  ($\mu$W)")
    ax[1].set_title("(b) Choosing the operating point")
    ax[1].legend(fontsize=8.4, loc="upper right")

    # -- (c) switching drive --
    ax[2].loglog(Q2, np.array(Psw_sc) * 1e6, "o-", ms=3.5,
                 color=COLORS["tcmt"], label=r"cascaded, self-consistent")
    ax[2].loglog(Q2, np.array(Psw_lin) * 1e6, ls="--", color=COLORS["tcmt"],
                 lw=1.1, alpha=0.7, label="cascaded, linear-cavity equivalent")
    ax[2].axhline(power_from_n(ppln, switching_photons(ppln, ppln.U_native))
                  * 1e6, color=COLORS["second"], ls="-.", label="native")
    ax[2].set_xlabel(r"SH quality factor $Q_2$")
    ax[2].set_ylabel(r"drive for a one-linewidth shift  ($\mu$W)")
    ax[2].set_title("(c) All-optical switching drive")
    ax[2].legend(fontsize=8.0, loc="upper right")
    fig.tight_layout()
    _save(fig, "ring_thresholds")
    return xs, Pb


def _drive_for_photons(p: RingParams, n_target, Domega, d1x=0.0):
    """Input power needed to reach n_target at fixed delta1 (self-consistent)."""
    lo, hi = 1e-9, 1.0
    if n_selfconsistent(p, hi, d1x * p.kappa1, Domega) < n_target:
        return np.inf
    for _ in range(120):
        mid = np.sqrt(lo * hi)
        if n_selfconsistent(p, mid, d1x * p.kappa1, Domega) > n_target:
            hi = mid
        else:
            lo = mid
    return np.sqrt(lo * hi)


def fig_bistability(dps):
    p = dps["ppln"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # -- (a) hysteresis: n vs P for cascade (defocusing) and native --
    #  The +Domega and -Domega branches give identical |n| curves (only the
    #  sign of the shift, and hence of the required pump detuning, flips), so
    #  only one is drawn; the optimal Domega is shown alongside.
    cases = [
        (+2.0, True, COLORS["tcmt"], r"cascade $\Delta\omega=+2\gamma_2$"),
        (+3.3, True, COLORS["highlight"],
         r"cascade $\Delta\omega=+3.3\gamma_2$ (optimum)"),
    ]
    P = np.logspace(-8.0, -1.5, 1100)
    for xs, casc, col, lab in cases:
        D0 = xs * p.gamma2
        U0 = float(U_casc(p, D0)) + p.U_native
        kc0 = float(kappa_casc(p, D0))
        d1 = np.sign(U0) * 1.45 * bistability_delta_min(p, U0, kc0)
        lo, hi = [], []
        for PP in P:
            r = cw_photon_numbers(p, PP, d1, D0, cascade=casc)
            lo.append(r.min())
            hi.append(r.max())
        ax[0].loglog(P * 1e6, hi, color=col, label=lab)
        ax[0].loglog(P * 1e6, lo, color=col, lw=1.0, alpha=0.55)
    # native reference at its own optimal detuning
    Un = p.U_native
    d1n = np.sign(Un) * 1.45 * np.sqrt(3.0) * p.kappa1
    lo, hi = [], []
    for PP in P:
        r = cw_photon_numbers(p, PP, d1n, 0.0, cascade=False)
        lo.append(r.min())
        hi.append(r.max())
    ax[0].loglog(P * 1e6, hi, color=COLORS["second"], ls="--",
                 label="native Kerr only")
    ax[0].loglog(P * 1e6, lo, color=COLORS["second"], ls="--", lw=1.0,
                 alpha=0.55)
    ax[0].set_xlabel(r"input power  ($\mu$W)")
    ax[0].set_ylabel(r"intracavity photons $n$")
    ax[0].set_title("(a) Bistability turns on 2.5 decades earlier")
    ax[0].legend(fontsize=8.4, loc="upper left")

    # -- (b) transmission lineshapes at fixed power, both signs --
    P0 = 30e-6
    d1 = np.linspace(-16, 16, 1600) * p.kappa1
    for xs, col, lab in ((+2.0, COLORS["tcmt"], r"$\Delta\omega=+2\gamma_2$"),
                         (-2.0, COLORS["highlight"],
                          r"$\Delta\omega=-2\gamma_2$"),
                         (None, COLORS["second"], "native only")):
        Tup, Tdn = [], []
        for d in d1:
            if xs is None:
                r = cw_photon_numbers(p, P0, d, 0.0, cascade=False)
                Tv = [transmission(p, n, P0, d, 0.0, cascade=False) for n in r]
            else:
                D0 = xs * p.gamma2
                r = cw_photon_numbers(p, P0, d, D0)
                Tv = [transmission(p, n, P0, d, D0) for n in r]
            Tup.append(max(Tv))
            Tdn.append(min(Tv))
        ls = "--" if xs is None else "-"
        ax[1].plot(d1 / p.kappa1, Tdn, ls, color=col, label=lab)
        ax[1].plot(d1 / p.kappa1, Tup, ls, color=col, lw=0.9, alpha=0.5)
    ax[1].set_xlabel(r"pump detuning $\delta_1/\kappa_1$")
    ax[1].set_ylabel("transmission")
    ax[1].set_title(rf"(b) Lineshape tilt reverses, $P_{{\rm in}}={P0*1e6:.0f}\,\mu$W")
    ax[1].legend(fontsize=8.6, loc="lower left")
    fig.tight_layout()
    _save(fig, "ring_bistability")


def fig_stability(dps):
    p = dps["ppln"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    xs = np.linspace(0.6, 8.0, 26)
    n_c, sh, an = [], [], []
    for x in xs:
        D = x * p.gamma2
        nc, s = stability_ceiling(p, D)
        n_c.append(nc)
        sh.append(s)
        an.append(saturation_photons(p, D))
    n_c = np.array(n_c)
    an = np.array(an)

    ax[0].semilogy(xs, n_c, "o-", color=COLORS["tcmt"], ms=4,
                   label=r"numerical: max Re$\,\lambda(J)>0$")
    ax[0].semilogy(xs, an, "--", color=COLORS["analytic"],
                   label=r"analytic $n_{\rm sat}=\Delta\omega^2/2g^2$")
    ax[0].set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax[0].set_ylabel(r"critical photon number $n_c$")
    ax[0].set_title("(a) Dynamic range: where CW operation ends")
    ax[0].legend(fontsize=8.8)

    ax[1].plot(xs, np.array(sh), "o-", color=COLORS["tcmt"], ms=4,
               label=r"max $|\delta\omega_1|/\kappa_1$ at the ceiling")
    ax[1].axhline(np.sqrt(3), color=COLORS["wrong"], ls="--",
                  label=r"bistability needs $\sqrt{3}\,\kappa_1$")
    ax[1].axhline(1.0, color=COLORS["second"], ls=":",
                  label=r"switching needs $1\,\kappa_1$")
    ax[1].set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax[1].set_ylabel(r"available shift  $|\delta\omega_1|/\kappa_1$")
    ax[1].set_title("(b) Headroom above the useful thresholds")
    ax[1].set_yscale("log")
    ax[1].legend(fontsize=8.6, loc="lower right")
    fig.tight_layout()
    _save(fig, "ring_stability")
    return xs, n_c, an, sh


def fig_validation(dps):
    p = dps["ppln"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # -- (a) full TCMT vs exact CW cubic, transmission lineshape --
    P0 = 30e-6
    D0 = 2.0 * p.gamma2
    d1 = np.linspace(-10, 4, 60) * p.kappa1
    T_num, T_an, nb_num, nb_an = [], [], [], []
    for d in d1:
        n, nb, T = tcmt_steady(p, P0, d, D0)
        T_num.append(T)
        nb_num.append(nb)
        r = cw_photon_numbers(p, P0, d, D0)
        # the branch the upward integration from vacuum lands on: lowest
        T_an.append(transmission(p, r.min(), P0, d, D0))
        Deff = D0 - 2 * d
        nb_an.append(p.g_rate**2 * r.min()**2 / (Deff**2 + p.kappa2**2))
    ax[0].plot(d1 / p.kappa1, T_an, "-", color=COLORS["analytic"],
               label="closed form (exact CW)")
    ax[0].plot(d1 / p.kappa1, T_num, marker="o", ms=4, mfc="white", mew=1.1,
               ls="none", color=COLORS["tcmt"], label="TCMT integration")
    ax[0].set_xlabel(r"pump detuning $\delta_1/\kappa_1$")
    ax[0].set_ylabel("transmission")
    ax[0].set_title("(a) Two-mode TCMT vs closed form")
    ax[0].legend(fontsize=9)

    # -- (b) true small-signal modulation transfer function --
    P0 = 20e-6
    Om = np.logspace(-1.3, 2.3, 320) * p.kappa1
    Hc = modulation_transfer(p, P0, 0.0, D0, Om, cascade=True)
    Hn = modulation_transfer(p, P0, 0.0, D0, Om, cascade=False)
    ax[1].loglog(Om / TWOPI / 1e6, Hc, color=COLORS["tcmt"],
                 label="cascade: SH dynamics retained")
    ax[1].loglog(Om / TWOPI / 1e6, Hn, color=COLORS["second"], ls="--",
                 label="matched instantaneous Kerr")
    ax[1].axhline(1 / np.sqrt(2), color="k", lw=0.8, ls=":")
    i_c = int(np.argmin(np.abs(Hc - 1 / np.sqrt(2))))
    i_n = int(np.argmin(np.abs(Hn - 1 / np.sqrt(2))))
    lam = jacobian_poles(p, P0, 0.0, D0)
    lam = sorted(lam, key=lambda z: abs(z))
    ax[1].text(0.97, 0.62,
               rf"$f_{{3\rm dB}}$: {Om[i_c]/TWOPI/1e6:.0f} vs "
               rf"{Om[i_n]/TWOPI/1e6:.0f} MHz" "\n"
               rf"the SH costs {100*(1-Om[i_c]/Om[i_n]):.0f}% of bandwidth"
               "\n\n"
               rf"poles: FM $|\lambda|={abs(lam[0])/p.kappa1:.1f}\kappa_1$,"
               "\n" rf"SH $|\lambda|={abs(lam[-1])/p.kappa1:.0f}\kappa_1$",
               transform=ax[1].transAxes, ha="right", va="top", fontsize=8.2,
               color="0.25")
    ax[1].set_xlabel(r"modulation frequency $\Omega/2\pi$  (MHz)")
    ax[1].set_ylabel(r"$|\delta n(\Omega)/\delta n(0)|$")
    ax[1].set_title("(b) Small-signal response: cavity-limited")
    ax[1].set_ylim(2e-3, 3.0)
    ax[1].legend(fontsize=8.4, loc="lower left")
    fig.tight_layout()
    _save(fig, "ring_validation")
    return np.array(T_num), np.array(T_an)


def fig_quantum_comb(dps):
    ppln, lu = dps["ppln"], dps["lu"]
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # -- (a) anharmonicity: dispersive cascade vs resonant chi2 --
    x = np.logspace(-0.5, 1.5, 300)
    for q, col, lab in ((ppln, COLORS["tcmt"], "target ring"),
                        (lu, COLORS["highlight"], "Lu 2020 device")):
        D = x * q.gamma2
        ax[0].loglog(x, anharmonicity_fom(q, D), color=col,
                     label=f"cascaded, {lab}")
        ax[0].axhline(resonant_fom(q), color=col, ls="--",
                      label=rf"resonant $g/\kappa_1$, {lab}")
    ax[0].axhline(1.0, color="k", ls=":", lw=1.2)
    ax[0].text(0.35, 1.4, "photon blockade", fontsize=9)
    ax[0].set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax[0].set_ylabel(r"single-photon anharmonicity")
    ax[0].set_title("(a) The cascade is not a route to blockade")
    ax[0].legend(fontsize=7.8, loc="lower left")

    # -- (b) which SH mode mediates which comb interaction: FSR walk-off --
    #  SFG of FM lines (m_p+q1, m_p+q2) lands on SH mode 2 m_p + j,
    #  j = q1+q2, whose detuning walks off as j (Omega_2 - Omega_1) with
    #  Omega_2/Omega_1 = ng1/ng2.  Only j = 0 (which includes the pump SPM
    #  and the degenerate FWM p+p -> s+i) stays resonant.
    j = np.arange(-6, 7)
    D0 = 2.0 * ppln.gamma2
    for ng2, col, lab in ((ppln.ng2, COLORS["tcmt"],
                           rf"$n_{{g,2}}={ppln.ng2}$ (as fabricated)"),
                          (ppln.ng1 * (1.0 + 6e-4), COLORS["second"],
                           r"$|\Delta n_g|/n_g=6\times10^{-4}$"
                           " (dispersion-engineered)")):
        q = ppln.replace(ng2=ng2)
        dOmega = q.FSR * (q.ng1 / ng2 - 1.0)      # Omega_2 - Omega_1
        Dj = D0 - j * dOmega
        ratio = np.abs(U_casc(q, Dj) / U_casc(q, D0))
        ax[1].semilogy(j, ratio, "o-", ms=5, color=col, label=lab)
        if ng2 == ppln.ng2:
            ax[1].annotate(rf"$|\Omega_2-\Omega_1|/2\pi={abs(dOmega)/TWOPI/1e9:.0f}$ GHz"
                           "\n" rf"vs $\gamma_2/2\pi={q.gamma2/TWOPI/1e9:.1f}$ GHz",
                           xy=(1.4, 2e-2), fontsize=8.6, color="0.3")
    ax[1].set_xlabel(r"SH mode index $j=q_1+q_2$ relative to $2m_p$")
    ax[1].set_ylabel(r"cascaded coupling / its $j=0$ value")
    ax[1].set_title("(b) Only one SH mode is resonant")
    ax[1].legend(fontsize=8.0, loc="lower center")
    ax[1].set_ylim(1e-2, 4.0)
    fig.tight_layout()
    _save(fig, "ring_quantum_comb")


def thermal_shift_full(p: RingParams, th: ThermalModel, P_in, Domega):
    """Thermal resonance shift including the heat the cascade itself makes.

    Linear load: the fraction of the dropped power that is absorbed.
    Cascaded load: the FM power funnelled into the SH and dissipated there,
    P_casc = hbar w1 * 2 kappa_casc n^2, which is EVEN in Domega.
    """
    D = np.atleast_1d(np.asarray(Domega, dtype=float))
    n = np.array([n_selfconsistent(p, P_in, 0.0, d) for d in D])
    kc = np.asarray(kappa_casc(p, D), dtype=float)
    P_lin = HBAR * p.omega1 * 2.0 * p.kappa1_0 * n          # FM intrinsic loss
    P_casc = HBAR * p.omega1 * 2.0 * kc * n**2              # into the SH channel
    # Only the INTRINSIC part of the SH decay stays on chip; the external
    # part leaves through the SH port.  Of what stays, only the absorptive
    # fraction becomes heat -- the rest is scattered out.
    P_casc_int = P_casc * (p.kappa2_0 / p.kappa2)
    P_abs = th.absorbed_fraction * (P_lin + P_casc_int)
    out = abs(th.dw_dT(p)) * P_abs * th.R_th
    return out if np.ndim(Domega) else float(out[0])


def fig_thermal(dps):
    p = dps["ppln"]
    th = ThermalModel()
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.1))

    # -- (a) magnitudes vs power --
    P = np.logspace(-7.5, -3.5, 300)
    D0 = 2.0 * p.gamma2
    n = np.array([n_selfconsistent(p, PP, 0.0, D0) for PP in P])
    n_lin = np.array([n_from_power(p, PP) for PP in P])
    kerr = abs(float(U_casc(p, D0))) * n
    nat = p.U_native * n_lin
    therm = np.array([thermal_shift_full(p, th, PP, D0) for PP in P])
    ax[0].loglog(P * 1e6, kerr / p.kappa1, color=COLORS["tcmt"],
                 label="cascaded Kerr")
    ax[0].loglog(P * 1e6, nat / p.kappa1, color=COLORS["second"], ls="--",
                 label="native Kerr")
    ax[0].loglog(P * 1e6, therm / p.kappa1, color=COLORS["wrong"], ls="-.",
                 label=r"thermal ($f_{\rm abs}=0.3$)")
    th_lo = ThermalModel(absorbed_fraction=0.05)
    ax[0].loglog(P * 1e6,
                 np.array([thermal_shift_full(p, th_lo, PP, D0) for PP in P])
                 / p.kappa1, color=COLORS["wrong"], ls=":", lw=1.2,
                 label=r"thermal ($f_{\rm abs}=0.05$)")
    ax[0].axhline(np.sqrt(3), color="k", ls="-", lw=0.9)
    ax[0].text(2e-2, 2.3, r"bistability threshold $\sqrt{3}\kappa_1$",
               fontsize=8.4)
    ax[0].set_xlabel(r"input power  ($\mu$W)")
    ax[0].set_ylabel(r"$|\delta\omega_1|/\kappa_1$")
    ax[0].set_title("(a) Thermal background is not negligible")
    ax[0].legend(fontsize=8.2, loc="lower right")

    # -- (b) the discriminator: antisymmetric part in Domega --
    x = np.linspace(-6, 6, 601)          # symmetric grid, odd length
    D = x * p.gamma2
    P0 = 20e-6
    n0 = np.array([n_selfconsistent(p, P0, 0.0, d) for d in D])
    kerr = -np.asarray(U_casc(p, D), dtype=float) * n0     # shift = -U n
    therm = -thermal_shift_full(p, th, P0, D)              # always red
    total = kerr + therm
    anti = 0.5 * (total - total[::-1])
    ax[1].plot(x, kerr / p.kappa1, color=COLORS["tcmt"], lw=3.0, alpha=0.35,
               label=r"cascaded Kerr (odd in $\Delta\omega$)")
    ax[1].plot(x, therm / p.kappa1, color=COLORS["wrong"], ls="-.",
               label=r"thermal (even, always red)")
    ax[1].plot(x, total / p.kappa1, color="k", ls=":",
               label="measured total")
    ax[1].plot(x, anti / p.kappa1, color=COLORS["highlight"], lw=1.4,
               label=r"antisymmetric part of the total")
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax[1].set_ylabel(r"$\delta\omega_1/\kappa_1$")
    ax[1].set_title(rf"(b) Isolating the cascade at $P_{{\rm in}}={P0*1e6:.0f}\,\mu$W")
    ax[1].legend(fontsize=8.0, loc="lower right")
    resid = np.max(np.abs(anti - kerr)) / np.max(np.abs(kerr))
    fig.tight_layout()
    _save(fig, "ring_thermal")
    return resid



# =====================================================================
def summary(dps):
    print("\n" + "=" * 78)
    print("MICRORING CASCADE: SUMMARY TABLE  (Domega = 2 gamma2)")
    print("=" * 78)
    hdr = (f"{'design':<30}{'g/2pi':>9}{'Q2':>9}{'R':>9}"
           f"{'U_c/2pi':>11}{'P_bist':>11}{'P_MI':>10}")
    print(hdr)
    print(f"{'':<30}{'[MHz]':>9}{'':>9}{'':>9}{'[Hz/ph]':>11}"
          f"{'[uW]':>11}{'[uW]':>10}")
    print("-" * 78)
    rows = {}
    for key in ("lu", "ppln", "modal"):
        q = dps[key]
        D = 2.0 * q.gamma2
        U = float(U_casc(q, D))
        kc = float(kappa_casc(q, D))
        Ut = U + q.U_native
        Pb = bistability_threshold_exact(q, D)[0]
        Pm = mi_threshold_bogoliubov(q, D, 10e6 * TWOPI, q_max=40,
                                     d1_grid=81)[0]
        print(f"{q.label:<30}{q.g_rate/TWOPI/1e6:>9.2f}{q.Q2_load:>9.0f}"
              f"{float(R_ratio(q, D)):>9.0f}{U/TWOPI:>11.0f}"
              f"{Pb*1e6:>11.2f}{Pm*1e6:>10.2f}")
        rows[key] = dict(R=float(R_ratio(q, D)), U=U, Pb=Pb, Pm=Pm, kc=kc)
    # native reference on the target ring
    q = dps["ppln"]
    Pb_n = bistability_threshold_power(q, q.U_native, 0.0)  # exact: no Deff dep.
    Pm_n = (q.kappa1 / q.U_native) * HBAR * q.omega1 * q.kappa1
    print("-" * 78)
    print(f"{'native Kerr only (same ring)':<30}{'-':>9}{'-':>9}{1:>9.0f}"
          f"{q.U_native/TWOPI:>11.2f}{Pb_n*1e6:>11.1f}{Pm_n*1e6:>10.1f}")
    print("=" * 78)

    q = dps["ppln"]
    D = 2.0 * q.gamma2
    U = float(U_casc(q, D)) + q.U_native
    kc = float(kappa_casc(q, D))
    nc, shift = stability_ceiling(q, D)
    print("\nTarget ring, further numbers:")
    print(f"    n2 equivalent            : {float(n2_casc_eff(q, D)):.3e} m^2/W"
          f"  ({float(n2_casc_eff(q, D))/q.n2_kerr:+.0f} x native)")
    print(f"    chi3 equivalent          : {float(chi3_casc_eff(q, D)):.3e}"
          f" m^2/V^2")
    print(f"    switching photons        : {switching_photons(q, U):.3e}"
          f"   (native {switching_photons(q, q.U_native):.3e})")
    n_sw = switching_photons(q, U)
    print(f"    switching drive: linear-equiv {power_from_n(q, n_sw)*1e6:.2f} uW, "
          f"self-consistent at delta1=0 {_drive_for_photons(q, n_sw, D)*1e6:.2f} uW")
    print(f"    intracavity switch energy: "
          f"{HBAR*q.omega1*switching_photons(q, U)*1e15:.2f} fJ"
          f"   (native {HBAR*q.omega1*switching_photons(q, q.U_native)*1e15:.1f} fJ)")
    print(f"    stability ceiling n_c    : {nc:.3e} photons"
          f"  -> max shift {shift:.1f} kappa1")
    print(f"    analytic n_sat           : {saturation_photons(q, D):.3e}")
    print(f"    ceiling input power      : {power_from_n(q, nc)*1e6:.1f} uW")
    print(f"    dispersive anharmonicity : {float(anharmonicity_fom(q, D)):.2e}"
          f"   (resonant g/kappa1 = {resonant_fom(q):.2e})")
    print(f"    cascaded Kerr bandwidth  : {np.sqrt(D**2+q.kappa2**2)/TWOPI/1e9:.2f} GHz"
          f"   (FM cavity {q.kappa1/TWOPI/1e6:.0f} MHz -> cavity-limited)")
    th = ThermalModel()
    xs = np.linspace(0.9, 12.0, 600)
    Pb, Pm = threshold_vs_detuning(q, xs)          # mi=False: Pm is None (Bogoliubov scan is slow)
    ib = int(np.argmin(Pb))
    print(f"    optimal Domega (bistab.) : {xs[ib]:.2f} gamma2 -> "
          f"{Pb[ib]*1e6:.1f} uW  (at 2 gamma2: {Pb[np.argmin(abs(xs-2))]*1e6:.1f} uW)")
    if Pm is not None:
        im = int(np.argmin(Pm))
        print(f"    optimal Domega (paramet.): {xs[im]:.2f} gamma2 -> "
              f"{Pm[im]*1e6:.2f} uW")
    n20 = n_selfconsistent(q, 20e-6, 0.0, D)
    print(f"    n at 20 uW (self-cons.)  : {n20:.3e}"
          f"   (linear formula {n_from_power(q, 20e-6):.3e})")
    print(f"    cascaded loss kc*n/kappa1: "
          f"{float(kappa_casc(q, D))*n20/q.kappa1:.2f}"
          f"   -> linewidth broadening factor {1+float(kappa_casc(q, D))*n20/q.kappa1:.2f}")
    print(f"    thermal shift at 20 uW   : "
          f"{thermal_shift_full(q, th, 20e-6, D)/q.kappa1:.2f} kappa1"
          f"   (cascade {abs(float(U_casc(q, D)))*n20/q.kappa1:.2f} kappa1)")


def main():
    dps = design_points()
    print("Generating microring figures ...")
    fig_R_landscape(dps)
    fig_sign_flip(dps)
    xs_th, Pb_th = fig_thresholds(dps)
    fig_bistability(dps)
    xs, n_c, an, sh = fig_stability(dps)
    T_num, T_an = fig_validation(dps)
    fig_quantum_comb(dps)
    resid_anti = fig_thermal(dps)

    print("\nValidation residuals:")
    err = np.max(np.abs(T_num - T_an))
    print(f"    max |T_TCMT - T_closedform| = {err:.3e}")
    ratio = np.array(n_c) / np.array(an)
    print(f"    antisym. extraction recovers Kerr to {resid_anti:.2e} rel.")
    qnw = dps["ppln"].replace(ng2=dps["ppln"].ng1)
    D_ = 2 * qnw.gamma2
    d1_ = -3.77 * qnw.kappa1
    Deff_ = D_ - 2 * d1_
    U_ = float(U_casc(qnw, Deff_)) + qnw.U_native
    kc_ = float(kappa_casc(qnw, Deff_))
    aa = -3 * (U_**2 + kc_**2); bb = 4 * (U_ * d1_ - qnw.kappa1 * kc_)
    cc = -(d1_**2 + qnw.kappa1**2)
    rr = np.roots([aa, bb, cc]); rr = rr[np.abs(rr.imag) < 1e-9].real
    n_an = rr[rr > 0].min()
    lo_, hi_ = n_an * 1e-3, n_an * 1e3
    for _ in range(200):
        m_ = np.sqrt(lo_ * hi_)
        if mi_gain(qnw, m_, d1_, D_, 1, 0.0) > 0: hi_ = m_
        else: lo_ = m_
    print(f"    4x4 Bogoliubov vs analytic 2x2 reduction (zero walk-off): "
          f"{abs(np.sqrt(lo_*hi_) - n_an)/n_an:.2e} rel.")
    print(f"    numerical ceiling / analytic n_sat: min {ratio.min():.3f}, "
          f"max {ratio.max():.3f}, median {np.median(ratio):.3f}")

    summary(dps)


if __name__ == "__main__":
    main()
