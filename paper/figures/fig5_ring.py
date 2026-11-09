"""Fig. 5 (draft): microrings.

(a) R = |U_casc|/U_nat versus Q2 at a cold mismatch of 2 gamma2 for the
    poled ring (design B of Supplement 1, Table S2) with g scaled from the
    fitted rate of Lu 2020 (upper solid) or from its measured SHG efficiency
    (lower solid), and for the unpoled ring (design C, dashed), all with
    Q2_int = 2 Q2.  Stars: design A (the device of Lu 2020) with the fitted
    g (filled) and with g from its efficiency (open).  Diamond: the AlN ring
    of Cui 2025, from its measured extremes of the self-phase-modulation
    tuning (1 -+ 3.7 at |Domega| = gamma2/2) scaled to Domega = 2 gamma2 with
    the Lorentzian of U_casc (factor 8/17).
(b) Cascaded shift of the pump resonance of design A versus the pump-frame
    Domega at 5 uW input, with the pump locked to the shifted resonance:
    main-text Eq. (7) from the measured efficiency (solid) and from the
    efficiency implied by the fitted g (dashed); dots are the locked-pump
    steady state with the nonlinear loss included.

Model: cascaded_chi3/ring_params.py (photon units, kappa = half-width).
"""

import numpy as np
import matplotlib.pyplot as plt

from paperstyle import apply, WIDTH, BLUE, RED, GREEN, GREY, panel_label, save
import ring_params as rp

TW = 2 * np.pi
HBAR = 1.054571817e-34


def shift_eq7(p, eta, P, Dw):
    """Main-text Eq. (7): first-order cascaded shift [rad/s] of the pump
    resonance, from the SHG efficiency eta [1/W] and input power P [W]."""
    k1, k2 = p.kappa1, p.kappa2
    cpl = (k1 / (2 * p.kappa1_ex)) * (k2 / (2 * p.kappa2_ex))
    return eta * P * cpl * k1 * k2 * Dw / (Dw**2 + k2**2)


def shift_locked(p, g, P, Dw):
    """Cascaded shift -U_casc n [rad/s] with the pump locked to the shifted
    resonance; n solves n (k1 + k_casc n)^2 = 2 k1e P / (hbar w1)."""
    k1, k2 = p.kappa1, p.kappa2
    U = -2 * g**2 * Dw / (Dw**2 + k2**2)
    kc = 2 * g**2 * k2 / (Dw**2 + k2**2)
    drive = 2 * p.kappa1_ex * P / (HBAR * p.omega1)
    n = drive / k1**2
    for _ in range(200):
        n = drive / (k1 + kc * n) ** 2
    return -U * n


def main():
    apply()
    dps = rp.design_points()
    A, B, Cc = dps["lu"], dps["ppln"], dps["modal"]
    g_fit = A.g_rate
    g_eta = rp.g_from_eta_norm(A, A.eta_norm_ref)
    scale = g_eta / g_fit
    print(f"design A: g_fit/2pi = {g_fit/TW/1e6:.3f} MHz, "
          f"g_eta/2pi = {g_eta/TW/1e6:.3f} MHz")

    fig = plt.figure(figsize=(WIDTH, 2.05))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    # (a) R vs Q2
    Q2 = np.logspace(4, 6.1, 200)

    def R_curve(q0, gfac=1.0):
        out = []
        for Q in Q2:
            q = q0.replace(Q2_load=Q, Q2_int=2 * Q, g=q0.g_rate * gfac)
            out.append(float(rp.R_ratio(q, 2 * q.gamma2)))
        return np.array(out)

    R_hi, R_lo = R_curve(B), R_curve(B, scale)
    ax_a.fill_between(Q2, R_lo, R_hi, color=BLUE, alpha=0.15, lw=0)
    ax_a.loglog(Q2, R_hi, "-", color=BLUE)
    ax_a.loglog(Q2, R_lo, "-", color=BLUE, lw=0.7)
    ax_a.loglog(Q2, R_curve(Cc), "--", color=BLUE)
    RA_fit = float(rp.R_ratio(A, 2 * A.gamma2))
    RA_eta = float(rp.R_ratio(A.replace(g=g_eta), 2 * A.gamma2))
    ax_a.plot([A.Q2_load] * 2, [RA_eta, RA_fit], "-", color=RED, lw=0.6)
    ax_a.plot([A.Q2_load], [RA_fit], "*", color=RED, ms=7, mec="k", mew=0.4)
    ax_a.plot([A.Q2_load], [RA_eta], "*", color="white", ms=7, mec=RED,
              mew=0.7)
    R_aln = 3.7 * 8 / 17
    ax_a.plot([1.83e5], [R_aln], "D", color=GREEN, ms=4, mec="k", mew=0.4)
    print(f"R: A {RA_eta:.0f}-{RA_fit:.0f}, B {R_lo[np.argmin(abs(Q2-3e5))]:.0f}"
          f"-{R_hi[np.argmin(abs(Q2-3e5))]:.0f}, AlN {R_aln:.2f}")
    ax_a.axhline(1.0, color="0.75", lw=0.5, zorder=0)
    ax_a.set_xlabel(r"$Q_2$")
    ax_a.set_ylabel(r"$R$")
    ax_a.set_ylim(0.3, 3e4)
    panel_label(ax_a, "(a)")

    # (b) predicted shift of the device of Lu 2020 at 5 uW
    P0 = 5e-6
    x = np.linspace(-6, 6, 601)
    Dw = x * A.gamma2
    eta_fit = rp.eta_norm_from_g(A, g_fit)
    ax_b.axhline(0.0, color="0.75", lw=0.5, zorder=0)
    ax_b.axvline(0.0, color="0.75", lw=0.5, zorder=0)
    ax_b.plot(x, shift_eq7(A, A.eta_norm_ref, P0, Dw) / TW / 1e6, "-",
              color=BLUE)
    ax_b.plot(x, shift_eq7(A, eta_fit, P0, Dw) / TW / 1e6, "--", color=BLUE)
    xd = np.linspace(-6, 6, 25)
    for g in (g_eta, g_fit):
        ax_b.plot(xd, shift_locked(A, g, P0, xd * A.gamma2) / TW / 1e6, "o",
                  color=BLUE, ms=2.2, mfc="white", mew=0.6)
    k1 = A.kappa1
    s2 = shift_eq7(A, A.eta_norm_ref, P0, 2 * A.gamma2)
    print(f"design A at 5 uW, Domega = 2 gamma2: shift/k1 = {s2/k1:.3f}, "
          f"{s2/TW/1e6:.1f} MHz (kappa1/2pi = {k1/TW/1e6:.0f} MHz); "
          f"fitted-g ratio {eta_fit/A.eta_norm_ref:.2f}; "
          f"native {A.U_native*rp.n_from_power(A, P0)/TW/1e6:.3f} MHz")
    ax_b.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax_b.set_ylabel(r"cascaded shift $\delta\omega_{\mathrm{casc}}/2\pi$ (MHz)")
    ax_b.set_xlim(-6, 6)
    panel_label(ax_b, "(b)")

    fig.tight_layout(w_pad=1.0)
    save(fig, "fig5_ring")


if __name__ == "__main__":
    main()
