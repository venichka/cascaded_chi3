"""Fig. S1 of Supplement 1: checks of the cascaded susceptibility.

(a) Cascaded coefficient alpha_casc of the pump mode, real (Kerr) and
    imaginary (loss) parts, from steady-state integration of the two-mode
    equations (markers) and from the closed form
    alpha_casc = -beta1* beta / (Domega + i gamma2/2) (lines).
    Units: beta = gamma2 = 1 (N1 = N2), so alpha_casc is in units of
    |beta|^2 / gamma2.  gamma1 = 10 gamma2, gamma1^e = gamma1/2,
    delta1 = 0, alpha3 = 0, weak drive.
(b) Cascaded susceptibility of a one-dimensional two-layer Bloch cell,
    computed from the full driven second-harmonic field (finite
    differences, no modal expansion) and divided by the single-mode closed
    form, versus Domega / 2 omega1.  Dashed line: the single-mode pole
    without linearisation, 1 / (1 - Domega / 4 omega1).

Run:  conda run -n ml_base python fig_s1_checks.py
Writes fig_s1_checks.pdf and .png next to this script.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.linalg import eigh
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "cascaded_chi3"))
from params import TCMTParams          # noqa: E402
from tcmt import steady_state          # noqa: E402

# ------------------------------------------------------------------ style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Palatino", "Palatino Linotype", "TeX Gyre Pagella",
                   "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.1,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})
BLUE, RED, GREY = "#1f4e9c", "#b2352b", "0.45"


# --------------------------------------------------------- panel (a) data
def panel_a():
    Dw = np.linspace(-6.0, 6.0, 25)
    base = dict(gamma1=10.0, gamma1_rad=5.0, alpha3=0.0, delta1=0.0)
    s_plus = 0.01
    re_num, im_num = np.empty_like(Dw), np.empty_like(Dw)
    for i, d in enumerate(Dw):
        p = TCMTParams(Domega=float(d), **base)
        a1, a2 = steady_state(p, s_plus, t_settle=200.0)
        # i conj(a1) a2 = i alpha_casc |a1|^2 a1  with beta = 1
        alpha = np.conj(a1) * a2 / (abs(a1) ** 2 * a1)
        re_num[i], im_num[i] = alpha.real, alpha.imag
    Dw_fine = np.linspace(-6.0, 6.0, 801)
    alpha_cf = -1.0 / (Dw_fine + 0.5j)
    return Dw, re_num, im_num, Dw_fine, alpha_cf


# --------------------------------------------------------- panel (b) data
L, NG = 1.0, 1200
Z = np.arange(NG) * L / NG
DZ = L / NG
LAYER_A = Z < 0.45 * L
EPS_A, EPS_B = 2.2 ** 2, 1.45 ** 2
CHI = np.where(LAYER_A, 1.0, 0.0)


def lap_bloch(k):
    D = np.zeros((NG, NG), complex)
    i = np.arange(NG)
    D[i, i] = 2.0 / DZ ** 2
    D[i[:-1], i[:-1] + 1] = -1.0 / DZ ** 2
    D[i[1:], i[1:] - 1] = -1.0 / DZ ** 2
    ph = np.exp(1j * k * L)
    D[NG - 1, 0] = -ph / DZ ** 2
    D[0, NG - 1] = -np.conj(ph) / DZ ** 2
    return D


def modes(k, eps, nmax=8):
    w2, V = eigh(lap_bloch(k), np.diag(eps), subset_by_index=[0, nmax - 1])
    return np.sqrt(w2.clip(0)), V


def integ(f):
    return np.sum(f) * DZ


def sh_eps(d):
    """SH permittivity, detuned by a scaling d that sets Domega."""
    return np.where(LAYER_A, EPS_A * (1 + d), EPS_B * (1 + 0.3 * d))


def panel_b():
    k1 = 0.37 * np.pi / L
    eps1 = np.where(LAYER_A, EPS_A, EPS_B)
    w_all, V1 = modes(k1, eps1)
    w1, f1 = w_all[0], V1[:, 0]
    W = 2 * w1
    w0, _ = modes(2 * k1, sh_eps(0.0))
    n_sh = int(np.argmin(abs(w0 - W)))

    def mismatch(d, target):
        w, _ = modes(2 * k1, sh_eps(d))
        return (W - w[n_sh]) / W - target

    targets = np.array([-2e-2, -1e-2, -5e-3, -2e-3, -1e-3,
                        1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
    ratio = []
    for tgt in targets:
        d = brentq(mismatch, -0.3, 0.3, args=(tgt,), xtol=1e-13)
        eps2 = sh_eps(d)
        w, V = modes(2 * k1, eps2)
        w2, f2 = w[n_sh], V[:, n_sh]
        Dw = W - w2
        N2 = integ(eps2 * abs(f2) ** 2)
        M = integ(CHI * np.conj(f2) * f1 ** 2)
        I4 = integ(abs(f1) ** 4)
        A = sp.csr_matrix(lap_bloch(2 * k1) - W ** 2 * np.diag(eps2))
        E2 = spl.spsolve(A, W ** 2 * CHI * f1 ** 2)
        P_casc = 2 * CHI * E2 * np.conj(f1)
        chi3_full = integ(np.conj(f1) * P_casc) / (3 * I4)
        chi3_cf = -(2 * w1 / 3) * abs(M) ** 2 / (N2 * Dw * I4)
        ratio.append(chi3_full.real / chi3_cf)
    return targets, np.array(ratio)


def main():
    Dw, re_num, im_num, Dw_fine, alpha_cf = panel_a()
    x, ratio = panel_b()
    for xi, ri in zip(x, ratio):
        print(f"Dw/2w1 = {xi:+.0e}   full/closed form = {ri:.6f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.25, 2.05))

    ax1.axhline(0, color="0.75", lw=0.5, zorder=0)
    ax1.axvline(0, color="0.75", lw=0.5, zorder=0)
    ax1.plot(Dw_fine, alpha_cf.real, "-", color=BLUE)
    ax1.plot(Dw_fine, alpha_cf.imag, "--", color=RED)
    ax1.plot(Dw, re_num, "o", ms=3.0, mfc="white", mec=BLUE, mew=0.8)
    ax1.plot(Dw, im_num, "s", ms=2.8, mfc="white", mec=RED, mew=0.8)
    ax1.set_xlabel(r"$\Delta\omega/\gamma_2$")
    ax1.set_ylabel(r"$\alpha_{\mathrm{casc}}\,\gamma_2/|\beta|^2$")
    ax1.set_xlim(-6, 6)
    ax1.set_ylim(-1.15, 2.15)
    ax1.text(-5.6, 0.32, r"Re", color=BLUE)
    ax1.text(0.55, 1.75, r"Im", color=RED)
    ax1.text(0.03, 0.96, "(a)", transform=ax1.transAxes, va="top")

    xs = np.linspace(-0.022, 0.022, 200)
    ax2.axhline(1, color="0.75", lw=0.5, zorder=0)
    ax2.plot(100 * xs, 1 / (1 - xs / 2), "--", color=GREY)
    ax2.plot(100 * x, ratio, "o", ms=3.0, mfc="white", mec=BLUE, mew=0.8)
    ax2.set_xlabel(r"$10^{2}\,\Delta\omega/2\omega_1$")
    ax2.set_ylabel(r"$\chi^{\mathrm{full}}_{\mathrm{casc}}/\chi_{\mathrm{casc}}$")
    ax2.set_xlim(-2.2, 2.2)
    ax2.text(0.03, 0.96, "(b)", transform=ax2.transAxes, va="top")

    fig.tight_layout(w_pad=1.2)
    out = HERE / "fig_s1_checks"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)
    print(f"wrote {out}.pdf/.png")


if __name__ == "__main__":
    main()
