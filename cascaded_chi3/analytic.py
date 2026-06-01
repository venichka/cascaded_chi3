"""Closed-form expressions from Section 5 of the note (dimensionless,
beta = gamma2 = 1, so |beta|^2 = 1).

Adiabatic elimination of a2 (set da2/dtau = 0) gives

    a2_adiabatic = i * a1^2 / (1/2 - i*Domega)

and the resulting effective FM equation contains a complex coefficient on
|a1|^2 * a1:

    kappa_casc = -1 / (1/2 - i*Domega)
               = -(1/2 + i*Domega) / (Domega^2 + 1/4).

Writing the cascaded contribution as (i*alpha_eff_casc + r) * |a1|^2 * a1
with real alpha_eff_casc (the cascaded Kerr) and real r (cascaded loss):

    alpha_eff_casc = -Domega / (Domega^2 + 1/4)
    r              = -(1/2)  / (Domega^2 + 1/4)   (i.e. two-photon-like loss)

Far off-resonance (|Domega| >> 1/2):
    alpha_eff_casc -> -1 / Domega        (the often-quoted simple form)
"""

import numpy as np


def alpha3_eff_casc_full(Domega):
    """Full adiabatic-elimination cascaded Kerr coefficient (real)."""
    return -Domega / (np.asarray(Domega) ** 2 + 0.25)


def alpha3_eff_casc_simple(Domega):
    """Leading |Domega| >> 1/2 form quoted in the note."""
    return -1.0 / np.asarray(Domega)


def gamma_2pa_eff(Domega):
    """Cascaded two-photon-like loss (positive coefficient on |a1|^2)."""
    return 0.5 / (np.asarray(Domega) ** 2 + 0.25)


def kappa_casc(Domega):
    """Complex coefficient on |a1|^2 * a1 from adiabatic SH elimination."""
    return -1.0 / (0.5 - 1j * np.asarray(Domega))


def a2_adiabatic(a1, Domega):
    """Adiabatic-elimination steady-state SH amplitude given a1."""
    return 1j * a1 * a1 / (0.5 - 1j * Domega)


def cw_kerr_only_intracavity(p, s_plus):
    """Real positive roots of the static Kerr-only model.

    The cascaded coefficient is set by the *pump-frame* SH detuning
    Delta_SH = 2*delta1 + Domega, NOT by Domega alone (the latter is the
    SH detuning at delta1 = 0).

    Solves |a1|^2 * [(g1 + g_eff*|a1|^2)^2 + (delta1 + a_eff*|a1|^2)^2]
                = gamma1_rad * |s_plus|^2
    with a_eff = alpha3 + alpha_eff_casc(Delta_SH), g_eff = gamma_2pa_eff(Delta_SH),
    g1 = gamma1/2.  Returns sorted ndarray of nonneg real |a1|^2 roots.
    """
    Delta_SH = 2 * p.delta1 + p.Domega
    # beta = 0 means the SH cavity is decoupled -- the cascade adds no
    # contribution, only the bare Kerr alpha_3 acts on the FM. This is
    # the "native" architecture.
    beta_sq = float(getattr(p, "beta", 1.0)) ** 2
    a_eff = p.alpha3 + beta_sq * alpha3_eff_casc_full(Delta_SH)
    g_eff = beta_sq * gamma_2pa_eff(Delta_SH)
    g1 = 0.5 * p.gamma1
    rhs = p.gamma1_rad * abs(s_plus) ** 2

    # cubic in x = |a1|^2:
    c3 = g_eff**2 + a_eff**2
    c2 = 2 * (g1 * g_eff + p.delta1 * a_eff)
    c1 = g1**2 + p.delta1**2
    coeffs = [c3, c2, c1, -rhs]

    roots = np.roots(coeffs)
    real = roots[np.abs(roots.imag) < 1e-8 * (1 + np.abs(roots))].real
    real = real[real > -1e-12]
    real = np.clip(real, 0.0, None)
    return np.sort(real)
