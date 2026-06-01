"""Linear stability of the trivial fixed point of the TCMT.

Around (a1_FP, a2_FP) = (a1_FP, i*beta*a1_FP^2/(gamma2/2 - i*Delta_SH)),
small perturbations (delta_a1, delta_a2) obey a 4x4 linear system in
the real-packed state (Re a1, Re a2, Im a1, Im a2). The cascaded
parametric mixing comes from the term i*beta^* a1^* a2 in the FM
equation, which couples delta_a1 to delta_a1^* via the SH amplitude
a2_FP -- the standard parametric instability mechanism.

This module returns the eigenvalues of the Jacobian as a function of
the saturation parameter sat = |beta a1_FP|/|Delta_SH|, at fixed
gamma_1/gamma_2, gamma_1_rad/gamma_1, Delta_SH/gamma_2, delta_p/gamma_2.
The threshold is the smallest sat for which max(Re lambda) > 0.
"""

from __future__ import annotations

import numpy as np

from analytic import a2_adiabatic
from params import TCMTParams


def jacobian(p: TCMTParams, a1_fp: complex, a2_fp: complex) -> np.ndarray:
    """4x4 real Jacobian of the TCMT around (a1_fp, a2_fp).

    State y = [Re a1, Re a2, Im a1, Im a2]. Returns dF/dy where F is
    the RHS of tcmt._rhs at constant pump (with the constant +sqrt(g1r)*s
    dropped, since the Jacobian is unaffected).
    """
    x1, y1 = a1_fp.real, a1_fp.imag
    x2, y2 = a2_fp.real, a2_fp.imag
    g1 = 0.5 * p.gamma1
    a3 = p.alpha3
    n1 = x1**2 + y1**2  # |a1|^2

    # FM RHS in real form (from tcmt._rhs):
    # da1 = (i*delta1 - g1)*a1 + i*conj(a1)*a2 + i*alpha3*|a1|^2 a1
    # Real part : -g1*x1 - delta1*y1 - (x1*y2 - y1*x2) - a3*n1*y1
    # Imag part : -g1*y1 + delta1*x1 + (x1*x2 + y1*y2) + a3*n1*x1

    # Partial derivatives w.r.t. (x1, x2, y1, y2)
    # F1 = -g1*x1 - delta1*y1 - x1*y2 + y1*x2 - a3*(x1^2+y1^2)*y1
    # F2 = -0.5*x2 - (2*delta1+Domega)*y2 - 2*x1*y1                 (SH real)
    # F3 = -g1*y1 + delta1*x1 + x1*x2 + y1*y2 + a3*(x1^2+y1^2)*x1
    # F4 = -0.5*y2 + (2*delta1+Domega)*x2 + x1^2 - y1^2             (SH imag)

    d1 = p.delta1
    Dsh = 2 * p.delta1 + p.Domega

    J = np.zeros((4, 4))
    # dF1/d.
    J[0, 0] = -g1 - y2 - 2 * a3 * x1 * y1
    J[0, 1] = y1
    J[0, 2] = -d1 + x2 - a3 * (n1 + 2 * y1**2)
    J[0, 3] = -x1
    # dF2/d.
    J[1, 0] = -2 * y1
    J[1, 1] = -0.5
    J[1, 2] = -2 * x1
    J[1, 3] = -Dsh
    # dF3/d.
    J[2, 0] = d1 + x2 + a3 * (n1 + 2 * x1**2)
    J[2, 1] = x1
    J[2, 2] = -g1 + y2 + 2 * a3 * x1 * y1
    J[2, 3] = y1
    # dF4/d.
    J[3, 0] = 2 * x1
    J[3, 1] = Dsh
    J[3, 2] = -2 * y1
    J[3, 3] = -0.5
    return J


def stability_vs_sat(p: TCMTParams, sat_grid):
    """For each sat parameter in sat_grid, build the trivial FP from the
    cubic, build the Jacobian, return max(Re eigenvalue).

    sat = |beta a1| / |Delta_SH| with beta = 1, gamma_2 = 1.
    Pick a1_FP real positive at this magnitude (the gauge in which the
    pump is real-positive at steady state, on FM resonance).
    """
    Dsh = 2 * p.delta1 + p.Domega
    out = np.empty_like(sat_grid)
    for i, sat in enumerate(sat_grid):
        a1_fp = abs(Dsh) * sat + 0j  # real-positive, |a1| = sat * |Delta_SH|
        a2_fp = a2_adiabatic(a1_fp, Dsh)
        J = jacobian(p, a1_fp, a2_fp)
        eigs = np.linalg.eigvals(J)
        out[i] = np.max(eigs.real)
    return out


def find_threshold(p: TCMTParams,
                   sat_lo: float = 0.5, sat_hi: float = 3.0,
                   n: int = 400):
    """Smallest sat in [sat_lo, sat_hi] where max(Re lambda) crosses 0."""
    sat = np.linspace(sat_lo, sat_hi, n)
    re_max = stability_vs_sat(p, sat)
    crossings = np.where(np.diff(np.sign(re_max)) > 0)[0]
    if not crossings.size:
        return None, sat, re_max
    i = crossings[0]
    sat_thr = np.interp(0.0, [re_max[i], re_max[i + 1]], [sat[i], sat[i + 1]])
    return sat_thr, sat, re_max
