"""Two-mode TCMT for cascaded chi^(2):chi^(2) -> chi^(3)_eff.

Equations (dimensionless; beta = gamma2 = 1; rotating frames at omega_p
for the FM and 2*omega_p for the SH):

    da1/dtau = (i*delta1 - 0.5*gamma1) * a1
               + i * conj(a1) * a2
               + i * alpha3 * |a1|^2 * a1
               + sqrt(gamma1_rad) * s_plus
    da2/dtau = (i*(2*delta1 + Domega) - 0.5) * a2
               + i * a1**2

with delta1 = (omega_p - omega1)/gamma2 the pump detuning from the FM
resonance and Domega = (2*omega1 - omega2)/gamma2 the SH detuning. The
input-output relation at the FM port is s_minus = -s_plus + sqrt(gamma1_rad)*a1.

scipy.solve_ivp does not handle complex state, so the integrator is fed a
real-packed state y = [Re(a1), Re(a2), Im(a1), Im(a2)].
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from params import TCMTParams


def _rhs(t, y, p: TCMTParams, sp):
    a1 = y[0] + 1j * y[2]
    a2 = y[1] + 1j * y[3]
    s = sp(t) if callable(sp) else sp
    Domega_t = p.Domega(t) if callable(p.Domega) else p.Domega
    delta1_t = p.delta1(t) if callable(p.delta1) else p.delta1

    da1 = (
        (1j * delta1_t - 0.5 * p.gamma1) * a1
        + 1j * p.beta * np.conj(a1) * a2
        + 1j * p.alpha3 * (a1.real**2 + a1.imag**2) * a1
        + np.sqrt(p.gamma1_rad) * s
    )
    da2 = (
        (1j * (2 * delta1_t + Domega_t) - 0.5) * a2
        + 1j * p.beta * a1 * a1
    )
    return [da1.real, da2.real, da1.imag, da2.imag]


def _pack(a1: complex, a2: complex):
    return np.array([a1.real, a2.real, a1.imag, a2.imag], dtype=float)


def _unpack(sol_y: np.ndarray) -> np.ndarray:
    a = np.empty((2, sol_y.shape[1]), dtype=complex)
    a[0] = sol_y[0] + 1j * sol_y[2]
    a[1] = sol_y[1] + 1j * sol_y[3]
    return a


def integrate(
    p: TCMTParams,
    sp,
    t_span,
    y0=(0 + 0j, 0 + 0j),
    t_eval=None,
    method: str = "DOP853",
    rtol: float = 1e-10,
    atol: float = 1e-13,
    max_step=np.inf,
):
    """Integrate the TCMT. sp is constant, scalar callable sp(t), or 0.

    Returns the scipy OdeResult with an extra attribute `.a` of shape
    (2, len(t)) holding complex (a1, a2) trajectories.
    """
    a1_0, a2_0 = y0
    y0_real = _pack(complex(a1_0), complex(a2_0))
    sol = solve_ivp(
        _rhs,
        t_span,
        y0_real,
        args=(p, sp),
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    sol.a = _unpack(sol.y)
    return sol


def steady_state(
    p: TCMTParams,
    sp_const: complex,
    t_settle: float | None = None,
    **kwargs,
):
    """Integrate at constant pump until transients die; return (a1, a2)."""
    if t_settle is None:
        # Slowest mode is whichever of the two has the smaller decay rate.
        rate = min(p.gamma1, 1.0)
        t_settle = 30.0 / rate
    sol = integrate(p, sp_const, (0.0, t_settle), **kwargs)
    return sol.a[0, -1], sol.a[1, -1]


def attractor(
    p: TCMTParams,
    sp_const: complex,
    t_settle: float = 250.0,
    t_avg: float = 250.0,
    n_samples: int = 1200,
    **kwargs,
):
    """Integrate the full TCMT to the attractor, then time-average over
    the last t_avg interval.

    Works regardless of whether the trivial fixed point is stable
    (regime A: returns the fixed-point value) or has lost stability to
    parametric oscillation (regime B: returns the limit-cycle time
    average).  In either case the returned values are well-defined
    experimentally observable quantities.

    Returns a dict with:
        n1_mean, n1_std  : <|a_1|^2> and its std over t_avg
        n2_mean, n2_std  : <|a_2|^2> and its std
        T_mean, T_std    : <|s_-/s_+|^2> and its std (FM-port transmission)
        a1_last          : final-sample complex a_1 (snapshot, not avg)
        is_limit_cycle   : bool, True if std/mean exceeds 1% (heuristic)
    """
    max_step = kwargs.pop("max_step", 0.5)
    t_eval = np.linspace(t_settle, t_settle + t_avg, n_samples)
    sol = integrate(p, sp_const, (0.0, t_settle + t_avg), t_eval=t_eval,
                    max_step=max_step, **kwargs)
    a1 = sol.a[0]
    a2 = sol.a[1]
    n1 = np.abs(a1) ** 2
    n2 = np.abs(a2) ** 2
    s_minus = -sp_const + np.sqrt(p.gamma1_rad) * a1
    sp_sq = abs(complex(sp_const)) ** 2
    T_t = np.abs(s_minus) ** 2 / max(sp_sq, 1e-300)
    return {
        "n1_mean": float(np.mean(n1)),
        "n1_std":  float(np.std(n1)),
        "n2_mean": float(np.mean(n2)),
        "n2_std":  float(np.std(n2)),
        "T_mean":  float(np.mean(T_t)),
        "T_std":   float(np.std(T_t)),
        "a1_last": complex(a1[-1]),
        "is_limit_cycle": float(np.std(n1)) > 0.01 * max(float(np.mean(n1)), 1e-12),
    }
