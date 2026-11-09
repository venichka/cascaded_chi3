"""Bistability threshold with cascaded loss versus the native Kerr response.

Steady state of Eq. (1) with a complex Kerr coefficient alpha = alpha' + i alpha''
(Supplement 1, Section 4.C):

    G |s|^2 = y [ (delta1 + alpha' y)^2 + (gamma1/2 + alpha'' y)^2 ],   y = |a1|^2.

For each R the native coefficient is alpha3 = |Re alpha_casc| / R with
alpha_casc = -1/(Domega + i/2) (beta = gamma2 = 1) at Domega = -2 (same sign as
alpha3).  The bistability threshold is the smallest input power at which the
response becomes bistable, minimised over the pump detuning.  Prints the
reduction of the threshold relative to the native one and compares it with 1+R.

Run:  conda run -n ml_base python bistability_penalty.py
"""

import numpy as np

G1, G1E = 1.0, 0.5          # the penalty does not depend on gamma1 (scaling)


def threshold(alpha):
    a1, a2 = alpha.real, alpha.imag
    g = G1 / 2
    rho = a2 / abs(a1)
    dmin = g / np.tan(np.radians(30) - np.arctan(rho))
    best = np.inf
    for d in -np.sign(a1) * np.linspace(dmin * (1 + 1e-7), dmin * 3, 6000):
        r = np.roots([3 * abs(alpha) ** 2, 4 * (a1 * d + a2 * g), d * d + g * g])
        n = r[np.abs(r.imag) < 1e-12].real
        if len(n) < 2:
            continue
        P = n * ((d + a1 * n) ** 2 + (g + a2 * n) ** 2) / G1E
        best = min(best, P.max())
    return best


def main():
    ac = -1.0 / (-2.0 + 0.5j)
    print(" R    P_nat/P_casc   1+R    (1+R)/(P_nat/P_casc)")
    for R in (5, 10, 20, 50, 100, 1e4):
        a3 = abs(ac.real) / R
        ratio = threshold(a3 + 0j) / threshold(a3 + ac)
        print(f"{R:6g}  {ratio:10.3f}  {1 + R:7.0f}   {(1 + R) / ratio:6.2f}")


if __name__ == "__main__":
    main()
