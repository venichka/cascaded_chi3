"""Bulk LiNbO3 cascade at 1.064 um versus the Z-scan of DeSalvo et al. (1996).

The bound-wave expression chi_casc = -2|chi2|^2 / [3(n2^2 - n1^2)]
(Supplement 1, Eq. bulk) for the e-polarised d33 interaction, converted to an
intensity index n2_I = 3 chi / (4 eps0 c n1^2), added to the native index
without the cascade, 3.9e-19 m^2/W (Phillips et al. 2011), should reproduce the
measured 0.91e-19 m^2/W.

Indices: Zelmon et al. (1997), congruent LiNbO3, extraordinary, 21 C.
Run:  conda run -n ml_base python bulk_desalvo_check.py
"""

import numpy as np

EPS0, C0 = 8.8541878128e-12, 2.99792458e8


def n_e_zelmon(lam_um):
    """Extraordinary index of congruent LiNbO3 (Zelmon 1997)."""
    l2 = lam_um ** 2
    A = (2.9804, 0.5981, 8.9543)
    B = (0.02047, 0.0666, 416.08)
    return np.sqrt(1 + sum(a * l2 / (l2 - b) for a, b in zip(A, B)))


def main():
    d33 = 25.2e-12                       # Shoji et al. 1997, 1.064 um
    chi2 = 2 * d33
    n1, n2 = n_e_zelmon(1.064), n_e_zelmon(0.532)
    chi_casc = -2 * chi2 ** 2 / (3 * (n2 ** 2 - n1 ** 2))
    n2I_casc = 3 * chi_casc / (4 * EPS0 * C0 * n1 ** 2)
    n2I_nat = 3.9e-19
    print(f"n_e(1.064) = {n1:.4f}, n_e(0.532) = {n2:.4f}")
    print(f"chi_casc = {chi_casc:.3e} m^2/V^2, n2_I(casc) = {n2I_casc:.3e} m^2/W")
    print(f"native + cascade = {n2I_nat + n2I_casc:.3e} m^2/W "
          f"(DeSalvo 1996: 0.91e-19 m^2/W)")


if __name__ == "__main__":
    main()
