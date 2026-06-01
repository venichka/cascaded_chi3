"""
Architecture-comparison simulation toolbox.

Builds on the TCMT integrator from cascaded_chi3_simulation.py. The
fundamental model is the same; the new content here is dedicated
infrastructure for:

  (a) Architecture comparison A vs A+B (native-Kerr vs cascade)
  (b) Application observables (SPM, bistability, FWM, sat absorber,
      AM modulator, sensor) with quantitative R-scaling

All units are dimensionless: gamma_2 = 1 sets the time scale. The
"comparison" is internal between alpha_3_native and alpha_3_casc,
both expressed in the same units.

Boyd conventions (matching the LaTeX note):
    P^(2)(2w) = eps_0 chi^(2) E^2                       (SHG; degeneracy 1)
    P^(2)(w)_DFG = 2 eps_0 chi^(2) E(2w) E^*(w)         (DFG; degeneracy 2)
    P^(3)(w) = 3 eps_0 chi^(3) |E|^2 E                  (SPM; degeneracy 3)
Modal couplings:
    beta = (eps_0 w_1 / 4) * Integral[chi^(2)_ijk f_2*_i f_1_j f_1_k]
    alpha_3 = (3 eps_0 w_1 / 8) * Integral[chi^(3)_ijkl f_1*_i f_1_j f_1_k f_1*_l]
"""

import os
import sys
import numpy as np
from dataclasses import dataclass

# Pull the validated TCMT integrator from the previous simulation
sys.path.insert(0, '/home/claude/sim')
from cascaded_chi3_simulation import (
    tcmt_rhs,
    steady_state_cw,
    alpha3_casc_full,
    alpha3_loss_casc,
)

# ============================================================
#         MATERIAL PARAMETERS AND R-RATIO HELPERS
# ============================================================

@dataclass
class Material:
    """Bulk material parameters for cascade enhancement estimates.

    chi2: |chi^(2)| in m/V
    chi3_native: |chi^(3)_native| in m^2/V^2
    n_FM, n_SH: refractive indices at FM and SH wavelengths
    name: display name
    """
    name: str
    chi2: float           # m/V
    chi3_native: float    # m^2/V^2
    n_FM: float
    n_SH: float

    @property
    def n2sq(self):
        return self.n_SH ** 2


# Representative materials (numbers from the LaTeX note table)
MATERIALS = {
    'LiNbO3':   Material('LiNbO$_3$',  5e-11,  2e-21,    2.21, 2.26),
    'GaP':      Material('GaP',        1e-10,  1e-19,    3.05, 3.10),
    'AlGaAs':   Material('AlGaAs',     1e-10,  1e-18,    3.30, 3.50),
    'AlN':      Material('AlN',        4e-12,  2e-21,    2.05, 2.10),
    'Si':       Material('Si',         0.0,    4e-18,    3.48, 3.45),       # chi2=0
    'Si3N4':    Material('Si$_3$N$_4$',0.0,    2.5e-21,  2.00, 2.05),       # chi2=0
}


def cascade_ratio_R(mat, Q2, eta_geom=0.3, Domega_over_gamma2=2.0, omega1=1.22e15):
    """Return the cascaded enhancement ratio R = |chi^(3)_casc| / |chi^(3)_native|.

    Uses the order-of-magnitude form
        |chi^(3)_casc| ~ omega_1 / (6 n_2^2 |Domega + i gamma_2/2|) * |chi^(2)|^2 * eta_geom

    with Domega = (Domega/gamma_2) * gamma_2 and gamma_2 = omega_2 / Q_2 ~ 2 omega_1 / Q_2.
    """
    if mat.chi2 == 0:
        return np.nan  # No cascade possible

    gamma_2 = 2.0 * omega1 / Q2
    Domega = Domega_over_gamma2 * gamma_2
    denom = np.sqrt(Domega**2 + (gamma_2 / 2.0) ** 2)
    chi3_casc = omega1 * mat.chi2 ** 2 * eta_geom / (6.0 * mat.n2sq * denom)
    return chi3_casc / mat.chi3_native


def crossover_Q2(mat, eta_geom=0.3):
    """Q_2 at which |chi^(3)_casc| = |chi^(3)_native| with Domega=2 gamma_2.

    From R = 1: 1 = omega_1 / (6 n_2^2 sqrt((Domega)^2 + (gamma_2/2)^2)) * |chi^(2)|^2 eta / chi3_native.
    With Domega = 2 gamma_2 and gamma_2 = 2 omega_1 / Q_2:
        sqrt(4 + 0.25) gamma_2 = ... -> Q_2_crossover ~ 12 sqrt(4.25) n_2^2 chi3_n / (|chi^(2)|^2 eta)
    """
    if mat.chi2 == 0:
        return np.inf
    # From eq.: R = 1 => Q_2 = 12 * sqrt(4.25) * n_2^2 * chi3_n / (|chi^(2)|^2 * eta)
    return 12.0 * np.sqrt(4.25) * mat.n2sq * mat.chi3_native / (mat.chi2 ** 2 * eta_geom)


# ============================================================
#         OUTPUT DIRECTORY & STYLING
# ============================================================
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 11.5,
    'axes.labelsize': 11,
    'legend.fontsize': 9.5,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 110,
    'savefig.dpi': 160,
    'savefig.bbox': 'tight',
    'mathtext.fontset': 'cm',
})

OUT = '/home/claude/sim_comparison/figures'
os.makedirs(OUT, exist_ok=True)


def adjacent_save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path)
    print(f'  saved {path}')
    return path
