"""Default dimensionless parameters for the cascaded-chi3 TCMT model.

All quantities below are dimensionless. The unit of frequency is gamma2
(the SH energy decay rate, equal to the FWHM linewidth of the SH
resonance). The unit of mode amplitude is gamma2/|beta|, so that beta = 1
in these units. Time unit is 1/gamma2.

Physical-unit reference (LiNbO3 telecom-band metasurface, from Section 9
of the note) is exposed only for annotations.
"""

from dataclasses import dataclass, replace


@dataclass
class TCMTParams:
    """Dimensionless TCMT parameters.

    Domega      : SH detuning (2*omega1 - omega2) / gamma2.
    gamma1      : FM energy decay rate / gamma2 (typ. > 1).
    gamma1_rad  : radiative part of gamma1 (port coupling, <= gamma1).
    alpha3      : bare Kerr coefficient, alpha3 * gamma2 / |beta_ref|^2.
                  Default 0 isolates the cascaded contribution.
    delta1      : pump detuning from FM resonance, (omega_p - omega1) / gamma2.
    beta        : SH coupling strength, dimensionless (default 1.0).
                  Set to 0 to disable cascading -- the FM equation reduces
                  to a single resonance with only the bare Kerr alpha3.
                  This gives the "native-Kerr" architecture (A) under the
                  same TCMT framework as the cascaded architecture (B).
    """
    Domega: float = 2.0
    gamma1: float = 10.0
    gamma1_rad: float = 5.0
    alpha3: float = 0.0
    delta1: float = 0.0
    beta: float = 1.0


# ---- helper for parameter sweeps (returns a new dataclass) ----------
def with_(p: TCMTParams, **kw) -> TCMTParams:
    return replace(p, **kw)


# ---- physical-unit reference (annotations only) --------------------
LAMBDA1_REAL = 1.55e-6                 # m
OMEGA1_REAL  = 2 * 3.141592653589793 * 299_792_458.0 / LAMBDA1_REAL  # rad/s
Q2_REF       = 4_000
GAMMA2_REAL  = 2 * OMEGA1_REAL / Q2_REF              # rad/s
N1_REAL, N2_REAL = 2.21, 2.26                        # LiNbO3 (kept as legacy ref)
CHI2_REAL    = 5e-11                                 # m/V (=2 d33 in LiNbO3)
ETA_GEOM     = 0.3                                   # geometric overlap
CHI3_NATIVE  = 2e-21                                 # m^2/V^2, LiNbO3 native


# ---- Host-material database for cross-host comparisons -------------
# Values at telecom (lambda_1 = 1.55 um). |chi^(2)| = 2 d_eff. n is the
# refractive index at the FM band (used for both n_1 and n_2 here as an
# approximation -- the few-percent dispersion is captured in N1/N2 only
# for LiNbO3 above). Inversion-symmetric hosts (Si, Si3N4) have chi^(2)
# = 0 by selection rule; cascade is structurally unavailable and we
# tag them accordingly with cascade_available=False.
HOSTS = {
    "LiNbO3":  dict(chi2=5.0e-11, chi3=2.0e-21, n=2.21, color="C0",
                    cascade_available=True),
    "GaP":     dict(chi2=1.0e-10, chi3=1.0e-19, n=3.05, color="C2",
                    cascade_available=True),
    "AlGaAs":  dict(chi2=1.0e-10, chi3=1.0e-18, n=3.30, color="C1",
                    cascade_available=True),
    "AlN":     dict(chi2=1.0e-11, chi3=2.0e-21, n=2.10, color="C4",
                    cascade_available=True),
    "Si":      dict(chi2=0.0,     chi3=4.0e-18, n=3.50, color="C5",
                    cascade_available=False),
    "Si3N4":   dict(chi2=0.0,     chi3=2.5e-21, n=2.00, color="C6",
                    cascade_available=False),
    "GaAs":    dict(chi2=3.5e-10, chi3=1.4e-18, n=3.50, color="C3",
                    cascade_available=True),
}
