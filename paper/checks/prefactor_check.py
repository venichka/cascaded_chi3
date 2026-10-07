"""
Independent numerical checks of the cascaded-chi3 prefactor (2026-10-06).

Units: c = eps0 = 1.  Boyd amplitude convention (real field = E e^{-iwt} + c.c.):
    P(2w)  = chi2 E1^2          (SHG)
    P(w)   = 2 chi2 E2 E1^*     (DFG)
    P(w)   = 3 chi3 |E1|^2 E1   (SPM)
Normal modes are real-frequency Bloch modes of a 1-D layered cell, with
N_n = int eps_n |f_n|^2 (eps0 = 1).  Claimed (corrected) results:
    beta     = w1 M / N2,                     M = int chi2 f2^* f1 f1
    alpha3   = (3 w1 / 2) chi3 int|f1|^4 / N1
    chi3_2R  = -(2 w1/3) |M|^2 / (N2 Dw int|f1|^4)
    chi3_1R  = -(2 w1/3) |Mt|^2 / (N2 Dw V_cell),  Mt = int chi2 f2^* e e
The note's boxed forms are 1/4 of chi3_2R and chi3_1R.

Test A  beta: exact driven SH field vs single-mode amplitude -beta/Dw.
Test B  full 2R chain: exact SH field -> exact DFG polarisation -> chi3_eff
        defined exactly as the note does (match the |a|^2 a coefficient,
        i.e. <f1, P_casc> = 3 chi3_eff <f1,|f1|^2 f1>), vs chi3_2R.
Test C  1R chain: prescribed uniform FM field, exact SH, cell-integrated DFG
        polarisation matched to 3 chi3 |E|^2 E V_cell, vs chi3_1R.
Test D  alpha3 coefficient: exact self-consistent nonlinear eigenfrequency
        shift of the FM mode with a native local Kerr chi3 (iterated linear
        eigenproblems, no projection formula), vs dw = -alpha3 |a|^2 with
        alpha3 = (3 w1/2) chi3 int|f1|^4 / N1.
Test E  uniform-medium limit vs the bulk bound-wave result.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.linalg import eigh
from scipy.optimize import brentq

L, Ng = 1.0, 1200
z = np.arange(Ng) * L / Ng
dz = L / Ng
layerA = z < 0.45 * L
epsA, epsB = 2.2**2, 1.45**2
chi = np.where(layerA, 1.0, 0.0)          # chi2 only in layer A


def lap_bloch(k):
    """-d^2/dz^2 with Bloch condition f(z+L) = e^{ikL} f(z) (dense, Hermitian)."""
    D = np.zeros((Ng, Ng), complex)
    i = np.arange(Ng)
    D[i, i] = 2.0 / dz**2
    D[i[:-1], i[:-1] + 1] = -1.0 / dz**2
    D[i[1:], i[1:] - 1] = -1.0 / dz**2
    ph = np.exp(1j * k * L)
    D[Ng - 1, 0] = -ph / dz**2
    D[0, Ng - 1] = -np.conj(ph) / dz**2
    return D


def modes(k, eps, nmax=8):
    w2, V = eigh(lap_bloch(k), np.diag(eps), subset_by_index=[0, nmax - 1])
    return np.sqrt(w2.clip(0)), V


def integ(f):
    return np.sum(f) * dz


def sh_eps(d):
    return np.where(layerA, epsA * (1 + d), epsB * (1 + 0.3 * d))


def solve_sh(k2, W, eps2, src):
    """Exact driven SH field: (-d^2 - W^2 eps2) E2 = W^2 src."""
    A = sp.csr_matrix(lap_bloch(k2) - W**2 * np.diag(eps2))
    return spl.spsolve(A, W**2 * src)


print("=" * 72)
# ---------------------------------------------------------------- 2R setup
k1 = 0.37 * np.pi / L
eps1 = np.where(layerA, epsA, epsB)
w_all, V1 = modes(k1, eps1)
w1, f1 = w_all[0], V1[:, 0]
N1 = integ(eps1 * abs(f1)**2)
W = 2 * w1
w0, _ = modes(2 * k1, sh_eps(0.0))
n_sh = int(np.argmin(abs(w0 - W)))


def mismatch(d, target):
    w, _ = modes(2 * k1, sh_eps(d))
    return (W - w[n_sh]) / W - target


print("Tests A/B  (doubly resonant, non-uniform layered cell)")
for target in (2e-2, 5e-3, 1e-3):
    d = brentq(mismatch, -0.3, 0.3, args=(target,), xtol=1e-12)
    eps2 = sh_eps(d)
    w, V = modes(2 * k1, eps2)
    w2, f2 = w[n_sh], V[:, n_sh]
    Dw = W - w2
    N2 = integ(eps2 * abs(f2)**2)
    M = integ(chi * np.conj(f2) * f1**2)
    I4 = integ(abs(f1)**4)
    E2 = solve_sh(2 * k1, W, eps2, chi * f1**2)
    # Test A
    a2_exact = integ(eps2 * np.conj(f2) * E2) / N2
    a2_claim = -(w1 * M / N2) / Dw
    # Test B: exact DFG polarisation, chi3_eff by matching the projected cubic term
    P_casc = 2 * chi * E2 * np.conj(f1)
    chi3_exact = integ(np.conj(f1) * P_casc) / (3 * I4)
    chi3_claim = -(2 * w1 / 3) * abs(M)**2 / (N2 * Dw * I4)
    chi3_note = chi3_claim / 4
    print(f"  Dw/W={target:6.0e}  A: a2_exact/a2_claim={abs(a2_exact/a2_claim):.5f}"
          f"   B: chi3_exact/claim={chi3_exact.real/chi3_claim:.5f}"
          f"  /note={chi3_exact.real/chi3_note:.4f}")

# ---------------------------------------------------------------- Test C (1R)
print("Test C  (singly resonant: prescribed uniform FM field, Bloch k = 0)")
wc, _ = modes(0.0, sh_eps(0.0))
n_c = 2   # a non-trivial SH band at k = 0 with non-zero overlap


def mm_c(d, w1c, target):
    w, _ = modes(0.0, sh_eps(d))
    return (2 * w1c - w[n_c]) / (2 * w1c) - target


w1c = 0.5 * wc[n_c] * 1.0       # pump such that 2 w1c sits on band n_c at d = 0
Eloc = 1.0
for target in (2e-2, 5e-3, 1e-3):
    d = brentq(mm_c, -0.3, 0.3, args=(w1c, target), xtol=1e-12)
    eps2 = sh_eps(d)
    w, V = modes(0.0, eps2)
    w2, f2 = w[n_c], V[:, n_c]
    Wc = 2 * w1c
    Dw = Wc - w2
    N2 = integ(eps2 * abs(f2)**2)
    Mt = integ(chi * np.conj(f2))                 # e = unit vector, uniform
    E2 = solve_sh(0.0, Wc, eps2, chi * Eloc**2)
    P_cell = integ(2 * chi * E2 * np.conj(Eloc))
    chi3_exact = P_cell / (3 * abs(Eloc)**2 * Eloc * L)
    chi3_claim = -(2 * w1c / 3) * abs(Mt)**2 / (N2 * Dw * L)
    print(f"  Dw/W={target:6.0e}  chi3_exact/claim={chi3_exact.real/chi3_claim:.5f}"
          f"  /note={chi3_exact.real/(chi3_claim/4):.4f}"
          f"   (|Mt|^2 = {abs(Mt)**2:.3e}, band {n_c})")

# ---------------------------------------------------------------- Test D (alpha3)
print("Test D  (native Kerr: exact self-consistent nonlinear eigenfrequency)")
chi3n = 1.0
K = lap_bloch(k1)
for amp in (0.05, 0.1, 0.2):
    # mode amplitude a: field E = a f1 with f1 normalised to N1 at the start;
    # iterate eps_nl = eps1 + 3 chi3 |E|^2 (Boyd SPM), keeping the projection
    # of E on f1 equal to a.
    E = amp * f1
    for it in range(40):
        epsnl = eps1 + 3 * chi3n * abs(E)**2 * layerA
        wn, Vn = eigh(K, np.diag(epsnl), subset_by_index=[0, 0])
        g = Vn[:, 0]
        c = integ(eps1 * np.conj(f1) * g) / N1
        E = amp * g / c
    dw_exact = np.sqrt(wn[0]) - w1
    alpha3 = (3 * w1 / 2) * chi3n * integ(abs(f1)**4 * layerA) / N1
    dw_claim = -alpha3 * amp**2
    print(f"  |a|={amp:.0e}: dw_exact/dw_claim = {dw_exact/dw_claim:.5f}"
          f"   (note's alpha3 = 3w1/8 would give {dw_exact/(dw_claim/4):.4f})")

# ---------------------------------------------------------------- Test E
print("Test E  (uniform-medium limit vs bulk bound-wave)")
n1, n2, c2 = 2.21, 2.26, 1.0
w1u = 1.0
w2u = 2 * w1u * n1 / n2
Dwu = 2 * w1u - w2u
f01, f02 = 1 / n1, 1 / n2                      # V = 1, N = 1
M = c2 * f02 * f01**2
I4 = f01**4
chi3_claim = -(2 * w1u / 3) * M**2 / (Dwu * I4)
chi3_exact = -2 * c2**2 / (3 * (n2**2 - n1**2))
print(f"  exact/claim = {chi3_exact/chi3_claim:.4f}   exact/note = {chi3_exact/(chi3_claim/4):.4f}"
      f"   (residual = first order in n2-n1)")
print("=" * 72)
