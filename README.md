# Cascaded χ² : χ² → χ³ in a dielectric metasurface

Theory notes on a doubly-resonant dielectric metasurface — a qBIC at the fundamental ω₁ and a high-Q BIC (or under-coupled qBIC) near 2ω₁ — designed to synthesise an effective Kerr response on ω₁ by cascading two χ² processes (SHG followed by DFG) instead of relying on the host material's native χ³. Because the SH channel is engineered not to radiate, the second harmonic can be adiabatically eliminated, leaving a Kerr-like term on the fundamental alone.

The derivation starts from Boyd-convention nonlinear polarisations and Brillouin-normalised modes, goes through both a polarisation-tensor route (χ²·G·χ²) and a temporal coupled-mode route, and ends with a single closed-form expression for χ³_casc(ω₁). A Python TCMT simulation validates the analytical predictions at the same operating points, and a sister note works out what changes if only the SH resonance is present.
