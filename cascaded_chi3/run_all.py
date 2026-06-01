"""Regenerate every figure in figures/.

Run via:
    conda run -n ml_base python run_all.py
"""

from __future__ import annotations

import time
from importlib import import_module

MODULES = [
    "fig_validation",
    "fig_sign_flip",
    "fig_sum_rule",
    "fig_bistability",
    "fig_saturation",
    "fig_pulse_bandwidth",
    "fig_manley_rowe",
    "fig_bulk_vs_metasurface",
    "fig_transmission_lineshape",
    "fig_overlap_eta",
    "fig_n2_scaling",
    "fig_parametric_onset",
    "fig_cascaded_vs_native",
    "fig_phase_shift_per_pump",
    "fig_arch_response",
    "fig_loss_budget",
    "fig_phase_shift_per_pump_2d",
    "fig_fwm_efficiency",
    "fig_am_modulator",
    "fig_sensor_sensitivity",
    "fig_switching_energy",
    "fig_xpm",
    "fig_optical_limiter",
    "fig_saturable_absorber",
]


def main():
    for name in MODULES:
        t0 = time.time()
        print(f"=== {name} ===")
        mod = import_module(name)
        mod.main()
        print(f"    ({time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
