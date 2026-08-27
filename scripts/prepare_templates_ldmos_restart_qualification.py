#!/usr/bin/env python3
"""Generate Templates/LDMOS exact-mesh stage-1.5 qualification decks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence


def solver(method: str = "newton") -> dict[str, Any]:
    if method == "frozen_state":
        return {"method": method, "impact_ionization": "none"}
    return {
        "method": "newton",
        "max_iter": 100,
        "reltol": 1.0e-7,
        "abstol": 1.0e-12,
        "damping_psi": 0.35,
        "warm_start": True,
        "carrier_statistics": {"model": "fermi_dirac"},
        "bandgap_narrowing": {
            "model": "old_slotboom",
            "fermi_statistics_correction": True,
        },
        "quasi_fermi_update_limit_V": 0.1,
        "line_search_mode": "block_filter",
        "residual_filter_gamma": 1.0e-4,
        "residual_filter_envelope_factor": 2.0,
        "continuity_row_scaling": {
            "enabled": True,
            "flux_fraction": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale": 1.0e-18,
            "min_weight": 1.0e-12,
            "max_weight": 1.0e12,
        },
    }


def contacts(gate_bias_V: float, drain_bias_V: float) -> list[dict[str, Any]]:
    return [
        {
            "name": "gate", "type": "metal_gate", "bias": gate_bias_V,
            "flatband_voltage": 0.0,
        },
        {"name": "drain", "type": "ohmic", "bias": drain_bias_V},
        {"name": "source", "type": "ohmic", "bias": 0.0},
        {"name": "substrate", "type": "ohmic", "bias": 0.0},
    ]


def deck(name: str,
         initial_state: str,
         output_state: str,
         gate_bias_V: float,
         drain_bias_V: float,
         frozen: bool = False,
         terminal_balance: bool = False) -> dict[str, Any]:
    sweep: dict[str, Any] = {
        "mode": "iv",
        "contact": "drain",
        "current_contact": "drain",
        "start": drain_bias_V,
        "stop": drain_bias_V,
        "step": 1.0,
        "initial_state_file": initial_state,
        "write_state_file": output_state,
        "write_vtk": False,
    }
    if frozen:
        sweep["frozen_state_compute_current"] = False
    if terminal_balance:
        sweep["diagnostics"] = {"terminal_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": f"{name}_terminal_balance.csv",
        }}
    return {
        "_comment": (
            "Templates/LDMOS stage-1.5 solver qualification only. "
            "The zero-flatband metal-gate mapping is provisional; stage 2 owns "
            "the PolySi work-function control."
        ),
        "simulation_type": "dc_sweep",
        "mesh_file": "../vela_exact_topology/mesh.json",
        "node_doping_file": "../vela_exact_topology/doping.csv",
        "output_csv": f"{name}.csv",
        "scaling": {"mode": "unit_scaling"},
        "contacts": contacts(gate_bias_V, drain_bias_V),
        "solver": solver("frozen_state" if frozen else "newton"),
        "sweep": sweep,
    }


def prepare(stage1_dir: Path) -> dict[str, Path]:
    output = stage1_dir / "qualification"
    output.mkdir(parents=True, exist_ok=True)
    specifications = {
        "frozen_eq_0v": (
            "sentaurus_eq_0v_state.csv", "frozen_eq_0v_roundtrip.csv", 0.0, 0.0,
            True, False),
        "reclose_eq_0v": (
            "sentaurus_eq_0v_state.csv", "reclose_eq_0v_state.csv", 0.0, 0.0,
            False, False),
        "reclose_eq_0v_repeat": (
            "reclose_eq_0v_state.csv", "reclose_eq_0v_repeat_state.csv", 0.0, 0.0,
            False, False),
        "reclose_idvg_vg0_vd0p1": (
            "sentaurus_idvg_vg0_vd0p1_state.csv", "reclose_idvg_vg0_vd0p1_state.csv",
            0.0, 0.1, False, True),
        "reclose_idvg_vg0_vd0p1_repeat": (
            "reclose_idvg_vg0_vd0p1_state.csv",
            "reclose_idvg_vg0_vd0p1_repeat_state.csv",
            0.0, 0.1, False, True),
        "reclose_idvg_vg0_vd0p1_repeat2": (
            "reclose_idvg_vg0_vd0p1_repeat_state.csv",
            "reclose_idvg_vg0_vd0p1_repeat2_state.csv",
            0.0, 0.1, False, True),
        "reclose_idvd_vg4_vd0p1": (
            "sentaurus_idvd_vg4_vd0p1_state.csv", "reclose_idvd_vg4_vd0p1_state.csv",
            4.0, 0.1, False, True),
        "reclose_idvd_vg4_vd0p1_repeat": (
            "reclose_idvd_vg4_vd0p1_state.csv",
            "reclose_idvd_vg4_vd0p1_repeat_state.csv",
            4.0, 0.1, False, True),
    }
    written: dict[str, Path] = {}
    for name, args in specifications.items():
        path = output / f"{name}.json"
        path.write_text(json.dumps(deck(name, *args), indent=2) + "\n", encoding="utf-8")
        written[name] = path
    return written


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    written = prepare(args.stage1_dir)
    print(json.dumps({"written": [str(path) for path in written.values()]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
