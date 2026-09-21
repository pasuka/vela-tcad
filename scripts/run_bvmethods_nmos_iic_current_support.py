#!/usr/bin/env python3
"""Replay accepted BVmethods NMOS states with one IIC current support."""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_BASE = (
    RUN_ROOT
    / "vela_validation/btbt_e2_iic_operator_controls_20260804/"
      "eparallel_edge_scalar/postprocess_only/simulation.json"
)
DEFAULT_INITIAL = (
    RUN_ROOT
    / "vela_validation/btbt_e2_adaptive_0_7_20260804/segment_6p0_6p5/"
      "states/accepted_state_bias_6p400000.csv"
)
DEFAULT_OUTPUT = (
    RUN_ROOT / "vela_validation/btbt_e2_iic_cell_vector_20260805"
)


def absolute(path: Path) -> str:
    return str(path.resolve())


def configure_outputs(config: dict[str, Any], output: Path) -> None:
    config["output_csv"] = absolute(output / "sweep.csv")
    sweep = config["sweep"]
    sweep["write_state_file"] = absolute(output / "last_state.csv")
    sweep["write_state_every_point_prefix"] = absolute(
        output / "states" / "accepted_state"
    )
    sweep["vtk_prefix"] = absolute(output / "vtk" / "state")
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["qf_bounds"]["csv_file"] = absolute(output / "qf_bounds.csv")
    diagnostics["sg_avalanche_edges"] = {
        "enabled": True,
        "csv_file": absolute(output / "sg_avalanche_edges.csv"),
    }
    diagnostics["terminal_current_method_compare"] = {
        "enabled": True,
        "contacts": ["drain", "source", "substrate"],
        "csv_file": absolute(output / "terminal_current_method_compare.csv"),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": ["drain", "source", "substrate"],
        "csv_file": absolute(output / "continuity_balance.csv"),
    }
    carrier = config["solver"]["carrier_row_convergence"]
    carrier["diagnostic_csv"] = absolute(output / "carrier_row_convergence.csv")
    carrier["trace_csv"] = absolute(output / "carrier_row_trace.csv")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--initial-state", type=Path, default=DEFAULT_INITIAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--biases", nargs="+", type=float, default=[6.4])
    parser.add_argument(
        "--current-approximation",
        default="nodal_vector_current_reconstructed",
    )
    parser.add_argument("--current-magnitude-mode", default="edge_scalar_abs")
    parser.add_argument(
        "--eparallel-field-recovery",
        choices=("edge_adjacent_cells", "nodal_vertex_star"),
    )
    parser.add_argument(
        "--source-volume-policy",
        choices=("genius_truncated", "genius_conservative", "edge_half_box", "edge_box"),
    )
    parser.add_argument(
        "--source-mapping-mode",
        choices=(
            "node_F_node_alpha_node_G",
            "edge_F_edge_alpha_edge_G_to_node",
            "cell_F_cell_alpha_cell_G_to_node",
            "nodal_eparallel_p1",
        ),
    )
    parser.add_argument(
        "--btbt-source-integration",
        choices=("semiconductor_cell_lumped", "transport_node_lumped"),
    )
    parser.add_argument("--taun", type=float)
    parser.add_argument("--taup", type=float)
    parser.add_argument(
        "--node-volume-policy",
        choices=("barycentric", "mixed_voronoi"),
    )
    parser.add_argument(
        "--high-field-gradient-discretization",
        choices=("edge_projection", "transport_cell_vector"),
    )
    parser.add_argument("--runner", type=Path)
    args = parser.parse_args()

    if args.runner is None:
        executable = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
        args.runner = REPO / "build-release" / executable

    config = copy.deepcopy(json.loads(args.base.read_text(encoding="utf-8-sig")))
    config["materials_file"] = absolute(
        REPO / "reference_tcad/bvmethods_sentaurus2018/vela/materials_sentaurus2018.json"
    )
    config["solver"]["bandgap_narrowing"] = {
        "model": "old_slotboom",
        "fermi_statistics_correction": True,
    }
    impact = config["solver"]["impact_ionization"]
    impact.update({
        "model": "van_overstraeten",
        "coupling_mode": "postprocess_only",
        "driving_force": "eparallel",
        "generation": "current_density",
        "current_approximation": args.current_approximation,
        "current_magnitude_mode": args.current_magnitude_mode,
    })
    if args.eparallel_field_recovery is not None:
        impact["eparallel_field_recovery"] = args.eparallel_field_recovery
    if args.source_volume_policy is not None:
        impact["source_volume_policy"] = args.source_volume_policy
    if args.source_mapping_mode is not None:
        impact["source_mapping_mode"] = args.source_mapping_mode
    if args.btbt_source_integration is not None:
        config["solver"].setdefault("band_to_band", {})[
            "source_integration"
        ] = args.btbt_source_integration
    if args.taun is not None:
        config["solver"]["taun"] = args.taun
    if args.taup is not None:
        config["solver"]["taup"] = args.taup
    if args.node_volume_policy is not None:
        config.setdefault("mesh_geometry", {})[
            "node_volume_policy"
        ] = args.node_volume_policy
    if args.high_field_gradient_discretization is not None:
        config["solver"].setdefault("mobility", {})[
            "high_field_gradient_discretization"
        ] = args.high_field_gradient_discretization
    sweep = config["sweep"]
    sweep.update({
        "start": args.biases[0],
        "stop": args.biases[-1],
        "bias_points": args.biases,
        "initial_state_file": absolute(args.initial_state),
        "stop_on_failure": True,
    })
    args.output.mkdir(parents=True, exist_ok=True)
    configure_outputs(config, args.output)
    config["_validation_case"] = {
        "purpose": "BVmethods NMOS current-type IIC current-support closure",
        "requested_biases_V": args.biases,
        "sentaurus_reference_BV_V": 6.377494277837012,
    }
    config_path = args.output / "simulation.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    completed = subprocess.run(
        [str(args.runner.resolve()), "--config", str(config_path.resolve())],
        cwd=args.output,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (args.output / "run.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        print("\n".join(completed.stdout.splitlines()[-40:]))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
