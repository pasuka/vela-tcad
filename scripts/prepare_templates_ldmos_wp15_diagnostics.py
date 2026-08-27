#!/usr/bin/env python3
"""Prepare exact-mesh WP1.5 Newton/Gummel qualification variants.

The generated decks differ only in solver globalization controls.  They keep
the Templates/LDMOS G4 state, mesh, contacts, and physics contract fixed so a
failed or successful reclose can be attributed to solver behavior.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def diagnostic_deck(base: dict[str, Any], output_dir: Path, name: str) -> dict[str, Any]:
    deck = copy.deepcopy(base)
    deck["_comment"] = (
        "Templates/LDMOS WP1.5 exact-mesh solver qualification; physics and "
        "initial state are held fixed."
    )
    deck["output_csv"] = str((output_dir / f"{name}.csv").resolve())
    solver = deck["solver"]
    solver["diagnostics"] = True
    sweep = deck["sweep"]
    sweep["write_state_file"] = str((output_dir / f"{name}_state.csv").resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((output_dir / f"{name}_newton_history.csv").resolve()),
        "attempts_csv_file": str((output_dir / f"{name}_newton_attempts.csv").resolve()),
        "iterations_csv_file": str((output_dir / f"{name}_newton_iterations.csv").resolve()),
        "rejected_state_directory": str((output_dir / f"{name}_rejected_states").resolve()),
    }
    diagnostics["terminal_balance"]["csv_file"] = str(
        (output_dir / f"{name}_terminal_balance.csv").resolve()
    )
    diagnostics["srh_balance"]["csv_file"] = str(
        (output_dir / f"{name}_srh_balance.csv").resolve()
    )
    return deck


def prepare(base_path: Path, output_dir: Path) -> dict[str, str]:
    base = read_json(base_path)
    variants: dict[str, dict[str, Any]] = {}

    variants["baseline_block_filter"] = diagnostic_deck(
        base, output_dir, "baseline_block_filter"
    )

    merit = diagnostic_deck(base, output_dir, "merit")
    merit["solver"]["line_search_mode"] = "merit"
    variants["merit"] = merit

    merit_repeat = diagnostic_deck(base, output_dir, "merit_repeat")
    merit_repeat["solver"]["line_search_mode"] = "merit"
    merit_repeat["sweep"]["initial_state_file"] = str(
        (output_dir / "merit_state.csv").resolve()
    )
    variants["merit_repeat"] = merit_repeat

    merit_drain_ramp = diagnostic_deck(base, output_dir, "merit_drain_0p1")
    merit_drain_ramp["solver"]["line_search_mode"] = "merit"
    merit_drain_ramp["sweep"].update(
        {
            "bias_points": [0.0, 0.1],
            "start": 0.0,
            "stop": 0.1,
            "step": 0.1,
            "initial_state_file": str((output_dir / "merit_repeat_state.csv").resolve()),
        }
    )
    variants["merit_drain_0p1"] = merit_drain_ramp

    no_qf_cap = diagnostic_deck(base, output_dir, "block_filter_no_qf_cap")
    no_qf_cap["solver"]["quasi_fermi_update_limit_V"] = 0.0
    variants["block_filter_no_qf_cap"] = no_qf_cap

    tight_qf_cap = diagnostic_deck(base, output_dir, "block_filter_qf_0p025")
    tight_qf_cap["solver"]["quasi_fermi_update_limit_V"] = 0.025
    variants["block_filter_qf_0p025"] = tight_qf_cap

    damped = diagnostic_deck(base, output_dir, "block_filter_damping_0p25")
    damped["solver"]["damping_factor"] = 0.25
    variants["block_filter_damping_0p25"] = damped

    uniform = diagnostic_deck(base, output_dir, "uniform_trust_region")
    uniform["solver"].update(
        {
            "quasi_fermi_update_limit_mode": "uniform_trust_region",
            "quasi_fermi_trust_region_shrink_factor": 0.5,
            "quasi_fermi_trust_region_min_multiplier": 0.0625,
        }
    )
    variants["uniform_trust_region"] = uniform

    contact_majority = diagnostic_deck(base, output_dir, "merit_contact_majority")
    contact_majority["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_majority"}
    )
    variants["merit_contact_majority"] = contact_majority

    contact_basin = diagnostic_deck(base, output_dir, "merit_contact_basin")
    contact_basin["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_basin"}
    )
    variants["merit_contact_basin"] = contact_basin

    contact_basin_last_bias_repeat = diagnostic_deck(
        base, output_dir, "merit_contact_basin_last_bias_repeat"
    )
    contact_basin_last_bias_repeat["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_basin"}
    )
    last_bias = float(base["sweep"]["bias_points"][-1])
    contact_basin_last_bias_repeat["sweep"].update(
        {
            "bias_points": [last_bias],
            "start": last_bias,
            "stop": last_bias,
            "initial_state_file": str(
                (output_dir / "merit_contact_basin_state.csv").resolve()
            ),
        }
    )
    variants["merit_contact_basin_last_bias_repeat"] = (
        contact_basin_last_bias_repeat
    )

    contact_basin_last_bias_guard = diagnostic_deck(
        base, output_dir, "merit_contact_basin_last_bias_guard"
    )
    contact_basin_last_bias_guard["solver"].update(
        {
            "line_search_mode": "merit",
            "quasi_fermi_reference": "contact_basin",
            "max_iter": 0,
        }
    )
    contact_basin_last_bias_guard["sweep"].update(
        {
            "bias_points": [last_bias],
            "start": last_bias,
            "stop": last_bias,
            "initial_state_file": str(
                (output_dir / "merit_contact_basin_last_bias_repeat_state.csv").resolve()
            ),
        }
    )
    variants["merit_contact_basin_last_bias_guard"] = contact_basin_last_bias_guard

    contact_basin_unscaled = diagnostic_deck(
        base, output_dir, "merit_contact_basin_no_row_scaling"
    )
    contact_basin_unscaled["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_basin"}
    )
    contact_basin_unscaled["solver"].setdefault("continuity_row_scaling", {})[
        "enabled"
    ] = False
    variants["merit_contact_basin_no_row_scaling"] = contact_basin_unscaled

    contact_basin_unscaled_repeat = diagnostic_deck(
        base, output_dir, "merit_contact_basin_no_row_scaling_repeat"
    )
    contact_basin_unscaled_repeat["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_basin"}
    )
    contact_basin_unscaled_repeat["solver"].setdefault(
        "continuity_row_scaling", {}
    )["enabled"] = False
    contact_basin_unscaled_repeat["sweep"].update(
        {
            "bias_points": [last_bias],
            "start": last_bias,
            "stop": last_bias,
            "initial_state_file": str(
                (output_dir / "merit_contact_basin_no_row_scaling_state.csv").resolve()
            ),
        }
    )
    variants["merit_contact_basin_no_row_scaling_repeat"] = (
        contact_basin_unscaled_repeat
    )

    contact_basin_repeat = diagnostic_deck(
        base, output_dir, "merit_contact_basin_repeat"
    )
    contact_basin_repeat["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_basin"}
    )
    contact_basin_repeat["sweep"]["initial_state_file"] = str(
        (output_dir / "merit_contact_basin_state.csv").resolve()
    )
    variants["merit_contact_basin_repeat"] = contact_basin_repeat

    contact_basin_reclose_guard = diagnostic_deck(
        base, output_dir, "merit_contact_basin_reclose_guard"
    )
    contact_basin_reclose_guard["solver"].update(
        {
            "line_search_mode": "merit",
            "quasi_fermi_reference": "contact_basin",
            "max_iter": 0,
        }
    )
    contact_basin_reclose_guard["sweep"]["initial_state_file"] = str(
        (output_dir / "merit_contact_basin_repeat_state.csv").resolve()
    )
    variants["merit_contact_basin_reclose_guard"] = contact_basin_reclose_guard

    contact_basin_drain = diagnostic_deck(
        base, output_dir, "merit_contact_basin_drain_0p1"
    )
    contact_basin_drain["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_reference": "contact_basin"}
    )
    contact_basin_drain["sweep"].update(
        {
            "bias_points": [0.0, 0.1],
            "start": 0.0,
            "stop": 0.1,
            "step": 0.1,
            "initial_state_file": str(
                (output_dir / "merit_contact_basin_repeat_state.csv").resolve()
            ),
        }
    )
    variants["merit_contact_basin_drain_0p1"] = contact_basin_drain

    merit_tight_qf = diagnostic_deck(base, output_dir, "merit_qf_0p025")
    merit_tight_qf["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_update_limit_V": 0.025}
    )
    variants["merit_qf_0p025"] = merit_tight_qf

    merit_qf_near_frozen = diagnostic_deck(base, output_dir, "merit_qf_1e_6")
    merit_qf_near_frozen["solver"].update(
        {"line_search_mode": "merit", "quasi_fermi_update_limit_V": 1.0e-6}
    )
    variants["merit_qf_1e_6"] = merit_qf_near_frozen

    qualified_floor = diagnostic_deck(base, output_dir, "merit_stall_1e_6")
    qualified_floor["solver"].update(
        {"line_search_mode": "merit", "stall_residual_floor": 1.0e-6}
    )
    variants["merit_stall_1e_6"] = qualified_floor

    qualified_floor_repeat = diagnostic_deck(
        base, output_dir, "merit_stall_1e_6_repeat"
    )
    qualified_floor_repeat["solver"].update(
        {"line_search_mode": "merit", "stall_residual_floor": 1.0e-6}
    )
    qualified_floor_repeat["sweep"]["initial_state_file"] = str(
        (output_dir / "merit_stall_1e_6_state.csv").resolve()
    )
    variants["merit_stall_1e_6_repeat"] = qualified_floor_repeat

    poisson_only = diagnostic_deck(base, output_dir, "poisson_only")
    poisson_only["solver"]["method"] = "poisson_only"
    variants["poisson_only"] = poisson_only

    poisson_handoff = diagnostic_deck(base, output_dir, "poisson_handoff_newton")
    poisson_handoff["solver"].update(
        {"method": "newton", "line_search_mode": "merit"}
    )
    poisson_handoff["sweep"]["initial_state_file"] = str(
        (output_dir / "poisson_only_state.csv").resolve()
    )
    variants["poisson_handoff_newton"] = poisson_handoff

    poisson_handoff_repeat = diagnostic_deck(
        base, output_dir, "poisson_handoff_newton_repeat"
    )
    poisson_handoff_repeat["solver"].update(
        {"method": "newton", "line_search_mode": "merit"}
    )
    poisson_handoff_repeat["sweep"]["initial_state_file"] = str(
        (output_dir / "poisson_handoff_newton_state.csv").resolve()
    )
    variants["poisson_handoff_newton_repeat"] = poisson_handoff_repeat

    poisson_handoff_frozen_qf = diagnostic_deck(
        base, output_dir, "poisson_handoff_newton_qf_1e_12"
    )
    poisson_handoff_frozen_qf["solver"].update(
        {
            "method": "newton",
            "line_search_mode": "merit",
            "quasi_fermi_update_limit_V": 1.0e-12,
        }
    )
    poisson_handoff_frozen_qf["sweep"]["initial_state_file"] = str(
        (output_dir / "poisson_only_state.csv").resolve()
    )
    variants["poisson_handoff_newton_qf_1e_12"] = poisson_handoff_frozen_qf

    poisson_handoff_frozen_qf_repeat = diagnostic_deck(
        base, output_dir, "poisson_handoff_newton_qf_1e_12_repeat"
    )
    poisson_handoff_frozen_qf_repeat["solver"].update(
        {
            "method": "newton",
            "line_search_mode": "merit",
            "quasi_fermi_update_limit_V": 1.0e-12,
        }
    )
    poisson_handoff_frozen_qf_repeat["sweep"]["initial_state_file"] = str(
        (output_dir / "poisson_handoff_newton_qf_1e_12_state.csv").resolve()
    )
    variants["poisson_handoff_newton_qf_1e_12_repeat"] = (
        poisson_handoff_frozen_qf_repeat
    )

    poisson_handoff_frozen_qf_guard = diagnostic_deck(
        base, output_dir, "poisson_handoff_newton_qf_1e_12_guard"
    )
    poisson_handoff_frozen_qf_guard["solver"].update(
        {
            "method": "newton",
            "line_search_mode": "merit",
            "quasi_fermi_update_limit_V": 1.0e-12,
            "max_iter": 0,
        }
    )
    poisson_handoff_frozen_qf_guard["sweep"]["initial_state_file"] = str(
        (output_dir / "poisson_handoff_newton_qf_1e_12_repeat_state.csv").resolve()
    )
    variants["poisson_handoff_newton_qf_1e_12_guard"] = (
        poisson_handoff_frozen_qf_guard
    )

    poisson_handoff_drain = diagnostic_deck(
        base, output_dir, "poisson_handoff_drain_0p1"
    )
    poisson_handoff_drain["solver"].update(
        {"method": "newton", "line_search_mode": "merit"}
    )
    poisson_handoff_drain["sweep"].update(
        {
            "bias_points": [0.0, 0.1],
            "start": 0.0,
            "stop": 0.1,
            "step": 0.1,
            "initial_state_file": str(
                (output_dir / "poisson_handoff_newton_repeat_state.csv").resolve()
            ),
        }
    )
    variants["poisson_handoff_drain_0p1"] = poisson_handoff_drain

    gummel = diagnostic_deck(base, output_dir, "gummel")
    gummel["solver"]["method"] = "gummel"
    variants["gummel"] = gummel

    hybrid = diagnostic_deck(base, output_dir, "gummel_newton")
    hybrid["solver"].update(
        {
            "method": "gummel_newton",
            "handoff": {
                "gummel_max_iter": 25,
                "newton_max_iter": 100,
                "require_gummel_convergence": False,
                "fallback": "off",
            },
        }
    )
    variants["gummel_newton"] = hybrid

    paths: dict[str, str] = {}
    for name, deck in variants.items():
        path = output_dir / f"{name}.json"
        write_json(path, deck)
        paths[name] = str(path.resolve())

    write_json(
        output_dir / "manifest.json",
        {
            "schema": "vela.templates_ldmos.wp15_solver_matrix.v1",
            "base_config": str(base_path.resolve()),
            "invariants": ["mesh", "doping", "materials", "contacts", "initial_state"],
            "variants": paths,
        },
    )
    return paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    paths = prepare(args.base_config, args.output_dir)
    print(json.dumps(paths, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
