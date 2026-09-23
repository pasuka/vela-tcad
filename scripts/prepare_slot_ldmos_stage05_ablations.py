#!/usr/bin/env python3
"""Prepare deterministic single-point Slot-LDMOS Stage-05 ablations."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = (
    REPO
    / "build-release/reference_tcad/slot_ldmos_sentaurus2022/run01"
    / "vela_ready_poly_nontransport"
)
TARGET_INNER_V = 0.018374398259206587
CASES = {
    "residual_off": {
        "coupling_mode": "postprocess_only",
        "source_jacobian": "finite_difference",
        "gummel_max_iter": 50,
    },
    "jacobian_finite_difference": {
        "coupling_mode": "self_consistent",
        "source_jacobian": "finite_difference",
        "gummel_max_iter": 50,
    },
    "jacobian_frozen": {
        "coupling_mode": "self_consistent",
        "source_jacobian": "frozen",
        "gummel_max_iter": 50,
    },
    "handoff_direct_newton": {
        "coupling_mode": "self_consistent",
        "source_jacobian": "finite_difference",
        "gummel_max_iter": 0,
    },
}


def build_case(base: dict[str, Any], name: str) -> dict[str, Any]:
    if name not in CASES:
        raise ValueError(f"unknown Stage-05 ablation {name!r}")
    settings = CASES[name]
    document = copy.deepcopy(base)
    root = f"outputs/ablations/stage05_first_probe/{name}"
    document["_comment"] = (
        "Stage-05 first-probe ablation at fixed inner drain voltage "
        f"{TARGET_INNER_V:.17g} V; no external-circuit root finder."
    )
    document["output_csv"] = f"{root}/iv.csv"

    solver = document["solver"]
    impact = solver["impact_ionization"]
    impact["coupling_mode"] = settings["coupling_mode"]
    impact["source_jacobian"] = settings["source_jacobian"]
    solver["handoff"]["gummel_max_iter"] = settings["gummel_max_iter"]

    sweep = document["sweep"]
    sweep.pop("external_circuit", None)
    sweep.pop("boundary_control", None)
    sweep["start"] = TARGET_INNER_V
    sweep["stop"] = TARGET_INNER_V
    sweep["bias_points"] = [TARGET_INNER_V]
    sweep["write_state_file"] = f"{root}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{root}/states/state"
    sweep["diagnostics"] = {
        "release_bv_config_audit": {
            "enabled": True,
            "csv_file": f"{root}/avalanche_summary.csv",
            "summary_file": f"{root}/avalanche_summary.md",
        },
        "newton_history": {
            "enabled": True,
            "csv_file": f"{root}/newton_history.csv",
            "attempts_csv_file": f"{root}/newton_attempts.csv",
            "iterations_csv_file": f"{root}/newton_iterations.csv",
        },
    }
    document["_stage05_ablation"] = {
        "case": name,
        "fixed_inner_voltage_V": TARGET_INNER_V,
        **settings,
    }
    return document


def prepare(bundle: Path) -> dict[str, Any]:
    bundle = bundle.resolve()
    base_path = bundle / "simulation_05_avalanche_on_60v.json"
    base = json.loads(base_path.read_text(encoding="utf-8"))
    manifest_cases: list[dict[str, Any]] = []
    for name in CASES:
        document = build_case(base, name)
        filename = f"simulation_ablation_stage05_{name}.json"
        (bundle / filename).write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        output_root = (
            bundle / "outputs" / "ablations" / "stage05_first_probe" / name
        )
        (output_root / "states").mkdir(parents=True, exist_ok=True)
        manifest_cases.append(
            {
                "name": name,
                "config": filename,
                **CASES[name],
            }
        )
    manifest = {
        "schema": "vela.slot_ldmos.stage05_first_probe_ablations.v1",
        "base_config": base_path.name,
        "fixed_inner_voltage_V": TARGET_INNER_V,
        "cases": manifest_cases,
        "comparisons": {
            "residual": ["residual_off", "jacobian_finite_difference"],
            "jacobian": ["jacobian_frozen", "jacobian_finite_difference"],
            "handoff": ["handoff_direct_newton", "jacobian_finite_difference"],
        },
    }
    (bundle / "stage05_first_probe_ablations_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    args = parser.parse_args()
    print(json.dumps(prepare(args.bundle), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
