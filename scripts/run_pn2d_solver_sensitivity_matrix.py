#!/usr/bin/env python3
"""Generate PN2D low-bias solver sensitivity configs."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]


CASES: list[dict[str, Any]] = [
    {"name": "baseline", "solver": {}},
    {"name": "newton_damping_0p5", "solver": {"damping_factor": 0.5}},
    {"name": "newton_max_update_0p05", "solver": {"max_update": 0.05}},
    {
        "name": "branch_guard_0p05",
        "continuation": {
            "branch_acceptance": {
                "psi_phin_jump": True,
                "max_psi_phin_jump_V": 0.05,
            },
        },
    },
    {
        "name": "bank_rose_like_damped",
        "solver": {
            "damping_factor": 0.5,
            "line_search": True,
            "max_update": 0.05,
        },
    },
    {
        "name": "qf_hard_limit_0p0259",
        "solver": {"quasi_fermi_update_limit_V": 0.0259},
    },
    {
        "name": "bgn_contact_ni_match_offset_0p0197",
        "solver": {
            "bandgap_narrowing": {
                "model": "old_slotboom",
                "offset_eV": 0.0197,
            },
        },
    },
    {
        "name": "sentaurus_split_ni_slotboom",
        "config": {
            "materials_file": str(
                REPO_ROOT
                / "reference_tcad"
                / "pn2d_sentaurus2018"
                / "source"
                / "pn2d_sentaurus2018_iv_materials.json"
            ),
        },
        "solver": {
            "bandgap_narrowing": "slotboom",
        },
    },
]


REVERSE_LOW_SWEEP = {
    "start": 0,
    "stop": -5,
    "step": -0.05,
    "write_vtk": True,
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--bias-window", choices=["reverse_low"], required=True)
    return parser.parse_args(argv)


def merge_nested(target: dict[str, Any], update: dict[str, Any]) -> None:
    for key, value in update.items():
        if (
            isinstance(value, dict)
            and isinstance(target.get(key), dict)
        ):
            merge_nested(target[key], value)
        else:
            target[key] = copy.deepcopy(value)


def resolve_base_input_paths(config: dict[str, Any], base_dir: Path) -> None:
    for key in ("mesh_file", "materials_file", "node_doping_file"):
        value = config.get(key)
        if not isinstance(value, str):
            continue
        path = Path(value)
        if not path.is_absolute():
            path = base_dir / path
        config[key] = str(path.resolve())


def config_for_case(base_config: dict[str, Any],
                    case: dict[str, Any],
                    case_dir: Path) -> dict[str, Any]:
    config = copy.deepcopy(base_config)
    config["output_csv"] = str(case_dir / "iv.csv")

    sweep = config.setdefault("sweep", {})
    if not isinstance(sweep, dict):
        raise TypeError("base config sweep must be a JSON object")
    sweep.update(REVERSE_LOW_SWEEP)
    sweep["vtk_prefix"] = str(case_dir / "vtk" / case["name"])

    if "config" in case:
        merge_nested(config, case["config"])
    if "solver" in case:
        solver = config.setdefault("solver", {})
        if not isinstance(solver, dict):
            raise TypeError("base config solver must be a JSON object")
        merge_nested(solver, case["solver"])
    if "continuation" in case:
        continuation = sweep.setdefault("continuation", {})
        if not isinstance(continuation, dict):
            raise TypeError("base config sweep.continuation must be a JSON object")
        merge_nested(continuation, case["continuation"])
    return config


def status_for_case(case: dict[str, Any], dry_run: bool) -> tuple[str, str | None]:
    if dry_run:
        return "generated", "dry-run config generation only"
    return "generated", "config generation only; execute the generated config with vela_example_runner"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    base_config_path = args.base_config.resolve()
    base_config = json.loads(base_config_path.read_text(encoding="utf-8"))
    resolve_base_input_paths(base_config, base_config_path.parent)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_cases: list[dict[str, Any]] = []
    for case in CASES:
        case_dir = out_dir / case["name"]
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "vtk").mkdir(parents=True, exist_ok=True)

        config = config_for_case(base_config, case, case_dir)
        config_path = case_dir / "simulation.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

        status, reason = status_for_case(case, args.dry_run)
        manifest_case: dict[str, Any] = {
            "name": case["name"],
            "config": str(config_path),
            "out_dir": str(case_dir),
            "status": status,
        }
        if reason is not None:
            manifest_case["reason"] = reason
        manifest_cases.append(manifest_case)

    manifest = {
        "bias_window": args.bias_window,
        "dry_run": args.dry_run,
        "base_config": str(base_config_path),
        "cases": manifest_cases,
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
