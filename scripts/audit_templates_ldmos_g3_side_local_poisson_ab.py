#!/usr/bin/env python3
"""Run a default-off region-local AverageBox Poisson A/B at the G3 gm endpoints."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.audit_templates_ldmos_averagebox_node4492 import (
    coefficient_edge,
    edge_key,
    parse_debug_block,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    if os.name == "nt":
        env["PATH"] = os.pathsep.join(
            [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin", env.get("PATH", "")]
        )
    env["VELA_LINEAR_SOLVER"] = "sparselu"
    return env


def run_config(runner: Path, config: dict[str, Any], path: Path) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(path.resolve())],
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (path.parent / f"{path.stem}.stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (path.parent / f"{path.stem}.stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    try:
        result = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"{config['simulation_type']} produced no JSON: "
            f"{completed.stderr or completed.stdout}"
        ) from error
    result["return_code"] = completed.returncode
    return result


def material_permittivities(materials: dict[str, Any]) -> dict[str, float]:
    return {
        str(material["name"]).lower(): float(material["eps_r"])
        for material in materials["materials"]
    }


def region_averagebox_poisson_profile(
    mesh: dict[str, Any],
    coefficients: dict[int, list[float]],
    eps_r_by_material: dict[str, float],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    points = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    regions = {int(region["id"]): region for region in mesh["regions"]}
    weighted_couple_um: dict[tuple[int, int], float] = defaultdict(float)
    adjacent_eps: dict[tuple[int, int], list[float]] = defaultdict(list)
    missing: list[int] = []

    for cell in mesh["triangles"]:
        cell_id = int(cell["id"])
        if cell_id not in coefficients:
            missing.append(cell_id)
            continue
        nodes = [int(value) for value in cell["node_ids"]]
        if len(nodes) != 3 or len(coefficients[cell_id]) != 3:
            raise ValueError(f"cell {cell_id} is not a Tri3 debug record")
        region = regions[int(cell["region_id"])]
        material = str(region["material"]).lower()
        if material not in eps_r_by_material:
            raise ValueError(f"missing eps_r for material {material}")
        eps_r = eps_r_by_material[material]
        for coefficient_index in range(3):
            node0, node1, _ = coefficient_edge(nodes, coefficient_index)
            key = edge_key(node0, node1)
            length_um = math.dist(points[node0], points[node1])
            local_couple_um = coefficients[cell_id][coefficient_index] * length_um
            weighted_couple_um[key] += eps_r * local_couple_um
            adjacent_eps[key].append(eps_r)
    if missing:
        raise ValueError(f"debug oracle lacks mesh cells: {missing[:20]}")

    rows: list[dict[str, Any]] = []
    interface_edges = 0
    for key in sorted(weighted_couple_um):
        eps_values = adjacent_eps[key]
        average_eps = sum(eps_values) / len(eps_values)
        effective_couple_um = weighted_couple_um[key] / average_eps
        if len({round(value, 14) for value in eps_values}) > 1:
            interface_edges += 1
        rows.append({
            "node0": key[0],
            "node1": key[1],
            "couple_m": effective_couple_um * 1.0e-6,
        })
    return rows, {
        "records": len(rows),
        "interface_edges": interface_edges,
        "zero_couples": sum(float(row["couple_m"]) == 0.0 for row in rows),
        "definition": "sum(eps_r_cell*averagebox_local_couple)/arithmetic_mean(eps_r_edge_cells)",
    }


def state_delta(left: Path, right: Path) -> dict[str, float]:
    a = {int(row["node_id"]): row for row in read_csv(left)}
    b = {int(row["node_id"]): row for row in read_csv(right)}
    result: dict[str, float] = {}
    for name in ("psi", "phin", "phip"):
        values = [float(b[node][name]) - float(a[node][name]) for node in sorted(a)]
        result[f"{name}_l2_V"] = math.sqrt(sum(value * value for value in values))
        result[f"{name}_maximum_abs_V"] = max(map(abs, values), default=0.0)
    return result


def endpoint(
    runner: Path,
    baseline_config_path: Path,
    initial_state: Path,
    profile: Path,
    profile_records: int,
    output: Path,
    poisson_charge_volume_policy: str,
) -> dict[str, Any]:
    baseline = json.loads(baseline_config_path.read_text(encoding="utf-8"))
    baseline.setdefault("discretization", {})[
        "poisson_charge_volume_policy"
    ] = poisson_charge_volume_policy
    baseline["state_file"] = str(initial_state.resolve())
    baseline_state = output / "baseline_state.csv"
    baseline["output_state_file"] = str(baseline_state.resolve())
    baseline_result = run_config(runner, baseline, output / "baseline.json")

    candidate = deepcopy(baseline)
    candidate_state = output / "candidate_state.csv"
    candidate["output_state_file"] = str(candidate_state.resolve())
    geometry = candidate.setdefault("mesh_geometry", {})
    geometry["poisson_couple_profile"] = "templates_ldmos_region_averagebox"
    geometry["external_averagebox_poisson_couples_file"] = str(profile.resolve())
    geometry["external_averagebox_poisson_expected_edges"] = profile_records
    candidate_result = run_config(runner, candidate, output / "candidate.json")

    jvp_result: dict[str, Any] | None = None
    if candidate_result.get("converged") and candidate_state.exists():
        jvp_csv = output / "candidate_jvp.csv"
        jvp = deepcopy(candidate)
        jvp["simulation_type"] = "newton_jvp_probe"
        jvp["state_file"] = str(candidate_state.resolve())
        jvp["output_csv"] = str(jvp_csv.resolve())
        jvp.pop("output_state_file", None)
        jvp["directions"] = [
            {
                "name": "interface_psi_ring1",
                "mode": "psi",
                "node_ids": [3721, 3974, 3973, 4091, 3727, 4021, 3771],
                "adjacent_cell_rings": 1,
                "exclude_contacts": False,
                "amplitude_V": 1.0e-6,
            },
            {
                "name": "interface_psi_minus_phin_ring1",
                "mode": "psi_minus_phin",
                "node_ids": [3721, 3974, 3973, 4091, 3727, 4021, 3771],
                "adjacent_cell_rings": 1,
                "exclude_contacts": False,
                "amplitude_V": 1.0e-6,
            },
        ]
        jvp_result = run_config(runner, jvp, output / "candidate_jvp.json")
        jvp_rows = read_csv(jvp_csv)
        jvp_result["maximum_relative_error"] = max(
            float(row["relative_error"]) for row in jvp_rows
        )

    result: dict[str, Any] = {
        "bias_V": next(
            float(contact["bias"])
            for contact in baseline["contacts"] if contact["name"].lower() == "gate"
        ),
        "baseline": baseline_result,
        "candidate": candidate_result,
        "candidate_jacobian": jvp_result,
    }
    if baseline_state.exists() and candidate_state.exists():
        result["candidate_minus_baseline_state"] = state_delta(
            baseline_state, candidate_state
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--measure-coefficients", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, action="append", required=True)
    parser.add_argument("--initial-state", type=Path, action="append", required=True)
    parser.add_argument("--sentaurus-current", type=float, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--poisson-charge-volume-policy",
        choices=("global", "material_local"),
        default="global",
    )
    args = parser.parse_args()
    if not (
        len(args.baseline_config) == len(args.initial_state)
        == len(args.sentaurus_current) == 2
    ):
        raise ValueError("exactly two aligned endpoint inputs are required")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    first = json.loads(args.baseline_config[0].read_text(encoding="utf-8"))
    mesh = json.loads(Path(first["mesh_file"]).read_text(encoding="utf-8"))
    materials = json.loads(Path(first["materials_file"]).read_text(encoding="utf-8"))
    profile_rows, profile_metadata = region_averagebox_poisson_profile(
        mesh,
        parse_debug_block(args.measure_coefficients.resolve(), "Coefficients"),
        material_permittivities(materials),
    )
    profile = output / "poisson_couples.csv"
    write_csv(profile, profile_rows)

    endpoints = [
        endpoint(
            args.runner.resolve(), args.baseline_config[index],
            args.initial_state[index], profile, len(profile_rows),
            output / f"endpoint_{index}",
            args.poisson_charge_volume_policy,
        )
        for index in range(2)
    ]
    delta_v = endpoints[1]["bias_V"] - endpoints[0]["bias_V"]
    baseline_currents = [
        float(item["baseline"]["contact_currents_A_per_um"]["drain"])
        for item in endpoints
    ]
    candidate_currents = [
        float(item["candidate"]["contact_currents_A_per_um"]["drain"])
        for item in endpoints
    ]
    sentaurus_gm = (args.sentaurus_current[1] - args.sentaurus_current[0]) / delta_v
    baseline_gm = (baseline_currents[1] - baseline_currents[0]) / delta_v
    candidate_gm = (candidate_currents[1] - candidate_currents[0]) / delta_v
    summary = {
        "schema": "vela.templates_ldmos.g3_side_local_poisson_ab.v1",
        "contract": {
            "carrier_transport": "qualified external AverageBox",
            "baseline_poisson": "mesh default barycentric/cotangent-fallback",
            "candidate_poisson": "region-local material-weighted AverageBox couples",
            "poisson_charge_volume_policy":
                args.poisson_charge_volume_policy,
            "initial_state": "independent converged Vela endpoint",
            "ialmob": "disabled",
            "predictor": "disabled",
            "production_default_changed": False,
        },
        "profile": profile_metadata,
        "endpoints": endpoints,
        "gm_A_per_um_V": {
            "sentaurus": sentaurus_gm,
            "baseline": baseline_gm,
            "candidate": candidate_gm,
        },
        "gm_ratio_to_sentaurus": {
            "baseline": baseline_gm / sentaurus_gm,
            "candidate": candidate_gm / sentaurus_gm,
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "profile": profile_metadata,
        "converged": [item["candidate"]["converged"] for item in endpoints],
        "gm_ratio_to_sentaurus": summary["gm_ratio_to_sentaurus"],
        "summary": str((output / "summary.json").resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
