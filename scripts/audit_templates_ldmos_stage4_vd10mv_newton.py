#!/usr/bin/env python3
"""Audit the failed Templates/LDMOS Vg=4 V, Vd=10 mV Newton transfer.

The audit replays the frozen initial, best and final states with the exact D5
physics contract.  It records residual/term/row diagnostics and checks analytic
Jacobian-vector products on each potential DOF in the one-ring patch around the
frozen electron/Poisson hot spots.  This is intentionally diagnostic-only: it
does not change the production sweep or enable predictor/IALMob.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


STATE_LABELS = ("initial", "best", "final")
FROZEN_HOTSPOTS = (675, 682, 683, 687, 695, 699, 701, 702, 707, 713, 729, 738)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def status_from_stdout(stdout: str) -> dict[str, Any]:
    for line in reversed(stdout.splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                value = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(value, dict):
                return value
    raise RuntimeError(f"runner emitted no JSON status:\n{stdout}")


def run_probe(runner: Path, config_path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    config_path.with_suffix(".log").write_text(completed.stdout, encoding="utf-8")
    status = status_from_stdout(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(json.dumps({
            "return_code": completed.returncode,
            "config": str(config_path),
            "status": status,
        }, indent=2))
    return status


def set_contact_bias(config: dict[str, Any], name: str, bias: float) -> None:
    for contact in config["contacts"]:
        if contact["name"] == name:
            contact["bias"] = bias
            return
    raise ValueError(f"contact {name!r} not present")


def cells_from_mesh(mesh: dict[str, Any]) -> list[list[int]]:
    cells = mesh.get("cells", mesh.get("triangles"))
    if not isinstance(cells, list):
        raise ValueError("mesh JSON has no cells or triangles array")
    result: list[list[int]] = []
    for cell in cells:
        nodes = cell.get("node_ids", cell.get("nodes"))
        if isinstance(nodes, list):
            result.append([int(node) for node in nodes])
    return result


def contact_nodes(mesh: dict[str, Any]) -> set[int]:
    result: set[int] = set()
    for contact in mesh.get("contacts", []):
        nodes = contact.get("node_ids", contact.get("nodes", []))
        result.update(int(node) for node in nodes)
    return result


def one_ring(cells: list[list[int]], seeds: set[int]) -> set[int]:
    result = set(seeds)
    for nodes in cells:
        if seeds.intersection(nodes):
            result.update(nodes)
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def summarize_terms(path: Path) -> dict[str, Any]:
    rows = read_rows(path)
    electron = [float(row["electron_residual"]) for row in rows]
    ranked = sorted(rows, key=lambda row: abs(float(row["electron_residual"])), reverse=True)
    return {
        "electron_residual_l2": norm(electron),
        "electron_residual_max_abs": max(abs(value) for value in electron),
        "top_electron_rows": [{
            "node_id": int(row["node_id"]),
            "residual": float(row["electron_residual"]),
            "flux": float(row["electron_flux"]),
            "flux_abs_sum": float(row["electron_flux_abs_sum"]),
            "recombination": float(row["electron_recombination"]),
        } for row in ranked[:12]],
    }


def summarize_carrier_rows(path: Path) -> dict[str, Any]:
    rows = read_rows(path)
    ranked = sorted(rows, key=lambda row: abs(float(row["electron_residual"])), reverse=True)
    return {
        "top_electron_rows": [{
            "node_id": int(row["node_id"]),
            "residual": float(row["electron_residual"]),
            "diagonal": float(row["electron_diagonal"]),
            "row_abs_sum": float(row["electron_row_abs_sum"]),
            "row_l2_norm": float(row["electron_row_l2_norm"]),
            "raw_delta_phin_V": float(row["raw_delta_phin_V"]),
            "capped_delta_phin_V": float(row["capped_delta_phin_V"]),
        } for row in ranked[:12]],
    }


def summarize_jvp(path: Path, contacts: set[int]) -> dict[str, Any]:
    rows = read_rows(path)
    ranked = sorted(rows, key=lambda row: float(row["phin_relative_error"]), reverse=True)
    return {
        "max_relative_error": max(float(row["relative_error"]) for row in rows),
        "max_phin_relative_error": max(float(row["phin_relative_error"]) for row in rows),
        "directions_over_1e_4": sum(float(row["phin_relative_error"]) > 1.0e-4 for row in rows),
        "top_directions": [{
            "direction": row["direction"],
            "node_id": int(row["direction"].split("_")[-1]),
            "is_contact": int(row["direction"].split("_")[-1]) in contacts,
            "relative_error": float(row["relative_error"]),
            "phin_relative_error": float(row["phin_relative_error"]),
            "max_error_node": int(row["max_error_node"]),
            "analytic_at_max_error": float(row["analytic_at_max_error"]),
            "finite_difference_at_max_error": float(row["finite_difference_at_max_error"]),
        } for row in ranked[:16]],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--states-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    runner = args.runner.resolve()
    base_path = args.base_config.resolve()
    states_dir = args.states_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)

    base = load_json(base_path)
    mesh_path = Path(base["mesh_file"])
    mesh = load_json(mesh_path)
    cells = cells_from_mesh(mesh)
    contacts = contact_nodes(mesh)
    patch_nodes = sorted(one_ring(cells, set(FROZEN_HOTSPOTS)))
    free_patch_nodes = [node for node in patch_nodes if node not in contacts]
    write_json(output / "input_manifest.json", {
        "base_config": str(base_path),
        "base_config_sha256": sha256(base_path),
        "mesh_file": str(mesh_path),
        "mesh_sha256": sha256(mesh_path),
        "frozen_hotspots": FROZEN_HOTSPOTS,
        "one_ring_nodes": patch_nodes,
        "free_one_ring_nodes": free_patch_nodes,
        "contact_one_ring_nodes": sorted(set(patch_nodes).intersection(contacts)),
    })

    summaries: dict[str, Any] = {}
    for label in STATE_LABELS:
        state = states_dir / f"attempt_2_bias_0p010000_{label}.csv"
        if not state.is_file():
            raise FileNotFoundError(state)
        state_out = output / label
        state_out.mkdir()
        common = json.loads(json.dumps(base))
        common.pop("sweep", None)
        common["state_file"] = str(state)
        set_contact_bias(common, "drain", 0.01)

        cases: list[tuple[str, str, dict[str, Any]]] = []
        for name, simulation_type in (
            ("residual", "newton_residual_probe"),
            ("carrier_terms", "newton_carrier_term_probe"),
            ("carrier_rows", "newton_carrier_row_probe"),
            ("newton_step", "newton_step_probe"),
        ):
            config = json.loads(json.dumps(common))
            config["simulation_type"] = simulation_type
            config["output_csv"] = str(state_out / f"{name}.csv")
            if name == "carrier_terms":
                config["carrier_term_probe"] = {"solved_equation_terms": True}
            cases.append((name, simulation_type, config))

        jvp = json.loads(json.dumps(common))
        jvp["simulation_type"] = "newton_jvp_probe"
        jvp["output_csv"] = str(state_out / "jvp.csv")
        jvp["directions"] = [{
            "name": f"psi_node_{node}",
            "mode": "psi",
            "amplitude_V": 1.0e-7,
            "node_ids": [node],
            "node_index_base": 0,
            "exclude_contacts": False,
        } for node in patch_nodes]
        cases.append(("jvp", "newton_jvp_probe", jvp))

        statuses: dict[str, Any] = {}
        for name, _, config in cases:
            config_path = state_out / f"{name}.json"
            write_json(config_path, config)
            statuses[name] = run_probe(runner, config_path)

        summaries[label] = {
            "state_file": str(state),
            "state_sha256": sha256(state),
            "statuses": statuses,
            "terms": summarize_terms(state_out / "carrier_terms.csv"),
            "carrier_rows": summarize_carrier_rows(state_out / "carrier_rows.csv"),
            "jvp": summarize_jvp(state_out / "jvp.csv", contacts),
        }

    result = {
        "schema": "vela.templates_ldmos.stage4_vd10mv_newton_audit.v1",
        "contract": {
            "gate_bias_V": 4.0,
            "drain_bias_V": 0.01,
            "predictor_enabled": False,
            "ialmob_enabled": False,
            "state_sequence": ["initial_iteration_0", "best_iteration_3", "final_iteration_40"],
        },
        "patch": {
            "frozen_hotspots": FROZEN_HOTSPOTS,
            "one_ring_node_count": len(patch_nodes),
            "free_one_ring_node_count": len(free_patch_nodes),
        },
        "states": summaries,
    }
    write_json(output / "summary.json", result)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
