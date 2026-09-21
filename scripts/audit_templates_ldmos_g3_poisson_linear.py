#!/usr/bin/env python3
"""Audit WP1.5 Poisson/Jacobian scaling at the two frozen G3 failures."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


FOCUS_NODE = 4601
STAGING_RELATIVE = Path(
    "reference_staging/templates_ldmos_sentaurus2022/"
    "phase01_original_20260826_02/stage1_v4/phase23_t2022_contract_v5_wp15"
)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def run_probe(runner: Path, config_path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    status: dict[str, Any] = {}
    for line in reversed(completed.stdout.splitlines()):
        if line.lstrip().startswith("{"):
            status = json.loads(line)
            break
    if completed.returncode != 0:
        raise RuntimeError(
            f"probe failed for {config_path} ({completed.returncode}):\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    status["returncode"] = completed.returncode
    return status


def archive_baseline_drain_state(
    runner: Path,
    base_config_path: Path,
    output: Path,
    archived_state: Path,
) -> None:
    reproduction = output / "baseline_drain_reproduction"
    reproduction.mkdir(parents=True, exist_ok=True)
    config = json.loads(base_config_path.read_text(encoding="utf-8"))
    config["solver"]["linear_equilibration"] = {"mode": "off"}
    config["output_csv"] = str((reproduction / "curve.csv").resolve())
    config["sweep"]["write_state_file"] = str(
        (reproduction / "final_state.csv").resolve()
    )
    history = config["sweep"]["diagnostics"]["newton_history"]
    history["csv_file"] = "newton_history.csv"
    history["attempts_csv_file"] = "newton_attempts.csv"
    history["iterations_csv_file"] = "newton_iterations.csv"
    history["rejected_state_directory"] = "rejected_states"
    config_path = reproduction / "config.json"
    write_json(config_path, config)
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 1:
        raise RuntimeError(
            "baseline drain reproduction must retain the registered failure; "
            f"got return code {completed.returncode}:\n{completed.stdout}\n"
            f"{completed.stderr}"
        )
    candidates = sorted(
        (reproduction / "rejected_states").glob("*0p002145_best.csv")
    )
    if len(candidates) != 1:
        raise RuntimeError(
            "baseline drain reproduction did not produce one frozen 0.002145 V best state"
        )
    archived_state.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(candidates[0], archived_state)


def read_focus_rows(path: Path, nodes: set[int]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            node = int(row["node_id"])
            if node not in nodes:
                continue
            result.append(
                {
                    "node_id": node,
                    "x": float(row["x"]),
                    "y": float(row["y"]),
                    "poisson_residual": float(row["poisson_residual"]),
                    "diagonal": float(row["diagonal"]),
                    "poisson_row_l2_norm": float(row["poisson_row_l2_norm"]),
                    "poisson_column_l2_norm": float(row["poisson_column_l2_norm"]),
                    "full_row_l2_norm": float(row["full_row_l2_norm"]),
                    "full_column_l2_norm": float(row["full_column_l2_norm"]),
                    "raw_delta_psi_V": float(row["raw_delta_psi_V"]),
                    "equilibrated_delta_psi_V": float(
                        row["equilibrated_delta_psi_V"]
                    ),
                }
            )
    return sorted(result, key=lambda row: row["node_id"])


def case_config(
    base: dict[str, Any],
    state: Path,
    output_csv: Path,
    contact_biases: dict[str, float],
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    config["simulation_type"] = "newton_poisson_linear_probe"
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output_csv.resolve())
    config["focus_node"] = FOCUS_NODE
    config.pop("sweep", None)
    for contact in config.get("contacts", []):
        name = str(contact.get("name", ""))
        if name in contact_biases:
            contact["bias"] = contact_biases[name]
    return config


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    repo = args.repo.resolve()
    runner = (args.runner or repo / "build/vela_example_runner.exe").resolve()
    staging = repo / STAGING_RELATIVE
    output = (
        args.output.resolve()
        if args.output
        else repo / "reference_staging/templates_ldmos_g3_node4601_audit_20260829"
    )
    output.mkdir(parents=True, exist_ok=True)

    inputs = output / "frozen_inputs"
    drain_input = inputs / "drain_vd0p002145_baseline_best.csv"
    seed_input = inputs / "idvg_vg0_vd0p1_baseline_best.csv"
    drain_config = staging / "g3_drain_prebias_sentaurus_path.json"
    if not drain_input.is_file():
        archive_baseline_drain_state(
            runner, drain_config, output, drain_input
        )
    if not seed_input.is_file():
        seed_input.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(
            staging / (
                "g3_idvg_seed_rejected_states/attempt_1_bias_0p000000_best.csv"
            ),
            seed_input,
        )

    cases = {
        "drain_vd0p002145_best": {
            "config": drain_config,
            "state": drain_input,
            "contact_biases": {"drain": 0.00214466666666667},
        },
        "idvg_vg0_vd0p1_best": {
            "config": staging / "g3_idvg_seed.json",
            "state": seed_input,
            "contact_biases": {"gate": 0.0, "drain": 0.1},
        },
    }

    report: dict[str, Any] = {
        "schema": "vela.templates_ldmos.g3_poisson_linear_audit.v1",
        "focus_node": FOCUS_NODE,
        "scope": "fixed-state linear diagnostic; no production solver change",
        "cases": {},
    }
    for name, paths in cases.items():
        base = json.loads(paths["config"].read_text(encoding="utf-8"))
        case_dir = output / name
        csv_path = case_dir / "poisson_linear.csv"
        config_path = case_dir / "probe.json"
        write_json(
            config_path,
            case_config(
                base, paths["state"], csv_path, paths["contact_biases"]
            ),
        )
        status = run_probe(runner, config_path)
        patch_nodes = {int(node) for node in status["focus_patch_nodes"]}
        report["cases"][name] = {
            "source_config": str(paths["config"].resolve()),
            "source_state": str(paths["state"].resolve()),
            "probe_config": str(config_path.resolve()),
            "probe_csv": str(csv_path.resolve()),
            "status": status,
            "focus_patch_rows": read_focus_rows(csv_path, patch_nodes),
        }

    write_json(output / "summary.json", report)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
