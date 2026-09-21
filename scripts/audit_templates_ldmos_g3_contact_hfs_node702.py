#!/usr/bin/env python3
"""Audit the G3 contact-HFS Newton row and Jacobian around node 702."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from copy import deepcopy
from pathlib import Path


FOCUS_NODE = 702
NEIGHBOUR_NODES = [659, 701, 702, 703, 5545, 5546]


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    ucrt = r"D:\msys64\ucrt64\bin"
    env["PATH"] = ucrt + os.pathsep + env.get("PATH", "")
    return env


def run_probe(runner: Path, config: dict, config_path: Path) -> dict:
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"probe failed ({config['simulation_type']}): {completed.stderr or completed.stdout}")
    return json.loads(completed.stdout)


def read_focus_row(path: Path) -> dict[str, str]:
    for row in csv.DictReader(path.open(newline="", encoding="utf-8")):
        if int(row["node_id"]) == FOCUS_NODE:
            return row
    raise ValueError(f"node {FOCUS_NODE} missing from {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--state-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gate-bias", type=float, default=2.0 / 3.0)
    parser.add_argument("--newton-iterations", type=Path)
    parser.add_argument(
        "--disable-high-field",
        action="store_true",
        help="replace the mobility block by the constant low-field control",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    base = deepcopy(json.loads(args.base_config.read_text(encoding="utf-8")))
    base.pop("sweep", None)
    base["state_file"] = str(args.state_file.resolve())
    if args.disable_high_field:
        base["solver"]["mobility"] = {"model": "constant"}
    for contact in base["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = args.gate_bias

    outputs: dict[str, dict] = {}
    for name, simulation_type in (
        ("carrier_row", "newton_carrier_row_probe"),
        ("carrier_terms", "newton_carrier_term_probe"),
        ("edge_mobility", "edge_mobility_probe"),
    ):
        config = deepcopy(base)
        config["simulation_type"] = simulation_type
        config["output_csv"] = str((args.output_dir / f"{name}.csv").resolve())
        outputs[name] = run_probe(
            args.runner, config, args.output_dir / f"{name}.json")

    jvp = deepcopy(base)
    jvp["simulation_type"] = "newton_jvp_probe"
    jvp["output_csv"] = str((args.output_dir / "jvp.csv").resolve())
    jvp["directions"] = [
        {
            "name": f"psi_node_{node}",
            "mode": "psi",
            "node_ids": [node],
            "exclude_contacts": False,
            "amplitude_V": 1.0e-6,
        }
        for node in NEIGHBOUR_NODES
    ] + [
        {
            "name": f"phin_node_{FOCUS_NODE}_{amplitude:.0e}",
            "mode": "phin",
            "node_ids": [FOCUS_NODE],
            "exclude_contacts": False,
            "amplitude_V": amplitude,
        }
        for amplitude in (
            1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7,
            1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12,
        )
    ]
    outputs["jvp"] = run_probe(
        args.runner, jvp, args.output_dir / "jvp.json")

    mesh = json.loads(Path(base["mesh_file"]).read_text(encoding="utf-8"))
    cells = [
        cell for cell in mesh["triangles"]
        if FOCUS_NODE in cell["node_ids"]
    ]
    contact_membership = {
        contact["name"]: sorted(
            set(contact.get("node_ids", [])) & set(NEIGHBOUR_NODES))
        for contact in mesh["contacts"]
        if set(contact.get("node_ids", [])) & set(NEIGHBOUR_NODES)
    }

    row = read_focus_row(args.output_dir / "carrier_row.csv")
    terms = read_focus_row(args.output_dir / "carrier_terms.csv")
    jvp_rows = list(csv.DictReader(
        (args.output_dir / "jvp.csv").open(newline="", encoding="utf-8")))

    block_filter = None
    if args.newton_iterations is not None:
        iteration_rows = [
            item for item in csv.DictReader(
                args.newton_iterations.open(newline="", encoding="utf-8"))
            if item["attempt_id"] == "5"
        ]
        block_filter = {
            "iterations": len(iteration_rows) - 1,
            "accepted_full_steps": sum(
                item["event"] == "accepted_iteration"
                and float(item["damping"]) == 1.0
                and item["line_search_attempts"] == "1"
                and item["line_search_accepted"] == "1"
                for item in iteration_rows
            ),
            "line_search_rejections": sum(
                item["line_search_accepted"] == "0" for item in iteration_rows),
        }

    summary = {
        "schema": "vela.templates_ldmos.g3_contact_hfs_node702_audit.v1",
        "focus_node": FOCUS_NODE,
        "focus_node_is_contact": any(
            FOCUS_NODE in contact.get("node_ids", [])
            for contact in mesh["contacts"]),
        "adjacent_cells": [
            {"cell_id": cell["id"], "node_ids": cell["node_ids"]}
            for cell in cells
        ],
        "contact_nodes_in_patch": contact_membership,
        "carrier_row": row,
        "carrier_terms": terms,
        "jvp": jvp_rows,
        "block_filter": block_filter,
        "probe_status": outputs,
    }
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "pass",
        "summary": str(summary_path.resolve()),
        "focus_electron_residual": float(row["electron_residual"]),
        "max_jvp_relative_error": max(float(item["relative_error"]) for item in jvp_rows),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
