#!/usr/bin/env python3
"""Run the WP3 T3 HFS-gradient A/B on two frozen LDMOS G3 SSS states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any


FIXED_NODES = (3721, 3974, 3973, 4091, 3727, 4021, 3771)
VARIANTS = ("edge_projection", "transport_cell_vector")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty table")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    return env


def run_probe(runner: Path, config: dict[str, Any], config_path: Path) -> dict[str, Any]:
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (config_path.parent / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (config_path.parent / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"probe failed for {config_path}: {completed.stderr or completed.stdout}"
        )
    lines = completed.stdout.strip().splitlines()
    return json.loads(lines[-1]) if lines else {"return_code": completed.returncode}


def mesh_node_sets(mesh_path: Path) -> dict[str, set[int]]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    material_by_region = {
        int(region["id"]): str(region.get("material", "")).lower()
        for region in mesh["regions"]
    }
    silicon_regions = {
        region_id
        for region_id, material in material_by_region.items()
        if material in {"si", "silicon"}
    }
    node_regions: dict[int, set[int]] = defaultdict(set)
    neighbors: dict[int, set[int]] = defaultdict(set)
    silicon_nodes: set[int] = set()
    for cell in mesh["triangles"]:
        region_id = int(cell["region_id"])
        nodes = [int(node) for node in cell["node_ids"]]
        for node in nodes:
            node_regions[node].add(region_id)
            if region_id in silicon_regions:
                silicon_nodes.add(node)
            neighbors[node].update(other for other in nodes if other != node)
    contact_nodes = {
        int(node)
        for contact in mesh["contacts"]
        for node in contact.get("node_ids", [])
    }
    direct_interface = {
        node
        for node in silicon_nodes
        if any(region not in silicon_regions for region in node_regions[node])
    }
    interface_band = set(direct_interface)
    for node in direct_interface:
        interface_band.update(neighbors[node] & silicon_nodes)
    free_silicon = silicon_nodes - contact_nodes
    return {
        "silicon": silicon_nodes,
        "contact": contact_nodes,
        "free_silicon": free_silicon,
        "direct_interface": direct_interface,
        "interface_band": interface_band,
        "free_direct_interface": direct_interface & free_silicon,
        "free_interface_band": interface_band & free_silicon,
    }


def residual_map(rows: list[dict[str, str]]) -> dict[int, float]:
    return {int(row["node_id"]): float(row["electron_residual"]) for row in rows}


def norm(values: list[float]) -> dict[str, float]:
    return {
        "l2": math.sqrt(sum(value * value for value in values)),
        "maximum_abs": max((abs(value) for value in values), default=0.0),
    }


def norm_on_nodes(residuals: dict[int, float], nodes: set[int] | tuple[int, ...]) -> dict[str, Any]:
    ordered = sorted(nodes)
    missing = [node for node in ordered if node not in residuals]
    if missing:
        raise ValueError(f"residual probe is missing nodes: {missing}")
    return {"node_count": len(ordered), **norm([residuals[node] for node in ordered])}


def ratio(candidate: float, baseline: float) -> float:
    if baseline == 0.0:
        return 0.0 if candidate == 0.0 else math.inf
    return candidate / baseline


def no_worse(candidate: float, baseline: float) -> bool:
    return candidate <= baseline * (1.0 + 1.0e-12) + 1.0e-30


def compare_metric(candidate: dict[str, Any], baseline: dict[str, Any]) -> dict[str, Any]:
    return {
        "baseline": baseline,
        "candidate": candidate,
        "ratio": {
            "l2": ratio(candidate["l2"], baseline["l2"]),
            "maximum_abs": ratio(candidate["maximum_abs"], baseline["maximum_abs"]),
        },
        "not_worse": (
            no_worse(candidate["l2"], baseline["l2"])
            and no_worse(candidate["maximum_abs"], baseline["maximum_abs"])
        ),
    }


def top_nodes(residuals: dict[int, float], eligible: set[int], count: int = 7) -> list[int]:
    return sorted(eligible, key=lambda node: (-abs(residuals[node]), node))[:count]


def analyze_endpoint(
    baseline_rows: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
    node_sets: dict[str, set[int]],
) -> dict[str, Any]:
    baseline = residual_map(baseline_rows)
    candidate = residual_map(candidate_rows)
    eligible = node_sets["silicon"]
    baseline_top = top_nodes(baseline, eligible)
    candidate_top = top_nodes(candidate, eligible)

    fixed_baseline = norm_on_nodes(baseline, FIXED_NODES)
    fixed_candidate = norm_on_nodes(candidate, FIXED_NODES)
    per_node = []
    for node in FIXED_NODES:
        baseline_abs = abs(baseline[node])
        candidate_abs = abs(candidate[node])
        per_node.append({
            "node_id": node,
            "baseline": baseline[node],
            "candidate": candidate[node],
            "absolute_ratio": ratio(candidate_abs, baseline_abs),
            "passes_half_gate": candidate_abs <= 0.5 * baseline_abs,
        })
    fixed_comparison = compare_metric(fixed_candidate, fixed_baseline)
    fixed_comparison["per_node"] = per_node
    fixed_comparison["passes_half_gate"] = (
        fixed_comparison["ratio"]["l2"] <= 0.5
        and fixed_comparison["ratio"]["maximum_abs"] <= 0.5
        and all(item["passes_half_gate"] for item in per_node)
    )

    candidate_hotspot_comparison = compare_metric(
        norm_on_nodes(candidate, set(candidate_top)),
        norm_on_nodes(baseline, set(candidate_top)),
    )
    interface_comparison = compare_metric(
        norm_on_nodes(candidate, node_sets["interface_band"]),
        norm_on_nodes(baseline, node_sets["interface_band"]),
    )
    silicon_comparison = compare_metric(
        norm_on_nodes(candidate, eligible),
        norm_on_nodes(baseline, eligible),
    )
    migration_pass = all(
        item["not_worse"]
        for item in (candidate_hotspot_comparison, interface_comparison, silicon_comparison)
    )
    return {
        "fixed_seven": fixed_comparison,
        "candidate_top_seven": {
            "node_ids": candidate_top,
            **candidate_hotspot_comparison,
        },
        "interface_band": interface_comparison,
        "all_silicon": silicon_comparison,
        "hotspot_migration": {
            "baseline_top_seven": baseline_top,
            "candidate_top_seven": candidate_top,
            "overlap": sorted(set(baseline_top) & set(candidate_top)),
            "removed": sorted(set(baseline_top) - set(candidate_top)),
            "added": sorted(set(candidate_top) - set(baseline_top)),
        },
        "passes_migration_guard": migration_pass,
        "passes_t3_gate": fixed_comparison["passes_half_gate"] and migration_pass,
    }


def parse_endpoint(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("endpoint must be LABEL=CONFIG")
    label, config = value.split("=", 1)
    return label, Path(config)


def validate_contract(config: dict[str, Any]) -> None:
    if config.get("simulation_type") != "newton_carrier_term_probe":
        raise ValueError("T3 requires the production newton_carrier_term_probe")
    if "predictor" in config.get("sweep", {}):
        raise ValueError("T3 requires predictor disabled")
    if "ialmob" in json.dumps(config["solver"].get("mobility", {})).lower():
        raise ValueError("T3 requires IALMob disabled")
    if config["solver"]["mobility"].get("high_field_driving_force") != "quasi_fermi_gradient":
        raise ValueError("T3 requires quasi-Fermi-gradient HFS driving force")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--endpoint", type=parse_endpoint, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.endpoint) != 2:
        raise ValueError("T3 requires exactly two endpoints")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs: dict[str, Any] = {}
    rows_by_endpoint: dict[str, dict[str, list[dict[str, str]]]] = {}
    statuses: dict[str, dict[str, Any]] = {}
    artifacts: dict[str, dict[str, Any]] = {}
    common_mesh: Path | None = None
    node_sets: dict[str, set[int]] | None = None

    for label, config_path in args.endpoint:
        base = json.loads(config_path.resolve().read_text(encoding="utf-8"))
        validate_contract(base)
        mesh_path = Path(base["mesh_file"]).resolve()
        if common_mesh is None:
            common_mesh = mesh_path
            node_sets = mesh_node_sets(mesh_path)
        elif mesh_path != common_mesh:
            raise ValueError("both T3 endpoints must use the same mesh")
        inputs[label] = {
            "base_config": str(config_path.resolve()),
            "state_file": str(Path(base["state_file"]).resolve()),
            "gate_bias_V": next(
                float(contact["bias"])
                for contact in base["contacts"]
                if str(contact["name"]).lower() == "gate"
            ),
            "prior_baseline_csv": str(Path(base["output_csv"]).resolve()),
        }
        rows_by_endpoint[label] = {}
        statuses[label] = {}
        artifacts[label] = {}
        for variant in VARIANTS:
            run_dir = output / "probes" / label / variant
            run_dir.mkdir(parents=True, exist_ok=True)
            config = deepcopy(base)
            config["solver"]["mobility"]["high_field_gradient_discretization"] = variant
            csv_path = run_dir / "carrier_terms.csv"
            config["output_csv"] = str(csv_path.resolve())
            statuses[label][variant] = run_probe(
                args.runner, config, run_dir / "config.json"
            )
            rows_by_endpoint[label][variant] = read_csv(csv_path)
            artifacts[label][variant] = {
                "config": str((run_dir / "config.json").resolve()),
                "carrier_terms_csv": str(csv_path.resolve()),
                "carrier_terms_sha256": sha256(csv_path),
            }
        prior_baseline = Path(base["output_csv"]).resolve()
        artifacts[label]["baseline_reproduction"] = {
            "prior_csv_sha256": sha256(prior_baseline),
            "rerun_csv_sha256": artifacts[label]["edge_projection"]["carrier_terms_sha256"],
            "byte_identical": (
                sha256(prior_baseline)
                == artifacts[label]["edge_projection"]["carrier_terms_sha256"]
            ),
        }

    assert node_sets is not None
    endpoints = {
        label: analyze_endpoint(
            rows_by_endpoint[label]["edge_projection"],
            rows_by_endpoint[label]["transport_cell_vector"],
            node_sets,
        )
        for label, _ in args.endpoint
    }
    report = {
        "schema": "vela.templates_ldmos.g3_hfs_gradient_ab.v1",
        "status": "complete",
        "contract": {
            "task": "WP3 T3",
            "mode": "read_only_fixed_complete_sentaurus_state",
            "baseline": "edge_projection",
            "candidate": "transport_cell_vector",
            "probe": "newton_carrier_term_probe",
            "fixed_nodes": list(FIXED_NODES),
            "interface_band_definition": "all silicon nodes on a Si/non-Si interface plus their silicon one-ring",
            "all_silicon_definition": "all nodes incident to at least one silicon triangle, including production boundary rows",
            "main_gate": "both endpoints fixed-seven L2, max and every point <= 0.5x baseline",
            "migration_guard": "candidate-own top-7, complete interface band and all free-silicon L2/max do not worsen",
            "reclose": False,
            "curve_sweep": False,
            "ialmob": False,
            "predictor": False,
            "production_default_changed": False,
        },
        "inputs": inputs,
        "mesh_sets": {name: len(nodes) for name, nodes in node_sets.items()},
        "probe_status": statuses,
        "artifacts": artifacts,
        "endpoints": endpoints,
        "verdict": {
            "passes_main_gate": all(
                item["fixed_seven"]["passes_half_gate"] for item in endpoints.values()
            ),
            "passes_migration_guard": all(
                item["passes_migration_guard"] for item in endpoints.values()
            ),
        },
    }
    report["verdict"]["passes_t3_gate"] = (
        report["verdict"]["passes_main_gate"]
        and report["verdict"]["passes_migration_guard"]
    )
    analysis_dir = output / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    summary = analysis_dir / "summary.json"
    fixed_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    for label, endpoint in endpoints.items():
        for item in endpoint["fixed_seven"]["per_node"]:
            fixed_rows.append({"endpoint": label, **item})
        for name in ("fixed_seven", "candidate_top_seven", "interface_band", "all_silicon"):
            metric = endpoint[name]
            metric_rows.append({
                "endpoint": label,
                "set": name,
                "node_count": metric["candidate"]["node_count"],
                "baseline_l2": metric["baseline"]["l2"],
                "candidate_l2": metric["candidate"]["l2"],
                "l2_ratio": metric["ratio"]["l2"],
                "baseline_maximum_abs": metric["baseline"]["maximum_abs"],
                "candidate_maximum_abs": metric["candidate"]["maximum_abs"],
                "maximum_abs_ratio": metric["ratio"]["maximum_abs"],
                "not_worse": metric["not_worse"],
            })
    write_csv(analysis_dir / "fixed_seven.csv", fixed_rows)
    write_csv(analysis_dir / "metric_ratios.csv", metric_rows)
    report["analysis_artifacts"] = {
        "fixed_seven": str((analysis_dir / "fixed_seven.csv").resolve()),
        "metric_ratios": str((analysis_dir / "metric_ratios.csv").resolve()),
        "summary": str(summary.resolve()),
    }
    summary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": str(summary), **report["verdict"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
