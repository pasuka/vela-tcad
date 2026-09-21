#!/usr/bin/env python3
"""Evaluate the frozen WP3 T5 H1 magnitude-average control algebraically."""

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
Q_C = 1.602176634e-19


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


def run_probe(runner: Path, config: dict[str, Any], config_path: Path) -> dict[str, Any]:
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    env = dict(os.environ)
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        text=True,
        capture_output=True,
        env=env,
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


def finite(value: str, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {label}")
    return result


def triangle_gradient(
    nodes: dict[int, tuple[float, float]], values: dict[int, float], tri: list[int]
) -> tuple[float, float, float]:
    n0, n1, n2 = tri
    x0, y0 = nodes[n0]
    x1, y1 = nodes[n1]
    x2, y2 = nodes[n2]
    det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
    if abs(det) <= 1e-300:
        raise ValueError(f"degenerate triangle {tri}")
    dv10 = values[n1] - values[n0]
    dv20 = values[n2] - values[n0]
    gx = (dv10 * (y2 - y0) - dv20 * (y1 - y0)) / det
    gy = ((x1 - x0) * dv20 - (x2 - x0) * dv10) / det
    return gx, gy, abs(det)


def geometry(
    mesh_path: Path, state_path: Path, carrier_contact_names: set[str]
) -> dict[str, Any]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = {
        int(row["id"]): (float(row["x"]), float(row["y"]))
        for row in mesh["nodes"]
    }
    state_rows = read_csv(state_path)
    phin = {int(row["node_id"]): finite(row["phin"], "phin") for row in state_rows}
    phip = {int(row["node_id"]): finite(row["phip"], "phip") for row in state_rows}
    material = {
        int(row["id"]): str(row.get("material", "")).lower()
        for row in mesh["regions"]
    }
    silicon_regions = {
        region for region, name in material.items() if name in {"si", "silicon"}
    }
    contact_nodes = {
        int(node)
        for contact in mesh["contacts"]
        for node in contact.get("node_ids", [])
    }
    carrier_dirichlet_nodes = {
        int(node)
        for contact in mesh["contacts"]
        if str(contact.get("name", "")).lower() in carrier_contact_names
        for node in contact.get("node_ids", [])
    }
    silicon_nodes: set[int] = set()
    node_regions: dict[int, set[int]] = defaultdict(set)
    neighbors: dict[int, set[int]] = defaultdict(set)
    edge_cells: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for cell in mesh["triangles"]:
        region = int(cell["region_id"])
        tri = [int(node) for node in cell["node_ids"]]
        for node in tri:
            node_regions[node].add(region)
            neighbors[node].update(other for other in tri if other != node)
        if region not in silicon_regions:
            continue
        silicon_nodes.update(tri)
        ngx, ngy, area2 = triangle_gradient(nodes, phin, tri)
        pgx, pgy, _ = triangle_gradient(nodes, phip, tri)
        item = {
            "area2": area2,
            "electron_gradient_norm_V_per_um": math.hypot(ngx, ngy),
            "hole_gradient_norm_V_per_um": math.hypot(pgx, pgy),
            "touches_contact": bool(set(tri) & contact_nodes),
        }
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            edge_cells[tuple(sorted((a, b)))].append(item)
    direct_interface = {
        node
        for node in silicon_nodes
        if any(region not in silicon_regions for region in node_regions[node])
    }
    interface_band = set(direct_interface)
    for node in direct_interface:
        interface_band.update(neighbors[node] & silicon_nodes)
    return {
        "edge_cells": edge_cells,
        "contact_nodes": contact_nodes,
        "carrier_dirichlet_nodes": carrier_dirichlet_nodes,
        "silicon_nodes": silicon_nodes,
        "direct_interface": direct_interface,
        "interface_band": interface_band,
    }


def magnitude_average_drive(
    cells: list[dict[str, Any]], carrier: str
) -> tuple[float, bool]:
    if not cells:
        return 0.0, False
    total_area = sum(cell["area2"] for cell in cells)
    field = sum(
        cell["area2"] * cell[f"{carrier}_gradient_norm_V_per_um"] for cell in cells
    ) / total_area
    return field * 1e6, any(cell["touches_contact"] for cell in cells)


def field_limited_mobility(mu0: float, field: float, vsat: float, beta: float) -> float:
    if mu0 <= 0.0:
        return 0.0
    if field <= 0.0:
        return mu0
    ratio = mu0 * abs(field) / vsat
    return mu0 / (1.0 + ratio**beta) ** (1.0 / beta)


def residual_map(rows: list[dict[str, str]]) -> dict[int, float]:
    return {int(row["node_id"]): finite(row["electron_residual"], "residual") for row in rows}


def norm(values: list[float]) -> dict[str, float]:
    return {
        "l2": math.sqrt(sum(value * value for value in values)),
        "maximum_abs": max((abs(value) for value in values), default=0.0),
    }


def ratio(candidate: float, baseline: float) -> float:
    if baseline == 0.0:
        return 0.0 if candidate == 0.0 else math.inf
    return candidate / baseline


def no_worse(candidate: float, baseline: float) -> bool:
    return candidate <= baseline * (1.0 + 1e-12) + 1e-30


def node_metric(values: dict[int, float], nodes: set[int] | tuple[int, ...]) -> dict[str, float | int]:
    return {"node_count": len(nodes), **norm([values[node] for node in sorted(nodes)])}


def top_nodes(values: dict[int, float], eligible: set[int], count: int = 7) -> list[int]:
    return sorted(eligible, key=lambda node: (-abs(values[node]), node))[:count]


def compare_metrics(candidate: dict[str, Any], baseline: dict[str, Any]) -> dict[str, Any]:
    return {
        "baseline": baseline,
        "candidate": candidate,
        "ratio": {
            "l2": ratio(candidate["l2"], baseline["l2"]),
            "maximum_abs": ratio(candidate["maximum_abs"], baseline["maximum_abs"]),
        },
        "not_worse": no_worse(candidate["l2"], baseline["l2"])
        and no_worse(candidate["maximum_abs"], baseline["maximum_abs"]),
    }


def evaluate_gate(
    baseline: dict[int, float], candidate: dict[int, float], sets: dict[str, Any]
) -> dict[str, Any]:
    per_node = [
        {
            "node_id": node,
            "baseline": baseline[node],
            "candidate": candidate[node],
            "absolute_ratio": ratio(abs(candidate[node]), abs(baseline[node])),
            "passes_half_gate": abs(candidate[node]) <= 0.5 * abs(baseline[node]),
        }
        for node in FIXED_NODES
    ]
    fixed = compare_metrics(
        node_metric(candidate, FIXED_NODES), node_metric(baseline, FIXED_NODES)
    )
    fixed["per_node"] = per_node
    fixed["passes_half_gate"] = (
        fixed["ratio"]["l2"] <= 0.5
        and fixed["ratio"]["maximum_abs"] <= 0.5
        and all(row["passes_half_gate"] for row in per_node)
    )
    baseline_top = top_nodes(baseline, sets["silicon_nodes"])
    candidate_top = top_nodes(candidate, sets["silicon_nodes"])
    candidate_hotspot = compare_metrics(
        node_metric(candidate, set(candidate_top)), node_metric(baseline, set(candidate_top))
    )
    interface = compare_metrics(
        node_metric(candidate, sets["interface_band"]),
        node_metric(baseline, sets["interface_band"]),
    )
    silicon = compare_metrics(
        node_metric(candidate, sets["silicon_nodes"]),
        node_metric(baseline, sets["silicon_nodes"]),
    )
    migration = all(row["not_worse"] for row in (candidate_hotspot, interface, silicon))
    return {
        "fixed_seven": fixed,
        "candidate_top_seven": {"node_ids": candidate_top, **candidate_hotspot},
        "interface_band": interface,
        "all_silicon": silicon,
        "hotspot_migration": {
            "baseline_top_seven": baseline_top,
            "candidate_top_seven": candidate_top,
            "overlap": sorted(set(baseline_top) & set(candidate_top)),
            "removed": sorted(set(baseline_top) - set(candidate_top)),
            "added": sorted(set(candidate_top) - set(baseline_top)),
        },
        "passes_migration_guard": migration,
        "passes_t5_gate": fixed["passes_half_gate"] and migration,
    }


def reconstruct_transport(
    sg_rows: list[dict[str, str]], node_count: int, contact_nodes: set[int]
) -> tuple[dict[int, float], dict[int, float]]:
    sums = {node: 0.0 for node in range(node_count)}
    abs_sums = {node: 0.0 for node in range(node_count)}
    for row in sg_rows:
        n0 = int(row["node0"])
        n1 = int(row["node1"])
        flux = finite(row["electron_flux"], "edge flux")
        sums[n0] += flux
        sums[n1] -= flux
        abs_sums[n0] += abs(flux)
        abs_sums[n1] += abs(flux)
    # The production carrier-term diagnostic applies the Dirichlet carrier-row
    # replacement after assembly.  Reproduce that exact row semantics here.
    for node in contact_nodes:
        sums[node] = 0.0
        abs_sums[node] = 0.0
    return sums, abs_sums


def reconstruction_error(
    rows: list[dict[str, str]], reconstructed: dict[int, float], physical_scale: float
) -> dict[str, Any]:
    diffs = []
    targets = []
    max_node = -1
    max_abs = -1.0
    for row in rows:
        node = int(row["node_id"])
        target = finite(row["electron_flux"], "carrier transport")
        diff = reconstructed[node] - target
        diffs.append(diff)
        targets.append(target)
        if abs(diff) > max_abs:
            max_abs = abs(diff)
            max_node = node
    l2 = math.sqrt(sum(value * value for value in diffs))
    target_l2 = math.sqrt(sum(value * value for value in targets))
    return {
        "relative_l2": l2 / max(target_l2, 1e-300),
        "maximum_absolute_raw": max_abs,
        "maximum_absolute_A_per_um": max_abs * physical_scale,
        "maximum_error_node": max_node,
        "passes": (
            l2 / max(target_l2, 1e-300) <= 1e-10
            and max_abs * physical_scale <= 1e-18
        ),
    }


def particle_flux_physical_scale(sg_rows: list[dict[str, str]]) -> dict[str, float]:
    ratios = []
    for row in sg_rows:
        raw = finite(row["electron_flux"], "raw edge flux")
        if raw == 0.0:
            continue
        physical = (
            finite(row["electron_particle_line_flux_per_m_s"], "particle line flux")
            * Q_C
            * 1e-6
        )
        ratios.append(physical / raw)
    if not ratios:
        raise ValueError("no nonzero electron edge flux for physical scaling")
    reference = ratios[0]
    maximum_relative_spread = max(
        abs(value - reference) / max(abs(value), abs(reference), 1e-300)
        for value in ratios
    )
    if maximum_relative_spread > 1e-12:
        raise ValueError(
            f"particle-line physical scale is not constant: {maximum_relative_spread}"
        )
    return {
        "A_per_um_per_raw": reference,
        "maximum_relative_spread": maximum_relative_spread,
    }


def parse_endpoint(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("endpoint must be LABEL=CONFIG")
    label, path = value.split("=", 1)
    return label, Path(path)


def validate_contract(config: dict[str, Any]) -> None:
    if config.get("simulation_type") != "newton_carrier_term_probe":
        raise ValueError("T5 requires newton_carrier_term_probe endpoint inputs")
    mobility = config["solver"]["mobility"]
    if mobility.get("model") != "constant_field":
        raise ValueError("T5 algebraic control is frozen to constant_field HFS")
    if mobility.get("high_field_driving_force") != "quasi_fermi_gradient":
        raise ValueError("T5 requires GradQF HFS")
    if not mobility.get("contact_electric_field_fallback", False):
        raise ValueError("T5 requires the qualified contact fallback")
    if "predictor" in config.get("sweep", {}):
        raise ValueError("T5 requires predictor disabled")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--t3-root", type=Path, required=True)
    parser.add_argument("--endpoint", type=parse_endpoint, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.endpoint) != 2:
        raise ValueError("T5 requires exactly two endpoints")

    protocol_path = args.protocol.resolve()
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    candidates = protocol.get("candidate_set", [])
    if protocol.get("status") != "frozen_before_result_calculation":
        raise ValueError("T5 protocol is not frozen")
    if [row.get("id") for row in candidates] != ["transport_cell_magnitude_average"]:
        raise ValueError("T5 protocol candidate set differs from the frozen single variant")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    t3_root = args.t3_root.resolve()
    endpoint_results: dict[str, Any] = {}
    probe_status: dict[str, Any] = {}
    artifact_map: dict[str, Any] = {}
    edge_output_rows: list[dict[str, Any]] = []
    node_output_rows: list[dict[str, Any]] = []
    reconstruction_rows: list[dict[str, Any]] = []

    for label, config_path_arg in args.endpoint:
        config_path = config_path_arg.resolve()
        base = json.loads(config_path.read_text(encoding="utf-8"))
        validate_contract(base)
        registered = next(
            row for row in protocol["endpoints"] if row["label"] == label
        )
        if sha256(config_path) != registered["config_sha256"]:
            raise ValueError(f"endpoint config hash mismatch for {label}")
        mesh_path = Path(base["mesh_file"]).resolve()
        state_path = Path(base["state_file"]).resolve()
        carrier_contact_names = {
            str(contact["name"]).lower()
            for contact in base["contacts"]
            if str(contact.get("type", "")).lower() == "ohmic"
        }
        sets = geometry(mesh_path, state_path, carrier_contact_names)
        carrier_rows: dict[str, list[dict[str, str]]] = {}
        sg_rows: dict[str, list[dict[str, str]]] = {}
        probe_status[label] = {}
        artifact_map[label] = {}

        for variant in VARIANTS:
            carrier_path = t3_root / "probes" / label / variant / "carrier_terms.csv"
            carrier_config_path = t3_root / "probes" / label / variant / "config.json"
            carrier_config = json.loads(carrier_config_path.read_text(encoding="utf-8"))
            if Path(carrier_config["state_file"]).resolve() != state_path:
                raise ValueError(f"T3 carrier state mismatch for {label}/{variant}")
            carrier_rows[variant] = read_csv(carrier_path)

            run_dir = output / "probes" / label / variant
            run_dir.mkdir(parents=True, exist_ok=True)
            sg_config = deepcopy(base)
            sg_config["simulation_type"] = "sg_edge_flux_probe"
            sg_config["solver"]["mobility"]["high_field_gradient_discretization"] = variant
            sg_path = run_dir / "sg_edges.csv"
            sg_config["output_csv"] = str(sg_path)
            probe_status[label][variant] = run_probe(
                args.runner, sg_config, run_dir / "sg_config.json"
            )
            sg_rows[variant] = read_csv(sg_path)
            artifact_map[label][variant] = {
                "carrier_terms": str(carrier_path.resolve()),
                "carrier_terms_sha256": sha256(carrier_path),
                "sg_edges": str(sg_path.resolve()),
                "sg_edges_sha256": sha256(sg_path),
                "sg_config": str((run_dir / "sg_config.json").resolve()),
            }

        if len(sg_rows[VARIANTS[0]]) != len(sg_rows[VARIANTS[1]]):
            raise ValueError(f"SG edge-count mismatch for {label}")
        node_count = len(carrier_rows[VARIANTS[0]])
        reconstruction: dict[str, Any] = {}
        for variant in VARIANTS:
            physical_scale = particle_flux_physical_scale(sg_rows[variant])
            sums, _ = reconstruct_transport(
                sg_rows[variant], node_count, sets["carrier_dirichlet_nodes"]
            )
            error = reconstruction_error(
                carrier_rows[variant], sums, physical_scale["A_per_um_per_raw"]
            )
            error["physical_scale"] = physical_scale
            reconstruction[variant] = error
            reconstruction_rows.append({"endpoint": label, "variant": variant, **error})
        if not all(row["passes"] for row in reconstruction.values()):
            raise RuntimeError(f"T5 extreme reconstruction gate failed for {label}")

        edge_by_variant = {
            variant: {int(row["edge_id"]): row for row in rows}
            for variant, rows in sg_rows.items()
        }
        if set(edge_by_variant[VARIANTS[0]]) != set(edge_by_variant[VARIANTS[1]]):
            raise ValueError(f"SG edge identity mismatch for {label}")
        mobility = base["solver"]["mobility"]
        beta = float(mobility["electron_high_field_beta"])
        vsat = float(mobility["electron_saturation_velocity_m_s"]) * 1e-2
        candidate_flux: dict[int, float] = {}
        max_formula_rel = 0.0
        max_dual_candidate_rel = 0.0
        for edge_id in sorted(edge_by_variant[VARIANTS[0]]):
            edge = edge_by_variant[VARIANTS[0]][edge_id]
            global_edge = edge_by_variant[VARIANTS[1]][edge_id]
            for key in ("node0", "node1", "length_m", "couple_m", "psi0_V", "psi1_V", "phin0_V", "phin1_V", "electron_bernoulli_argument", "electron_generalized_einstein_factor"):
                if edge[key] != global_edge[key]:
                    raise ValueError(f"non-mobility SG field changed at {label} edge {edge_id}: {key}")
            n0 = int(edge["node0"])
            n1 = int(edge["node1"])
            cells = sets["edge_cells"].get(tuple(sorted((n0, n1))), [])
            drive, contact_fallback = magnitude_average_drive(cells, "electron")
            if contact_fallback:
                drive = finite(edge["electron_mobility_field_V_m"], "contact drive")
                if edge["electron_mobility_field_V_m"] != global_edge["electron_mobility_field_V_m"]:
                    raise ValueError(f"contact fallback differs between extremes at {label} edge {edge_id}")
            mu0 = 0.1417 if cells else 0.0
            candidate_mu = field_limited_mobility(mu0, drive, vsat, beta)
            candidates_from_extreme = []
            for variant, row in ((VARIANTS[0], edge), (VARIANTS[1], global_edge)):
                observed_mu = finite(row["electron_mobility_m2_V_s"], "mobility")
                observed_drive = finite(row["electron_mobility_field_V_m"], "drive")
                expected_mu = field_limited_mobility(mu0, observed_drive, vsat, beta)
                formula_rel = abs(observed_mu - expected_mu) / max(abs(observed_mu), abs(expected_mu), 1e-300)
                max_formula_rel = max(max_formula_rel, formula_rel)
                observed_flux = finite(row["electron_flux"], "edge flux")
                candidates_from_extreme.append(
                    observed_flux * candidate_mu / observed_mu if observed_mu > 0.0 else 0.0
                )
            candidate_edge_flux = candidates_from_extreme[0]
            dual_rel = abs(candidates_from_extreme[0] - candidates_from_extreme[1]) / max(
                abs(candidates_from_extreme[0]), abs(candidates_from_extreme[1]), 1e-300
            )
            max_dual_candidate_rel = max(max_dual_candidate_rel, dual_rel)
            candidate_flux[edge_id] = candidate_edge_flux
            edge_output_rows.append({
                "endpoint": label,
                "edge_id": edge_id,
                "node0": n0,
                "node1": n1,
                "adjacent_silicon_cells": len(cells),
                "contact_fallback": contact_fallback,
                "edge_projection_drive_V_m": edge["electron_mobility_field_V_m"],
                "transport_cell_vector_drive_V_m": global_edge["electron_mobility_field_V_m"],
                "magnitude_average_drive_V_m": drive,
                "candidate_mobility_m2_V_s": candidate_mu,
                "candidate_flux": candidate_edge_flux,
                "candidate_flux_dual_relative_error": dual_rel,
            })
        if max_formula_rel > 1e-12:
            raise RuntimeError(f"HFS formula reconstruction failed for {label}: {max_formula_rel}")
        if max_dual_candidate_rel > 1e-12:
            raise RuntimeError(f"dual-extreme candidate reconstruction failed for {label}: {max_dual_candidate_rel}")

        baseline_values = residual_map(carrier_rows[VARIANTS[0]])
        candidate_values = dict(baseline_values)
        baseline_edges = edge_by_variant[VARIANTS[0]]
        delta = {node: 0.0 for node in range(node_count)}
        for edge_id, flux in candidate_flux.items():
            row = baseline_edges[edge_id]
            edge_delta = flux - finite(row["electron_flux"], "baseline flux")
            n0 = int(row["node0"])
            n1 = int(row["node1"])
            delta[n0] += edge_delta
            delta[n1] -= edge_delta
        for node in range(node_count):
            if node not in sets["carrier_dirichlet_nodes"]:
                candidate_values[node] += delta[node]
            node_output_rows.append({
                "endpoint": label,
                "node_id": node,
                "baseline_residual": baseline_values[node],
                "transport_delta": 0.0 if node in sets["carrier_dirichlet_nodes"] else delta[node],
                "candidate_residual": candidate_values[node],
                "is_contact": node in sets["contact_nodes"],
                "is_carrier_dirichlet": node in sets["carrier_dirichlet_nodes"],
            })
        gate = evaluate_gate(baseline_values, candidate_values, sets)
        endpoint_results[label] = {
            "extreme_reconstruction": reconstruction,
            "hfs_formula_max_relative_error": max_formula_rel,
            "dual_extreme_candidate_max_relative_error": max_dual_candidate_rel,
            "gate": gate,
        }

    verdict = {
        "passes_main_gate": all(
            row["gate"]["fixed_seven"]["passes_half_gate"]
            for row in endpoint_results.values()
        ),
        "passes_migration_guard": all(
            row["gate"]["passes_migration_guard"] for row in endpoint_results.values()
        ),
    }
    verdict["passes_t5_gate"] = verdict["passes_main_gate"] and verdict["passes_migration_guard"]
    verdict["next_action"] = (
        "nominate_for_cpp_and_t6_only" if verdict["passes_t5_gate"]
        else "stop_t5_h1_targeted_control_no_t6"
    )
    analysis = output / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    write_csv(analysis / "candidate_edges.csv", edge_output_rows)
    write_csv(analysis / "candidate_nodes.csv", node_output_rows)
    write_csv(analysis / "extreme_reconstruction.csv", reconstruction_rows)
    report = {
        "schema": "vela.templates_ldmos.g3_wp3_t5_hfs_targeted.v1",
        "status": "complete",
        "protocol": str(protocol_path),
        "protocol_sha256": sha256(protocol_path),
        "candidate": "transport_cell_magnitude_average",
        "contract": {
            "mode": "read_only_fixed_complete_sentaurus_state_algebraic_reconstruction",
            "reclose": False,
            "curve_sweep": False,
            "production_default_changed": False,
            "candidate_count": 1,
        },
        "probe_status": probe_status,
        "artifacts": artifact_map,
        "endpoints": endpoint_results,
        "verdict": verdict,
        "analysis_artifacts": {
            "candidate_edges": str((analysis / "candidate_edges.csv").resolve()),
            "candidate_nodes": str((analysis / "candidate_nodes.csv").resolve()),
            "extreme_reconstruction": str((analysis / "extreme_reconstruction.csv").resolve()),
        },
    }
    summary = analysis / "summary.json"
    report["analysis_artifacts"]["summary"] = str(summary.resolve())
    summary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": str(summary), **verdict}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
