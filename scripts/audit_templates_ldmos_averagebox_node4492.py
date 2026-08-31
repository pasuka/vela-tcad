#!/usr/bin/env python3
"""Audit direct Sentaurus AverageBox geometry at LDMOS node 4492.

The script regenerates read-only Vela SG and first-step probes from the frozen
SSS state, maps the grid-numbered Sentaurus Measure/Coefficients debug file to
the exact imported triangles, and reassembles the target one-ring without
changing the production mesh or solver policy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any


DEBUG_LINE = re.compile(
    r"^\s*(?P<grd>\d+)\s+(?P<des>-?\d+)\s+(?P<type>\d+)\s+"
    r"(?P<values>[^#]+?)\s*$"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_debug_block(path: Path, name: str) -> dict[int, list[float]]:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf"\n\s*{re.escape(name)}\s*\{{.*?\n(?P<body>.*?)\n\s*\}}",
        text,
        re.DOTALL,
    )
    if match is None:
        raise ValueError(f"{name} block not found in {path}")
    result: dict[int, list[float]] = {}
    for line in match.group("body").splitlines():
        parsed = DEBUG_LINE.match(line)
        if parsed is None or int(parsed.group("type")) != 2:
            continue
        design = int(parsed.group("des"))
        if design >= 0:
            result[design] = [
                float(value) for value in parsed.group("values").split()
            ]
    return result


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def coefficient_edge(nodes: list[int], coefficient_index: int) -> tuple[int, int, int]:
    """Return (edge node 0, edge node 1, opposite local vertex).

    MeasureCoefficients.debug follows the TDR triangle-edge ordering, not the
    local-vertex ordering used by Measure.  For a triangle [v0, v1, v2], its
    three coefficient entries correspond to edges [v2-v0, v1-v2, v0-v1].
    The mapping is independently qualified below by comparing ordinary acute
    cells against their analytic cotangent coefficients.
    """
    mapping = ((2, 0, 1), (1, 2, 0), (0, 1, 2))
    i, j, opposite = mapping[coefficient_index]
    return nodes[i], nodes[j], opposite


def angle_degrees(
    center: tuple[float, float], left: tuple[float, float], right: tuple[float, float]
) -> float:
    u = (left[0] - center[0], left[1] - center[1])
    v = (right[0] - center[0], right[1] - center[1])
    scale = math.hypot(*u) * math.hypot(*v)
    cosine = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / scale))
    return math.degrees(math.acos(cosine))


def local_vela_couple(
    a: tuple[float, float], b: tuple[float, float], opposite: tuple[float, float]
) -> tuple[float, str]:
    edge = (b[0] - a[0], b[1] - a[1])
    length = math.hypot(*edge)
    u = (a[0] - opposite[0], a[1] - opposite[1])
    v = (b[0] - opposite[0], b[1] - opposite[1])
    cross = abs(u[0] * v[1] - u[1] * v[0])
    cotangent = (u[0] * v[0] + u[1] * v[1]) / cross
    if cotangent >= 0.0:
        return 0.5 * cotangent * length, "cotangent"
    area = 0.5 * cross
    return area / (3.0 * length), "positive_barycentric_fallback"


def means(n0: float, n1: float) -> dict[str, float]:
    result = {
        "arithmetic": 0.5 * (n0 + n1),
        "geometric": math.sqrt(n0 * n1),
        "harmonic": 2.0 * n0 * n1 / (n0 + n1),
    }
    result["logarithmic"] = (
        n0 if n0 == n1 else (n1 - n0) / math.log(n1 / n0)
    )
    return result


def norm(values: list[float]) -> dict[str, float]:
    return {
        "l1": sum(abs(value) for value in values),
        "l2": math.sqrt(sum(value * value for value in values)),
        "maximum_abs": max((abs(value) for value in values), default=0.0),
    }


def read_scalar_field(root: Path, name: str) -> dict[int, float]:
    rows = read_csv(root / "fields" / f"{name}_region0.csv")
    return {int(row["node_id"]): float(row["component0"]) for row in rows}


def first_step_alignment(
    vela_rows: list[dict[str, Any]],
    initial_root: Path,
    newton_root: Path,
    poststep_root: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    initial = {
        name: read_scalar_field(initial_root, name)
        for name in ("ElectrostaticPotential", "eQuasiFermiPotential",
                     "hQuasiFermiPotential")
    }
    post = {
        name: read_scalar_field(poststep_root, name)
        for name in ("ElectrostaticPotential", "eQuasiFermiPotential",
                     "hQuasiFermiPotential", "eDensity", "hDensity")
    }
    newton = {
        name: read_scalar_field(newton_root, name)
        for name in ("NewtonStepElectrostaticPotentialUpdate",
                     "NewtonStepEDensityUpdate", "NewtonStepHDensityUpdate",
                     "PoissonRhs", "eContinuityRhs", "hContinuityRhs")
    }
    rows = []
    for vela in vela_rows:
        node = int(vela["node_id"])
        rows.append({
            "node_id": node,
            "is_target": bool(vela["is_target"]),
            "vela_delta_psi_V": float(vela["baseline_delta_psi_V"]),
            "vela_delta_phin_V": float(vela["baseline_delta_phin_V"]),
            "sentaurus_delta_psi_V": (
                post["ElectrostaticPotential"][node]
                - initial["ElectrostaticPotential"][node]
            ),
            "sentaurus_newton_step_psi_update_V":
                newton["NewtonStepElectrostaticPotentialUpdate"][node],
            "sentaurus_delta_phin_V": (
                post["eQuasiFermiPotential"][node]
                - initial["eQuasiFermiPotential"][node]
            ),
            "sentaurus_delta_phip_V": (
                post["hQuasiFermiPotential"][node]
                - initial["hQuasiFermiPotential"][node]
            ),
            "sentaurus_newton_step_e_density_update_cm3":
                newton["NewtonStepEDensityUpdate"][node],
            "sentaurus_newton_step_h_density_update_cm3":
                newton["NewtonStepHDensityUpdate"][node],
            "sentaurus_post_e_density_cm3": post["eDensity"][node],
            "sentaurus_post_h_density_cm3": post["hDensity"][node],
            "sentaurus_poisson_rhs_C": newton["PoissonRhs"][node],
            "sentaurus_e_continuity_rhs_C": newton["eContinuityRhs"][node],
            "sentaurus_h_continuity_rhs_C": newton["hContinuityRhs"][node],
        })
    target = next(row for row in rows if row["is_target"])
    summary = {
        "target": target,
        "one_ring_max_abs_sentaurus_delta_psi_V": max(
            abs(row["sentaurus_delta_psi_V"]) for row in rows
        ),
        "one_ring_max_abs_sentaurus_delta_phin_V": max(
            abs(row["sentaurus_delta_phin_V"]) for row in rows
        ),
        "one_ring_max_abs_sentaurus_delta_phip_V": max(
            abs(row["sentaurus_delta_phip_V"]) for row in rows
        ),
    }
    return rows, summary


def run_probe(
    runner: Path, base: dict[str, Any], output: Path, simulation_type: str,
) -> tuple[Path, dict[str, Any]]:
    csv_path = output / f"{simulation_type}.csv"
    config_path = output / f"{simulation_type}.json"
    config = json.loads(json.dumps(base))
    config["simulation_type"] = simulation_type
    config["output_csv"] = str(csv_path.resolve())
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    status = {
        "return_code": completed.returncode,
        "stdout": completed.stdout,
        "config": str(config_path.resolve()),
        "output_csv": str(csv_path.resolve()),
    }
    if completed.returncode != 0 or not csv_path.is_file():
        raise RuntimeError(json.dumps(status, indent=2))
    return csv_path, status


def audit(
    mesh: dict[str, Any], sg_rows: list[dict[str, str]],
    carrier_rows: list[dict[str, str]], step_rows: list[dict[str, str]],
    measures: dict[int, list[float]], coefficients: dict[int, list[float]],
    target: int,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    points = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    cells = {
        int(cell["id"]): [int(value) for value in cell["node_ids"]]
        for cell in mesh["triangles"]
    }
    incident = {cell: nodes for cell, nodes in cells.items() if target in nodes}
    if not incident:
        raise ValueError(f"target node {target} has no incident triangles")
    missing = sorted(set(incident) - set(coefficients))
    if missing:
        raise ValueError(f"AverageBox debug lacks target cells: {missing}")

    cell_output: list[dict[str, Any]] = []
    sent_couple_um: dict[tuple[int, int], float] = defaultdict(float)
    vela_couple_um: dict[tuple[int, int], float] = defaultdict(float)
    target_measure_um2 = 0.0
    for cell, nodes in sorted(incident.items()):
        local_points = [points[node] for node in nodes]
        angles = [
            angle_degrees(local_points[i], local_points[(i + 1) % 3],
                          local_points[(i + 2) % 3])
            for i in range(3)
        ]
        target_local = nodes.index(target)
        target_measure_um2 += measures[cell][target_local]
        for coefficient_index in range(3):
            a, b, opposite = coefficient_edge(nodes, coefficient_index)
            key = edge_key(a, b)
            length_um = math.dist(points[a], points[b])
            sent_local = coefficients[cell][coefficient_index] * length_um
            vela_local, policy = local_vela_couple(
                points[a], points[b], points[nodes[opposite]]
            )
            sent_couple_um[key] += sent_local
            vela_couple_um[key] += vela_local
            cell_output.append({
                "cell_id": cell,
                "coefficient_index": coefficient_index,
                "local_opposite": opposite,
                "opposite_node": nodes[opposite],
                "edge_node0": a,
                "edge_node1": b,
                "edge_incident_to_target": target in key,
                "maximum_angle_deg": max(angles),
                "target_local_measure_um2": measures[cell][target_local],
                "averagebox_coefficient": coefficients[cell][coefficient_index],
                "averagebox_local_couple_um": sent_local,
                "vela_local_couple_um": vela_local,
                "vela_policy": policy,
            })

    sg_by_key = {edge_key(int(row["node0"]), int(row["node1"])): row for row in sg_rows}
    carrier_by_node = {int(row["node_id"]): row for row in carrier_rows}
    step_by_node = {int(row["node_id"]): row for row in step_rows}
    one_ring = {target}
    for nodes in incident.values():
        one_ring.update(nodes)
    variants = ["averagebox_current_sg", "averagebox_arithmetic",
                "averagebox_geometric", "averagebox_logarithmic",
                "averagebox_harmonic"]
    delta_by_variant: dict[str, dict[int, float]] = {
        name: defaultdict(float) for name in variants
    }
    edge_output: list[dict[str, Any]] = []
    reproduction_errors = []
    density_errors: dict[str, list[float]] = defaultdict(list)
    for neighbor in sorted(one_ring - {target}):
        key = edge_key(target, neighbor)
        row = sg_by_key[key]
        reported_couple = float(row["couple_m"])
        reconstructed = vela_couple_um[key] * 1.0e-6
        reproduction_errors.append(abs(reconstructed - reported_couple))
        sent_couple = sent_couple_um[key] * 1.0e-6
        per_couple = float(row["electron_scaled_flux_per_couple_m"])
        baseline_flux = float(row["electron_flux"])
        current_candidate = per_couple * sent_couple
        target_sign = 1.0 if int(row["node0"]) == target else -1.0
        physical_density_flux = float(
            row["electron_particle_flux_density_per_couple_m2_s"]
        )
        qf_gradient = (
            (float(row["phin1_V"]) - float(row["phin0_V"]))
            / float(row["length_m"])
        )
        mobility = float(row["electron_mobility_m2_V_s"])
        effective_density = (
            abs(physical_density_flux / (mobility * qf_gradient))
            if mobility != 0.0 and qf_gradient != 0.0 else math.nan
        )
        simple = means(
            float(row["electron_density0_m3"]),
            float(row["electron_density1_m3"]),
        )
        candidate_fluxes = {"averagebox_current_sg": current_candidate}
        for name, density in simple.items():
            candidate_fluxes[f"averagebox_{name}"] = (
                current_candidate * density / effective_density
                if effective_density > 0.0 else current_candidate
            )
            density_errors[name].append(abs(density / effective_density - 1.0))
        for name, candidate in candidate_fluxes.items():
            delta = candidate - baseline_flux
            node0, node1 = int(row["node0"]), int(row["node1"])
            delta_by_variant[name][node0] += delta
            delta_by_variant[name][node1] -= delta
        edge_output.append({
            "edge_id": int(row["edge_id"]),
            "node0": int(row["node0"]),
            "node1": int(row["node1"]),
            "target_sign": target_sign,
            "reported_vela_couple_m": reported_couple,
            "reconstructed_vela_couple_m": reconstructed,
            "sentaurus_averagebox_couple_m": sent_couple,
            "averagebox_over_vela_couple": (
                sent_couple / reported_couple if reported_couple else math.inf
            ),
            "baseline_electron_flux": baseline_flux,
            "averagebox_current_sg_flux": current_candidate,
            "electron_effective_sg_density_m3": effective_density,
            **{f"{name}_density_m3": value for name, value in simple.items()},
        })

    baseline_ring = []
    variant_ring: dict[str, list[float]] = {name: [] for name in variants}
    node_rows = []
    for node in sorted(one_ring):
        raw = carrier_by_node[node]
        baseline = float(raw["electron_residual"])
        baseline_ring.append(baseline)
        record: dict[str, Any] = {
            "node_id": node,
            "is_target": node == target,
            "baseline_electron_residual": baseline,
            "baseline_delta_psi_V": float(step_by_node[node]["delta_psi_V"]),
            "baseline_delta_phin_V": float(step_by_node[node]["delta_phin_V"]),
            "baseline_delta_phip_V": float(step_by_node[node]["delta_phip_V"]),
            "baseline_trial_phin_residual": float(
                step_by_node[node]["trial_phin_residual"]
            ),
        }
        for name in variants:
            candidate = baseline + delta_by_variant[name][node]
            record[f"{name}_electron_residual"] = candidate
            variant_ring[name].append(candidate)
        node_rows.append(record)

    baseline_norm = norm(baseline_ring)
    variants_summary = {}
    for name in variants:
        candidate_norm = norm(variant_ring[name])
        target_candidate = next(
            row[f"{name}_electron_residual"] for row in node_rows
            if row["is_target"]
        )
        variants_summary[name] = {
            "one_ring_residual": candidate_norm,
            "one_ring_l2_over_baseline": (
                candidate_norm["l2"] / baseline_norm["l2"]
            ),
            "target_residual": target_candidate,
            "target_abs_over_baseline": (
                abs(target_candidate)
                / abs(float(carrier_by_node[target]["electron_residual"]))
            ),
        }
    direct = variants_summary["averagebox_current_sg"]
    gate = (
        direct["one_ring_l2_over_baseline"] <= 0.5
        and direct["target_abs_over_baseline"] <= 0.5
    )
    summary = {
        "schema": "vela.templates_ldmos.averagebox_node4492_audit.v1",
        "target_node": target,
        "incident_cells": sorted(incident),
        "one_ring_nodes": sorted(one_ring),
        "sentaurus_averagebox_target_measure_um2": target_measure_um2,
        "maximum_absolute_vela_couple_reproduction_error_m": max(
            reproduction_errors, default=0.0
        ),
        "baseline_one_ring_residual": baseline_norm,
        "variants": variants_summary,
        "simple_density_mean_median_relative_errors": {
            name: sorted(values)[len(values) // 2]
            for name, values in density_errors.items()
        },
        "baseline_first_step_target": {
            "delta_psi_V": float(step_by_node[target]["delta_psi_V"]),
            "delta_phin_V": float(step_by_node[target]["delta_phin_V"]),
            "delta_phip_V": float(step_by_node[target]["delta_phip_V"]),
            "phin_residual": float(step_by_node[target]["phin_residual"]),
            "trial_phin_residual": float(step_by_node[target]["trial_phin_residual"]),
        },
        "averagebox_first_step_implementation_gate": {
            "requires_target_abs_ratio_at_most": 0.5,
            "requires_one_ring_l2_ratio_at_most": 0.5,
            "passed": gate,
        },
    }
    return summary, cell_output, edge_output + node_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--carrier-terms", type=Path, required=True)
    parser.add_argument("--measure-coefficients", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-initial-import", type=Path)
    parser.add_argument("--sentaurus-newton-import", type=Path)
    parser.add_argument("--sentaurus-poststep-import", type=Path)
    parser.add_argument("--target-node", type=int, default=4492)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base = json.loads(args.base_config.read_text(encoding="utf-8"))
    sg_path, sg_status = run_probe(args.runner.resolve(), base, output, "sg_edge_flux_probe")
    step_path, step_status = run_probe(args.runner.resolve(), base, output, "newton_step_probe")
    summary, cells, combined = audit(
        json.loads(args.mesh.read_text(encoding="utf-8")),
        read_csv(sg_path), read_csv(args.carrier_terms), read_csv(step_path),
        parse_debug_block(args.measure_coefficients, "Measure"),
        parse_debug_block(args.measure_coefficients, "Coefficients"),
        args.target_node,
    )
    edges = [row for row in combined if "edge_id" in row]
    nodes = [row for row in combined if "node_id" in row]
    write_csv(output / "cell_local_coefficients.csv", cells)
    write_csv(output / "target_edges.csv", edges)
    write_csv(output / "one_ring_nodes.csv", nodes)
    sentaurus_roots = (
        args.sentaurus_initial_import,
        args.sentaurus_newton_import,
        args.sentaurus_poststep_import,
    )
    if any(path is not None for path in sentaurus_roots):
        if not all(path is not None for path in sentaurus_roots):
            raise ValueError("all three Sentaurus import directories are required")
        alignment, alignment_summary = first_step_alignment(
            nodes,
            args.sentaurus_initial_import.resolve(),
            args.sentaurus_newton_import.resolve(),
            args.sentaurus_poststep_import.resolve(),
        )
        write_csv(output / "first_step_alignment.csv", alignment)
        summary["first_step_alignment"] = alignment_summary
    summary["probe_status"] = {"sg": sg_status, "first_step": step_status}
    summary["artifacts"] = {
        "cell_local_coefficients": str((output / "cell_local_coefficients.csv").resolve()),
        "target_edges": str((output / "target_edges.csv").resolve()),
        "one_ring_nodes": str((output / "one_ring_nodes.csv").resolve()),
    }
    if "first_step_alignment" in summary:
        summary["artifacts"]["first_step_alignment"] = str(
            (output / "first_step_alignment.csv").resolve()
        )
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str((output / "summary.json").resolve()),
        "direct_averagebox": summary["variants"]["averagebox_current_sg"],
        "gate": summary["averagebox_first_step_implementation_gate"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
