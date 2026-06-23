#!/usr/bin/env python3
"""Prepare a coupled QF-gradient PN2D BV predictor restart state.

This diagnostic helper starts from a converged Vela state, selects active SG
edges around requested support classes, blends their endpoint electron/hole QF
potentials toward the Sentaurus endpoint pattern, and reconstructs carriers from
Vela inferred-ni. It intentionally leaves psi unchanged so the experiment
isolates carrier-continuity/QF-gradient branch movement.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_sg_flux_forms as fluxforms

STATE_FIELDS = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
NODE_FIELDS = [
    "node_id",
    "support_class",
    "target_weight",
    "active_edge_ids",
    "phin_baseline_V",
    "phin_predictor_V",
    "phin_sentaurus_V",
    "phip_baseline_V",
    "phip_predictor_V",
    "phip_sentaurus_V",
    "electron_density_baseline_m3",
    "electron_density_predictor_m3",
    "hole_density_baseline_m3",
    "hole_density_predictor_m3",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--state-csv-name", default="coupled_qf_predictor_state.csv")
    parser.add_argument("--support-class", action="append", dest="support_classes")
    parser.add_argument("--active-axis", choices=["x", "y"], default="x")
    parser.add_argument("--active-source-relative-threshold", type=float, default=0.0)
    parser.add_argument("--blend", type=float, default=1.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def selected_support(path: Path, selected_classes: set[str]) -> dict[int, str]:
    selected: dict[int, str] = {}
    for row in read_csv_rows(path):
        support_class = row.get("support_class", "")
        if support_class in selected_classes:
            selected[int(row["node_id"])] = support_class
    return selected


def parse_support_classes(raw: list[str] | None) -> set[str]:
    if not raw:
        return {"false_negative"}
    classes: set[str] = set()
    for item in raw:
        classes.update(part.strip() for part in item.split(",") if part.strip())
    if not classes:
        raise SystemExit("--support-class must select at least one support class")
    return classes


def inferred_ni(density: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if density <= 0.0:
        return 0.0
    if carrier == "electron":
        return density / fluxforms.sgdiag.limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return density / fluxforms.sgdiag.limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def density_from_qf(ni: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if ni <= 0.0:
        return 0.0
    if carrier == "electron":
        return ni * fluxforms.sgdiag.limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return ni * fluxforms.sgdiag.limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def collect_active_edges(
    support_nodes: dict[int, str],
    active_by_node: dict[int, list[dict[str, Any]]],
    mesh_edges: dict[int, dict[str, Any]],
) -> tuple[dict[int, dict[str, Any]], dict[int, set[int]]]:
    edges: dict[int, dict[str, Any]] = {}
    node_edges: dict[int, set[int]] = {}
    for support_node in support_nodes:
        for item in active_by_node.get(support_node, []):
            edge_id = int(item["edge_id"])
            mesh_edge = mesh_edges.get(edge_id)
            if mesh_edge is None:
                continue
            merged = {**item, **mesh_edge}
            edges[edge_id] = merged
            for endpoint in (int(mesh_edge["node0"]), int(mesh_edge["node1"])):
                node_edges.setdefault(endpoint, set()).add(edge_id)
    return edges, node_edges


def build_predictor(
    vela: dict[str, list[float]],
    sentaurus: dict[str, list[float]],
    target_nodes: set[int],
    blend: float,
    vt: float,
) -> dict[str, list[float]]:
    predictor = {key: list(values) for key, values in vela.items()}
    for node_id in target_nodes:
        psi = vela["psi"][node_id]
        ni_e = inferred_ni(vela["n"][node_id], psi, vela["phin"][node_id], vt, "electron")
        ni_h = inferred_ni(vela["p"][node_id], psi, vela["phip"][node_id], vt, "hole")
        phin = vela["phin"][node_id] + blend * (sentaurus["phin"][node_id] - vela["phin"][node_id])
        phip = vela["phip"][node_id] + blend * (sentaurus["phip"][node_id] - vela["phip"][node_id])
        predictor["phin"][node_id] = phin
        predictor["phip"][node_id] = phip
        predictor["n"][node_id] = density_from_qf(ni_e, psi, phin, vt, "electron")
        predictor["p"][node_id] = density_from_qf(ni_h, psi, phip, vt, "hole")
    return predictor


def write_state(path: Path, state: dict[str, list[float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATE_FIELDS)
        writer.writeheader()
        for node_id in range(len(state["psi"])):
            writer.writerow({
                "node_id": node_id,
                "psi": state["psi"][node_id],
                "phin": state["phin"][node_id],
                "phip": state["phip"][node_id],
                "electrons_m3": state["n"][node_id],
                "holes_m3": state["p"][node_id],
            })


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not math.isfinite(args.blend):
        raise SystemExit("--blend must be finite")
    support_classes = parse_support_classes(args.support_classes)
    nodes, triangles, _ = fluxforms.sgdiag.read_mesh(args.mesh)
    node_count = len(nodes)
    mesh_edges = {int(edge["edge_id"]): edge for edge in fluxforms.sgdiag.build_edges(nodes, triangles)}
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    vela = fluxforms.load_vela_state(fluxfac.discover_vtk(args.vela_vtk_root, args.bias), node_count)
    sentaurus = fluxforms.load_sentaurus_state(args.sentaurus_dir, node_count)
    support_nodes = selected_support(args.support_csv, support_classes)
    active_by_node = fluxfac.load_active_cxx_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    active_edges, node_edges = collect_active_edges(support_nodes, active_by_node, mesh_edges)
    target_nodes = set(node_edges)
    predictor = build_predictor(vela, sentaurus, target_nodes, args.blend, vt)

    edge_gradient_n: list[float] = []
    edge_gradient_p: list[float] = []
    for edge in active_edges.values():
        i = int(edge["node0"])
        j = int(edge["node1"])
        base_n = vela["phin"][j] - vela["phin"][i]
        pred_n = predictor["phin"][j] - predictor["phin"][i]
        base_p = vela["phip"][j] - vela["phip"][i]
        pred_p = predictor["phip"][j] - predictor["phip"][i]
        edge_gradient_n.append(pred_n - base_n)
        edge_gradient_p.append(pred_p - base_p)

    node_rows: list[dict[str, Any]] = []
    for node_id in sorted(target_nodes):
        edge_ids = sorted(node_edges.get(node_id, set()))
        node_rows.append({
            "node_id": node_id,
            "support_class": support_nodes.get(node_id) or "edge_endpoint",
            "target_weight": 1.0,
            "active_edge_ids": ";".join(str(edge_id) for edge_id in edge_ids),
            "phin_baseline_V": vela["phin"][node_id],
            "phin_predictor_V": predictor["phin"][node_id],
            "phin_sentaurus_V": sentaurus["phin"][node_id],
            "phip_baseline_V": vela["phip"][node_id],
            "phip_predictor_V": predictor["phip"][node_id],
            "phip_sentaurus_V": sentaurus["phip"][node_id],
            "electron_density_baseline_m3": vela["n"][node_id],
            "electron_density_predictor_m3": predictor["n"][node_id],
            "hole_density_baseline_m3": vela["p"][node_id],
            "hole_density_predictor_m3": predictor["p"][node_id],
        })

    summary = {
        "bias_V": args.bias,
        "blend": args.blend,
        "support_classes": sorted(support_classes),
        "selected_support_count": len(support_nodes),
        "active_edge_count": len(active_edges),
        "target_node_count": len(target_nodes),
        "electron_qf_gradient_delta_median_V": median(edge_gradient_n),
        "hole_qf_gradient_delta_median_V": median(edge_gradient_p),
        "electron_density_factor_median": median([
            predictor["n"][node_id] / vela["n"][node_id]
            for node_id in target_nodes if vela["n"][node_id] != 0.0
        ]),
        "hole_density_factor_median": median([
            predictor["p"][node_id] / vela["p"][node_id]
            for node_id in target_nodes if vela["p"][node_id] != 0.0
        ]),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_state(args.out_dir / args.state_csv_name, predictor)
    write_csv(args.out_dir / "coupled_qf_predictor_nodes.csv", node_rows, NODE_FIELDS)
    (args.out_dir / "coupled_qf_predictor_state_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
