#!/usr/bin/env python3
"""Decompose electron continuity row terms on selected BV transition edges."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_flux_forms as fluxforms
import diagnose_pn2d_bv_transition_flux_replay as fluxreplay


EDGE_FIELDS = [
    "bias_V",
    "state",
    "start_node",
    "target_contact",
    "step",
    "edge_id",
    "node_from",
    "node_to",
    "edge_mid_x_um",
    "length_m",
    "couple_m",
    "electron_flux_integral_s_inv",
    "d_flux_from_s_inv_per_V",
    "d_flux_to_s_inv_per_V",
]


NODE_FIELDS = [
    "bias_V",
    "state",
    "node_id",
    "x_um",
    "y_um",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "selected_edge_count",
    "selected_edge_ids",
    "node_volume_m2",
    "vela_phin_V",
    "sentaurus_phin_V",
    "state_phin_V",
    "vela_to_sentaurus_phin_delta_V",
    "state_phin_minus_sentaurus_V",
    "electron_residual_contribution_s_inv",
    "sentaurus_electron_residual_contribution_s_inv",
    "electron_residual_delta_vs_sentaurus_s_inv",
    "d_residual_d_self_phin_s_inv_per_V",
    "d_residual_d_neighbor_phin_sum_s_inv_per_V",
    "local_newton_step_phin_V",
    "local_delta_newton_step_phin_V",
    "local_delta_step_alignment_with_sentaurus",
    "vela_avalanche_node_integral_s_inv",
    "vela_srh_node_integral_s_inv",
    "sentaurus_generation_node_integral_s_inv",
    "sentaurus_srh_node_integral_s_inv",
    "electron_full_residual_s_inv",
    "sentaurus_electron_full_residual_s_inv",
    "electron_full_residual_delta_vs_sentaurus_s_inv",
    "source_delta_vs_sentaurus_s_inv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--start-node", type=int)
    parser.add_argument("--target-contact", default="")
    parser.add_argument("--edge-ids", default="")
    parser.add_argument("--x-min-um", type=float)
    parser.add_argument("--x-max-um", type=float)
    parser.add_argument("--delta-v", type=float, default=1.0e-4)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def avg2(a: float, b: float) -> float:
    return 0.5 * (a + b)


def safe_div(numerator: float, denominator: float) -> float | None:
    if denominator == 0.0 or not math.isfinite(denominator):
        return None
    return numerator / denominator


def sign_alignment(step: float | None, target: float) -> str:
    if step is None or step == 0.0 or target == 0.0:
        return "neutral"
    return "toward_sentaurus" if step * target > 0.0 else "away_from_sentaurus"


def load_optional_sentaurus_scalar(root: Path, name: str) -> dict[int, float]:
    try:
        return fluxreplay.load_sentaurus_scalar(root, name)
    except FileNotFoundError:
        return {}


def sentaurus_scalar_or_zero(values: dict[int, float], node_id: int) -> float:
    return values.get(node_id, 0.0)


def add_node_term(
    terms: dict[tuple[str, int], dict[str, Any]],
    state: str,
    node_id: int,
    edge_id: int,
    contribution: float,
    d_self: float,
    d_neighbor: float,
) -> None:
    key = (state, node_id)
    row = terms.setdefault(key, {
        "state": state,
        "node_id": node_id,
        "selected_edge_ids": set(),
        "electron_residual_contribution_s_inv": 0.0,
        "d_residual_d_self_phin_s_inv_per_V": 0.0,
        "d_residual_d_neighbor_phin_sum_s_inv_per_V": 0.0,
    })
    row["selected_edge_ids"].add(edge_id)
    row["electron_residual_contribution_s_inv"] += contribution
    row["d_residual_d_self_phin_s_inv_per_V"] += d_self
    row["d_residual_d_neighbor_phin_sum_s_inv_per_V"] += d_neighbor


def make_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    nodes, triangles, _ = cont.sgdiag.read_mesh(args.mesh)
    volumes = cont.sgdiag.node_volumes(nodes, triangles)
    vt = cont.sgdiag.K_B_OVER_Q * args.temperature_k
    vtk_path = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    scalars = cont.sgdiag.parse_vtk_scalars(vtk_path)
    required_len = max(nodes) + 1 if nodes else 0
    fluxreplay.require_vtk_scalars(
        scalars,
        ("Potential", "ElectronQuasiFermi", "ElectronMobility"),
        required_len,
        vtk_path,
    )
    for optional_name in ("AvalancheGeneration", "SRHRecombination"):
        if optional_name in scalars:
            fluxreplay.require_vtk_scalars(scalars, (optional_name,), required_len, vtk_path)
    sentaurus_nodes = cont.load_sentaurus_nodes(args.sentaurus_dir)
    sent_nearest = fluxreplay.nearest_cache(sentaurus_nodes, nodes)
    sent_psi = fluxreplay.load_sentaurus_scalar(args.sentaurus_dir, "ElectrostaticPotential")
    sent_phin = fluxreplay.load_sentaurus_scalar(args.sentaurus_dir, "eQuasiFermiPotential")
    sent_generation = load_optional_sentaurus_scalar(args.sentaurus_dir, "ImpactIonization")
    sent_srh = load_optional_sentaurus_scalar(args.sentaurus_dir, "srhRecombination")
    vela_avalanche = scalars.get("AvalancheGeneration", [0.0 for _ in nodes])
    vela_srh = scalars.get("SRHRecombination", [0.0 for _ in nodes])
    ni = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    states = fluxreplay.make_state_vectors(nodes, scalars, sent_nearest, sent_psi, sent_phin)
    selected = fluxreplay.selected_path_rows(args)
    fluxreplay.require_selected_path_rows(selected)

    edge_rows: list[dict[str, Any]] = []
    node_terms: dict[tuple[str, int], dict[str, Any]] = {}
    selected_nodes: set[int] = set()
    for path_row in selected:
        i = int(path_row["node_from"])
        j = int(path_row["node_to"])
        selected_nodes.update((i, j))
        edge_id = int(path_row["edge_id"])
        length_m = float(path_row.get("length_m") or 0.0)
        if length_m <= 0.0:
            length_m = math.hypot(
                (nodes[j]["x_um"] - nodes[i]["x_um"]) * 1.0e-6,
                (nodes[j]["y_um"] - nodes[i]["y_um"]) * 1.0e-6,
            )
        couple_m = float(path_row.get("couple_m") or 1.0)
        mobility = avg2(scalars["ElectronMobility"][i], scalars["ElectronMobility"][j])
        for state_name, state in states.items():
            psi = state["psi"]
            phin = state["phin"]
            assert isinstance(psi, list)
            assert isinstance(phin, list)
            flux = fluxreplay.flux_value(i, j, length_m, mobility, psi, phin, ni, vt)
            integral = flux * couple_m
            d_from = fluxreplay.derivative_at(
                i, i, j, length_m, couple_m, mobility, psi, phin, ni, vt, args.delta_v)
            d_to = fluxreplay.derivative_at(
                j, i, j, length_m, couple_m, mobility, psi, phin, ni, vt, args.delta_v)
            edge_rows.append({
                "bias_V": args.bias,
                "state": state_name,
                "start_node": path_row.get("start_node", ""),
                "target_contact": path_row.get("target_contact", ""),
                "step": path_row.get("step", ""),
                "edge_id": edge_id,
                "node_from": i,
                "node_to": j,
                "edge_mid_x_um": avg2(nodes[i]["x_um"], nodes[j]["x_um"]),
                "length_m": length_m,
                "couple_m": couple_m,
                "electron_flux_integral_s_inv": integral,
                "d_flux_from_s_inv_per_V": d_from,
                "d_flux_to_s_inv_per_V": d_to,
            })
            add_node_term(node_terms, state_name, i, edge_id, integral, d_from, d_to)
            add_node_term(node_terms, state_name, j, edge_id, -integral, -d_to, -d_from)

    sentaurus_by_node = {
        node_id: node_terms.get(("sentaurus", node_id), {}).get("electron_residual_contribution_s_inv", 0.0)
        for node_id in selected_nodes
    }
    sentaurus_full_by_node: dict[int, float] = {}
    source_terms_by_node: dict[int, dict[str, float]] = {}
    for node_id in selected_nodes:
        sent_id, _ = sent_nearest[node_id]
        volume = volumes[node_id]
        vela_avalanche_integral = vela_avalanche[node_id] * volume
        vela_srh_integral = vela_srh[node_id] * volume
        sent_generation_integral = sentaurus_scalar_or_zero(sent_generation, sent_id) * 1.0e6 * volume
        sent_srh_integral = sentaurus_scalar_or_zero(sent_srh, sent_id) * 1.0e6 * volume
        source_terms_by_node[node_id] = {
            "vela_avalanche": vela_avalanche_integral,
            "vela_srh": vela_srh_integral,
            "sentaurus_generation": sent_generation_integral,
            "sentaurus_srh": sent_srh_integral,
        }
        sentaurus_full_by_node[node_id] = (
            sentaurus_by_node.get(node_id, 0.0) + sent_srh_integral - sent_generation_integral)
    node_rows: list[dict[str, Any]] = []
    for (state_name, node_id), term in sorted(node_terms.items(), key=lambda item: (item[0][1], item[0][0])):
        sent_id, sent_dist = sent_nearest[node_id]
        state = states[state_name]
        phin = state["phin"]
        assert isinstance(phin, list)
        residual = float(term["electron_residual_contribution_s_inv"])
        delta = residual - sentaurus_by_node.get(node_id, 0.0)
        d_self = float(term["d_residual_d_self_phin_s_inv_per_V"])
        local_step = safe_div(-residual, d_self)
        local_delta_step = safe_div(-delta, d_self)
        target_delta = sent_phin[sent_id] - scalars["ElectronQuasiFermi"][node_id]
        edge_ids = sorted(term["selected_edge_ids"])
        source_terms = source_terms_by_node[node_id]
        if state_name == "sentaurus":
            full_residual = sentaurus_full_by_node[node_id]
        else:
            full_residual = residual + source_terms["vela_srh"] - source_terms["vela_avalanche"]
        full_delta = full_residual - sentaurus_full_by_node[node_id]
        source_delta = (
            source_terms["vela_srh"] - source_terms["vela_avalanche"]
            - source_terms["sentaurus_srh"] + source_terms["sentaurus_generation"])
        node_rows.append({
            "bias_V": args.bias,
            "state": state_name,
            "node_id": node_id,
            "x_um": nodes[node_id]["x_um"],
            "y_um": nodes[node_id]["y_um"],
            "sentaurus_node_id": sent_id,
            "sentaurus_distance_um": sent_dist,
            "selected_edge_count": len(edge_ids),
            "selected_edge_ids": ";".join(str(edge_id) for edge_id in edge_ids),
            "node_volume_m2": volumes[node_id],
            "vela_phin_V": scalars["ElectronQuasiFermi"][node_id],
            "sentaurus_phin_V": sent_phin[sent_id],
            "state_phin_V": phin[node_id],
            "vela_to_sentaurus_phin_delta_V": target_delta,
            "state_phin_minus_sentaurus_V": phin[node_id] - sent_phin[sent_id],
            "electron_residual_contribution_s_inv": residual,
            "sentaurus_electron_residual_contribution_s_inv": sentaurus_by_node.get(node_id, 0.0),
            "electron_residual_delta_vs_sentaurus_s_inv": delta,
            "d_residual_d_self_phin_s_inv_per_V": d_self,
            "d_residual_d_neighbor_phin_sum_s_inv_per_V": term["d_residual_d_neighbor_phin_sum_s_inv_per_V"],
            "local_newton_step_phin_V": local_step,
            "local_delta_newton_step_phin_V": local_delta_step,
            "local_delta_step_alignment_with_sentaurus": sign_alignment(local_delta_step, target_delta),
            "vela_avalanche_node_integral_s_inv": source_terms["vela_avalanche"],
            "vela_srh_node_integral_s_inv": source_terms["vela_srh"],
            "sentaurus_generation_node_integral_s_inv": source_terms["sentaurus_generation"],
            "sentaurus_srh_node_integral_s_inv": source_terms["sentaurus_srh"],
            "electron_full_residual_s_inv": full_residual,
            "sentaurus_electron_full_residual_s_inv": sentaurus_full_by_node[node_id],
            "electron_full_residual_delta_vs_sentaurus_s_inv": full_delta,
            "source_delta_vs_sentaurus_s_inv": source_delta,
        })
    return edge_rows, node_rows


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def summarize(node_rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_state: dict[str, list[float]] = defaultdict(list)
    full_by_state: dict[str, list[float]] = defaultdict(list)
    source_delta_values: list[float] = []
    alignments: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in node_rows:
        state = str(row["state"])
        if state == "sentaurus":
            continue
        value = optional_float(row.get("electron_residual_delta_vs_sentaurus_s_inv"))
        if value is not None:
            by_state[state].append(abs(value))
        full_value = optional_float(row.get("electron_full_residual_delta_vs_sentaurus_s_inv"))
        if full_value is not None:
            full_by_state[state].append(abs(full_value))
        source_delta = optional_float(row.get("source_delta_vs_sentaurus_s_inv"))
        if source_delta is not None:
            source_delta_values.append(abs(source_delta))
        alignments[state][str(row.get("local_delta_step_alignment_with_sentaurus", "neutral"))] += 1
    median_by_state = {
        state: statistics.median(values)
        for state, values in by_state.items()
        if values
    }
    best = min(median_by_state, key=median_by_state.get) if median_by_state else ""
    full_median_by_state = {
        state: statistics.median(values)
        for state, values in full_by_state.items()
        if values
    }
    full_best = min(full_median_by_state, key=full_median_by_state.get) if full_median_by_state else ""
    return {
        "node_count": len({int(row["node_id"]) for row in node_rows}),
        "row_count": len(node_rows),
        "median_abs_residual_delta_vs_sentaurus_s_inv": median_by_state,
        "best_residual_match_state": best,
        "median_abs_full_residual_delta_vs_sentaurus_s_inv": full_median_by_state,
        "best_full_residual_match_state": full_best,
        "median_abs_source_delta_vs_sentaurus_s_inv": (
            statistics.median(source_delta_values) if source_delta_values else None
        ),
        "local_delta_step_alignment_counts": {
            state: dict(counts) for state, counts in alignments.items()
        },
    }


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    args = parse_args()
    edge_rows, node_rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "transition_row_residual_edges.csv", edge_rows, EDGE_FIELDS)
    write_csv(args.out_dir / "transition_row_residual_nodes.csv", node_rows, NODE_FIELDS)
    summary = summarize(node_rows)
    summary.update({
        "bias_V": args.bias,
        "path_edge_csv": str(args.path_edge_csv),
        "x_min_um": args.x_min_um,
        "x_max_um": args.x_max_um,
        "target_contact": args.target_contact,
        "start_node": args.start_node,
    })
    (args.out_dir / "transition_row_residual_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
