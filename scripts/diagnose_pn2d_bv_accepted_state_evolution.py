#!/usr/bin/env python3
"""Summarize accepted-state carrier branch movement from a VTK/IV sweep."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


STATE_FIELDS = [
    "step_index",
    "bias_V",
    "vtk",
    "current_total_A_per_um",
    "iterations",
    "accepted_step",
    "attempted_step",
    "retry_count",
    "electron_density_median_cm3",
    "electron_density_p95_cm3",
    "electron_density_max_cm3",
    "hole_density_median_cm3",
    "electron_density_jump_median_dex",
    "electron_density_jump_p95_abs_dex",
    "electron_density_jump_max_dex",
    "top_electron_jump_node_id",
    "top_electron_jump_x_um",
    "top_electron_jump_y_um",
    "top_electron_jump_before_cm3",
    "top_electron_jump_after_cm3",
    "top_electron_qf_before_V",
    "top_electron_qf_after_V",
    "top_potential_before_V",
    "top_potential_after_V",
]


EDGE_FIELDS = [
    "step_index",
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "x0_um",
    "y0_um",
    "x1_um",
    "y1_um",
    "current_total_A_per_um",
    "electron_qf_drop_V",
    "electron_qf_drop_delta_V",
    "potential_drop_V",
    "potential_drop_delta_V",
    "electron_density_avg_cm3",
    "electron_density_avg_jump_dex",
    "electron_density_node0_cm3",
    "electron_density_node1_cm3",
    "electron_density_node0_jump_dex",
    "electron_density_node1_jump_dex",
]


def main() -> int:
    args = parse_args()
    nodes, triangles, _ = sgdiag.read_mesh(args.mesh)
    edges = {int(edge["edge_id"]): edge for edge in sgdiag.build_edges(nodes, triangles)}
    edge_ids = parse_int_list(args.edge_ids)
    for edge_id in edge_ids:
        if edge_id not in edges:
            raise SystemExit(f"missing edge id: {edge_id}")
    states = select_states(discover_vtks(args.vtk_root), args.bias_min, args.bias_max)
    if not states:
        raise SystemExit("no VTK states selected")
    iv = load_iv(args.iv_csv)

    state_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    previous: dict[str, Any] | None = None
    top_jump: dict[str, Any] | None = None
    for state in states:
        scalars = sgdiag.parse_vtk_scalars(state["path"])
        require_scalars(scalars)
        iv_row = nearest_iv_row(iv, state["bias"])
        row = state_row(state, scalars, previous, nodes, iv_row)
        state_rows.append(row)
        if row.get("top_electron_jump_node_id", "") != "":
            jump_value = float(row["electron_density_jump_max_dex"])
            if top_jump is None or jump_value > float(top_jump["electron_density_jump_max_dex"]):
                top_jump = row
        for edge_id in edge_ids:
            edge_rows.append(edge_row(state, scalars, previous, edges[edge_id], nodes, iv_row))
        previous = scalars

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "accepted_state_evolution.csv", state_rows, STATE_FIELDS)
    write_csv(args.out_dir / "accepted_state_edge_evolution.csv", edge_rows, EDGE_FIELDS)
    summary = {
        "state_count": len(state_rows),
        "edge_row_count": len(edge_rows),
        "bias_min": args.bias_min,
        "bias_max": args.bias_max,
        "edge_ids": edge_ids,
        "max_electron_density_jump": summarize_top_jump(top_jump),
    }
    (args.out_dir / "accepted_state_evolution_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vtk-root", type=Path, required=True)
    parser.add_argument("--iv-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias-min", type=float, required=True)
    parser.add_argument("--bias-max", type=float, required=True)
    parser.add_argument("--edge-ids", required=True)
    return parser.parse_args()


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 10)


def discover_vtks(root: Path) -> list[dict[str, Any]]:
    pattern = re.compile(r"_(?P<step>\d+)_(?P<bias>[-+0-9.]+)V\.vtk$")
    states: list[dict[str, Any]] = []
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if not match:
            continue
        states.append({
            "step_index": int(match.group("step")),
            "bias": float(match.group("bias")),
            "path": path,
        })
    states.sort(key=lambda item: item["step_index"])
    return states


def select_states(states: list[dict[str, Any]], bias_min: float, bias_max: float) -> list[dict[str, Any]]:
    lo = min(bias_min, bias_max)
    hi = max(bias_min, bias_max)
    return [state for state in states if lo - 1.0e-12 <= float(state["bias"]) <= hi + 1.0e-12]


def load_iv(path: Path) -> list[tuple[float, dict[str, str]]]:
    rows: list[tuple[float, dict[str, str]]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if "converged" in row and str(row.get("converged", "")).strip() not in {"1", "true", "True"}:
                continue
            try:
                rows.append((float(row["bias_V"]), row))
            except (KeyError, ValueError):
                continue
    return rows


def nearest_iv_row(rows: list[tuple[float, dict[str, str]]], bias: float) -> dict[str, str] | None:
    if not rows:
        return None
    best_bias, best_row = min(rows, key=lambda item: abs(item[0] - bias))
    if abs(best_bias - bias) <= 1.0e-3:
        return best_row
    return None


def require_scalars(scalars: dict[str, list[float]]) -> None:
    missing = [name for name in ["Electrons", "Holes", "ElectronQuasiFermi", "Potential"] if name not in scalars]
    if missing:
        raise RuntimeError(f"VTK is missing required scalars: {', '.join(missing)}")


def state_row(
    state: dict[str, Any],
    scalars: dict[str, list[float]],
    previous: dict[str, list[float]] | None,
    nodes: dict[int, dict[str, float]],
    iv_row: dict[str, str] | None,
) -> dict[str, Any]:
    electrons = np.asarray(scalars["Electrons"], dtype=float)
    holes = np.asarray(scalars["Holes"], dtype=float)
    row: dict[str, Any] = {
        "step_index": state["step_index"],
        "bias_V": state["bias"],
        "vtk": str(state["path"]),
        **iv_values(iv_row),
        "electron_density_median_cm3": float(np.nanmedian(electrons)),
        "electron_density_p95_cm3": float(np.nanpercentile(electrons, 95)),
        "electron_density_max_cm3": float(np.nanmax(electrons)),
        "hole_density_median_cm3": float(np.nanmedian(holes)),
        "electron_density_jump_median_dex": "",
        "electron_density_jump_p95_abs_dex": "",
        "electron_density_jump_max_dex": "",
        "top_electron_jump_node_id": "",
    }
    if previous is None:
        return row
    jump = signed_log10_jump(np.asarray(previous["Electrons"], dtype=float), electrons)
    abs_jump = np.abs(jump)
    finite = np.isfinite(jump)
    top_index = int(np.nanargmax(abs_jump)) if np.any(finite) else -1
    row.update({
        "electron_density_jump_median_dex": float(np.nanmedian(jump[finite])),
        "electron_density_jump_p95_abs_dex": float(np.nanpercentile(abs_jump[finite], 95)),
        "electron_density_jump_max_dex": float(jump[top_index]),
        "top_electron_jump_node_id": top_index,
        "top_electron_jump_x_um": nodes[top_index]["x_um"],
        "top_electron_jump_y_um": nodes[top_index]["y_um"],
        "top_electron_jump_before_cm3": float(previous["Electrons"][top_index]),
        "top_electron_jump_after_cm3": float(electrons[top_index]),
        "top_electron_qf_before_V": float(previous["ElectronQuasiFermi"][top_index]),
        "top_electron_qf_after_V": float(scalars["ElectronQuasiFermi"][top_index]),
        "top_potential_before_V": float(previous["Potential"][top_index]),
        "top_potential_after_V": float(scalars["Potential"][top_index]),
    })
    return row


def edge_row(
    state: dict[str, Any],
    scalars: dict[str, list[float]],
    previous: dict[str, list[float]] | None,
    edge: dict[str, Any],
    nodes: dict[int, dict[str, float]],
    iv_row: dict[str, str] | None,
) -> dict[str, Any]:
    i = int(edge["node0"])
    j = int(edge["node1"])
    phin = scalars["ElectronQuasiFermi"]
    psi = scalars["Potential"]
    electrons = scalars["Electrons"]
    density_avg = 0.5 * (electrons[i] + electrons[j])
    row: dict[str, Any] = {
        "step_index": state["step_index"],
        "bias_V": state["bias"],
        "edge_id": int(edge["edge_id"]),
        "node0": i,
        "node1": j,
        "x0_um": nodes[i]["x_um"],
        "y0_um": nodes[i]["y_um"],
        "x1_um": nodes[j]["x_um"],
        "y1_um": nodes[j]["y_um"],
        "current_total_A_per_um": (iv_row or {}).get("current_total_A_per_um", ""),
        "electron_qf_drop_V": phin[j] - phin[i],
        "potential_drop_V": psi[j] - psi[i],
        "electron_density_avg_cm3": density_avg,
        "electron_density_node0_cm3": electrons[i],
        "electron_density_node1_cm3": electrons[j],
        "electron_qf_drop_delta_V": "",
        "potential_drop_delta_V": "",
        "electron_density_avg_jump_dex": "",
        "electron_density_node0_jump_dex": "",
        "electron_density_node1_jump_dex": "",
    }
    if previous is None:
        return row
    previous_drop = previous["ElectronQuasiFermi"][j] - previous["ElectronQuasiFermi"][i]
    previous_potential_drop = previous["Potential"][j] - previous["Potential"][i]
    previous_avg = 0.5 * (previous["Electrons"][i] + previous["Electrons"][j])
    row.update({
        "electron_qf_drop_delta_V": (phin[j] - phin[i]) - previous_drop,
        "potential_drop_delta_V": (psi[j] - psi[i]) - previous_potential_drop,
        "electron_density_avg_jump_dex": log10_ratio(density_avg, previous_avg),
        "electron_density_node0_jump_dex": log10_ratio(electrons[i], previous["Electrons"][i]),
        "electron_density_node1_jump_dex": log10_ratio(electrons[j], previous["Electrons"][j]),
    })
    return row


def iv_values(row: dict[str, str] | None) -> dict[str, Any]:
    if row is None:
        return {
            "current_total_A_per_um": "",
            "iterations": "",
            "accepted_step": "",
            "attempted_step": "",
            "retry_count": "",
        }
    diagnostics = parse_step_diagnostics(row.get("step_diagnostics", ""))
    return {
        "current_total_A_per_um": row.get("current_total_A_per_um", ""),
        "iterations": row.get("iterations", ""),
        "accepted_step": diagnostics.get("accepted_step", ""),
        "attempted_step": diagnostics.get("attempted_step", ""),
        "retry_count": diagnostics.get("retry_count", ""),
    }


def parse_step_diagnostics(raw: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for part in raw.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def signed_log10_jump(before: np.ndarray, after: np.ndarray) -> np.ndarray:
    eps = 1.0e-300
    result = np.full_like(after, math.nan, dtype=float)
    mask = np.isfinite(before) & np.isfinite(after)
    result[mask] = np.log10(np.maximum(after[mask], eps)) - np.log10(np.maximum(before[mask], eps))
    return result


def log10_ratio(after: float, before: float) -> float | None:
    if after <= 0.0 or before <= 0.0:
        return None
    return math.log10(after / before)


def summarize_top_jump(row: dict[str, Any] | None) -> dict[str, Any]:
    if row is None:
        return {}
    return {
        "bias_V": row["bias_V"],
        "step_index": row["step_index"],
        "jump_dex": row["electron_density_jump_max_dex"],
        "node_id": row["top_electron_jump_node_id"],
        "x_um": row.get("top_electron_jump_x_um", ""),
        "y_um": row.get("top_electron_jump_y_um", ""),
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


if __name__ == "__main__":
    raise SystemExit(main())
