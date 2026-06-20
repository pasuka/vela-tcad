#!/usr/bin/env python3
"""Track transition-node qF overshoot across accepted Vela continuation states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


NODE_FIELDS = [
    "step_index",
    "bias_V",
    "vtk",
    "sentaurus_bias_V",
    "sentaurus_dir",
    "node_id",
    "x_um",
    "y_um",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "vela_psi_V",
    "sentaurus_psi_V",
    "delta_psi_V",
    "vela_phin_V",
    "sentaurus_phin_V",
    "delta_phin_V",
    "abs_delta_phin_V",
    "vela_psi_minus_phin_V",
    "sentaurus_psi_minus_phin_V",
    "delta_psi_minus_phin_V",
    "vela_electron_density_cm3",
    "sentaurus_electron_density_cm3",
    "log10_vela_over_sentaurus_electron_density",
    "phin_mismatch_class",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias-min", type=float, required=True)
    parser.add_argument("--bias-max", type=float, required=True)
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--phin-threshold-v", type=float, default=0.01)
    parser.add_argument("--max-accepted-gap-v", type=float, default=0.1)
    return parser.parse_args()


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def bias_token(bias: float) -> str:
    return f"{bias:g}".replace("-", "m").replace("+", "").replace(".", "p")


def parse_sentaurus_bias(path: Path) -> float | None:
    match = re.search(r"sentaurus_(?P<bias>[-+]?\d+(?:\.\d+)?|m\d+(?:p\d+)?)v$", path.name)
    if not match:
        return None
    token = match.group("bias")
    if token.startswith("m"):
        return -float(token[1:].replace("p", "."))
    return float(token.replace("p", "."))


def discover_sentaurus_dirs(root: Path) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for path in sorted(root.iterdir()):
        if not path.is_dir() or not (path / "nodes.csv").exists():
            continue
        bias = parse_sentaurus_bias(path)
        if bias is not None:
            result.append({"bias": bias, "path": path})
    result.sort(key=lambda item: float(item["bias"]))
    return result


def discover_vtks(root: Path) -> list[dict[str, Any]]:
    pattern = re.compile(r"_(?P<step>\d+)_(?P<bias>[-+0-9.]+)V\.vtk$")
    result: list[dict[str, Any]] = []
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if not match:
            continue
        result.append({
            "step_index": int(match.group("step")),
            "bias": float(match.group("bias")),
            "path": path,
        })
    result.sort(key=lambda item: int(item["step_index"]))
    return result


def select_bias_window(states: list[dict[str, Any]], bias_min: float, bias_max: float) -> list[dict[str, Any]]:
    lo = min(bias_min, bias_max)
    hi = max(bias_min, bias_max)
    return [state for state in states if lo - 1.0e-12 <= float(state["bias"]) <= hi + 1.0e-12]


def nearest_sentaurus(sentaurus: list[dict[str, Any]], bias: float) -> dict[str, Any]:
    if not sentaurus:
        raise FileNotFoundError("no Sentaurus exports found")
    return min(sentaurus, key=lambda item: abs(float(item["bias"]) - bias))


def load_scalar(root: Path, name: str) -> dict[int, float]:
    raw = cont.load_components(root, name)
    return {node_id: values[0] for node_id, values in raw.items()}


def load_optional_scalar(root: Path, name: str) -> dict[int, float]:
    try:
        return load_scalar(root, name)
    except FileNotFoundError:
        return {}


def log10_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate <= 0.0 or reference <= 0.0:
        return None
    return math.log10(candidate / reference)


def class_for_delta(delta: float, threshold: float) -> str:
    if abs(delta) >= threshold:
        return "large"
    return "small"


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, _, _ = sgdiag.read_mesh(args.mesh)
    selected_nodes = parse_int_list(args.nodes)
    for node_id in selected_nodes:
        if node_id not in nodes:
            raise SystemExit(f"missing mesh node: {node_id}")
    states = select_bias_window(discover_vtks(args.vela_vtk_root), args.bias_min, args.bias_max)
    if not states:
        raise SystemExit("no Vela VTK states selected")
    sentaurus = discover_sentaurus_dirs(args.sentaurus_root)

    rows: list[dict[str, Any]] = []
    sent_cache: dict[Path, dict[str, Any]] = {}
    for state in states:
        scalars = sgdiag.parse_vtk_scalars(state["path"])
        for name in ("Potential", "ElectronQuasiFermi", "Electrons"):
            if name not in scalars:
                raise RuntimeError(f"VTK is missing required scalar {name}: {state['path']}")
        sent = nearest_sentaurus(sentaurus, float(state["bias"]))
        sent_path = sent["path"]
        if sent_path not in sent_cache:
            sent_nodes = cont.load_sentaurus_nodes(sent_path)
            sent_cache[sent_path] = {
                "nodes": sent_nodes,
                "psi": load_scalar(sent_path, "ElectrostaticPotential"),
                "phin": load_scalar(sent_path, "eQuasiFermiPotential"),
                "density": load_optional_scalar(sent_path, "eDensity"),
            }
        sent_data = sent_cache[sent_path]
        sent_nodes = sent_data["nodes"]
        for node_id in selected_nodes:
            node = nodes[node_id]
            sent_id, sent_dist = cont.nearest_node(sent_nodes, node["x_um"], node["y_um"])
            vela_psi = scalars["Potential"][node_id]
            vela_phin = scalars["ElectronQuasiFermi"][node_id]
            sent_psi = sent_data["psi"][sent_id]
            sent_phin = sent_data["phin"][sent_id]
            vela_density_cm3 = scalars["Electrons"][node_id] * 1.0e-6
            sent_density_cm3 = sent_data["density"].get(sent_id)
            delta_phin = vela_phin - sent_phin
            rows.append({
                "step_index": state["step_index"],
                "bias_V": state["bias"],
                "vtk": str(state["path"]),
                "sentaurus_bias_V": sent["bias"],
                "sentaurus_dir": str(sent_path),
                "node_id": node_id,
                "x_um": node["x_um"],
                "y_um": node["y_um"],
                "sentaurus_node_id": sent_id,
                "sentaurus_distance_um": sent_dist,
                "vela_psi_V": vela_psi,
                "sentaurus_psi_V": sent_psi,
                "delta_psi_V": vela_psi - sent_psi,
                "vela_phin_V": vela_phin,
                "sentaurus_phin_V": sent_phin,
                "delta_phin_V": delta_phin,
                "abs_delta_phin_V": abs(delta_phin),
                "vela_psi_minus_phin_V": vela_psi - vela_phin,
                "sentaurus_psi_minus_phin_V": sent_psi - sent_phin,
                "delta_psi_minus_phin_V": (vela_psi - vela_phin) - (sent_psi - sent_phin),
                "vela_electron_density_cm3": vela_density_cm3,
                "sentaurus_electron_density_cm3": sent_density_cm3,
                "log10_vela_over_sentaurus_electron_density": log10_ratio(vela_density_cm3, sent_density_cm3),
                "phin_mismatch_class": class_for_delta(delta_phin, args.phin_threshold_v),
            })
    return rows


def finite_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def sign_alternations(rows: list[dict[str, Any]], threshold: float) -> int:
    selected = [
        row for row in sorted(rows, key=lambda item: float(item["x_um"]))
        if abs(float(row["delta_phin_V"])) >= threshold
    ]
    signs = [1 if float(row["delta_phin_V"]) > 0.0 else -1 for row in selected]
    return sum(1 for before, after in zip(signs, signs[1:]) if before != after)


def summarize(rows: list[dict[str, Any]], max_gap_v: float, threshold_v: float) -> dict[str, Any]:
    by_step: dict[int, list[dict[str, Any]]] = {}
    for row in rows:
        by_step.setdefault(int(row["step_index"]), []).append(row)
    state_summaries: list[dict[str, Any]] = []
    for step, step_rows in sorted(by_step.items()):
        abs_phin = [abs(float(row["delta_phin_V"])) for row in step_rows]
        density_logs = [
            value for value in (finite_float(row["log10_vela_over_sentaurus_electron_density"]) for row in step_rows)
            if value is not None
        ]
        state_summaries.append({
            "step_index": step,
            "bias_V": float(step_rows[0]["bias_V"]),
            "sentaurus_bias_V": float(step_rows[0]["sentaurus_bias_V"]),
            "node_count": len(step_rows),
            "median_abs_delta_phin_V": statistics.median(abs_phin),
            "max_abs_delta_phin_V": max(abs_phin),
            "median_log10_vela_over_sentaurus_electron_density": (
                statistics.median(density_logs) if density_logs else None
            ),
            "phin_sign_alternations": sign_alternations(step_rows, threshold_v),
        })
    first_large = next(
        (item for item in state_summaries if item["max_abs_delta_phin_V"] >= threshold_v),
        None,
    )
    first_alternating_large = next(
        (
            item for item in state_summaries
            if item["max_abs_delta_phin_V"] >= threshold_v
            and int(item["phin_sign_alternations"]) > 0
        ),
        None,
    )
    gaps = [
        abs(float(after["bias_V"]) - float(before["bias_V"]))
        for before, after in zip(state_summaries, state_summaries[1:])
    ]
    largest_gap = max(gaps) if gaps else 0.0
    gap_before_first = None
    if first_large is not None:
        idx = state_summaries.index(first_large)
        if idx > 0:
            gap_before_first = abs(
                float(first_large["bias_V"]) - float(state_summaries[idx - 1]["bias_V"]))
    gap_before_first_alternating = None
    if first_alternating_large is not None:
        idx = state_summaries.index(first_alternating_large)
        if idx > 0:
            gap_before_first_alternating = abs(
                float(first_alternating_large["bias_V"]) - float(state_summaries[idx - 1]["bias_V"]))
    return {
        "state_count": len(state_summaries),
        "node_count": len({int(row["node_id"]) for row in rows}),
        "row_count": len(rows),
        "phin_threshold_V": threshold_v,
        "largest_accepted_gap_V": largest_gap,
        "gap_before_first_large_phin_mismatch_V": gap_before_first,
        "gap_before_first_alternating_large_phin_mismatch_V": gap_before_first_alternating,
        "needs_intermediate_restart": (
            (gap_before_first is not None and gap_before_first > max_gap_v)
            or (gap_before_first_alternating is not None and gap_before_first_alternating > max_gap_v)
        ),
        "first_large_phin_mismatch_bias_V": (
            first_large["bias_V"] if first_large is not None else None
        ),
        "first_large_phin_mismatch": first_large,
        "first_alternating_large_phin_mismatch_bias_V": (
            first_alternating_large["bias_V"] if first_alternating_large is not None else None
        ),
        "first_alternating_large_phin_mismatch": first_alternating_large,
        "states": state_summaries,
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
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "transition_continuation_nodes.csv", rows)
    summary = summarize(rows, args.max_accepted_gap_v, args.phin_threshold_v)
    summary.update({
        "bias_min": args.bias_min,
        "bias_max": args.bias_max,
        "nodes": parse_int_list(args.nodes),
        "max_accepted_gap_V": args.max_accepted_gap_v,
    })
    (args.out_dir / "transition_continuation_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
