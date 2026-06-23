#!/usr/bin/env python3
"""Probe whether absolute QF/psi state alignment explains PN2D BV source gaps."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


EDGE_FIELDS = [
    "bias_V",
    "edge_id",
    "relation",
    "node0",
    "node1",
    "electron_state_factor",
    "hole_state_factor",
    "flux_weighted_state_factor",
    "vela_source_density_m3_s",
    "sentaurus_generation_avg_m3_s",
    "original_log10_vela_over_sentaurus_generation",
    "state_scaled_source_density_m3_s",
    "state_scaled_log10_vela_over_sentaurus_generation",
    "source_gap_recovered_decades",
    "mean_delta_psi_minus_phin_V",
    "mean_delta_phip_minus_psi_V",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--continuity-edges", type=Path, required=True)
    parser.add_argument("--continuity-nodes", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="", help="Comma-separated bias filter; default uses all biases.")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def parse_biases(raw: str) -> set[float] | None:
    if not raw.strip():
        return None
    return {bias_key(float(item.strip())) for item in raw.split(",") if item.strip()}


def bias_key(value: float) -> float:
    return round(value, 12)


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def positive_ratio(candidate: Any, reference: Any) -> float | None:
    cand = optional_float(candidate)
    ref = optional_float(reference)
    if cand is None or ref is None or cand <= 0.0 or ref <= 0.0:
        return None
    return cand / ref


def geometric_pair(left: float | None, right: float | None) -> float | None:
    values = [value for value in (left, right) if value is not None and value > 0.0]
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    return math.sqrt(values[0] * values[1])


def log10_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate <= 0.0 or reference <= 0.0:
        return None
    return math.log10(candidate / reference)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def node_map(rows: list[dict[str, str]]) -> dict[tuple[float, int], dict[str, str]]:
    result: dict[tuple[float, int], dict[str, str]] = {}
    for row in rows:
        result[(bias_key(float(row["bias_V"])), int(row["node_id"]))] = row
    return result


def row_for_edge(edge: dict[str, str], nodes: dict[tuple[float, int], dict[str, str]]) -> dict[str, Any] | None:
    bias = bias_key(float(edge["bias_V"]))
    node0 = nodes.get((bias, int(edge["node0"])))
    node1 = nodes.get((bias, int(edge["node1"])))
    if node0 is None or node1 is None:
        return None

    electron_factor = geometric_pair(
        positive_ratio(node0.get("vela_electron_density_from_sentaurus_qf_cm3"), node0.get("vela_electron_density_cm3")),
        positive_ratio(node1.get("vela_electron_density_from_sentaurus_qf_cm3"), node1.get("vela_electron_density_cm3")),
    )
    hole_factor = geometric_pair(
        positive_ratio(node0.get("vela_hole_density_from_sentaurus_qf_cm3"), node0.get("vela_hole_density_cm3")),
        positive_ratio(node1.get("vela_hole_density_from_sentaurus_qf_cm3"), node1.get("vela_hole_density_cm3")),
    )
    if electron_factor is None and hole_factor is None:
        return None

    electron_flux = optional_float(edge.get("vela_electron_flux_abs_m2_s")) or 0.0
    hole_flux = optional_float(edge.get("vela_hole_flux_abs_m2_s")) or 0.0
    weighted_terms = []
    if electron_factor is not None and electron_flux > 0.0:
        weighted_terms.append((electron_flux, electron_factor))
    if hole_factor is not None and hole_flux > 0.0:
        weighted_terms.append((hole_flux, hole_factor))
    if weighted_terms:
        flux_weighted_factor = sum(weight * factor for weight, factor in weighted_terms) / sum(
            weight for weight, _ in weighted_terms)
    else:
        available = [factor for factor in (electron_factor, hole_factor) if factor is not None]
        flux_weighted_factor = statistics.fmean(available)

    source = optional_float(edge.get("vela_source_density_m3_s"))
    sentaurus = optional_float(edge.get("sentaurus_generation_avg_m3_s"))
    scaled_source = source * flux_weighted_factor if source is not None else None
    original_log = log10_ratio(source, sentaurus)
    scaled_log = log10_ratio(scaled_source, sentaurus)
    recovered = scaled_log - original_log if original_log is not None and scaled_log is not None else None

    delta_e = [optional_float(node.get("delta_psi_minus_phin_V")) for node in (node0, node1)]
    delta_h = [optional_float(node.get("delta_phip_minus_psi_V")) for node in (node0, node1)]
    delta_e = [value for value in delta_e if value is not None]
    delta_h = [value for value in delta_h if value is not None]

    return {
        "bias_V": float(edge["bias_V"]),
        "edge_id": int(edge["edge_id"]),
        "relation": edge.get("relation", ""),
        "node0": int(edge["node0"]),
        "node1": int(edge["node1"]),
        "electron_state_factor": electron_factor,
        "hole_state_factor": hole_factor,
        "flux_weighted_state_factor": flux_weighted_factor,
        "vela_source_density_m3_s": source,
        "sentaurus_generation_avg_m3_s": sentaurus,
        "original_log10_vela_over_sentaurus_generation": original_log,
        "state_scaled_source_density_m3_s": scaled_source,
        "state_scaled_log10_vela_over_sentaurus_generation": scaled_log,
        "source_gap_recovered_decades": recovered,
        "mean_delta_psi_minus_phin_V": statistics.fmean(delta_e) if delta_e else None,
        "mean_delta_phip_minus_psi_V": statistics.fmean(delta_h) if delta_h else None,
    }


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(bias_key(float(row["bias_V"])), []).append(row)
    bias_summaries = []
    for bias, items in sorted(groups.items()):
        recovered = [float(row["source_gap_recovered_decades"]) for row in items if row.get("source_gap_recovered_decades") is not None]
        original = [float(row["original_log10_vela_over_sentaurus_generation"]) for row in items if row.get("original_log10_vela_over_sentaurus_generation") is not None]
        scaled = [float(row["state_scaled_log10_vela_over_sentaurus_generation"]) for row in items if row.get("state_scaled_log10_vela_over_sentaurus_generation") is not None]
        active_items = [
            row for row in items
            if row.get("original_log10_vela_over_sentaurus_generation") is not None
            and float(row["original_log10_vela_over_sentaurus_generation"]) > -20.0
        ]
        active_recovered = [
            float(row["source_gap_recovered_decades"])
            for row in active_items
            if row.get("source_gap_recovered_decades") is not None
        ]
        active_original = [
            float(row["original_log10_vela_over_sentaurus_generation"])
            for row in active_items
            if row.get("original_log10_vela_over_sentaurus_generation") is not None
        ]
        active_scaled = [
            float(row["state_scaled_log10_vela_over_sentaurus_generation"])
            for row in active_items
            if row.get("state_scaled_log10_vela_over_sentaurus_generation") is not None
        ]
        focus = next((row for row in items if row.get("relation") == "focus"), None)
        bias_summaries.append({
            "bias_V": bias,
            "edge_count": len(items),
            "median_original_log10_vela_over_sentaurus_generation": statistics.median(original) if original else None,
            "median_state_scaled_log10_vela_over_sentaurus_generation": statistics.median(scaled) if scaled else None,
            "median_source_gap_recovered_decades": statistics.median(recovered) if recovered else None,
            "active_edge_count": len(active_items),
            "active_median_original_log10_vela_over_sentaurus_generation": statistics.median(active_original) if active_original else None,
            "active_median_state_scaled_log10_vela_over_sentaurus_generation": statistics.median(active_scaled) if active_scaled else None,
            "active_median_source_gap_recovered_decades": statistics.median(active_recovered) if active_recovered else None,
            "focus_edge": focus,
        })
    return {"biases": bias_summaries, "edge_rows": len(rows)}


def main() -> int:
    args = parse_args()
    selected_biases = parse_biases(args.biases)
    nodes = node_map(read_csv_rows(args.continuity_nodes))
    rows = []
    for edge in read_csv_rows(args.continuity_edges):
        bias = bias_key(float(edge["bias_V"]))
        if selected_biases is not None and bias not in selected_biases:
            continue
        item = row_for_edge(edge, nodes)
        if item is not None:
            rows.append(item)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "absolute_state_feedback_edges.csv", rows, EDGE_FIELDS)
    (args.out_dir / "absolute_state_feedback_summary.json").write_text(
        json.dumps(clean_json(summarize(rows)), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())