#!/usr/bin/env python3
"""Merge active-overlap BV residual, state, flux, and mobility diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any


BASE_FIELDS = [
    "node_id",
    "x_m",
    "y_m",
    "support_class",
    "electron_required_source_s_inv",
    "hole_required_source_s_inv",
    "psi_V",
    "phin_V",
    "phip_V",
    "electron_density_m3",
    "hole_density_m3",
    "ni_eff_m3",
    "node_volume_m2",
    "all_combined_source_s_inv",
    "sentaurus_edgeavg_current_scaled_source_s_inv",
    "hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined",
    "hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined",
    "hybrid_e2cxx_h_edgeavg_current_residual_l2",
    "sentaurus_edgeavg_electron_current_scaled_source_s_inv",
    "sentaurus_edgeavg_hole_current_scaled_source_s_inv",
]

EDGE_FIELDS = [
    "edge_current_row_count",
    "cxx_particle_flux_m2_s_mean",
    "sentaurus_edgeavg_particle_flux_m2_s_mean",
    "sentaurus_support_particle_flux_m2_s_mean",
    "sentaurus_edgeavg_electron_current_fraction_mean",
    "sentaurus_edgeavg_hole_current_fraction_mean",
    "cxx_generation_density_m3_s_mean",
    "sentaurus_edgeavg_generation_m3_s_mean",
    "sentaurus_support_generation_m3_s_mean",
    "cxx_generation_over_sentaurus_edgeavg_generation_mean",
    "cxx_particle_flux_over_sentaurus_edgeavg_current_mean",
]

MOBILITY_FIELDS = [
    "electron_final_over_sentaurus_mobility",
    "hole_final_over_sentaurus_mobility",
    "electric_field_edge_abs_vela_over_sentaurus",
    "electron_drive_vela_over_sentaurus_qf",
    "hole_drive_vela_over_sentaurus_qf",
]

MIXED_FIELDS = [
    "mixed_electron_density_geommean_over_sentaurus",
    "mixed_hole_density_geommean_over_sentaurus",
    "mixed_particle_flux_over_sentaurus",
    "mixed_generation_over_sentaurus",
]

RESIDUAL_PROXY_FIELDS = [
    "sentaurus_state_electron_residual_over_impact",
    "sentaurus_state_hole_residual_over_impact",
    "shifted_vela_qf_electron_residual_over_impact",
    "shifted_vela_qf_hole_residual_over_impact",
]

CSV_FIELDS = BASE_FIELDS + EDGE_FIELDS + MOBILITY_FIELDS + MIXED_FIELDS + RESIDUAL_PROXY_FIELDS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-node-csv", type=Path, required=True)
    parser.add_argument("--edge-replay-csv", type=Path, required=True)
    parser.add_argument("--edge-current-csv", type=Path, required=True)
    parser.add_argument("--mobility-csv", type=Path, default=None)
    parser.add_argument("--mixed-state-csv", type=Path, default=None)
    parser.add_argument("--residual-proxy-csv", type=Path, default=None)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--support-class", default="overlap")
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def read_csv(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with open_path(path, "r", newline="") as handle:
        return list(csv.DictReader(handle))


def finite_float(raw: Any, default: float = 0.0) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if denominator == 0.0:
        return None
    return numerator / denominator


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None:
        return True
    if row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def nonimpact_terms(row: dict[str, str], carrier: str) -> float:
    return (
        finite_float(row.get(f"{carrier}_flux"))
        + finite_float(row.get(f"{carrier}_recombination"))
        + finite_float(row.get(f"{carrier}_gauge"))
        + finite_float(row.get(f"{carrier}_boundary"))
    )


def rows_by_node(rows: list[dict[str, str]], variant: str | None, bias: float | None) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        if variant is not None and row.get("variant") not in (None, "", variant):
            continue
        if not matches_bias(row, bias):
            continue
        node_id = row.get("node_id")
        if node_id not in (None, ""):
            result[node_id] = row
    return result


def mean(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def aggregate_edge_current(rows: list[dict[str, str]], bias: float | None) -> dict[str, dict[str, Any]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        if not matches_bias(row, bias):
            continue
        node_id = row.get("support_node_id")
        if node_id not in (None, ""):
            grouped.setdefault(node_id, []).append(row)

    result: dict[str, dict[str, Any]] = {}
    for node_id, items in grouped.items():
        edge_particle = [finite_float(row.get("sentaurus_edgeavg_particle_flux_m2_s")) for row in items]
        electron_fractions: list[float] = []
        hole_fractions: list[float] = []
        for row in items:
            particle = finite_float(row.get("sentaurus_edgeavg_particle_flux_m2_s"))
            electron_fraction = safe_ratio(
                finite_float(row.get("sentaurus_edgeavg_electron_current_flux_m2_s")),
                particle,
            )
            hole_fraction = safe_ratio(
                finite_float(row.get("sentaurus_edgeavg_hole_current_flux_m2_s")),
                particle,
            )
            if electron_fraction is not None:
                electron_fractions.append(electron_fraction)
            if hole_fraction is not None:
                hole_fractions.append(hole_fraction)
        result[node_id] = {
            "edge_current_row_count": len(items),
            "cxx_particle_flux_m2_s_mean": mean([
                finite_float(row.get("cxx_particle_flux_m2_s")) for row in items
            ]),
            "sentaurus_edgeavg_particle_flux_m2_s_mean": mean(edge_particle),
            "sentaurus_support_particle_flux_m2_s_mean": mean([
                finite_float(row.get("sentaurus_support_particle_flux_m2_s")) for row in items
            ]),
            "sentaurus_edgeavg_electron_current_fraction_mean": mean(electron_fractions),
            "sentaurus_edgeavg_hole_current_fraction_mean": mean(hole_fractions),
            "cxx_generation_density_m3_s_mean": mean([
                finite_float(row.get("cxx_generation_density_m3_s")) for row in items
            ]),
            "sentaurus_edgeavg_generation_m3_s_mean": mean([
                finite_float(row.get("sentaurus_edgeavg_generation_m3_s")) for row in items
            ]),
            "sentaurus_support_generation_m3_s_mean": mean([
                finite_float(row.get("sentaurus_support_generation_m3_s")) for row in items
            ]),
            "cxx_generation_over_sentaurus_edgeavg_generation_mean": mean([
                finite_float(row.get("cxx_generation_over_sentaurus_edgeavg_generation")) for row in items
            ]),
            "cxx_particle_flux_over_sentaurus_edgeavg_current_mean": mean([
                finite_float(row.get("cxx_particle_flux_over_sentaurus_edgeavg_current")) for row in items
            ]),
        }
    return result


def mobility_by_node(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {
        row["node_id"]: row
        for row in rows
        if row.get("node_id") not in (None, "")
    }


def residual_proxy_by_node(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    variants = {
        "sentaurus_state_sentaurus_node_source": "sentaurus_state",
        "shifted_vela_qf_sentaurus_node_source": "shifted_vela_qf",
    }
    for row in rows:
        node_id = row.get("node_id")
        prefix = variants.get(row.get("variant", ""))
        if node_id in (None, "") or prefix is None:
            continue
        item = result.setdefault(node_id, {})
        item[f"{prefix}_electron_residual_over_impact"] = optional_float(
            row.get("electron_residual_over_impact")
        )
        item[f"{prefix}_hole_residual_over_impact"] = optional_float(
            row.get("hole_residual_over_impact")
        )
    return result


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    replay = rows_by_node(read_csv(args.edge_replay_csv), args.variant, args.bias)
    edge_current = aggregate_edge_current(read_csv(args.edge_current_csv), args.bias)
    mobility = mobility_by_node(read_csv(args.mobility_csv))
    mixed = rows_by_node(read_csv(args.mixed_state_csv), args.variant, None)
    residual_proxy = residual_proxy_by_node(read_csv(args.residual_proxy_csv))

    rows: list[dict[str, Any]] = []
    for exact in read_csv(args.exact_node_csv):
        if exact.get("variant") != args.variant:
            continue
        if not matches_bias(exact, args.bias):
            continue
        if exact.get("support_class") != args.support_class:
            continue
        node_id = exact.get("node_id", "")
        replay_row = replay.get(node_id, {})
        mobility_row = mobility.get(node_id, {})
        mixed_row = mixed.get(node_id, {})
        item: dict[str, Any] = {
            "node_id": node_id,
            "x_m": exact.get("x", exact.get("x_um", "")),
            "y_m": exact.get("y", exact.get("y_um", "")),
            "support_class": exact.get("support_class", ""),
            "electron_required_source_s_inv": nonimpact_terms(exact, "electron"),
            "hole_required_source_s_inv": nonimpact_terms(exact, "hole"),
            "psi_V": exact.get("psi_V", ""),
            "phin_V": exact.get("phin_V", ""),
            "phip_V": exact.get("phip_V", ""),
            "electron_density_m3": exact.get("electron_density_m3", ""),
            "hole_density_m3": exact.get("hole_density_m3", ""),
            "ni_eff_m3": exact.get("ni_eff_m3", ""),
            "node_volume_m2": exact.get("node_volume_m2", ""),
            "all_combined_source_s_inv": replay_row.get("all_combined_source", ""),
            "sentaurus_edgeavg_current_scaled_source_s_inv": replay_row.get(
                "sentaurus_edgeavg_current_scaled_source", ""
            ),
            "hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined": replay_row.get(
                "hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined", ""
            ),
            "hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined": replay_row.get(
                "hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined", ""
            ),
            "hybrid_e2cxx_h_edgeavg_current_residual_l2": replay_row.get(
                "hybrid_e2cxx_h_edgeavg_current_residual_l2", ""
            ),
            "sentaurus_edgeavg_electron_current_scaled_source_s_inv": replay_row.get(
                "sentaurus_edgeavg_electron_current_scaled_source", ""
            ),
            "sentaurus_edgeavg_hole_current_scaled_source_s_inv": replay_row.get(
                "sentaurus_edgeavg_hole_current_scaled_source", ""
            ),
        }
        item.update(edge_current.get(node_id, {}))
        for key in MOBILITY_FIELDS:
            item[key] = mobility_row.get(key, "")
        item.update({
            "mixed_electron_density_geommean_over_sentaurus": mixed_row.get(
                "electron_density_geommean_over_sentaurus", ""
            ),
            "mixed_hole_density_geommean_over_sentaurus": mixed_row.get(
                "hole_density_geommean_over_sentaurus", ""
            ),
            "mixed_particle_flux_over_sentaurus": mixed_row.get(
                "particle_flux_over_sentaurus", ""
            ),
            "mixed_generation_over_sentaurus": mixed_row.get("generation_over_sentaurus", ""),
        })
        item.update(residual_proxy.get(node_id, {}))
        rows.append(item)
    rows.sort(key=lambda row: finite_float(row.get("y_m")))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    with open_path(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def summary(rows: list[dict[str, Any]], args: argparse.Namespace) -> dict[str, Any]:
    max_row = max(
        rows,
        key=lambda row: finite_float(row.get("hybrid_e2cxx_h_edgeavg_current_residual_l2")),
        default={},
    )
    y_values = [finite_float(row.get("y_m")) for row in rows]
    return clean_json({
        "variant": args.variant,
        "bias_V": args.bias,
        "support_class": args.support_class,
        "row_count": len(rows),
        "y_min_m": min(y_values) if y_values else None,
        "y_max_m": max(y_values) if y_values else None,
        "max_hybrid_residual_node_id": max_row.get("node_id", ""),
        "max_hybrid_residual_l2": optional_float(
            max_row.get("hybrid_e2cxx_h_edgeavg_current_residual_l2")
        ),
    })


def main() -> int:
    args = parse_args()
    rows = build_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_overlap_outlier_table.csv", rows)
    payload = summary(rows, args)
    write_text(
        args.out_dir / "active_overlap_outlier_summary.json",
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
