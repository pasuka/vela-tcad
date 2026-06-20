#!/usr/bin/env python3
"""Compare residual-proxy continuity units against exact carrier-row units."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any


Q = 1.602176634e-19

CSV_FIELDS = [
    "node_id",
    "x_m",
    "y_m",
    "support_class",
    "node_volume_m2",
    "inferred_continuity_scale_from_srh",
    "inverse_inferred_continuity_scale_from_srh",
    "exact_all_combined_source",
    "exact_sentaurus_edgeavg_current_scaled_source",
    "exact_electron_required_source",
    "exact_hole_required_source",
    "proxy_variant",
    "transport_model",
    "proxy_impact_source_s_inv",
    "proxy_electron_required_s_inv",
    "proxy_hole_required_s_inv",
    "proxy_electron_residual_over_impact",
    "proxy_hole_residual_over_impact",
    "exact_edgeavg_current_over_proxy_impact_scale",
    "exact_all_combined_over_proxy_impact_scale",
    "exact_electron_required_over_proxy_required_scale",
    "exact_hole_required_over_proxy_required_scale",
    "electron_required_scale_over_edgeavg_scale",
    "hole_required_scale_over_edgeavg_scale",
    "edgeavg_scale_over_inverse_srh_scale",
    "edgeavg_scale_over_node_volume",
    "edgeavg_scale_over_q_node_volume",
]

SUMMARY_METRICS = [
    "exact_edgeavg_current_over_proxy_impact_scale",
    "exact_all_combined_over_proxy_impact_scale",
    "exact_electron_required_over_proxy_required_scale",
    "exact_hole_required_over_proxy_required_scale",
    "electron_required_scale_over_edgeavg_scale",
    "hole_required_scale_over_edgeavg_scale",
    "edgeavg_scale_over_inverse_srh_scale",
    "edgeavg_scale_over_node_volume",
    "edgeavg_scale_over_q_node_volume",
    "proxy_electron_residual_over_impact",
    "proxy_hole_residual_over_impact",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-node-csv", type=Path, required=True)
    parser.add_argument("--edge-replay-csv", type=Path, required=True)
    parser.add_argument("--residual-proxy-csv", type=Path, required=True)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument(
        "--proxy-variant",
        default="probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source",
    )
    parser.add_argument("--transport-model", default="qf_old_slotboom_ni")
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


def read_csv(path: Path) -> list[dict[str, str]]:
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


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator in (None, 0.0):
        return None
    return numerator / denominator


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None or row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def nonimpact_terms(row: dict[str, str], carrier: str) -> float:
    return (
        finite_float(row.get(f"{carrier}_flux"))
        + finite_float(row.get(f"{carrier}_recombination"))
        + finite_float(row.get(f"{carrier}_gauge"))
        + finite_float(row.get(f"{carrier}_boundary"))
    )


def replay_by_node(rows: list[dict[str, str]], variant: str, bias: float | None) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        if row.get("variant") not in (None, "", variant):
            continue
        if not matches_bias(row, bias):
            continue
        node_id = row.get("node_id")
        if node_id not in (None, ""):
            result[node_id] = row
    return result


def proxy_by_node(
    rows: list[dict[str, str]],
    variant: str,
    transport_model: str,
    support_class: str,
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        if row.get("variant") != variant:
            continue
        if row.get("transport_model") != transport_model:
            continue
        if row.get("support_class") not in (None, "", support_class):
            continue
        node_id = row.get("node_id")
        if node_id not in (None, ""):
            result[node_id] = row
    return result


def required_from_proxy(row: dict[str, str], carrier: str) -> float:
    transport = finite_float(row.get(f"{carrier}_transport_s_inv"))
    srh = finite_float(row.get("srh_source_s_inv"))
    explicit = optional_float(row.get(f"{carrier}_residual_proxy_s_inv"))
    impact = finite_float(row.get("impact_source_s_inv"))
    return impact + explicit if explicit is not None else transport + srh


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    replay_rows = replay_by_node(read_csv(args.edge_replay_csv), args.variant, args.bias)
    proxy_rows = proxy_by_node(
        read_csv(args.residual_proxy_csv),
        args.proxy_variant,
        args.transport_model,
        args.support_class,
    )
    rows: list[dict[str, Any]] = []
    for exact in read_csv(args.exact_node_csv):
        if exact.get("variant") != args.variant:
            continue
        if exact.get("support_class") != args.support_class:
            continue
        if not matches_bias(exact, args.bias):
            continue
        node_id = exact.get("node_id", "")
        replay = replay_rows.get(node_id)
        proxy = proxy_rows.get(node_id)
        if replay is None or proxy is None:
            continue

        proxy_impact = finite_float(proxy.get("impact_source_s_inv"))
        proxy_e_required = required_from_proxy(proxy, "electron")
        proxy_h_required = required_from_proxy(proxy, "hole")
        exact_all = optional_float(replay.get("all_combined_source"))
        exact_edgeavg = optional_float(replay.get("sentaurus_edgeavg_current_scaled_source"))
        exact_e_required = optional_float(replay.get("electron_required_source"))
        if exact_e_required is None:
            exact_e_required = nonimpact_terms(exact, "electron")
        exact_h_required = optional_float(replay.get("hole_required_source"))
        if exact_h_required is None:
            exact_h_required = nonimpact_terms(exact, "hole")

        edgeavg_scale = safe_ratio(exact_edgeavg, proxy_impact)
        all_scale = safe_ratio(exact_all, proxy_impact)
        e_required_scale = safe_ratio(exact_e_required, proxy_e_required)
        h_required_scale = safe_ratio(exact_h_required, proxy_h_required)
        continuity_scale = optional_float(exact.get("inferred_continuity_scale_from_srh"))
        inverse_continuity_scale = safe_ratio(1.0, continuity_scale)
        node_volume = optional_float(exact.get("node_volume_m2"))

        rows.append({
            "node_id": node_id,
            "x_m": exact.get("x", exact.get("x_um", "")),
            "y_m": exact.get("y", exact.get("y_um", "")),
            "support_class": exact.get("support_class", ""),
            "node_volume_m2": node_volume,
            "inferred_continuity_scale_from_srh": continuity_scale,
            "inverse_inferred_continuity_scale_from_srh": inverse_continuity_scale,
            "exact_all_combined_source": exact_all,
            "exact_sentaurus_edgeavg_current_scaled_source": exact_edgeavg,
            "exact_electron_required_source": exact_e_required,
            "exact_hole_required_source": exact_h_required,
            "proxy_variant": args.proxy_variant,
            "transport_model": args.transport_model,
            "proxy_impact_source_s_inv": proxy_impact,
            "proxy_electron_required_s_inv": proxy_e_required,
            "proxy_hole_required_s_inv": proxy_h_required,
            "proxy_electron_residual_over_impact": safe_ratio(proxy_e_required - proxy_impact, proxy_impact),
            "proxy_hole_residual_over_impact": safe_ratio(proxy_h_required - proxy_impact, proxy_impact),
            "exact_edgeavg_current_over_proxy_impact_scale": edgeavg_scale,
            "exact_all_combined_over_proxy_impact_scale": all_scale,
            "exact_electron_required_over_proxy_required_scale": e_required_scale,
            "exact_hole_required_over_proxy_required_scale": h_required_scale,
            "electron_required_scale_over_edgeavg_scale": safe_ratio(e_required_scale, edgeavg_scale),
            "hole_required_scale_over_edgeavg_scale": safe_ratio(h_required_scale, edgeavg_scale),
            "edgeavg_scale_over_inverse_srh_scale": safe_ratio(edgeavg_scale, inverse_continuity_scale),
            "edgeavg_scale_over_node_volume": safe_ratio(edgeavg_scale, node_volume),
            "edgeavg_scale_over_q_node_volume": safe_ratio(edgeavg_scale, Q * node_volume if node_volume is not None else None),
        })
    rows.sort(key=lambda row: finite_float(row.get("y_m")))
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def metric_summary(rows: list[dict[str, Any]], key: str) -> dict[str, float | None]:
    values = finite_values(rows, key)
    if not values:
        return {"median": None, "mean": None, "min": None, "max": None}
    return {
        "median": statistics.median(values),
        "mean": statistics.fmean(values),
        "min": min(values),
        "max": max(values),
    }


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    with open_path(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> int:
    args = parse_args()
    rows = build_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "continuity_scaling_decomposition_nodes.csv", rows)
    payload = clean_json({
        "variant": args.variant,
        "proxy_variant": args.proxy_variant,
        "transport_model": args.transport_model,
        "bias_V": args.bias,
        "support_class": args.support_class,
        "row_count": len(rows),
        "metrics": {
            key: metric_summary(rows, key)
            for key in SUMMARY_METRICS
        },
    })
    write_text(
        args.out_dir / "continuity_scaling_decomposition_summary.json",
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
