#!/usr/bin/env python3
"""Summarize active-edge impact feedback ownership hypotheses."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


CSV_FIELDS = [
    "support_class",
    "active_edge_generation_over_sentaurus_median",
    "active_edge_generation_over_sentaurus_mean",
    "active_endpoint_area_fraction_summed",
    "active_endpoint_area_fraction_median",
    "endpoint_feedback_over_sentaurus",
    "full_active_edge_feedback_over_sentaurus",
    "electron_required_impact_scale_median",
    "hole_required_impact_scale_median",
    "endpoint_feedback_times_electron_required_scale",
    "endpoint_feedback_times_hole_required_scale",
    "full_active_edge_feedback_times_electron_required_scale",
    "full_active_edge_feedback_times_hole_required_scale",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--active-edge-replay-summary", type=Path, required=True)
    parser.add_argument("--source-geometry-summary", type=Path, required=True)
    parser.add_argument("--impact-scale-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--variant", default="predictor")
    parser.add_argument("--state", default="predictor")
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def optional_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def product(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return a * b


def class_names(replay: dict[str, Any], geometry: dict[str, Any], scale: dict[str, Any], state: str) -> list[str]:
    names = set((replay.get("support_classes") or {}).keys())
    names.update((geometry.get("support_classes") or {}).keys())
    scale_classes = (((scale.get("impact_scale_sensitivity") or {})
        .get("states") or {}).get(state) or {}).get("support_classes") or {}
    names.update(scale_classes.keys())
    return sorted(str(name) for name in names)


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    replay = load_json(args.active_edge_replay_summary)
    geometry = load_json(args.source_geometry_summary)
    scale = load_json(args.impact_scale_summary)
    scale_classes = (((scale.get("impact_scale_sensitivity") or {})
        .get("states") or {}).get(args.state) or {}).get("support_classes") or {}
    rows: list[dict[str, Any]] = []
    for name in class_names(replay, geometry, scale, args.state):
        replay_variant = (((replay.get("support_classes") or {}).get(name) or {})
            .get("variants") or {}).get(args.variant) or {}
        geometry_class = (geometry.get("support_classes") or {}).get(name) or {}
        scale_class = scale_classes.get(name) or {}

        generation_median = optional_float(replay_variant.get("generation_over_sentaurus_median"))
        generation_mean = optional_float(replay_variant.get("generation_over_sentaurus_mean"))
        area_summed = optional_float(geometry_class.get("summed_active_endpoint_area_fraction"))
        area_median = optional_float(geometry_class.get("cxx_active_endpoint_area_fraction_median"))
        electron_scale = optional_float(scale_class.get("electron_required_impact_scale_median"))
        hole_scale = optional_float(scale_class.get("hole_required_impact_scale_median"))
        endpoint_feedback = product(generation_median, area_summed)
        full_feedback = generation_median
        rows.append({
            "support_class": name,
            "active_edge_generation_over_sentaurus_median": generation_median,
            "active_edge_generation_over_sentaurus_mean": generation_mean,
            "active_endpoint_area_fraction_summed": area_summed,
            "active_endpoint_area_fraction_median": area_median,
            "endpoint_feedback_over_sentaurus": endpoint_feedback,
            "full_active_edge_feedback_over_sentaurus": full_feedback,
            "electron_required_impact_scale_median": electron_scale,
            "hole_required_impact_scale_median": hole_scale,
            "endpoint_feedback_times_electron_required_scale": product(endpoint_feedback, electron_scale),
            "endpoint_feedback_times_hole_required_scale": product(endpoint_feedback, hole_scale),
            "full_active_edge_feedback_times_electron_required_scale": product(full_feedback, electron_scale),
            "full_active_edge_feedback_times_hole_required_scale": product(full_feedback, hole_scale),
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "impact_feedback_ownership_summary.csv", rows)
    support_summary = {}
    for row in rows:
        full_feedback = optional_float(row.get("full_active_edge_feedback_over_sentaurus"))
        endpoint_feedback = optional_float(row.get("endpoint_feedback_over_sentaurus"))
        support_summary[row["support_class"]] = {
            "endpoint_feedback_over_sentaurus": endpoint_feedback,
            "full_active_edge_feedback_over_sentaurus": full_feedback,
            "endpoint_feedback_closes_unit_feedback": (
                endpoint_feedback is not None and 0.8 <= endpoint_feedback <= 1.25
            ),
            "full_active_edge_closes_unit_feedback": (
                full_feedback is not None and 0.8 <= full_feedback <= 1.25
            ),
        }
    payload = clean_json({
        "variant": args.variant,
        "state": args.state,
        "row_count": len(rows),
        "support_classes": support_summary,
    })
    (args.out_dir / "impact_feedback_ownership_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())