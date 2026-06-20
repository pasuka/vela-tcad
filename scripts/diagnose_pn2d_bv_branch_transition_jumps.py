#!/usr/bin/env python3
"""Summarize same-bias BV errors and adjacent-state density jumps."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

import compare_pn2d_bv_multibias_fields as field_compare


DENSITY_QUANTITIES = {"electron_density", "hole_density"}


@dataclass(frozen=True)
class Point:
    bias: float
    label: str
    signed_dir: Path
    sentaurus_dir: Path
    vtk_path: Path


def main() -> int:
    args = parse_args()
    quantities = field_compare.parse_csv_list(args.quantities)
    regions = field_compare.parse_csv_list(args.regions)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    points = load_points(args)
    candidate_curve = (
        field_compare.load_curve_points(args.curve_candidate)
        if args.curve_candidate is not None
        else []
    )

    summary_rows = same_bias_rows(points, quantities, regions, candidate_curve, args.current_floor)
    jump_rows = adjacent_jump_rows(points, quantities, regions)

    write_csv(args.out_dir / "branch_transition_summary.csv", summary_rows)
    write_csv(args.out_dir / "branch_transition_adjacent_jumps.csv", jump_rows)
    report = {
        "point_count": len(points),
        "points": [
            {
                "bias_V": point.bias,
                "label": point.label,
                "signed_dir": str(point.signed_dir),
                "sentaurus_dir": str(point.sentaurus_dir),
                "vela_vtk": str(point.vtk_path),
            }
            for point in points
        ],
        "summary_row_count": len(summary_rows),
        "adjacent_jump_row_count": len(jump_rows),
        "quantities": quantities,
        "regions": regions,
    }
    (args.out_dir / "branch_transition_summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signed-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--points",
        required=True,
        help=(
            "Comma-separated entries 'bias:signed_dir_label' or "
            "'bias:signed_dir_label:vtk_path'."
        ),
    )
    parser.add_argument("--curve-candidate", type=Path)
    parser.add_argument("--current-floor", type=float, default=1.0e-25)
    parser.add_argument("--quantities", default="electron_density,hole_density")
    parser.add_argument("--regions", default="all,junction")
    return parser.parse_args()


def load_points(args: argparse.Namespace) -> list[Point]:
    points: list[Point] = []
    for raw in field_compare.parse_csv_list(args.points):
        parts = raw.split(":")
        if len(parts) not in {2, 3}:
            raise SystemExit(f"invalid point entry: {raw}")
        bias = float(parts[0])
        label = parts[1]
        signed_dir = args.signed_root / label
        metadata_path = signed_dir / "signed_field_partition_summary.json"
        if not metadata_path.exists():
            raise SystemExit(f"missing signed summary metadata: {metadata_path}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        sentaurus_dir = resolve_path(Path(str(metadata.get("sentaurus_dir", ""))))
        if len(parts) == 3:
            vtk_path = resolve_path(Path(parts[2]), base=args.signed_root)
        else:
            vtk_path = resolve_path(Path(str(metadata.get("vela_vtk", ""))))
        if not sentaurus_dir.exists():
            raise SystemExit(f"missing Sentaurus export: {sentaurus_dir}")
        if not vtk_path.exists():
            raise SystemExit(f"missing Vela VTK: {vtk_path}")
        points.append(Point(bias, label, signed_dir, sentaurus_dir, vtk_path))
    return points


def resolve_path(path: Path, base: Path | None = None) -> Path:
    if path.is_absolute():
        return path
    if base is not None and (base / path).exists():
        return base / path
    return Path.cwd() / path


def same_bias_rows(
    points: list[Point],
    quantities: list[str],
    regions: list[str],
    candidate_curve: list[tuple[float, float]],
    current_floor: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for point in points:
        summary_by_key = load_signed_summary(point.signed_dir)
        sentaurus_current = load_sentaurus_current(point.sentaurus_dir)
        vela_current = field_compare.interpolate_curve(candidate_curve, point.bias)
        current_log_ratio = current_log10_ratio(sentaurus_current, vela_current, current_floor)
        current_status = current_status_for(sentaurus_current, vela_current, current_floor)

        for quantity in quantities:
            for region in regions:
                source = summary_by_key.get((quantity, region), {})
                row = {
                    "bias_V": point.bias,
                    "label": point.label,
                    "quantity": quantity,
                    "region": region,
                    "status": source.get("status", "missing_summary"),
                    "sentaurus_current_A": sentaurus_current,
                    "vela_current_A_per_um": vela_current,
                    "current_log10_abs_ratio": current_log_ratio,
                    "current_abs_log10_error": abs(current_log_ratio)
                    if math.isfinite(current_log_ratio)
                    else math.nan,
                    "current_status": current_status,
                }
                for key in [
                    "signed_log10_median_vela_over_sentaurus",
                    "signed_log10_p95_abs",
                    "signed_delta_median",
                    "signed_delta_p95_abs",
                    "sentaurus_median",
                    "vela_median",
                    "vela_over_sentaurus_median",
                    "direction",
                    "count",
                ]:
                    row[key] = source.get(key, "")
                rows.append(row)
    return rows


def load_signed_summary(signed_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = signed_dir / "signed_field_partition_summary.csv"
    if not path.exists():
        raise SystemExit(f"missing signed partition CSV: {path}")
    with path.open(newline="") as handle:
        return {
            (row.get("quantity", ""), row.get("region", "")): row
            for row in csv.DictReader(handle)
        }


def load_sentaurus_current(sentaurus_dir: Path) -> float | None:
    fields_dir = sentaurus_dir / "fields"
    candidates = [
        fields_dir / "ContactCurrentFlux_region2.csv",
        fields_dir / "ContactCurrent_region2.csv",
        sentaurus_dir / "ContactCurrentFlux_region2.csv",
    ]
    for path in candidates:
        if not path.exists():
            continue
        total = 0.0
        found = False
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            component = next(
                (name for name in (reader.fieldnames or []) if name != "node_id"),
                None,
            )
            if component is None:
                continue
            for row in reader:
                value = row.get(component, "")
                if value == "":
                    continue
                total += float(value)
                found = True
        if found:
            return total
    return None


def current_status_for(reference: float | None, candidate: float | None, floor: float) -> str:
    if reference is None:
        return "missing_reference"
    if candidate is None:
        return "missing_candidate"
    if abs(reference) < floor or abs(candidate) < floor:
        return "below_current_floor"
    if reference == 0.0 or candidate == 0.0:
        return "zero_current"
    return "ok"


def current_log10_ratio(reference: float | None, candidate: float | None, floor: float) -> float:
    if current_status_for(reference, candidate, floor) != "ok":
        return math.nan
    assert reference is not None
    assert candidate is not None
    return math.log10(abs(candidate) / abs(reference))


def adjacent_jump_rows(
    points: list[Point],
    quantities: list[str],
    regions: list[str],
) -> list[dict[str, Any]]:
    if len(points) < 2:
        return []
    sx, sy, _ = field_compare.load_sentaurus_mesh(points[0].sentaurus_dir)
    node_ids = field_compare.load_sentaurus_node_ids(points[0].sentaurus_dir)
    masks = {
        "all": np.ones(len(node_ids), dtype=bool),
        "centerline": centerline_mask(sy),
        **field_compare.load_region_masks(points[0].sentaurus_dir, sx, sy),
    }

    state_cache = {point.label: load_vela_state(point, sx, sy) for point in points}
    rows: list[dict[str, Any]] = []
    for before, after in zip(points, points[1:]):
        before_state = state_cache[before.label]
        after_state = state_cache[after.label]
        for quantity in quantities:
            if quantity not in DENSITY_QUANTITIES:
                continue
            before_values = before_state.get(quantity)
            after_values = after_state.get(quantity)
            if before_values is None or after_values is None:
                continue
            signed_jump = signed_log10_jump(before_values, after_values)
            abs_jump = np.abs(signed_jump)
            for region in regions:
                mask = masks.get(region)
                if mask is None:
                    raise SystemExit(f"unknown region: {region}")
                rows.append(
                    summarize_jump(
                        before,
                        after,
                        quantity,
                        region,
                        mask,
                        node_ids,
                        sx,
                        sy,
                        before_values,
                        after_values,
                        signed_jump,
                        abs_jump,
                        before_state,
                        after_state,
                    )
                )
    return rows


def centerline_mask(y_um: np.ndarray) -> np.ndarray:
    indices = field_compare.centerline_indices(y_um)
    mask = np.zeros(len(y_um), dtype=bool)
    mask[indices] = True
    return mask


def load_vela_state(point: Point, sx: np.ndarray, sy: np.ndarray) -> dict[str, np.ndarray]:
    vtk = field_compare.parse_vtk(point.vtk_path)
    vx = np.asarray(vtk["x"], dtype=float)
    vy = np.asarray(vtk["y"], dtype=float)
    result: dict[str, np.ndarray] = {}
    for quantity, spec in field_compare.FIELD_SPECS.items():
        raw = vtk["scalars"].get(spec["vela"])
        if raw is None:
            continue
        result[quantity] = field_compare.nearest_values(
            vx,
            vy,
            np.asarray(raw, dtype=float) * float(spec["scale"]),
            sx,
            sy,
        )
    return result


def signed_log10_jump(before: np.ndarray, after: np.ndarray) -> np.ndarray:
    eps = 1.0e-300
    mask = np.isfinite(before) & np.isfinite(after)
    result = np.full_like(before, math.nan, dtype=float)
    result[mask] = (
        np.log10(np.maximum(np.abs(after[mask]), eps))
        - np.log10(np.maximum(np.abs(before[mask]), eps))
    )
    return result


def summarize_jump(
    before: Point,
    after: Point,
    quantity: str,
    region: str,
    mask: np.ndarray,
    node_ids: list[int],
    x_um: np.ndarray,
    y_um: np.ndarray,
    before_values: np.ndarray,
    after_values: np.ndarray,
    signed_jump: np.ndarray,
    abs_jump: np.ndarray,
    before_state: dict[str, np.ndarray],
    after_state: dict[str, np.ndarray],
) -> dict[str, Any]:
    finite = mask & np.isfinite(signed_jump)
    if not np.any(finite):
        return {
            "from_bias_V": before.bias,
            "to_bias_V": after.bias,
            "quantity": quantity,
            "region": region,
            "status": "empty",
            "count": 0,
        }
    indices = np.where(finite)[0]
    top_index = int(indices[np.argmax(abs_jump[finite])])
    row = {
        "from_bias_V": before.bias,
        "to_bias_V": after.bias,
        "from_label": before.label,
        "to_label": after.label,
        "quantity": quantity,
        "region": region,
        "status": "ok",
        "count": int(np.sum(finite)),
        "signed_log10_jump_median": float(np.nanmedian(signed_jump[finite])),
        "signed_log10_jump_p95_abs": float(np.nanpercentile(abs_jump[finite], 95)),
        "top_node_id": node_ids[top_index],
        "top_x_um": float(x_um[top_index]),
        "top_y_um": float(y_um[top_index]),
        "top_value_before": float(before_values[top_index]),
        "top_value_after": float(after_values[top_index]),
        "top_signed_log10_jump": float(signed_jump[top_index]),
        "top_abs_log10_jump": float(abs_jump[top_index]),
        "top_ratio_after_over_before": 10.0 ** float(signed_jump[top_index]),
    }
    for field in ["potential", "electric_field"]:
        before_field = before_state.get(field)
        after_field = after_state.get(field)
        if before_field is not None and after_field is not None:
            row[f"top_{field}_before"] = float(before_field[top_index])
            row[f"top_{field}_after"] = float(after_field[top_index])
            row[f"top_{field}_delta"] = float(after_field[top_index] - before_field[top_index])
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
