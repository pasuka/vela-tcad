#!/usr/bin/env python3
"""Compare thresholded PN2D BV avalanche-generation support sets."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


ROW_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "sentaurus_avalanche_cm3_s",
    "vela_avalanche_cm3_s",
    "sentaurus_active",
    "vela_active",
    "support_class",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--percentile", type=float, default=99.0)
    parser.add_argument(
        "--absolute-threshold",
        type=float,
        default=None,
        help="Use this shared cm^-3 s^-1 threshold instead of per-side percentile thresholds.",
    )
    parser.add_argument(
        "--vela-scale",
        type=float,
        default=1.0e-6,
        help="Scale Vela VTK AvalancheGeneration into Sentaurus cm^-3 s^-1 units.",
    )
    parser.add_argument(
        "--match-tolerance-um",
        type=float,
        default=1.0e-4,
        help="Coordinate quantization tolerance for matching Sentaurus CSV nodes to ASCII VTK points.",
    )
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_sentaurus_nodes(path: Path) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    for row in read_csv_rows(path):
        nodes[int(row["id"])] = (float(row["x_um"]), float(row["y_um"]))
    return nodes


def load_sentaurus_avalanche(export_dir: Path) -> dict[int, float]:
    for name in ["ImpactIonization", "AvalancheGeneration"]:
        path = export_dir / "fields" / f"{name}_region0.csv"
        if path.exists():
            return {
                int(row["node_id"]): float(row["component0"])
                for row in read_csv_rows(path)
                if row.get("component0") not in (None, "")
            }
    raise FileNotFoundError(f"missing Sentaurus avalanche field under {export_dir / 'fields'}")


def parse_vtk(path: Path) -> tuple[list[tuple[float, float]], dict[str, list[float]]]:
    lines = path.read_text().splitlines()
    points: list[tuple[float, float]] = []
    scalars: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "POINTS":
            count = int(parts[1])
            i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < count * 3:
                values.extend(float(item) for item in lines[i].split())
                i += 1
            for index in range(0, count * 3, 3):
                points.append((values[index] * 1.0e6, values[index + 1] * 1.0e6))
            continue
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values = []
            while i < len(lines) and len(values) < len(points):
                next_parts = lines[i].split()
                if not next_parts:
                    i += 1
                    continue
                if next_parts[0] in {"SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA"}:
                    break
                values.extend(float(item) for item in next_parts)
                i += 1
            scalars[name] = values[:len(points)]
            continue
        i += 1
    return points, scalars


def coord_key(point: tuple[float, float], tolerance_um: float) -> tuple[int, int]:
    if tolerance_um <= 0.0:
        return (hash(round(point[0], 12)), hash(round(point[1], 12)))
    return (
        round(point[0] / tolerance_um),
        round(point[1] / tolerance_um),
    )


def percentile(values: list[float], pct: float) -> float:
    clean = sorted(abs(value) for value in values if math.isfinite(value))
    if not clean:
        return 0.0
    if pct <= 0.0:
        return clean[0]
    if pct >= 100.0:
        return clean[-1]
    idx = (len(clean) - 1) * pct / 100.0
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return clean[lo]
    return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)


def support_class(sentaurus_active: bool, vela_active: bool) -> str:
    if sentaurus_active and vela_active:
        return "overlap"
    if sentaurus_active:
        return "false_negative"
    if vela_active:
        return "false_positive"
    return "inactive"


def peak(rows: list[dict[str, Any]], key: str) -> dict[str, Any] | None:
    if not rows:
        return None
    return max(rows, key=lambda row: abs(float(row[key])))


def distance_um(lhs: dict[str, Any] | None, rhs: dict[str, Any] | None) -> float | None:
    if lhs is None or rhs is None:
        return None
    return math.hypot(float(lhs["x_um"]) - float(rhs["x_um"]), float(lhs["y_um"]) - float(rhs["y_um"]))


def make_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    nodes = load_sentaurus_nodes(args.sentaurus_dir / "nodes.csv")
    sentaurus_values = load_sentaurus_avalanche(args.sentaurus_dir)
    vtk_points, vtk_scalars = parse_vtk(args.vela_vtk)
    if "AvalancheGeneration" not in vtk_scalars:
        raise RuntimeError(f"VTK is missing AvalancheGeneration scalar: {args.vela_vtk}")
    vela_by_coord = {
        coord_key(point, args.match_tolerance_um): value * args.vela_scale
        for point, value in zip(vtk_points, vtk_scalars["AvalancheGeneration"])
    }
    matched: list[dict[str, Any]] = []
    missing = 0
    for node_id, point in nodes.items():
        value = vela_by_coord.get(coord_key(point, args.match_tolerance_um))
        if value is None:
            missing += 1
            continue
        if node_id not in sentaurus_values:
            continue
        matched.append({
            "node_id": node_id,
            "x_um": point[0],
            "y_um": point[1],
            "sentaurus_avalanche_cm3_s": sentaurus_values[node_id],
            "vela_avalanche_cm3_s": value,
        })
    s_threshold = (
        args.absolute_threshold
        if args.absolute_threshold is not None
        else percentile([row["sentaurus_avalanche_cm3_s"] for row in matched], args.percentile)
    )
    v_threshold = (
        args.absolute_threshold
        if args.absolute_threshold is not None
        else percentile([row["vela_avalanche_cm3_s"] for row in matched], args.percentile)
    )
    for row in matched:
        s_active = abs(float(row["sentaurus_avalanche_cm3_s"])) >= float(s_threshold)
        v_active = abs(float(row["vela_avalanche_cm3_s"])) >= float(v_threshold)
        row["sentaurus_active"] = int(s_active)
        row["vela_active"] = int(v_active)
        row["support_class"] = support_class(s_active, v_active)
    summary = summarize(matched, missing, args.percentile, s_threshold, v_threshold)
    summary["match_tolerance_um"] = args.match_tolerance_um
    return matched, summary


def sum_abs(rows: list[dict[str, Any]], key: str, predicate: str | None = None) -> float:
    selected = rows if predicate is None else [row for row in rows if row["support_class"] == predicate]
    return sum(abs(float(row[key])) for row in selected)


def summarize(
    rows: list[dict[str, Any]],
    missing_vela_points: int,
    pct: float,
    sentaurus_threshold: float,
    vela_threshold: float,
) -> dict[str, Any]:
    overlap = [row for row in rows if row["support_class"] == "overlap"]
    false_positive = [row for row in rows if row["support_class"] == "false_positive"]
    false_negative = [row for row in rows if row["support_class"] == "false_negative"]
    sentaurus_active = [row for row in rows if int(row["sentaurus_active"])]
    vela_active = [row for row in rows if int(row["vela_active"])]
    union_count = len(overlap) + len(false_positive) + len(false_negative)
    s_peak = peak(rows, "sentaurus_avalanche_cm3_s")
    v_peak = peak(rows, "vela_avalanche_cm3_s")
    return {
        "matched_points": len(rows),
        "missing_vela_points": missing_vela_points,
        "percentile": pct,
        "sentaurus_threshold_cm3_s": sentaurus_threshold,
        "vela_threshold_cm3_s": vela_threshold,
        "sentaurus_active_count": len(sentaurus_active),
        "vela_active_count": len(vela_active),
        "overlap_count": len(overlap),
        "false_positive_count": len(false_positive),
        "false_negative_count": len(false_negative),
        "union_count": union_count,
        "jaccard": (len(overlap) / union_count) if union_count else 1.0,
        "sentaurus_active_sum_cm3_s": sum_abs(sentaurus_active, "sentaurus_avalanche_cm3_s"),
        "vela_active_sum_cm3_s": sum_abs(vela_active, "vela_avalanche_cm3_s"),
        "overlap_sentaurus_sum_cm3_s": sum_abs(overlap, "sentaurus_avalanche_cm3_s"),
        "overlap_vela_sum_cm3_s": sum_abs(overlap, "vela_avalanche_cm3_s"),
        "false_positive_vela_sum_cm3_s": sum_abs(false_positive, "vela_avalanche_cm3_s"),
        "false_negative_sentaurus_sum_cm3_s": sum_abs(false_negative, "sentaurus_avalanche_cm3_s"),
        "sentaurus_peak": s_peak,
        "vela_peak": v_peak,
        "peak_separation_um": distance_um(s_peak, v_peak),
    }


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in ROW_FIELDS})


def main() -> int:
    args = parse_args()
    rows, summary = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.out_dir / "thresholded_avalanche_support_nodes.csv", rows)
    (args.out_dir / "thresholded_avalanche_support_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
