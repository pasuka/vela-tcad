#!/usr/bin/env python3
"""Join PN2D BV p99 avalanche support nodes with local Sentaurus/Vela factors."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Q = 1.602176634e-19


CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "sentaurus_active",
    "vela_active",
    "sentaurus_source_cm3_s",
    "vela_source_cm3_s",
    "source_ratio_vela_over_sentaurus",
    "sentaurus_electric_field_V_cm",
    "vela_electric_field_V_cm",
    "electric_field_ratio_vela_over_sentaurus",
    "vela_electron_high_field_drive_V_cm",
    "vela_hole_high_field_drive_V_cm",
    "sentaurus_electron_mobility_cm2_V_s",
    "vela_electron_mobility_cm2_V_s",
    "electron_mobility_ratio_vela_over_sentaurus",
    "sentaurus_hole_mobility_cm2_V_s",
    "vela_hole_mobility_cm2_V_s",
    "hole_mobility_ratio_vela_over_sentaurus",
    "sentaurus_electron_density_cm3",
    "vela_electron_density_cm3",
    "electron_density_ratio_vela_over_sentaurus",
    "sentaurus_hole_density_cm3",
    "vela_hole_density_cm3",
    "hole_density_ratio_vela_over_sentaurus",
    "sentaurus_electron_current_density_A_cm2",
    "sentaurus_hole_current_density_A_cm2",
    "sentaurus_current_density_sum_abs_A_cm2",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--match-tolerance-um", type=float, default=1.0e-4)
    parser.add_argument("--vela-source-scale", type=float, default=1.0e-6)
    parser.add_argument("--vela-field-scale-to-v-cm", type=float, default=1.0)
    parser.add_argument("--vela-mobility-scale-to-cm2-v-s", type=float, default=1.0e4)
    parser.add_argument("--vela-density-scale-to-cm3", type=float, default=1.0e-6)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def coord_key(x_um: float, y_um: float, tolerance_um: float) -> tuple[int, int]:
    return (round(x_um / tolerance_um), round(y_um / tolerance_um))


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


def load_sentaurus_field(export_dir: Path, name: str) -> dict[int, float]:
    path = export_dir / "fields" / f"{name}_region0.csv"
    if not path.exists():
        return {}
    values: dict[int, float] = {}
    for row in read_csv_rows(path):
        value = optional_float(row.get("component0"))
        if value is not None:
            values[int(row["node_id"])] = value
    return values


def vtk_lookup(
    points: list[tuple[float, float]],
    scalars: dict[str, list[float]],
    tolerance_um: float,
) -> dict[tuple[int, int], dict[str, float]]:
    result: dict[tuple[int, int], dict[str, float]] = {}
    for index, (x_um, y_um) in enumerate(points):
        item: dict[str, float] = {}
        for name, values in scalars.items():
            if index < len(values):
                item[name] = values[index]
        result[coord_key(x_um, y_um, tolerance_um)] = item
    return result


def field_value(fields: dict[str, dict[int, float]], name: str, node_id: int) -> float | None:
    return fields.get(name, {}).get(node_id)


def vtk_value(vtk_item: dict[str, float] | None, name: str, scale: float = 1.0) -> float | None:
    if vtk_item is None or name not in vtk_item:
        return None
    return vtk_item[name] * scale


def make_factor_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    support_rows = [
        row for row in read_csv_rows(args.support_csv)
        if row.get("support_class") and row.get("support_class") != "inactive"
    ]
    fields = {
        name: load_sentaurus_field(args.sentaurus_dir, name)
        for name in [
            "ImpactIonization",
            "ElectricField",
            "eMobility",
            "hMobility",
            "eDensity",
            "hDensity",
            "eCurrentDensity",
            "hCurrentDensity",
        ]
    }
    points, scalars = parse_vtk(args.vela_vtk)
    vtk_by_coord = vtk_lookup(points, scalars, args.match_tolerance_um)
    rows: list[dict[str, Any]] = []
    for support in support_rows:
        node_id = int(support["node_id"])
        x_um = float(support["x_um"])
        y_um = float(support["y_um"])
        vtk_item = vtk_by_coord.get(coord_key(x_um, y_um, args.match_tolerance_um))
        sentaurus_source = optional_float(support.get("sentaurus_avalanche_cm3_s"))
        if sentaurus_source is None:
            sentaurus_source = field_value(fields, "ImpactIonization", node_id)
        vela_source = optional_float(support.get("vela_avalanche_cm3_s"))
        if vela_source is None:
            vela_source = vtk_value(vtk_item, "AvalancheGeneration", args.vela_source_scale)
        s_efield = field_value(fields, "ElectricField", node_id)
        v_efield = vtk_value(vtk_item, "ElectricField", args.vela_field_scale_to_v_cm)
        s_emob = field_value(fields, "eMobility", node_id)
        v_emob = vtk_value(vtk_item, "ElectronMobility", args.vela_mobility_scale_to_cm2_v_s)
        s_hmob = field_value(fields, "hMobility", node_id)
        v_hmob = vtk_value(vtk_item, "HoleMobility", args.vela_mobility_scale_to_cm2_v_s)
        s_edens = field_value(fields, "eDensity", node_id)
        v_edens = vtk_value(vtk_item, "Electrons", args.vela_density_scale_to_cm3)
        s_hdens = field_value(fields, "hDensity", node_id)
        v_hdens = vtk_value(vtk_item, "Holes", args.vela_density_scale_to_cm3)
        s_jn = field_value(fields, "eCurrentDensity", node_id)
        s_jp = field_value(fields, "hCurrentDensity", node_id)
        rows.append({
            "node_id": node_id,
            "x_um": x_um,
            "y_um": y_um,
            "support_class": support["support_class"],
            "sentaurus_active": support.get("sentaurus_active", ""),
            "vela_active": support.get("vela_active", ""),
            "sentaurus_source_cm3_s": sentaurus_source,
            "vela_source_cm3_s": vela_source,
            "source_ratio_vela_over_sentaurus": ratio(vela_source, sentaurus_source),
            "sentaurus_electric_field_V_cm": s_efield,
            "vela_electric_field_V_cm": v_efield,
            "electric_field_ratio_vela_over_sentaurus": ratio(v_efield, s_efield),
            "vela_electron_high_field_drive_V_cm": vtk_value(
                vtk_item, "ElectronHighFieldDrive", args.vela_field_scale_to_v_cm),
            "vela_hole_high_field_drive_V_cm": vtk_value(
                vtk_item, "HoleHighFieldDrive", args.vela_field_scale_to_v_cm),
            "sentaurus_electron_mobility_cm2_V_s": s_emob,
            "vela_electron_mobility_cm2_V_s": v_emob,
            "electron_mobility_ratio_vela_over_sentaurus": ratio(v_emob, s_emob),
            "sentaurus_hole_mobility_cm2_V_s": s_hmob,
            "vela_hole_mobility_cm2_V_s": v_hmob,
            "hole_mobility_ratio_vela_over_sentaurus": ratio(v_hmob, s_hmob),
            "sentaurus_electron_density_cm3": s_edens,
            "vela_electron_density_cm3": v_edens,
            "electron_density_ratio_vela_over_sentaurus": ratio(v_edens, s_edens),
            "sentaurus_hole_density_cm3": s_hdens,
            "vela_hole_density_cm3": v_hdens,
            "hole_density_ratio_vela_over_sentaurus": ratio(v_hdens, s_hdens),
            "sentaurus_electron_current_density_A_cm2": s_jn,
            "sentaurus_hole_current_density_A_cm2": s_jp,
            "sentaurus_current_density_sum_abs_A_cm2": (
                abs(s_jn or 0.0) + abs(s_jp or 0.0)
            ),
        })
    return rows


def median(values: list[float]) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def mean(values: list[float]) -> float | None:
    clean = [value for value in values if math.isfinite(value)]
    if not clean:
        return None
    return sum(clean) / len(clean)


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for cls in sorted({str(row["support_class"]) for row in rows}):
        subset = [row for row in rows if row["support_class"] == cls]
        source_ratios = [
            float(row["source_ratio_vela_over_sentaurus"])
            for row in subset
            if row.get("source_ratio_vela_over_sentaurus") is not None
        ]
        field_ratios = [
            float(row["electric_field_ratio_vela_over_sentaurus"])
            for row in subset
            if row.get("electric_field_ratio_vela_over_sentaurus") is not None
        ]
        electron_mobility_ratios = [
            float(row["electron_mobility_ratio_vela_over_sentaurus"])
            for row in subset
            if row.get("electron_mobility_ratio_vela_over_sentaurus") is not None
        ]
        hole_mobility_ratios = [
            float(row["hole_mobility_ratio_vela_over_sentaurus"])
            for row in subset
            if row.get("hole_mobility_ratio_vela_over_sentaurus") is not None
        ]
        electron_density_ratios = [
            float(row["electron_density_ratio_vela_over_sentaurus"])
            for row in subset
            if row.get("electron_density_ratio_vela_over_sentaurus") is not None
        ]
        hole_density_ratios = [
            float(row["hole_density_ratio_vela_over_sentaurus"])
            for row in subset
            if row.get("hole_density_ratio_vela_over_sentaurus") is not None
        ]
        sentaurus_current_density = [
            float(row["sentaurus_current_density_sum_abs_A_cm2"])
            for row in subset
            if row.get("sentaurus_current_density_sum_abs_A_cm2") is not None
        ]
        electron_drive_ratios = [
            float(row["vela_electron_high_field_drive_V_cm"]) / float(row["sentaurus_electric_field_V_cm"])
            for row in subset
            if row.get("vela_electron_high_field_drive_V_cm") is not None
            and row.get("sentaurus_electric_field_V_cm") not in (None, 0.0)
        ]
        hole_drive_ratios = [
            float(row["vela_hole_high_field_drive_V_cm"]) / float(row["sentaurus_electric_field_V_cm"])
            for row in subset
            if row.get("vela_hole_high_field_drive_V_cm") is not None
            and row.get("sentaurus_electric_field_V_cm") not in (None, 0.0)
        ]
        result[cls] = {
            "count": len(subset),
            "source_ratio_median": median(source_ratios),
            "source_ratio_mean": mean(source_ratios),
            "electric_field_ratio_median": median(field_ratios),
            "electron_mobility_ratio_median": median(electron_mobility_ratios),
            "hole_mobility_ratio_median": median(hole_mobility_ratios),
            "electron_density_ratio_median": median(electron_density_ratios),
            "hole_density_ratio_median": median(hole_density_ratios),
            "electron_drive_over_sentaurus_field_median": median(electron_drive_ratios),
            "hole_drive_over_sentaurus_field_median": median(hole_drive_ratios),
            "sentaurus_current_density_sum_abs_median_A_cm2": median(sentaurus_current_density),
            "sentaurus_source_sum_cm3_s": sum(
                abs(float(row["sentaurus_source_cm3_s"] or 0.0)) for row in subset),
            "vela_source_sum_cm3_s": sum(
                abs(float(row["vela_source_cm3_s"] or 0.0)) for row in subset),
        }
    return result


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = make_factor_rows(args)
    summary = {
        "row_count": len(rows),
        "support_classes": class_summary(rows),
        "match_tolerance_um": args.match_tolerance_um,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "support_local_factors.csv", rows)
    (args.out_dir / "support_local_factors_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
