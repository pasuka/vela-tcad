#!/usr/bin/env python3
"""Compare PN2D 0V electric-field magnitude from Vela potential vs Sentaurus."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


REPORT_NAME = "pn2d_0v_electric_field_comparison.json"
MARKDOWN_NAME = "pn2d_0v_electric_field_comparison.md"
NODE_CSV_NAME = "pn2d_0v_electric_field_nodes.csv"
UNIT_CANDIDATES = {
    "V_per_m": 1.0,
    "V_per_cm": 1.0e-2,
    "V_per_um": 1.0e-6,
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def field_files(root: Path, field: str) -> list[Path]:
    exact = root / f"{field}_region0.csv"
    if exact.is_file():
        return [exact]
    return sorted(root.glob(f"{field}_region*.csv"))


def read_sentaurus_scalar(fields_root: Path, field: str) -> dict[int, float]:
    values: dict[int, float] = {}
    paths = field_files(fields_root, field)
    if not paths:
        raise FileNotFoundError(f"missing Sentaurus field export: {field}")
    for path in paths:
        for row in read_csv_rows(path):
            values[int(float(row["node_id"]))] = float(row["component0"])
    return values


def read_nodes(path: Path) -> list[dict[str, float]]:
    if not path.is_file():
        return []
    nodes: list[dict[str, float]] = []
    for row in read_csv_rows(path):
        nodes.append({
            "id": int(float(row["id"])),
            "x_um": float(row.get("x_um", row.get("x", 0.0))),
            "y_um": float(row.get("y_um", row.get("y", 0.0))),
        })
    return sorted(nodes, key=lambda item: item["id"])


def read_doping(path: Path) -> dict[int, float]:
    if not path.is_file():
        return {}
    doping: dict[int, float] = {}
    for row in read_csv_rows(path):
        node_id = int(float(row["node_id"]))
        donors = float(row.get("donors_cm3", 0.0) or 0.0)
        acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0)
        doping[node_id] = donors - acceptors
    return doping


def read_contacts(reference_root: Path, mesh: dict[str, Any]) -> dict[str, list[int]]:
    contacts_csv = reference_root / "contacts.csv"
    if contacts_csv.is_file():
        contacts: dict[str, list[int]] = {}
        for row in read_csv_rows(contacts_csv):
            raw = row.get("node_ids", "")
            contacts[row["name"]] = [int(item) for item in raw.split(";") if item.strip()]
        return contacts
    contacts = {}
    for entry in mesh.get("contacts", []):
        contacts[str(entry.get("name", ""))] = [int(node) for node in entry.get("node_ids", [])]
    return contacts


def parse_vtk(path: Path) -> tuple[list[tuple[float, float, float]], dict[str, list[float]]]:
    lines = path.read_text().splitlines()
    points: list[tuple[float, float, float]] = []
    fields: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "POINTS":
            npoints = int(parts[1])
            i += 1
            raw: list[float] = []
            while i < len(lines) and len(raw) < 3 * npoints:
                raw.extend(float(item) for item in lines[i].split())
                i += 1
            points = [
                (raw[3 * idx], raw[3 * idx + 1], raw[3 * idx + 2])
                for idx in range(npoints)
            ]
            continue
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < len(points):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                head = line.split()[0]
                if head in {"SCALARS", "CELL_DATA", "POINT_DATA", "VECTORS"}:
                    break
                values.extend(float(item) for item in line.split())
                i += 1
            fields[name] = values[:len(points)]
            continue
        i += 1
    return points, fields


def mesh_triangles(mesh: dict[str, Any]) -> list[list[int]]:
    triangles: list[list[int]] = []
    for entry in mesh.get("triangles", []):
        raw = entry.get("node_ids", entry.get("nodes", entry))
        if isinstance(raw, dict):
            continue
        ids = [int(node) for node in raw]
        if len(ids) == 3:
            triangles.append(ids)
    return triangles


def triangle_gradient(points: list[tuple[float, float, float]],
                      potential: list[float],
                      ids: list[int]) -> float | None:
    i0, i1, i2 = ids
    if max(ids) >= len(points) or max(ids) >= len(potential):
        return None
    x0, y0, _ = points[i0]
    x1, y1, _ = points[i1]
    x2, y2, _ = points[i2]
    p0, p1, p2 = potential[i0], potential[i1], potential[i2]
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if abs(det) <= 1.0e-300:
        return None
    dpsi_dx = ((p1 - p0) * (y2 - y0) - (p2 - p0) * (y1 - y0)) / det
    dpsi_dy = ((x1 - x0) * (p2 - p0) - (x2 - x0) * (p1 - p0)) / det
    return math.hypot(dpsi_dx, dpsi_dy)


def derive_nodal_field_v_per_m(points: list[tuple[float, float, float]],
                               potential: list[float],
                               triangles: list[list[int]]) -> dict[int, float]:
    accum: dict[int, list[float]] = {}
    for ids in triangles:
        magnitude = triangle_gradient(points, potential, ids)
        if magnitude is None:
            continue
        for node_id in ids:
            accum.setdefault(node_id, []).append(magnitude)
    return {node_id: sum(values) / len(values) for node_id, values in accum.items() if values}


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = int(round((len(ordered) - 1) * p))
    return ordered[index]


def value_stats(values: list[float]) -> dict[str, float | None]:
    return {
        "points": len(values),
        "min": min(values) if values else None,
        "mean": sum(values) / len(values) if values else None,
        "median": statistics.median(values) if values else None,
        "p95": quantile(values, 0.95),
        "max": max(values) if values else None,
    }


def compare_values(reference: dict[int, float],
                   candidate: dict[int, float],
                   groups: dict[str, set[int]]) -> dict[str, Any]:
    common_all = sorted(set(reference) & set(candidate))
    out: dict[str, Any] = {}
    for name, ids in groups.items():
        common = [node for node in common_all if node in ids]
        diffs = [candidate[node] - reference[node] for node in common]
        abs_diffs = [abs(value) for value in diffs]
        rel_diffs = [
            abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
            for node in common
        ]
        worst = None
        if common:
            worst_node = max(common, key=lambda node: abs(candidate[node] - reference[node]))
            worst = {
                "node_id": worst_node,
                "reference": reference[worst_node],
                "candidate": candidate[worst_node],
                "abs_diff": abs(candidate[worst_node] - reference[worst_node]),
                "rel_diff": abs(candidate[worst_node] - reference[worst_node])
                / max(abs(reference[worst_node]), 1.0e-300),
            }
        out[name] = {
            "points": len(common),
            "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
            "median_abs_diff": statistics.median(abs_diffs) if abs_diffs else None,
            "p95_abs_diff": quantile(abs_diffs, 0.95),
            "max_abs_diff": max(abs_diffs) if abs_diffs else None,
            "mean_rel_diff": sum(rel_diffs) / len(rel_diffs) if rel_diffs else None,
            "median_rel_diff": statistics.median(rel_diffs) if rel_diffs else None,
            "p95_rel_diff": quantile(rel_diffs, 0.95),
            "max_rel_diff": max(rel_diffs) if rel_diffs else None,
            "worst": worst,
        }
    return out


def selected_node_groups(nodes: list[dict[str, float]],
                         contacts: dict[str, list[int]],
                         doping: dict[int, float],
                         compared_ids: set[int]) -> dict[str, set[int]]:
    groups: dict[str, set[int]] = {"all": set(compared_ids)}
    if nodes:
        ys = [node["y_um"] for node in nodes]
        y_mid = 0.5 * (min(ys) + max(ys))
        y_span = max(ys) - min(ys)
        y_halfwidth = max(y_span * 0.02, 0.01)
        groups["centerline"] = {
            node["id"] for node in nodes if abs(node["y_um"] - y_mid) <= y_halfwidth
        }
        xs = [node["x_um"] for node in nodes]
        x_mid = 0.5 * (min(xs) + max(xs))
        x_halfwidth = max((max(xs) - min(xs)) * 0.02, 0.02)
        junction = {node_id for node_id, value in doping.items() if abs(value) <= 1.0e-30}
        groups["junction_near"] = junction or {
            node["id"] for node in nodes if abs(node["x_um"] - x_mid) <= x_halfwidth
        }
    for name, ids in contacts.items():
        groups[f"contact:{name}"] = set(ids)
    return groups


def median_ratio(reference: dict[int, float], candidate: dict[int, float]) -> float | None:
    ratios: list[float] = []
    for node in sorted(set(reference) & set(candidate)):
        ref = reference[node]
        if abs(ref) <= 1.0e-300:
            continue
        ratios.append(candidate[node] / ref)
    return statistics.median(ratios) if ratios else None


def select_unit(reference: dict[int, float],
                base_v_per_m: dict[int, float],
                groups: dict[str, set[int]]) -> tuple[str, dict[str, Any]]:
    reports: dict[str, Any] = {}
    common = sorted(set(reference) & set(base_v_per_m))
    reference_values = [abs(reference[node]) for node in common]
    reference_median = statistics.median(reference_values) if reference_values else None
    reference_mean = sum(reference_values) / len(reference_values) if reference_values else None
    for name, multiplier in UNIT_CANDIDATES.items():
        candidate = {node: value * multiplier for node, value in base_v_per_m.items()}
        candidate_values = [abs(candidate[node]) for node in common]
        candidate_median = statistics.median(candidate_values) if candidate_values else None
        candidate_mean = sum(candidate_values) / len(candidate_values) if candidate_values else None
        entry = {
            "multiplier_from_V_per_m": multiplier,
            "median_candidate_over_reference": median_ratio(reference, candidate),
            "median_magnitude_ratio": (
                candidate_median / reference_median
                if candidate_median is not None and reference_median not in (None, 0.0)
                else None
            ),
            "mean_magnitude_ratio": (
                candidate_mean / reference_mean
                if candidate_mean is not None and reference_mean not in (None, 0.0)
                else None
            ),
            "stats": compare_values(reference, candidate, groups),
        }
        reports[name] = entry

    def score(name: str) -> tuple[float, float, float]:
        ratio = reports[name]["median_magnitude_ratio"]
        mean_ratio = reports[name]["mean_magnitude_ratio"]
        all_stats = reports[name]["stats"]["all"]
        ratio_score = abs(math.log10(abs(ratio))) if ratio not in (None, 0.0) else math.inf
        mean_ratio_score = abs(math.log10(abs(mean_ratio))) if mean_ratio not in (None, 0.0) else math.inf
        median_abs = all_stats["median_abs_diff"] if all_stats["median_abs_diff"] is not None else math.inf
        return (ratio_score, mean_ratio_score, median_abs)

    return min(reports, key=score), reports


def top_differences(reference: dict[int, float],
                    candidate: dict[int, float],
                    limit: int = 10) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for node in sorted(set(reference) & set(candidate)):
        abs_diff = abs(candidate[node] - reference[node])
        rows.append({
            "node_id": node,
            "sentaurus": reference[node],
            "vela": candidate[node],
            "abs_diff": abs_diff,
            "rel_diff": abs_diff / max(abs(reference[node]), 1.0e-300),
        })
    return sorted(rows, key=lambda item: float(item["abs_diff"]), reverse=True)[:limit]


def find_default_vtk(reference_root: Path) -> Path:
    candidates = sorted((reference_root / "reports").glob("**/*.vtk"))
    if not candidates:
        raise FileNotFoundError("no VTK file found; pass --vtk")
    preferred = [path for path in candidates if "0v" in path.name.lower()]
    return preferred[0] if preferred else candidates[0]


def write_node_csv(path: Path,
                   reference: dict[int, float],
                   candidate: dict[int, float],
                   unit: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["node_id", "unit", "sentaurus_E", "vela_E", "abs_diff", "rel_diff"],
        )
        writer.writeheader()
        for node in sorted(set(reference) & set(candidate)):
            abs_diff = abs(candidate[node] - reference[node])
            writer.writerow({
                "node_id": node,
                "unit": unit,
                "sentaurus_E": f"{reference[node]:.17g}",
                "vela_E": f"{candidate[node]:.17g}",
                "abs_diff": f"{abs_diff:.17g}",
                "rel_diff": f"{abs_diff / max(abs(reference[node]), 1.0e-300):.17g}",
            })


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Electric Field Comparison",
        "",
        f"Status: {report.get('status')}",
        f"Source: {report.get('source')}",
        f"Best unit: {report.get('unit_selection', {}).get('best_unit')}",
        "",
        "Vela electric field is derived from the linear triangle gradient of VTK `Potential`; "
        "Sentaurus uses `ElectricField_region0.csv`.",
        "",
        "| Group | Points | Mean abs | Median abs | P95 abs | Max abs |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    stats = report.get("stats", {}).get("best_unit", {})
    for name, item in stats.items():
        lines.append(
            f"| {name} | {item.get('points')} | {item.get('mean_abs_diff')} | "
            f"{item.get('median_abs_diff')} | {item.get('p95_abs_diff')} | {item.get('max_abs_diff')} |"
        )
    lines.extend(["", "Worst nodes:"])
    for item in report.get("top_differences", []):
        lines.append(
            f"- node {item.get('node_id')}: Sentaurus={item.get('sentaurus')}, "
            f"Vela={item.get('vela')}, abs_diff={item.get('abs_diff')}"
        )
    if report.get("error"):
        lines.extend(["", f"Error: {report['error']}"])
    return "\n".join(lines) + "\n"


def build_report(reference_root: Path, vtk_path: Path) -> tuple[dict[str, Any], dict[int, float], str]:
    fields_root = reference_root / "sim_fields" / "0v" / "fields"
    sentaurus = read_sentaurus_scalar(fields_root, "ElectricField")
    mesh = read_json(reference_root / "vela" / "mesh.json")
    points, vtk_fields = parse_vtk(vtk_path)
    if "Potential" not in vtk_fields:
        raise RuntimeError("missing Vela VTK scalar field: Potential")
    triangles = mesh_triangles(mesh)
    if not triangles:
        raise RuntimeError("mesh.json does not contain triangles")
    vela_v_per_m = derive_nodal_field_v_per_m(points, vtk_fields["Potential"], triangles)
    nodes = read_nodes(reference_root / "nodes.csv")
    doping = read_doping(reference_root / "doping.csv")
    contacts = read_contacts(reference_root, mesh)
    compared_ids = set(sentaurus) & set(vela_v_per_m)
    groups = selected_node_groups(nodes, contacts, doping, compared_ids)
    best_unit, unit_reports = select_unit(sentaurus, vela_v_per_m, groups)
    best_candidate = {
        node: value * UNIT_CANDIDATES[best_unit] for node, value in vela_v_per_m.items()
    }
    best_stats = unit_reports[best_unit]["stats"]
    report = {
        "status": "pass" if compared_ids else "fail",
        "source": "derived_from_potential_gradient",
        "reference_field": "ElectricField",
        "vtk_field": "Potential",
        "vtk_path": str(vtk_path),
        "points": {
            "sentaurus": len(sentaurus),
            "vela": len(vela_v_per_m),
            "compared": len(compared_ids),
            "missing_vela_nodes": len(set(sentaurus) - set(vela_v_per_m)),
            "missing_sentaurus_nodes": len(set(vela_v_per_m) - set(sentaurus)),
        },
        "value_stats": {
            "sentaurus": value_stats([sentaurus[node] for node in sorted(compared_ids)]),
            "vela_V_per_m": value_stats([vela_v_per_m[node] for node in sorted(compared_ids)]),
            "vela_best_unit": value_stats([best_candidate[node] for node in sorted(compared_ids)]),
        },
        "unit_selection": {
            "best_unit": best_unit,
            "candidates": unit_reports,
        },
        "stats": {
            "best_unit": best_stats,
        },
        "top_differences": top_differences(sentaurus, best_candidate),
        "notes": [
            "Vela ElectricField is not read as a VTK scalar here; it is derived from -grad(Potential) magnitude.",
            "Magnitude comparison is sign-insensitive; it targets distribution and unit scaling first.",
        ],
    }
    return report, best_candidate, best_unit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--vtk", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME
    node_csv_path = args.output_dir / NODE_CSV_NAME
    try:
        vtk_path = args.vtk or find_default_vtk(args.reference_root)
        report, best_candidate, best_unit = build_report(args.reference_root, vtk_path)
        sentaurus = read_sentaurus_scalar(args.reference_root / "sim_fields" / "0v" / "fields", "ElectricField")
        write_node_csv(node_csv_path, sentaurus, best_candidate, best_unit)
        report["node_csv"] = str(node_csv_path)
        write_json(report_path, report)
        md_path.parent.mkdir(parents=True, exist_ok=True)
        md_path.write_text(markdown_report(report))
        return 0 if report["status"] == "pass" else 1
    except Exception as exc:
        report = {"status": "error", "error": str(exc)}
        write_json(report_path, report)
        md_path.parent.mkdir(parents=True, exist_ok=True)
        md_path.write_text(markdown_report(report))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
