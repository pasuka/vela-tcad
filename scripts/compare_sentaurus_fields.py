#!/usr/bin/env python3
"""Compare exported Sentaurus/Vela node field CSVs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def parse_field_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def field_files(root: Path, field: str) -> list[Path]:
    exact = root / f"{field}_region0.csv"
    if exact.is_file():
        return [exact]
    return sorted(root.glob(f"{field}_region*.csv"))


def load_field_csv(path: Path) -> dict[int, list[float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "node_id" not in reader.fieldnames:
            raise ValueError(f"{path} is missing node_id column")
        component_columns = sorted(
            (name for name in reader.fieldnames if name.startswith("component")),
            key=lambda name: int(name.removeprefix("component") or "0"),
        )
        if not component_columns:
            component_columns = [name for name in reader.fieldnames if name != "node_id"]
        values: dict[int, list[float]] = {}
        for row in reader:
            node_id = int(float(row["node_id"]))
            values[node_id] = [float(row[name]) for name in component_columns]
        return values


def compare_values(reference: dict[int, list[float]],
                   candidate: dict[int, list[float]]) -> dict[str, object]:
    common_nodes = sorted(set(reference) & set(candidate))
    max_abs = 0.0
    max_rel = 0.0
    points = 0
    worst: dict[str, object] | None = None
    for node_id in common_nodes:
        ref_components = reference[node_id]
        cand_components = candidate[node_id]
        for component, (ref_value, cand_value) in enumerate(zip(ref_components, cand_components)):
            abs_diff = abs(cand_value - ref_value)
            rel_diff = abs_diff / max(abs(ref_value), 1.0e-300)
            points += 1
            if abs_diff > max_abs:
                max_abs = abs_diff
                worst = {
                    "node_id": node_id,
                    "component": component,
                    "reference": ref_value,
                    "candidate": cand_value,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_diff,
                }
            max_rel = max(max_rel, rel_diff)
    return {
        "points_compared": points,
        "reference_points": len(reference),
        "candidate_points": len(candidate),
        "missing_reference_nodes": len(set(candidate) - set(reference)),
        "missing_candidate_nodes": len(set(reference) - set(candidate)),
        "max_abs_diff": max_abs,
        "max_rel_diff": max_rel,
        "worst": worst,
    }


def compare_field(reference_root: Path,
                  candidate_root: Path,
                  field: str) -> dict[str, object]:
    reference_files = field_files(reference_root, field)
    candidate_files = field_files(candidate_root, field)
    if not reference_files:
        return {"status": "missing_reference", "points_compared": 0}
    if not candidate_files:
        return {"status": "missing_candidate", "points_compared": 0}
    by_region: dict[str, object] = {}
    merged_reference: dict[int, list[float]] = {}
    merged_candidate: dict[int, list[float]] = {}
    for reference_file in reference_files:
        suffix = reference_file.name.removeprefix(field + "_")
        candidate_file = candidate_root / f"{field}_{suffix}"
        if not candidate_file.is_file():
            by_region[suffix] = {"status": "missing_candidate", "points_compared": 0}
            continue
        reference_values = load_field_csv(reference_file)
        candidate_values = load_field_csv(candidate_file)
        region_report = compare_values(reference_values, candidate_values)
        region_report["status"] = "pass"
        by_region[suffix] = region_report
        merged_reference.update(reference_values)
        merged_candidate.update(candidate_values)
    merged_report = compare_values(merged_reference, merged_candidate)
    merged_report["regions"] = by_region
    merged_report["status"] = "pass"
    return merged_report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-fields", required=True, type=Path)
    parser.add_argument("--candidate-fields", required=True, type=Path)
    parser.add_argument("--fields", required=True, help="Comma-separated field names")
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--max-abs-diff", type=float, default=math.inf)
    parser.add_argument("--max-rel-diff", type=float, default=math.inf)
    args = parser.parse_args()

    report: dict[str, object] = {"status": "pass", "fields": {}}
    failed = False
    for field in parse_field_list(args.fields):
        field_report = compare_field(args.reference_fields, args.candidate_fields, field)
        status = field_report.get("status")
        if status != "pass":
            failed = True
        elif (
            float(field_report.get("max_abs_diff", 0.0)) > args.max_abs_diff
            or float(field_report.get("max_rel_diff", 0.0)) > args.max_rel_diff
            or int(field_report.get("missing_candidate_nodes", 0)) > 0
            or int(field_report.get("missing_reference_nodes", 0)) > 0
        ):
            field_report["status"] = "fail"
            failed = True
        report["fields"][field] = field_report
    if failed:
        report["status"] = "fail"
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
