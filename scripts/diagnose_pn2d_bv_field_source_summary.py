#!/usr/bin/env python3
"""Summarize PN2D BV electric-field and avalanche/source parity."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


FIELDS = [
    "bias_V",
    "quantity",
    "points",
    "sentaurus_max",
    "vela_max",
    "log10_vela_over_sentaurus_max",
    "sentaurus_p95",
    "vela_p95",
    "log10_vela_over_sentaurus_p95",
    "sentaurus_sum_proxy",
    "vela_sum_proxy",
    "log10_vela_over_sentaurus_sum_proxy",
    "field_compare_status",
    "classification",
]


QUANTITIES = {
    "electric_field": {
        "sentaurus": ["ElectricField"],
        "vela": ["ElectricField"],
        "vela_scale": 1.0,
    },
    "avalanche_generation": {
        "sentaurus": ["ImpactIonization", "AvalancheGeneration"],
        "vela": ["AvalancheGeneration", "ImpactIonization"],
        "vela_scale": 1.0e-6,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--field-compare", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", required=True)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def discover_sentaurus_dir(root: Path, bias: float) -> Path | None:
    for candidate in [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias:g}v",
    ]:
        if (candidate / "nodes.csv").exists():
            return candidate
    return None


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    result: dict[float, Path] = {}
    mtimes: dict[float, float] = {}
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if not match:
            continue
        key = bias_key(float(match.group("bias")))
        mtime = path.stat().st_mtime
        if key not in result or mtime >= mtimes[key]:
            result[key] = path
            mtimes[key] = mtime
    return result


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_scalar_csv(export_dir: Path, names: list[str]) -> list[float]:
    for name in names:
        path = export_dir / "fields" / f"{name}_region0.csv"
        if not path.exists():
            continue
        values: list[float] = []
        for row in read_csv_rows(path):
            raw = row.get("component0")
            if raw not in (None, ""):
                values.append(float(raw))
        return values
    return []


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    npoints = 0
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "POINTS":
            npoints = int(parts[1])
            break
    fields: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < npoints:
                text = lines[i].strip()
                if not text:
                    i += 1
                    continue
                head = text.split()[0]
                if head in {"SCALARS", "VECTORS", "CELL_DATA", "POINT_DATA"}:
                    break
                values.extend(float(item) for item in text.split())
                i += 1
            fields[name] = values[:npoints]
            continue
        i += 1
    return fields


def pick_vtk_field(fields: dict[str, list[float]], names: list[str]) -> list[float]:
    for name in names:
        if name in fields:
            return fields[name]
    return []


def scaled(values: list[float], factor: float) -> list[float]:
    return [value * factor for value in values]


def pct(values: list[float], percentile: float) -> float | None:
    clean = sorted(abs(value) for value in values if math.isfinite(value))
    if not clean:
        return None
    idx = (len(clean) - 1) * percentile / 100.0
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return clean[lo]
    return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)


def log10_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate == 0.0 or reference == 0.0:
        return None
    return math.log10(abs(candidate) / abs(reference))


def stats(values: list[float]) -> dict[str, float | None]:
    clean = [abs(value) for value in values if math.isfinite(value)]
    if not clean:
        return {"points": 0, "max": None, "p95": None, "sum_proxy": None}
    return {
        "points": len(clean),
        "max": max(clean),
        "p95": pct(clean, 95),
        "sum_proxy": sum(clean),
    }


def load_field_compare_status(path: Path) -> dict[tuple[float, str], str]:
    result: dict[tuple[float, str], str] = {}
    for row in read_csv_rows(path):
        try:
            key = (bias_key(float(row["bias_V"])), row["quantity"])
        except (KeyError, ValueError):
            continue
        result[key] = row.get("avalanche_status") or row.get("status") or ""
    return result


def classify(quantity: str, field_status: str, log_sum: float | None) -> str:
    if quantity == "avalanche_generation" and field_status:
        return "threshold_sensitive_source_check"
    if log_sum is not None and abs(log_sum) > 0.3:
        return "source_or_field_sum_mismatch"
    return "field_source_summary_available"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    biases = parse_biases(args.biases)
    vela_vtks = discover_vela_vtks(args.vela_vtk_root)
    field_status = load_field_compare_status(args.field_compare)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        sentaurus_dir = discover_sentaurus_dir(args.sentaurus_root, bias)
        vtk_path = vela_vtks.get(bias_key(bias))
        vtk_fields = parse_vtk_scalars(vtk_path) if vtk_path is not None else {}
        for quantity, spec in QUANTITIES.items():
            s_values = load_scalar_csv(sentaurus_dir, spec["sentaurus"]) if sentaurus_dir else []
            v_values = scaled(pick_vtk_field(vtk_fields, spec["vela"]), float(spec["vela_scale"]))
            s_stats = stats(s_values)
            v_stats = stats(v_values)
            log_sum = log10_ratio(v_stats["sum_proxy"], s_stats["sum_proxy"])
            status = field_status.get((bias_key(bias), quantity), "")
            rows.append({
                "bias_V": bias,
                "quantity": quantity,
                "points": min(int(s_stats["points"] or 0), int(v_stats["points"] or 0)),
                "sentaurus_max": s_stats["max"],
                "vela_max": v_stats["max"],
                "log10_vela_over_sentaurus_max": log10_ratio(v_stats["max"], s_stats["max"]),
                "sentaurus_p95": s_stats["p95"],
                "vela_p95": v_stats["p95"],
                "log10_vela_over_sentaurus_p95": log10_ratio(v_stats["p95"], s_stats["p95"]),
                "sentaurus_sum_proxy": s_stats["sum_proxy"],
                "vela_sum_proxy": v_stats["sum_proxy"],
                "log10_vela_over_sentaurus_sum_proxy": log_sum,
                "field_compare_status": status,
                "classification": classify(quantity, status, log_sum),
            })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "field_source_summary.csv", rows)
    (args.out_dir / "summary.json").write_text(
        json.dumps({"biases": biases, "rows": rows}, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
