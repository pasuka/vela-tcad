#!/usr/bin/env python3
"""Cross-check Sentaurus PLT terminal current against TDR ContactCurrentFlux."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sentaurus-root", type=Path, help="Root containing sentaurus_<bias>v exports.")
    group.add_argument("--export-dir", type=Path, help="One Sentaurus TDR export directory.")
    parser.add_argument("--sentaurus-plt", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--plt-bias-column", default=None)
    parser.add_argument("--plt-current-column", default=None)
    parser.add_argument("--max-relative-difference", type=float, default=0.01)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument(
        "--fail-on-mismatch",
        action="store_true",
        help="Return non-zero when PLT and ContactCurrentFlux differ beyond the tolerance.",
    )
    return parser.parse_args()


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def optional_float(raw: str | None) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_quoted_list(text: str, key: str) -> list[str]:
    match = re.search(rf"{re.escape(key)}\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def parse_plt_values_block(text: str, column_count: int) -> list[list[float]]:
    match = re.search(r"Values\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        match = re.search(r"Data\s*\{(.*?)\}", text, re.DOTALL | re.IGNORECASE)
    if not match:
        return []
    numbers = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
            match.group(1),
        )
    ]
    if column_count <= 0 or len(numbers) % column_count != 0:
        return []
    return [
        numbers[index:index + column_count]
        for index in range(0, len(numbers), column_count)
    ]


def load_plt_current(
    path: Path,
    bias: float,
    bias_column: str,
    current_column: str,
) -> float:
    if not path.exists():
        raise FileNotFoundError(f"Sentaurus PLT not found: {path}")
    text = path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if bias_column not in datasets:
        raise ValueError(f"PLT bias column not found: {bias_column}")
    if current_column not in datasets:
        raise ValueError(f"PLT current column not found: {current_column}")
    bias_index = datasets.index(bias_column)
    current_index = datasets.index(current_column)
    points: dict[float, float] = {}
    for row in parse_plt_values_block(text, len(datasets)):
        points[bias_key(row[bias_index])] = row[current_index]
    key = bias_key(bias)
    if key in points:
        return points[key]
    if not points:
        raise ValueError(f"PLT contains no numeric data rows: {path}")
    nearest = min(points, key=lambda item: abs(item - key))
    if abs(nearest - key) > 1.0e-9:
        raise ValueError(f"PLT bias not found: {bias:g}")
    return points[nearest]


def discover_export_dir(root: Path | None, export_dir: Path | None, bias: float) -> Path:
    if export_dir is not None:
        return export_dir
    assert root is not None
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "field_manifest.json").exists():
            return candidate
    raise FileNotFoundError(
        "No Sentaurus export directory with field_manifest.json found for "
        f"bias {bias:g} under {root}"
    )


def load_contact_current_flux(export_dir: Path, contact: str) -> dict[str, Any]:
    manifest_path = export_dir / "field_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"field_manifest.json not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    datasets = manifest.get("datasets") or manifest.get("fields") or (
        manifest if isinstance(manifest, list) else []
    )
    for item in datasets:
        if item.get("name") != "ContactCurrentFlux":
            continue
        if str(item.get("region_name", "")).lower() != contact.lower():
            continue
        region = item.get("region")
        path = export_dir / "fields" / f"ContactCurrentFlux_region{region}.csv"
        if not path.exists():
            raise FileNotFoundError(f"ContactCurrentFlux CSV not found: {path}")
        for row in load_csv(path):
            value = optional_float(row.get("component0"))
            if value is not None:
                return {"region": region, "value_A": value}
        raise ValueError(f"ContactCurrentFlux has no finite component0 value: {path}")
    raise ValueError(f"ContactCurrentFlux for contact {contact!r} not found in {manifest_path}")


def relative_difference(lhs: float, rhs: float) -> float:
    return abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0e-300)


def log10_ratio(lhs: float, rhs: float) -> float | None:
    if lhs == 0.0 or rhs == 0.0:
        return None
    return math.log10(abs(lhs) / abs(rhs))


def build_summary(args: argparse.Namespace) -> dict[str, Any]:
    contact = args.contact
    bias_column = args.plt_bias_column or f"{contact} OuterVoltage"
    current_column = args.plt_current_column or f"{contact} TotalCurrent"
    export_dir = discover_export_dir(args.sentaurus_root, args.export_dir, args.bias)
    flux = load_contact_current_flux(export_dir, contact)
    plt_current = load_plt_current(
        args.sentaurus_plt,
        args.bias,
        bias_column,
        current_column,
    )
    flux_current = float(flux["value_A"])
    rel = relative_difference(flux_current, plt_current)
    status = (
        "sentaurus_plt_contact_flux_consistent"
        if rel <= args.max_relative_difference
        else "sentaurus_plt_contact_flux_mismatch"
    )
    return {
        "bias_V": args.bias,
        "contact": contact,
        "export_dir": str(export_dir),
        "sentaurus_plt": str(args.sentaurus_plt),
        "plt_bias_column": bias_column,
        "plt_current_column": current_column,
        "contact_current_flux_region": flux["region"],
        "contact_current_flux_A": flux_current,
        "plt_current_A": plt_current,
        "relative_difference": rel,
        "log10_contact_flux_over_plt": log10_ratio(flux_current, plt_current),
        "max_relative_difference": args.max_relative_difference,
        "status": status,
    }


def main() -> int:
    args = parse_args()
    try:
        summary = build_summary(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    if args.fail_on_mismatch and summary["status"].endswith("_mismatch"):
        print(
            "Sentaurus PLT terminal current and TDR ContactCurrentFlux differ "
            f"by relative {summary['relative_difference']:.6g}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
