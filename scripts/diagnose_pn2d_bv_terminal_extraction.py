#!/usr/bin/env python3
"""Run PN2D BV terminal-current extraction diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_DECK = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "vela" / "simulation_bv.json"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
DEFAULT_OUT_DIR = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reports" / "terminal_extraction"
TARGET_BIASES = [-18.24, -19.24, -20.0]


METHOD_COLUMNS = {
    "sgflux": "I_sgflux_A_per_um",
    "residual": "I_residual_A_per_um",
    "sgflux_with_qf_floor": "I_sgflux_with_qf_floor_A_per_um",
    "sg_avalanche_source_integral_total": "sg_avalanche_source_integral_total",
}


def float_or_none(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def absolutize(value: str, base: Path) -> str:
    p = Path(value)
    if p.is_absolute():
        return str(p)
    return str((base / p).resolve())


def prepare_deck(base_deck: Path, out_dir: Path, start: float, stop: float, step: float,
                 contacts: list[str]) -> Path:
    base_deck = base_deck.resolve()
    deck_dir = base_deck.parent
    data = json.loads(base_deck.read_text(encoding="utf-8"))
    for key in ["mesh_file", "node_doping_file", "materials_file"]:
        if key in data and isinstance(data[key], str):
            data[key] = absolutize(data[key], deck_dir)

    curve_csv = out_dir / "terminal_extraction_curve.csv"
    compare_csv = out_dir / "terminal_current_method_compare.csv"
    data["output_csv"] = str(curve_csv.resolve())
    sweep = dict(data.get("sweep", {}))
    sweep.update({
        "mode": "bv_reverse",
        "start": start,
        "stop": stop,
        "step": step,
        "max_step": abs(step),
        "min_step": min(abs(step), 1.0e-10),
        "write_vtk": False,
        "contact": sweep.get("contact", "Anode"),
        "current_contact": sweep.get("current_contact", "Anode"),
    })
    diagnostics = dict(sweep.get("diagnostics", {}))
    diagnostics["terminal_current_method_compare"] = {
        "enabled": True,
        "contacts": contacts,
        "csv_file": str(compare_csv.resolve()),
    }
    diagnostics["contact_current_qf_floor"] = {
        "enabled": True,
        "contacts": contacts,
    }
    sweep["diagnostics"] = diagnostics
    data["sweep"] = sweep

    out_dir.mkdir(parents=True, exist_ok=True)
    deck_path = out_dir / "simulation_terminal_extraction.json"
    deck_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return deck_path


def points_for_contact(rows: list[dict[str, str]], contact: str, column: str) -> list[tuple[float, float]]:
    points: dict[float, float] = {}
    for row in rows:
        if row.get("contact") != contact:
            continue
        bias = float_or_none(row.get("bias_V"))
        value = float_or_none(row.get(column))
        if bias is None or value is None or value == 0.0:
            continue
        points[bias] = value
    return sorted(points.items())


def value_at(points: list[tuple[float, float]], bias: float, log_abs: bool = True) -> float | None:
    if not points:
        return None
    ordered = sorted(points)
    for b, value in ordered:
        if abs(b - bias) <= 1.0e-9:
            return value
    for (b0, v0), (b1, v1) in zip(ordered, ordered[1:]):
        lo = min(b0, b1)
        hi = max(b0, b1)
        if lo <= bias <= hi and b0 != b1:
            t = (bias - b0) / (b1 - b0)
            if log_abs and v0 != 0.0 and v1 != 0.0:
                sign = 1.0 if abs(v1) >= abs(v0) else (1.0 if v0 > 0.0 else -1.0)
                logv = (1.0 - t) * math.log10(abs(v0)) + t * math.log10(abs(v1))
                return sign * (10.0 ** logv)
            return (1.0 - t) * v0 + t * v1
    return None


def safe_ratio(a: float | None, b: float | None) -> float | None:
    if a is None or b is None or b == 0.0:
        return None
    return abs(a) / abs(b)


def classify(rows: list[dict[str, Any]]) -> str:
    by_method = {row["method"]: row for row in rows if row.get("bias_from_V") == TARGET_BIASES[0] and row.get("bias_to_V") == TARGET_BIASES[-1]}
    sg = by_method.get("sgflux", {}).get("growth_ratio")
    residual = by_method.get("residual", {}).get("growth_ratio")
    if residual is not None and residual >= 1.5 and sg is not None and 0.8 <= sg <= 1.2:
        return "extraction_artifact"
    return "physics_magnitude_gap"


def build_growth_rows(compare_rows: list[dict[str, str]], contact: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    value_rows: list[dict[str, Any]] = []
    growth_rows: list[dict[str, Any]] = []
    for method, column in METHOD_COLUMNS.items():
        points = points_for_contact(compare_rows, contact, column)
        values = {bias: value_at(points, bias, log_abs=(method != "sg_avalanche_source_integral_total")) for bias in TARGET_BIASES}
        for bias in TARGET_BIASES:
            value_rows.append({"bias_V": bias, "contact": contact, "method": method, "value": values[bias]})
        for left, right in [(TARGET_BIASES[0], TARGET_BIASES[1]), (TARGET_BIASES[1], TARGET_BIASES[2]), (TARGET_BIASES[0], TARGET_BIASES[2])]:
            growth_rows.append({
                "contact": contact,
                "method": method,
                "bias_from_V": left,
                "bias_to_V": right,
                "growth_ratio": safe_ratio(values[right], values[left]),
            })
    return value_rows, growth_rows


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# PN2D BV terminal extraction diagnostic",
        "",
        f"classification: `{payload['classification']}`",
        "",
        "## Growth ratios",
        "",
        "| method | from V | to V | growth |",
        "| --- | ---: | ---: | ---: |",
    ]
    for row in payload["growth_rows"]:
        lines.append(
            f"| {row['method']} | {row['bias_from_V']} | {row['bias_to_V']} | {row['growth_ratio']} |"
        )
    lines += [
        "",
        "## Artifacts",
        "",
        f"- compare_csv: `{payload['compare_csv']}`",
        f"- curve_csv: `{payload['curve_csv']}`",
        f"- deck: `{payload['deck']}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def self_test() -> int:
    rows = []
    for bias, sg, residual, source in [
        (-18.24, 1.0, 1.0, 1.0),
        (-19.24, 1.01, 1.6, 1.2),
        (-20.0, 1.02, 2.4, 1.4),
    ]:
        rows.append({
            "bias_V": str(bias),
            "contact": "Anode",
            "I_sgflux_A_per_um": str(sg),
            "I_residual_A_per_um": str(residual),
            "I_sgflux_with_qf_floor_A_per_um": str(sg),
            "sg_avalanche_source_integral_total": str(source),
        })
    _, growth = build_growth_rows(rows, "Anode")
    assert classify(growth) == "extraction_artifact"
    print("self-test passed")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--base-deck", type=Path, default=DEFAULT_DECK)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--start", type=float, default=0.0)
    parser.add_argument("--stop", type=float, default=-20.0)
    parser.add_argument("--step", type=float, default=-0.25)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--contacts", default="Anode,Cathode")
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        return self_test()
    out_dir = args.out_dir.resolve()
    contacts = [item.strip() for item in args.contacts.split(",") if item.strip()]
    deck = prepare_deck(args.base_deck, out_dir, args.start, args.stop, args.step, contacts)
    if not args.skip_run:
        subprocess.run([str(args.runner.resolve()), "--config", str(deck)], check=True)
    compare_csv = out_dir / "terminal_current_method_compare.csv"
    compare_rows = load_csv(compare_csv)
    value_rows, growth_rows = build_growth_rows(compare_rows, args.contact)
    classification = classify(growth_rows)
    payload = {
        "classification": classification,
        "contact": args.contact,
        "target_biases_V": TARGET_BIASES,
        "deck": str(deck),
        "curve_csv": str((out_dir / "terminal_extraction_curve.csv").resolve()),
        "compare_csv": str(compare_csv.resolve()),
        "value_rows": value_rows,
        "growth_rows": growth_rows,
    }
    write_csv(out_dir / "terminal_extraction_values.csv", value_rows)
    write_csv(out_dir / "terminal_extraction_growth.csv", growth_rows)
    (out_dir / "terminal_extraction_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_markdown(out_dir / "terminal_extraction_report.md", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
