#!/usr/bin/env python3
"""Follow-up diagnostics for PN2D BV reference coverage and M-1 consistency."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_IONIZATION_DIR = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reports" / "ionization_magnitude"
DEFAULT_OUT_DIR = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reports" / "ionization_followup"
DEFAULT_REFERENCE_CURVE = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reference_curves" / "pn2d_sentaurus2018_bv_reference.csv"
TARGET_BIASES = [round(-16.0 - 0.25 * i, 12) for i in range(17)]
Q = 1.602176634e-19


FIELD_ALIASES = {
    "ElectrostaticPotential": ["ElectrostaticPotential", "Potential"],
    "ElectricField": ["ElectricField"],
    "eDensity": ["eDensity", "Electrons"],
    "hDensity": ["hDensity", "Holes"],
    "eCurrent": ["eCurrent", "eCurrentDensity"],
    "hCurrent": ["hCurrent", "hCurrentDensity"],
    "TotalCurrent": ["TotalCurrent", "TotalCurrentDensity"],
    "ContactCurrentFlux": ["ContactCurrentFlux"],
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
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def nearest_row(rows: list[dict[str, Any]], bias: float, tolerance: float = 1.0e-6) -> dict[str, Any] | None:
    best = None
    best_delta = float("inf")
    for row in rows:
        b = float_or_none(str(row.get("bias_V")))
        if b is None:
            continue
        delta = abs(b - bias)
        if delta < best_delta:
            best = row
            best_delta = delta
    return best if best is not None and best_delta <= tolerance else None


def curve_points(path: Path) -> list[tuple[float, float]]:
    points: dict[float, float] = {}
    for row in load_csv(path):
        bias = float_or_none(row.get("bias_V"))
        value = float_or_none(row.get("current_total_A_per_um"))
        if value is None:
            value = float_or_none(row.get("current_total"))
        if bias is None or value is None or value == 0.0:
            continue
        previous = points.get(bias)
        if previous is None or abs(value) > abs(previous):
            points[bias] = value
    return sorted(points.items())


def value_at(points: list[tuple[float, float]], bias: float, log_abs: bool = True) -> float | None:
    ordered = sorted(points)
    for b, value in ordered:
        if abs(b - bias) <= 1.0e-9:
            return value
    for (b0, v0), (b1, v1) in zip(ordered, ordered[1:]):
        if min(b0, b1) <= bias <= max(b0, b1) and b0 != b1:
            t = (bias - b0) / (b1 - b0)
            if log_abs and v0 != 0.0 and v1 != 0.0:
                sign = 1.0 if v1 >= 0.0 else -1.0
                logv = (1.0 - t) * math.log10(abs(v0)) + t * math.log10(abs(v1))
                return sign * (10.0 ** logv)
            return (1.0 - t) * v0 + t * v1
    return None


def safe_ratio(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or den == 0.0:
        return None
    return abs(num) / abs(den)


def growth_minus_one(points: list[tuple[float, float]], left: float, right: float) -> float | None:
    left_value = value_at(points, left)
    right_value = value_at(points, right)
    if left_value is None or right_value is None or left_value == 0.0:
        return None
    return abs(right_value) / abs(left_value) - 1.0


def field_names_from_dir(fields_dir: Path) -> set[str]:
    names: set[str] = set()
    if fields_dir.exists():
        for path in fields_dir.glob("*.csv"):
            name = path.stem
            if "_region" in name:
                name = name.split("_region", 1)[0]
            names.add(name)
    return names


def metadata_bias(metadata: dict[str, Any]) -> float | None:
    values: list[float] = []
    for field in metadata.get("fields", []):
        if field.get("name") == "ContactExternalVoltage":
            raw = field.get("raw_values") or []
            if raw:
                values.append(float(raw[0]))
    negatives = [value for value in values if value < 0.0]
    if negatives:
        return min(negatives)
    return min(values) if values else None


def is_bv_field_export(metadata_path: Path, metadata: dict[str, Any], bias: float | None) -> bool:
    text = (str(metadata_path).lower() + " " + str(metadata.get("source", "")).lower())
    if "pn2d" not in text:
        return False
    if "bv" in text:
        return True
    return bias is not None and bias < 0.0 and "iv" not in text


def canonical_field_presence(names: set[str]) -> dict[str, bool]:
    result: dict[str, bool] = {}
    for canonical, aliases in FIELD_ALIASES.items():
        result[canonical] = any(alias in names for alias in aliases)
    return result


def scan_reference_coverage(search_root: Path, reference_curve: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    exports: list[dict[str, Any]] = []
    raw_files: list[dict[str, Any]] = []
    for metadata_path in search_root.rglob("metadata.json"):
        if ".git" in metadata_path.parts:
            continue
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        bias = metadata_bias(metadata)
        if not is_bv_field_export(metadata_path, metadata, bias):
            continue
        fields_dir = metadata_path.parent / "fields"
        names = field_names_from_dir(fields_dir)
        presence = canonical_field_presence(names)
        exports.append({
            "bias_V": bias,
            "metadata": str(metadata_path),
            "fields_dir": str(fields_dir) if fields_dir.exists() else None,
            "source": metadata.get("source"),
            "field_names": sorted(names),
            **{f"has_{key}": value for key, value in presence.items()},
        })
    for root, dirs, files in os.walk(search_root):
        dirs[:] = [d for d in dirs if d != ".git"]
        for name in files:
            lower = name.lower()
            if not lower.endswith((".tdr", ".plt", ".csv")):
                continue
            if "pn2d" not in lower and "sentaurus2018" not in lower:
                continue
            if "bv" not in lower and "multibias" not in lower and lower != "pn2d_sentaurus2018_bv_reference.csv":
                continue
            path = Path(root) / name
            try:
                size = path.stat().st_size
            except OSError:
                size = None
            raw_files.append({"path": str(path), "suffix": path.suffix.lower(), "size_bytes": size})
    current_points = []
    if reference_curve.exists():
        for bias, current in curve_points(reference_curve):
            if min(TARGET_BIASES) - 1.0e-9 <= bias <= max(TARGET_BIASES) + 1.0e-9:
                current_points.append({"bias_V": bias, "current_total_A_per_um": current})
    rows_by_bias: dict[float, dict[str, Any]] = {}
    for export in exports:
        bias = export.get("bias_V")
        if bias is None:
            continue
        rounded = round(float(bias), 9)
        existing = rows_by_bias.get(rounded)
        if existing is None:
            rows_by_bias[rounded] = export
        else:
            existing_fields = len(existing.get("field_names", []))
            export_fields = len(export.get("field_names", []))
            if export_fields > existing_fields:
                rows_by_bias[rounded] = export
    coverage_rows = []
    for bias in TARGET_BIASES:
        export = rows_by_bias.get(round(bias, 9))
        current = value_at([(row["bias_V"], row["current_total_A_per_um"]) for row in current_points], bias)
        row = {
            "bias_V": bias,
            "has_sentaurus_field_export": export is not None,
            "has_reference_current": current is not None,
            "reference_current_total_A_per_um": current,
        }
        if export is not None:
            row.update(export)
        coverage_rows.append(row)
    summary = {
        "field_export_count": len(exports),
        "raw_file_count": len(raw_files),
        "reference_curve": str(reference_curve),
        "available_field_biases_V": [row["bias_V"] for row in coverage_rows if row["has_sentaurus_field_export"]],
        "missing_field_biases_V": [row["bias_V"] for row in coverage_rows if not row["has_sentaurus_field_export"]],
        "current_bias_count_in_window": len(current_points),
        "only_minus20_field_available": [row["bias_V"] for row in coverage_rows if row["has_sentaurus_field_export"]] == [-20.0],
        "raw_files": raw_files,
        "all_field_exports": exports,
    }
    return coverage_rows, summary


def classify_m_consistency(script_m: float | None, terminal_m: float | None, source_m: float | None) -> str:
    if script_m is not None and terminal_m is not None:
        scale = max(abs(terminal_m), abs(source_m or 0.0), 1.0e-30)
        if abs(script_m) > 3.0 * scale:
            return "script_integral_overestimate"
    if source_m is not None and terminal_m is not None:
        if abs(source_m) > 3.0 * max(abs(terminal_m), 1.0e-30):
            return "physical_M_not_reaching_contact"
    return "consistent_low_M"


def load_terminal_contact_points(path: Path, contact: str, column: str) -> list[tuple[float, float]]:
    points: dict[float, float] = {}
    if not path.exists():
        return []
    for row in load_csv(path):
        if row.get("contact") != contact:
            continue
        bias = float_or_none(row.get("bias_V"))
        value = float_or_none(row.get(column))
        if bias is None or value is None or value == 0.0:
            continue
        points[bias] = value
    return sorted(points.items())


def load_source_points(path: Path) -> list[tuple[float, float]]:
    if not path.exists():
        return []
    rows: dict[float, float] = {}
    for row in load_csv(path):
        bias = float_or_none(row.get("bias_V"))
        source = float_or_none(row.get("sg_avalanche_source_integral_total"))
        if bias is None or source is None:
            continue
        rows[bias] = source
    return sorted(rows.items())


def build_m_consistency(ionization_dir: Path, reference_curve: Path, contact: str) -> dict[str, Any]:
    report_json = json.loads((ionization_dir / "ionization_magnitude_report.json").read_text(encoding="utf-8"))
    points_rows = report_json.get("rows", [])
    focus_bias = -20.0
    low_seed_bias = -16.0
    knee_seed_bias = -18.2432189285
    focus_row = nearest_row(points_rows, focus_bias, tolerance=1.0e-6) or {}
    script_integral = float_or_none(str(focus_row.get("alpha_integral_vela_field")))
    script_m = float_or_none(str(focus_row.get("M_minus_1_vela_field_estimate")))
    script_m_source = "reported_M_minus_1_vela_field_estimate"
    if script_m is None and script_integral is not None:
        script_m = math.expm1(min(script_integral, 700.0))
        script_m_source = "reconstructed_from_alpha_integral_for_audit_only"

    terminal_points = load_terminal_contact_points(
        ionization_dir / "terminal_current_method_compare.csv", contact, "I_sgflux_A_per_um")
    residual_points = load_terminal_contact_points(
        ionization_dir / "terminal_current_method_compare.csv", contact, "I_residual_A_per_um")
    source_points_raw = load_source_points(ionization_dir / "terminal_current_method_compare.csv")
    source_current_points = [(bias, source * Q * 1.0e-6) for bias, source in source_points_raw]
    reference_points = curve_points(reference_curve)

    i_low = abs(value_at(terminal_points, low_seed_bias) or 0.0)
    i_focus = abs(value_at(terminal_points, focus_bias) or 0.0)
    i_knee = abs(value_at(terminal_points, knee_seed_bias) or 0.0)
    source_focus = abs(value_at(source_current_points, focus_bias, log_abs=False) or 0.0)
    source_low = abs(value_at(source_current_points, low_seed_bias, log_abs=False) or 0.0)
    source_knee = abs(value_at(source_current_points, knee_seed_bias, log_abs=False) or 0.0)

    terminal_m_low_to_focus = None if i_low == 0.0 else i_focus / i_low - 1.0
    terminal_m_knee_to_focus = None if i_knee == 0.0 else i_focus / i_knee - 1.0
    source_m_total_vs_low = None if i_low == 0.0 else source_focus / i_low
    source_m_delta_vs_low = None if i_low == 0.0 else (source_focus - source_low) / i_low
    source_m_delta_vs_knee = None if i_knee == 0.0 else (source_focus - source_knee) / i_knee
    residual_delta = None
    residual_low = abs(value_at(residual_points, low_seed_bias) or 0.0)
    residual_focus = abs(value_at(residual_points, focus_bias) or 0.0)
    if residual_low != 0.0:
        residual_delta = residual_focus / residual_low - 1.0

    reference_m_low_to_focus = growth_minus_one(reference_points, low_seed_bias, focus_bias)
    reference_m_knee_to_focus = growth_minus_one(reference_points, knee_seed_bias, focus_bias)
    classification = classify_m_consistency(script_m, terminal_m_knee_to_focus, source_m_total_vs_low)
    audit = {
        "integration_path": "junction-normal line band through the x=1um metallurgical junction, not a carrier streamline",
        "carrier_counting": "uses alpha_n + alpha_p, which double-counts for a single-carrier multiplication estimate",
        "driving_field_used_by_script": "|grad psi| from electrostatic potential profile",
        "deck_driving_force": "quasi_fermi_gradient / current_density edge-current source",
        "formula_used": "M_minus_1 = exp(integral(alpha_n + alpha_p) dl) - 1",
        "verdict": "the line integral is reconstructed only as the audited legacy scalar for this deck; it is not a terminal-current M-1 estimate",
    }
    continuity_balance_path = ionization_dir / "continuity_balance.csv"
    continuity_balance_rows = load_csv(continuity_balance_path) if continuity_balance_path.exists() else []
    continuity_note = (
        f"present with {len(continuity_balance_rows)} rows; terminal_current_method_compare residual current is used for the scalar terminal-current M-1 cross-check"
        if continuity_balance_path.exists()
        else "not present in the existing ionization_magnitude run; terminal_current_method_compare residual current is used as the continuity-conservative terminal check"
    )
    return {
        "classification": classification,
        "audit": audit,
        "focus_bias_V": focus_bias,
        "low_seed_bias_V": low_seed_bias,
        "knee_seed_bias_V": knee_seed_bias,
        "script_line_integral": script_integral,
        "script_line_integral_M_minus_1": script_m,
        "script_line_integral_M_minus_1_source": script_m_source,
        "terminal_current_M_minus_1_low_to_focus": terminal_m_low_to_focus,
        "terminal_current_M_minus_1_knee_to_focus": terminal_m_knee_to_focus,
        "terminal_residual_M_minus_1_low_to_focus": residual_delta,
        "source_current_A_per_um_at_focus": source_focus,
        "source_current_A_per_um_at_low_seed": source_low,
        "source_current_A_per_um_at_knee_seed": source_knee,
        "source_total_M_minus_1_eff_vs_low_seed": source_m_total_vs_low,
        "source_delta_M_minus_1_eff_vs_low_seed": source_m_delta_vs_low,
        "source_delta_M_minus_1_eff_vs_knee_seed": source_m_delta_vs_knee,
        "reference_current_M_minus_1_low_to_focus": reference_m_low_to_focus,
        "reference_current_M_minus_1_knee_to_focus": reference_m_knee_to_focus,
        "terminal_current_A_per_um_at_low_seed": i_low,
        "terminal_current_A_per_um_at_knee_seed": i_knee,
        "terminal_current_A_per_um_at_focus": i_focus,
        "continuity_balance_csv": str(continuity_balance_path) if continuity_balance_path.exists() else None,
        "continuity_balance_rows": len(continuity_balance_rows),
        "continuity_balance_note": continuity_note,
    }


def write_coverage_markdown(path: Path, rows: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    lines = [
        "# PN2D BV Sentaurus Reference Coverage",
        "",
        f"Only -20V field available: `{summary['only_minus20_field_available']}`",
        "",
        "| bias V | field export | reference current | potential | eDensity | hDensity | current field | fields dir |",
        "| ---: | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        current_field = bool(row.get("has_eCurrent") or row.get("has_hCurrent") or row.get("has_TotalCurrent") or row.get("has_ContactCurrentFlux"))
        lines.append(
            f"| {row['bias_V']} | {row['has_sentaurus_field_export']} | {row['has_reference_current']} | "
            f"{row.get('has_ElectrostaticPotential')} | {row.get('has_eDensity')} | {row.get('has_hDensity')} | "
            f"{current_field} | {row.get('fields_dir')} |"
        )
    lines += [
        "",
        "## Missing Field Biases",
        "",
        f"`{summary['missing_field_biases_V']}`",
        "",
        "## Raw File Inventory",
        "",
        f"Raw files matching pn2d BV/multibias/reference criteria: `{summary['raw_file_count']}`",
    ]
    for item in summary["raw_files"][:80]:
        lines.append(f"- `{item['path']}` ({item['suffix']}, {item['size_bytes']} bytes)")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_m_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# PN2D BV M-1 Consistency",
        "",
        f"classification: `{payload['classification']}`",
        "",
        "## Audit",
        "",
    ]
    for key, value in payload["audit"].items():
        lines.append(f"- {key}: {value}")
    lines += [
        "",
        "## Numeric Cross-Checks",
        "",
        "| metric | value |",
        "| --- | ---: |",
    ]
    for key in [
        "script_line_integral",
        "script_line_integral_M_minus_1",
        "terminal_current_M_minus_1_low_to_focus",
        "terminal_current_M_minus_1_knee_to_focus",
        "terminal_residual_M_minus_1_low_to_focus",
        "source_total_M_minus_1_eff_vs_low_seed",
        "source_delta_M_minus_1_eff_vs_low_seed",
        "source_delta_M_minus_1_eff_vs_knee_seed",
        "reference_current_M_minus_1_low_to_focus",
        "reference_current_M_minus_1_knee_to_focus",
        "source_current_A_per_um_at_focus",
        "terminal_current_A_per_um_at_focus",
    ]:
        lines.append(f"| {key} | {payload.get(key)} |")
    lines += [
        "",
        "## Continuity Balance",
        "",
        f"continuity_balance_csv: `{payload.get('continuity_balance_csv')}`",
        "",
        payload.get("continuity_balance_note", ""),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def self_test() -> int:
    assert classify_m_consistency(6.9, 0.5, 0.5) == "script_integral_overestimate"
    assert classify_m_consistency(0.6, 0.1, 1.0) == "physical_M_not_reaching_contact"
    assert classify_m_consistency(0.4, 0.35, 0.3) == "consistent_low_M"
    assert safe_ratio(2.0, 4.0) == 0.5
    print("self-test passed")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--search-root", type=Path, default=REPO)
    parser.add_argument("--ionization-dir", type=Path, default=DEFAULT_IONIZATION_DIR)
    parser.add_argument("--reference-curve", type=Path, default=DEFAULT_REFERENCE_CURVE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--contact", default="Anode")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        return self_test()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    coverage_rows, coverage_summary = scan_reference_coverage(args.search_root.resolve(), args.reference_curve.resolve())
    coverage_payload = {"rows": coverage_rows, "summary": coverage_summary}
    coverage_json = out_dir / "reference_field_coverage.json"
    coverage_md = out_dir / "reference_field_coverage.md"
    coverage_csv = out_dir / "reference_field_coverage.csv"
    coverage_json.write_text(json.dumps(coverage_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(coverage_csv, coverage_rows)
    write_coverage_markdown(coverage_md, coverage_rows, coverage_summary)

    consistency = build_m_consistency(args.ionization_dir.resolve(), args.reference_curve.resolve(), args.contact)
    consistency_json = out_dir / "m_minus_1_consistency.json"
    consistency_md = out_dir / "m_minus_1_consistency.md"
    consistency_json.write_text(json.dumps(consistency, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_m_markdown(consistency_md, consistency)

    print(json.dumps({
        "reference_field_coverage_json": str(coverage_json),
        "reference_field_coverage_md": str(coverage_md),
        "m_minus_1_consistency_json": str(consistency_json),
        "m_minus_1_consistency_md": str(consistency_md),
        "only_minus20_field_available": coverage_summary["only_minus20_field_available"],
        "m_minus_1_classification": consistency["classification"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())