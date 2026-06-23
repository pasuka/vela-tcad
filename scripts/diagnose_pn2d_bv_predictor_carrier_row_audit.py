#!/usr/bin/env python3
"""Audit PN2D BV predictor active-endpoint carrier rows and residual terms."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag

NODE_FIELDS = [
    "state",
    "node_id",
    "support_class",
    "initial_psi_minus_phin_shift_V",
    "initial_phip_minus_psi_shift_V",
    "electron_residual",
    "hole_residual",
    "electron_diagonal",
    "hole_diagonal",
    "electron_offdiag_fraction",
    "hole_offdiag_fraction",
    "electron_raw_delta_phin_V",
    "hole_raw_delta_phip_V",
    "electron_raw_delta_rollback_fraction",
    "hole_raw_delta_rollback_fraction",
    "electron_flux",
    "electron_recombination",
    "electron_impact",
    "electron_term_sum",
    "hole_flux",
    "hole_recombination",
    "hole_impact",
    "hole_term_sum",
    "impact_combined_source",
]

IMPACT_SCALE_NODE_FIELDS = [
    "state",
    "node_id",
    "support_class",
    "impact_scale",
    "electron_impact",
    "electron_term_sum",
    "electron_adjusted_impact",
    "electron_adjusted_residual",
    "electron_required_impact_scale",
    "hole_impact",
    "hole_term_sum",
    "hole_adjusted_impact",
    "hole_adjusted_residual",
    "hole_required_impact_scale",
]

STATE_FIELDS = {
    "psi": "ElectrostaticPotential",
    "phin": "eQuasiFermiPotential",
    "phip": "hQuasiFermiPotential",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--carrier-row-csv", action="append", default=[], metavar="STATE=CSV")
    parser.add_argument("--carrier-term-csv", action="append", default=[], metavar="STATE=CSV")
    parser.add_argument("--base-config", type=Path)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--state-csv", action="append", default=[], metavar="STATE=CSV")
    parser.add_argument("--state-vtk", action="append", default=[], metavar="STATE=VTK")
    parser.add_argument("--step-csv", type=Path, help="Optional step probe CSV used to synthesize trial state.")
    parser.add_argument("--include-trial-from-step", action="store_true")
    parser.add_argument("--bias", type=float, default=-20.0)
    parser.add_argument("--bias-contact", default="Anode")
    parser.add_argument("--ground-contact", default="Cathode")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument(
        "--impact-scale",
        action="append",
        type=float,
        default=[],
        help="Impact-source multiplier to evaluate against carrier term rows.",
    )
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def safe_ratio(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or abs(den) < 1.0e-300:
        return None
    return num / den


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def mean(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def clean_json(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


def parse_named_paths(items: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for raw in items:
        if "=" not in raw:
            raise RuntimeError(f"expected STATE=PATH, got {raw!r}")
        name, value = raw.split("=", 1)
        name = name.strip()
        if not name:
            raise RuntimeError("state name must be non-empty")
        if name in result:
            raise RuntimeError(f"duplicate state name {name}")
        result[name] = Path(value)
    return result


def support_metadata(path: Path) -> dict[int, dict[str, Any]]:
    meta: dict[int, dict[str, Any]] = {}
    for row in read_csv_rows(path):
        node_id = int(row["node_id"])
        phin_base = finite_float(row.get("phin_baseline_V"))
        phin_pred = finite_float(row.get("phin_predictor_V"))
        phip_base = finite_float(row.get("phip_baseline_V"))
        phip_pred = finite_float(row.get("phip_predictor_V"))
        if phin_base is None or phin_pred is None or phip_base is None or phip_pred is None:
            continue
        meta[node_id] = {
            "support_class": row.get("support_class") or "unknown",
            "initial_psi_minus_phin_shift_V": -(phin_pred - phin_base),
            "initial_phip_minus_psi_shift_V": phip_pred - phip_base,
        }
    if not meta:
        raise RuntimeError("support csv has no predictor metadata rows")
    return meta


def row_by_node(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["node_id"]): row for row in read_csv_rows(path)}


def offdiag_fraction(row: dict[str, str], carrier: str) -> float | None:
    offdiag = finite_float(row.get(f"{carrier}_offdiag_abs_sum"))
    total = finite_float(row.get(f"{carrier}_row_abs_sum"))
    return safe_ratio(offdiag, total)


def state_node_rows(
    state: str,
    row_rows: dict[int, dict[str, str]],
    term_rows: dict[int, dict[str, str]],
    meta: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for node_id, node_meta in sorted(meta.items()):
        row = row_rows.get(node_id, {})
        term = term_rows.get(node_id, {})
        initial_e = finite_float(node_meta["initial_psi_minus_phin_shift_V"])
        initial_h = finite_float(node_meta["initial_phip_minus_psi_shift_V"])
        raw_e = finite_float(row.get("raw_delta_phin_V"))
        raw_h = finite_float(row.get("raw_delta_phip_V"))
        rows.append({
            "state": state,
            "node_id": node_id,
            "support_class": node_meta["support_class"],
            "initial_psi_minus_phin_shift_V": initial_e,
            "initial_phip_minus_psi_shift_V": initial_h,
            "electron_residual": finite_float(row.get("electron_residual")),
            "hole_residual": finite_float(row.get("hole_residual")),
            "electron_diagonal": finite_float(row.get("electron_diagonal")),
            "hole_diagonal": finite_float(row.get("hole_diagonal")),
            "electron_offdiag_fraction": offdiag_fraction(row, "electron"),
            "hole_offdiag_fraction": offdiag_fraction(row, "hole"),
            "electron_raw_delta_phin_V": raw_e,
            "hole_raw_delta_phip_V": raw_h,
            "electron_raw_delta_rollback_fraction": safe_ratio(raw_e, initial_e),
            "hole_raw_delta_rollback_fraction": safe_ratio(-raw_h if raw_h is not None else None, initial_h),
            "electron_flux": finite_float(term.get("electron_flux")),
            "electron_recombination": finite_float(term.get("electron_recombination")),
            "electron_impact": finite_float(term.get("electron_impact")),
            "electron_term_sum": finite_float(term.get("electron_term_sum")),
            "hole_flux": finite_float(term.get("hole_flux")),
            "hole_recombination": finite_float(term.get("hole_recombination")),
            "hole_impact": finite_float(term.get("hole_impact")),
            "hole_term_sum": finite_float(term.get("hole_term_sum")),
            "impact_combined_source": finite_float(term.get("impact_combined_source")),
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values = []
    for row in rows:
        value = finite_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    classes = sorted({row["support_class"] for row in rows})
    result: dict[str, Any] = {"count": len(rows), "support_classes": {}}
    metrics = [
        "electron_residual", "hole_residual",
        "electron_diagonal", "hole_diagonal",
        "electron_offdiag_fraction", "hole_offdiag_fraction",
        "electron_raw_delta_phin_V", "hole_raw_delta_phip_V",
        "electron_raw_delta_rollback_fraction", "hole_raw_delta_rollback_fraction",
        "electron_flux", "electron_recombination", "electron_impact", "electron_term_sum",
        "hole_flux", "hole_recombination", "hole_impact", "hole_term_sum",
        "impact_combined_source",
    ]
    for support_class in classes:
        subset = [row for row in rows if row["support_class"] == support_class]
        item: dict[str, Any] = {"count": len(subset)}
        for metric in metrics:
            values = finite_values(subset, metric)
            if metric.endswith("_V"):
                prefix = metric[:-2]
                item[f"{prefix}_median_V"] = median(values)
                item[f"{prefix}_mean_V"] = mean(values)
            else:
                item[f"{metric}_median"] = median(values)
                item[f"{metric}_mean"] = mean(values)
            if metric.endswith("residual"):
                item[f"{metric}_abs_median"] = median([abs(v) for v in values])
        result["support_classes"][support_class] = item
    return result


def impact_scale_label(scale: float) -> str:
    return f"{scale:.12g}"


def adjusted_residual(term_sum: float | None, impact: float | None, scale: float) -> float | None:
    if term_sum is None or impact is None:
        return None
    return term_sum + (scale - 1.0) * impact


def required_impact_scale(term_sum: float | None, impact: float | None) -> float | None:
    if term_sum is None or impact is None or abs(impact) < 1.0e-300:
        return None
    return 1.0 - term_sum / impact


def impact_scale_rows(rows: list[dict[str, Any]], scales: list[float]) -> list[dict[str, Any]]:
    unique_scales = sorted({scale for scale in scales if math.isfinite(scale)})
    result: list[dict[str, Any]] = []
    for row in rows:
        electron_term = finite_float(row.get("electron_term_sum"))
        electron_impact = finite_float(row.get("electron_impact"))
        hole_term = finite_float(row.get("hole_term_sum"))
        hole_impact = finite_float(row.get("hole_impact"))
        electron_required = required_impact_scale(electron_term, electron_impact)
        hole_required = required_impact_scale(hole_term, hole_impact)
        for scale in unique_scales:
            result.append({
                "state": row["state"],
                "node_id": row["node_id"],
                "support_class": row["support_class"],
                "impact_scale": impact_scale_label(scale),
                "electron_impact": electron_impact,
                "electron_term_sum": electron_term,
                "electron_adjusted_impact": electron_impact * scale if electron_impact is not None else None,
                "electron_adjusted_residual": adjusted_residual(electron_term, electron_impact, scale),
                "electron_required_impact_scale": electron_required,
                "hole_impact": hole_impact,
                "hole_term_sum": hole_term,
                "hole_adjusted_impact": hole_impact * scale if hole_impact is not None else None,
                "hole_adjusted_residual": adjusted_residual(hole_term, hole_impact, scale),
                "hole_required_impact_scale": hole_required,
            })
    return result


def summarize_impact_scale_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    states: dict[str, Any] = {}
    for state in sorted({str(row["state"]) for row in rows}):
        state_rows_subset = [row for row in rows if str(row["state"]) == state]
        state_item: dict[str, Any] = {"support_classes": {}}
        for support_class in sorted({str(row["support_class"]) for row in state_rows_subset}):
            class_rows = [row for row in state_rows_subset if str(row["support_class"]) == support_class]
            class_item: dict[str, Any] = {
                "count": len({int(row["node_id"]) for row in class_rows}),
                "electron_required_impact_scale_median": median(
                    finite_values(class_rows, "electron_required_impact_scale")),
                "hole_required_impact_scale_median": median(
                    finite_values(class_rows, "hole_required_impact_scale")),
                "scales": [],
            }
            for scale in sorted({str(row["impact_scale"]) for row in class_rows}, key=float):
                scale_rows_subset = [row for row in class_rows if str(row["impact_scale"]) == scale]
                electron_values = finite_values(scale_rows_subset, "electron_adjusted_residual")
                hole_values = finite_values(scale_rows_subset, "hole_adjusted_residual")
                class_item["scales"].append({
                    "impact_scale": scale,
                    "electron_adjusted_residual_median": median(electron_values),
                    "electron_adjusted_residual_abs_median": median([abs(v) for v in electron_values]),
                    "hole_adjusted_residual_median": median(hole_values),
                    "hole_adjusted_residual_abs_median": median([abs(v) for v in hole_values]),
                })
            state_item["support_classes"][support_class] = class_item
        states[state] = state_item
    return {"states": states}

def compare_states(left: list[dict[str, Any]], right: list[dict[str, Any]]) -> dict[str, Any]:
    left_by_node = {int(row["node_id"]): row for row in left}
    rows: list[dict[str, Any]] = []
    delta_metrics = [
        "electron_residual", "hole_residual", "electron_raw_delta_phin_V", "hole_raw_delta_phip_V",
        "electron_flux", "electron_recombination", "electron_impact", "electron_term_sum",
        "hole_flux", "hole_recombination", "hole_impact", "hole_term_sum", "impact_combined_source",
    ]
    for row in right:
        node_id = int(row["node_id"])
        if node_id not in left_by_node:
            continue
        item = {"node_id": node_id, "support_class": row["support_class"]}
        for metric in delta_metrics:
            item[f"{metric}_delta"] = (
                finite_float(row.get(metric)) - finite_float(left_by_node[node_id].get(metric))
                if finite_float(row.get(metric)) is not None and finite_float(left_by_node[node_id].get(metric)) is not None
                else None
            )
        rows.append(item)
    classes = sorted({row["support_class"] for row in rows})
    result: dict[str, Any] = {"support_classes": {}}
    for support_class in classes:
        subset = [row for row in rows if row["support_class"] == support_class]
        item: dict[str, Any] = {"count": len(subset)}
        for metric in delta_metrics:
            values = finite_values(subset, f"{metric}_delta")
            item[f"{metric}_delta_median"] = median(values)
            item[f"{metric}_delta_mean"] = mean(values)
        result["support_classes"][support_class] = item
    return result


def state_csv_to_fields(state_csv: Path, fields_dir: Path) -> None:
    rows = read_csv_rows(state_csv)
    values = {field: {} for field in STATE_FIELDS}
    for row in rows:
        node_id = int(row["node_id"])
        values["psi"][node_id] = float(row["psi"])
        values["phin"][node_id] = float(row["phin"])
        values["phip"][node_id] = float(row["phip"])
    write_fields(values, fields_dir)


def step_csv_trial_to_fields(step_csv: Path, fields_dir: Path) -> None:
    values = {field: {} for field in STATE_FIELDS}
    for row in read_csv_rows(step_csv):
        node_id = int(row["node_id"])
        values["psi"][node_id] = float(row["trial_psi"])
        values["phin"][node_id] = float(row["trial_phin"])
        values["phip"][node_id] = float(row["trial_phip"])
    write_fields(values, fields_dir)


def vtk_to_fields(vtk: Path, fields_dir: Path) -> None:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    values = {
        "psi": {idx: value for idx, value in enumerate(scalars["Potential"])},
        "phin": {idx: value for idx, value in enumerate(scalars["ElectronQuasiFermi"])},
        "phip": {idx: value for idx, value in enumerate(scalars["HoleQuasiFermi"])},
    }
    write_fields(values, fields_dir)


def write_fields(values: dict[str, dict[int, float]], fields_dir: Path) -> None:
    fields_dir.mkdir(parents=True, exist_ok=True)
    for source_field, output_name in STATE_FIELDS.items():
        path = fields_dir / f"{output_name}_region0.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["node_id", "component0"])
            for node_id in range(max(values[source_field]) + 1):
                writer.writerow([node_id, f"{values[source_field][node_id]:.17g}"])


def resolve_config_path(base_config: Path, raw: str) -> str:
    path = Path(raw)
    if path.is_absolute():
        return str(path)
    return str((base_config.parent / path).resolve())


def make_probe_config(base_config: Path, output_csv: Path, fields_dir: Path, simulation_type: str,
                      bias: float, bias_contact: str, ground_contact: str) -> dict[str, Any]:
    cfg = json.loads(base_config.read_text(encoding="utf-8-sig"))
    cfg["simulation_type"] = simulation_type
    cfg["output_csv"] = str(output_csv.resolve())
    cfg["state_fields_dir"] = str(fields_dir.resolve())
    cfg.pop("sweep", None)
    for key in ["mesh_file", "node_doping_file", "materials_file"]:
        if key in cfg:
            cfg[key] = resolve_config_path(base_config, cfg[key])
    for contact in cfg.get("contacts", []):
        name = contact.get("name")
        if name == bias_contact:
            contact["bias"] = bias
        elif name == ground_contact:
            contact["bias"] = 0.0
    return cfg


def run_runner(runner: Path, config: Path) -> dict[str, Any]:
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve())],
        cwd=config.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"runner failed for {config}: {result.stderr.strip() or result.stdout.strip()}")
    return json.loads(result.stdout)


def ensure_probe_csvs(args: argparse.Namespace) -> tuple[dict[str, Path], dict[str, Path], dict[str, Any]]:
    row_csvs = parse_named_paths(args.carrier_row_csv)
    term_csvs = parse_named_paths(args.carrier_term_csv)
    if row_csvs and term_csvs:
        return row_csvs, term_csvs, {}
    if args.base_config is None or args.runner is None:
        raise RuntimeError("--base-config and --runner are required when probe CSVs are not supplied")
    state_csvs = parse_named_paths(args.state_csv)
    state_vtks = parse_named_paths(args.state_vtk)
    state_inputs: dict[str, tuple[str, Path]] = {}
    for name, path in state_csvs.items():
        state_inputs[name] = ("csv", path)
    for name, path in state_vtks.items():
        if name in state_inputs:
            raise RuntimeError(f"duplicate state input {name}")
        state_inputs[name] = ("vtk", path)
    if args.include_trial_from_step:
        if args.step_csv is None:
            raise RuntimeError("--step-csv is required with --include-trial-from-step")
        state_inputs["trial"] = ("step_trial", args.step_csv)
    if not state_inputs:
        raise RuntimeError("no state inputs supplied")
    statuses: dict[str, Any] = {}
    configs_dir = args.out_dir / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)
    for state, (kind, path) in state_inputs.items():
        fields_dir = args.out_dir / "state_fields" / state
        if kind == "csv":
            state_csv_to_fields(path, fields_dir)
        elif kind == "vtk":
            vtk_to_fields(path, fields_dir)
        elif kind == "step_trial":
            step_csv_trial_to_fields(path, fields_dir)
        row_csv = args.out_dir / "probes" / f"{state}_carrier_rows.csv"
        term_csv = args.out_dir / "probes" / f"{state}_carrier_terms.csv"
        row_csv.parent.mkdir(parents=True, exist_ok=True)
        row_config = configs_dir / f"{state}_carrier_row_probe.json"
        term_config = configs_dir / f"{state}_carrier_term_probe.json"
        row_config.write_text(json.dumps(make_probe_config(
            args.base_config, row_csv, fields_dir, "newton_carrier_row_probe",
            args.bias, args.bias_contact, args.ground_contact), indent=2) + "\n", encoding="utf-8")
        term_config.write_text(json.dumps(make_probe_config(
            args.base_config, term_csv, fields_dir, "newton_carrier_term_probe",
            args.bias, args.bias_contact, args.ground_contact), indent=2) + "\n", encoding="utf-8")
        row_csvs[state] = row_csv
        term_csvs[state] = term_csv
        statuses[state] = {"row_config": str(row_config), "term_config": str(term_config)}
        if not args.prepare_only:
            statuses[state]["row_status"] = run_runner(args.runner, row_config)
            statuses[state]["term_status"] = run_runner(args.runner, term_config)
    if statuses:
        (args.out_dir / "predictor_carrier_row_probe_status.json").write_text(
            json.dumps(clean_json(statuses), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return row_csvs, term_csvs, statuses


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    row_csvs, term_csvs, statuses = ensure_probe_csvs(args)
    if args.prepare_only:
        print(json.dumps(clean_json(statuses), sort_keys=True))
        return 0
    meta = support_metadata(args.support_csv)
    states = [state for state in row_csvs if state in term_csvs]
    state_rows: dict[str, list[dict[str, Any]]] = {}
    all_rows: list[dict[str, Any]] = []
    for state in states:
        rows = state_node_rows(state, row_by_node(row_csvs[state]), row_by_node(term_csvs[state]), meta)
        state_rows[state] = rows
        all_rows.extend(rows)
    write_csv(args.out_dir / "predictor_carrier_row_audit_nodes.csv", all_rows, NODE_FIELDS)
    scale_rows = impact_scale_rows(all_rows, args.impact_scale) if args.impact_scale else []
    if scale_rows:
        write_csv(
            args.out_dir / "predictor_carrier_impact_scale_nodes.csv",
            scale_rows,
            IMPACT_SCALE_NODE_FIELDS,
        )
    comparisons: dict[str, Any] = {}
    for left, right in [("baseline", "predictor"), ("predictor", "trial"), ("baseline", "trial")]:
        if left in state_rows and right in state_rows:
            comparisons[f"{right}_minus_{left}"] = compare_states(state_rows[left], state_rows[right])
    summary = {
        "states": {state: summarize_rows(rows) for state, rows in state_rows.items()},
        "comparisons": comparisons,
        "probe_status": statuses,
    }
    if scale_rows:
        summary["impact_scale_sensitivity"] = summarize_impact_scale_rows(scale_rows)
    (args.out_dir / "predictor_carrier_row_audit_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
