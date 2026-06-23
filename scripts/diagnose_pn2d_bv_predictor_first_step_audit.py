#!/usr/bin/env python3
"""Audit the first Newton update from a PN2D BV predictor state.

The audit can either analyze existing newton_step_probe/newton_block_step_probe
CSV files or prepare and run those probes from a base runner deck plus a
DCSweep-compatible state CSV.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Any

NODE_FIELDS = [
    "mode",
    "node_id",
    "support_class",
    "initial_psi_minus_phin_shift_V",
    "initial_phip_minus_psi_shift_V",
    "delta_psi_minus_phin_V",
    "delta_phip_minus_psi_V",
    "electron_rollback_fraction",
    "hole_rollback_fraction",
    "delta_psi_V",
    "delta_phin_V",
    "delta_phip_V",
    "phin_residual",
    "phip_residual",
    "trial_phin_residual",
    "trial_phip_residual",
]

SCALAR_FIELDS = {
    "psi": "ElectrostaticPotential",
    "phin": "eQuasiFermiPotential",
    "phip": "hQuasiFermiPotential",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--step-csv", type=Path)
    parser.add_argument("--block-step-csv", type=Path)
    parser.add_argument("--base-config", type=Path)
    parser.add_argument("--state-csv", type=Path)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--bias", type=float, default=-20.0)
    parser.add_argument("--bias-contact", default="Anode")
    parser.add_argument("--ground-contact", default="Cathode")
    parser.add_argument("--block-modes", default="poisson_only,carrier_only")
    parser.add_argument("--prepare-only", action="store_true")
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


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = finite_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def clean_json(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


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
        raise RuntimeError(
            "support csv must contain predictor metadata columns phin_baseline_V, "
            "phin_predictor_V, phip_baseline_V, and phip_predictor_V")
    return meta


def normalized_step_rows(path: Path, mode: str) -> list[dict[str, str]]:
    rows = []
    for row in read_csv_rows(path):
        item = dict(row)
        item["mode"] = mode
        rows.append(item)
    return rows


def audit_rows_for_mode(rows: list[dict[str, str]], meta: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in rows:
        node_id = int(row["node_id"])
        if node_id not in meta:
            continue
        node_meta = meta[node_id]
        initial_e = finite_float(node_meta["initial_psi_minus_phin_shift_V"])
        initial_h = finite_float(node_meta["initial_phip_minus_psi_shift_V"])
        delta_e = finite_float(row.get("delta_psi_minus_phin_V"))
        delta_h = finite_float(row.get("delta_phip_minus_psi_V"))
        out.append({
            "mode": row.get("mode", "full_newton"),
            "node_id": node_id,
            "support_class": node_meta["support_class"],
            "initial_psi_minus_phin_shift_V": initial_e,
            "initial_phip_minus_psi_shift_V": initial_h,
            "delta_psi_minus_phin_V": delta_e,
            "delta_phip_minus_psi_V": delta_h,
            "electron_rollback_fraction": safe_ratio(-(delta_e or 0.0), initial_e),
            "hole_rollback_fraction": safe_ratio(-(delta_h or 0.0), initial_h),
            "delta_psi_V": finite_float(row.get("delta_psi_V")),
            "delta_phin_V": finite_float(row.get("delta_phin_V")),
            "delta_phip_V": finite_float(row.get("delta_phip_V")),
            "phin_residual": finite_float(row.get("phin_residual")),
            "phip_residual": finite_float(row.get("phip_residual")),
            "trial_phin_residual": finite_float(row.get("trial_phin_residual")),
            "trial_phip_residual": finite_float(row.get("trial_phip_residual")),
        })
    return out


def summarize_mode(rows: list[dict[str, Any]]) -> dict[str, Any]:
    classes = sorted({row["support_class"] for row in rows})
    result: dict[str, Any] = {"count": len(rows), "support_classes": {}}
    for support_class in classes:
        subset = [row for row in rows if row["support_class"] == support_class]
        result["support_classes"][support_class] = {
            "count": len(subset),
            "electron_rollback_fraction_median": median(finite_values(subset, "electron_rollback_fraction")),
            "electron_rollback_fraction_mean": mean(finite_values(subset, "electron_rollback_fraction")),
            "hole_rollback_fraction_median": median(finite_values(subset, "hole_rollback_fraction")),
            "hole_rollback_fraction_mean": mean(finite_values(subset, "hole_rollback_fraction")),
            "delta_psi_minus_phin_median_V": median(finite_values(subset, "delta_psi_minus_phin_V")),
            "delta_phip_minus_psi_median_V": median(finite_values(subset, "delta_phip_minus_psi_V")),
            "delta_psi_median_V": median(finite_values(subset, "delta_psi_V")),
            "delta_phin_median_V": median(finite_values(subset, "delta_phin_V")),
            "delta_phip_median_V": median(finite_values(subset, "delta_phip_V")),
            "phin_residual_abs_median": median([abs(v) for v in finite_values(subset, "phin_residual")]),
            "phip_residual_abs_median": median([abs(v) for v in finite_values(subset, "phip_residual")]),
            "trial_phin_residual_abs_median": median([abs(v) for v in finite_values(subset, "trial_phin_residual")]),
            "trial_phip_residual_abs_median": median([abs(v) for v in finite_values(subset, "trial_phip_residual")]),
        }
    return result


def state_csv_to_fields(state_csv: Path, fields_dir: Path) -> None:
    rows = read_csv_rows(state_csv)
    values = {field: {} for field in SCALAR_FIELDS}
    for row in rows:
        node_id = int(row["node_id"])
        values["psi"][node_id] = float(row["psi"])
        values["phin"][node_id] = float(row["phin"])
        values["phip"][node_id] = float(row["phip"])
    fields_dir.mkdir(parents=True, exist_ok=True)
    for source_field, output_name in SCALAR_FIELDS.items():
        output = fields_dir / f"{output_name}_region0.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["node_id", "component0"])
            for node_id in range(max(values[source_field]) + 1):
                writer.writerow([node_id, f"{values[source_field][node_id]:.17g}"])


def resolve_config_path(base_config: Path, raw: str) -> str:
    path = Path(raw)
    if path.is_absolute():
        return str(path)
    return str((base_config.parent / path).resolve())


def make_probe_config(
    base_config: Path,
    output_csv: Path,
    fields_dir: Path,
    simulation_type: str,
    bias: float,
    bias_contact: str,
    ground_contact: str,
    block_modes: list[str] | None = None,
) -> dict[str, Any]:
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
    if block_modes is not None:
        cfg["block_modes"] = block_modes
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
        raise RuntimeError(
            f"runner failed for {config}: {result.stderr.strip() or result.stdout.strip()}")
    return json.loads(result.stdout)


def ensure_probe_outputs(args: argparse.Namespace) -> tuple[Path, Path, dict[str, Any]]:
    if args.step_csv and args.block_step_csv:
        return args.step_csv, args.block_step_csv, {}
    missing = [name for name, value in [
        ("--base-config", args.base_config),
        ("--state-csv", args.state_csv),
    ] if value is None]
    if missing:
        raise RuntimeError("missing required probe inputs: " + ", ".join(missing))
    if args.runner is None and not args.prepare_only:
        raise RuntimeError("--runner is required unless --prepare-only or existing CSVs are supplied")

    fields_dir = args.out_dir / "state_fields"
    configs_dir = args.out_dir / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)
    state_csv_to_fields(args.state_csv, fields_dir)
    block_modes = [item.strip() for item in args.block_modes.split(",") if item.strip()]
    step_csv = args.out_dir / "newton_step_probe.csv"
    block_csv = args.out_dir / "newton_block_step_probe.csv"
    step_config = configs_dir / "newton_step_probe.json"
    block_config = configs_dir / "newton_block_step_probe.json"
    step_config.write_text(json.dumps(make_probe_config(
        args.base_config, step_csv, fields_dir, "newton_step_probe",
        args.bias, args.bias_contact, args.ground_contact), indent=2) + "\n", encoding="utf-8")
    block_config.write_text(json.dumps(make_probe_config(
        args.base_config, block_csv, fields_dir, "newton_block_step_probe",
        args.bias, args.bias_contact, args.ground_contact, block_modes), indent=2) + "\n", encoding="utf-8")
    statuses: dict[str, Any] = {"prepared_configs": [str(step_config), str(block_config)]}
    if args.prepare_only:
        return step_csv, block_csv, statuses
    statuses["newton_step_probe"] = run_runner(args.runner, step_config)
    statuses["newton_block_step_probe"] = run_runner(args.runner, block_config)
    (args.out_dir / "predictor_first_step_probe_status.json").write_text(
        json.dumps(clean_json(statuses), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return step_csv, block_csv, statuses


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    step_csv, block_csv, statuses = ensure_probe_outputs(args)
    if args.prepare_only:
        print(json.dumps(clean_json(statuses), sort_keys=True))
        return 0
    meta = support_metadata(args.support_csv)
    rows_by_mode: dict[str, list[dict[str, Any]]] = {}
    rows_by_mode["full_newton"] = audit_rows_for_mode(normalized_step_rows(step_csv, "full_newton"), meta)
    for row in read_csv_rows(block_csv):
        mode = row.get("mode", "")
        rows_by_mode.setdefault(mode, []).extend(audit_rows_for_mode([row], meta))
    all_rows = [row for rows in rows_by_mode.values() for row in rows]
    write_csv(args.out_dir / "predictor_first_step_audit_nodes.csv", all_rows, NODE_FIELDS)
    summary = {
        "step_csv": str(step_csv),
        "block_step_csv": str(block_csv),
        "support_node_count": len(meta),
        "modes": {mode: summarize_mode(rows) for mode, rows in sorted(rows_by_mode.items())},
        "probe_status": statuses,
    }
    (args.out_dir / "predictor_first_step_audit_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
