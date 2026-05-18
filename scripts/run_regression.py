#!/usr/bin/env python3
"""Run the Vela engineering examples as a lightweight regression suite."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

EXAMPLES = [
    {
        "name": "pn_diode_iv",
        "source": "pn_diode",
        "config": Path("examples/pn_diode/simulation_iv.json"),
        "expected": [Path("outputs/pn_iv.csv"), Path("outputs/pn_sweep_0000_0V.vtk")],
        "checks": ["csv_converged", "finite_outputs", "iv_trend", "dc_sweep_regression"],
    },
    {
        "name": "pn_diode_cv",
        "source": "pn_diode",
        "config": Path("examples/pn_diode/simulation_cv.json"),
        "expected": [Path("outputs/pn_cv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "pn_diode_bv",
        "source": "pn_diode",
        "config": Path("examples/pn_diode/simulation_bv.json"),
        "expected": [Path("outputs/pn_bv.csv"), Path("outputs/pn_bv_sweep_0000_0V.vtk")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "moscap",
        "source": "moscap",
        "config": Path("examples/moscap/simulation.json"),
        "expected": [Path("outputs/moscap_poisson.vtk")],
        "checks": ["finite_outputs", "moscap_interface"],
    },
    {
        "name": "nmos2d",
        "source": "nmos2d",
        "config": Path("examples/nmos2d/simulation.json"),
        "expected": [Path("outputs/nmos_poisson.vtk")],
        "checks": ["finite_outputs"],
    },
    {
        "name": "nmos2d_dd_iv",
        "source": "nmos2d_dd",
        "config": Path("examples/nmos2d_dd/simulation_iv.json"),
        "expected": [
            Path("outputs/nmos2d_dd_iv.csv"),
            Path("outputs/nmos2d_dd_sweep_0000_0V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression", "mos_trends"],
    },
    {
        "name": "nmos2d_dd_idvg",
        "source": "nmos2d_dd",
        "config": Path("examples/nmos2d_dd/simulation_idvg.json"),
        "expected": [Path("outputs/nmos2d_dd_idvg.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "nmos2d_dd_cv",
        "source": "nmos2d_dd",
        "config": Path("examples/nmos2d_dd/simulation_cv.json"),
        "expected": [Path("outputs/nmos2d_dd_cv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "nmos2d_mos_dd_iv",
        "source": "nmos2d_mos_dd",
        "config": Path("examples/nmos2d_mos_dd/simulation_iv.json"),
        "expected": [
            Path("outputs/nmos2d_mos_dd_iv.csv"),
            Path("outputs/nmos2d_mos_dd_iv_sweep_0000_0V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression", "mos_trends"],
    },
    {
        "name": "nmos2d_mos_dd_idvg",
        "source": "nmos2d_mos_dd",
        "config": Path("examples/nmos2d_mos_dd/simulation_idvg.json"),
        "expected": [Path("outputs/nmos2d_mos_dd_idvg.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "nmos2d_mos_dd_cv",
        "source": "nmos2d_mos_dd",
        "config": Path("examples/nmos2d_mos_dd/simulation_cv.json"),
        "expected": [Path("outputs/nmos2d_mos_dd_cv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "nmos2d_mos_dd_bv",
        "source": "nmos2d_mos_dd",
        "config": Path("examples/nmos2d_mos_dd/simulation_bv.json"),
        "expected": [Path("outputs/nmos2d_mos_dd_bv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "pmos2d_dd_iv",
        "source": "pmos2d_dd",
        "config": Path("examples/pmos2d_dd/simulation_iv.json"),
        "expected": [
            Path("outputs/pmos2d_dd_iv.csv"),
            Path("outputs/pmos2d_dd_sweep_0000_0V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression", "mos_trends"],
    },
    {
        "name": "pmos2d_dd_idvg",
        "source": "pmos2d_dd",
        "config": Path("examples/pmos2d_dd/simulation_idvg.json"),
        "expected": [Path("outputs/pmos2d_dd_idvg.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "pmos2d_dd_cv",
        "source": "pmos2d_dd",
        "config": Path("examples/pmos2d_dd/simulation_cv.json"),
        "expected": [Path("outputs/pmos2d_dd_cv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "pmos2d_mos_dd_iv",
        "source": "pmos2d_mos_dd",
        "config": Path("examples/pmos2d_mos_dd/simulation_iv.json"),
        "expected": [
            Path("outputs/pmos2d_mos_dd_iv.csv"),
            Path("outputs/pmos2d_mos_dd_iv_sweep_0000_0V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression", "mos_trends"],
    },
    {
        "name": "pmos2d_mos_dd_idvg",
        "source": "pmos2d_mos_dd",
        "config": Path("examples/pmos2d_mos_dd/simulation_idvg.json"),
        "expected": [Path("outputs/pmos2d_mos_dd_idvg.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "pmos2d_mos_dd_cv",
        "source": "pmos2d_mos_dd",
        "config": Path("examples/pmos2d_mos_dd/simulation_cv.json"),
        "expected": [Path("outputs/pmos2d_mos_dd_cv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "pmos2d_mos_dd_bv",
        "source": "pmos2d_mos_dd",
        "config": Path("examples/pmos2d_mos_dd/simulation_bv.json"),
        "expected": [Path("outputs/pmos2d_mos_dd_bv.csv")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "ldmos2d_poisson",
        "source": "ldmos2d",
        "config": Path("examples/ldmos2d/simulation_iv.json"),
        "expected": [Path("outputs/ldmos2d_reverse_poisson.vtk")],
        "checks": ["finite_outputs"],
    },
    {
        "name": "igbt2d_poisson",
        "source": "igbt2d",
        "config": Path("examples/igbt2d/simulation_poisson.json"),
        "expected": [Path("outputs/igbt2d_poisson.vtk")],
        "checks": ["finite_outputs"],
    },
    {
        "name": "igbt2d_iv",
        "source": "igbt2d",
        "config": Path("examples/igbt2d/simulation_iv.json"),
        "expected": [Path("outputs/igbt2d_iv.csv"), Path("outputs/igbt2d_iv_sweep_0000_0V.vtk")],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression"],
    },
    {
        "name": "schottky_diode_2d_iv",
        "source": "schottky_diode_2d",
        "config": Path("examples/schottky_diode_2d/simulation_iv.json"),
        "expected": [
            Path("outputs/schottky_iv.csv"),
            Path("outputs/schottky_sweep_0000_-0.1V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression",
                   "schottky_iv_trend"],
    },
]



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--runner", type=Path, default=None, help="Path to vela_example_runner")
    parser.add_argument("--workdir", type=Path, default=None, help="Directory for copied examples and outputs")
    return parser.parse_args()


def default_runner(repo: Path) -> Path:
    candidates = [repo / "build" / "vela_example_runner", repo / "build" / "src" / "tools" / "vela_example_runner"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def copy_example(repo: Path, workdir: Path, name: str, source: str | None = None) -> Path:
    src = repo / "examples" / (source or name)
    dst = workdir / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    (dst / "outputs").mkdir(parents=True, exist_ok=True)
    return dst


def has_nan_or_inf(path: Path) -> bool:
    text = path.read_text(errors="ignore")
    for token in text.replace(",", " ").split():
        try:
            value = float(token)
        except ValueError:
            lowered = token.lower()
            if lowered in {"nan", "+nan", "-nan", "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"}:
                return True
            continue
        if not math.isfinite(value):
            return True
    return False


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def read_vtk_scalar(path: Path, field: str) -> list[float]:
    lines = path.read_text().splitlines()
    marker = f"SCALARS {field} double 1"
    for idx, line in enumerate(lines):
        if line.strip() == marker:
            start = idx + 2
            values = []
            for raw in lines[start:]:
                stripped = raw.strip()
                if not stripped or stripped.startswith("SCALARS") or stripped.startswith("POINT_DATA"):
                    break
                try:
                    values.append(float(stripped))
                except ValueError:
                    break
            return values
    raise AssertionError(f"VTK field '{field}' not found in {path}")


def node_index(mesh: dict[str, Any], x: float, y: float) -> int:
    best = min(mesh["nodes"], key=lambda n: abs(n["x"] - x) + abs(n["y"] - y))
    if abs(best["x"] - x) > 1e-15 or abs(best["y"] - y) > 1e-15:
        raise AssertionError(f"No mesh node near ({x}, {y})")
    return int(best["id"])


def output_csv_path(example_dir: Path, cfg: dict[str, Any]) -> Path:
    if "output_csv" in cfg:
        return example_dir / cfg["output_csv"]
    return example_dir / cfg.get("sweep", {}).get("csv_file", "dc_sweep.csv")


def dc_sweep_regression_options(cfg: dict[str, Any]) -> tuple[dict[str, Any], str, bool, bool]:
    regression = cfg.get("regression", {})
    reg = regression.get("dc_sweep", {})
    mode = normalize_curve_mode(str(cfg.get("sweep", {}).get("mode", "iv")))
    declared = bool(regression.get("declared_converged", True))
    allow_final_bv = bool(reg.get("allow_nonconverged_final_bv_point", False)) and mode == "bv_reverse"
    return reg, mode, declared, allow_final_bv


def nonconverged_row_allowed(
    *,
    declared_converged: bool,
    allow_final_bv_nonconverged: bool,
    row_index: int,
    row_count: int,
) -> bool:
    return (not declared_converged) or (allow_final_bv_nonconverged and row_index == row_count)


def disallowed_nonconverged_rows(
    rows: list[dict[str, str]],
    *,
    declared_converged: bool,
    allow_final_bv_nonconverged: bool,
) -> list[tuple[int, dict[str, str]]]:
    bad: list[tuple[int, dict[str, str]]] = []
    for row_index, row in enumerate(rows, start=1):
        if row.get("converged") == "1":
            continue
        if nonconverged_row_allowed(
            declared_converged=declared_converged,
            allow_final_bv_nonconverged=allow_final_bv_nonconverged,
            row_index=row_index,
            row_count=len(rows),
        ):
            continue
        bad.append((row_index, row))
    return bad


def nonzero_runner_exit_allowed(cfg: dict[str, Any]) -> bool:
    _reg, _mode, declared, allow_final_bv = dc_sweep_regression_options(cfg)
    return (not declared) or allow_final_bv


def check_csv_converged(example_dir: Path) -> None:
    cfg = json.loads((example_dir / "simulation.json").read_text())
    rows = read_csv(output_csv_path(example_dir, cfg))
    if not rows:
        raise AssertionError(f"{example_dir.name}: DC sweep CSV contains no sweep rows")
    _reg, _mode, declared, allow_final_bv = dc_sweep_regression_options(cfg)

    bad = disallowed_nonconverged_rows(
        rows,
        declared_converged=declared,
        allow_final_bv_nonconverged=allow_final_bv,
    )
    if bad:
        diagnostics = [
            {
                "example": example_dir.name,
                "row_index": row_index,
                "field": "converged",
                "actual": row.get("converged"),
                "bias_V": row.get("bias_V"),
                "failure_reason": row.get("failure_reason", ""),
                "validation_diagnostics": row.get("validation_diagnostics", ""),
            }
            for row_index, row in bad
        ]
        raise AssertionError(f"Non-converged declared sweep rows: {diagnostics}")


def expected_sweep_voltages(sweep: dict[str, Any]) -> list[float]:
    start = float(sweep["start"])
    stop = float(sweep["stop"])
    step = float(sweep["step"])
    if step == 0.0:
        raise AssertionError("PN sweep step must be non-zero")
    direction = 1.0 if step > 0.0 else -1.0
    if (stop - start) * step < 0.0:
        raise AssertionError("PN sweep step does not move start toward stop")

    values = [start]
    voltage = start + step
    while direction * (voltage - stop) <= 1.0e-12:
        values.append(voltage)
        voltage += step
    if abs(values[-1] - stop) > 1.0e-12:
        values.append(stop)
    return values



def normalize_curve_mode(mode: str) -> str:
    normalized = mode.lower().replace("-", "_")
    if normalized in {"", "iv"}:
        return "iv"
    if normalized in {"cv", "cv_quasistatic"}:
        return "cv_quasistatic"
    if normalized in {"bv", "bv_reverse", "reverse_breakdown"}:
        return "bv_reverse"
    return normalized


def cv_charge_columns(sweep: dict[str, Any], row: dict[str, str]) -> tuple[str, str]:
    charge_cfg = sweep.get("terminal_charge", {})
    per_meter = bool(charge_cfg.get("per_meter", sweep.get("charge_per_meter", True)))
    preferred = ("charge_C_per_m", "capacitance_F_per_m") if per_meter else ("charge_C", "capacitance_F")
    alternate = ("charge_C", "capacitance_F") if per_meter else ("charge_C_per_m", "capacitance_F_per_m")
    if all(column in row for column in preferred):
        return preferred
    if all(column in row for column in alternate):
        return alternate
    return preferred


def curve_column(row: dict[str, str], modern: str, legacy: str | None = None) -> str:
    if modern in row:
        return row[modern]
    if legacy is not None and legacy in row:
        return row[legacy]
    raise KeyError(modern)


def parse_step_diagnostics(value: str) -> dict[str, str]:
    diagnostics: dict[str, str] = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Malformed step diagnostic item {item!r}")
        key, raw = item.split("=", maxsplit=1)
        diagnostics[key] = raw
    return diagnostics


def check_dc_sweep_regression(example_dir: Path) -> dict[str, Any]:
    cfg = json.loads((example_dir / "simulation.json").read_text())
    example = example_dir.name
    rows = read_csv(output_csv_path(example_dir, cfg))
    if not rows:
        raise AssertionError(f"{example}: DC sweep CSV contains no sweep rows")

    sweep_cfg = cfg.get("sweep", {})
    reg, mode, declared_converged, allow_final_bv_nonconverged = dc_sweep_regression_options(cfg)

    expected_rows = reg.get("expected_rows")
    if expected_rows is None and "sweep" in cfg:
        expected_rows = len(expected_sweep_voltages(cfg["sweep"]))
    if expected_rows is not None and len(rows) != int(expected_rows):
        raise AssertionError(
            f"{example}: field=rows actual={len(rows)} expected={expected_rows}")

    max_abs_attempted = float(reg.get("max_abs_attempted_step", math.inf))
    max_abs_accepted = float(reg.get("max_abs_accepted_step", math.inf))
    max_retry = int(reg.get("max_retry_count", 0))
    require_monotone_abs_current = bool(reg.get("require_monotone_abs_current", False))
    current_abs_tolerance = float(reg.get("current_monotone_abs_tolerance", 1.0e-30))
    current_rel_tolerance = float(reg.get("current_monotone_rel_tolerance", 1.0e-12))
    require_monotone_max_field = bool(reg.get("require_monotone_max_field", False))
    max_field_abs_tolerance = float(reg.get("max_field_monotone_abs_tolerance", 1.0e-9))
    max_field_rel_tolerance = float(reg.get("max_field_monotone_rel_tolerance", 1.0e-12))
    min_converged_rows = reg.get("min_converged_rows")
    min_max_field = reg.get("min_max_electric_field_V_per_m")
    max_max_field = reg.get("max_max_electric_field_V_per_m")
    allow_zero_capacitance = bool(reg.get("allow_zero_capacitance", False))
    expected_zero_capacitance_rows = reg.get("expected_zero_capacitance_rows")
    min_nonzero_capacitance_rows = int(reg.get(
        "min_nonzero_capacitance_rows",
        0 if expected_zero_capacitance_rows is not None else 1))
    modern_required = ("mode", "bias_contact", "bias_V", "current_contact",
                       "current_electron", "current_hole", "current_total",
                       "converged", "iterations", "step_diagnostics")
    missing_modern = [column for column in modern_required if column not in rows[0]]
    if missing_modern:
        raise AssertionError(
            f"{example}: DC sweep CSV is missing required curve schema column(s): "
            + ", ".join(missing_modern))

    if mode == "cv_quasistatic":
        for column in cv_charge_columns(sweep_cfg, rows[0]):
            if column not in rows[0]:
                raise AssertionError(f"{example}: CV sweep row 1 is missing column '{column}'")
    if mode == "bv_reverse":
        for column in ("max_electric_field_V_per_m", "current_jump_ratio",
                       "breakdown_detected", "breakdown_voltage", "criterion",
                       "last_stable_bias", "failed_bias", "failure_reason"):
            if column not in rows[0]:
                raise AssertionError(f"{example}: BV sweep row 1 is missing column '{column}'")

    max_attempted_seen = 0.0
    max_accepted_seen = 0.0
    max_retry_seen = 0
    converged_rows = 0
    final_current_total = 0.0
    current_values: list[tuple[int, float]] = []
    max_fields: list[tuple[int, float]] = []
    max_electric_field_seen = 0.0
    breakdown_detected = False
    nonzero_capacitance_rows = 0
    zero_capacitance_rows = 0
    for row_index, row in enumerate(rows, start=1):
        label = f"{example}: DC sweep"
        converged = row.get("converged") == "1"
        if converged:
            converged_rows += 1
        current_total = parse_finite_float(row, "current_total", label, row_index)
        final_current_total = current_total
        if converged:
            current_values.append((row_index, abs(current_total)))
        try:
            diagnostics = parse_step_diagnostics(row["step_diagnostics"])
            attempted = abs(float(diagnostics["attempted_step"]))
            accepted = abs(float(diagnostics["accepted_step"]))
            retry = int(diagnostics["retry_count"])
        except (KeyError, ValueError) as exc:
            raise AssertionError(
                f"{example}: row {row_index} field=step_diagnostics actual={row.get('step_diagnostics')!r} "
                f"is unparseable") from exc
        max_attempted_seen = max(max_attempted_seen, attempted)
        max_accepted_seen = max(max_accepted_seen, accepted)
        max_retry_seen = max(max_retry_seen, retry)
        if mode == "cv_quasistatic":
            charge_column, capacitance_column = cv_charge_columns(sweep_cfg, row)
            _ = parse_finite_float(row, charge_column, f"{example}: CV sweep", row_index)
            capacitance = parse_finite_float(row, capacitance_column, f"{example}: CV sweep", row_index)
            if row_index > 1:
                if abs(capacitance) > 0.0:
                    nonzero_capacitance_rows += 1
                else:
                    zero_capacitance_rows += 1
                    if not allow_zero_capacitance:
                        raise AssertionError(
                            f"{example}: row {row_index} field={capacitance_column} actual={capacitance} "
                            "is zero differential capacitance")
        if mode == "bv_reverse":
            max_field = parse_finite_float(row, "max_electric_field_V_per_m", f"{example}: BV sweep", row_index)
            max_electric_field_seen = max(max_electric_field_seen, max_field)
            if max_field < 0.0:
                raise AssertionError(
                    f"{example}: row {row_index} field=max_electric_field_V_per_m actual={max_field} is negative")
            if converged:
                if min_max_field is not None and max_field < float(min_max_field):
                    raise AssertionError(
                        f"{example}: row {row_index} field=max_electric_field_V_per_m actual={max_field} "
                        f"is below regression minimum {min_max_field}")
                if max_max_field is not None and max_field > float(max_max_field):
                    raise AssertionError(
                        f"{example}: row {row_index} field=max_electric_field_V_per_m actual={max_field} "
                        f"exceeds regression maximum {max_max_field}")
                max_fields.append((row_index, max_field))
            try:
                jump = float(row["current_jump_ratio"])
            except (KeyError, ValueError) as exc:
                raise AssertionError(
                    f"{example}: row {row_index} field=current_jump_ratio actual={row.get('current_jump_ratio')!r} "
                    "is invalid") from exc
            if math.isnan(jump) or jump < 0.0:
                raise AssertionError(
                    f"{example}: row {row_index} field=current_jump_ratio actual={jump} is negative or NaN")
            row_breakdown = row.get("breakdown_detected") == "1"
            breakdown_detected = breakdown_detected or row_breakdown
            if row_breakdown and not row.get("criterion"):
                raise AssertionError(
                    f"{example}: row {row_index} field=criterion actual={row.get('criterion')!r} "
                    "is empty after breakdown_detected=1")
            if row_breakdown and row.get("criterion") == "last_stable_before_nonconvergence":
                last_stable = parse_finite_float(row, "last_stable_bias", f"{example}: BV sweep", row_index)
                failed_bias = parse_finite_float(row, "failed_bias", f"{example}: BV sweep", row_index)
                if not row.get("failure_reason"):
                    raise AssertionError(
                        f"{example}: row {row_index} field=failure_reason actual={row.get('failure_reason')!r} "
                        "is empty for non-convergence breakdown")
                stop = float(sweep_cfg["stop"])
                direction = 1.0 if stop >= last_stable else -1.0
                if direction * (failed_bias - last_stable) <= 0.0:
                    raise AssertionError(
                        f"{example}: row {row_index} field=failed_bias actual={failed_bias} is not beyond "
                        f"last_stable_bias {last_stable} toward stop {stop}")
        if attempted > max_abs_attempted + 1.0e-12:
            raise AssertionError(
                f"{example}: row {row_index} field=attempted_step actual={attempted} "
                f"exceeds regression limit {max_abs_attempted}")
        if accepted > max_abs_accepted + 1.0e-12:
            raise AssertionError(
                f"{example}: row {row_index} field=accepted_step actual={accepted} "
                f"exceeds regression limit {max_abs_accepted}")
        if retry > max_retry:
            raise AssertionError(
                f"{example}: row {row_index} field=retry_count actual={retry} "
                f"exceeds regression limit {max_retry}")
        if not converged and not nonconverged_row_allowed(
            declared_converged=declared_converged,
            allow_final_bv_nonconverged=allow_final_bv_nonconverged,
            row_index=row_index,
            row_count=len(rows),
        ):
            raise AssertionError(
                f"{example}: row {row_index} field=converged actual={row.get('converged')!r} "
                "is non-converged and not allowed by regression.dc_sweep")

    if min_converged_rows is not None and converged_rows < int(min_converged_rows):
        raise AssertionError(
            f"{example}: field=converged_rows actual={converged_rows} expected_at_least={min_converged_rows}")

    current_trend_checked = False
    if require_monotone_abs_current:
        if not current_values:
            raise AssertionError(
                f"{example}: require_monotone_abs_current requested but no converged current values were read")
        for idx in range(1, len(current_values)):
            previous_row, previous = current_values[idx - 1]
            current_row, current = current_values[idx]
            tolerance = current_abs_tolerance + current_rel_tolerance * max(abs(previous), abs(current))
            if current + tolerance < previous:
                values = [value for _, value in current_values]
                raise AssertionError(
                    f"{example}: row {current_row} field=abs(current_total) actual={current} decreased "
                    f"below row {previous_row} value {previous} beyond tolerance {tolerance}; "
                    f"converged values: {values}")
        current_trend_checked = True

    max_field_trend_checked = False
    if mode == "bv_reverse" and require_monotone_max_field:
        if not max_fields:
            raise AssertionError(
                f"{example}: require_monotone_max_field requested but no converged field values were read")
        for idx in range(1, len(max_fields)):
            previous_row, previous = max_fields[idx - 1]
            current_row, current = max_fields[idx]
            tolerance = max_field_abs_tolerance + max_field_rel_tolerance * max(abs(previous), abs(current), 1.0)
            if current + tolerance < previous:
                max_field_values = [field for _, field in max_fields]
                raise AssertionError(
                    f"{example}: row {current_row} field=max_electric_field_V_per_m actual={current} "
                    f"decreased below row {previous_row} value {previous} beyond tolerance {tolerance}; "
                    f"converged values: {max_field_values}")
        max_field_trend_checked = True

    if mode == "cv_quasistatic":
        if nonzero_capacitance_rows < min_nonzero_capacitance_rows:
            raise AssertionError(
                f"{example}: field=nonzero_capacitance_rows actual={nonzero_capacitance_rows} "
                f"expected_at_least={min_nonzero_capacitance_rows}")
        if (expected_zero_capacitance_rows is not None
                and zero_capacitance_rows != int(expected_zero_capacitance_rows)):
            raise AssertionError(
                f"{example}: field=zero_capacitance_rows actual={zero_capacitance_rows} "
                f"expected={expected_zero_capacitance_rows}")

    result = {
        "rows": len(rows),
        "converged_rows": converged_rows,
        "max_attempted_step": max_attempted_seen,
        "max_accepted_step": max_accepted_seen,
        "max_retry_count_seen": max_retry_seen,
        "final_current_total": final_current_total,
        "current_trend_checked": current_trend_checked,
        # Retain the historic keys for downstream consumers while adding the
        # clearer P0-stability summary field names above.
        "max_abs_attempted_step": max_attempted_seen,
        "max_abs_accepted_step": max_accepted_seen,
        "max_retry_count": max_retry_seen,
    }
    if mode == "bv_reverse":
        result["max_electric_field_seen"] = max_electric_field_seen
        result["breakdown_detected"] = breakdown_detected
        result["max_electric_field_V_per_m"] = [field for _, field in max_fields]
        result["max_field_trend_checked"] = max_field_trend_checked
    if mode == "cv_quasistatic":
        result["nonzero_capacitance_rows"] = nonzero_capacitance_rows
    return result

def check_schottky_iv_trend(example_dir: Path) -> dict[str, Any]:
    cfg = json.loads((example_dir / "simulation.json").read_text())
    rows = read_csv(output_csv_path(example_dir, cfg))
    if not rows:
        raise AssertionError("Schottky IV CSV contains no rows")
    voltages = [
        parse_finite_float(row, "bias_V", "Schottky IV", idx)
        for idx, row in enumerate(rows, start=1)
    ]
    currents = [
        parse_finite_float(row, "current_total", "Schottky IV", idx)
        for idx, row in enumerate(rows, start=1)
    ]
    # The Dirichlet-barrier prototype is not calibrated, but the engineering
    # smoke regression checks that forward bias gives at least as much |I| as
    # zero / reverse bias.  We grab the most-forward and most-reverse bias
    # rows from the converged set rather than assuming any particular order.
    converged = [
        (v, i) for v, i in zip(voltages, currents)
        if math.isfinite(v) and math.isfinite(i)
    ]
    if not converged:
        raise AssertionError("Schottky IV CSV contains no finite converged rows")
    converged.sort(key=lambda item: item[0])
    reverse_abs = abs(converged[0][1])
    forward_abs = abs(converged[-1][1])
    if forward_abs + 1.0e-30 < reverse_abs:
        raise AssertionError(
            "Schottky IV forward |I| should not be smaller than reverse |I|; "
            f"reverse={reverse_abs} forward={forward_abs} "
            f"all voltages={voltages} currents={currents}")
    return {
        "voltages": voltages,
        "currents": currents,
        "reverse_abs": reverse_abs,
        "forward_abs": forward_abs,
    }


def check_iv_trend(example_dir: Path) -> dict[str, Any]:

    cfg = json.loads((example_dir / "simulation.json").read_text())
    expected = expected_sweep_voltages(cfg["sweep"])
    rows = read_csv(example_dir / "outputs" / "pn_iv.csv")
    voltages = [float(curve_column(row, "bias_V", "voltage")) for row in rows]
    if len(voltages) != len(expected):
        raise AssertionError(
            f"PN sweep wrote {len(voltages)} rows, expected {len(expected)}: {voltages}"
        )
    for actual, want in zip(voltages, expected):
        if abs(actual - want) > 1.0e-9:
            raise AssertionError(
                f"PN sweep voltage {actual} does not match expected {want}; all voltages: {voltages}"
            )
    if abs(voltages[-1] - float(cfg["sweep"]["stop"])) > 1.0e-9:
        raise AssertionError(
            f"PN sweep ended at {voltages[-1]} V instead of configured stop "
            f"{cfg['sweep']['stop']} V"
        )

    currents = [abs(float(curve_column(row, "current_total", "total_current"))) for row in rows]
    if currents[-1] + 1e-40 < currents[0]:
        raise AssertionError(f"PN diode current trend is not forward-increasing: {currents}")
    return {"voltages": voltages, "expected_voltages": expected}


def contact_node_ids(mesh: dict[str, Any], name: str) -> list[int]:
    for contact in mesh["contacts"]:
        if contact["name"] == name:
            return [int(node_id) for node_id in contact["node_ids"]]
    raise AssertionError(f"Contact '{name}' not found in mesh")


def contact_biases(cfg: dict[str, Any]) -> dict[str, float]:
    return {contact["name"]: float(contact["bias"]) for contact in cfg["contacts"]}


def average_potential(psi: list[float], node_ids: list[int]) -> float:
    return sum(psi[node_id] for node_id in node_ids) / len(node_ids)


def check_moscap_interface(example_dir: Path) -> dict[str, float]:
    cfg = json.loads((example_dir / "simulation.json").read_text())
    reg = cfg["regression"]
    mesh = json.loads((example_dir / cfg["mesh_file"]).read_text())
    psi = read_vtk_scalar(example_dir / cfg["output_vtk"], "potential_V")
    x = float(reg["x_probe"])
    y_if = float(reg["interface_y"])
    y_si = float(reg["silicon_probe_y"])
    y_ox = float(reg["oxide_probe_y"])

    biases = contact_biases(cfg)
    gate_bias = biases["gate"]
    body_bias = biases["body"]
    expected_bias = gate_bias - body_bias
    contact_tol = float(reg.get("contact_bias_tol", 1.0e-8))
    min_probe_drop = float(reg.get("min_probe_potential_drop", 1.0e-3))

    gate_potential = average_potential(psi, contact_node_ids(mesh, "gate"))
    body_potential = average_potential(psi, contact_node_ids(mesh, "body"))
    contact_drop = gate_potential - body_potential
    if abs(contact_drop - expected_bias) > contact_tol:
        raise AssertionError(
            f"MOSCAP contact potential drop {contact_drop} V does not match configured "
            f"gate/body bias {expected_bias} V"
        )

    i_if_a = node_index(mesh, x, y_if)
    i_if_b = node_index(mesh, x, y_if)
    jump = abs(psi[i_if_a] - psi[i_if_b])
    if jump > float(reg["potential_jump_tol"]):
        raise AssertionError(f"MOSCAP interface potential jump {jump} exceeds tolerance")
    eps_si = 11.7
    eps_ox = 3.9
    si_probe = psi[node_index(mesh, x, y_si)]
    ox_probe = psi[node_index(mesh, x, y_ox)]
    probe_drop = abs(ox_probe - si_probe)
    if probe_drop < min_probe_drop:
        raise AssertionError(f"MOSCAP probe potential drop {probe_drop} V is below {min_probe_drop} V")
    field_si = abs((psi[i_if_a] - si_probe) / (y_if - y_si))
    field_ox = abs((ox_probe - psi[i_if_b]) / (y_ox - y_if))
    disp_si = eps_si * field_si
    disp_ox = eps_ox * field_ox
    denom = max(abs(disp_si), abs(disp_ox), 1e-300)
    rel = abs(disp_si - disp_ox) / denom
    if rel > float(reg["displacement_rel_tol"]):
        raise AssertionError(f"MOSCAP displacement mismatch {rel} exceeds tolerance")
    return {
        "potential_jump": jump,
        "displacement_rel_error": rel,
        "contact_potential_drop": contact_drop,
        "probe_potential_drop": probe_drop,
    }



def set_contact_bias(cfg: dict[str, Any], name: str, bias: float) -> None:
    for contact in cfg["contacts"]:
        if contact["name"] == name:
            contact["bias"] = bias
            return
    raise AssertionError(f"Contact '{name}' not found in config")


def parse_finite_float(row: dict[str, str], column: str, label: str, row_index: int) -> float:
    try:
        value = float(row[column])
    except KeyError as exc:
        raise AssertionError(f"{label} row {row_index} is missing column '{column}'") from exc
    except ValueError as exc:
        raise AssertionError(
            f"{label} row {row_index} has non-numeric {column} value {row.get(column)!r}") from exc
    if not math.isfinite(value):
        raise AssertionError(f"{label} row {row_index} has non-finite {column} value {value}")
    return value


def validate_sweep_voltages(
    actual: list[float],
    expected: list[float],
    label: str,
    tolerance: float = 1.0e-9,
) -> None:
    if len(actual) != len(expected):
        raise AssertionError(f"{label} wrote {len(actual)} rows, expected {len(expected)}: {actual}")
    for actual_value, expected_value in zip(actual, expected):
        if abs(actual_value - expected_value) > tolerance:
            raise AssertionError(
                f"{label} voltage {actual_value} does not match expected {expected_value}; "
                f"all voltages: {actual}")


def assert_monotone_non_decreasing(values: list[float], label: str, tolerance: float = 1.0e-18) -> None:
    for value in values:
        if not math.isfinite(value):
            raise AssertionError(f"{label} contains non-finite value: {values}")
    for left, right in zip(values, values[1:]):
        if right + tolerance < left:
            raise AssertionError(f"{label} is not monotone non-decreasing: {values}")


def check_mos_trends(example_dir: Path, runner: Path) -> dict[str, Any]:
    cfg_path = example_dir / "simulation.json"
    cfg = json.loads(cfg_path.read_text())
    reg = cfg.get("regression", {}).get("mos", {})
    device = reg.get("device", example_dir.name.split("2d", maxsplit=1)[0]).lower()
    if device not in {"nmos", "pmos"}:
        raise AssertionError(f"Unsupported MOS regression device '{device}'")

    rows = read_csv(output_csv_path(example_dir, cfg))
    if not rows:
        raise AssertionError("MOS Id-Vd CSV contains no rows")
    voltages = [
        parse_finite_float(row, "bias_V", "MOS Id-Vd", idx)
        for idx, row in enumerate(rows, start=1)
    ]
    currents = [
        parse_finite_float(row, "current_total", "MOS Id-Vd", idx)
        for idx, row in enumerate(rows, start=1)
    ]
    abs_currents = [abs(value) for value in currents]
    expected_voltages = expected_sweep_voltages(cfg["sweep"])
    validate_sweep_voltages(voltages, expected_voltages, "MOS Id-Vd")
    assert_monotone_non_decreasing(abs_currents, f"{device.upper()} |Id|-Vd trend")

    polarity = float(reg.get("drain_current_sign", 1.0))
    on_current = currents[-1]
    if polarity * on_current <= 0.0:
        raise AssertionError(
            f"{device.upper()} drain current sign {on_current} does not match expected polarity {polarity}")

    idvg = reg.get("idvg", {})
    if not idvg:
        return {"idvd_currents": currents, "idvd_voltages": voltages}

    idvg_cfg = json.loads(cfg_path.read_text())
    drain_bias = float(idvg.get("drain_bias", cfg["sweep"]["stop"]))
    gate_start = float(idvg["gate_start"])
    gate_stop = float(idvg["gate_stop"])
    gate_step = float(idvg["gate_step"])
    set_contact_bias(idvg_cfg, "drain", drain_bias)
    idvg_cfg["output_csv"] = idvg.get("output_csv", f"outputs/{device}_idvg.csv")
    idvg_cfg["sweep"] = dict(idvg_cfg["sweep"])
    idvg_cfg["sweep"].update({
        "contact": "gate",
        "start": gate_start,
        "stop": gate_stop,
        "step": gate_step,
        "current_contact": "drain",
        "write_vtk": False,
    })
    idvg_cfg_path = example_dir / f"{device}_idvg_regression.json"
    idvg_cfg_path.write_text(json.dumps(idvg_cfg, indent=2) + "\n")
    proc = subprocess.run(
        [str(runner), "--config", str(idvg_cfg_path)],
        cwd=example_dir,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise AssertionError(f"{device.upper()} Id-Vg regression run failed: {proc.stderr.strip()}")

    idvg_rows = read_csv(output_csv_path(example_dir, idvg_cfg))
    if not idvg_rows:
        raise AssertionError(f"{device.upper()} Id-Vg CSV contains no rows")
    idvg_gate = [
        parse_finite_float(row, "bias_V", f"{device.upper()} Id-Vg", idx)
        for idx, row in enumerate(idvg_rows, start=1)
    ]
    idvg_currents = [
        parse_finite_float(row, "current_total", f"{device.upper()} Id-Vg", idx)
        for idx, row in enumerate(idvg_rows, start=1)
    ]
    expected_idvg_gate = expected_sweep_voltages(idvg_cfg["sweep"])
    validate_sweep_voltages(idvg_gate, expected_idvg_gate, f"{device.upper()} Id-Vg")
    idvg_abs = [abs(value) for value in idvg_currents]
    assert_monotone_non_decreasing(idvg_abs, f"{device.upper()} |Id|-Vg trend")
    if idvg_abs[-1] <= idvg_abs[0]:
        raise AssertionError(
            f"{device.upper()} Id-Vg on-current did not increase: {idvg_abs}")

    return {
        "device": device,
        "idvd_voltages": voltages,
        "idvd_currents": currents,
        "idvg_gate_voltages": idvg_gate,
        "idvg_currents": idvg_currents,
    }


def run_example(runner: Path, repo: Path, workdir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    name = spec["name"]
    source = spec.get("source")
    config_source = repo / spec["config"]
    source_dir = repo / "examples" / (source or name)
    if config_source.name != "simulation.json" and (source_dir / "simulation.json").exists():
        raise AssertionError(
            f"named-deck example {name!r} must not ship a source simulation.json; "
            "keep simulation.json generated from the selected regression deck")
    example_dir = copy_example(repo, workdir, name, source)
    config = example_dir / config_source.name
    canonical_config = example_dir / "simulation.json"
    if config.resolve() != canonical_config.resolve():
        shutil.copyfile(config, canonical_config)
    run_cfg = json.loads(canonical_config.read_text())
    allow_nonzero_runner_exit = nonzero_runner_exit_allowed(run_cfg)
    proc = subprocess.run([str(runner), "--config", str(config)], cwd=example_dir, text=True, capture_output=True)
    result: dict[str, Any] = {
        "name": name,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
        "passed": False,
        "checks": {},
    }
    try:
        if proc.returncode != 0 and not allow_nonzero_runner_exit:
            raise AssertionError(f"runner exited with {proc.returncode}: {proc.stderr.strip()}")
        if proc.returncode != 0:
            result["runner_nonzero_exit_allowed"] = True
        for rel in spec["expected"]:
            out = example_dir / rel
            if not out.exists() or out.stat().st_size == 0:
                raise AssertionError(f"missing or empty output file: {rel}")
        if "finite_outputs" in spec["checks"]:
            for rel in spec["expected"]:
                if has_nan_or_inf(example_dir / rel):
                    raise AssertionError(f"NaN/Inf detected in {rel}")
            result["checks"]["finite_outputs"] = True
        if "csv_converged" in spec["checks"]:
            check_csv_converged(example_dir)
            result["checks"]["csv_converged"] = True
        if "dc_sweep_regression" in spec["checks"]:
            result["checks"]["dc_sweep_regression"] = check_dc_sweep_regression(example_dir)
        if "iv_trend" in spec["checks"]:
            result["checks"]["iv_trend"] = check_iv_trend(example_dir)
        if "moscap_interface" in spec["checks"]:
            result["checks"]["moscap_interface"] = check_moscap_interface(example_dir)
        if "mos_trends" in spec["checks"]:
            result["checks"]["mos_trends"] = check_mos_trends(example_dir, runner)
        if "schottky_iv_trend" in spec["checks"]:
            result["checks"]["schottky_iv_trend"] = check_schottky_iv_trend(example_dir)
        result["passed"] = True

    except Exception as ex:  # noqa: BLE001 - regression summary should capture all failures.
        result["error"] = str(ex)
    return result


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    runner = (args.runner or default_runner(repo)).resolve()
    workdir = (args.workdir or (repo / "build" / "regression_output")).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    if not runner.exists():
        print(f"Regression runner not found: {runner}", file=sys.stderr)
        return 2
    results = [run_example(runner, repo, workdir, spec) for spec in EXAMPLES]
    summary = {"runner": str(runner), "workdir": str(workdir), "passed": all(r["passed"] for r in results), "examples": results}
    (workdir / "regression_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
