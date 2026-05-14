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
        "name": "pn_diode",
        "config": Path("examples/pn_diode/simulation.json"),
        "expected": [Path("outputs/pn_iv.csv"), Path("outputs/pn_sweep_0000_0V.vtk")],
        "checks": ["csv_converged", "finite_outputs", "iv_trend", "dc_sweep_regression"],
    },
    {
        "name": "moscap",
        "config": Path("examples/moscap/simulation.json"),
        "expected": [Path("outputs/moscap_poisson.vtk")],
        "checks": ["finite_outputs", "moscap_interface"],
    },
    {
        "name": "nmos2d",
        "config": Path("examples/nmos2d/simulation.json"),
        "expected": [Path("outputs/nmos_poisson.vtk")],
        "checks": ["finite_outputs"],
    },
    {
        "name": "nmos2d_dd",
        "config": Path("examples/nmos2d_dd/simulation.json"),
        "expected": [
            Path("outputs/nmos2d_dd_iv.csv"),
            Path("outputs/nmos2d_dd_sweep_0000_0V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression", "mos_trends"],
    },
    {
        "name": "pmos2d_dd",
        "config": Path("examples/pmos2d_dd/simulation.json"),
        "expected": [
            Path("outputs/pmos2d_dd_iv.csv"),
            Path("outputs/pmos2d_dd_sweep_0000_0V.vtk"),
        ],
        "checks": ["csv_converged", "finite_outputs", "dc_sweep_regression", "mos_trends"],
    },
    {
        "name": "ldmos2d_poisson",
        "config": Path("examples/ldmos2d_poisson/simulation.json"),
        "expected": [Path("outputs/ldmos2d_poisson.vtk")],
        "checks": ["finite_outputs"],
    },
    {
        "name": "igbt2d_poisson",
        "config": Path("examples/igbt2d_poisson/simulation.json"),
        "expected": [Path("outputs/igbt2d_poisson.vtk")],
        "checks": ["finite_outputs"],
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


def copy_example(repo: Path, workdir: Path, name: str) -> Path:
    src = repo / "examples" / name
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


def check_csv_converged(example_dir: Path) -> None:
    cfg = json.loads((example_dir / "simulation.json").read_text())
    rows = read_csv(output_csv_path(example_dir, cfg))
    if not rows:
        raise AssertionError("DC sweep CSV contains no sweep rows")
    declared = bool(cfg.get("regression", {}).get("declared_converged", True))
    if declared:
        bad = [row for row in rows if row.get("converged") != "1"]
        if bad:
            raise AssertionError(f"Non-converged declared sweep rows: {bad}")


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
    reg = cfg.get("regression", {}).get("dc_sweep", {})
    rows = read_csv(output_csv_path(example_dir, cfg))
    if not rows:
        raise AssertionError("DC sweep CSV contains no sweep rows")

    expected_rows = reg.get("expected_rows")
    if expected_rows is None and "sweep" in cfg:
        expected_rows = len(expected_sweep_voltages(cfg["sweep"]))
    if expected_rows is not None and len(rows) != int(expected_rows):
        raise AssertionError(f"DC sweep wrote {len(rows)} rows, expected {expected_rows}")

    max_abs_attempted = float(reg.get("max_abs_attempted_step", math.inf))
    max_abs_accepted = float(reg.get("max_abs_accepted_step", math.inf))
    max_retry = int(reg.get("max_retry_count", 0))
    modern_required = ("mode", "bias_contact", "bias_V", "current_contact",
                       "current_electron", "current_hole", "current_total",
                       "converged", "iterations", "step_diagnostics")
    missing_modern = [column for column in modern_required if column not in rows[0]]
    if missing_modern:
        raise AssertionError(
            "DC sweep CSV is missing required curve schema column(s): "
            + ", ".join(missing_modern))

    mode = cfg.get("sweep", {}).get("mode", "iv")
    if mode == "cv_quasistatic":
        for column in ("charge_C_per_m", "capacitance_F_per_m"):
            if column not in rows[0]:
                raise AssertionError(f"CV sweep CSV is missing required column '{column}'")
    if mode == "bv_reverse":
        for column in ("max_electric_field_V_per_m", "current_jump_ratio",
                       "breakdown_detected", "breakdown_voltage", "criterion"):
            if column not in rows[0]:
                raise AssertionError(f"BV sweep CSV is missing required column '{column}'")

    max_attempted_seen = 0.0
    max_accepted_seen = 0.0
    max_retry_seen = 0
    for row_index, row in enumerate(rows, start=1):
        try:
            diagnostics = parse_step_diagnostics(row["step_diagnostics"])
            attempted = abs(float(diagnostics["attempted_step"]))
            accepted = abs(float(diagnostics["accepted_step"]))
            retry = int(diagnostics["retry_count"])
        except (KeyError, ValueError) as exc:
            raise AssertionError(
                f"DC sweep row {row_index} has unparseable step diagnostics: {row}") from exc
        max_attempted_seen = max(max_attempted_seen, attempted)
        max_accepted_seen = max(max_accepted_seen, accepted)
        max_retry_seen = max(max_retry_seen, retry)
        if mode == "cv_quasistatic":
            _ = parse_finite_float(row, "charge_C_per_m", "CV sweep", row_index)
            capacitance = parse_finite_float(row, "capacitance_F_per_m", "CV sweep", row_index)
            if row_index > 1 and abs(capacitance) <= 0.0:
                raise AssertionError(f"CV sweep row {row_index} has zero differential capacitance")
        if mode == "bv_reverse":
            max_field = parse_finite_float(row, "max_electric_field_V_per_m", "BV sweep", row_index)
            if max_field < 0.0:
                raise AssertionError(f"BV sweep row {row_index} has negative max electric field")
            jump = parse_finite_float(row, "current_jump_ratio", "BV sweep", row_index)
            if jump < 0.0:
                raise AssertionError(f"BV sweep row {row_index} has negative current jump ratio")
            if row.get("breakdown_detected") == "1" and not row.get("criterion"):
                raise AssertionError(f"BV sweep row {row_index} detected breakdown without a criterion")
        if attempted > max_abs_attempted + 1.0e-12:
            raise AssertionError(
                f"attempted_step {attempted} exceeds regression limit {max_abs_attempted}")
        if accepted > max_abs_accepted + 1.0e-12:
            raise AssertionError(
                f"accepted_step {accepted} exceeds regression limit {max_abs_accepted}")
        if retry > max_retry:
            raise AssertionError(f"retry_count {retry} exceeds regression limit {max_retry}")

    return {
        "rows": len(rows),
        "max_abs_attempted_step": max_attempted_seen,
        "max_abs_accepted_step": max_accepted_seen,
        "max_retry_count": max_retry_seen,
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
    example_dir = copy_example(repo, workdir, name)
    config = example_dir / "simulation.json"
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
        if proc.returncode != 0:
            raise AssertionError(f"runner exited with {proc.returncode}: {proc.stderr.strip()}")
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
