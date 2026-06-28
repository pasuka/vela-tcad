#!/usr/bin/env python3
"""Run PN2D BV sweeps for VanOverstraeten parameter-set candidates."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
E_CHARGE_C = 1.602176634e-19
DEPTH_M = 1.0e-6

CASES = [
    ("BV-P0-default", "default"),
    ("BV-P1-fit-A-only", "sentaurus_fit_A_only"),
    ("BV-P2-fit-A-B", "sentaurus_fit_A_B"),
    ("BV-P3-fit-A-B-switch", "sentaurus_fit_A_B_switch"),
]

DEFAULT_BIASES = [0.0, -1.0, -5.0, -10.0, -16.0, -18.0, -20.0]
WINDOWS = {
    "full_semiconductor": None,
    "junction_window": (0.7, 1.3),
    "center_window": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [
            {str(key).strip(): str(value).strip() for key, value in row.items()}
            for row in csv.DictReader(handle, skipinitialspace=True)
        ]


def write_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def fmt(value: Any) -> str:
    number = finite_float(value)
    if number is None:
        return ""
    return f"{number:.12g}"


def parse_biases(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 9)


def signed_bias_token(bias: float) -> str:
    text = f"{bias:g}".replace("+", "").replace(".", "p")
    return text


def default_runner() -> Path:
    exe = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
    build = REPO / "build" / exe
    if build.exists():
        return build
    return REPO / "build-release" / exe


def resolve_path(base: Path, raw: Any) -> str:
    if raw in (None, ""):
        return ""
    path = Path(str(raw))
    if not path.is_absolute():
        path = base.parent / path
    return str(path.resolve())


def build_case_config(
    base_config: Path,
    case_dir: Path,
    case_id: str,
    parameter_set: str,
    biases: list[float],
) -> Path:
    config = read_json(base_config)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in config:
            config[key] = resolve_path(base_config, config[key])
    config["output_csv"] = f"{case_id}.csv"

    solver = config.setdefault("solver", {})
    impact = solver.setdefault("impact_ionization", {})
    if isinstance(impact, str):
        impact = {"model": impact}
        solver["impact_ionization"] = impact
    impact["model"] = "van_overstraeten"
    impact["parameter_set"] = parameter_set

    sweep = config.setdefault("sweep", {})
    sweep["mode"] = "bv_reverse"
    sweep["bias_points"] = biases
    sweep["start"] = biases[0]
    sweep["stop"] = biases[-1]
    sweep["step"] = -abs(float(sweep.get("step", -5.0)) or 5.0)
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = "vtk/dc_sweep"
    sweep["write_state_file"] = f"{case_id}_last_state.csv"
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["sg_avalanche_edges"] = {
        "enabled": True,
        "csv_file": str((case_dir / f"{case_id}_sg_avalanche_edges.csv").resolve()),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": [sweep.get("current_contact", sweep.get("contact", "Anode"))],
        "csv_file": str((case_dir / f"{case_id}_continuity_balance.csv").resolve()),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((case_dir / f"{case_id}_newton_history.csv").resolve()),
    }

    config["_parameter_matrix_metadata"] = {
        "case_id": case_id,
        "parameter_set": parameter_set,
        "bias_sign_convention": "PN2D BV uses sweep.contact bias_V directly; negative Anode bias is reverse bias in this deck.",
    }

    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "vtk").mkdir(parents=True, exist_ok=True)
    config_path = case_dir / f"{case_id}.json"
    write_json(config_path, config)
    return config_path


def run_case(runner: Path, config: Path, case_dir: Path) -> int:
    print(f":: running {config.name}", flush=True)
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve())],
        cwd=case_dir,
        check=False,
    )
    return result.returncode


class VtkData:
    def __init__(self) -> None:
        self.points: list[tuple[float, float]] = []
        self.cells: list[tuple[int, int, int]] = []
        self.fields: dict[str, list[float]] = {}


def parse_vtk(path: Path) -> VtkData:
    data = VtkData()
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "POINTS":
            count = int(parts[1])
            i += 1
            for _ in range(count):
                xyz = [float(value) for value in lines[i].split()[:3]]
                data.points.append((xyz[0], xyz[1]))
                i += 1
            continue
        if len(parts) >= 3 and parts[0] == "CELLS":
            count = int(parts[1])
            i += 1
            for _ in range(count):
                cell = [int(value) for value in lines[i].split()]
                if len(cell) >= 4 and cell[0] == 3:
                    data.cells.append((cell[1], cell[2], cell[3]))
                i += 1
            continue
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA", "CELL_TYPES",
                }:
                    break
                values.extend(float(value) for value in next_parts)
                i += 1
            data.fields[name] = values
            continue
        if len(parts) >= 3 and parts[0] == "VECTORS":
            name = parts[1]
            i += 1
            values = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA", "CELL_TYPES",
                }:
                    break
                comps = [float(value) for value in next_parts[:3]]
                values.append(math.sqrt(sum(value * value for value in comps)))
                i += 1
            data.fields[name] = values
            continue
        i += 1
    return data


def discover_vtks(case_dir: Path) -> dict[float, Path]:
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    found: dict[float, Path] = {}
    for path in sorted(case_dir.rglob("*.vtk")):
        match = pattern.search(path.name)
        if match:
            found[bias_key(float(match.group("bias")))] = path
    return found


def field_values(data: VtkData, names: list[str]) -> list[float]:
    for name in names:
        values = data.fields.get(name)
        if values:
            return values
    return []


def max_with_coord(data: VtkData, values: list[float], scale: float = 1.0) -> tuple[float | None, float | None, float | None]:
    if not values or not data.points:
        return None, None, None
    indexed = [(abs(value * scale), i) for i, value in enumerate(values) if i < len(data.points)]
    if not indexed:
        return None, None, None
    value, index = max(indexed, key=lambda item: item[0])
    x_m, y_m = data.points[index]
    return value, x_m * 1.0e6, y_m * 1.0e6


def triangle_area_m2(points: list[tuple[float, float]], cell: tuple[int, int, int]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = (points[cell[0]], points[cell[1]], points[cell[2]])
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) * 0.5


def centroid_x_um(points: list[tuple[float, float]], cell: tuple[int, int, int]) -> float:
    return sum(points[node][0] for node in cell) / 3.0 * 1.0e6


def integrate_generation(
    points: list[tuple[float, float]],
    cells: list[tuple[int, int, int]],
    values_m3_s: list[float],
    window: tuple[float, float] | None,
) -> float:
    if not points or not cells or not values_m3_s:
        return 0.0
    total = 0.0
    for cell in cells:
        if max(cell) >= len(values_m3_s):
            continue
        cx = centroid_x_um(points, cell)
        if window is not None and not (window[0] <= cx <= window[1]):
            continue
        avg = sum(values_m3_s[node] for node in cell) / 3.0
        total += avg * triangle_area_m2(points, cell) * DEPTH_M
    return total


def parse_step_used(row: dict[str, str]) -> float | None:
    text = row.get("step_diagnostics", "")
    match = re.search(r"accepted_step=([^;]+)", text)
    return finite_float(match.group(1)) if match else None


def load_reference_current(path: Path) -> dict[float, float]:
    refs: dict[float, float] = {}
    for row in read_rows(path):
        bias = finite_float(row.get("bias_V"))
        current = finite_float(row.get("current_total") or row.get("sentaurus_current_A"))
        if bias is not None and current is not None:
            refs.setdefault(bias_key(bias), current)
    return refs


def load_sentaurus_generation_integrals(root: Path) -> dict[float, dict[str, float]]:
    result: dict[float, dict[str, float]] = {}
    if not root.exists():
        return result
    pattern = re.compile(r"sentaurus_(?P<bias>[-+0-9p.]+)v$")
    for case_dir in root.iterdir():
        if not case_dir.is_dir():
            continue
        match = pattern.match(case_dir.name)
        if not match:
            continue
        raw_bias = match.group("bias").replace("p", ".")
        try:
            bias = float(raw_bias)
        except ValueError:
            continue
        nodes = load_sentaurus_nodes(case_dir / "nodes.csv")
        cells = load_sentaurus_elements(case_dir / "elements.csv")
        values_cm3_s = load_sentaurus_scalar(case_dir / "fields" / "ImpactIonization_region0.csv")
        values_m3_s = [value * 1.0e6 for value in values_cm3_s]
        result[bias_key(bias)] = {
            name: integrate_generation(nodes, cells, values_m3_s, window)
            for name, window in WINDOWS.items()
        }
    return result


def load_sentaurus_nodes(path: Path) -> list[tuple[float, float]]:
    rows = read_rows(path)
    if not rows:
        return []
    points: list[tuple[float, float]] = []
    for row in rows:
        x_um = finite_float(row.get("x_um")) or 0.0
        y_um = finite_float(row.get("y_um")) or 0.0
        points.append((x_um * 1.0e-6, y_um * 1.0e-6))
    return points


def load_sentaurus_elements(path: Path) -> list[tuple[int, int, int]]:
    cells: list[tuple[int, int, int]] = []
    for row in read_rows(path):
        try:
            cells.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
        except (KeyError, ValueError):
            continue
    return cells


def load_sentaurus_scalar(path: Path) -> list[float]:
    rows = read_rows(path)
    if not rows:
        return []
    indexed: list[tuple[int, float]] = []
    for row in rows:
        node = int(row.get("node_id", len(indexed)))
        value = finite_float(row.get("component0") or row.get("value")) or 0.0
        indexed.append((node, value))
    count = max((node for node, _ in indexed), default=-1) + 1
    values = [0.0] * count
    for node, value in indexed:
        values[node] = value
    return values


def summarize_vtk(data: VtkData) -> dict[str, Any]:
    out: dict[str, Any] = {}
    field_specs = {
        "max_ElectricField_V_per_cm": (["ElectricFieldVector", "ElectricField"], 1.0),
        "max_Fn_V_per_cm": (["ElectronHighFieldDrive", "ElectronImpactIonizationDrive"], 1.0),
        "max_Fp_V_per_cm": (["HoleHighFieldDrive", "HoleImpactIonizationDrive"], 1.0),
        "max_eAlphaAvalanche_cm_inv": (["ElectronAlphaAvalanche"], 0.01),
        "max_hAlphaAvalanche_cm_inv": (["HoleAlphaAvalanche"], 0.01),
        "max_AvalancheGeneration_m3_per_s": (["AvalancheGeneration"], 1.0),
    }
    for label, (names, scale) in field_specs.items():
        value, x_um, y_um = max_with_coord(data, field_values(data, names), scale)
        out[label] = value
        out[f"{label}_x_um"] = x_um
        out[f"{label}_y_um"] = y_um

    generation = field_values(data, ["AvalancheGeneration"])
    for name, window in WINDOWS.items():
        integral = integrate_generation(data.points, data.cells, generation, window)
        out[f"integral_AvalancheGeneration_{name}_per_s_per_um"] = integral
        out[f"q_integral_AvalancheGeneration_{name}_A_per_um"] = integral * E_CHARGE_C
    return out


def nearest_vtk(vtks: dict[float, Path], bias: float) -> Path | None:
    return vtks.get(bias_key(bias))


def collect_case_rows(
    case_root: Path,
    case_id: str,
    parameter_set: str,
    biases: list[float],
    reference_current: dict[float, float],
    sentaurus_generation: dict[float, dict[str, float]],
    run_return_code: int,
) -> list[dict[str, Any]]:
    case_dir = case_root / case_id
    sweep_csv = case_dir / f"{case_id}.csv"
    curve_rows = {bias_key(finite_float(row.get("bias_V")) or 0.0): row for row in read_rows(sweep_csv)}
    vtks = discover_vtks(case_dir)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        row = curve_rows.get(bias_key(bias), {})
        vtk_path = nearest_vtk(vtks, bias)
        vtk_summary = summarize_vtk(parse_vtk(vtk_path)) if vtk_path else {}
        sent_q = sentaurus_generation.get(bias_key(bias), {})
        item: dict[str, Any] = {
            "case_id": case_id,
            "parameter_set": parameter_set,
            "bias_V": bias,
            "bias_sign_convention": "negative Anode bias is reverse bias for this PN2D deck",
            "run_return_code": run_return_code,
            "convergence_status": "converged" if row.get("converged") == "1" else "not_converged",
            "newton_iterations": row.get("newton_iterations", ""),
            "iterations": row.get("iterations", ""),
            "voltage_step_used_V": parse_step_used(row),
            "terminal_electron_current_A": finite_float(row.get("current_electron")),
            "terminal_hole_current_A": finite_float(row.get("current_hole")),
            "terminal_total_current_A": finite_float(row.get("current_total")),
            "terminal_electron_current_A_per_um": finite_float(row.get("current_electron_A_per_um")),
            "terminal_hole_current_A_per_um": finite_float(row.get("current_hole_A_per_um")),
            "terminal_total_current_A_per_um": finite_float(row.get("current_total_A_per_um")),
            "sentaurus_terminal_total_current_A": reference_current.get(bias_key(bias)),
            "vtk_file": str(vtk_path) if vtk_path else "",
        }
        item.update(vtk_summary)
        for name in WINDOWS:
            sent_integral = sent_q.get(name)
            item[f"sentaurus_integral_AvalancheGeneration_{name}_per_s_per_um"] = sent_integral
            item[f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um"] = (
                sent_integral * E_CHARGE_C if sent_integral is not None else None
            )
        rows.append(item)
    return rows


def log_score(rows: list[dict[str, Any]], value_col: str, ref_col: str) -> float | None:
    scores: list[float] = []
    for row in rows:
        value = finite_float(row.get(value_col))
        ref = finite_float(row.get(ref_col))
        if value is None or ref is None or value == 0.0 or ref == 0.0:
            continue
        scores.append(abs(math.log10(abs(value) / abs(ref))))
    return median(scores)


def median(values: list[float]) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def ratio_median(rows: list[dict[str, Any]], value_col: str, ref_col: str) -> float | None:
    ratios = []
    for row in rows:
        value = finite_float(row.get(value_col))
        ref = finite_float(row.get(ref_col))
        if value is not None and ref is not None and ref != 0.0:
            ratios.append(abs(value) / abs(ref))
    return median(ratios)


def build_summary_rows(by_bias: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case_id, parameter_set in CASES:
        subset = [row for row in by_bias if row["case_id"] == case_id]
        rows.append({
            "case_id": case_id,
            "parameter_set": parameter_set,
            "bias_count": len(subset),
            "converged_count": sum(1 for row in subset if row["convergence_status"] == "converged"),
            "median_current_log10_abs_error": log_score(
                subset, "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A"),
            "median_qAvalanche_full_log10_abs_error": log_score(
                subset,
                "q_integral_AvalancheGeneration_full_semiconductor_A_per_um",
                "sentaurus_q_integral_AvalancheGeneration_full_semiconductor_A_per_um",
            ),
            "median_qAvalanche_junction_log10_abs_error": log_score(
                subset,
                "q_integral_AvalancheGeneration_junction_window_A_per_um",
                "sentaurus_q_integral_AvalancheGeneration_junction_window_A_per_um",
            ),
            "median_qAvalanche_left_shoulder_log10_abs_error": log_score(
                subset,
                "q_integral_AvalancheGeneration_left_shoulder_A_per_um",
                "sentaurus_q_integral_AvalancheGeneration_left_shoulder_A_per_um",
            ),
            "median_qAvalanche_right_shoulder_log10_abs_error": log_score(
                subset,
                "q_integral_AvalancheGeneration_right_shoulder_A_per_um",
                "sentaurus_q_integral_AvalancheGeneration_right_shoulder_A_per_um",
            ),
            "median_current_abs_ratio_to_sentaurus": ratio_median(
                subset, "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A"),
            "median_qAvalanche_full_abs_ratio_to_sentaurus": ratio_median(
                subset,
                "q_integral_AvalancheGeneration_full_semiconductor_A_per_um",
                "sentaurus_q_integral_AvalancheGeneration_full_semiconductor_A_per_um",
            ),
        })
    return rows


def best_case(summary_rows: list[dict[str, Any]], score_col: str) -> str:
    best = ""
    best_score = math.inf
    for row in summary_rows:
        score = finite_float(row.get(score_col))
        if score is not None and score < best_score:
            best = str(row["parameter_set"])
            best_score = score
    return best or "insufficient_reference_data"


def best_complete_case(summary_rows: list[dict[str, Any]], score_col: str) -> str:
    complete = [
        row for row in summary_rows
        if row.get("bias_count") == row.get("converged_count")
    ]
    return best_case(complete or summary_rows, score_col)


def row_for(summary_rows: list[dict[str, Any]], parameter_set: str) -> dict[str, Any]:
    for row in summary_rows:
        if row["parameter_set"] == parameter_set:
            return row
    return {}


def improves(a: Any, b: Any, margin: float = 0.1) -> bool:
    av = finite_float(a)
    bv = finite_float(b)
    return av is not None and bv is not None and av + margin < bv


def write_markdown(path: Path, summary_rows: list[dict[str, Any]], by_bias: list[dict[str, Any]]) -> None:
    default = row_for(summary_rows, "default")
    a_only = row_for(summary_rows, "sentaurus_fit_A_only")
    ab = row_for(summary_rows, "sentaurus_fit_A_B")
    ab_switch = row_for(summary_rows, "sentaurus_fit_A_B_switch")
    best_current_partial = best_case(summary_rows, "median_current_log10_abs_error")
    best_q_partial = best_case(summary_rows, "median_qAvalanche_full_log10_abs_error")
    best_current = best_complete_case(summary_rows, "median_current_log10_abs_error")
    best_q = best_complete_case(summary_rows, "median_qAvalanche_full_log10_abs_error")

    a_only_current = improves(
        a_only.get("median_current_log10_abs_error"),
        default.get("median_current_log10_abs_error"),
    )
    a_only_q = improves(
        a_only.get("median_qAvalanche_full_log10_abs_error"),
        default.get("median_qAvalanche_full_log10_abs_error"),
    )
    ab_shoulders = improves(
        ab.get("median_qAvalanche_left_shoulder_log10_abs_error"),
        a_only.get("median_qAvalanche_left_shoulder_log10_abs_error"),
        0.05,
    ) or improves(
        ab.get("median_qAvalanche_right_shoulder_log10_abs_error"),
        a_only.get("median_qAvalanche_right_shoulder_log10_abs_error"),
        0.05,
    )
    switch_center = improves(
        ab_switch.get("median_qAvalanche_junction_log10_abs_error"),
        ab.get("median_qAvalanche_junction_log10_abs_error"),
        0.05,
    )
    switch_shoulder_worse = (
        improves(ab.get("median_qAvalanche_left_shoulder_log10_abs_error"),
                 ab_switch.get("median_qAvalanche_left_shoulder_log10_abs_error"), 0.05)
        or improves(ab.get("median_qAvalanche_right_shoulder_log10_abs_error"),
                    ab_switch.get("median_qAvalanche_right_shoulder_log10_abs_error"), 0.05)
    )
    ab_ratio = finite_float(ab.get("median_qAvalanche_full_abs_ratio_to_sentaurus"))
    ab_still_low = ab_ratio is not None and ab_ratio < 0.5
    ab_incomplete = ab.get("bias_count") != ab.get("converged_count")

    lines = [
        "# VanOverstraeten BV Parameter Matrix",
        "",
        "- bias sign convention: negative Anode bias is reverse bias for this PN2D deck.",
        "- parameter_set values are written into each generated case config; default remains the built-in material parameter set.",
        "- q integral columns use q * integral(AvalancheGeneration dV) over a 1 um depth.",
        "",
        "## Summary",
        "",
        "| case | parameter_set | converged/biases | current score | qAvalanche full score | qAvalanche full ratio |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['case_id']} | {row['parameter_set']} | "
            f"{row['converged_count']}/{row['bias_count']} | "
            f"{fmt(row.get('median_current_log10_abs_error'))} | "
            f"{fmt(row.get('median_qAvalanche_full_log10_abs_error'))} | "
            f"{fmt(row.get('median_qAvalanche_full_abs_ratio_to_sentaurus'))} |"
        )
    lines.extend([
        "",
        "## Answers",
        "",
        f"1. Closest complete terminal current curve: {best_current}. Partial converged-point winner: {best_current_partial}.",
        f"2. Closest complete q integral AvalancheGeneration dV: {best_q}. Partial converged-point winner: {best_q_partial}.",
        f"3. fit_A_only mainly fixes current scale: {'yes' if a_only_current and not a_only_q else 'mixed/no'}.",
        f"4. fit_A_B improves shoulder source integrals: {'yes' if ab_shoulders else 'no or insufficient reference data'}.",
        f"5. fit_A_B_switch improves center high field but worsens other regions: {'yes' if switch_center and switch_shoulder_worse else 'no or not clearly'}.",
    ])
    if best_current == "sentaurus_fit_A_B" and best_q == "sentaurus_fit_A_B":
        lines.append("6. Remaining difference is now more likely source mapping or current density than alpha parameters alone.")
        lines.append("7. sentaurus_fit_A_B is the recommended temporary calibration parameter_set; enter the next stage source mapping check.")
    elif ab_incomplete:
        lines.append("6. sentaurus_fit_A_B improves low-bias terminal-current scale but does not complete the BV sweep; Newton path and source/current feedback are now prominent.")
        lines.append("7. Do not promote sentaurus_fit_A_B as a BV calibration yet; inspect source mapping/current density feedback and continuation stability first.")
    elif ab_still_low:
        lines.append("6. Remaining difference is still significantly low after fit_A_B; continue checking cutoff/smoothing/RefDens or alpha calculation location.")
        lines.append("7. Do not promote sentaurus_fit_A_B beyond temporary diagnostics yet.")
    else:
        lines.append("6. Remaining difference is mixed across alpha parameters, source mapping, current density, cutoff/smoothing/RefDens, and Newton path.")
        lines.append("7. Use the by-bias CSV to decide whether sentaurus_fit_A_B is good enough as a temporary calibration set.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append(f"- by-bias rows: {path.parent / 'vanoverstraeten_bv_parameter_matrix_by_bias.csv'}")
    lines.append(f"- summary rows: {path.parent / 'vanoverstraeten_bv_parameter_matrix_summary.csv'}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def output_fields() -> list[str]:
    base = [
        "case_id", "parameter_set", "bias_V", "bias_sign_convention", "run_return_code",
        "convergence_status", "newton_iterations", "iterations", "voltage_step_used_V",
        "terminal_electron_current_A", "terminal_hole_current_A", "terminal_total_current_A",
        "terminal_electron_current_A_per_um", "terminal_hole_current_A_per_um",
        "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A",
        "max_ElectricField_V_per_cm", "max_ElectricField_V_per_cm_x_um",
        "max_ElectricField_V_per_cm_y_um",
        "max_Fn_V_per_cm", "max_Fn_V_per_cm_x_um", "max_Fn_V_per_cm_y_um",
        "max_Fp_V_per_cm", "max_Fp_V_per_cm_x_um", "max_Fp_V_per_cm_y_um",
        "max_eAlphaAvalanche_cm_inv", "max_eAlphaAvalanche_cm_inv_x_um",
        "max_eAlphaAvalanche_cm_inv_y_um",
        "max_hAlphaAvalanche_cm_inv", "max_hAlphaAvalanche_cm_inv_x_um",
        "max_hAlphaAvalanche_cm_inv_y_um",
        "max_AvalancheGeneration_m3_per_s", "max_AvalancheGeneration_m3_per_s_x_um",
        "max_AvalancheGeneration_m3_per_s_y_um",
    ]
    for name in WINDOWS:
        base.append(f"integral_AvalancheGeneration_{name}_per_s_per_um")
        base.append(f"q_integral_AvalancheGeneration_{name}_A_per_um")
        base.append(f"sentaurus_integral_AvalancheGeneration_{name}_per_s_per_um")
        base.append(f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um")
    base.append("vtk_file")
    return base


def summary_fields() -> list[str]:
    return [
        "case_id",
        "parameter_set",
        "bias_count",
        "converged_count",
        "median_current_log10_abs_error",
        "median_qAvalanche_full_log10_abs_error",
        "median_qAvalanche_junction_log10_abs_error",
        "median_qAvalanche_left_shoulder_log10_abs_error",
        "median_qAvalanche_right_shoulder_log10_abs_error",
        "median_current_abs_ratio_to_sentaurus",
        "median_qAvalanche_full_abs_ratio_to_sentaurus",
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare/simulation_coarse_previous_full20_aligned.json",
    )
    parser.add_argument("--runner", type=Path, default=default_runner())
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO / "build/diagnostics",
    )
    parser.add_argument(
        "--reference-current-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_vm_vector/reference_curves/pn2d_sentaurus2018_coarse7x3_bv_reference.csv",
    )
    parser.add_argument(
        "--sentaurus-multibias-dir",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare/sentaurus_multibias",
    )
    parser.add_argument(
        "--bias-points",
        default=",".join(f"{bias:g}" for bias in DEFAULT_BIASES),
    )
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = parse_biases(args.bias_points)
    if not biases:
        raise ValueError("--bias-points must contain at least one bias")
    if biases[0] != 0.0:
        biases.insert(0, 0.0)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    case_root = args.out_dir / "vanoverstraeten_bv_parameter_matrix_cases"
    if args.skip_run:
        test_case_root = args.out_dir / "cases"
        if test_case_root.exists():
            case_root = test_case_root

    reference_current = load_reference_current(args.reference_current_csv)
    sentaurus_generation = load_sentaurus_generation_integrals(args.sentaurus_multibias_dir)

    all_rows: list[dict[str, Any]] = []
    for case_id, parameter_set in CASES:
        case_dir = case_root / case_id
        config_path = build_case_config(args.base_config, case_dir, case_id, parameter_set, biases)
        run_return_code = 0
        if not args.skip_run:
            run_return_code = run_case(args.runner, config_path, case_dir)
        all_rows.extend(
            collect_case_rows(
                case_root,
                case_id,
                parameter_set,
                biases,
                reference_current,
                sentaurus_generation,
                run_return_code,
            )
        )

    by_bias_path = args.out_dir / "vanoverstraeten_bv_parameter_matrix_by_bias.csv"
    summary_path = args.out_dir / "vanoverstraeten_bv_parameter_matrix_summary.csv"
    md_path = args.out_dir / "vanoverstraeten_bv_parameter_matrix.md"
    summary_rows = build_summary_rows(all_rows)
    write_rows(by_bias_path, all_rows, output_fields())
    write_rows(summary_path, summary_rows, summary_fields())
    write_markdown(md_path, summary_rows, all_rows)
    print(by_bias_path)
    print(summary_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
