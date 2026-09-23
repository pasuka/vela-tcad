#!/usr/bin/env python3
"""Run and close the two-stage Charon/Sentaurus Schottky reference workflow."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "reference_tcad" / "schottky_charon_sentaurus2018"
REFERENCE = FIXTURE / "sentaurus_forward_reference.csv"


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(text)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def converged_rows(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    return [row for row in rows if row.get("converged", "1") == "1"]


def _absolute(path: Path) -> str:
    return str(path.resolve())


def materialize(output_root: Path) -> tuple[Path, Path]:
    output_root = output_root.resolve()
    stage_a_dir = output_root / "stage_a"
    stage_b_dir = output_root / "stage_b"
    stage_a_dir.mkdir(parents=True, exist_ok=True)
    stage_b_dir.mkdir(parents=True, exist_ok=True)

    stage_a = json.loads(
        (FIXTURE / "vela" / "simulation_iv.json").read_text(encoding="utf-8"))
    stage_b = json.loads(
        (FIXTURE / "vela" / "simulation_iv_arclength.json").read_text(
            encoding="utf-8"))
    for config in (stage_a, stage_b):
        config["mesh_file"] = _absolute(FIXTURE / "vela" / "mesh.json")
        config["node_doping_file"] = _absolute(FIXTURE / "vela" / "doping.csv")
        config["solver"]["global_continuity_closure"] = {
            # Retain the global audit columns, but do not make their near-zero
            # source normalization a nonlinear acceptance gate.  Terminal KCL
            # is checked below with an absolute-or-relative hybrid rule.
            "mode": "report", "tolerance": 0.01, "source_floor": 1.0e-14}

    stage_a["output_csv"] = _absolute(stage_a_dir / "sweep.csv")
    stage_a["sweep"]["write_state_every_point_prefix"] = _absolute(
        stage_a_dir / "states" / "schottky")
    stage_a["sweep"]["write_state_file"] = _absolute(stage_a_dir / "last_state.h5")

    stage_b["output_csv"] = _absolute(stage_b_dir / "sweep.csv")
    stage_b["sweep"]["initial_state_file"] = _absolute(
        stage_a_dir / "states" / "schottky_bias_0p820000.h5")
    stage_b["sweep"]["write_state_every_point_prefix"] = _absolute(
        stage_b_dir / "states" / "schottky")
    stage_b["sweep"]["write_state_file"] = _absolute(stage_b_dir / "last_state.h5")

    path_a = stage_a_dir / "simulation.json"
    path_b = stage_b_dir / "simulation.json"
    write_text_lf(path_a, json.dumps(stage_a, indent=2) + "\n")
    write_text_lf(path_b, json.dumps(stage_b, indent=2) + "\n")
    return path_a, path_b


def run_stage(runner: Path, config: Path) -> None:
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve())],
        cwd=str(ROOT), text=True, capture_output=True, check=False)
    write_text_lf(config.parent / "runner.stdout.log", result.stdout)
    write_text_lf(config.parent / "runner.stderr.log", result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {config} with {result.returncode}: {result.stderr[-1000:]}")


def merge_curves(stage_a_csv: Path, stage_b_csv: Path, output: Path) -> list[dict[str, float]]:
    a = converged_rows(stage_a_csv)
    b = converged_rows(stage_b_csv)
    points: list[dict[str, float]] = []
    for row in a:
        points.append({
            "bias_V": float(row["bias_V"]),
            "current_total_A_per_um": float(row["current_total_A_per_um"]),
        })
    last_a = points[-1]["bias_V"]
    for row in b:
        bias = float(row["bias_V"])
        if bias > last_a + 1.0e-12:
            points.append({
                "bias_V": bias,
                "current_total_A_per_um": float(row["current_total_A_per_um"]),
            })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["bias_V", "current_total_A_per_um"])
        writer.writeheader()
        writer.writerows(points)
    return points


def _log_interpolate(points: list[dict[str, float]], bias: float) -> float:
    ordered = sorted(points, key=lambda point: point["bias_V"])
    for left, right in zip(ordered, ordered[1:]):
        if left["bias_V"] <= bias <= right["bias_V"]:
            if right["bias_V"] == left["bias_V"]:
                return abs(right["current_total_A_per_um"])
            fraction = ((bias - left["bias_V"])
                        / (right["bias_V"] - left["bias_V"]))
            current0 = left["current_total_A_per_um"]
            current1 = right["current_total_A_per_um"]
            if current0 != 0.0 and current1 != 0.0 and math.copysign(
                    1.0, current0) == math.copysign(1.0, current1):
                y0 = math.log10(abs(current0))
                y1 = math.log10(abs(current1))
                return 10.0 ** (y0 + fraction * (y1 - y0))
            # Around the equilibrium sign crossing, the reference comparator
            # deliberately falls back to signed linear interpolation.
            return abs(current0 + fraction * (current1 - current0))
    raise ValueError(f"candidate curve does not bracket {bias:g} V")


def compare_curves(candidate: list[dict[str, float]]) -> dict[str, Any]:
    reference = [
        {"bias_V": float(row["bias_V"]),
         "current_total_A_per_um": float(row["current_total_A_per_um"])}
        for row in read_csv(REFERENCE)
        if float(row["bias_V"]) > 0.0
    ]
    comparisons = []
    for row in reference:
        vela = _log_interpolate(candidate, row["bias_V"])
        sentaurus = abs(row["current_total_A_per_um"])
        comparisons.append({
            "bias_V": row["bias_V"],
            "sentaurus_current_A_per_um": sentaurus,
            "vela_current_A_per_um": vela,
            "absolute_log10_error_dex": abs(math.log10(vela / sentaurus)),
            "relative_error": abs(vela - sentaurus) / sentaurus,
        })
    positive = [point for point in candidate if point["bias_V"] >= 0.01]
    monotonic = all(
        abs(right["current_total_A_per_um"])
        >= abs(left["current_total_A_per_um"]) * (1.0 - 1.0e-10)
        for left, right in zip(positive, positive[1:]))
    at_one = next(row for row in comparisons if abs(row["bias_V"] - 1.0) < 1.0e-12)
    maximum = max(row["absolute_log10_error_dex"] for row in comparisons)
    return {
        "status": "pass" if len(comparisons) == 24 and monotonic and maximum <= 0.5 else "fail",
        "points_compared": len(comparisons),
        "bias_range_V": [comparisons[0]["bias_V"], comparisons[-1]["bias_V"]],
        "trend_match": monotonic,
        "maximum_log10_current_error_dex": maximum,
        "threshold_dex": 0.5,
        "one_volt": at_one,
        "points": comparisons,
    }


def audit_stage_b(rows: list[dict[str, str]]) -> dict[str, Any]:
    biases = [float(row["bias_V"]) for row in rows]
    currents = [abs(float(row["current_total_A_per_um"])) for row in rows]
    backsteps = sum(right < left - 1.0e-12 for left, right in zip(biases, biases[1:]))
    decreases = sum(right < left * (1.0 - 1.0e-10)
                    for left, right in zip(currents, currents[1:]))
    return {
        "points": len(rows),
        "last_bias_V": biases[-1],
        "voltage_backsteps": backsteps,
        "current_decreases": decreases,
        "status": "pass" if biases[-1] >= 1.0 and backsteps == 0 and decreases == 0 else "fail",
    }


def audit_kcl(rows: list[dict[str, str]]) -> dict[str, Any]:
    requested = (0.0, 0.4, 1.0)
    points = []
    for target in requested:
        row = min(rows, key=lambda item: abs(float(item["bias_V"]) - target))
        actual = float(row["bias_V"])
        current = float(row["current_total_A_per_um"])
        electron = float(row["current_electron_A_per_um"])
        hole = float(row["current_hole_A_per_um"])
        # Vela stores conventional terminal current as electron current minus
        # hole-particle current; preserve that sign convention in the KCL sum.
        residual = abs(current - electron + hole)
        relative = residual / max(abs(current), 1.0e-300)
        terminal_kcl = residual <= 1.0e-14 or relative <= 0.01
        forward_sign = target == 0.0 or current > 0.0
        points.append({
            "requested_bias_V": target,
            "actual_bias_V": actual,
            "current_total_A_per_um": current,
            "electron_current_A_per_um": electron,
            "hole_current_A_per_um": hole,
            "terminal_kcl_absolute_residual_A_per_um": residual,
            "terminal_kcl_relative_residual": relative,
            "global_electron_continuity_ratio_diagnostic": float(
                row["global_electron_continuity_closure_ratio"]),
            "global_hole_continuity_ratio_diagnostic": float(
                row["global_hole_continuity_closure_ratio"]),
            "forward_current_sign_ok": forward_sign,
            "status": "pass" if terminal_kcl and forward_sign else "fail",
        })
    return {
        "rule": "terminal |It-Ie+Ih| <= 1e-14 A/um OR relative residual <= 1%; forward current positive",
        "points": points,
        "status": "pass" if all(point["status"] == "pass" for point in points) else "fail",
    }


def closeout(output_root: Path) -> dict[str, Any]:
    stage_a_csv = output_root / "stage_a" / "sweep.csv"
    stage_b_csv = output_root / "stage_b" / "sweep.csv"
    candidate_path = output_root / "vela_schottky_iv_combined.csv"
    candidate = merge_curves(stage_a_csv, stage_b_csv, candidate_path)
    stage_a_rows = converged_rows(stage_a_csv)
    stage_b_rows = converged_rows(stage_b_csv)
    curve = compare_curves(candidate)
    continuation = audit_stage_b(stage_b_rows)
    kcl = audit_kcl(stage_a_rows + stage_b_rows)
    result = {
        "schema": "vela.schottky.reference_workflow.v1",
        "status": "pass" if all(
            item["status"] == "pass" for item in (curve, continuation, kcl)) else "fail",
        "physics": ["Poisson", "Electron", "Hole", "SRH", "thermionic Robin"],
        "excluded_physics": [
            "barrier lowering", "tunnelling", "high-field mobility",
            "self heating", "AC", "series resistance"],
        "stage_a": {"points": len(stage_a_rows), "last_bias_V": float(stage_a_rows[-1]["bias_V"])},
        "stage_b": continuation,
        "curve_acceptance": curve,
        "three_bias_boundary_kcl": kcl,
        "artifacts": {
            "stage_a_sha256": sha256(stage_a_csv),
            "stage_b_sha256": sha256(stage_b_csv),
            "combined_sha256": sha256(candidate_path),
        },
    }
    write_text_lf(
        output_root / "validation.json", json.dumps(result, indent=2) + "\n")
    write_markdown(output_root / "validation.md", result)
    return result


def write_markdown(path: Path, result: dict[str, Any]) -> None:
    curve = result["curve_acceptance"]
    one = curve["one_volt"]
    lines = [
        "# Simple Schottky reference closeout", "",
        f"Status: **{result['status'].upper()}**", "",
        f"Compared {curve['points_compared']} Sentaurus points; maximum error "
        f"is `{curve['maximum_log10_current_error_dex']:.6g} dex`.", "",
        f"At 1 V, Vela is `{one['vela_current_A_per_um']:.9g} A/um` and "
        f"Sentaurus is `{one['sentaurus_current_A_per_um']:.9g} A/um` "
        f"(`{100.0 * one['relative_error']:.4f}%`).", "",
        "| Requested bias | Actual bias | Terminal KCL abs / relative | Total current | Result |",
        "| ---: | ---: | ---: | ---: | --- |",
    ]
    for point in result["three_bias_boundary_kcl"]["points"]:
        lines.append(
            f"| {point['requested_bias_V']:.3g} V | {point['actual_bias_V']:.9g} V | "
            f"{point['terminal_kcl_absolute_residual_A_per_um']:.3g} A/um / "
            f"{point['terminal_kcl_relative_residual']:.3g} | "
            f"{point['current_total_A_per_um']:.6g} A/um | {point['status']} |")
    lines.extend(["", "No optional Schottky physics is enabled.", ""])
    write_text_lf(path, "\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--publish-json", type=Path)
    parser.add_argument("--publish-md", type=Path)
    args = parser.parse_args()
    path_a, path_b = materialize(args.output_root)
    if args.execute:
        if args.runner is None or not args.runner.is_file():
            parser.error("--execute requires an existing --runner")
        run_stage(args.runner, path_a)
        run_stage(args.runner, path_b)
    result = closeout(args.output_root.resolve())
    if args.publish_json is not None:
        args.publish_json.parent.mkdir(parents=True, exist_ok=True)
        write_text_lf(args.publish_json, json.dumps(result, indent=2) + "\n")
    if args.publish_md is not None:
        args.publish_md.parent.mkdir(parents=True, exist_ok=True)
        write_markdown(args.publish_md, result)
    print(json.dumps({"status": result["status"], "output": str(args.output_root.resolve())}))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
