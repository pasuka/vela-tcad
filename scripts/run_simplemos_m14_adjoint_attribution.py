#!/usr/bin/env python3
"""Run SimpleMOS M14 frozen-state and adjoint terminal-current attribution."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable


REPO = Path(__file__).resolve().parents[1]
M10 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m10_fixed_state_replay")
M8 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
      / "m8_original_physics/vela")
M8_DEEP = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
           / "m8_deep_off_diagnostics/vela")
M13 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m13_spatial_attribution")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m14_adjoint_attribution")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "adjoint_attribution")
SPECTRUM = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "terminal_sensitivity/m12_error_spectrum_points.csv")
M13_SUMMARY = (REPO / "reference_tcad/simplemos_sentaurus2022"
               / "spatial_attribution/m13_state_summary.csv")
RUNNER = REPO / "build-release/vela_example_runner.exe"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def drain_tag(value: float) -> str:
    return "0p05" if value == 0.05 else "1"


def gate_tag(value: float) -> str:
    return {0.0: "0", 0.05: "0p05", 0.8: "0p8", 2.5: "2p5"}[value]


def state_tag(device: str, drain: float, gate: float) -> str:
    return f"{device}_vd_{drain_tag(drain)}_vg_{gate_tag(gate)}"


def vela_state(device: str, drain: float, gate: float) -> Path:
    branch = M8 / device / "workflow" / f"vd_{drain_tag(drain)}"
    if gate == 0.0:
        return branch / "10_drain_ramp_accepted_state.csv"
    if gate == 0.05:
        return (M8_DEEP / f"{device}_vd_{drain_tag(drain)}"
                / "vg_0p05/accepted_state.csv")
    if gate == 0.8:
        return (M13 / "self_consistent"
                / f"{device}_vd_{drain_tag(drain)}_vg_0p8"
                / "accepted_state.csv")
    if gate == 2.5:
        return branch / "20_gate_sweep_accepted_state.csv"
    raise ValueError(gate)


def run_config(config: Path, runner: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run([str(runner), "--config", str(config)],
                               cwd=REPO, text=True, capture_output=True,
                               env=env, check=False)
    config.with_suffix(".stdout.txt").write_text(completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"M14 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def write_hybrid(path: Path, base_rows: list[dict[str, str]],
                 alt_rows: list[dict[str, str]], fields: tuple[str, ...],
                 mask: Callable[[int], bool]) -> None:
    columns = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
    alt = {int(row["node_id"]): row for row in alt_rows}
    rows = []
    for source in base_rows:
        node = int(source["node_id"])
        row = {column: source[column] for column in columns}
        if mask(node):
            for field in fields:
                row[field] = alt[node][field]
        rows.append(row)
    write_csv(path, rows)


def vector_from_residual(path: Path) -> list[float]:
    rows = read_csv(path)
    values: list[float] = []
    for column in ("psi_residual", "phin_residual", "phip_residual"):
        values.extend(float(row[column]) for row in rows)
    return values


def dot(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise ValueError("vector size mismatch")
    return sum(a * b for a, b in zip(left, right))


def linear_fit(rows: list[dict[str, Any]]) -> dict[str, float]:
    xs = [float(row["gate_voltage_V"]) for row in rows]
    ys = [float(row["signed_log10_ratio_dex"]) for row in rows]
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    denominator = sum((x - xbar) ** 2 for x in xs)
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denominator
    intercept = ybar - slope * xbar
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    return {"slope_dex_per_V": slope, "intercept_dex": intercept,
            "rmse_dex": math.sqrt(sum(r * r for r in residuals) / len(residuals))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable", type=Path, default=PORTABLE)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    portable = args.portable.resolve()
    output.mkdir(parents=True, exist_ok=True)
    portable.mkdir(parents=True, exist_ok=True)

    m13_by_state = {row["state"]: row for row in read_csv(M13_SUMMARY)}
    cases: list[dict[str, Any]] = []
    jobs: list[tuple[Path, str, Path]] = []
    for device in ("n17", "n21"):
        for drain in (0.05, 1.0):
            for gate in (0.0, 0.05, 0.8, 2.5):
                tag = state_tag(device, drain, gate)
                case_dir = output / tag
                case_dir.mkdir(parents=True, exist_ok=True)
                source_config = M10 / "replay" / tag / "vela_drive_edge_mobility.json"
                config = read_json(source_config)
                base_state = vela_state(device, drain, gate)
                sent_state = M10 / "replay" / tag / "sentaurus_state_for_vela.csv"
                for path in (base_state, sent_state, Path(config["mesh_file"])):
                    if not path.is_file():
                        raise FileNotFoundError(path)

                adjoint_csv = case_dir / "adjoint.csv"
                adjoint_config = case_dir / "adjoint.json"
                adjoint_deck = dict(config)
                adjoint_deck.update({
                    "simulation_type": "terminal_current_adjoint_probe",
                    "state_file": str(base_state.resolve()),
                    "output_csv": str(adjoint_csv.resolve()),
                    "contact": "drain",
                    "m14": {"read_only": True, "role": "baseline_adjoint"},
                })
                write_json(adjoint_config, adjoint_deck)
                jobs.append((adjoint_config, "adjoint", adjoint_csv))

                mesh = read_json(Path(config["mesh_file"]))
                coordinates = {int(node["id"]): (float(node["x"]), float(node["y"]))
                               for node in mesh["nodes"]}
                base_rows = read_csv(base_state)
                sent_rows = read_csv(sent_state)
                variants: list[tuple[str, tuple[str, ...], Callable[[int], bool]]] = [
                    ("field_psi", ("psi",), lambda _node: True),
                    ("field_phin", ("phin",), lambda _node: True),
                    ("field_phip", ("phip",), lambda _node: True),
                    ("field_all", ("psi", "phin", "phip"), lambda _node: True),
                    ("region_body", ("psi", "phin", "phip"),
                     lambda node, c=coordinates: c[node][0] > 0.25),
                    ("region_source", ("psi", "phin", "phip"),
                     lambda node, c=coordinates: c[node][0] <= 0.25 and c[node][1] < -0.3),
                    ("region_drain", ("psi", "phin", "phip"),
                     lambda node, c=coordinates: c[node][0] <= 0.25 and c[node][1] > 0.3),
                    ("region_channel", ("psi", "phin", "phip"),
                     lambda node, c=coordinates: c[node][0] <= 0.25 and -0.3 <= c[node][1] <= 0.3),
                ]
                all_variants = [("baseline", tuple(), lambda _node: False)] + variants
                variant_metadata = []
                for name, fields, mask in all_variants:
                    state_path = case_dir / f"{name}_state.csv"
                    write_hybrid(state_path, base_rows, sent_rows, fields, mask)
                    residual_csv = case_dir / f"{name}_residual.csv"
                    probe_config = case_dir / f"{name}_functional.json"
                    deck = dict(config)
                    deck.update({
                        "simulation_type": "terminal_current_functional_probe",
                        "state_file": str(state_path.resolve()),
                        "contact": "drain",
                        "residual_output_csv": str(residual_csv.resolve()),
                        "m14": {"read_only": True, "variant": name},
                    })
                    write_json(probe_config, deck)
                    jobs.append((probe_config, name, residual_csv))
                    variant_metadata.append({"name": name, "state": state_path,
                                             "config": probe_config,
                                             "residual": residual_csv})
                cases.append({"state": tag, "device": device, "drain": drain,
                              "gate": gate, "dir": case_dir,
                              "adjoint_config": adjoint_config,
                              "adjoint_csv": adjoint_csv,
                              "base_state": base_state, "sent_state": sent_state,
                              "variants": variant_metadata,
                              "m13": m13_by_state[tag]})

    selected = [(config, role, expected) for config, role, expected in jobs
                if args.force or not expected.exists() or not config.with_suffix(".stdout.txt").exists()]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, args.runner.resolve()):
                   (config, role) for config, role, _expected in selected}
        for future in as_completed(futures):
            config, role = futures[future]
            future.result()
            print(f"completed {config.parent.name}: {role}", flush=True)

    substitution_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    for case in cases:
        adjoint_status = read_json(case["adjoint_config"].with_suffix(".stdout.txt"))
        adjoint_data = read_csv(case["adjoint_csv"])
        gradient: list[float] = []
        adjoint: list[float] = []
        for column in ("dI_dpsi_scaled", "dI_dphin_scaled", "dI_dphip_scaled"):
            gradient.extend(float(row[column]) for row in adjoint_data)
        for column in ("lambda_poisson", "lambda_electron", "lambda_hole"):
            adjoint.extend(float(row[column]) for row in adjoint_data)
        base_state = {int(row["node_id"]): row for row in read_csv(case["base_state"])}
        sent_state = {int(row["node_id"]): row for row in read_csv(case["sent_state"])}
        potential_scale = float(adjoint_status["potential_scale_V"])
        block_delta: dict[str, list[float]] = {}
        for field in ("psi", "phin", "phip"):
            block_delta[field] = [
                (float(sent_state[node][field]) - float(base_state[node][field]))
                / potential_scale for node in sorted(base_state)]
        full_delta = block_delta["psi"] + block_delta["phin"] + block_delta["phip"]

        by_variant: dict[str, dict[str, Any]] = {}
        for variant in case["variants"]:
            status = read_json(variant["config"].with_suffix(".stdout.txt"))
            residual = vector_from_residual(variant["residual"])
            by_variant[variant["name"]] = {"status": status, "residual": residual}
        baseline = by_variant["baseline"]
        base_current = float(baseline["status"]["current_A_per_um"])
        if not math.isclose(base_current, float(adjoint_status["current_A_per_um"]),
                            rel_tol=1e-11, abs_tol=1e-30):
            raise AssertionError(f"adjoint/current mismatch for {case['state']}")

        for variant in case["variants"]:
            name = variant["name"]
            if name == "baseline":
                continue
            trial = by_variant[name]
            delta_residual = [a - b for a, b in zip(trial["residual"], baseline["residual"])]
            if name == "field_psi":
                delta = block_delta["psi"] + [0.0] * (2 * len(block_delta["psi"]))
            elif name == "field_phin":
                delta = ([0.0] * len(block_delta["psi"]) + block_delta["phin"]
                         + [0.0] * len(block_delta["psi"]))
            elif name == "field_phip":
                delta = [0.0] * (2 * len(block_delta["psi"])) + block_delta["phip"]
            else:
                hybrid_rows = {int(row["node_id"]): row
                               for row in read_csv(variant["state"])}
                delta = []
                for field in ("psi", "phin", "phip"):
                    delta.extend((float(hybrid_rows[node][field]) -
                                  float(base_state[node][field])) / potential_scale
                                 for node in sorted(base_state))
            direct_linear = dot(gradient, delta)
            direct_exact = float(trial["status"]["current_A_per_um"]) - base_current
            relaxation = -dot(adjoint, delta_residual)
            substitution_rows.append({
                "state": case["state"], "device": case["device"],
                "drain_voltage_V": case["drain"], "gate_voltage_V": case["gate"],
                "variant": name, "baseline_current_A_per_um": base_current,
                "direct_linear_current_change_A_per_um": direct_linear,
                "direct_exact_current_change_A_per_um": direct_exact,
                "adjoint_relaxation_current_change_A_per_um": relaxation,
                "two_layer_sum_A_per_um": direct_exact + relaxation,
                "direct_linear_relative_to_current": direct_linear / max(abs(base_current), 1e-300),
                "direct_exact_relative_to_current": direct_exact / max(abs(base_current), 1e-300),
                "adjoint_relaxation_relative_to_current": relaxation / max(abs(base_current), 1e-300),
            })
        all_row = next(row for row in substitution_rows
                       if row["state"] == case["state"] and row["variant"] == "field_all")
        state_rows.append({
            "state": case["state"], "device": case["device"],
            "drain_voltage_V": case["drain"], "gate_voltage_V": case["gate"],
            "terminal_signed_log10_error_dex": case["m13"]["terminal_signed_log10_error_dex"],
            "baseline_current_A_per_um": base_current,
            "adjoint_relative_residual": adjoint_status["adjoint_relative_residual"],
            "all_field_direct_exact_relative_to_current": all_row["direct_exact_relative_to_current"],
            "all_field_adjoint_relaxation_relative_to_current": all_row["adjoint_relaxation_relative_to_current"],
            "all_field_two_layer_relative_to_current":
                float(all_row["two_layer_sum_A_per_um"]) / max(abs(base_current), 1e-300),
        })

    write_csv(portable / "m14_state_summary.csv", state_rows)
    write_csv(portable / "m14_substitution_summary.csv", substitution_rows)

    spectrum_rows = [row for row in read_csv(SPECTRUM)
                     if row["device"] in {"n21", "n22", "n23", "n24"}
                     and math.isclose(float(row["drain_voltage_V"]), 0.05)
                     and float(row["gate_voltage_V"]) <= 0.95 + 1e-12]
    if len(spectrum_rows) != 80:
        raise AssertionError(f"expected 80 weak-region points, got {len(spectrum_rows)}")
    write_csv(portable / "m14_weak_region_80_points.csv", spectrum_rows)
    fits = {device: linear_fit([row for row in spectrum_rows if row["device"] == device])
            for device in ("n21", "n22", "n23", "n24")}
    abs_errors = [abs(float(row["signed_log10_ratio_dex"])) for row in spectrum_rows]
    low_drain_rows = [row for row in substitution_rows
                      if float(row["drain_voltage_V"]) == 0.05]
    field_means = {}
    for variant in ("field_psi", "field_phin", "field_phip", "field_all"):
        selected_rows = [row for row in low_drain_rows if row["variant"] == variant]
        field_means[variant] = {
            "mean_abs_direct_exact_relative_to_current":
                sum(abs(float(row["direct_exact_relative_to_current"])) for row in selected_rows)
                / len(selected_rows),
            "mean_abs_adjoint_relaxation_relative_to_current":
                sum(abs(float(row["adjoint_relaxation_relative_to_current"])) for row in selected_rows)
                / len(selected_rows),
        }
    region_means = {}
    for variant in ("region_source", "region_channel", "region_drain", "region_body"):
        selected_rows = [row for row in low_drain_rows if row["variant"] == variant]
        region_means[variant] = {
            "mean_abs_direct_exact_relative_to_current":
                sum(abs(float(row["direct_exact_relative_to_current"])) for row in selected_rows)
                / len(selected_rows),
            "mean_abs_adjoint_relaxation_relative_to_current":
                sum(abs(float(row["adjoint_relaxation_relative_to_current"])) for row in selected_rows)
                / len(selected_rows),
        }
    report = {
        "schema": "vela.simplemos.sdevice.m14_adjoint_attribution.report.v1",
        "status": "complete_with_declared_field_coverage",
        "execution": {
            "weak_region_terminal_point_count": len(spectrum_rows),
            "paired_spatial_state_count": len(cases),
            "adjoint_solve_count": len(cases),
            "substitution_variant_count_per_state": 8,
            "new_sentaurus_execution": False,
            "default_model_changed": False,
            "cpp_changed": True,
        },
        "coverage_boundary": {
            "terminal_error_spectrum": "complete 80/80",
            "paired_nodal_field_attribution": "16 representative states only",
            "missing_for_80_point_spatial_extension":
                "Sentaurus nodal field exports for n21-n24 at every Vg=0..0.95 point",
        },
        "adjoint_quality": {
            "maximum_relative_residual": max(float(row["adjoint_relative_residual"])
                                             for row in state_rows),
        },
        "weak_region": {
            "all_vela_above_sentaurus": all(float(row["signed_log10_ratio_dex"]) > 0
                                             for row in spectrum_rows),
            "mean_signed_error_dex": sum(float(row["signed_log10_ratio_dex"])
                                          for row in spectrum_rows) / len(spectrum_rows),
            "maximum_absolute_error_dex": max(abs_errors),
            "linear_fits_by_device": fits,
        },
        "low_drain_representative_attribution": {
            "field_means": field_means,
            "region_means": region_means,
        },
        "interpretation_limits": [
            "The adjoint is local to the converged Vela state and is a first-order self-consistent sensitivity.",
            "Field and region substitutions can be large; exact frozen current changes are reported separately from the adjoint response.",
            "Carrier densities are derived from psi and quasi-Fermi potentials and are not substituted as independent Newton unknowns.",
            "Mobility and SRH evidence is linked from M10-M12 rather than represented as independent state-vector coordinates.",
        ],
    }
    write_json(portable / "m14_adjoint_attribution_report.json", report)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
