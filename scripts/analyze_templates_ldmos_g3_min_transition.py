#!/usr/bin/env python3
"""Align Templates/LDMOS G3 fixed-state transport and the first Newton step.

The report deliberately separates three contracts:

* immutable 0.0231559774221138 V state operator replay;
* Sentaurus contact/current fields versus Vela SG/contact-cut integration;
* the 0.0231559774221138 -> 0.0316416579413075 V first Newton step.

Sentaurus NewtonPlot exports may be supplied for both the history-aware path and
the Save/Load no-predictor control.  Generated probe decks and large CSV files
belong in an ignored reference_staging directory.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import subprocess
import statistics
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
SOURCE_BIAS_V = 0.0231559774221138
TARGET_BIAS_V = 0.0316416579413075


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def finite(values: Iterable[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def distribution(values: Iterable[float]) -> dict[str, float | int]:
    samples = sorted(finite(values))
    if not samples:
        return {"count": 0, "maximum": math.nan, "median": math.nan, "l2": math.nan}
    return {
        "count": len(samples),
        "maximum": samples[-1],
        "median": statistics.median(samples),
        "l2": math.sqrt(sum(value * value for value in samples)),
    }


def signed_ratio(candidate: float, reference: float) -> dict[str, float]:
    result = {
        "candidate": candidate,
        "reference": reference,
        "absolute_difference": candidate - reference,
    }
    if reference != 0.0:
        result["signed_ratio"] = candidate / reference
        result["magnitude_error_dex"] = abs(
            math.log10(max(abs(candidate), 1.0e-300))
            - math.log10(max(abs(reference), 1.0e-300))
        )
    return result


def source_case(bias: float) -> dict[str, Any]:
    return {
        "group": "templates_ldmos_g3",
        "bias_kind": "drain",
        "bias_V": bias,
        "gate_bias_V": 0.0,
        "drain_bias_V": bias,
    }


def probe_suite(
    base: Any,
    runner: Path,
    state_file: Path,
    output_dir: Path,
    sentaurus_export: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    fixed_case = source_case(SOURCE_BIAS_V)
    artifacts: dict[str, str] = {}
    statuses: dict[str, Any] = {}
    rows: dict[str, list[dict[str, str]]] = {}
    for probe in (
        "sg_edge_flux_probe",
        "edge_mobility_probe",
        "newton_residual_probe",
        "newton_carrier_term_probe",
    ):
        csv_path, status = base.run_probe(
            runner, fixed_case, output_dir, probe, state_file
        )
        label = probe.removesuffix("_probe")
        artifacts[label] = str(csv_path.resolve())
        statuses[label] = status
        rows[label] = read_csv(csv_path)

    step_csv, step_status = base.run_probe(
        runner,
        source_case(TARGET_BIAS_V),
        output_dir,
        "newton_step_probe",
        state_file,
    )
    artifacts["newton_step"] = str(step_csv.resolve())
    statuses["newton_step"] = step_status
    rows["newton_step"] = read_csv(step_csv)

    sg_rows = rows["sg_edge_flux"]
    mobility_rows = rows["edge_mobility"]
    edge_comparison = output_dir / "sentaurus_nodal_current_vs_vela_sg_edges.csv"
    edge_metrics = base.edge_transport_metrics(
        sentaurus_export, sg_rows, mobility_rows, edge_comparison
    )
    artifacts["edge_comparison"] = str(edge_comparison.resolve())

    return {
        "state_file": str(state_file.resolve()),
        "artifacts": artifacts,
        "statuses": statuses,
        "sentaurus_projection_vs_vela_sg": edge_metrics,
        "endpoint_density": base.endpoint_density_metrics(sentaurus_export, sg_rows),
        "contact_cut_currents_A_per_um": contact_cut_currents(sentaurus_export, sg_rows),
        "residual_extrema": residual_extrema(rows["newton_residual"]),
        "carrier_term_extrema": carrier_term_extrema(rows["newton_carrier_term"]),
        "first_step_extrema": step_extrema(rows["newton_step"]),
    }


def contact_nodes(export_dir: Path) -> dict[str, set[int]]:
    result: dict[str, set[int]] = {}
    for row in read_csv(export_dir / "contacts.csv"):
        result[row["name"]] = {
            int(value) for value in row["node_ids"].split(";") if value
        }
    return result


def contact_cut_currents(
    export_dir: Path, edge_rows: list[dict[str, str]]
) -> dict[str, dict[str, float | int]]:
    elementary_charge = 1.602176634e-19
    result: dict[str, dict[str, float | int]] = {}
    for name, nodes in contact_nodes(export_dir).items():
        electron_particles = 0.0
        hole_particles = 0.0
        count = 0
        for row in edge_rows:
            at0 = int(row["node0"]) in nodes
            at1 = int(row["node1"]) in nodes
            if at0 == at1:
                continue
            sign = 1.0 if at0 else -1.0
            electron_particles += sign * float(
                row["electron_particle_line_flux_per_m_s"]
            )
            hole_particles += sign * float(row["hole_particle_line_flux_per_m_s"])
            count += 1
        electron = -elementary_charge * electron_particles * 1.0e-6
        hole = elementary_charge * hole_particles * 1.0e-6
        result[name] = {
            "crossing_edge_count": count,
            "electron_A_per_um": electron,
            "hole_A_per_um": hole,
            "total_A_per_um": electron + hole,
        }
    return result


def contact_scalars(export_dir: Path, field_name: str) -> dict[str, float]:
    manifest = json.loads(
        (export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    result: dict[str, float] = {}
    for field in manifest["fields"]:
        if field["name"] != field_name:
            continue
        rows = read_csv(export_dir / "fields" / field["csv_file"])
        if rows:
            result[field["region_name"]] = float(rows[0]["component0"])
    return result


def residual_extrema(rows: list[dict[str, str]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for block, column in (
        ("poisson", "psi_residual"),
        ("electron", "phin_residual"),
        ("hole", "phip_residual"),
    ):
        ranked = sorted(rows, key=lambda row: abs(float(row[column])), reverse=True)
        top = ranked[0]
        result[block] = {
            "absolute": distribution(abs(float(row[column])) for row in rows),
            "signed_sum": sum(float(row[column]) for row in rows),
            "top": {
                "node_id": int(top["node_id"]),
                # Probe node coordinates preserve the unit-scaling mesh's um
                # coordinates; unlike SG edge rows they are not exported in m.
                "x_um": float(top["x"]),
                "y_um": float(top["y"]),
                "value": float(top[column]),
            },
        }
    return result


def carrier_term_extrema(rows: list[dict[str, str]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for carrier in ("electron", "hole"):
        columns = [
            f"{carrier}_flux",
            f"{carrier}_recombination",
            f"{carrier}_residual",
        ]
        result[carrier] = {
            column.removeprefix(f"{carrier}_"): {
                "absolute": distribution(abs(float(row[column])) for row in rows),
                "signed_sum": sum(float(row[column]) for row in rows),
            }
            for column in columns
        }
    return result


def step_extrema(rows: list[dict[str, str]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for block, delta, residual, trial in (
        ("poisson", "delta_psi_V", "psi_residual", "trial_psi_residual"),
        ("electron", "delta_phin_V", "phin_residual", "trial_phin_residual"),
        ("hole", "delta_phip_V", "phip_residual", "trial_phip_residual"),
    ):
        ranked = sorted(rows, key=lambda row: abs(float(row[delta])), reverse=True)
        top = ranked[0]
        result[block] = {
            "delta_absolute": distribution(abs(float(row[delta])) for row in rows),
            "residual_absolute": distribution(abs(float(row[residual])) for row in rows),
            "trial_residual_absolute": distribution(abs(float(row[trial])) for row in rows),
            "top_delta": {
                "node_id": int(top["node_id"]),
                "x_um": float(top["x"]),
                "y_um": float(top["y"]),
                "value_V": float(top[delta]),
            },
        }
    return result


def state_delta(
    left: Path, right: Path, node_ids: set[int] | None = None
) -> dict[str, Any]:
    left_rows = {int(row["node_id"]): row for row in read_csv(left)}
    right_rows = {int(row["node_id"]): row for row in read_csv(right)}
    common = sorted(left_rows.keys() & right_rows.keys())
    if node_ids is not None:
        common = [node for node in common if node in node_ids]
    result: dict[str, Any] = {"common_nodes": len(common)}
    for column in ("psi", "phin", "phip"):
        result[f"{column}_absolute_difference_V"] = distribution(
            abs(float(left_rows[node][column]) - float(right_rows[node][column]))
            for node in common
        )
    for column in ("electrons_m3", "holes_m3"):
        errors = []
        for node in common:
            left_value = max(float(left_rows[node][column]), 1.0)
            right_value = max(float(right_rows[node][column]), 1.0)
            errors.append(abs(math.log10(left_value) - math.log10(right_value)))
        result[f"{column}_absolute_difference_dex"] = distribution(errors)
    return result


def field(export_dir: Path, name: str, region: int = 0) -> dict[int, float]:
    path = export_dir / "fields" / f"{name}_region{region}.csv"
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(path)
    }


def sentaurus_newton_export(export_dir: Path) -> dict[str, Any]:
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(export_dir / "nodes.csv")
    }
    result: dict[str, Any] = {}
    for key, name, unit in (
        ("poisson_rhs", "PoissonRhs", "C"),
        ("electron_rhs", "eContinuityRhs", "A"),
        ("hole_rhs", "hContinuityRhs", "A"),
        ("psi_update", "NewtonStepElectrostaticPotentialUpdate", "V"),
        ("electron_density_update", "NewtonStepEDensityUpdate", "cm^-3"),
        ("hole_density_update", "NewtonStepHDensityUpdate", "cm^-3"),
    ):
        values = field(export_dir, name)
        top_node = max(values, key=lambda node: abs(values[node]))
        result[key] = {
            "unit": unit,
            "absolute": distribution(abs(value) for value in values.values()),
            "signed_sum": sum(values.values()),
            "top": {
                "node_id": top_node,
                "x_um": nodes[top_node][0],
                "y_um": nodes[top_node][1],
                "value": values[top_node],
            },
        }
    return result


def sentaurus_newton_state_delta(left: Path, right: Path) -> dict[str, Any]:
    left_psi = field(left, "ElectrostaticPotential")
    right_psi = field(right, "ElectrostaticPotential")
    common = sorted(left_psi.keys() & right_psi.keys())
    result: dict[str, Any] = {
        "common_silicon_nodes": len(common),
        "psi_absolute_difference_V": distribution(
            abs(left_psi[node] - right_psi[node]) for node in common
        ),
    }
    for label, name in (("electron", "eDensity"), ("hole", "hDensity")):
        left_density = field(left, name)
        right_density = field(right, name)
        result[f"{label}_density_absolute_difference_dex"] = distribution(
            abs(
                math.log10(max(abs(left_density[node]), 1.0))
                - math.log10(max(abs(right_density[node]), 1.0))
            )
            for node in common
        )
    return result


def constant_high_field_control(baseline: Path, output: Path) -> Path:
    """Materialize the closest current Vela proxy for HFS without DopingDep.

    Unit-scaling decks consume cm/um TCAD internal values even though legacy
    key spellings contain SI-looking suffixes.  Setting CT mu_min equal to the
    material mobility makes its doping factor identically constant, leaving
    only the high-field limiter active.
    """
    config = json.loads(baseline.read_text(encoding="utf-8"))
    config["solver"]["mobility"] = {
        "model": "caughey_thomas_field",
        "doping_concentration_basis": "net_doping",
        "high_field_driving_force": "quasi_fermi_gradient",
        "high_field_gradient_discretization": "edge_projection",
        "electron_mu_min_m2_V_s": 1417.0,
        "electron_nref_m3": 1.0e17,
        "electron_alpha": 1.0,
        "electron_saturation_velocity_m_s": 1.07e7,
        "electron_high_field_beta": 1.109,
        "hole_mu_min_m2_V_s": 470.5,
        "hole_nref_m3": 1.0e17,
        "hole_alpha": 1.0,
        "hole_saturation_velocity_m_s": 8.37e6,
        "hole_high_field_beta": 1.213,
    }
    output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return output


def redirect_csv_files(value: Any, output_dir: Path, prefix: str = "") -> None:
    """Keep a generated control's mutable artifacts in its ignored run dir."""
    if isinstance(value, dict):
        for key, child in value.items():
            label = f"{prefix}_{key}" if prefix else key
            if key in {"csv_file", "summary_file", "write_state_file"}:
                suffix = Path(str(child)).suffix or ".csv"
                value[key] = str((output_dir / f"{label}{suffix}").resolve())
            else:
                redirect_csv_files(child, output_dir, label)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            redirect_csv_files(child, output_dir, f"{prefix}_{index}")


def run_secant_sweep_control(
    base: Any, runner: Path, baseline: Path, output_dir: Path
) -> dict[str, Any]:
    """Run the unchanged G3 deck with only the existing secant predictor enabled."""
    output_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(baseline.read_text(encoding="utf-8"))
    config["_comment"] = (
        "G3 single-factor control: baseline physics with Vela secant sweep predictor."
    )
    config["output_csv"] = str((output_dir / "iv.csv").resolve())
    config["sweep"].setdefault("continuation", {})["predictor"] = {
        "mode": "secant",
        "fields": ["psi", "phin", "phip"],
        "max_extrapolation_ratio": 2.0,
    }
    redirect_csv_files(config["sweep"], output_dir, "sweep")
    config_path = output_dir / "g3_secant_predictor_control.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    log_path = (output_dir / "g3_secant_predictor_control.log").resolve()
    completed = subprocess.run(
        [
            str(runner),
            "--config",
            str(config_path.resolve()),
            "--log",
            str(log_path),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=base.runner_environment(),
        check=False,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    csv_path = Path(config["output_csv"])
    rows = read_csv(csv_path) if csv_path.exists() else []
    status: dict[str, Any] = {}
    if completed.stdout.strip():
        try:
            status = json.loads(completed.stdout.strip().splitlines()[-1])
        except json.JSONDecodeError:
            status = {"stdout_last_line": completed.stdout.strip().splitlines()[-1]}
    return {
        "return_code": completed.returncode,
        "status": status,
        "accepted_biases_V": [
            float(row["bias_V"]) for row in rows if row.get("converged") == "1"
        ],
        "attempted_biases_V": [float(row["bias_V"]) for row in rows],
        "curve_rows": len(rows),
        "config": str(config_path.resolve()),
        "output_csv": str(csv_path.resolve()),
        "log": str(log_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--sentaurus-export", type=Path, required=True)
    parser.add_argument("--sentaurus-state", type=Path, required=True)
    parser.add_argument("--vela-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-constant-high-field-control", action="store_true")
    parser.add_argument("--run-secant-sweep-control", action="store_true")
    parser.add_argument(
        "--sentaurus-newton-export",
        action="append",
        default=[],
        metavar="LABEL=DIR",
    )
    args = parser.parse_args()

    base = load_module(
        "templates_ldmos_formula_base",
        REPO / "scripts/run_transportmodels_sentaurus_formula_replay.py",
    )
    base.IDVG_CONFIG = args.baseline_config.resolve()
    # The reusable TransportModels helpers target that template's silicon
    # region 3.  Templates/LDMOS has its silicon material in region 0.
    transportmodels_field = base.field
    base.field = lambda export_dir, name, _region, components=1: transportmodels_field(
        export_dir, name, 0, components
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sentaurus = probe_suite(
        base,
        args.runner.resolve(),
        args.sentaurus_state.resolve(),
        args.output_dir / "sentaurus_state",
        args.sentaurus_export.resolve(),
    )
    vela = probe_suite(
        base,
        args.runner.resolve(),
        args.vela_state.resolve(),
        args.output_dir / "vela_state",
        args.sentaurus_export.resolve(),
    )
    controls: dict[str, Any] = {}
    if args.run_constant_high_field_control:
        control_config = constant_high_field_control(
            args.baseline_config.resolve(),
            args.output_dir / "constant_high_field_control_baseline.json",
        )
        base.IDVG_CONFIG = control_config.resolve()
        controls["sentaurus_state_constant_high_field"] = probe_suite(
            base,
            args.runner.resolve(),
            args.sentaurus_state.resolve(),
            args.output_dir / "sentaurus_state_constant_high_field",
            args.sentaurus_export.resolve(),
        )
    if args.run_secant_sweep_control:
        controls["vela_secant_sweep"] = run_secant_sweep_control(
            base,
            args.runner.resolve(),
            args.baseline_config.resolve(),
            args.output_dir / "vela_secant_sweep",
        )
    sentaurus_contact = contact_scalars(
        args.sentaurus_export, "ContactCurrentFlux"
    )
    terminal_alignment = {}
    for contact, reference in sentaurus_contact.items():
        if contact in sentaurus["contact_cut_currents_A_per_um"]:
            terminal_alignment[contact] = signed_ratio(
                float(
                    sentaurus["contact_cut_currents_A_per_um"][contact][
                        "total_A_per_um"
                    ]
                ),
                reference,
            )

    newton_exports: dict[str, Any] = {}
    newton_export_paths: dict[str, Path] = {}
    for encoded in args.sentaurus_newton_export:
        if "=" not in encoded:
            raise ValueError("--sentaurus-newton-export requires LABEL=DIR")
        label, path = encoded.split("=", 1)
        newton_export_paths[label] = Path(path)
        newton_exports[label] = sentaurus_newton_export(Path(path))

    predictor_state_effect: dict[str, Any] = {}
    if {"history_n0", "no_predictor_n0"} <= newton_export_paths.keys():
        predictor_state_effect["iteration_0"] = sentaurus_newton_state_delta(
            newton_export_paths["history_n0"],
            newton_export_paths["no_predictor_n0"],
        )
    if {"history_n1", "no_predictor_n1"} <= newton_export_paths.keys():
        predictor_state_effect["iteration_1"] = sentaurus_newton_state_delta(
            newton_export_paths["history_n1"],
            newton_export_paths["no_predictor_n1"],
        )

    report = {
        "schema": "vela.templates_ldmos_g3_min_transition_alignment.v1",
        "status": "complete",
        "biases_V": {"source": SOURCE_BIAS_V, "target": TARGET_BIAS_V},
        "contracts": {
            "fixed_state": "immutable state replay through production Vela operators",
            "sentaurus_edge_reference": "nodal current density projected onto Vela primal edges; diagnostic, not a discrete identity",
            "vela_terminal": "signed SG contact-cut integration per 1 um depth",
            "newton": "first assembled residual and production Newton step at target contact bias",
        },
        "sentaurus_tdr_contact_current_A": sentaurus_contact,
        "sentaurus_state_terminal_alignment": terminal_alignment,
        "sentaurus_state": sentaurus,
        "vela_last_accepted_state": vela,
        "sentaurus_vs_vela_state": {
            "all_nodes": state_delta(args.sentaurus_state, args.vela_state),
            "silicon_nodes": state_delta(
                args.sentaurus_state,
                args.vela_state,
                set(field(args.sentaurus_export.resolve(), "eDensity")),
            ),
        },
        "sentaurus_newton_exports": newton_exports,
        "sentaurus_predictor_state_effect": predictor_state_effect,
        "controls": controls,
    }
    summary = args.output_dir / "summary.json"
    summary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": str(summary.resolve())}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
