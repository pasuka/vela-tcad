#!/usr/bin/env python3
"""Audit Templates/LDMOS Fermi and OldSlotboom state mapping at fixed bias."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable


BIAS_V = 0.0231559774221138
HEAVY_DOPING_M3 = 1.0e25
REPRESENTATIVE_NODES = (173, 4573, 4670)
THERMAL_VOLTAGE_300K = 0.025851999786435


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def distribution(values: Iterable[float]) -> dict[str, float | int]:
    samples = sorted(value for value in values if math.isfinite(value))
    if not samples:
        return {"count": 0, "median": math.nan, "maximum": math.nan,
                "l2": math.nan}
    return {
        "count": len(samples),
        "median": statistics.median(samples),
        "maximum": samples[-1],
        "l2": math.sqrt(sum(value * value for value in samples)),
    }


def scalar_field(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(path)
    }


def density_error(row: dict[str, str], carrier: str) -> float:
    supplied = max(float(row[f"input_{carrier}_density_m3"]), 1.0)
    reconstructed = max(
        float(row[f"reconstructed_{carrier}_density_m3"]), 1.0)
    return abs(math.log10(supplied) - math.log10(reconstructed))


def summarize_variant(
    rows: list[dict[str, str]],
    silicon_nodes: set[int],
) -> dict[str, Any]:
    silicon = [row for row in rows if int(row["node_id"]) in silicon_nodes]
    free = [row for row in silicon if row["constrained"] == "0"]
    heavy_n = [
        row for row in silicon
        if float(row["net_doping_m3"]) >= HEAVY_DOPING_M3
    ]
    heavy_p = [
        row for row in silicon
        if float(row["net_doping_m3"]) <= -HEAVY_DOPING_M3
    ]
    heavy_majority_density_errors = (
        [density_error(row, "electron") for row in heavy_n] +
        [density_error(row, "hole") for row in heavy_p]
    )
    heavy_majority_qf_errors = (
        [abs(float(row["electron_qf_mapping_error_V"])) for row in heavy_n] +
        [abs(float(row["hole_qf_mapping_error_V"])) for row in heavy_p]
    )
    return {
        "silicon_nodes": len(silicon),
        "free_silicon_nodes": len(free),
        "electron_density_abs_error_dex": distribution(
            density_error(row, "electron") for row in silicon),
        "hole_density_abs_error_dex": distribution(
            density_error(row, "hole") for row in silicon),
        "electron_qf_mapping_abs_error_V": distribution(
            abs(float(row["electron_qf_mapping_error_V"])) for row in silicon),
        "hole_qf_mapping_abs_error_V": distribution(
            abs(float(row["hole_qf_mapping_error_V"])) for row in silicon),
        "heavy_majority_density_abs_error_dex": distribution(
            heavy_majority_density_errors),
        "heavy_majority_qf_mapping_abs_error_V": distribution(
            heavy_majority_qf_errors),
        "electron_qf_mapping_signed_error_V": distribution(
            float(row["electron_qf_mapping_error_V"]) for row in silicon),
        "hole_qf_mapping_signed_error_V": distribution(
            float(row["hole_qf_mapping_error_V"]) for row in silicon),
        "free_poisson_residual_abs": distribution(
            abs(float(row["production_residual"])) for row in free),
        "representative_nodes": {
            str(node): next(
                (row for row in rows if int(row["node_id"]) == node), None)
            for node in REPRESENTATIVE_NODES
        },
    }


def write_reference_aligned_state(
    source: Path, destination: Path, old_slotboom_deg0_eV: float,
) -> dict[str, float]:
    """Apply a diagnostic half-dEg0 QF shift; this is not a production transform."""
    rows = read_csv(source)
    if not rows:
        raise ValueError("Sentaurus state is empty")
    electron_offset = -0.5 * old_slotboom_deg0_eV
    hole_offset = 0.5 * old_slotboom_deg0_eV
    fieldnames = list(rows[0])
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            adjusted = dict(row)
            adjusted["phin"] = format(
                float(row["phin"]) + electron_offset, ".17g")
            adjusted["phip"] = format(
                float(row["phip"]) + hole_offset, ".17g")
            writer.writerow(adjusted)
    return {
        "old_slotboom_dEg0_eV": old_slotboom_deg0_eV,
        "electron_qf_offset_V": electron_offset,
        "hole_qf_offset_V": hole_offset,
    }


def variants(base_solver: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    production = deepcopy(base_solver)
    no_correction = deepcopy(production)
    no_correction["bandgap_narrowing"]["fermi_statistics_correction"] = False

    correction_only = deepcopy(production)
    correction_only["bandgap_narrowing"]["coefficient_eV"] = 0.0
    correction_only["bandgap_narrowing"]["offset_eV"] = 0.0

    fermi_no_bgn = deepcopy(production)
    fermi_no_bgn["bandgap_narrowing"] = {"model": "none"}

    boltzmann_oldslotboom = deepcopy(no_correction)
    boltzmann_oldslotboom["carrier_statistics"] = {"model": "boltzmann"}

    boltzmann_no_bgn = deepcopy(fermi_no_bgn)
    boltzmann_no_bgn["carrier_statistics"] = {"model": "boltzmann"}

    return [
        ("fermi_oldslotboom_correction", production),
        ("fermi_oldslotboom_no_correction", no_correction),
        ("fermi_correction_only", correction_only),
        ("fermi_no_bgn", fermi_no_bgn),
        ("boltzmann_oldslotboom", boltzmann_oldslotboom),
        ("boltzmann_no_bgn", boltzmann_no_bgn),
    ]


def probe_config(
    baseline: dict[str, Any], solver: dict[str, Any], state: Path,
    output: Path,
) -> dict[str, Any]:
    config = deepcopy(baseline)
    config["simulation_type"] = "newton_poisson_term_probe"
    config["solver"] = solver
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output.resolve())
    config.pop("sweep", None)
    config.pop("output_vtk", None)
    for contact in config.get("contacts", []):
        if contact["name"] == "drain":
            contact["bias"] = BIAS_V
        elif contact["name"] in {"gate", "source", "substrate"}:
            contact["bias"] = 0.0
    return config


def run_probe(runner: Path, config: Path, log: Path) -> dict[str, Any]:
    environment = os.environ.copy()
    if os.name == "nt":
        prefixes = [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin"]
        environment["PATH"] = os.pathsep.join(prefixes + [environment.get("PATH", "")])
    completed = subprocess.run(
        [str(runner), "--config", str(config), "--log", str(log)],
        text=True, capture_output=True, env=environment, check=False)
    config.with_suffix(".stdout.txt").write_text(
        completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(
        completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)
    return json.loads(completed.stdout.strip().splitlines()[-1])


def add_bgn_field_comparison(
    summaries: dict[str, Any], rows_by_variant: dict[str, list[dict[str, str]]],
    sentaurus_bgn: dict[int, float], silicon_nodes: set[int],
) -> None:
    base_ni = {
        int(row["node_id"]): float(row["ni_eff_m3"])
        for row in rows_by_variant["fermi_no_bgn"]
    }
    for name, rows in rows_by_variant.items():
        errors: list[float] = []
        inferred: dict[str, Any] = {}
        for row in rows:
            node = int(row["node_id"])
            ni = float(row["ni_eff_m3"])
            ni0 = base_ni.get(node, 0.0)
            if node in silicon_nodes and node in sentaurus_bgn and ni > 0.0 and ni0 > 0.0:
                apparent_bgn = 2.0 * THERMAL_VOLTAGE_300K * math.log(ni / ni0)
                errors.append(abs(apparent_bgn - sentaurus_bgn[node]))
                if node in REPRESENTATIVE_NODES:
                    inferred[str(node)] = {
                        "vela_apparent_bgn_eV": apparent_bgn,
                        "sentaurus_bandgap_narrowing_eV": sentaurus_bgn[node],
                        "absolute_error_eV": abs(apparent_bgn - sentaurus_bgn[node]),
                    }
        summaries[name]["sentaurus_bgn_abs_error_eV"] = distribution(errors)
        summaries[name]["representative_bgn"] = inferred


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--sentaurus-state", type=Path, required=True)
    parser.add_argument("--sentaurus-export", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--sentaurus-old-slotboom-deg0-eV", type=float, required=True,
        help="Sentaurus dEg0 used to test the band-centre QF gauge transform")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    silicon_nodes = set(scalar_field(
        args.sentaurus_export / "fields" / "DopingConcentration_region0.csv"))
    sentaurus_bgn = scalar_field(
        args.sentaurus_export / "fields" / "BandgapNarrowing_region0.csv")

    summaries: dict[str, Any] = {}
    rows_by_variant: dict[str, list[dict[str, str]]] = {}
    artifacts: dict[str, Any] = {}
    for name, solver in variants(baseline["solver"]):
        csv_path = output_dir / f"{name}.csv"
        config_path = output_dir / f"{name}.json"
        config_path.write_text(json.dumps(probe_config(
            baseline, solver, args.sentaurus_state, csv_path), indent=2) + "\n",
            encoding="utf-8")
        status = run_probe(
            args.runner.resolve(), config_path.resolve(),
            (output_dir / f"{name}.log").resolve())
        rows = read_csv(csv_path)
        rows_by_variant[name] = rows
        summaries[name] = summarize_variant(rows, silicon_nodes)
        artifacts[name] = {
            "config": str(config_path), "csv": str(csv_path), "status": status}

    add_bgn_field_comparison(
        summaries, rows_by_variant, sentaurus_bgn, silicon_nodes)

    aligned_state = output_dir / "sentaurus_state_reference_aligned.csv"
    reference_transform = write_reference_aligned_state(
        args.sentaurus_state, aligned_state,
        args.sentaurus_old_slotboom_deg0_eV)
    aligned_name = "fermi_oldslotboom_correction_reference_aligned"
    aligned_csv = output_dir / f"{aligned_name}.csv"
    aligned_config = output_dir / f"{aligned_name}.json"
    aligned_config.write_text(json.dumps(probe_config(
        baseline, baseline["solver"], aligned_state, aligned_csv), indent=2) + "\n",
        encoding="utf-8")
    aligned_status = run_probe(
        args.runner.resolve(), aligned_config.resolve(),
        (output_dir / f"{aligned_name}.log").resolve())
    aligned_rows = read_csv(aligned_csv)
    reference_control = summarize_variant(aligned_rows, silicon_nodes)
    artifacts[aligned_name] = {
        "config": str(aligned_config), "csv": str(aligned_csv),
        "state": str(aligned_state), "status": aligned_status}
    ranking = sorted(summaries, key=lambda name: (
        summaries[name]["heavy_majority_qf_mapping_abs_error_V"]["median"],
        summaries[name]["heavy_majority_density_abs_error_dex"]["median"],
    ))
    result = {
        "schema": "vela.templates_ldmos.fermi_bgn_mapping_audit.v1",
        "status": "complete",
        "bias_V": BIAS_V,
        "state_file": str(args.sentaurus_state.resolve()),
        "silicon_nodes": len(silicon_nodes),
        "heavy_doping_threshold_m3": HEAVY_DOPING_M3,
        "ranking": ranking,
        "variants": summaries,
        "reference_energy_control": {
            "transform": reference_transform,
            "metrics": reference_control,
        },
        "artifacts": artifacts,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": str(summary_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
