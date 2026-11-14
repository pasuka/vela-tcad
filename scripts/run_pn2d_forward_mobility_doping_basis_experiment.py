#!/usr/bin/env python3
"""Run and summarize PN2D mobility doping-basis candidates."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

from generate_pn2d_config import render_named_template


CANDIDATES = (
    "net_doping",
    "total_impurity",
    "cell_reconstructed_total_impurity",
)
ANCHORS = (1, 2, 5, 10, 15, 20)


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def read_vtk_scalar(path: Path, name: str) -> list[float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    marker = f"SCALARS {name} "
    start = next(i for i, line in enumerate(lines) if line.startswith(marker))
    values: list[float] = []
    for line in lines[start + 2 :]:
        if line.startswith(("SCALARS ", "VECTORS ", "FIELD ")):
            break
        if line.strip():
            values.append(float(line))
    return values


def rms(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values) / len(values))


def sentaurus_scalar(fields_root: Path, bias: int, name: str) -> list[float]:
    values = {
        int(row["node_id"]): float(row["component0"])
        for row in csv_rows(
            fields_root / f"{bias}v" / "fields" / f"{name}_region0.csv"
        )
    }
    return [values[node] for node in sorted(values)]


def sentaurus_current(fields_root: Path, bias: int) -> float:
    path = (
        fields_root
        / f"{bias}v"
        / "fields"
        / "ContactCurrentFlux_region2.csv"
    )
    return abs(float(csv_rows(path)[0]["component0"]))


def run_candidate(
    runner: Path,
    baseline: dict,
    output_root: Path,
    basis: str,
) -> tuple[Path, dict]:
    candidate_dir = output_root / basis
    candidate_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(json.dumps(baseline))
    config["solver"]["mobility"]["doping_concentration_basis"] = basis
    prefix = candidate_dir / f"vela_{basis}"
    config["output_csv"] = str(prefix.with_suffix(".csv").resolve())
    config["sweep"]["vtk_prefix"] = str(prefix.resolve())
    config["sweep"]["write_state_file"] = str(
        (candidate_dir / f"vela_{basis}_last_state.csv").resolve()
    )
    config["sweep"]["diagnostics"]["newton_history"]["csv_file"] = str(
        (candidate_dir / f"vela_{basis}_newton_history.csv").resolve()
    )
    config_path = candidate_dir / f"simulation_{basis}.json"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    (candidate_dir / "runner.stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (candidate_dir / "runner.stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(
            f"{basis} failed with exit code {completed.returncode}: "
            f"{completed.stderr[-1000:]}"
        )
    return candidate_dir, config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--baseline-config", type=Path)
    source.add_argument(
        "--template",
        action="store_true",
        help="Use the versioned pn2d_iv production template.",
    )
    parser.add_argument("--mesh-file", type=Path)
    parser.add_argument("--node-doping-file", type=Path)
    parser.add_argument("--materials-file", type=Path)
    parser.add_argument("--sentaurus-fields", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.template:
        missing = [
            name for name, value in (
                ("--mesh-file", args.mesh_file),
                ("--node-doping-file", args.node_doping_file),
                ("--materials-file", args.materials_file),
            ) if value is None
        ]
        if missing:
            parser.error(f"--template requires {', '.join(missing)}")
        baseline, _ = render_named_template(
            "pn2d_iv",
            {
                "solver_profile": "legacy",
                "mesh_file": str(args.mesh_file.resolve()),
                "node_doping_file": str(args.node_doping_file.resolve()),
                "materials_file": str(args.materials_file.resolve()),
            },
            allow_absolute_paths=True,
        )
    else:
        baseline = json.loads(
            args.baseline_config.read_text(encoding="utf-8-sig")
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    selected_rows = []
    node_rows = []

    sent_e = sentaurus_scalar(args.sentaurus_fields, 20, "eMobility")
    sent_h = sentaurus_scalar(args.sentaurus_fields, 20, "hMobility")
    sent_psi = sentaurus_scalar(
        args.sentaurus_fields, 20, "ElectrostaticPotential"
    )
    sent_phin = sentaurus_scalar(
        args.sentaurus_fields, 20, "eQuasiFermiPotential"
    )
    sent_phip = sentaurus_scalar(
        args.sentaurus_fields, 20, "hQuasiFermiPotential"
    )

    for basis in CANDIDATES:
        candidate_dir, _ = run_candidate(
            args.runner.resolve(), baseline, args.out_dir.resolve(), basis
        )
        iv_path = candidate_dir / f"vela_{basis}.csv"
        iv = {
            round(float(row["bias_V"]), 9): row
            for row in csv_rows(iv_path)
            if row["converged"] == "1"
        }
        vtk = candidate_dir / f"vela_{basis}_0200_20V.vtk"
        e_mu = read_vtk_scalar(vtk, "ElectronMobility")
        h_mu = read_vtk_scalar(vtk, "HoleMobility")
        psi = read_vtk_scalar(vtk, "Potential")
        phin = read_vtk_scalar(vtk, "ElectronQuasiFermi")
        phip = read_vtk_scalar(vtk, "HoleQuasiFermi")
        junction = range(7, 16)
        e_rmse = rms([e_mu[i] - sent_e[i] for i in junction])
        h_rmse = rms([h_mu[i] - sent_h[i] for i in junction])
        psi_rmse = rms([a - b for a, b in zip(psi, sent_psi)])
        common_rmse = rms(
            [
                ((psi[i] - sent_psi[i]) + (phin[i] - sent_phin[i]) +
                 (phip[i] - sent_phip[i])) / 3.0
                for i in range(len(psi))
            ]
        )
        anchor_errors = []
        for bias in ANCHORS:
            sent_current = sentaurus_current(args.sentaurus_fields, bias)
            vela_current = abs(float(iv[float(bias)]["current_total_A_per_um"]))
            rel = (vela_current - sent_current) / sent_current
            anchor_errors.append(abs(rel))
            selected_rows.append(
                {
                    "basis": basis,
                    "bias_V": bias,
                    "sentaurus_A_per_um": sent_current,
                    "vela_A_per_um": vela_current,
                    "relative_difference_percent": 100.0 * rel,
                }
            )
        for node in junction:
            node_rows.append(
                {
                    "basis": basis,
                    "node_id": node,
                    "sentaurus_e_mobility_cm2_V_s": sent_e[node],
                    "vela_e_mobility_cm2_V_s": e_mu[node],
                    "sentaurus_h_mobility_cm2_V_s": sent_h[node],
                    "vela_h_mobility_cm2_V_s": h_mu[node],
                }
            )
        summaries.append(
            {
                "basis": basis,
                "converged_points": len(iv),
                "junction_nodes_7_15_e_mobility_rmse_cm2_V_s": e_rmse,
                "junction_nodes_7_15_h_mobility_rmse_cm2_V_s": h_rmse,
                "exact_anchor_median_abs_current_error_percent": (
                    100.0 * sorted(anchor_errors)[len(anchor_errors) // 2]
                ),
                "current_error_20V_percent": selected_rows[-1][
                    "relative_difference_percent"
                ],
                "psi_error_20V_rms_V": psi_rmse,
                "common_mode_error_20V_rms_V": common_rmse,
            }
        )

    def write_csv(path: Path, values: list[dict]) -> None:
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=values[0].keys())
            writer.writeheader()
            writer.writerows(values)

    write_csv(args.out_dir / "candidate_summary.csv", summaries)
    write_csv(args.out_dir / "exact_anchor_current_comparison.csv", selected_rows)
    write_csv(args.out_dir / "junction_mobility_comparison_20V.csv", node_rows)
    (args.out_dir / "candidate_summary.json").write_text(
        json.dumps({"candidates": summaries}, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
