#!/usr/bin/env python3
"""Regression coverage for the VanOverstraeten A=2 B-scale full BV matrix runner."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_vanoverstraeten_A2_B_scale_full_bv_matrix.py"


class VanOverstraetenA2BScaleFullBvMatrixTest(unittest.TestCase):
    def test_skip_run_summarizes_existing_b_scale_case_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vo_a2_b_scale_full_bv_") as td:
            tmp = Path(td)
            out_dir = tmp / "matrix"
            base_config = tmp / "base_bv.json"
            base_config.write_text(
                json.dumps({
                    "simulation_type": "dc_sweep",
                    "mesh_file": "mesh.json",
                    "node_doping_file": "doping.csv",
                    "output_csv": "base.csv",
                    "solver": {
                        "impact_ionization": {
                            "model": "van_overstraeten",
                            "parameter_set": "default",
                        },
                    },
                    "sweep": {
                        "mode": "bv_reverse",
                        "contact": "anode",
                        "current_contact": "anode",
                        "start": 0.0,
                        "stop": -20.0,
                        "step": -5.0,
                        "write_vtk": True,
                        "vtk_prefix": "vtk/dc_sweep",
                    },
                }, indent=2),
                encoding="utf-8",
            )

            reference_csv = tmp / "sentaurus_current.csv"
            write_rows(reference_csv, ["bias_V", "current_total"], [
                {"bias_V": "0", "current_total": "0"},
                {"bias_V": "-5", "current_total": "3e-6"},
            ])

            sentaurus_dir = tmp / "sentaurus_multibias"
            for bias_token, generation in [("0", 1.0e10), ("-5", 2.0e20)]:
                sent_dir = sentaurus_dir / f"sentaurus_{bias_token}v"
                (sent_dir / "fields").mkdir(parents=True)
                write_rows(sent_dir / "nodes.csv", ["id", "x_um", "y_um"], [
                    {"id": "0", "x_um": "0.75", "y_um": "0"},
                    {"id": "1", "x_um": "1.0", "y_um": "0"},
                    {"id": "2", "x_um": "1.25", "y_um": "0"},
                ])
                write_rows(sent_dir / "elements.csv", ["id", "node0", "node1", "node2"], [
                    {"id": "0", "node0": "0", "node1": "1", "node2": "2"},
                ])
                write_rows(
                    sent_dir / "fields" / "ImpactIonization_region0.csv",
                    ["node_id", "component0"],
                    [
                        {"node_id": "0", "component0": str(generation)},
                        {"node_id": "1", "component0": str(2.0 * generation)},
                        {"node_id": "2", "component0": str(generation)},
                    ],
                )

            cases = [
                ("BV-A2-B0p85", "0.85", 2.5, "1"),
                ("BV-A2-B0p90", "0.90", 2.0, "1"),
                ("BV-A2-B0p95", "0.95", 1.5, "1"),
                ("BV-A2-B1p00", "1.00", 1.0, "1"),
                ("BV-A2-B1p05", "1.05", 0.8, "1"),
                ("BV-A2-B1p10", "1.10", 0.6, "0"),
            ]
            for case_id, _, scale, converged in cases:
                case_dir = out_dir / "cases" / case_id
                (case_dir / "vtk").mkdir(parents=True)
                write_rows(
                    case_dir / f"{case_id}.csv",
                    [
                        "mode",
                        "bias_contact",
                        "bias_V",
                        "current_contact",
                        "current_electron",
                        "current_hole",
                        "current_total",
                        "converged",
                        "iterations",
                        "newton_iterations",
                        "step_diagnostics",
                        "current_total_A_per_um",
                        "current_electron_A_per_um",
                        "current_hole_A_per_um",
                    ],
                    [
                        make_curve_row("0", scale, "1"),
                        make_curve_row("-5", scale, converged),
                    ],
                )
                write_vtk(case_dir / "vtk" / "dc_sweep_0000_0V.vtk", scale)
                write_vtk(case_dir / "vtk" / "dc_sweep_0001_-5V.vtk", scale)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--base-config",
                    str(base_config),
                    "--out-dir",
                    str(out_dir),
                    "--reference-current-csv",
                    str(reference_csv),
                    "--sentaurus-multibias-dir",
                    str(sentaurus_dir),
                    "--bias-points",
                    "0,-5",
                    "--skip-run",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            out_csv = out_dir / "vanoverstraeten_A2_B_scale_full_bv_matrix.csv"
            summary_md = out_dir / "vanoverstraeten_A2_B_scale_full_bv_matrix_summary.md"
            self.assertTrue(out_csv.exists())
            self.assertTrue(summary_md.exists())

            with out_csv.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 12)
            first = rows[0]
            for column in [
                "case_id",
                "A_scale",
                "B_scale",
                "parameter_set",
                "bias_sign_convention",
                "convergence_status",
                "terminal_total_current_ratio_to_sentaurus",
                "actual_voltage_step_V",
                "max_Fn_V_per_cm",
                "max_eAlphaAvalanche_cm_inv",
                "q_integral_AvalancheGeneration_full_semiconductor_A_per_um",
                "q_integral_AvalancheGeneration_junction_window_ratio_to_sentaurus",
            ]:
                self.assertIn(column, first)
            self.assertTrue(all(row["A_scale"] == "2" for row in rows))
            self.assertIn("1.00", {row["B_scale"] for row in rows})

            text = summary_md.read_text(encoding="utf-8")
            self.assertIn("negative Anode bias is reverse bias", text)
            self.assertIn("B_scale<1", text)
            self.assertIn("B_low_scale", text)
            self.assertIn("sentaurus_fit_A_B", text)


def make_curve_row(bias: str, scale: float, converged: str) -> dict[str, str]:
    return {
        "mode": "bv_reverse",
        "bias_contact": "anode",
        "bias_V": bias,
        "current_contact": "anode",
        "current_electron": str(scale),
        "current_hole": str(scale * 0.5),
        "current_total": str(scale * 1.5),
        "converged": converged,
        "iterations": "2",
        "newton_iterations": "3",
        "step_diagnostics": f"attempted_step={bias};accepted_step={bias};retry_count=0",
        "current_total_A_per_um": str(scale * 1.0e-6),
        "current_electron_A_per_um": str(scale * 7.0e-7),
        "current_hole_A_per_um": str(scale * 3.0e-7),
    }


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_vtk(path: Path, scale: float) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("test\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        handle.write("POINTS 3 double\n")
        handle.write("7.5e-7 0 0\n1.0e-6 0 0\n1.25e-6 0 0\n")
        handle.write("CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n")
        handle.write("POINT_DATA 3\n")
        handle.write("VECTORS ElectricFieldVector double\n")
        handle.write(f"{scale * 1.0e5} 0 0\n{scale * 1.2e5} 0 0\n{scale * 0.9e5} 0 0\n")
        for name, values in {
            "ElectronHighFieldDrive": [scale * 2.0e5, scale * 2.2e5, scale * 1.8e5],
            "HoleHighFieldDrive": [scale * 1.5e5, scale * 1.7e5, scale * 1.4e5],
            "ElectronAlphaAvalanche": [scale * 1.0e6, scale * 1.2e6, scale * 0.8e6],
            "HoleAlphaAvalanche": [scale * 0.9e6, scale * 1.1e6, scale * 0.7e6],
            "AvalancheGeneration": [scale * 1.0e20, scale * 2.0e20, scale * 1.0e20],
        }.items():
            handle.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            for value in values:
                handle.write(f"{value}\n")


if __name__ == "__main__":
    unittest.main()