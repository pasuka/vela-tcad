#!/usr/bin/env python3
"""Regression coverage for the A2/B1.05 spatial validation runner."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_vanoverstraeten_A2_B105_spatial_validation.py"


class VanOverstraetenA2B105SpatialValidationTest(unittest.TestCase):
    def test_skip_run_writes_spatial_validation_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vo_a2_b105_spatial_") as td:
            tmp = Path(td)
            out_dir = tmp / "diagnostics"
            base_config = tmp / "base_bv.json"
            base_config.write_text(
                json.dumps({
                    "mesh_file": "mesh.json",
                    "node_doping_file": "doping.csv",
                    "materials_file": "materials.json",
                    "solver": {"impact_ionization": {"model": "van_overstraeten"}},
                    "sweep": {
                        "mode": "bv_reverse",
                        "contact": "Anode",
                        "current_contact": "Anode",
                        "bias_points": [0, -5],
                        "write_vtk": True,
                        "vtk_prefix": "vtk/dc_sweep",
                    },
                }),
                encoding="utf-8",
            )
            reference_csv = tmp / "sentaurus_current.csv"
            write_rows(reference_csv, ["bias_V", "current_total"], [
                {"bias_V": "0", "current_total": "0"},
                {"bias_V": "-5", "current_total": "1e-6"},
            ])
            sentaurus_dir = tmp / "sentaurus_multibias"
            write_sentaurus_case(sentaurus_dir / "sentaurus_0v", 1.0)
            write_sentaurus_case(sentaurus_dir / "sentaurus_-5v", 2.0)

            case_dir = out_dir / "vanoverstraeten_A2_B105_spatial_validation_case" / "BV-A2-B1p05-spatial"
            (case_dir / "vtk").mkdir(parents=True)
            write_rows(
                case_dir / "BV-A2-B1p05-spatial.csv",
                [
                    "bias_V", "current_electron", "current_hole", "current_total",
                    "current_electron_A_per_um", "current_hole_A_per_um",
                    "current_total_A_per_um", "converged", "iterations",
                    "newton_iterations", "step_diagnostics",
                ],
                [
                    make_curve_row("0", 1.0),
                    make_curve_row("-5", 1.1),
                ],
            )
            write_vtk(case_dir / "vtk" / "dc_sweep_0000_0V.vtk", 1.0)
            write_vtk(case_dir / "vtk" / "dc_sweep_0001_-5V.vtk", 2.0)

            b_matrix_csv = tmp / "b_matrix.csv"
            write_rows(
                b_matrix_csv,
                [
                    "B_scale", "bias_V", "terminal_total_current_ratio_to_sentaurus",
                    "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus",
                    "q_integral_AvalancheGeneration_junction_window_ratio_to_sentaurus",
                    "q_integral_AvalancheGeneration_center_window_ratio_to_sentaurus",
                    "q_integral_AvalancheGeneration_left_shoulder_ratio_to_sentaurus",
                    "q_integral_AvalancheGeneration_right_shoulder_ratio_to_sentaurus",
                ],
                [
                    {"B_scale": "1.00", "bias_V": "-5.0",
                     "terminal_total_current_ratio_to_sentaurus": "0.9",
                     "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus": "1.5",
                     "q_integral_AvalancheGeneration_junction_window_ratio_to_sentaurus": "1.4",
                     "q_integral_AvalancheGeneration_center_window_ratio_to_sentaurus": "1.3",
                     "q_integral_AvalancheGeneration_left_shoulder_ratio_to_sentaurus": "1.2",
                     "q_integral_AvalancheGeneration_right_shoulder_ratio_to_sentaurus": "1.1"},
                ],
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--base-config", str(base_config),
                    "--out-dir", str(out_dir),
                    "--reference-current-csv", str(reference_csv),
                    "--sentaurus-multibias-dir", str(sentaurus_dir),
                    "--b-scale-matrix-csv", str(b_matrix_csv),
                    "--bias-points", "0,-5",
                    "--skip-run",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            out_csv = out_dir / "vanoverstraeten_A2_B105_spatial_validation.csv"
            summary_md = out_dir / "vanoverstraeten_A2_B105_spatial_validation_summary.md"
            self.assertTrue(out_csv.exists())
            self.assertTrue(summary_md.exists())
            with out_csv.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            row = rows[1]
            for column in [
                "A_scale",
                "B_scale",
                "qG_full_semiconductor_ratio_to_sentaurus",
                "qG_left_shoulder_ratio_to_sentaurus",
                "sentaurus_max_Gava_cm3_per_s",
                "max_Gava_distance_to_sentaurus_um",
                "A2_B100_qG_full_semiconductor_ratio_to_sentaurus",
            ]:
                self.assertIn(column, row)
            self.assertEqual(row["A_scale"], "2.0")
            self.assertEqual(row["B_scale"], "1.05")

            text = summary_md.read_text(encoding="utf-8")
            self.assertIn("negative Anode bias is reverse bias", text)
            self.assertIn("source mapping", text)
            self.assertIn("B_low_scale", text)
            self.assertIn("current density and continuity source feedback", text)


def make_curve_row(bias: str, scale: float) -> dict[str, str]:
    return {
        "bias_V": bias,
        "current_electron": str(scale * 0.7e-6),
        "current_hole": str(scale * 0.3e-6),
        "current_total": str(scale * 1.0e-6),
        "current_electron_A_per_um": str(scale * 0.7e-6),
        "current_hole_A_per_um": str(scale * 0.3e-6),
        "current_total_A_per_um": str(scale * 1.0e-6),
        "converged": "1",
        "iterations": "4",
        "newton_iterations": "5",
        "step_diagnostics": f"accepted_step={bias};retry_count=0",
    }


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_sentaurus_case(case_dir: Path, scale: float) -> None:
    (case_dir / "fields").mkdir(parents=True)
    write_rows(case_dir / "nodes.csv", ["id", "x_um", "y_um"], [
        {"id": "0", "x_um": "0.75", "y_um": "0"},
        {"id": "1", "x_um": "1.0", "y_um": "0"},
        {"id": "2", "x_um": "1.25", "y_um": "0.25"},
    ])
    write_rows(case_dir / "elements.csv", ["id", "node0", "node1", "node2"], [
        {"id": "0", "node0": "0", "node1": "1", "node2": "2"},
    ])
    for name, values in {
        "ImpactIonization": [scale * 1.0e14, scale * 2.0e14, scale * 1.0e14],
        "ElectricField": [scale * 1.0e5, scale * 2.0e5, scale * 1.0e5],
        "eAlphaAvalanche": [scale * 1.0e6, scale * 2.0e6, scale * 1.0e6],
        "hAlphaAvalanche": [scale * 0.8e6, scale * 1.6e6, scale * 0.8e6],
    }.items():
        write_rows(case_dir / "fields" / f"{name}_region0.csv", ["node_id", "component0"], [
            {"node_id": str(idx), "component0": str(value)}
            for idx, value in enumerate(values)
        ])


def write_vtk(path: Path, scale: float) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("test\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        handle.write("POINTS 3 double\n")
        handle.write("7.5e-7 0 0\n1.0e-6 0 0\n1.25e-6 2.5e-7 0\n")
        handle.write("CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n")
        handle.write("POINT_DATA 3\n")
        handle.write("VECTORS ElectricFieldVector double\n")
        handle.write(f"{scale * 1.0e5} 0 0\n{scale * 2.0e5} 0 0\n{scale * 1.0e5} 0 0\n")
        for name, values in {
            "ElectronHighFieldDrive": [scale * 2.0e5, scale * 2.1e5, scale * 1.9e5],
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
