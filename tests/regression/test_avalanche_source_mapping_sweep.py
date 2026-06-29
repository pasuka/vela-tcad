#!/usr/bin/env python3
"""Regression coverage for the avalanche source mapping sweep diagnostic."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_avalanche_source_mapping_sweep.py"


class AvalancheSourceMappingSweepTest(unittest.TestCase):
    def test_skip_run_outputs_mapping_modes_and_summary_answers(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_avalanche_mapping_") as td:
            tmp = Path(td)
            out_dir = tmp / "diagnostics"
            case_dir = out_dir / "case"
            (case_dir / "vtk").mkdir(parents=True)
            write_rows(
                case_dir / "sweep.csv",
                [
                    "bias_V", "current_total_A_per_um", "converged",
                    "iterations", "newton_iterations", "step_diagnostics",
                ],
                [
                    {"bias_V": "0", "current_total_A_per_um": "0", "converged": "1",
                     "iterations": "1", "newton_iterations": "1", "step_diagnostics": "accepted_step=0"},
                    {"bias_V": "-5", "current_total_A_per_um": "1.0e-6", "converged": "1",
                     "iterations": "2", "newton_iterations": "3", "step_diagnostics": "accepted_step=-5"},
                ],
            )
            write_vtk(case_dir / "vtk" / "dc_sweep_0000_0V.vtk", scale=1.0)
            write_vtk(case_dir / "vtk" / "dc_sweep_0001_-5V.vtk", scale=2.0)

            sentaurus_dir = tmp / "sentaurus_multibias"
            write_sentaurus_case(sentaurus_dir / "sentaurus_0v", scale=1.0)
            write_sentaurus_case(sentaurus_dir / "sentaurus_-5v", scale=2.0)
            reference_current = tmp / "sentaurus_current.csv"
            write_rows(reference_current, ["bias_V", "current_total"], [
                {"bias_V": "0", "current_total": "0"},
                {"bias_V": "-5", "current_total": "1.0e-6"},
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--case-dir", str(case_dir),
                    "--sweep-csv", str(case_dir / "sweep.csv"),
                    "--sentaurus-multibias-dir", str(sentaurus_dir),
                    "--reference-current-csv", str(reference_current),
                    "--out-dir", str(out_dir),
                    "--bias-points", "0,-5",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            out_csv = out_dir / "avalanche_source_mapping_sweep.csv"
            summary_md = out_dir / "avalanche_source_mapping_sweep_summary.md"
            self.assertTrue(out_csv.exists())
            self.assertTrue(summary_md.exists())
            with out_csv.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            modes = {row["mapping_mode"] for row in rows}
            self.assertEqual(modes, {
                "node_F_node_alpha_node_G",
                "edge_F_edge_alpha_edge_G_to_node",
                "cell_F_cell_alpha_cell_G_to_node",
                "cell_F_cell_alpha_cell_G_integral_only",
            })
            self.assertEqual(len(rows), 8)
            for column in [
                "bias_V", "mapping_mode", "full_qG_ratio", "junction_qG_ratio",
                "center_qG_ratio", "left_shoulder_qG_ratio", "right_shoulder_qG_ratio",
                "max_Gava_self", "max_Gava_sentaurus", "max_Gava_distance_um",
                "terminal_current_ratio", "convergence_status",
            ]:
                self.assertIn(column, rows[0])
            self.assertTrue(any(row["full_qG_ratio"] for row in rows))

            text = summary_md.read_text(encoding="utf-8")
            self.assertIn("max Gava position", text)
            self.assertIn("left/right shoulder", text)
            self.assertIn("cell integral only", text)
            self.assertIn("source mapping", text)
            self.assertIn("B_low/B_high", text)
            self.assertIn("diagnostic option", text)


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
        handle.write("POINTS 4 double\n")
        handle.write("7.5e-7 0 0\n1.0e-6 0 0\n1.25e-6 2.5e-7 0\n1.0e-6 5.0e-7 0\n")
        handle.write("CELLS 2 8\n3 0 1 3\n3 1 2 3\nCELL_TYPES 2\n5\n5\n")
        handle.write("POINT_DATA 4\n")
        for name, values in {
            "ElectronImpactIonizationDrive": [scale * 1.0e7, scale * 1.2e7, scale * 1.1e7, scale * 0.9e7],
            "HoleImpactIonizationDrive": [scale * 0.8e7, scale * 0.9e7, scale * 1.0e7, scale * 0.7e7],
            "ElectronAlphaAvalanche": [scale * 1.0e6, scale * 2.0e6, scale * 1.5e6, scale * 0.5e6],
            "HoleAlphaAvalanche": [scale * 0.8e6, scale * 1.0e6, scale * 1.2e6, scale * 0.4e6],
            "AvalancheGeneration": [scale * 1.0e20, scale * 2.0e20, scale * 1.5e20, scale * 0.5e20],
        }.items():
            handle.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            for value in values:
                handle.write(f"{value}\n")
        for name, values in {
            "ElectronCurrentDensityVector": [(1.0, 0, 0), (2.0, 0, 0), (1.5, 0, 0), (0.5, 0, 0)],
            "HoleCurrentDensityVector": [(0.8, 0, 0), (1.0, 0, 0), (1.2, 0, 0), (0.4, 0, 0)],
            "ElectricFieldVector": [(scale * 1.0e5, 0, 0)] * 4,
        }.items():
            handle.write(f"VECTORS {name} double\n")
            for x, y, z in values:
                handle.write(f"{scale * x} {scale * y} {scale * z}\n")


def write_sentaurus_case(case_dir: Path, scale: float) -> None:
    (case_dir / "fields").mkdir(parents=True)
    write_rows(case_dir / "nodes.csv", ["id", "x_um", "y_um"], [
        {"id": "0", "x_um": "0.75", "y_um": "0"},
        {"id": "1", "x_um": "1.0", "y_um": "0"},
        {"id": "2", "x_um": "1.25", "y_um": "0.25"},
        {"id": "3", "x_um": "1.0", "y_um": "0.5"},
    ])
    write_rows(case_dir / "elements.csv", ["id", "node0", "node1", "node2"], [
        {"id": "0", "node0": "0", "node1": "1", "node2": "3"},
        {"id": "1", "node0": "1", "node1": "2", "node2": "3"},
    ])
    write_rows(case_dir / "fields" / "ImpactIonization_region0.csv", ["node_id", "component0"], [
        {"node_id": "0", "component0": str(scale * 1.0e14)},
        {"node_id": "1", "component0": str(scale * 2.0e14)},
        {"node_id": "2", "component0": str(scale * 1.5e14)},
        {"node_id": "3", "component0": str(scale * 0.5e14)},
    ])


if __name__ == "__main__":
    unittest.main()
