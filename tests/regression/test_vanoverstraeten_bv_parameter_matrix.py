#!/usr/bin/env python3
"""Regression coverage for the VanOverstraeten BV parameter matrix runner."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_vanoverstraeten_bv_parameter_matrix.py"


class VanOverstraetenBvParameterMatrixTest(unittest.TestCase):
    def test_skip_run_summarizes_existing_case_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vo_bv_matrix_") as td:
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
            ])

            sentaurus_dir = tmp / "sentaurus_multibias"
            sent_0v = sentaurus_dir / "sentaurus_0v"
            (sent_0v / "fields").mkdir(parents=True)
            write_rows(sent_0v / "nodes.csv", ["id", "x_um", "y_um"], [
                {"id": "0", "x_um": "0.75", "y_um": "0"},
                {"id": "1", "x_um": "1.0", "y_um": "0"},
                {"id": "2", "x_um": "1.25", "y_um": "0"},
            ])
            write_rows(sent_0v / "elements.csv", ["id", "node0", "node1", "node2"], [
                {"id": "0", "node0": "0", "node1": "1", "node2": "2"},
            ])
            write_rows(sent_0v / "fields" / "ImpactIonization_region0.csv",
                       ["node_id", "component0"], [
                {"node_id": "0", "component0": "1e10"},
                {"node_id": "1", "component0": "2e10"},
                {"node_id": "2", "component0": "1e10"},
            ])

            cases = [
                ("BV-P0-default", "default", 1.0),
                ("BV-P1-fit-A-only", "sentaurus_fit_A_only", 2.0),
                ("BV-P2-fit-A-B", "sentaurus_fit_A_B", 4.0),
                ("BV-P3-fit-A-B-switch", "sentaurus_fit_A_B_switch", 3.0),
            ]
            for case_id, _, scale in cases:
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
                        "max_electric_field_V_per_cm",
                    ],
                    [{
                        "mode": "bv_reverse",
                        "bias_contact": "anode",
                        "bias_V": "0",
                        "current_contact": "anode",
                        "current_electron": str(scale),
                        "current_hole": str(scale * 0.5),
                        "current_total": str(scale * 1.5),
                        "converged": "1",
                        "iterations": "2",
                        "newton_iterations": "3",
                        "step_diagnostics": "attempted_step=0;accepted_step=0;retry_count=0",
                        "current_total_A_per_um": str(scale * 1e-6),
                        "current_electron_A_per_um": str(scale * 7e-7),
                        "current_hole_A_per_um": str(scale * 3e-7),
                        "max_electric_field_V_per_cm": str(scale * 1e5),
                    }],
                )
                write_vtk(case_dir / "vtk" / "dc_sweep_0000_0V.vtk", scale)

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
                    "0",
                    "--skip-run",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            by_bias = out_dir / "vanoverstraeten_bv_parameter_matrix_by_bias.csv"
            summary = out_dir / "vanoverstraeten_bv_parameter_matrix_summary.csv"
            md = out_dir / "vanoverstraeten_bv_parameter_matrix.md"
            self.assertTrue(by_bias.exists())
            self.assertTrue(summary.exists())
            self.assertTrue(md.exists())

            with by_bias.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 4)
            first = rows[0]
            for column in [
                "case_id",
                "parameter_set",
                "bias_V",
                "terminal_total_current_A_per_um",
                "max_Fn_V_per_cm",
                "max_eAlphaAvalanche_cm_inv",
                "integral_AvalancheGeneration_full_semiconductor_per_s_per_um",
                "q_integral_AvalancheGeneration_full_semiconductor_A_per_um",
                "voltage_step_used_V",
            ]:
                self.assertIn(column, first)

            text = md.read_text(encoding="utf-8")
            self.assertIn("bias sign convention", text)
            self.assertIn("sentaurus_fit_A_B", text)
            self.assertIn("source mapping", text)


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_vtk(path: Path, scale: float) -> None:
    values = {
        "ElectricFieldVector": [(scale * 1.0e5, 0.0, 0.0), (scale * 1.2e5, 0.0, 0.0), (scale * 0.9e5, 0.0, 0.0)],
        "ElectronHighFieldDrive": [scale * 2.0e5, scale * 2.2e5, scale * 1.8e5],
        "HoleHighFieldDrive": [scale * 1.5e5, scale * 1.7e5, scale * 1.4e5],
        "ElectronAlphaAvalanche": [scale * 1.0e6, scale * 1.2e6, scale * 0.8e6],
        "HoleAlphaAvalanche": [scale * 0.9e6, scale * 1.1e6, scale * 0.7e6],
        "AvalancheGeneration": [scale * 1.0e20, scale * 2.0e20, scale * 1.0e20],
    }
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("test\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        handle.write("POINTS 3 double\n")
        handle.write("7.5e-7 0 0\n1.0e-6 0 0\n1.25e-6 0 0\n")
        handle.write("CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n")
        handle.write("POINT_DATA 3\n")
        for name, field_values in values.items():
            if isinstance(field_values[0], tuple):
                handle.write(f"VECTORS {name} double\n")
                for x, y, z in field_values:
                    handle.write(f"{x} {y} {z}\n")
            else:
                handle.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
                for value in field_values:
                    handle.write(f"{value}\n")


if __name__ == "__main__":
    unittest.main()
