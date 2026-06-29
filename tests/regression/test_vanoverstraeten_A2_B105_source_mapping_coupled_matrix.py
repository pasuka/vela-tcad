#!/usr/bin/env python3
"""Regression coverage for the A2/B1.05 coupled source mapping matrix runner."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_vanoverstraeten_A2_B105_source_mapping_coupled_matrix.py"


class VanOverstraetenA2B105SourceMappingCoupledMatrixTest(unittest.TestCase):
    def test_skip_run_collects_all_source_mapping_modes(self) -> None:
        source_case = (
            REPO / "build/diagnostics/vanoverstraeten_A2_B105_spatial_validation_case"
            / "BV-A2-B1p05-spatial"
        )
        if not source_case.exists():
            self.skipTest("A2_B105 spatial validation case is not available")
        with tempfile.TemporaryDirectory(prefix="vela_source_mapping_matrix_") as td:
            tmp = Path(td)
            out_dir = tmp / "diagnostics"
            base_config = tmp / "base.json"
            base_config.write_text(
                json.dumps({
                    "mesh_file": "mesh.json",
                    "node_doping_file": "doping.csv",
                    "materials_file": "materials.json",
                    "solver": {"impact_ionization": {"model": "van_overstraeten"}},
                    "sweep": {"bias_points": [0, -5], "write_vtk": True},
                }),
                encoding="utf-8",
            )
            case_root = out_dir / "vanoverstraeten_A2_B105_source_mapping_coupled_matrix_cases"
            for case_id in ["A2B105-map-node", "A2B105-map-edge", "A2B105-map-cell"]:
                target = case_root / case_id
                shutil.copytree(source_case, target)
                old_csv = target / "BV-A2-B1p05-spatial.csv"
                new_csv = target / f"{case_id}.csv"
                old_csv.rename(new_csv)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--base-config", str(base_config),
                    "--out-dir", str(out_dir),
                    "--bias-points", "0,-5",
                    "--skip-run",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            out_csv = out_dir / "vanoverstraeten_A2_B105_source_mapping_coupled_matrix.csv"
            summary = out_dir / "vanoverstraeten_A2_B105_source_mapping_coupled_matrix_summary.md"
            self.assertTrue(out_csv.exists())
            self.assertTrue(summary.exists())
            with out_csv.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 6)
            self.assertEqual(
                {row["source_mapping_mode"] for row in rows},
                {
                    "node_F_node_alpha_node_G",
                    "edge_F_edge_alpha_edge_G_to_node",
                    "cell_F_cell_alpha_cell_G_to_node",
                },
            )
            for column in [
                "terminal_total_current_ratio_to_sentaurus",
                "qG_full_semiconductor_ratio_to_sentaurus",
                "qG_left_shoulder_ratio_to_sentaurus",
                "max_Gava_distance_to_sentaurus_um",
                "node_reproduce_A2_B105_current_ratio_delta",
            ]:
                self.assertIn(column, rows[0])
            text = summary.read_text(encoding="utf-8")
            self.assertIn("Max Gava position", text)
            self.assertIn("Left/right shoulder", text)
            self.assertIn("Full/junction/center", text)
            self.assertIn("node_F_node_alpha_node_G reproduction", text)
            self.assertIn("sentaurus_fit_A_B", text)


if __name__ == "__main__":
    unittest.main()
