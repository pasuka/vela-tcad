#!/usr/bin/env python3
"""Regression coverage for the current-density export-path audit."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_current_density_export_path.py"


class CurrentDensityExportPathAuditTest(unittest.TestCase):
    def test_audit_reports_vector_sum_and_units(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_current_density_export_") as td:
            tmp = Path(td)
            vtk = tmp / "probe.vtk"
            out_md = tmp / "audit.md"
            out_csv = tmp / "consistency.csv"
            vtk.write_text(
                """# vtk DataFile Version 3.0
synthetic
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 double
0 0 0
1e-6 2e-6 0
POINT_DATA 2
VECTORS ElectronCurrentDensityVector double
1 2 0
-3 4 0
VECTORS HoleCurrentDensityVector double
5 7 0
11 -13 0
VECTORS TotalCurrentDensityVector double
6 9 0
8 -9 0
""",
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--vtk",
                    str(vtk),
                    "--out-md",
                    str(out_md),
                    "--out-csv",
                    str(out_csv),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(out_md.exists())
            self.assertTrue(out_csv.exists())

            with out_csv.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["total_equals_electron_plus_hole"], "yes")
            self.assertAlmostEqual(float(rows[0]["relative_error_mag"]), 0.0)
            self.assertEqual(rows[1]["total_equals_electron_plus_hole"], "yes")

            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("data source: node recovery / least-squares gradient post-processing", summary)
            self.assertIn("output unit: A/cm^2", summary)
            self.assertIn("TotalCurrentDensity = ElectronCurrentDensity + HoleCurrentDensity: yes", summary)
            self.assertIn("A/m^2 -> A/cm^2: explicit /1e4 conversion", summary)


if __name__ == "__main__":
    unittest.main()
