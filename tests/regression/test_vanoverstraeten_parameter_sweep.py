#!/usr/bin/env python3
"""Regression coverage for the VanOverstraeten parameter-set sweep diagnostic."""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
BUILD_DIR = Path(os.environ.get("VELA_BUILD_DIR", REPO / "build"))
EXE = BUILD_DIR / (
    "vanoverstraeten_parameter_sweep.exe"
    if sys.platform == "win32"
    else "vanoverstraeten_parameter_sweep"
)


class VanOverstraetenParameterSweepTest(unittest.TestCase):
    def test_tool_compares_all_parameter_sets_by_region_and_branch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vo_parameter_sweep_") as td:
            tmp = Path(td)
            compare_csv = tmp / "coarse_node_field_compare_aligned.csv"
            elements_csv = tmp / "elements.csv"
            out_csv = tmp / "vanoverstraeten_parameter_sweep.csv"
            summary_md = tmp / "vanoverstraeten_parameter_sweep_summary.md"

            elements_csv.write_text(
                "\n".join([
                    "id,node0,node1,node2,region,material",
                    "0,0,1,2,R.Si,Si",
                    "1,1,3,2,R.Si,Si",
                    "2,1,4,3,R.Si,Si",
                    "3,4,5,3,R.Si,Si",
                ]),
                encoding="utf-8",
            )

            fieldnames = [
                "bias_V",
                "node_id",
                "x_um",
                "y_um",
                "quantity",
                "sentaurus_value",
                "vela_value_scaled_to_sentaurus_units",
            ]
            nodes = [
                (0, 0.75, 0.0),
                (1, 1.00, 0.0),
                (2, 0.75, 0.5),
                (3, 1.00, 0.5),
                (4, 1.25, 0.0),
                (5, 1.25, 0.5),
            ]
            rows: list[dict[str, str]] = []
            field_v_per_cm = 3.0e5
            alpha_e = 6.78391642452e7 * math.exp(-1.21718982697e6 / field_v_per_cm)
            alpha_h = 1.41230834668e8 * math.exp(-1.99067614831e6 / field_v_per_cm)
            for node_id, x_um, y_um in nodes:
                e_qf = 30.0 * x_um
                h_qf = -30.0 * x_um
                for quantity, sentaurus in [
                    ("electron_qf", e_qf),
                    ("hole_qf", h_qf),
                    ("electron_alpha_avalanche", alpha_e),
                    ("hole_alpha_avalanche", alpha_h),
                ]:
                    rows.append({
                        "bias_V": "-20",
                        "node_id": str(node_id),
                        "x_um": str(x_um),
                        "y_um": str(y_um),
                        "quantity": quantity,
                        "sentaurus_value": f"{sentaurus:.15g}",
                        "vela_value_scaled_to_sentaurus_units": "0",
                    })

            with compare_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)

            result = subprocess.run(
                [
                    str(EXE),
                    "--input-csv",
                    str(compare_csv),
                    "--elements-csv",
                    str(elements_csv),
                    "--output-csv",
                    str(out_csv),
                    "--summary-md",
                    str(summary_md),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with out_csv.open(encoding="utf-8") as handle:
                out_rows = list(csv.DictReader(handle))
            self.assertEqual(len(out_rows), 12)

            first = out_rows[0]
            for column in [
                "region_label",
                "carrier",
                "alpha_default_cm_inv",
                "alpha_fit_A_only_cm_inv",
                "alpha_fit_A_B_cm_inv",
                "alpha_fit_A_B_switch_cm_inv",
                "ratio_default",
                "ratio_fit_A_only",
                "ratio_fit_A_B",
                "ratio_fit_A_B_switch",
                "branch_default",
                "branch_fit_A_B_switch",
            ]:
                self.assertIn(column, first)

            electron_center = next(
                row for row in out_rows
                if row["carrier"] == "electron" and row["region_label"] == "center"
            )
            self.assertEqual(electron_center["branch_default"], "low")
            self.assertEqual(electron_center["branch_fit_A_B_switch"], "high")
            self.assertAlmostEqual(
                float(electron_center["ratio_fit_A_B_switch"]), 1.0, delta=1.0e-9)

            regions = {row["region_label"] for row in out_rows}
            self.assertTrue({"center", "left_shoulder", "right_shoulder"}.issubset(regions))

            summary = summary_md.read_text(encoding="utf-8")
            self.assertIn("center", summary)
            self.assertIn("left_shoulder", summary)
            self.assertIn("sentaurus_fit_A_B_switch", summary)
            self.assertIn("2.5e5~4e5 V/cm", summary)


if __name__ == "__main__":
    unittest.main()
