#!/usr/bin/env python3
"""Regression coverage for the VanOverstraeten A-scale IIC sweep."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_vanoverstraeten_A_scale_iic_sweep.py"


class VanOverstraetenAScaleIicSweepTest(unittest.TestCase):
    def test_a_scale_sweep_integrates_fixed_solution_and_marks_low_signal(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vo_a_scale_iic_") as td:
            tmp = Path(td)
            compare_csv = tmp / "coarse_node_field_compare_aligned.csv"
            elements_csv = tmp / "elements.csv"
            out_csv = tmp / "vanoverstraeten_A_scale_iic_sweep.csv"
            summary_md = tmp / "vanoverstraeten_A_scale_iic_sweep_summary.md"

            write_rows(
                elements_csv,
                ["id", "node0", "node1", "node2", "region", "material"],
                [
                    {"id": "0", "node0": "0", "node1": "1", "node2": "2", "region": "R.Si", "material": "Si"},
                    {"id": "1", "node0": "1", "node1": "3", "node2": "2", "region": "R.Si", "material": "Si"},
                ],
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
                (1, 1.25, 0.0),
                (2, 0.75, 0.5),
                (3, 1.25, 0.5),
            ]
            rows: list[dict[str, str]] = []
            for bias, sentaurus_gen in [("-5", "2.0e18"), ("-1", "1.0e-40")]:
                for node_id, x_um, y_um in nodes:
                    e_qf = 30.0 * x_um
                    h_qf = -30.0 * x_um
                    for quantity, sentaurus, vela in [
                        ("electron_qf", e_qf, e_qf),
                        ("hole_qf", h_qf, h_qf),
                        ("electron_current_total", "0", "2.0"),
                        ("hole_current_total", "0", "3.0"),
                        ("avalanche", sentaurus_gen, "0"),
                    ]:
                        rows.append({
                            "bias_V": bias,
                            "node_id": str(node_id),
                            "x_um": str(x_um),
                            "y_um": str(y_um),
                            "quantity": quantity,
                            "sentaurus_value": str(sentaurus),
                            "vela_value_scaled_to_sentaurus_units": str(vela),
                        })
            write_rows(compare_csv, fieldnames, rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input-csv",
                    str(compare_csv),
                    "--elements-csv",
                    str(elements_csv),
                    "--output-csv",
                    str(out_csv),
                    "--summary-md",
                    str(summary_md),
                    "--scales",
                    "1,2",
                    "--sentaurus-qg-floor",
                    "1e-30",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with out_csv.open(encoding="utf-8") as handle:
                out_rows = list(csv.DictReader(handle))
            self.assertEqual(len(out_rows), 4)

            first = out_rows[0]
            for column in [
                "bias_V",
                "A_scale",
                "qG_full_self_A_per_um",
                "qG_full_sentaurus_A_per_um",
                "qG_full_ratio",
                "qG_junction_ratio",
                "qG_center_ratio",
                "qG_left_shoulder_ratio",
                "qG_right_shoulder_ratio",
                "max_G_self",
                "max_G_sentaurus",
                "best_region_error_metric",
                "signal_status",
            ]:
                self.assertIn(column, first)

            rows_by_key = {(row["bias_V"], row["A_scale"]): row for row in out_rows}
            scale_1 = rows_by_key[("-5", "1")]
            scale_2 = rows_by_key[("-5", "2")]
            self.assertAlmostEqual(
                float(scale_2["qG_full_self_A_per_um"]),
                2.0 * float(scale_1["qG_full_self_A_per_um"]),
                delta=abs(float(scale_1["qG_full_self_A_per_um"])) * 1.0e-10,
            )
            self.assertAlmostEqual(
                float(scale_2["qG_full_ratio"]),
                2.0 * float(scale_1["qG_full_ratio"]),
                delta=abs(float(scale_1["qG_full_ratio"])) * 1.0e-10,
            )
            self.assertEqual(rows_by_key[("-1", "1")]["signal_status"], "low_signal")
            self.assertEqual(rows_by_key[("-1", "1")]["qG_full_ratio"], "")

            summary = summary_md.read_text(encoding="utf-8")
            self.assertIn("-5V", summary)
            self.assertIn("single A_scale", summary)
            self.assertIn("sentaurus_fit_A_B", summary)
            self.assertIn("low-signal", summary)


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
