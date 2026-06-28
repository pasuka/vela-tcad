#!/usr/bin/env python3
"""Regression coverage for the alpha-from-Sentaurus-F diagnostic tool."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EXE = REPO / "build" / ("alpha_from_sentaurus_f.exe" if sys.platform == "win32" else "alpha_from_sentaurus_f")


class AlphaFromSentaurusFDiagnosticTest(unittest.TestCase):
    def test_tool_reconstructs_qf_gradients_and_writes_ratio_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_alpha_from_s_f_") as td:
            tmp = Path(td)
            compare_csv = tmp / "coarse_node_field_compare_aligned.csv"
            elements_csv = tmp / "elements.csv"
            out_csv = tmp / "alpha_from_sentaurus_F.csv"
            summary_md = tmp / "alpha_from_sentaurus_F_summary.md"

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
            for bias in ("-5", "-10"):
                for node_id, x_um, y_um in nodes:
                    e_qf_s = 20.0 * x_um + 30.0 * y_um
                    h_qf_s = -40.0 * x_um + 10.0 * y_um
                    e_qf_self = 22.0 * x_um + 30.0 * y_um
                    h_qf_self = -40.0 * x_um + 8.0 * y_um
                    e_alpha_s = "23000" if x_um == 1.0 else "1e8"
                    h_alpha_s = "11000" if x_um == 1.0 else "1e8"
                    for quantity, sentaurus, self_value in [
                        ("electron_qf", e_qf_s, e_qf_self),
                        ("hole_qf", h_qf_s, h_qf_self),
                        ("electron_alpha_avalanche", e_alpha_s, "1"),
                        ("hole_alpha_avalanche", h_alpha_s, "1"),
                    ]:
                        rows.append({
                            "bias_V": bias,
                            "node_id": str(node_id),
                            "x_um": str(x_um),
                            "y_um": str(y_um),
                            "quantity": quantity,
                            "sentaurus_value": str(sentaurus),
                            "vela_value_scaled_to_sentaurus_units": str(self_value),
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
            self.assertIn("Fn_sentaurus_V_per_cm", first)
            self.assertIn("alpha_n_self_from_S_F_cm_inv", first)
            self.assertAlmostEqual(float(first["Fn_sentaurus_V_per_cm"]), 360555.127546, delta=1.0e-3)
            self.assertAlmostEqual(float(first["Fp_sentaurus_V_per_cm"]), 412310.562562, delta=1.0e-3)
            self.assertGreater(float(first["alpha_n_self_from_S_F_cm_inv"]), 0.0)

            summary = summary_md.read_text(encoding="utf-8")
            self.assertIn("-5V", summary)
            self.assertIn("-10V", summary)
            self.assertIn("\u5269\u4f59\u5dee\u5f02\u96c6\u4e2d\u5728\u4e2d\u4f4e\u573a\u80a9\u90e8", summary)


if __name__ == "__main__":
    unittest.main()
