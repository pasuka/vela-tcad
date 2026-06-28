#!/usr/bin/env python3
"""Regression coverage for the raw Van Overstraeten alpha comparison report."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class VanOverstraetenRawAlphaCompareTest(unittest.TestCase):
    def test_report_writes_required_columns_regions_and_conclusion(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_raw_alpha_compare_") as td:
            tmp = Path(td)
            input_csv = tmp / "coarse_node_field_compare_aligned.csv"
            output_csv = tmp / "vanoverstraeten_raw_alpha_compare.csv"
            summary = tmp / "vanoverstraeten_raw_alpha_summary.md"

            rows = [
                {
                    "bias_V": "-20",
                    "node_id": "1",
                    "x_um": "0.75",
                    "y_um": "0.25",
                    "Fn_sentaurus_V_per_cm": "200000",
                    "Fn_self_default_V_per_cm": "200000",
                    "Fn_self_raw_V_per_cm": "200000",
                    "alpha_n_sentaurus_cm_inv": "1000",
                    "alpha_n_self_default_cm_inv": "100",
                    "Fp_sentaurus_V_per_cm": "200000",
                    "Fp_self_default_V_per_cm": "200000",
                    "Fp_self_raw_V_per_cm": "200000",
                    "alpha_p_sentaurus_cm_inv": "100",
                    "alpha_p_self_default_cm_inv": "10",
                },
                {
                    "bias_V": "-20",
                    "node_id": "2",
                    "x_um": "1.00",
                    "y_um": "0.25",
                    "Fn_sentaurus_V_per_cm": "500000",
                    "Fn_self_default_V_per_cm": "500000",
                    "Fn_self_raw_V_per_cm": "500000",
                    "alpha_n_sentaurus_cm_inv": "9000",
                    "alpha_n_self_default_cm_inv": "8800",
                    "Fp_sentaurus_V_per_cm": "500000",
                    "Fp_self_default_V_per_cm": "500000",
                    "Fp_self_raw_V_per_cm": "500000",
                    "alpha_p_sentaurus_cm_inv": "2200",
                    "alpha_p_self_default_cm_inv": "2100",
                },
            ]
            with input_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_vanoverstraeten_raw_alpha.py"),
                    "--input-csv",
                    str(input_csv),
                    "--output-csv",
                    str(output_csv),
                    "--summary-md",
                    str(summary),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with output_csv.open(encoding="utf-8") as handle:
                out_rows = list(csv.DictReader(handle))
            self.assertEqual(len(out_rows), 2)
            for name in [
                "Fn_self_raw_V_per_cm",
                "alpha_n_self_raw_cm_inv",
                "alpha_n_default_ratio",
                "alpha_n_raw_ratio",
                "Fp_self_raw_V_per_cm",
                "alpha_p_self_raw_cm_inv",
                "alpha_p_default_ratio",
                "alpha_p_raw_ratio",
            ]:
                self.assertIn(name, out_rows[0])
            self.assertGreater(float(out_rows[0]["alpha_n_self_raw_cm_inv"]), 10.0)

            text = summary.read_text(encoding="utf-8")
            self.assertIn("\u7ed3\u5de6\u80a9\u90e8", text)
            self.assertIn(
                "alpha \u504f\u5c0f\u4e3b\u8981\u6765\u81ea cutoff/smoothing/RefDens \u6291\u5236",
                text,
            )

    def test_report_accepts_space_padded_long_node_compare_csv(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_raw_alpha_long_") as td:
            tmp = Path(td)
            input_csv = tmp / "coarse_node_field_compare_aligned.csv"
            output_csv = tmp / "raw.csv"
            summary = tmp / "summary.md"
            input_csv.write_text(
                "\n".join([
                    "bias_V,node_id,x_um,y_um,quantity                  ,sentaurus_value,vela_value_scaled_to_sentaurus_units",
                    "-5,7,0.75,0.25,electric_field            ,200000,200000",
                    "-5,7,0.75,0.25,electron_alpha_avalanche ,1000,100",
                    "-5,7,0.75,0.25,hole_alpha_avalanche     ,100,10",
                ]),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_vanoverstraeten_raw_alpha.py"),
                    "--input-csv",
                    str(input_csv),
                    "--output-csv",
                    str(output_csv),
                    "--summary-md",
                    str(summary),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with output_csv.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["Fn_self_raw_V_per_cm"], "200000")
            self.assertEqual(rows[0]["alpha_n_self_default_cm_inv"], "100")
            self.assertGreater(float(rows[0]["alpha_n_self_raw_cm_inv"]), 100.0)


if __name__ == "__main__":
    unittest.main()
