#!/usr/bin/env python3
"""Regression tests for current-density vector audit diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
Q = 1.602176634e-19


def write_compare_row(
    writer: csv.DictWriter,
    *,
    bias: float,
    node: int,
    x_um: float,
    y_um: float,
    quantity: str,
    sentaurus: float,
    self_value: float,
) -> None:
    ratio = "" if sentaurus == 0.0 else self_value / sentaurus
    diff = self_value - sentaurus
    writer.writerow({
        "bias_V": bias,
        "vela_bias_V": bias,
        "node_id": node,
        "x_um": x_um,
        "y_um": y_um,
        "quantity": quantity,
        "sentaurus_field": quantity,
        "sentaurus_value": sentaurus,
        "vela_field": quantity,
        "vela_value_scaled_to_sentaurus_units": self_value,
        "comparison_basis": "synthetic",
        "vela_over_sentaurus": ratio,
        "diff": diff,
        "abs_diff": abs(diff),
    })


class CurrentDensityVectorAuditTest(unittest.TestCase):
    def test_audit_reports_scale_sign_and_gava_closure(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_current_vector_audit_") as tmp:
            root = Path(tmp)
            compare = root / "coarse_node_field_compare_aligned.csv"
            out_csv = root / "audit.csv"
            out_summary = root / "audit.md"
            fields = [
                "bias_V", "vela_bias_V", "node_id", "x_um", "y_um", "quantity",
                "sentaurus_field", "sentaurus_value", "vela_field",
                "vela_value_scaled_to_sentaurus_units", "comparison_basis",
                "vela_over_sentaurus", "diff", "abs_diff",
            ]
            with compare.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                # Electron self current is a sign-flipped 2x version of Sentaurus.
                for node, x_um, jn_sent, jp_sent in [
                    (0, 0.80, (3.0, 4.0), (0.0, 2.0)),
                    (1, 1.00, (6.0, 8.0), (0.0, 4.0)),
                    (2, 1.20, (9.0, 12.0), (0.0, 6.0)),
                ]:
                    jn_self = (-2.0 * jn_sent[0], -2.0 * jn_sent[1])
                    jp_self = (3.0 * jp_sent[0], 3.0 * jp_sent[1])
                    total_sent = (jn_sent[0] + jp_sent[0], jn_sent[1] + jp_sent[1])
                    total_self = (jn_self[0] + jp_self[0], jn_self[1] + jp_self[1])
                    values = {
                        "electron_current_density": (jn_sent, jn_self),
                        "hole_current_density": (jp_sent, jp_self),
                        "total_current_density": (total_sent, total_self),
                        "electric_field": ((1.0, 0.0), (1.01, 0.0)),
                    }
                    for prefix, (sent, self_value) in values.items():
                        write_compare_row(writer, bias=-5.0, node=node, x_um=x_um, y_um=0.25,
                                          quantity=f"{prefix}_x", sentaurus=sent[0], self_value=self_value[0])
                        write_compare_row(writer, bias=-5.0, node=node, x_um=x_um, y_um=0.25,
                                          quantity=f"{prefix}_y", sentaurus=sent[1], self_value=self_value[1])
                        write_compare_row(writer, bias=-5.0, node=node, x_um=x_um, y_um=0.25,
                                          quantity=f"{prefix}_mag",
                                          sentaurus=(sent[0] ** 2 + sent[1] ** 2) ** 0.5,
                                          self_value=(self_value[0] ** 2 + self_value[1] ** 2) ** 0.5)
                    write_compare_row(writer, bias=-5.0, node=node, x_um=x_um, y_um=0.25,
                                      quantity="electron_alpha_avalanche", sentaurus=10.0, self_value=10.0)
                    write_compare_row(writer, bias=-5.0, node=node, x_um=x_um, y_um=0.25,
                                      quantity="hole_alpha_avalanche", sentaurus=5.0, self_value=5.0)
                    sent_g = (10.0 * (jn_sent[0] ** 2 + jn_sent[1] ** 2) ** 0.5
                              + 5.0 * (jp_sent[0] ** 2 + jp_sent[1] ** 2) ** 0.5) / Q * 1.0e-6
                    self_g = (10.0 * (jn_self[0] ** 2 + jn_self[1] ** 2) ** 0.5
                              + 5.0 * (jp_self[0] ** 2 + jp_self[1] ** 2) ** 0.5) / Q * 1.0e-6
                    write_compare_row(writer, bias=-5.0, node=node, x_um=x_um, y_um=0.25,
                                      quantity="avalanche", sentaurus=sent_g, self_value=self_g)

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "audit_current_density_vectors.py"),
                    "--node-compare-csv",
                    str(compare),
                    "--out-csv",
                    str(out_csv),
                    "--summary-md",
                    str(out_summary),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with out_csv.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            electron_full = next(
                row for row in rows
                if row["row_type"] == "vector_summary"
                and row["carrier"] == "electron"
                and row["window"] == "full"
            )
            self.assertEqual(electron_full["sign_flip_better"], "1")
            self.assertAlmostEqual(float(electron_full["fit_scale_self_to_sentaurus"]), -0.5)
            gava_self = next(row for row in rows if row["row_type"] == "gava_summary" and row["side"] == "self")
            self.assertLess(float(gava_self["median_gava_rec_rel_error_A_per_m2_assumption"]), 1.0e-12)
            summary = out_summary.read_text(encoding="utf-8")
            self.assertIn("CurrentDensity vector sign convention", summary)
            self.assertIn("self exported current density can reconstruct self AvalancheGeneration", summary)


if __name__ == "__main__":
    unittest.main()
