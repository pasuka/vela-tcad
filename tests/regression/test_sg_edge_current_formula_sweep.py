"""Regression coverage for SG edge current formula sweep diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_sg_edge_current_formula_sweep.py"
Q = 1.602176634e-19


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def compare_row(
    bias: float,
    node_id: int,
    x_um: float,
    y_um: float,
    quantity: str,
    sentaurus: float,
    self_value: float,
) -> list[object]:
    return [
        bias,
        bias,
        node_id,
        x_um,
        y_um,
        quantity,
        "",
        sentaurus,
        "",
        self_value,
        "",
        "",
        "",
        "",
    ]


class SGEdgeCurrentFormulaSweepTest(unittest.TestCase):
    def test_synthetic_qf_arithmetic_candidate_matches_sentaurus_reference(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sg_formula_sweep_") as td:
            tmp = Path(td)
            internal = tmp / "internal.csv"
            topology = tmp / "sg_edges.csv"
            node_compare = tmp / "node_compare.csv"
            out_csv = tmp / "sg_edge_current_formula_sweep.csv"
            out_md = tmp / "sg_edge_current_formula_sweep_summary.md"

            bias = -20.0
            length_um = 1.0
            length_cm = length_um * 1.0e-4
            mu_n = 1000.0
            mu_p = 400.0
            n = 2.0e10
            p = 3.0e10
            phin0 = 0.0
            phin1 = -0.05
            phip0 = 0.0
            phip1 = -0.02
            jn_ref = -Q * mu_n * n * (phin1 - phin0) / length_cm
            jp_ref = -Q * mu_p * p * (phip1 - phip0) / length_cm

            write_csv(
                internal,
                [
                    "source_location_type",
                    "source_entity_id",
                    "bias_V",
                    "x_um",
                    "y_um",
                    "Fn_used_V_per_cm",
                    "Fp_used_V_per_cm",
                    "alpha_n_used_cm_inv",
                    "alpha_p_used_cm_inv",
                    "Jn_mag_used_A_per_cm2",
                    "Jp_mag_used_A_per_cm2",
                ],
                [["edge", 0, bias, 0.75, 0.0, 500.0, 200.0, 1.0, 1.0, jn_ref, jp_ref]],
            )
            write_csv(
                topology,
                [
                    "point_index",
                    "bias_V",
                    "edge_id",
                    "node0",
                    "node1",
                    "x0_um",
                    "y0_um",
                    "x1_um",
                    "y1_um",
                    "edge_length_m",
                    "edge_area_proxy_m2",
                ],
                [[0, bias, 0, 0, 1, 0.25, 0.0, 1.25, 0.0, length_um * 1.0e-6, 1.0e-12]],
            )

            header = [
                "bias_V",
                "vela_bias_V",
                "node_id",
                "x_um",
                "y_um",
                "quantity",
                "sentaurus_field",
                "sentaurus_value",
                "vela_field",
                "vela_value_scaled_to_sentaurus_units",
                "comparison_basis",
                "vela_over_sentaurus",
                "diff",
                "abs_diff",
            ]
            rows: list[list[object]] = []
            for node_id, x_um, phin, phip in ((0, 0.25, phin0, phip0), (1, 1.25, phin1, phip1)):
                rows.extend([
                    compare_row(bias, node_id, x_um, 0.0, "potential", 0.0, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "electron_qf", phin, phin),
                    compare_row(bias, node_id, x_um, 0.0, "hole_qf", phip, phip),
                    compare_row(bias, node_id, x_um, 0.0, "electron_density", n, n),
                    compare_row(bias, node_id, x_um, 0.0, "hole_density", p, p),
                    compare_row(bias, node_id, x_um, 0.0, "electron_mobility", mu_n, mu_n),
                    compare_row(bias, node_id, x_um, 0.0, "hole_mobility", mu_p, mu_p),
                    compare_row(bias, node_id, x_um, 0.0, "electron_current_density_x", jn_ref, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "electron_current_density_y", 0.0, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "hole_current_density_x", jp_ref, 0.0),
                    compare_row(bias, node_id, x_um, 0.0, "hole_current_density_y", 0.0, 0.0),
                ])
            write_csv(node_compare, header, rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--internal-audit",
                    str(internal),
                    "--self-edge-topology",
                    str(topology),
                    "--node-compare",
                    str(node_compare),
                    "--biases",
                    "-20",
                    "--out-csv",
                    str(out_csv),
                    "--out-summary",
                    str(out_md),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with out_csv.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            self.assertAlmostEqual(float(row["Jn_qf_arithmetic_density_A_per_cm2"]), jn_ref)
            self.assertAlmostEqual(float(row["Jp_qf_arithmetic_density_A_per_cm2"]), jp_ref)
            self.assertAlmostEqual(float(row["Jn_qf_arithmetic_density_ratio_to_sent_vector_mag"]), 1.0)
            self.assertAlmostEqual(float(row["Jp_qf_arithmetic_density_ratio_to_sent_vector_mag"]), 1.0)
            self.assertEqual(row["Jn_signal_label"], "high_signal")
            self.assertEqual(row["window_label"], "left_shoulder")

            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("Best electron formula vs Sentaurus edge vector magnitude: qf_arithmetic_density", summary)
            self.assertIn("Current self implementation closest candidate: electron=qf_arithmetic_density", summary)
            self.assertIn("small-signal edges are marked", summary)


if __name__ == "__main__":
    unittest.main()
