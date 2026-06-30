"""Regression coverage for edge current decomposition diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_edge_current_decomposition.py"
Q = 1.602176634e-19


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def node_row(
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
        "scalar_or_magnitude",
        "",
        "",
        "",
    ]


class EdgeCurrentDecompositionAuditTest(unittest.TestCase):
    def test_synthetic_edge_decomposes_q_mu_density_field_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_edge_current_decomp_") as td:
            tmp = Path(td)
            internal = tmp / "internal.csv"
            topology = tmp / "sg_edges.csv"
            node_compare = tmp / "coarse_node_field_compare_aligned.csv"
            out_csv = tmp / "edge_current_decomposition_audit.csv"
            out_md = tmp / "edge_current_decomposition_audit_summary.md"

            bias = -20.0
            length_um = 1.0
            fn = 1.0e5
            fp = 2.0e5
            n_sent = 1.0e10
            p_sent = 2.0e10
            n_self = 1.0e9
            p_self = 1.0e10
            mu_n = 1000.0
            mu_p_sent = 500.0
            mu_p_self = 250.0
            jn_sent = Q * mu_n * n_sent * fn
            jp_sent = Q * mu_p_sent * p_sent * fp
            jn_self = Q * mu_n * n_self * fn
            jp_self = Q * mu_p_self * p_self * fp

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
                    "Gava_n_used_cm_minus3_s_minus1",
                    "Gava_p_used_cm_minus3_s_minus1",
                    "Gava_total_used_cm_minus3_s_minus1",
                    "source_weight_or_volume_cm2_for_2D",
                    "qG_contribution_A_per_um",
                ],
                [[
                    "edge",
                    0,
                    bias,
                    0.75,
                    0.0,
                    fn,
                    fp,
                    1.0,
                    1.0,
                    jn_self,
                    jp_self,
                    0.0,
                    0.0,
                    0.0,
                    1.0e-8,
                    0.0,
                ]],
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
            for node_id, x_um, electron_qf, hole_qf in (
                (0, 0.25, 0.0, 0.0),
                (1, 1.25, 10.0, 20.0),
            ):
                rows.extend([
                    node_row(bias, node_id, x_um, 0.0, "electron_density", n_sent, n_self),
                    node_row(bias, node_id, x_um, 0.0, "hole_density", p_sent, p_self),
                    node_row(bias, node_id, x_um, 0.0, "electron_mobility", mu_n, mu_n),
                    node_row(bias, node_id, x_um, 0.0, "hole_mobility", mu_p_sent, mu_p_self),
                    node_row(bias, node_id, x_um, 0.0, "electron_qf", electron_qf, electron_qf),
                    node_row(bias, node_id, x_um, 0.0, "hole_qf", hole_qf, hole_qf),
                    node_row(bias, node_id, x_um, 0.0, "electron_current_density_x", jn_sent, 0.0),
                    node_row(bias, node_id, x_um, 0.0, "electron_current_density_y", 0.0, 0.0),
                    node_row(bias, node_id, x_um, 0.0, "hole_current_density_x", jp_sent, 0.0),
                    node_row(bias, node_id, x_um, 0.0, "hole_current_density_y", 0.0, 0.0),
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
                rows_out = list(csv.DictReader(handle))
            self.assertEqual(len(rows_out), 1)
            row = rows_out[0]
            self.assertEqual(row["window_label"], "left_shoulder")
            self.assertAlmostEqual(float(row["Jn_self_closure_ratio"]), 1.0)
            self.assertAlmostEqual(float(row["Jp_self_closure_ratio"]), 1.0)
            self.assertAlmostEqual(float(row["Jn_ratio_self_over_sent"]), 0.1)
            self.assertAlmostEqual(float(row["Jp_ratio_self_over_sent"]), 0.25)
            self.assertEqual(row["dominant_Jn_difference_factor"], "carrier_density")
            self.assertEqual(row["dominant_Jp_difference_factor"], "carrier_density")

            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("self Jn/Jp can be reconstructed from self q*mu*n*F: yes", summary)
            self.assertIn("Sentaurus edge Jn/Jp can be reconstructed from Sentaurus q*mu*n*F: yes", summary)
            self.assertIn("self exported node current density used: no", summary)


if __name__ == "__main__":
    unittest.main()
