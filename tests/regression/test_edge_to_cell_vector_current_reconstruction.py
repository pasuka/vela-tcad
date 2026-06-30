"""Regression coverage for edge-to-cell vector current reconstruction diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "diagnose_edge_to_cell_vector_current_reconstruction.py"
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


class EdgeToCellVectorCurrentReconstructionTest(unittest.TestCase):
    def test_single_triangle_reconstructs_cell_vector_and_baseline_qg(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_edge_to_cell_vector_") as td:
            tmp = Path(td)
            internal = tmp / "internal.csv"
            edges = tmp / "sg_edges.csv"
            elements = tmp / "elements.csv"
            compare = tmp / "node_compare.csv"
            out_csv = tmp / "edge_to_cell_vector_current_reconstruction.csv"
            out_md = tmp / "edge_to_cell_vector_current_reconstruction_summary.md"

            bias = -20.0
            alpha_n = 2.0
            alpha_p = 3.0
            jn = 4.0
            jp = 5.0
            gava = (alpha_n * jn + alpha_p * jp) / Q
            area_cm2 = 0.5e-8
            qg = Q * gava * area_cm2 / 1.0e4

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
                    "Gava_total_used_cm_minus3_s_minus1",
                    "source_weight_or_volume_cm2_for_2D",
                    "qG_contribution_A_per_um",
                ],
                [
                    ["edge", 0, bias, 0.5, 0.0, 1.0, 1.0, alpha_n, alpha_p, jn, jp, gava, area_cm2, qg],
                    ["edge", 1, bias, 0.5, 0.5, 1.0, 1.0, alpha_n, alpha_p, 0.0, 0.0, 0.0, area_cm2, 0.0],
                    ["edge", 2, bias, 0.0, 0.5, 1.0, 1.0, alpha_n, alpha_p, 0.0, 0.0, 0.0, area_cm2, 0.0],
                ],
            )
            write_csv(
                edges,
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
                [
                    [0, bias, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0e-6, area_cm2 / 1.0e4],
                    [0, bias, 1, 1, 2, 1.0, 0.0, 0.0, 1.0, 2.0 ** 0.5 * 1.0e-6, area_cm2 / 1.0e4],
                    [0, bias, 2, 2, 0, 0.0, 1.0, 0.0, 0.0, 1.0e-6, area_cm2 / 1.0e4],
                ],
            )
            write_csv(elements, ["id", "node0", "node1", "node2", "region", "material"], [[0, 0, 1, 2, "R.Si", "Si"]])

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
            for node_id, x_um, y_um, phin in ((0, 0.0, 0.0, 0.0), (1, 1.0, 0.0, -1.0e-4), (2, 0.0, 1.0, 0.0)):
                rows.extend([
                    compare_row(bias, node_id, x_um, y_um, "electron_qf", phin, phin),
                    compare_row(bias, node_id, x_um, y_um, "hole_qf", phin, phin),
                    compare_row(bias, node_id, x_um, y_um, "electron_density", 1.0e10, 1.0e10),
                    compare_row(bias, node_id, x_um, y_um, "hole_density", 1.0e10, 1.0e10),
                    compare_row(bias, node_id, x_um, y_um, "electron_mobility", 1000.0, 1000.0),
                    compare_row(bias, node_id, x_um, y_um, "hole_mobility", 500.0, 500.0),
                    compare_row(bias, node_id, x_um, y_um, "electron_alpha_avalanche", alpha_n, alpha_n),
                    compare_row(bias, node_id, x_um, y_um, "hole_alpha_avalanche", alpha_p, alpha_p),
                    compare_row(bias, node_id, x_um, y_um, "electron_current_density_x", jn, 0.0),
                    compare_row(bias, node_id, x_um, y_um, "electron_current_density_y", 0.0, 0.0),
                    compare_row(bias, node_id, x_um, y_um, "hole_current_density_x", jp, 0.0),
                    compare_row(bias, node_id, x_um, y_um, "hole_current_density_y", 0.0, 0.0),
                    compare_row(bias, node_id, x_um, y_um, "avalanche", gava, 0.0),
                ])
            write_csv(compare, header, rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--internal-audit",
                    str(internal),
                    "--self-edge-topology",
                    str(edges),
                    "--elements-csv",
                    str(elements),
                    "--node-compare",
                    str(compare),
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
            baseline = next(row for row in rows_out if row["mode"] == "current_edge_scalar_existing")
            tangent = next(row for row in rows_out if row["mode"] == "current_cell_vector_ls_from_edge_tangent")
            self.assertAlmostEqual(float(baseline["qG_cell_contribution_A_per_um"]), qg)
            self.assertGreater(float(tangent["Jn_mag_A_per_cm2"]), 0.0)
            self.assertIn("current_cell_vector_ls_from_edge_tangent", out_md.read_text(encoding="utf-8"))
            self.assertIn("edge scalar baseline reproduces current internal qG: yes", out_md.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
