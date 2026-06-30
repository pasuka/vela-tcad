"""Regression tests for dual-face normal-flux vector current reconstruction."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "diagnose_dual_face_vector_current_gava_reconstruction.py"
Q = 1.602176634e-19


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def compare_row(bias: float, node_id: int, x_um: float, y_um: float, quantity: str, sent: float) -> list[object]:
    return [bias, bias, node_id, x_um, y_um, quantity, "", sent, "", sent, "", "", "", ""]


class DualFaceVectorCurrentGavaReconstructionTest(unittest.TestCase):
    def test_true_dual_faces_reconstruct_cell_vector_gava(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_dual_face_vector_") as td:
            tmp = Path(td)
            dual = tmp / "dual_faces.csv"
            internal = tmp / "internal.csv"
            elements = tmp / "elements.csv"
            compare = tmp / "compare.csv"
            out_csv = tmp / "dual_reconstruction.csv"
            out_md = tmp / "dual_reconstruction.md"

            bias = -20.0
            alpha_n = 2.0
            alpha_p = 3.0
            jn = 4.0
            jp = 5.0
            area_cm2 = 0.5e-8
            sent_gava = (alpha_n * jn + alpha_p * jp) / Q

            write_csv(
                dual,
                [
                    "dual_face_id",
                    "cell_id",
                    "owner_node_id",
                    "neighbor_node_id",
                    "associated_primal_edge_id",
                    "dual_face_vertex0_x_um",
                    "dual_face_vertex0_y_um",
                    "dual_face_vertex1_x_um",
                    "dual_face_vertex1_y_um",
                    "dual_face_mid_x_um",
                    "dual_face_mid_y_um",
                    "dual_face_normal_x",
                    "dual_face_normal_y",
                    "dual_face_length_cm",
                    "orientation_sign",
                    "owner_subcv_area_cm2",
                    "neighbor_subcv_area_cm2",
                    "material_region",
                    "boundary_type",
                    "dual_type",
                ],
                [
                    [0, 0, 0, 1, 10, "", "", "", "", 0.5, 0.0, 1.0, 0.0, 1.0e-5, 1, area_cm2 / 3, area_cm2 / 3, "R.Si", "interior_bulk", "median_dual"],
                    [1, 0, 1, 2, 11, "", "", "", "", 0.5, 0.5, 0.0, 1.0, 1.0e-5, 1, area_cm2 / 3, area_cm2 / 3, "R.Si", "interior_bulk", "median_dual"],
                    [2, 0, 2, 0, 12, "", "", "", "", 0.0, 0.5, -1.0, 0.0, 1.0e-5, -1, area_cm2 / 3, area_cm2 / 3, "R.Si", "interior_bulk", "median_dual"],
                ],
            )
            write_csv(
                internal,
                [
                    "source_location_type",
                    "source_entity_id",
                    "bias_V",
                    "Fn_used_V_per_cm",
                    "Fp_used_V_per_cm",
                    "alpha_n_used_cm_inv",
                    "alpha_p_used_cm_inv",
                    "Jn_mag_used_A_per_cm2",
                    "Jp_mag_used_A_per_cm2",
                    "source_weight_or_volume_cm2_for_2D",
                    "qG_contribution_A_per_um",
                ],
                [
                    ["edge", 10, bias, 1.0, 1.0, alpha_n, alpha_p, jn, jp, area_cm2 / 3, 0.0],
                    ["edge", 11, bias, 1.0, 1.0, alpha_n, alpha_p, 0.0, 0.0, area_cm2 / 3, 0.0],
                    ["edge", 12, bias, 1.0, 1.0, alpha_n, alpha_p, jn, jp, area_cm2 / 3, 0.0],
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
            for node_id, x_um, y_um in ((0, 0.0, 0.0), (1, 1.0, 0.0), (2, 0.0, 1.0)):
                rows.extend([
                    compare_row(bias, node_id, x_um, y_um, "electron_alpha_avalanche", alpha_n),
                    compare_row(bias, node_id, x_um, y_um, "hole_alpha_avalanche", alpha_p),
                    compare_row(bias, node_id, x_um, y_um, "avalanche", sent_gava),
                ])
            write_csv(compare, header, rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--dual-face-csv",
                    str(dual),
                    "--internal-audit",
                    str(internal),
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
            self.assertEqual(len(rows_out), 1)
            row = rows_out[0]
            self.assertEqual(row["reconstruction_status"], "ok")
            self.assertAlmostEqual(float(row["Jn_mag_A_per_cm2"]), jn)
            self.assertAlmostEqual(float(row["Jp_mag_A_per_cm2"]), jp)
            self.assertAlmostEqual(float(row["Gava_cell_cm_minus3_s_minus1"]), sent_gava)
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("Does true dual-face vector reconstruction improve -20V full qG beyond previous 0.246? yes", summary)


if __name__ == "__main__":
    unittest.main()


