"""Regression coverage for local dual-face flux decomposition."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "decompose_dual_face_flux_for_hotspot.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def compare_row(bias: float, node: int, x: float, y: float, quantity: str, self_value: float) -> list[object]:
    return [bias, bias, node, x, y, quantity, "", self_value, "", self_value, "", "", "", ""]


class DualFaceFluxDecompositionTest(unittest.TestCase):
    def test_decomposes_face_contributions_and_dominant_jp_face(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_dual_face_flux_") as td:
            tmp = Path(td)
            dual = tmp / "dual_faces.csv"
            internal = tmp / "internal.csv"
            compare = tmp / "node_compare.csv"
            out_csv = tmp / "face_flux.csv"
            out_md = tmp / "face_flux.md"
            bias = -20.0

            write_csv(
                dual,
                [
                    "dual_face_id",
                    "cell_id",
                    "owner_node_id",
                    "neighbor_node_id",
                    "associated_primal_edge_id",
                    "dual_face_mid_x_um",
                    "dual_face_mid_y_um",
                    "dual_face_normal_x",
                    "dual_face_normal_y",
                    "dual_face_length_cm",
                    "orientation_sign",
                    "boundary_type",
                    "dual_type",
                ],
                [
                    [10, 7, 0, 1, 100, 0.5, 0.0, 1.0, 0.0, 2.0, 1.0, "boundary", "median_dual"],
                    [11, 7, 1, 2, 101, 0.5, 0.5, 0.0, 1.0, 1.0, 1.0, "subcv_internal", "median_dual"],
                    [12, 7, 2, 0, "", 0.0, 0.5, -1.0, 0.0, 1.0, 1.0, "subcv_internal", "median_dual"],
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
                ],
                [
                    ["edge", 100, bias, 1.0, 2.0, 3.0, 4.0, 5.0, 20.0, 1.0],
                    ["edge", 101, bias, 1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 1.0],
                ],
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
            for node, x, y, phip, p, mu in (
                (0, 0.0, 0.0, 0.0, 1.0e10, 100.0),
                (1, 1.0, 0.0, 0.2, 2.0e10, 200.0),
                (2, 0.0, 1.0, 0.1, 3.0e10, 300.0),
            ):
                rows.extend([
                    compare_row(bias, node, x, y, "hole_qf", phip),
                    compare_row(bias, node, x, y, "hole_density", p),
                    compare_row(bias, node, x, y, "hole_mobility", mu),
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
                    "--node-compare",
                    str(compare),
                    "--bias",
                    str(bias),
                    "--cell-id",
                    "7",
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
            self.assertEqual(len(rows_out), 3)
            dominant = max(rows_out, key=lambda row: abs(float(row["Jp_face_scalar_A_per_cm2"] or 0.0)))
            self.assertEqual(dominant["face_id"], "10")
            self.assertNotEqual(dominant["Jp_vector_reconstruction_contribution_x_A_per_cm2"], "")
            self.assertNotEqual(dominant["residual_p_A_per_cm2"], "")
            self.assertEqual(dominant["associated_nodes"], "0;1")
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("1. Which face dominates reconstructed Jp magnitude?", summary)
            self.assertIn("face_id=10", summary)
            self.assertIn("4. Would changing source mapping from cell centroid to node average reduce the hotspot?", summary)


if __name__ == "__main__":
    unittest.main()
