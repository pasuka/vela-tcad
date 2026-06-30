"""Regression coverage for local dual-face LS sensitivity diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_self_max_cell_dual_face_ls_sensitivity.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def compare_row(bias: float, node: int, x: float, y: float, quantity: str, self_value: float) -> list[object]:
    return [bias, bias, node, x, y, quantity, "", self_value, "", self_value, "", "", "", ""]


class SelfMaxCellDualFaceLSSensitivityTest(unittest.TestCase):
    def test_exports_ls_system_missing_faces_and_sensitivity_variants(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_ls_sensitivity_") as td:
            tmp = Path(td)
            dual = tmp / "dual_faces.csv"
            internal = tmp / "internal.csv"
            compare = tmp / "node_compare.csv"
            cells = tmp / "dual_cell.csv"
            out_csv = tmp / "ls_sensitivity.csv"
            out_md = tmp / "ls_sensitivity.md"
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
                    [72, 12, 0, 1, "", 1.10, 0.39, 0.70710678, 0.70710678, 1.0, 1.0, "subcv_internal", "median_dual"],
                    [74, 12, 1, 2, 25, 1.10, 0.46, -0.89442719, 0.44721360, 2.0, 1.0, "boundary", "median_dual"],
                    [75, 12, 2, 1, 25, 1.10, 0.46, 0.89442719, -0.44721360, 2.0, -1.0, "boundary", "median_dual"],
                    [76, 12, 2, 0, 26, 1.04, 0.39, 0.44721360, -0.89442719, 1.0, 1.0, "interior_bulk", "median_dual"],
                    [77, 12, 0, 2, 26, 1.04, 0.39, -0.44721360, 0.89442719, 1.0, -1.0, "interior_bulk", "median_dual"],
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
                    "electron_raw_signed_flux_proxy",
                    "hole_raw_signed_flux_proxy",
                ],
                [
                    ["edge", 25, bias, 10.0, 20.0, 30.0, 40.0, 1.0e-11, 2.0e-7, 1.0, 1.0],
                    ["edge", 26, bias, 10.0, 20.0, 30.0, 40.0, 5.0e-13, 8.0e-9, 1.0, 1.0],
                ],
            )
            write_csv(
                cells,
                [
                    "bias_V",
                    "cell_id",
                    "centroid_x_um",
                    "centroid_y_um",
                    "reconstruction_status",
                    "Jp_x_A_per_cm2",
                    "Jp_y_A_per_cm2",
                    "Jp_mag_A_per_cm2",
                    "Jn_x_A_per_cm2",
                    "Jn_y_A_per_cm2",
                    "Jn_mag_A_per_cm2",
                    "alpha_n_cell_cm_inv",
                    "alpha_p_cell_cm_inv",
                    "Gava_cell_cm_minus3_s_minus1",
                    "Gava_sent_cell_cm_minus3_s_minus1",
                ],
                [
                    [bias, 12, 1.0833333333, 0.4166666667, "ok", -1.0e-7, -1.0e-7, 1.414e-7, -1.0e-11, -1.0e-11, 1.414e-11, 30.0, 40.0, 4.0e13, 2.0e13],
                    [bias, 13, 1.0833333333, 0.1666666667, "ok", -2.0e-8, 0.0, 2.0e-8, -1.0e-12, 0.0, 1.0e-12, 30.0, 40.0, 5.0e12, 5.0e12],
                ],
            )
            compare_header = [
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
            compare_rows: list[list[object]] = []
            for node, x, y in ((0, 1.0, 0.3), (1, 1.1, 0.4), (2, 1.2, 0.5)):
                compare_rows.extend([
                    compare_row(bias, node, x, y, "hole_qf", -10.0 + node),
                    compare_row(bias, node, x, y, "hole_density", 1.0e4 + node),
                    compare_row(bias, node, x, y, "hole_mobility", 25.0),
                ])
            write_csv(compare, compare_header, compare_rows)

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
                    "--self-dual-cell-csv",
                    str(cells),
                    "--bias",
                    str(bias),
                    "--cell-id",
                    "12",
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
                rows = list(csv.DictReader(handle))
            row_types = {row["row_type"] for row in rows}
            self.assertTrue({"ls_system", "face", "missing_flux_constraint", "variant"}.issubset(row_types))

            variants = {row["variant_name"] for row in rows if row["row_type"] == "variant"}
            self.assertIn("remove_faces_74_75", variants)
            self.assertIn("boundary_weight_0p5", variants)
            self.assertIn("tikhonov_1e-04_trace", variants)

            face74 = next(row for row in rows if row["row_type"] == "face" and row["face_id"] == "74" and row["carrier"] == "hole")
            self.assertEqual(face74["is_boundary"], "True")
            self.assertNotEqual(face74["residual_A_per_cm2"], "")

            missing = [row for row in rows if row["row_type"] == "missing_flux_constraint"]
            self.assertEqual(missing[0]["face_id"], "72")
            self.assertIn("no_associated_primal_edge", missing[0]["missing_reason"])

            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("1. Is the LS reconstruction ill-conditioned?", summary)
            self.assertIn("3. Are faces 74/75 necessary for the hotspot?", summary)
            self.assertIn("7. Should boundary cell dual-face vector current use a guarded reconstruction?", summary)


if __name__ == "__main__":
    unittest.main()
