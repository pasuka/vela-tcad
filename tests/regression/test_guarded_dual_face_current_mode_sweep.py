"""Regression coverage for guarded dual-face current magnitude diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "sweep_guarded_dual_face_current_modes.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


class GuardedDualFaceCurrentModeSweepTest(unittest.TestCase):
    def test_guarded_modes_reduce_flagged_boundary_hotspot(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_guarded_modes_") as td:
            tmp = Path(td)
            dual_cell = tmp / "dual_cell.csv"
            edge_cell = tmp / "edge_cell.csv"
            dual_faces = tmp / "dual_faces.csv"
            internal = tmp / "internal.csv"
            ls = tmp / "ls.csv"
            out_csv = tmp / "guarded.csv"
            out_md = tmp / "guarded.md"

            cell_header = [
                "bias_V",
                "cell_id",
                "centroid_x_um",
                "centroid_y_um",
                "reconstruction_status",
                "mode",
                "Jn_mag_A_per_cm2",
                "Jp_mag_A_per_cm2",
                "alpha_n_cell_cm_inv",
                "alpha_p_cell_cm_inv",
                "Gava_cell_cm_minus3_s_minus1",
                "qG_cell_contribution_A_per_um",
                "Gava_sent_cell_cm_minus3_s_minus1",
                "qG_sent_cell_contribution_A_per_um",
                "window_label",
            ]
            write_csv(
                dual_cell,
                cell_header,
                [
                    [-20, 12, 1.083333, 0.416667, "ok", "current_cell_vector_ls_from_true_dual_face_normals", 1.0e-11, 9.0e-7, 10.0, 10.0, 5.7e13, 5.7e-17, 3.0e13, 3.0e-17, "center"],
                    [-20, 13, 1.083333, 0.166667, "ok", "current_cell_vector_ls_from_true_dual_face_normals", 1.0e-11, 2.0e-7, 10.0, 10.0, 1.3e13, 1.3e-17, 4.0e13, 4.0e-17, "center"],
                    [-20, 14, 1.20, 0.20, "ok", "current_cell_vector_ls_from_true_dual_face_normals", 1.0e-11, 4.0e-7, 10.0, 10.0, 2.6e13, 2.6e-17, 2.0e13, 2.0e-17, "right_shoulder"],
                ],
            )
            write_csv(
                edge_cell,
                cell_header,
                [
                    [-20, 12, 1.083333, 0.416667, "ok", "current_edge_scalar_existing", 1.0e-11, 1.0e-7, 10.0, 10.0, 6.0e12, 6.0e-18, 3.0e13, 3.0e-17, "center"],
                    [-20, 13, 1.083333, 0.166667, "ok", "current_edge_scalar_existing", 1.0e-11, 2.0e-7, 10.0, 10.0, 1.3e13, 1.3e-17, 4.0e13, 4.0e-17, "center"],
                    [-20, 14, 1.20, 0.20, "ok", "current_edge_scalar_existing", 1.0e-11, 2.0e-7, 10.0, 10.0, 1.3e13, 1.3e-17, 2.0e13, 2.0e-17, "right_shoulder"],
                ],
            )
            write_csv(
                dual_faces,
                [
                    "dual_face_id",
                    "cell_id",
                    "associated_primal_edge_id",
                    "dual_face_normal_x",
                    "dual_face_normal_y",
                    "dual_face_length_cm",
                    "orientation_sign",
                    "boundary_type",
                    "dual_type",
                ],
                [
                    [74, 12, 25, -0.894, 0.447, 1.0, 1.0, "boundary", "median_dual"],
                    [75, 12, 25, 0.894, -0.447, 1.0, -1.0, "boundary", "median_dual"],
                    [76, 12, 26, 0.447, -0.894, 1.0, 1.0, "interior_bulk", "median_dual"],
                    [77, 12, 26, -0.447, 0.894, 1.0, -1.0, "interior_bulk", "median_dual"],
                    [80, 13, 27, 1.0, 0.0, 1.0, 1.0, "interior_bulk", "median_dual"],
                    [81, 13, 28, 0.0, 1.0, 1.0, 1.0, "interior_bulk", "median_dual"],
                    [82, 14, 29, 1.0, 0.0, 1.0, 1.0, "interior_bulk", "median_dual"],
                ],
            )
            write_csv(
                internal,
                [
                    "source_location_type",
                    "source_entity_id",
                    "bias_V",
                    "Jn_mag_used_A_per_cm2",
                    "Jp_mag_used_A_per_cm2",
                    "electron_raw_signed_flux_proxy",
                    "hole_raw_signed_flux_proxy",
                ],
                [
                    ["edge", 25, -20, 1.0e-11, 9.0e-7, 1.0, 1.0],
                    ["edge", 26, -20, 1.0e-11, 1.0e-8, 1.0, 1.0],
                    ["edge", 27, -20, 1.0e-11, 2.0e-7, 1.0, 1.0],
                    ["edge", 28, -20, 1.0e-11, 2.0e-7, 1.0, 1.0],
                    ["edge", 29, -20, 1.0e-11, 4.0e-7, 1.0, 1.0],
                ],
            )
            write_csv(
                ls,
                ["row_type", "bias_V", "cell_id", "variant_name", "ratio_to_baseline", "hotspot_would_remain_rank1"],
                [["variant", -20, 12, "remove_faces_74_75", 0.02, "False"]],
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--dual-cell-csv",
                    str(dual_cell),
                    "--edge-cell-csv",
                    str(edge_cell),
                    "--dual-face-csv",
                    str(dual_faces),
                    "--internal-audit",
                    str(internal),
                    "--ls-sensitivity-csv",
                    str(ls),
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
            modes = {row["current_magnitude_mode"] for row in rows}
            self.assertIn("dual_face_vector_mag", modes)
            self.assertIn("guarded_fallback_edge_scalar", modes)
            self.assertIn("guarded_ratio_limiter_2", modes)
            fallback = next(row for row in rows if row["current_magnitude_mode"] == "guarded_fallback_edge_scalar")
            self.assertEqual(fallback["cell12_remains_rank1"], "False")
            self.assertLess(float(fallback["cell12_Gava_cm_minus3_s_minus1"]), 1.0e13)
            limiter = next(row for row in rows if row["current_magnitude_mode"] == "guarded_ratio_limiter_2")
            self.assertLess(float(limiter["cell12_Jp_A_per_cm2"]), 3.0e-7)
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("1. Which guarded mode best reduces the artificial hotspot at cell 12?", summary)
            self.assertIn("6. Do not recommend unguarded dual_face_vector_mag as default.", summary)


if __name__ == "__main__":
    unittest.main()
