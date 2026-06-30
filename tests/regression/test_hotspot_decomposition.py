"""Regression coverage for top-Gava hotspot decomposition diagnostics."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "decompose_gava_hotspots.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


class HotspotDecompositionTest(unittest.TestCase):
    def test_lists_sentaurus_and_self_hotspots_with_counterpart_distance(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_hotspot_decomp_") as td:
            tmp = Path(td)
            sent_case = tmp / "sentaurus_multibias" / "sentaurus_-20v"
            fields = sent_case / "fields"
            write_csv(sent_case / "nodes.csv", ["id", "x_um", "y_um"], [[0, 1.0, 0.0], [1, 1.2, 0.0], [2, 1.0, 1.0]])
            write_csv(sent_case / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [[0, 0, 1, 2, "R.Si", "Si"]])
            write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [[0, 10.0], [1, 30.0], [2, 20.0]])
            write_csv(fields / "eAlphaAvalanche_region0.csv", ["node_id", "component0"], [[0, 2.0], [1, 3.0], [2, 2.0]])
            write_csv(fields / "hAlphaAvalanche_region0.csv", ["node_id", "component0"], [[0, 4.0], [1, 5.0], [2, 4.0]])
            write_csv(fields / "eCurrentDensity_region0.csv", ["node_id", "component0"], [[0, 1.0], [1, 2.0], [2, 1.0]])
            write_csv(fields / "hCurrentDensity_region0.csv", ["node_id", "component0"], [[0, 1.0], [1, 3.0], [2, 1.0]])
            write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0], [1, 0.1], [2, 0.0]])
            write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0], [1, 0.0], [2, 0.2]])

            dual = tmp / "dual.csv"
            write_csv(
                dual,
                [
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
                    "Fn_cell_V_per_cm",
                    "Fp_cell_V_per_cm",
                    "Gava_cell_cm_minus3_s_minus1",
                    "qG_cell_contribution_A_per_um",
                    "window_label",
                ],
                [[-20.0, 7, 1.1, 0.1, "ok", "current_cell_vector_ls_from_true_dual_face_normals", 4.0, 5.0, 2.0, 3.0, 100.0, 200.0, 40.0, 1.0e-18, "center"]],
            )

            out_csv = tmp / "hotspots.csv"
            out_md = tmp / "hotspots.md"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--sentaurus-multibias-dir",
                    str(tmp / "sentaurus_multibias"),
                    "--self-case-dir",
                    str(tmp / "missing_case"),
                    "--self-edge-csv",
                    str(tmp / "missing_edges.csv"),
                    "--self-dual-cell-csv",
                    str(dual),
                    "--biases",
                    "-20",
                    "--top-n",
                    "2",
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
            self.assertEqual({row["source"] for row in rows}, {"sentaurus_node", "self_dual_cell"})
            sent_top = next(row for row in rows if row["source"] == "sentaurus_node" and row["rank"] == "1")
            self.assertEqual(sent_top["x_um"], "1.2")
            self.assertEqual(sent_top["window_label"], "right_shoulder")
            self.assertNotEqual(sent_top["nearest_counterpart_distance_um"], "")
            self.assertIn("Gava_n_cm_minus3_s_minus1", sent_top)
            self.assertIn("qG_contribution_A_per_um", sent_top)
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("1. Is self hotspot shifted because F/alpha peak shifts?", summary)
            self.assertIn("4. At -20V, why does max remain at (1.08333,0.416667) um?", summary)


if __name__ == "__main__":
    unittest.main()
