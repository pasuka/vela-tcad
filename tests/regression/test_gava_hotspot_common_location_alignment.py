"""Regression coverage for common-location Gava hotspot alignment."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "align_gava_hotspots_common_locations.py"


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


class GavaHotspotCommonLocationAlignmentTest(unittest.TestCase):
    def test_outputs_cell_and_node_common_location_hotspots(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_common_hotspot_") as td:
            tmp = Path(td)
            sent = tmp / "sentaurus_multibias" / "sentaurus_-20v"
            fields = sent / "fields"
            write_csv(sent / "nodes.csv", ["id", "x_um", "y_um"], [[0, 1.0, 0.0], [1, 1.2, 0.0], [2, 1.0, 0.6], [3, 1.2, 0.6]])
            write_csv(sent / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [[0, 0, 1, 2, "R.Si", "Si"], [1, 1, 3, 2, "R.Si", "Si"]])
            write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [[0, 10.0], [1, 30.0], [2, 20.0], [3, 5.0]])
            write_csv(fields / "eAlphaAvalanche_region0.csv", ["node_id", "component0"], [[0, 2.0], [1, 3.0], [2, 2.0], [3, 1.0]])
            write_csv(fields / "hAlphaAvalanche_region0.csv", ["node_id", "component0"], [[0, 4.0], [1, 5.0], [2, 4.0], [3, 1.0]])
            write_csv(fields / "eCurrentDensity_region0.csv", ["node_id", "component0"], [[0, 1.0], [1, 2.0], [2, 1.0], [3, 1.0]])
            write_csv(fields / "hCurrentDensity_region0.csv", ["node_id", "component0"], [[0, 1.0], [1, 3.0], [2, 1.0], [3, 1.0]])
            write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0], [1, 0.1], [2, 0.0], [3, 0.1]])
            write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0], [1, 0.0], [2, 0.2], [3, 0.2]])

            dual = tmp / "dual.csv"
            write_csv(
                dual,
                [
                    "bias_V",
                    "cell_id",
                    "centroid_x_um",
                    "centroid_y_um",
                    "reconstruction_status",
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
                [
                    [-20.0, 0, 1.0666666667, 0.2, "ok", 1.0, 10.0, 2.0, 3.0, 100.0, 200.0, 100.0, 1.0e-18, "center"],
                    [-20.0, 1, 1.1333333333, 0.4, "ok", 1.0, 1.0, 1.0, 1.0, 50.0, 60.0, 10.0, 1.0e-19, "center"],
                ],
            )

            out_csv = tmp / "alignment.csv"
            out_md = tmp / "alignment.md"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--sentaurus-multibias-dir",
                    str(tmp / "sentaurus_multibias"),
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
            source_types = {row["source_type"] for row in rows}
            self.assertIn("sentaurus_cell_average", source_types)
            self.assertIn("self_dual_cell", source_types)
            self.assertIn("self_node_recovered", source_types)
            self.assertIn("sentaurus_node", source_types)
            self.assertTrue(any(row["comparison_case"] == "sentaurus_max_node_nearest_self_node" for row in rows))
            self.assertTrue(any(row["comparison_case"] == "self_max_cell_nearest_sentaurus_cell" for row in rows))
            first_self_node = next(row for row in rows if row["source_type"] == "self_node_recovered" and row["rank"] == "1")
            self.assertNotEqual(first_self_node["nearest_counterpart_distance_um"], "")
            self.assertNotEqual(first_self_node["nearest_counterpart_ratio"], "")
            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("1. Does hotspot mismatch persist after comparing cell-vs-cell", summary)
            self.assertIn("7. If hotspot mismatch persists", summary)


if __name__ == "__main__":
    unittest.main()
