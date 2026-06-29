#!/usr/bin/env python3
"""Regression coverage for internal edge source versus Sentaurus comparison."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "compare_internal_edge_source_vs_sentaurus.py"
Q = 1.602176634e-19


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


class InternalEdgeSourceVsSentaurusTest(unittest.TestCase):
    def test_synthetic_edge_comparison_closes_ratios_and_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_edge_sentaurus_") as td:
            tmp = Path(td)
            internal = tmp / "internal.csv"
            topology = tmp / "sg_edges.csv"
            sentaurus = tmp / "sentaurus_multibias"
            out_csv = tmp / "edge_compare.csv"
            out_md = tmp / "edge_compare.md"

            area0 = 2.0e-8
            area1 = 3.0e-8
            g0 = (2.0 * 4.0 + 3.0 * 5.0) / Q
            g1 = (4.0 * 6.0 + 7.0 * 8.0) / Q
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
                    "Gava_reconstructed_from_used_terms",
                    "Gava_closure_relative_error",
                    "source_weight_or_volume_cm2_for_2D",
                    "qG_contribution_A_per_um",
                ],
                [
                    ["edge", 0, -20.0, 0.75, 0.0, 1.0e5, 1.0e5, 2.0, 3.0, 4.0, 5.0, 2.0 * 4.0 / Q, 3.0 * 5.0 / Q, g0, g0, 0.0, area0, Q * g0 * area0 / 1.0e4],
                    ["edge", 1, -20.0, 1.20, 0.0, 2.0e5, 2.0e5, 4.0, 7.0, 6.0, 8.0, 4.0 * 6.0 / Q, 7.0 * 8.0 / Q, g1, g1, 0.0, area1, Q * g1 * area1 / 1.0e4],
                ],
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
                [
                    [0, -20.0, 0, 0, 1, 0.25, 0.0, 1.25, 0.0, 1.0e-6, area0 / 1.0e4],
                    [0, -20.0, 1, 1, 2, 1.0, 0.0, 1.4, 0.0, 4.0e-7, area1 / 1.0e4],
                ],
            )
            fields = sentaurus / "sentaurus_-20v" / "fields"
            write_csv(fields / "eQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0], [1, 10.0], [2, 18.0]])
            write_csv(fields / "hQuasiFermiPotential_region0.csv", ["node_id", "component0"], [[0, 0.0], [1, 10.0], [2, 18.0]])
            write_csv(fields / "eAlphaAvalanche_region0.csv", ["node_id", "component0"], [[0, 2.0], [1, 2.0], [2, 6.0]])
            write_csv(fields / "hAlphaAvalanche_region0.csv", ["node_id", "component0"], [[0, 3.0], [1, 3.0], [2, 11.0]])
            write_csv(fields / "eCurrentDensity_region0.csv", ["node_id", "component0", "component1"], [[0, 4.0, 0.0], [1, 4.0, 0.0], [2, 8.0, 0.0]])
            write_csv(fields / "hCurrentDensity_region0.csv", ["node_id", "component0", "component1"], [[0, 5.0, 0.0], [1, 5.0, 0.0], [2, 11.0, 0.0]])
            write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [[0, g0], [1, g0], [2, g1]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--internal-audit",
                    str(internal),
                    "--self-edge-topology",
                    str(topology),
                    "--sentaurus-multibias-dir",
                    str(sentaurus),
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
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["edge_id"], "0")
            self.assertEqual(rows[0]["node0"], "0")
            self.assertAlmostEqual(float(rows[0]["Fn_ratio_self_over_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(rows[0]["alpha_n_ratio_self_over_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(rows[0]["Jn_ratio_self_over_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(rows[0]["Gava_ratio_self_over_sentaurus_method_A"]), 1.0)
            self.assertIn("qG_self_edge_A_per_um", rows[0])

            summary = out_md.read_text(encoding="utf-8")
            self.assertIn("self internal edge qG explains solver qG: yes", summary)
            self.assertIn("self exported node current density used: no", summary)
            self.assertIn("max Gava edge distance at -20V", summary)
            self.assertIn("Sentaurus nodal averaging can explain Sentaurus nodal Gava", summary)


if __name__ == "__main__":
    unittest.main()
