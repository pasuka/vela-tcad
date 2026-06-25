#!/usr/bin/env python3
"""Regression tests for PN2D current-support replay diagnostics."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class CurrentSupportReplayTest(unittest.TestCase):
    def test_replay_keeps_alpha_fixed_while_swapping_current_support(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_current_support_replay_") as tmp:
            root = Path(tmp)
            overlap = root / "high_field_overlap_edges.csv"
            sg_edges = root / "sg_avalanche_edges.csv"
            manifest = root / "manifest.json"
            mesh = root / "mesh.json"
            sentaurus_root = root / "sentaurus"
            fields = sentaurus_root / "sentaurus_-20v" / "fields"
            fields.mkdir(parents=True)
            out_dir = root / "out"

            self.write_csv(
                overlap,
                [
                    "bias_V",
                    "edge_id",
                    "node0",
                    "node1",
                    "sent_high_field_p95",
                    "sent_high_h_current_p95",
                    "sent_impact_m3_s",
                    "edge_area_proxy_m2",
                    "sent_e_particle_flux_m2_s",
                    "sent_h_particle_flux_m2_s",
                    "vela_e_flux_proxy_m2_s",
                    "vela_h_flux_proxy_m2_s",
                ],
                [
                    [-20.0, 0, 0, 1, "true", "false", 61.0, 1.0, 11.0, 13.0, 5.0, 7.0],
                    [-20.0, 1, 1, 2, "false", "false", 61.0, 1.0, 11.0, 13.0, 5.0, 7.0],
                    [-20.0, 2, 2, 0, "false", "false", 61.0, 1.0, 11.0, 13.0, 5.0, 7.0],
                ],
            )

            electron_edge1 = 1.0 / math.sqrt(2.0)
            hole_edge1 = 12.0 / math.sqrt(2.0)
            self.write_csv(
                sg_edges,
                [
                    "bias_V",
                    "edge_id",
                    "node0",
                    "node1",
                    "electron_alpha_m_inv",
                    "hole_alpha_m_inv",
                    "electron_flux_proxy",
                    "hole_flux_proxy",
                    "electron_raw_flux_proxy",
                    "hole_raw_flux_proxy",
                    "electron_raw_signed_flux_proxy",
                    "hole_raw_signed_flux_proxy",
                ],
                [
                    [-20.0, 0, 0, 1, 2.0, 3.0, 3.0, 0.0, 3.0, 0.0, 3.0, 0.0],
                    [-20.0, 1, 1, 2, 2.0, 3.0, abs(electron_edge1), abs(hole_edge1), abs(electron_edge1), abs(hole_edge1), electron_edge1, hole_edge1],
                    [-20.0, 2, 2, 0, 2.0, 3.0, 4.0, 12.0, 4.0, 12.0, -4.0, -12.0],
                ],
            )
            manifest.write_text(
                json.dumps([{"bias_V": -20.0, "sg_edges_csv": str(sg_edges)}]),
                encoding="utf-8",
            )
            mesh.write_text(
                json.dumps({
                    "nodes": [
                        {"id": 0, "x": 0.0, "y": 0.0},
                        {"id": 1, "x": 1.0, "y": 0.0},
                        {"id": 2, "x": 0.0, "y": 1.0},
                    ],
                    "triangles": [{"id": 0, "node_ids": [0, 1, 2]}],
                }),
                encoding="utf-8",
            )
            e_current = 11.0 * 1.602176634e-19 / 1.0e4
            h_current = 13.0 * 1.602176634e-19 / 1.0e4
            self.write_csv(
                fields / "eCurrentDensity_region0.csv",
                ["node_id", "component0"],
                [[0, e_current], [1, e_current], [2, e_current]],
            )
            self.write_csv(
                fields / "hCurrentDensity_region0.csv",
                ["node_id", "component0"],
                [[0, h_current], [1, h_current], [2, h_current]],
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_current_support_replay.py"),
                    "--overlap-csv",
                    str(overlap),
                    "--sg-edge-manifest",
                    str(manifest),
                    "--sentaurus-root",
                    str(sentaurus_root),
                    "--mesh",
                    str(mesh),
                    "--bias-min",
                    "-20",
                    "--bias-max",
                    "-20",
                    "--out-dir",
                    str(out_dir),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self.read_csv(out_dir / "current_support_replay_edges.csv")
            edge0 = next(row for row in rows if row["edge_id"] == "0")
            self.assertAlmostEqual(float(edge0["genius_edge_gradqf_density_m3_s"]), 6.0)
            self.assertAlmostEqual(float(edge0["charon_hcurl_centroid_density_m3_s"]), 46.0)
            self.assertAlmostEqual(float(edge0["sentaurus_node_current_injected_density_m3_s"]), 61.0)
            self.assertAlmostEqual(float(edge0["sentaurus_edge_current_injected_density_m3_s"]), 61.0)
            self.assertAlmostEqual(float(edge0["sentaurus_edge_current_injected_over_sentaurus"]), 1.0)

            summary = json.loads(
                (out_dir / "current_support_replay_summary.json").read_text(encoding="utf-8")
            )
            high_field = summary["by_bias"]["-20"]["sentaurus_p95_high_field"]
            self.assertEqual(high_field["count"], 1)
            self.assertAlmostEqual(high_field["charon_hcurl_centroid_area_weighted_ratio"], 46.0 / 61.0)
            self.assertAlmostEqual(high_field["sentaurus_node_current_injected_area_weighted_ratio"], 1.0)
            self.assertAlmostEqual(high_field["sentaurus_edge_current_injected_area_weighted_ratio"], 1.0)

            summary_rows = self.read_csv(out_dir / "current_support_replay_summary.csv")
            self.assertEqual(len(summary_rows), 1)
            self.assertAlmostEqual(
                float(summary_rows[0]["sentaurus_node_current_injected_over_sentaurus"]), 1.0
            )

    @staticmethod
    def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    @staticmethod
    def read_csv(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))


if __name__ == "__main__":
    unittest.main()
