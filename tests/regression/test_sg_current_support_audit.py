#!/usr/bin/env python3
"""Regression tests for PN2D SG current-support audit diagnostics."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class SgCurrentSupportAuditTest(unittest.TestCase):
    def test_audit_splits_electron_hole_support_and_qf_flux_factors(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sg_support_audit_") as tmp:
            root = Path(tmp)
            overlap = root / "high_field_overlap_edges.csv"
            sg_edges = root / "sg_avalanche_edges.csv"
            manifest = root / "manifest.json"
            mesh = root / "mesh.json"
            state = root / "state.csv"
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
                ],
                [[-20.0, 0, 0, 1, "true", "false", 61.0, 2.0, 11.0, 13.0]],
            )
            self.write_csv(
                sg_edges,
                [
                    "bias_V",
                    "edge_id",
                    "node0",
                    "node1",
                    "edge_length_m",
                    "edge_area_proxy_m2",
                    "electron_alpha_m_inv",
                    "hole_alpha_m_inv",
                    "electron_mobility_m2_V_s",
                    "hole_mobility_m2_V_s",
                    "electron_raw_signed_flux_proxy",
                    "hole_raw_signed_flux_proxy",
                    "electron_flux_proxy",
                    "hole_flux_proxy",
                ],
                [[-20.0, 0, 0, 1, 2.0, 2.0, 2.0, 3.0, 0.25, 0.4, -5.0, 7.0, 5.0, 7.0]],
            )
            manifest.write_text(
                json.dumps([{"bias_V": -20.0, "sg_edges_csv": str(sg_edges), "state_csv": str(state)}]),
                encoding="utf-8",
            )
            mesh.write_text(
                json.dumps({
                    "nodes": [{"id": 0, "x": 0.0, "y": 0.0}, {"id": 1, "x": 2.0e6, "y": 0.0}],
                    "triangles": [],
                }),
                encoding="utf-8",
            )
            self.write_csv(
                state,
                ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"],
                [[0, 0.0, 1.0, 2.0, 10.0, 20.0], [1, 0.5, 5.0, 7.0, 30.0, 40.0]],
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_current_support_audit.py"),
                    "--overlap-csv",
                    str(overlap),
                    "--sg-edge-manifest",
                    str(manifest),
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
            rows = self.read_csv(out_dir / "sg_current_support_edge_audit.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["electron_flux_over_sentaurus"]), 5.0 / 11.0)
            self.assertAlmostEqual(float(row["hole_flux_over_sentaurus"]), 7.0 / 13.0)
            self.assertAlmostEqual(float(row["electron_linear_qf_flux_m2_s"]), 10.0)
            self.assertAlmostEqual(float(row["hole_linear_qf_flux_m2_s"]), 30.0)
            self.assertAlmostEqual(float(row["electron_sg_over_linear_qf"]), 0.5)
            self.assertAlmostEqual(float(row["hole_sg_over_linear_qf"]), 7.0 / 30.0)

            summary = self.read_csv(out_dir / "sg_current_support_summary.csv")
            self.assertEqual(len(summary), 1)
            item = summary[0]
            self.assertAlmostEqual(float(item["vela_e_source_integral"]), 20.0)
            self.assertAlmostEqual(float(item["vela_h_source_integral"]), 42.0)
            self.assertAlmostEqual(float(item["sent_e_source_integral"]), 44.0)
            self.assertAlmostEqual(float(item["sent_h_source_integral"]), 78.0)
            self.assertAlmostEqual(float(item["electron_missing_fraction"]), 24.0 / 60.0)
            self.assertAlmostEqual(float(item["hole_missing_fraction"]), 36.0 / 60.0)

    def test_manifest_config_without_materials_file_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sg_support_missing_materials_") as tmp:
            root = Path(tmp)
            overlap = root / "high_field_overlap_edges.csv"
            sg_edges = root / "sg_avalanche_edges.csv"
            mesh = root / "mesh.json"
            state = root / "state.csv"
            config = root / "simulation.json"
            manifest = root / "manifest.json"
            out_dir = root / "out"

            self.write_csv(
                overlap,
                [
                    "bias_V", "edge_id", "node0", "node1", "sent_high_field_p95",
                    "sent_high_h_current_p95", "sent_impact_m3_s", "edge_area_proxy_m2",
                    "sent_e_particle_flux_m2_s", "sent_h_particle_flux_m2_s",
                ],
                [[-20.0, 0, 0, 1, "true", "false", 1.0, 1.0, 1.0, 1.0]],
            )
            self.write_csv(
                sg_edges,
                [
                    "bias_V", "edge_id", "node0", "node1", "edge_length_m",
                    "edge_area_proxy_m2", "electron_alpha_m_inv", "hole_alpha_m_inv",
                    "electron_mobility_m2_V_s", "hole_mobility_m2_V_s",
                    "electron_raw_signed_flux_proxy", "hole_raw_signed_flux_proxy",
                    "electron_flux_proxy", "hole_flux_proxy",
                ],
                [[-20.0, 0, 0, 1, 1.0, 1.0, 1.0, 1.0, 0.25, 0.4, 1.0, 1.0, 1.0, 1.0]],
            )
            mesh.write_text(
                json.dumps({
                    "nodes": [{"id": 0, "x": 0.0, "y": 0.0}, {"id": 1, "x": 1.0e6, "y": 0.0}],
                    "triangles": [],
                }),
                encoding="utf-8",
            )
            self.write_csv(
                state,
                ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"],
                [[0, 0.0, 0.0, 0.0, 1.0, 1.0], [1, 0.0, 0.0, 0.0, 1.0, 1.0]],
            )
            config.write_text(json.dumps({"sweep": {"write_state_file": str(state)}}), encoding="utf-8")
            manifest.write_text(
                json.dumps([{"bias_V": -20.0, "sg_edges_csv": str(sg_edges), "config": str(config)}]),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_current_support_audit.py"),
                    "--overlap-csv", str(overlap),
                    "--sg-edge-manifest", str(manifest),
                    "--mesh", str(mesh),
                    "--bias-min", "-20",
                    "--bias-max", "-20",
                    "--out-dir", str(out_dir),
                ],
                text=True,
                capture_output=True,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("materials_file", result.stderr)
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
