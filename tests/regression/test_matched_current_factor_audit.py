#!/usr/bin/env python3
"""Regression tests for matched Vela/Sentaurus current-support factor audit."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class MatchedCurrentFactorAuditTest(unittest.TestCase):
    def test_decomposes_density_mobility_and_qf_gradient_ratios(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_matched_factor_audit_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            overlap = root / "overlap.csv"
            sg_edges = root / "sg_edges.csv"
            state = root / "last_state.csv"
            manifest = root / "manifest.json"
            sent_case = root / "sentaurus" / "sentaurus_-20v"
            sent_fields = sent_case / "fields"
            sent_fields.mkdir(parents=True)
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [{"id": 0, "x": 0.0, "y": 0.0}, {"id": 1, "x": 1.0e6, "y": 0.0}],
                "triangles": [],
            }), encoding="utf-8")
            self.write_csv(overlap, [
                "bias_V", "edge_id", "node0", "node1", "sent_high_field_p95", "edge_area_proxy_m2",
                "sent_e_particle_flux_m2_s", "sent_h_particle_flux_m2_s",
            ], [[-20.0, 0, 0, 1, "true", 2.0, 100.0, 80.0]])
            self.write_csv(sg_edges, [
                "bias_V", "edge_id", "node0", "node1", "edge_length_m", "edge_area_proxy_m2",
                "electron_mobility_m2_V_s", "hole_mobility_m2_V_s",
                "electron_alpha_m_inv", "hole_alpha_m_inv",
                "electron_raw_signed_flux_proxy", "hole_raw_signed_flux_proxy",
            ], [[-20.0, 0, 0, 1, 1.0, 2.0, 0.2, 0.3, 5.0, 7.0, 10.0, 12.0]])
            self.write_csv(state, ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"], [
                [0, 0.0, 0.0, 1.0, 20.0e6, 60.0e6],
                [1, 0.0, 4.0, 7.0, 20.0e6, 60.0e6],
            ])
            manifest.write_text(json.dumps([{
                "bias_V": -20.0,
                "sg_edges_csv": str(sg_edges),
                "state_csv": str(state),
            }]), encoding="utf-8")

            self.write_csv(sent_case / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0], [1, 1.0e6, 0.0]])
            self.write_field(sent_fields / "eDensity_region0.csv", {0: 10.0, 1: 10.0})
            self.write_field(sent_fields / "hDensity_region0.csv", {0: 30.0, 1: 30.0})
            self.write_field(sent_fields / "eMobility_region0.csv", {0: 1000.0, 1: 1000.0})
            self.write_field(sent_fields / "hMobility_region0.csv", {0: 2000.0, 1: 2000.0})
            self.write_field(sent_fields / "eQuasiFermiPotential_region0.csv", {0: 0.0, 1: 2.0})
            self.write_field(sent_fields / "hQuasiFermiPotential_region0.csv", {0: 1.0, 1: 4.0})

            result = subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "diagnose_pn2d_bv_matched_current_factor_audit.py"),
                "--overlap-csv", str(overlap),
                "--sg-edge-manifest", str(manifest),
                "--sentaurus-root", str(root / "sentaurus"),
                "--mesh", str(mesh),
                "--bias-min", "-20",
                "--bias-max", "-20",
                "--out-dir", str(out_dir),
            ], text=True, capture_output=True)

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self.read_csv(out_dir / "matched_current_factor_edges.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["electron_density_ratio_vela_over_sentaurus"]), 2.0)
            self.assertAlmostEqual(float(row["electron_mobility_ratio_vela_over_sentaurus"]), 2.0)
            self.assertAlmostEqual(float(row["electron_qf_grad_ratio_vela_over_sentaurus"]), 2.0)
            self.assertAlmostEqual(float(row["electron_linear_support_ratio_vela_over_sentaurus"]), 8.0)
            self.assertAlmostEqual(float(row["hole_linear_support_ratio_vela_over_sentaurus"]), 6.0)

            summary = self.read_csv(out_dir / "matched_current_factor_summary.csv")
            self.assertEqual(len(summary), 1)
            self.assertAlmostEqual(float(summary[0]["electron_linear_support_ratio_median"]), 8.0)
            self.assertAlmostEqual(float(summary[0]["electron_density_ratio_high_field_median"]), 2.0)

    @staticmethod
    def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    @classmethod
    def write_field(cls, path: Path, values: dict[int, float]) -> None:
        cls.write_csv(path, ["node_id", "component0"], [[node, value] for node, value in values.items()])

    @staticmethod
    def read_csv(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))


if __name__ == "__main__":
    unittest.main()
