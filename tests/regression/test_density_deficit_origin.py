#!/usr/bin/env python3
"""Regression tests for PN2D density-deficit origin audit."""

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
VT = 0.025851999786435535


class DensityDeficitOriginAuditTest(unittest.TestCase):
    def test_splits_density_ratio_into_level_and_ni_eff_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_density_origin_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            overlap = root / "overlap.csv"
            sg_edges = root / "sg_edges.csv"
            state = root / "state.csv"
            manifest = root / "manifest.json"
            sent_case = root / "sentaurus" / "sentaurus_-20v"
            fields = sent_case / "fields"
            fields.mkdir(parents=True)
            out_dir = root / "out"

            level_ratio = math.exp(-2.0)
            ni_ratio = 0.5
            sent_n = 100.0
            vela_n = sent_n * level_ratio * ni_ratio
            hole_level_ratio = math.exp(-1.0)
            hole_ni_ratio = 2.0
            sent_p = 80.0
            vela_p = sent_p * hole_level_ratio * hole_ni_ratio

            mesh.write_text(json.dumps({
                "nodes": [{"id": 0, "x": 0.0, "y": 0.0}, {"id": 1, "x": 1.0e6, "y": 0.0}],
                "triangles": [],
            }), encoding="utf-8")
            self.write_csv(overlap, ["bias_V", "edge_id", "node0", "node1", "sent_high_field_p95"], [[-20.0, 0, 0, 1, "true"]])
            self.write_csv(sg_edges, ["bias_V", "edge_id", "node0", "node1", "edge_length_m"], [[-20.0, 0, 0, 1, 1.0]])
            self.write_csv(state, ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"], [
                [0, 0.0, 2.0 * VT, -VT, vela_n * 1.0e6, vela_p * 1.0e6],
                [1, 0.0, 2.0 * VT, -VT, vela_n * 1.0e6, vela_p * 1.0e6],
            ])
            manifest.write_text(json.dumps([{"bias_V": -20.0, "sg_edges_csv": str(sg_edges), "state_csv": str(state)}]), encoding="utf-8")
            self.write_csv(sent_case / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0], [1, 1.0e6, 0.0]])
            self.write_field(fields / "ElectrostaticPotential_region0.csv", {0: 0.0, 1: 0.0})
            self.write_field(fields / "eQuasiFermiPotential_region0.csv", {0: 0.0, 1: 0.0})
            self.write_field(fields / "hQuasiFermiPotential_region0.csv", {0: 0.0, 1: 0.0})
            self.write_field(fields / "eDensity_region0.csv", {0: sent_n, 1: sent_n})
            self.write_field(fields / "hDensity_region0.csv", {0: sent_p, 1: sent_p})

            result = subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "diagnose_pn2d_bv_density_deficit_origin.py"),
                "--overlap-csv", str(overlap),
                "--sg-edge-manifest", str(manifest),
                "--sentaurus-root", str(root / "sentaurus"),
                "--mesh", str(mesh),
                "--bias-min", "-20",
                "--bias-max", "-20",
                "--out-dir", str(out_dir),
            ], text=True, capture_output=True)

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self.read_csv(out_dir / "density_deficit_origin_edges.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["electron_level_density_ratio"]), level_ratio)
            self.assertAlmostEqual(float(row["electron_inferred_ni_eff_ratio"]), ni_ratio)
            self.assertAlmostEqual(float(row["hole_level_density_ratio"]), hole_level_ratio)
            self.assertAlmostEqual(float(row["hole_inferred_ni_eff_ratio"]), hole_ni_ratio)

            summary = self.read_csv(out_dir / "density_deficit_origin_summary.csv")
            self.assertAlmostEqual(float(summary[0]["electron_level_density_ratio_high_field_median"]), level_ratio)
            self.assertAlmostEqual(float(summary[0]["electron_inferred_ni_eff_ratio_high_field_median"]), ni_ratio)

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
