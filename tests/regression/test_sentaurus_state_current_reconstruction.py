#!/usr/bin/env python3
"""Regression tests for Sentaurus-state current reconstruction diagnostics."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
Q_C = 1.602176634e-19


class SentaurusStateCurrentReconstructionTest(unittest.TestCase):
    def test_reconstructs_edge_linear_qf_current_and_effective_mobility(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sent_state_current_") as tmp:
            root = Path(tmp)
            case = root / "sentaurus_-20v"
            fields = case / "fields"
            fields.mkdir(parents=True)
            out_dir = root / "out"

            self.write_csv(case / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0e6, 0.0],
                [2, 0.0, 1.0e6],
            ])
            self.write_csv(case / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                [0, 0, 1, 2, "R.Si", "Si"],
            ])
            self.write_field(fields / "eDensity_region0.csv", {0: 10.0, 1: 30.0, 2: 10.0})
            self.write_field(fields / "hDensity_region0.csv", {0: 20.0, 1: 20.0, 2: 20.0})
            self.write_field(fields / "eMobility_region0.csv", {0: 2500.0, 1: 2500.0, 2: 2500.0})
            self.write_field(fields / "hMobility_region0.csv", {0: 4000.0, 1: 4000.0, 2: 4000.0})
            self.write_field(fields / "eQuasiFermiPotential_region0.csv", {0: 1.0, 1: 5.0, 2: 1.0})
            self.write_field(fields / "hQuasiFermiPotential_region0.csv", {0: 2.0, 1: 2.0, 2: 7.0})
            # Edge 0-1 electron linear flux = (2500 cm2/Vs -> .25 m2/Vs)
            # * (20 cm-3 -> 2e7 m-3) * 4 V / 1 m = 2e7 m-2 s-1.
            # Current density below is 4x that value in A/cm2.
            e_current_a_cm2 = 4.0 * 2.0e7 * Q_C / 1.0e4
            h_current_a_cm2 = 3.0 * 4.0e7 * Q_C / 1.0e4
            self.write_field(fields / "eCurrentDensity_region0.csv", {0: e_current_a_cm2, 1: e_current_a_cm2, 2: 0.0})
            self.write_field(fields / "hCurrentDensity_region0.csv", {0: h_current_a_cm2, 1: 0.0, 2: h_current_a_cm2})
            self.write_field(fields / "ElectricField_region0.csv", {0: 1.0, 1: 2.0, 2: 3.0})

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sentaurus_state_current_reconstruction.py"),
                    "--sentaurus-root",
                    str(root),
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
            rows = self.read_csv(out_dir / "sentaurus_state_current_edges.csv")
            edge01 = next(row for row in rows if {row["node0"], row["node1"]} == {"0", "1"})
            self.assertAlmostEqual(float(edge01["electron_linear_qf_flux_m2_s"]), 2.0e7)
            self.assertAlmostEqual(float(edge01["electron_current_over_linear_qf"]), 4.0)
            self.assertAlmostEqual(float(edge01["electron_mu_eff_over_mobility"]), 4.0)

            summary = self.read_csv(out_dir / "sentaurus_state_current_summary.csv")
            self.assertEqual(len(summary), 1)
            self.assertGreater(float(summary[0]["electron_current_over_linear_qf_median"]), 0.0)
            payload = json.loads((out_dir / "sentaurus_state_current_summary.json").read_text())
            self.assertEqual(payload["schema"], "pn2d.sentaurus_state_current_reconstruction.v1")

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
