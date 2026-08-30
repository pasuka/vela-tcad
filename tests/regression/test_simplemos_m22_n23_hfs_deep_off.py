#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m22_n23_hfs_deep_off_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM22EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def test_report_closes_targeted_matrix(self) -> None:
        report = json.loads((ROOT / "m22_n23_hfs_deep_off_report.json").read_text())
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["state_count"], 8)
        self.assertEqual(report["frozen_formula_evaluations"], 16)

    def test_key_state_separates_direct_and_self_consistent_hfs(self) -> None:
        report = json.loads((ROOT / "m22_n23_hfs_deep_off_report.json").read_text())
        key = report["key_state"]
        findings = report["findings"]
        self.assertAlmostEqual(key["cross_solver_error_improvement_dex"],
                               0.019019967339194324, places=14)
        self.assertLess(abs(findings["key_full_state_frozen_hfs_effect_dex"]),
                        1.0e-5)
        self.assertGreater(findings["key_max_per_edge_sg_cancellation_condition"],
                           1.0e14)
        self.assertEqual(findings["key_max_frozen_drain_cut_mobility_effect_dex"],
                         0.0)

    def test_direct_bias_lattice_is_exact(self) -> None:
        with (ROOT / "m22_summary.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual([float(row["gate_voltage_V"]) for row in rows],
                         [0.0, 0.05, 0.1, 0.15])
        for row in rows:
            self.assertEqual(float(row["vela_full_current_A_per_um"]),
                             float(row["rerun_full_current_A_per_um"]))
            self.assertEqual(float(row["vela_no_hfs_current_A_per_um"]),
                             float(row["rerun_no_hfs_current_A_per_um"]))

    def test_portable_artifacts_and_figure_match_hashes(self) -> None:
        for artifact in self.evidence["artifacts"].values():
            path = REPO / artifact["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), artifact["sha256"], path)
        for figure in self.evidence["figures"]:
            path = REPO / figure["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), figure["sha256"], path)


if __name__ == "__main__":
    unittest.main()
