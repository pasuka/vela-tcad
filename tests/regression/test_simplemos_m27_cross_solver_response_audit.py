#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "cross_solver_response_audit")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m27_cross_solver_response_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM27EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.report = json.loads((
            ROOT / "m27_cross_solver_response_audit_report.json")
            .read_text(encoding="utf-8"))

    def test_shared_topology_and_raw_data_are_complete(self) -> None:
        execution = self.report["execution"]
        acceptance = self.report["acceptance"]
        self.assertEqual(execution["sentaurus_state_count"], 2)
        self.assertEqual(execution["common_silicon_node_count"], 942)
        self.assertEqual(execution["common_silicon_edge_count"], 2691)
        self.assertTrue(acceptance["all_checks_pass"])
        self.assertEqual(acceptance["maximum_bias_error_V"], 0.0)
        self.assertEqual(acceptance["nonfinite_fraction"], 0.0)

    def test_terminal_hfs_response_agrees_but_baseline_gap_remains(self) -> None:
        response = self.report["terminal_response"]
        self.assertLess(response["response_relative_disagreement"], 0.05)
        self.assertGreater(response["response_ratio_sentaurus_per_vela"], 0.95)
        self.assertLess(response["response_ratio_sentaurus_per_vela"], 1.05)
        self.assertLess(
            response["hfs_induced_gap_change_fraction_of_full_gap"], 0.01)

    def test_quasi_fermi_response_is_spatially_aligned(self) -> None:
        stats = self.report["cross_solver_response"]["responsive_node_phin"]
        self.assertGreater(stats["pearson"], 0.995)
        self.assertGreater(stats["sign_agreement_fraction"], 0.99)
        self.assertGreater(stats["origin_slope_sentaurus_per_vela"], 0.95)
        self.assertLess(stats["origin_slope_sentaurus_per_vela"], 1.10)
        self.assertLess(stats["p95_absolute_difference"], 1e-4)

    def test_mobility_response_is_aligned_on_responsive_active_edges(self) -> None:
        stats = self.report["cross_solver_response"]["active_edge_mobility"]
        self.assertGreater(stats["count"], 100)
        self.assertGreater(stats["pearson"], 0.95)
        self.assertLess(stats["p95_absolute_difference"], 0.05)

    def test_upstream_response_and_drain_pinning_match_m26_pathway(self) -> None:
        zones = self.report["zones"]
        upstream = zones["upstream_operator_source"]
        drain = zones["drain_first_layer"]
        self.assertGreater(upstream["vela_phin_response_abs_max_V"], 2e-3)
        self.assertGreater(upstream["sentaurus_phin_response_abs_max_V"], 2e-3)
        self.assertLess(drain["vela_phin_response_abs_max_V"], 1e-12)
        self.assertLess(drain["sentaurus_phin_response_abs_max_V"], 1e-12)

    def test_frozen_artifacts_match_hashes(self) -> None:
        for artifact in self.evidence["artifacts"].values():
            path = REPO / artifact["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), artifact["sha256"], path)
        for artifact in self.evidence["figures"]:
            path = REPO / artifact["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), artifact["sha256"], path)


if __name__ == "__main__":
    unittest.main()
