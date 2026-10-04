#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "qf_reference_transform_audit")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m25_qf_reference_transform_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM25EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.report = json.loads((
            ROOT / "m25_qf_reference_transform_audit_report.json")
            .read_text(encoding="utf-8"))

    def test_drain_cut_does_not_cross_reference_basin(self) -> None:
        topology = self.report["topology"]
        self.assertEqual(topology["drain_cut_edge_count"], 7)
        self.assertEqual(
            topology["drain_cut_reference_transition_edge_count"], 0)
        self.assertEqual(topology["reference_transition_edge_count"], 51)

    def test_reference_transform_and_sg_log_close(self) -> None:
        closure = self.report["closure"]
        self.assertTrue(closure["all_acceptance_checks_pass"])
        self.assertLess(closure["maximum_endpoint_transform_error_V"], 1e-16)
        self.assertLess(closure[
            "maximum_physical_drop_closure_error_V"], 1e-16)
        self.assertLess(closure[
            "maximum_sg_log_state_closure_error"], 1e-12)
        self.assertLess(closure[
            "maximum_factorized_flux_relative_error"], 1e-12)

    def test_reference_jump_is_required(self) -> None:
        findings = self.report["findings"]
        self.assertAlmostEqual(findings[
            "maximum_naive_increment_only_error_on_reference_transition_V"],
            0.05, places=15)
        self.assertGreater(findings[
            "maximum_reference_jump_thermal_voltages"], 1.9)
        self.assertLess(findings[
            "maximum_full_no_hfs_drain_cut_qf_drop_delta_V"], 1e-15)

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
