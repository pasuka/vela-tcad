#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/contact_boundary_audit"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m24_contact_boundary_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM24EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.report = json.loads((ROOT / "m24_contact_boundary_audit_report.json")
                                .read_text(encoding="utf-8"))

    def test_contact_boundary_is_pinned_and_unchanged(self) -> None:
        closure = self.report["closure"]
        findings = self.report["findings"]
        self.assertLess(closure["maximum_contact_qf_bias_error_V"], 1.0e-12)
        self.assertEqual(findings[
            "maximum_full_no_hfs_contact_physical_state_delta"], 0.0)
        self.assertLess(findings[
            "maximum_full_no_hfs_contact_qf_increment_delta_V"], 1.0e-20)

    def test_first_layer_carries_self_consistent_response(self) -> None:
        findings = self.report["findings"]
        self.assertAlmostEqual(findings[
            "self_consistent_cut_current_effect_dex"], -0.099826495218532,
            places=13)
        self.assertGreater(findings[
            "maximum_first_layer_qf_increment_delta_V"], 1.0e-17)
        self.assertGreater(findings[
            "maximum_first_layer_contact_column_share"], 0.4)
        self.assertLess(self.report["closure"][
            "maximum_continuity_term_absolute_closure"], 1.0e-18)

    def test_current_extraction_difference_is_bounded(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure[
            "maximum_sg_cut_current_relative_disagreement"], 1.0e-3)
        self.assertLess(closure[
            "maximum_sg_cut_current_absolute_disagreement_A_per_um"], 1.0e-18)

    def test_frozen_artifacts_match_hashes(self) -> None:
        for artifact in self.evidence["artifacts"].values():
            path = REPO / artifact["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), artifact["sha256"], path)
        for group in ("upstream_inputs", "figures"):
            for artifact in self.evidence[group]:
                path = REPO / artifact["path"]
                self.assertTrue(path.is_file(), path)
                self.assertEqual(sha256(path), artifact["sha256"], path)


if __name__ == "__main__":
    unittest.main()
