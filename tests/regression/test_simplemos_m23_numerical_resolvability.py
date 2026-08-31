#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/numerical_resolvability"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m23_numerical_resolvability_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM23EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.report = json.loads((ROOT / "m23_numerical_resolvability_report.json")
                                .read_text(encoding="utf-8"))

    def test_factorized_production_flux_is_resolved(self) -> None:
        findings = self.report["findings"]
        self.assertEqual(self.report["execution"]["edge_evaluation_count"], 36)
        self.assertTrue(findings[
            "all_edges_resolved_in_production_representation"])
        self.assertLess(findings[
            "maximum_production_vs_decimal_relative_error"], 1.0e-12)
        self.assertLess(findings[
            "maximum_one_ulp_flux_relative_envelope"], 1.0e-12)

    def test_raw_condition_is_not_production_error_bound(self) -> None:
        findings = self.report["findings"]
        self.assertGreater(findings["maximum_raw_subtractive_sg_condition"],
                           1.0e15)
        self.assertEqual(findings["absolute_export_collapsed_edge_count"], 10)
        self.assertLess(findings[
            "maximum_production_vs_decimal_relative_error"], 1.0e-12)

    def test_solver_path_noise_is_below_hfs_response(self) -> None:
        spread = self.report["findings"]["solver_path_spread"]
        self.assertLess(spread["full"]["spread_dex"], 1.0e-4)
        self.assertLess(spread["no_hfs"]["spread_dex"], 1.0e-7)
        self.assertGreater(0.09977223082251109 /
                           spread["full"]["spread_dex"], 1000.0)

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
