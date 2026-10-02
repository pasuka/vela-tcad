#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "first_layer_feedback_audit")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m26_first_layer_feedback_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM26EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.report = json.loads((
            ROOT / "m26_first_layer_feedback_audit_report.json")
            .read_text(encoding="utf-8"))

    def test_current_topology_and_complete_one_ring_are_frozen(self) -> None:
        execution = self.report["execution"]
        self.assertEqual(
            execution["first_layer_nodes"], [1089, 1090, 1092, 1099, 1100])
        self.assertEqual(execution["first_layer_node_count"], 5)
        self.assertEqual(execution["perturbation_column_count"], 12)
        self.assertEqual(execution["central_perturbation_state_count"], 48)
        self.assertEqual(execution["read_only_probe_count"], 150)

    def test_adjoint_fd_cut_and_term_closure_pass(self) -> None:
        closure = self.report["closure"]
        self.assertTrue(closure["all_acceptance_checks_pass"])
        self.assertLess(closure["maximum_adjoint_relative_residual"], 1e-8)
        self.assertLess(
            closure["maximum_direct_gradient_relative_disagreement"], 1e-4)
        self.assertLess(
            closure["maximum_functional_cut_gradient_relative_disagreement"],
            1e-8)
        self.assertLess(closure["maximum_fd_term_closure_relative"], 1e-8)
        self.assertLess(closure["maximum_first_order_relative_error"], 0.1)
        self.assertLess(closure["maximum_cut_prediction_relative_error"], 1e-3)

    def test_relaxation_dominates_frozen_operator_response(self) -> None:
        findings = self.report["findings"]
        self.assertLess(
            findings["turn_on_direct_frozen_fraction_of_actual"], 1e-3)
        self.assertGreater(
            findings["turn_on_adjoint_relaxation_fraction_of_actual"], 0.8)
        self.assertLess(
            findings["turn_on_adjoint_relaxation_fraction_of_actual"], 1.2)
        self.assertGreater(
            findings["turn_off_adjoint_relaxation_fraction_of_actual"], 0.8)
        self.assertLess(
            findings["turn_off_adjoint_relaxation_fraction_of_actual"], 1.2)
        self.assertLess(
            findings["turn_on_first_layer_operator_absolute_fraction"], 1e-6)

    def test_first_layer_is_terminal_conduit_not_operator_origin(self) -> None:
        fractions = self.report["findings"][
            "first_layer_cut_current_fraction_by_node_no_hfs_linearization"]
        self.assertGreater(fractions["1092"], 0.8)
        self.assertGreater(fractions["1100"], 0.17)
        self.assertAlmostEqual(fractions["1090"], 0.0, places=15)
        self.assertLess(
            self.report["closure"]["maximum_cut_prediction_relative_error"],
            1e-3)

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
