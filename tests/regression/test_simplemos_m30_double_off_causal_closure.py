import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_historical_source_provenance,
)


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/double_off_causal_closure"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m30_double_off_causal_closure_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM30DoubleOffCausalClosureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(
            (ROOT / "m30_double_off_causal_closure_report.json").read_text())
        cls.double_off = cls.report["double_off_decomposition"]
        cls.interactions = {
            row["field"]: row for row in cls.report["nodal_interaction_summary"]}

    def test_matrix_and_probe_contract_is_complete(self) -> None:
        acceptance = self.report["acceptance"]
        self.assertTrue(acceptance["all_checks_pass"])
        self.assertEqual(acceptance["state_count"], 8)
        self.assertEqual(acceptance["state_operator_pair_count"], 32)
        self.assertEqual(acceptance["probe_run_count"], 96)
        self.assertEqual(acceptance["common_silicon_node_count"], 942)
        self.assertLessEqual(acceptance["maximum_decomposition_closure_dex"], 1e-12)

    def test_double_off_gap_has_exact_two_layer_decomposition(self) -> None:
        self.assertAlmostEqual(self.double_off["native_gap_dex"],
                               1.3212392857887345, places=14)
        self.assertAlmostEqual(
            self.double_off["common_operator_state_contribution_dex"],
            0.8825886389719867, places=14)
        self.assertAlmostEqual(
            self.double_off["sentaurus_state_vela_operator_contribution_dex"],
            0.43844371351571293, places=14)
        self.assertEqual(self.double_off["decomposition_closure_dex"], 0.0)

    def test_deep_off_contact_cut_conditioning_is_reported(self) -> None:
        self.assertEqual(self.double_off["vela_drain_condition"], 1.0)
        self.assertGreater(self.double_off["sentaurus_state_drain_condition"], 11.7)
        self.assertGreater(self.double_off["sentaurus_state_source_condition"], 4.2)

    def test_nodal_interaction_disagreement_is_carrier_specific(self) -> None:
        self.assertGreater(self.interactions["phin"]["pearson"], 0.98)
        self.assertLess(self.interactions["phip"]["pearson"], 0.56)
        self.assertGreater(self.interactions["phip"]["p95_absolute_difference"], 1.5e-3)

    def test_strict_reclosure_controls_do_not_claim_a_new_solution(self) -> None:
        controls = self.report["solver_controls"]
        self.assertEqual(len(controls), 6)
        self.assertFalse(any(row["converged"] for row in controls))
        self.assertTrue(all(row["failure_reason"] == "line_search_non_decrease"
                            for row in controls))
        self.assertEqual(sum(row["same_bias_reclosure"] for row in controls), 4)
        self.assertEqual(sum(row["continuation_path_repeat"] for row in controls), 2)

    def test_audit_does_not_change_production_defaults(self) -> None:
        execution = self.report["execution"]
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["cpp_changed"])
        self.assertFalse(execution["default_model_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"],
                             sha256(REPO / artifact["path"]))
        for relative, expected in evidence["source_hashes"].items():
            self.assertTrue((REPO / relative).is_file())
        assert_historical_source_provenance(self, EVIDENCE)


if __name__ == "__main__":
    unittest.main()
