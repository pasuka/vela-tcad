import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
REPORT = ROOT / "post_qf_rebaseline/m45_post_qf_rebaseline_report.json"
EVIDENCE = ROOT / "simplemos_m45_post_qf_rebaseline_evidence.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM45PostQfRebaselineTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def test_rebaseline_is_complete_without_default_changes(self) -> None:
        self.assertEqual(self.report["status"], "complete")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertEqual(self.evidence["status"], "frozen")
        self.assertFalse(self.report["execution"]["default_model_changed"])
        self.assertFalse(self.report["execution"]["new_sentaurus_execution"])
        self.assertFalse(self.report["execution"]["historical_contracts_rewritten"])

    def test_stage_impacts_are_explicit(self) -> None:
        stages = {row["stage"]: row for row in self.report["stage_impact"]}
        self.assertEqual(set(stages), {"M31", "M32", "M37", "M41", "M42"})
        self.assertEqual(stages["M31"]["rebaseline_status"], "complete")
        self.assertEqual(stages["M32"]["rebaseline_status"], "complete")
        self.assertEqual(stages["M37"]["impact"], "bitwise unchanged control")
        self.assertEqual(stages["M41"]["rebaseline_status"], "failed")
        self.assertEqual(stages["M42"]["rebaseline_status"], "failed")

    def test_operator_and_qf_closure_are_rebaselined(self) -> None:
        finding = self.report["findings"]
        self.assertLess(
            finding["m42_maximum_compensated_cut_reported_gap_dex"], 1e-12)
        self.assertLess(
            finding["m42_minimum_legacy_cut_reported_gap_dex"], 1e-12)
        self.assertLess(
            finding["m42_maximum_compensated_free_electron_residual_l2"],
            2e-20)
        self.assertGreater(
            finding["m42_legacy_to_compensated_residual_ratio"], 1e9)
        self.assertTrue(finding["legacy_operator_consistency_restored_by_m43"])
        self.assertTrue(
            finding["compensated_reference_replay_closure_restored_by_m44"])

    def test_retained_and_superseded_claims_are_separated(self) -> None:
        finding = self.report["findings"]
        self.assertFalse(finding["m31_qualitative_conclusion_changed"])
        self.assertFalse(finding["m32_qualitative_conclusion_changed"])
        self.assertTrue(finding["m37_bitwise_unchanged"])
        self.assertTrue(finding["m41_historical_contract_superseded"])
        self.assertTrue(
            finding["m42_historical_legacy_mismatch_hypothesis_superseded"])
        self.assertTrue(finding["srh_feedback_conclusion_retained"])
        self.assertTrue(finding["contact_reconstruction_null_retained"])
        self.assertTrue(finding["qf_reference_increment_requirement_retained"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        for artifact in self.evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(expected, sha256(REPO / relative))


if __name__ == "__main__":
    unittest.main()
