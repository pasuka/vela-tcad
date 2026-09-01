import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "equal_ni_flux_ablation/m41_equal_ni_flux_ablation_report.json")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m41_equal_ni_flux_ablation_evidence.json")
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m41_equal_ni_flux_ablation_contract_v1.json")
M43_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m43_sg_kernel_consistency_evidence.json")
M44_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m44_qf_coordinate_consistency_evidence.json")
M45_REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
              / "post_qf_rebaseline/m45_post_qf_rebaseline_report.json")
M45_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m45_post_qf_rebaseline_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM41EqualNiFluxAblationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        cls.contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        cls.m43_evidence = json.loads(
            M43_EVIDENCE.read_text(encoding="utf-8"))
        cls.m44_evidence = json.loads(
            M44_EVIDENCE.read_text(encoding="utf-8"))
        cls.m45_report = json.loads(M45_REPORT.read_text(encoding="utf-8"))
        cls.m45_evidence = json.loads(
            M45_EVIDENCE.read_text(encoding="utf-8"))

    def test_historical_contract_is_explicitly_superseded(self) -> None:
        self.assertEqual(self.report["status"], "failed")
        self.assertFalse(self.report["acceptance"]["all_checks_pass"])
        self.assertFalse(self.report["acceptance"]["legacy_m40_reproduced"])
        self.assertEqual(self.m45_report["status"], "complete")
        self.assertTrue(self.m45_report["findings"][
            "m41_historical_contract_superseded"])
        self.assertTrue(self.report["acceptance"]["all_workflows_converged"])
        self.assertFalse(self.report["execution"]["default_model_changed"])
        self.assertFalse(self.evidence["default_model_changed"])
        self.assertEqual(
            self.report["numerical_contract"]["legacy"],
            "legacy_factor_difference")
        self.assertEqual(
            self.report["numerical_contract"]["compensated"],
            "compensated_log_expm1")

    def test_frozen_cut_exposes_m40_mixed_numerical_semantics(self) -> None:
        rows = self.report["frozen_cut_reconciliation"]
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertLess(row["compensated_cut_vs_reported_dex"], 0.002)
            self.assertGreater(row["legacy_cut_vs_reported_dex"], 2.0)

    def test_consistent_compensation_improves_both_srh_branches(self) -> None:
        rows = {row["cell"]: row
                for row in self.report["self_consistent_response"]}
        srh_off = rows["bgn_off_srh_off"]
        srh_on = rows["bgn_off_srh_on"]
        self.assertGreater(srh_off["gap_improvement_dex"], 1.17)
        self.assertGreater(srh_off["compensated_gap_dex"], 0.73)
        self.assertGreater(srh_on["gap_improvement_dex"], 1.37)
        self.assertGreater(srh_on["compensated_gap_dex"], 0.64)

        result = self.report["causal_result"]
        self.assertFalse(result["hypothesis_supported"])
        self.assertFalse(result["m40_no_bgn_factorial_clean"])
        self.assertFalse(
            result["m40_srh_on_close_parity_survives_consistent_flux"])
        self.assertTrue(result["srh_off_anomaly_partially_closed"])

    def test_contract_preserves_legacy_default(self) -> None:
        self.assertEqual(
            self.contract["evaluation_modes"]["legacy"],
            "legacy_factor_difference")
        self.assertFalse(
            self.contract["fixed_invariants"]["production_default_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        self.assertEqual(self.evidence["status"], "failed")
        for artifact in self.evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in self.evidence["source_hashes"].items():
            current = sha256(REPO / relative)
            if expected == current:
                continue
            superseding = (self.m45_evidence, self.m44_evidence,
                           self.m43_evidence)
            self.assertTrue(
                any(evidence["status"] == "frozen" and
                    evidence["source_hashes"].get(relative) == current
                    for evidence in superseding),
                f"changed M41 source is not frozen by M43/M44/M45: {relative}")


if __name__ == "__main__":
    unittest.main()
