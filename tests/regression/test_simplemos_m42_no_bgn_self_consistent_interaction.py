import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "no_bgn_self_consistent_interaction")
REPORT = ROOT / "m42_no_bgn_self_consistent_interaction_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m42_no_bgn_self_consistent_interaction_evidence.json")
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m42_no_bgn_self_consistent_interaction_contract_v1.json")
M43_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m43_sg_kernel_consistency_evidence.json")
M44_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m44_qf_coordinate_consistency_evidence.json")
M45_REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
              / "post_qf_rebaseline/m45_post_qf_rebaseline_report.json")
M45_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m45_post_qf_rebaseline_evidence.json")
M46_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m46_full_matrix_requalification_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


class SimpleMosM42NoBgnSelfConsistentInteractionTest(unittest.TestCase):
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
        cls.m46_evidence = json.loads(
            M46_EVIDENCE.read_text(encoding="utf-8"))

    def test_matrix_is_complete_and_defaults_are_unchanged(self) -> None:
        self.assertEqual(self.report["status"], "failed")
        self.assertFalse(self.report["acceptance"]["all_checks_pass"])
        self.assertFalse(self.report["acceptance"][
            "legacy_cut_reported_mismatch_reproduced"])
        self.assertEqual(self.m45_report["status"], "complete")
        self.assertTrue(self.m45_report["findings"][
            "m42_historical_legacy_mismatch_hypothesis_superseded"])
        execution = self.report["execution"]
        self.assertEqual(execution["self_consistent_state_count"], 8)
        self.assertEqual(execution["frozen_operator_pair_count"], 32)
        self.assertEqual(execution["qf_substitution_count"], 14)
        self.assertFalse(execution["default_model_changed"])
        self.assertFalse(execution["new_sentaurus_execution"])

    def test_contact_reconstruction_is_exact_null_control(self) -> None:
        contact = rows(ROOT / "m42_contact_reconstruction_response.csv")
        self.assertEqual(len(contact), 4)
        for row in contact:
            for field in (
                    "operator_cut_log_shift_dex",
                    "reported_current_log_shift_dex",
                    "maximum_abs_psi_delta_V",
                    "maximum_abs_phin_delta_V",
                    "maximum_abs_phip_delta_V",
                    "maximum_abs_logn_delta_dex",
                    "maximum_abs_logp_delta_dex"):
                self.assertEqual(float(row[field]), 0.0)
        self.assertTrue(self.report["causal_interpretation"][
            "contact_reconstruction_is_null_for_simplemos"])

    def test_srh_is_feedback_and_electron_qf_is_the_state_path(self) -> None:
        finding = self.report["findings"]
        self.assertEqual(
            finding["maximum_srh_direct_frozen_cut_delta_A_per_um"], 0.0)
        self.assertGreater(
            finding["maximum_srh_recombination_operator_change_l2"], 1.9e-13)
        for fraction in finding[
                "electron_qf_density_fraction_of_srh_shift"].values():
            self.assertGreater(fraction, 0.999)
            self.assertLess(fraction, 1.001)
        interpretation = self.report["causal_interpretation"]
        self.assertTrue(
            interpretation["srh_terminal_effect_is_state_feedback_not_direct_operator"])
        self.assertTrue(
            interpretation["electron_qf_density_is_primary_srh_feedback_path"])

    def test_legacy_and_compensated_terminal_paths_are_consistent(self) -> None:
        finding = self.report["findings"]
        self.assertLess(finding["minimum_legacy_cut_reported_gap_dex"], 1e-12)
        self.assertLess(
            finding["maximum_compensated_cut_reported_gap_dex"], 1e-12)
        self.assertGreater(
            finding["legacy_to_compensated_residual_separation_min_ratio"],
            1e9)
        self.assertLess(
            finding["maximum_compensated_free_electron_residual_l2"], 2e-20)
        self.assertGreater(
            finding["maximum_transport_direct_frozen_cut_shift_abs_dex"],
            3.2)
        self.assertFalse(self.report["causal_interpretation"][
            "legacy_state_terminal_extraction_is_operator_inconsistent"])

    def test_reference_increment_is_required_in_deep_off_state(self) -> None:
        finding = self.report["findings"]
        self.assertGreater(
            finding["compensated_absolute_qf_roundtrip_loss_dex"], 2.65)
        self.assertTrue(self.report["causal_interpretation"][
            "qf_reference_increment_representation_is_required"])
        qf = {(row["mode"], row["variant"]): row for row in rows(
            ROOT / "m42_qf_state_substitution.csv")}
        baseline = abs(float(qf[("compensated", "srh_on_baseline")][
            "drain_cut_total_A_per_um"]))
        absolute_only = abs(float(qf[("compensated", "absolute_qf_only_roundtrip")][
            "drain_cut_total_A_per_um"]))
        self.assertGreater(baseline / absolute_only, 400.0)

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        self.assertEqual(self.evidence["status"], "failed")
        self.assertFalse(self.evidence["default_model_changed"])
        for artifact in self.evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in self.evidence["source_hashes"].items():
            current = sha256(REPO / relative)
            if expected == current:
                continue
            superseding = (self.m46_evidence, self.m45_evidence, self.m44_evidence,
                           self.m43_evidence)
            self.assertTrue(
                any(evidence["status"] == "frozen" and
                    evidence["source_hashes"].get(relative) == current
                    for evidence in superseding),
                f"changed M42 source is not frozen by M43/M44/M45/M46: {relative}")


if __name__ == "__main__":
    unittest.main()
