import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/measure_contact_qf_ablation"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m37_measure_contact_qf_ablation_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM37MeasureContactQfAblationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            ROOT / "m37_measure_contact_qf_ablation_report.json").read_text())

    def test_ablation_is_complete_and_defaults_are_unchanged(self) -> None:
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertFalse(self.report["execution"]["default_model_changed"])
        self.assertEqual(len(self.report["boundary_measure_ablation"]["variants"]), 3)

    def test_boundary_support_is_minor_part_of_measure_response(self) -> None:
        fraction = self.report["boundary_measure_ablation"][
            "boundary_scope_fraction_of_all_signed_log_response"]
        self.assertGreater(fraction, 0.0)
        self.assertLess(fraction, 0.05)

    def test_contact_cut_has_no_direct_hfs_mobility_response(self) -> None:
        contact = self.report["contact_sg_ablation"]
        self.assertEqual(contact["drain_cut_edge_count"], 7)
        self.assertEqual(contact["maximum_mobility_log10_response_dex"], 0.0)
        self.assertEqual(contact["signed_direct_electron_current_delta_A_per_um"], 0.0)

    def test_qf_relaxation_dominates_finite_hfs_response(self) -> None:
        feedback = self.report["quasi_fermi_feedback_ablation"]
        self.assertLess(feedback["direct_fraction_of_actual"], 1e-4)
        self.assertGreater(feedback["adjoint_relaxation_fraction_of_actual"], 0.9)
        self.assertLess(feedback["first_order_relative_error"], 0.1)
        self.assertGreater(feedback["electron_equation_signed_fraction"], 0.999)

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in evidence["source_hashes"].items():
            self.assertEqual(expected, sha256(REPO / relative))


if __name__ == "__main__":
    unittest.main()
