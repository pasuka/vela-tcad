import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/electron_qf_scaled_perturbation"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m32_electron_qf_scaled_perturbation_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM32ElectronQfScaledPerturbationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            ROOT / "m32_electron_qf_scaled_perturbation_report.json").read_text())
        cls.by_fraction = {row["fraction"]: row
                           for row in cls.report["fraction_response"]}

    def test_probe_matrix_is_complete(self) -> None:
        acceptance = self.report["acceptance"]
        self.assertTrue(acceptance["all_checks_pass"])
        self.assertEqual(acceptance["fraction_count"], 6)
        self.assertEqual(acceptance["functional_probe_count"], 6)
        self.assertEqual(acceptance["sg_probe_count"], 6)
        self.assertEqual(acceptance["cross_block_probe_count"], 6)
        self.assertEqual(acceptance["drain_cut_edge_count"], 7)
        self.assertEqual(acceptance["condition_estimates_skipped_count"], 6)

    def test_no_tested_nonzero_fraction_is_in_tangent_validity_range(self) -> None:
        self.assertIsNone(
            self.report["linearization"]["largest_tested_valid_fraction"])
        self.assertEqual(self.by_fraction[0.01]["exact_delta_current_A_per_um"], 0.0)
        self.assertLess(self.by_fraction[0.03]["exact_delta_current_A_per_um"], 0.0)
        self.assertGreater(self.by_fraction[0.1]["exact_delta_current_A_per_um"], 0.0)
        for fraction in (0.01, 0.03, 0.1, 0.3, 1.0):
            self.assertGreater(
                self.by_fraction[fraction]["tangent_relative_error"], 0.1)

    def test_terminal_functional_closes_to_sg_drain_cut(self) -> None:
        for row in self.by_fraction.values():
            self.assertLess(abs(row["functional_minus_cut_A_per_um"]), 1e-27)

    def test_baseline_subtracted_cross_block_response_scales(self) -> None:
        cross = {row["fraction"]: row
                 for row in self.report["cross_block_scaling"]}
        full = cross[1.0]["induced_full_raw_delta_psi_l2_V"]
        for fraction in (0.01, 0.03, 0.1, 0.3):
            ratio = cross[fraction]["induced_full_raw_delta_psi_l2_V"] / full
            self.assertAlmostEqual(ratio, fraction, delta=5e-4)
        self.assertLess(
            self.report["acceptance"]["maximum_schur_relative_closure"], 1e-18)

    def test_full_fraction_drain_cut_is_cancellation_sensitive(self) -> None:
        self.assertGreater(self.by_fraction[1.0]["drain_cut_condition"], 11.5)
        top_edges = self.report["top_transport_edges_at_full_fraction"]
        self.assertTrue(any(row["is_drain_cut"] for row in top_edges[:5]))

    def test_audit_preserves_production_defaults(self) -> None:
        execution = self.report["execution"]
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["cpp_changed"])
        self.assertFalse(execution["default_model_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in evidence["source_hashes"].items():
            self.assertEqual(expected, sha256(REPO / relative))


if __name__ == "__main__":
    unittest.main()
