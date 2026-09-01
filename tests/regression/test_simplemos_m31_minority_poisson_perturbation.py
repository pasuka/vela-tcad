import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_source_hashes_current_or_m44,
)


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/minority_poisson_perturbation"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m31_minority_poisson_perturbation_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM31MinorityPoissonPerturbationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(
            (ROOT / "m31_minority_poisson_perturbation_report.json").read_text())
        cls.states = {row["variant"]: row
                      for row in cls.report["state_variants"]}
        cls.cross = {row["variant"]: row
                     for row in cls.report["cross_block"]}

    def test_probe_contract_is_complete(self) -> None:
        acceptance = self.report["acceptance"]
        self.assertTrue(acceptance["all_checks_pass"])
        self.assertEqual(acceptance["state_variant_count"], 8)
        self.assertEqual(acceptance["cross_block_variant_count"], 3)
        self.assertEqual(acceptance["common_silicon_node_count"], 942)
        self.assertEqual(acceptance["condition_estimates_skipped_count"], 3)

    def test_minority_hole_direct_current_response_is_small(self) -> None:
        hole = self.report["minority_hole"]
        self.assertEqual(hole["top20_exact_delta_current_A_per_um"], 0.0)
        self.assertLess(hole["full_current_shift_dex"], 5e-4)
        self.assertLess(abs(hole["full_exact_delta_current_A_per_um"]), 3e-22)

    def test_electron_state_controls_frozen_current(self) -> None:
        self.assertGreater(
            self.states["electron_full_replace"][
                "absolute_log10_current_shift_dex"], 4.19)
        self.assertAlmostEqual(
            self.states["electron_full_replace"][
                "absolute_log10_current_shift_dex"],
            self.states["qf_full_replace"][
                "absolute_log10_current_shift_dex"], places=6)

    def test_poisson_floor_is_spatially_concentrated_but_not_direct_current(self) -> None:
        floor = self.report["poisson_floor"]
        self.assertGreater(floor["top10_l2_share"], 0.996)
        self.assertEqual(floor["top10_step_current_shift_dex"], 0.0)
        self.assertEqual(floor["full_step_current_shift_dex"], 0.0)

    def test_cross_blocks_are_finite_and_condition_estimates_are_skipped(self) -> None:
        for row in self.cross.values():
            self.assertFalse(row["condition_estimates_computed"])
            self.assertLess(row["schur_relative_closure"], 1e-20)
            self.assertGreater(row["psi_qfp_product_l2"], 0.0)
        self.assertGreater(
            self.cross["electron_full_replace"]["psi_qfp_product_l2"],
            self.cross["hole_full_replace"]["psi_qfp_product_l2"])

    def test_audit_does_not_change_production_defaults(self) -> None:
        execution = self.report["execution"]
        self.assertTrue(execution["cpp_changed"])
        self.assertTrue(execution["diagnostic_only_cpp_changed"])
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"],
                             sha256(REPO / artifact["path"]))
        assert_source_hashes_current_or_m44(self, evidence)


if __name__ == "__main__":
    unittest.main()
