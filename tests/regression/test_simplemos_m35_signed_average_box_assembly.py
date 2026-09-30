import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_historical_source_provenance,
)


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/signed_average_box_assembly"
M33 = REPO / "reference_tcad/simplemos_sentaurus2022/region_resolved_interface_assembly"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m35_signed_average_box_assembly_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM35SignedAverageBoxAssemblyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            ROOT / "m35_signed_average_box_assembly_report.json").read_text())
        cls.by_variant = {row["variant"]: row
                          for row in cls.report["factorial_response"]}
        cls.effects = {row["factor"]: row
                       for row in cls.report["factorial_main_effects"]}

    def test_factorial_matrix_is_complete_and_converged(self) -> None:
        acceptance = self.report["acceptance"]
        self.assertTrue(acceptance["all_checks_pass"])
        self.assertEqual(acceptance["variant_count"], 8)
        self.assertEqual(acceptance["converged_count"], 8)
        self.assertEqual(acceptance["baseline_replay_relative_error"], 0.0)

    def test_sentaurus_interface_measure_oracle_closes(self) -> None:
        oracle = self.report["sentaurus_measure_oracle"]
        self.assertEqual(oracle["sentaurus_interface_node_count"], 21)
        self.assertLess(oracle["maximum_sentaurus_measure_relative_error"], 3e-14)
        self.assertEqual(oracle["nonpositive_assembled_si_node_measure_count"], 0)

    def test_coherent_signed_geometry_is_best_corner(self) -> None:
        baseline = self.by_variant["p0_t0_v0"]
        all_on = self.by_variant["p1_t1_v1"]
        self.assertAlmostEqual(baseline["absolute_error_dex"], 0.10941868092427424)
        self.assertAlmostEqual(all_on["absolute_error_dex"], 0.08440560595033898)
        self.assertAlmostEqual(
            all_on["error_improvement_vs_baseline_dex"], 0.025013074973935262)
        self.assertEqual(self.report["best_variant"]["variant"], "p1_t1_v1")
        self.assertGreater(self.report["gap_closure"]["all_region_resolved_fraction"], 0.22)

    def test_signed_measure_reverses_barycentric_volume_response(self) -> None:
        m33 = json.loads((M33 / "m33_region_resolved_interface_assembly_report.json").read_text())
        old = {row["factor"]: row["main_effect_log_current_dex"]
               for row in m33["factorial_main_effects"]}
        new = self.effects["transport_signed_average_box_node_volume"][
            "main_effect_log_current_dex"]
        self.assertGreater(old["transport_node_volume"], 0.015)
        self.assertLess(new, -0.005)

    def test_default_model_is_unchanged(self) -> None:
        self.assertFalse(self.report["execution"]["default_model_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        assert_historical_source_provenance(self, EVIDENCE)


if __name__ == "__main__":
    unittest.main()
