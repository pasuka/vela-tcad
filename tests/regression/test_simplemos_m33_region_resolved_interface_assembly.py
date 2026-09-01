import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_source_hashes_current_or_m44,
)


REPO = Path(__file__).resolve().parents[2]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "region_resolved_interface_assembly")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m33_region_resolved_interface_assembly_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM33RegionResolvedInterfaceAssemblyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            ROOT / "m33_region_resolved_interface_assembly_report.json").read_text())
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

    def test_geometry_audit_quantifies_shared_interface_weighting(self) -> None:
        geometry = self.report["geometry"]
        self.assertEqual(geometry["si_sio2_interface_edge_count"], 20)
        self.assertEqual(geometry["si_sio2_interface_node_count"], 21)
        self.assertAlmostEqual(
            geometry["transport_couple_ratio_median"], 4.984083881213519)
        self.assertAlmostEqual(
            geometry["poisson_coefficient_ratio_median"], 1.4272691926338978)
        self.assertGreater(geometry["total_over_si_node_volume_max"], 11.28)

    def test_coherent_region_resolution_has_bounded_improvement(self) -> None:
        baseline = self.by_variant["p0_t0_v0"]
        all_on = self.by_variant["p1_t1_v1"]
        self.assertAlmostEqual(baseline["absolute_error_dex"], 0.10941868092427424)
        self.assertAlmostEqual(all_on["absolute_error_dex"], 0.10476463269852927)
        self.assertAlmostEqual(
            all_on["error_improvement_vs_baseline_dex"], 0.004654048225744972)
        self.assertGreater(
            self.report["gap_closure"]["all_region_resolved_fraction"], 0.04)
        self.assertLess(
            self.report["gap_closure"]["all_region_resolved_fraction"], 0.05)

    def test_transport_couple_and_volume_effects_partly_cancel(self) -> None:
        transport = self.effects["transport_edge_coupling"][
            "main_effect_log_current_dex"]
        volume = self.effects["transport_node_volume"][
            "main_effect_log_current_dex"]
        poisson = self.effects["poisson_edge_coupling"][
            "main_effect_log_current_dex"]
        self.assertLess(transport, -0.018)
        self.assertGreater(volume, 0.015)
        self.assertLess(poisson, 0.0)
        self.assertLess(abs(poisson), abs(transport))

    def test_lowest_error_corner_is_not_the_coherent_all_on_variant(self) -> None:
        best = self.report["best_variant"]
        self.assertEqual(best["variant"], "p1_t1_v0")
        self.assertFalse(best["transport_node_volume"])
        self.assertFalse(self.report["execution"]["default_model_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        assert_source_hashes_current_or_m44(self, evidence)


if __name__ == "__main__":
    unittest.main()
