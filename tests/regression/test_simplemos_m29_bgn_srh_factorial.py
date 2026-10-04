import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/bgn_srh_factorial"


class SimpleMosM29BgnSrhFactorialTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(
            (ROOT / "m29_bgn_srh_factorial_report.json").read_text())
        cls.cells = {row["cell"]: row for row in cls.report["matrix"]["cells"]}

    def test_execution_and_acceptance_are_frozen(self) -> None:
        execution = self.report["execution"]
        self.assertEqual(execution["new_sentaurus_state_count"], 3)
        self.assertEqual(execution["vela_state_count"], 4)
        self.assertEqual(execution["reused_m28_state_count"], 1)
        self.assertFalse(execution["default_model_changed"])
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertEqual(self.report["acceptance"]["common_silicon_node_count"], 942)

    def test_baseline_and_best_cell_are_distinct(self) -> None:
        self.assertAlmostEqual(
            self.report["matrix"]["baseline_absolute_gap_dex"],
            0.09039871358507992, places=14)
        self.assertEqual(self.report["matrix"]["best_gap_cell"],
                         "bgn_off_srh_on")
        self.assertLess(self.report["matrix"]["best_absolute_gap_dex"], 0.0024)

    def test_local_single_factor_responses_are_close(self) -> None:
        sent = self.report["factor_effects"]["sentaurus"]
        vela = self.report["factor_effects"]["vela"]
        self.assertLess(abs(vela["bgn_effect_srh_on_dex"]
                            - sent["bgn_effect_srh_on_dex"]), 0.089)
        self.assertLess(abs(vela["srh_effect_bgn_on_dex"]
                            - sent["srh_effect_bgn_on_dex"]), 0.016)

    def test_double_off_cell_exposes_interaction_mismatch(self) -> None:
        self.assertGreater(
            self.cells["bgn_off_srh_off"]["absolute_log10_ratio"], 1.32)
        sent = self.report["factor_effects"]["sentaurus"]
        vela = self.report["factor_effects"]["vela"]
        self.assertGreater(abs(vela["bgn_srh_interaction_dex"]
                               - sent["bgn_srh_interaction_dex"]), 1.30)

    def test_nodal_factor_response_patterns_remain_aligned(self) -> None:
        spatial = {(row["factor"], row["field"]): row
                   for row in self.report["spatial_effects"]}
        self.assertGreater(spatial[("bgn", "phin")]["pearson"], 0.99)
        self.assertGreater(spatial[("bgn", "logn")]["pearson"], 0.96)
        self.assertGreater(spatial[("srh", "phin")]["pearson"], 0.9999)
        self.assertGreater(spatial[("srh", "logn")]["pearson"], 0.9999)


if __name__ == "__main__":
    unittest.main()
