import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_historical_source_provenance,
)


REPO = Path(__file__).resolve().parents[2]
REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "bgn_chain_first_divergence/m39_bgn_chain_first_divergence_report.json")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m39_bgn_chain_first_divergence_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM39BgnChainFirstDivergenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def test_report_is_complete_and_all_gates_pass(self) -> None:
        self.assertEqual(self.report["status"], "complete")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertTrue(all(self.report["acceptance"].values()))

    def test_first_divergence_is_m29_model_selection(self) -> None:
        finding = self.report["first_divergence"]
        self.assertEqual(finding["stage"], "C-1_model_selection_contract")
        self.assertEqual(finding["declared"], "none")
        self.assertEqual(finding["observed"], "BennettWilson")
        audit = self.report["m29_contract_audit"]
        self.assertEqual(audit["bgn_off_label_effective_model"], "BennettWilson")
        self.assertFalse(audit["bgn_off_label_matches_expected"])
        self.assertEqual(audit["bgn_off_c1_status"],
                         "diagnostic_cross_model_not_parity")

    def test_m27_oldslotboom_numeric_chain_closes(self) -> None:
        chain = self.report["m27_qualified_oldslotboom_chain"]
        self.assertTrue(chain["mesh_id_identity"])
        self.assertEqual(chain["c0_max_relative_difference"], 0.0)
        self.assertEqual(chain["c1_best_convention"], "total")
        self.assertLess(chain["c1_total_impurity_max_abs_difference_eV"], 1e-12)
        self.assertGreater(chain["c1_abs_net_max_abs_difference_eV"], 0.02)
        self.assertLess(chain["c3_max_joint_fit_rms"], 1e-10)
        self.assertAlmostEqual(chain["c5_conduction_share"], 0.5, places=10)
        self.assertAlmostEqual(chain["c5_valence_share"], 0.5, places=10)
        self.assertAlmostEqual(chain["c5_gap_slope"], -1.0, places=10)

    def test_missing_m29_doping_is_fail_closed(self) -> None:
        audit = self.report["m29_contract_audit"]
        self.assertEqual(audit["native_c0_status"],
                         "not_scored_missing_fields")
        self.assertEqual(audit["downstream_doping_source"],
                         "m27_full_same_mesh_input_oracle")
        self.assertTrue(audit["oracle_mesh_identity"]["qualified"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        self.assertEqual(self.evidence["status"], "frozen")
        for artifact in self.evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        assert_historical_source_provenance(self, EVIDENCE)


if __name__ == "__main__":
    unittest.main()
