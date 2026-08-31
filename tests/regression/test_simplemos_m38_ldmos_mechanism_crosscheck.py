import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "ldmos_mechanism_crosscheck/m38_ldmos_mechanism_crosscheck_report.json")


class SimpleMosM38LdmosMechanismCrosscheckTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_snapshot_is_qualified_but_not_claimed_frozen(self) -> None:
        self.assertEqual(self.report["status"], "qualified_read_only_snapshot")
        self.assertTrue(self.report["execution"]["ldmos_worktree_dirty_observed"])
        self.assertFalse(self.report["execution"]["ldmos_files_modified"])
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])

    def test_shared_qf_feedback_and_distinct_contact_response(self) -> None:
        simple = self.report["simplemos"]
        ldmos = self.report["ldmos"]
        self.assertEqual(simple["contact_sg_direct_electron_delta_A_per_um"], 0.0)
        self.assertGreater(ldmos["contact_hfs_operator_improvement_dex"], 1.0)
        self.assertGreater(ldmos["feedback_fraction_of_self_consistent_error"], 0.9)
        self.assertGreater(ldmos["phin_fraction_of_feedback"], 0.99)

    def test_signed_cotangent_does_not_change_feedback_attribution(self) -> None:
        ldmos = self.report["ldmos"]
        for key in (
            "signed_cotangent_operator_change_dex",
            "signed_cotangent_self_consistent_change_dex",
            "signed_cotangent_feedback_change_dex",
        ):
            self.assertLess(abs(ldmos[key]), 1e-3)


if __name__ == "__main__":
    unittest.main()
