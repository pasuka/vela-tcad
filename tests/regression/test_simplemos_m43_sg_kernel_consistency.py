import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference_tcad/simplemos_sentaurus2022"


class TestSimplemosM43SgKernelConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            REFERENCE / "sg_kernel_consistency"
            / "m43_sg_kernel_consistency_report.json").read_text())
        cls.evidence = json.loads((
            REFERENCE
            / "simplemos_m43_sg_kernel_consistency_evidence.json").read_text())
        with (REFERENCE / "sg_kernel_consistency"
              / "m43_sg_kernel_consistency_matrix.csv").open(newline="") as f:
            cls.rows = list(csv.DictReader(f))

    def test_evidence_is_frozen_and_complete(self) -> None:
        self.assertEqual(self.report["status"], "complete")
        self.assertEqual(self.evidence["status"], "frozen")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])

    def test_four_frozen_states_cover_both_numerical_modes(self) -> None:
        self.assertEqual(len(self.rows), 4)
        self.assertEqual({row["mode"] for row in self.rows},
                         {"legacy", "compensated"})
        self.assertEqual({row["srh"] for row in self.rows},
                         {"srh_on", "srh_off"})

    def test_edge_cut_functional_and_extractor_are_consistent(self) -> None:
        for row in self.rows:
            self.assertLessEqual(
                float(row["cut_vs_terminal_functional_gap_dex"]), 0.005)
            self.assertLessEqual(
                float(row[
                    "terminal_functional_vs_contact_extractor_gap_dex"]),
                0.005)

    def test_default_numerical_mode_and_physics_are_unchanged(self) -> None:
        execution = self.report["execution"]
        self.assertFalse(execution["default_physics_model_changed"])
        self.assertFalse(
            execution["default_equal_ni_flux_evaluation_changed"])


if __name__ == "__main__":
    unittest.main()
