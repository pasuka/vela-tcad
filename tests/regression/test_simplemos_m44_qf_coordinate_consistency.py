import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference_tcad/simplemos_sentaurus2022"


class TestSimplemosM44QfCoordinateConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            REFERENCE / "qf_coordinate_consistency"
            / "m44_qf_coordinate_consistency_report.json").read_text())
        cls.evidence = json.loads((
            REFERENCE
            / "simplemos_m44_qf_coordinate_consistency_evidence.json").read_text())
        with (REFERENCE / "qf_coordinate_consistency"
              / "m44_qf_coordinate_consistency_matrix.csv").open(newline="") as f:
            cls.rows = list(csv.DictReader(f))
        with (REFERENCE / "qf_coordinate_consistency"
              / "m44_qf_coordinate_edge_matrix.csv").open(newline="") as f:
            cls.edge_rows = list(csv.DictReader(f))

    def test_evidence_is_frozen_and_complete(self) -> None:
        self.assertEqual(self.report["status"], "complete")
        self.assertEqual(self.evidence["status"], "frozen")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])

    def test_four_frozen_states_cover_modes_and_srh(self) -> None:
        self.assertEqual(len(self.rows), 4)
        self.assertEqual({row["mode"] for row in self.rows},
                         {"legacy", "compensated"})
        self.assertEqual({row["srh"] for row in self.rows},
                         {"srh_on", "srh_off"})
        self.assertGreater(len(self.edge_rows), 0)

    def test_terminal_paths_close_at_deep_off_precision(self) -> None:
        for row in self.rows:
            self.assertLessEqual(
                float(row["cut_vs_terminal_functional_gap_dex"]), 1.0e-12)
            self.assertLessEqual(
                float(row[
                    "terminal_functional_vs_contact_extractor_gap_dex"]),
                1.0e-12)

    def test_reference_coordinates_and_mobility_are_identical(self) -> None:
        for row in self.edge_rows:
            self.assertLessEqual(
                abs(float(row["phin0_relative_difference_V"])), 1.0e-30)
            self.assertLessEqual(
                abs(float(row["phin1_relative_difference_V"])), 1.0e-30)
            self.assertLessEqual(
                float(row["mobility_drive_relative_difference"]), 1.0e-13)
            self.assertLessEqual(
                float(row["mobility_relative_difference"]), 1.0e-13)
            self.assertLessEqual(
                float(row["total_current_relative_difference"]), 1.0e-10)

    def test_m44_does_not_change_default_physics_or_hfs(self) -> None:
        execution = self.report["execution"]
        self.assertFalse(execution["default_physics_model_changed"])
        self.assertFalse(execution["default_hfs_model_changed"])
        self.assertFalse(
            execution["default_equal_ni_flux_evaluation_changed"])
        self.assertTrue(
            self.report["findings"]["reference_increment_packing_was_root_cause"])
        self.assertFalse(
            self.report["findings"]["mobility_reconstruction_was_root_cause"])


if __name__ == "__main__":
    unittest.main()
