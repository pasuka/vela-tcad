import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
REPORT = ROOT / "full_matrix_requalification/m46_full_matrix_requalification_report.json"
EVIDENCE = ROOT / "simplemos_m46_full_matrix_requalification_evidence.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


class SimpleMosM46FullMatrixRequalificationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def test_complete_matrix_passes_without_default_changes(self) -> None:
        self.assertEqual(self.report["status"], "complete")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertEqual(self.report["execution"]["curve_count"], 16)
        self.assertEqual(self.report["execution"]["direct_bias_point_count"],
                         816)
        self.assertFalse(self.report["execution"]["default_model_changed"])
        self.assertFalse(self.report["execution"]["new_sentaurus_execution"])

    def test_cross_tcad_acceptance_is_unchanged(self) -> None:
        finding = self.report["findings"]
        self.assertEqual(finding["passing_curve_count"], 16)
        self.assertEqual(finding["trend_match_count"], 16)
        self.assertAlmostEqual(
            finding["maximum_absolute_log10_ratio_dex"],
            0.10941868092427424, places=15)
        self.assertAlmostEqual(
            finding["maximum_relative_error"],
            0.28652633602019956, places=15)
        self.assertAlmostEqual(
            finding["maximum_absolute_endpoint_log10_ratio_dex"],
            0.013174412341053722, places=15)

    def test_all_curves_and_points_are_identical_to_m8(self) -> None:
        finding = self.report["findings"]
        self.assertEqual(finding["bitwise_identical_case_csv_count"], 16)
        self.assertEqual(finding["maximum_log_shift_vs_historical_m8_dex"],
                         0.0)
        self.assertEqual(
            finding["maximum_absolute_current_delta_vs_historical_m8_A_per_um"],
            0.0)
        self.assertFalse(
            finding["m43_m44_changed_original_physics_production_curves"])
        point_rows = rows(
            ROOT / "full_matrix_requalification/m46_pointwise_delta.csv")
        self.assertEqual(len(point_rows), 816)
        self.assertTrue(all(float(row["absolute_log_current_shift_dex"]) == 0.0
                            for row in point_rows))

    def test_worst_point_remains_n23_low_drain_weak_inversion(self) -> None:
        worst = self.report["findings"]["worst_case"]
        self.assertEqual(worst["case"], "n23_vd_0p05")
        self.assertEqual(worst["gate_voltage_V"], 0.05)
        self.assertAlmostEqual(worst["m46_absolute_log10_ratio_dex"],
                               0.10941868092427424, places=15)

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        self.assertEqual(self.evidence["status"], "frozen")
        for artifact in self.evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(expected, sha256(REPO / relative))


if __name__ == "__main__":
    unittest.main()
