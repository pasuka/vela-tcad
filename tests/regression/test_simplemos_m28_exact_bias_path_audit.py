import csv
import json
import math
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/exact_bias_path_audit"


class SimpleMosM28ExactBiasPathAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((ROOT / "m28_exact_bias_path_audit_report.json").read_text())
        with (ROOT / "m28_terminal_route_ledger.csv").open(newline="") as stream:
            cls.routes = {row["route"]: row for row in csv.DictReader(stream)}

    def test_execution_and_acceptance_are_frozen(self) -> None:
        self.assertEqual(self.report["execution"]["new_sentaurus_state_count"], 6)
        self.assertFalse(self.report["execution"]["default_model_changed"])
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])

    def test_original_initial_route_exactly_reproduces_m8(self) -> None:
        baseline = self.routes["m8_currentplot"]
        matched = self.routes["original_initial_dozero"]
        for field in ("full_current_A_per_um", "no_hfs_current_A_per_um",
                      "hfs_response_A_per_um"):
            self.assertEqual(float(matched[field]), float(baseline[field]))

    def test_dozero_is_null_for_tiny_schedule(self) -> None:
        m27 = self.routes["m27_tiny_endpoint"]
        dozero = self.routes["tiny_dozero"]
        for variant in ("full_current_A_per_um", "no_hfs_current_A_per_um"):
            self.assertLess(abs(float(dozero[variant]) - float(m27[variant])), 2.1e-20)

    def test_hfs_response_is_route_stable(self) -> None:
        self.assertLess(
            self.report["terminal"]["maximum_hfs_response_relative_difference_from_m8"],
            0.015)

    def test_path_shift_is_common_mode(self) -> None:
        terminal = self.report["terminal"]
        self.assertGreater(terminal["maximum_common_mode_offset_A_per_um"], 2.5e-16)
        self.assertLess(terminal["maximum_differential_offset_A_per_um"], 4.0e-19)

    def test_state_response_remains_aligned(self) -> None:
        for route in self.report["state_response"]["routes"]:
            self.assertEqual(route["node_count"], 942)
            self.assertTrue(math.isclose(
                route["hfs_phin_response_pearson_vs_m27"], 1.0,
                rel_tol=0.0, abs_tol=1e-12))


if __name__ == "__main__":
    unittest.main()
