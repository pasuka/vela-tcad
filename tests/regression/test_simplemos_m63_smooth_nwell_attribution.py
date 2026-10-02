import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "smooth_nwell_attribution"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM63SmoothNwellAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.freeze = read_json(
            ROOT / "simplemos_m63_smooth_nwell_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m63_smooth_nwell_attribution_evidence.json")
        cls.report = read_json(OUTPUT / "m63_smooth_nwell_attribution_report.json")
        cls.points = rows(OUTPUT / "m63_smooth_window_point_ledger.csv")
        cls.curves = rows(OUTPUT / "m63_curve_shape_ledger.csv")
        cls.thresholds = rows(OUTPUT / "m63_constant_current_threshold_ledger.csv")
        cls.pairs = rows(OUTPUT / "m63_nwell_pair_attribution_ledger.csv")
        cls.mobility = rows(OUTPUT / "m63_m13_mobility_crosscheck_ledger.csv")

    def test_contract_and_upstream_hashes_are_frozen(self):
        contract = ROOT / "simplemos_m63_smooth_nwell_attribution_contract_v1.json"
        self.assertEqual(self.freeze["status"], "frozen_before_analysis")
        self.assertEqual(self.freeze["contract_sha256"], sha256(contract))
        for relative, expected in self.freeze["upstream_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)

    def test_matrix_is_complete(self):
        self.assertEqual(len(self.points), 128)
        self.assertEqual(len(self.curves), 16)
        self.assertGreaterEqual(len(self.thresholds), 48)
        self.assertEqual(len(self.pairs), 8)
        self.assertEqual(len(self.mobility), 2)

    def test_horizontal_shift_closes_curve_and_pair_growth(self):
        self.assertGreaterEqual(
            self.report["summary"]["minimum_curve_horizontal_shift_explained_fraction"],
            0.80)
        self.assertGreaterEqual(
            self.report["summary"]["minimum_pair_growth_horizontal_shift_explained_fraction"],
            0.80)
        self.assertTrue(self.report["summary"]["all_nwell_growth_positive"])

    def test_m13_mobility_chain_does_not_share_growth_sign(self):
        self.assertEqual(self.report["summary"]["mobility_sign_support_count"], 0)
        self.assertTrue(all(int(row["sign_supports_mobility_growth"]) == 0
                            for row in self.mobility))

    def test_result_is_scoped_and_frozen(self):
        self.assertEqual(self.report["status"], "accepted")
        self.assertEqual(self.report["classification"],
                         "threshold_like_horizontal_shift_dominant")
        self.assertIn("not_closed", self.report["causal_scope"])
        self.assertEqual(self.evidence["status"], "frozen")
        self.assertFalse(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
