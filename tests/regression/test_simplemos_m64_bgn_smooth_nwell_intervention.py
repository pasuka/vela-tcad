import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "bgn_smooth_nwell_intervention"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM64BgnSmoothNwellInterventionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.freeze = read_json(
            ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_evidence.json")
        cls.report = read_json(
            OUTPUT / "m64_bgn_smooth_nwell_intervention_report.json")
        cls.points = rows(OUTPUT / "m64_bgn_on_off_point_ledger.csv")
        cls.curves = rows(OUTPUT / "m64_bgn_on_off_curve_ledger.csv")
        cls.pairs = rows(OUTPUT / "m64_bgn_nwell_pair_ledger.csv")
        cls.cases = rows(OUTPUT / "m64_execution_case_ledger.csv")

    def test_contract_and_upstream_hashes_are_frozen(self):
        contract = ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_contract_v1.json"
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(self.freeze["contract_sha256"], sha256(contract))
        for relative, expected in self.freeze["upstream_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)

    def test_matrix_and_model_qualification_are_complete(self):
        self.assertEqual(len(self.points), 128)
        self.assertEqual(len(self.curves), 16)
        self.assertEqual(len(self.pairs), 8)
        self.assertEqual(len(self.cases), 32)
        sentaurus = [row for row in self.cases if row["solver"] == "sentaurus"]
        vela = [row for row in self.cases if row["solver"] == "vela"]
        self.assertEqual(len(sentaurus), 16)
        self.assertEqual(len(vela), 16)
        self.assertTrue(all(int(row["bgn_off_qualified"]) == 1
                            for row in self.cases))

    def test_difference_in_differences_is_an_exact_identity(self):
        tolerance = 1e-12
        self.assertLessEqual(max(abs(float(row["identity_residual_dex"]))
                                 for row in self.points), tolerance)
        self.assertLessEqual(max(abs(float(
            row["difference_in_differences_identity_residual_dex"]))
            for row in self.pairs), tolerance)

    def test_m63_bgn_on_anchor_is_replayed(self):
        self.assertLessEqual(max(abs(float(row["m63_anchor_replay_error_dex"]))
                                 for row in self.pairs), 1e-12)

    def test_result_is_scoped_and_frozen(self):
        self.assertEqual(self.report["status"], "accepted")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertIn(self.report["classification"], {
            "bgn_self_consistent_response_dominant",
            "bgn_material_but_not_dominant",
            "bgn_independent_threshold_shift_dominant",
            "bgn_changes_curve_regime_noncomparably",
        })
        self.assertEqual(self.evidence["status"], "frozen")
        self.assertTrue(self.evidence["new_sentaurus_execution"])
        self.assertTrue(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["production_reference_replaced"])
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
