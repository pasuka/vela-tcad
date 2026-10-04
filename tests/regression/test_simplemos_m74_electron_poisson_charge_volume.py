import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m74_electron_poisson_charge_volume_contract_v1.json"
FREEZE = ROOT / "simplemos_m74_electron_poisson_charge_volume_contract_freeze_v1.json"
EVIDENCE = ROOT / "simplemos_m74_electron_poisson_charge_volume_evidence.json"
REPORT = ROOT / "electron_poisson_charge_volume/m74_electron_poisson_charge_volume_report.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM74EvidenceTest(unittest.TestCase):
    def test_frozen_accepted_electron_only_self_consistent_ab(self):
        self.assertEqual(load(FREEZE)["contract_sha256"], digest(CONTRACT))
        report = load(REPORT)
        self.assertTrue(report["acceptance"]["all_checks_pass"])
        self.assertEqual(
            report["classification"],
            "electron_poisson_volume_material_but_not_dominant")
        self.assertEqual(report["execution"]["new_sentaurus_solves"], 0)
        self.assertEqual(report["execution"]["new_vela_candidate_workflows"], 16)
        self.assertFalse(report["execution"]["production_default_changed"])

        summary = report["summary"]
        self.assertEqual(summary["improved_pair_count"], 8)
        self.assertGreaterEqual(summary["median_pair_closure_fraction"], 0.3)
        self.assertLess(summary["median_pair_closure_fraction"], 0.8)
        self.assertGreater(
            summary["high_nwell_median_candidate_absolute_error_dex"],
            summary["high_nwell_median_baseline_absolute_error_dex"])
        self.assertGreater(
            summary["maximum_curve_candidate_over_baseline_abs_log_shift_dex"],
            3.0)

        for name, count in (
            ("m74_case_ledger.csv", 16),
            ("m74_pair_ledger.csv", 8),
            ("m74_curve_response_ledger.csv", 274),
        ):
            path = ROOT / "electron_poisson_charge_volume" / name
            with path.open(newline="", encoding="utf-8") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), count)

    def test_evidence_hashes(self):
        evidence = load(EVIDENCE)
        self.assertEqual(evidence["status"], "frozen")
        self.assertFalse(evidence["new_sentaurus_execution"])
        self.assertTrue(evidence["new_vela_self_consistent_execution"])
        self.assertFalse(evidence["production_reference_replaced"])
        for relative, expected in evidence["artifacts"].items():
            self.assertEqual(digest(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
