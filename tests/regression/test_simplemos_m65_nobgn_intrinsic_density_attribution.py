import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_contract_freeze.json"
EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
REPORT = ROOT / "nobgn_intrinsic_density_attribution/m65_nobgn_intrinsic_density_attribution_report.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM65EvidenceTest(unittest.TestCase):
    def test_contract_was_frozen_before_execution(self):
        freeze = load(FREEZE)
        self.assertEqual(freeze["status"], "frozen_before_execution")
        self.assertEqual(freeze["contract_sha256"], digest(CONTRACT))

    def test_is_accepted_single_axis_evidence(self):
        report = load(REPORT)
        evidence = load(EVIDENCE)
        self.assertEqual(report["status"], "accepted")
        self.assertTrue(report["acceptance"]["all_checks_pass"])
        self.assertEqual(evidence["status"], "frozen")
        self.assertEqual(report["execution"]["single_axis"], "copied Vela Si.ni")
        self.assertFalse(report["execution"]["production_default_changed"])
        self.assertFalse(evidence["production_reference_replaced"])

    def test_pair_and_state_ledgers_are_complete(self):
        pair_path = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
        state_path = ROOT / "nobgn_intrinsic_density_attribution/m65_state_attribution_ledger.csv"
        with pair_path.open(newline="", encoding="utf-8") as stream:
            pairs = list(csv.DictReader(stream))
        with state_path.open(newline="", encoding="utf-8") as stream:
            states = list(csv.DictReader(stream))
        self.assertEqual(len(pairs), 8)
        self.assertEqual(len(states), 16)
        self.assertTrue(all(float(row["sentaurus_ni_relative_spread"]) <= 1e-12
                            for row in states))
        self.assertTrue(all(float(row["boltzmann_identity_p95_residual_dex"]) <= 5e-4
                            for row in states))

    def test_frozen_artifact_hashes_match(self):
        evidence = load(EVIDENCE)
        for relative, expected in evidence["artifacts"].items():
            self.assertEqual(digest(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
