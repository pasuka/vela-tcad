import csv
import hashlib
import json
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m72_native_gate_reaction_contract_v3.json"
FREEZE = ROOT / "simplemos_m72_native_gate_reaction_contract_freeze_v3.json"
EVIDENCE = ROOT / "simplemos_m72_native_gate_reaction_evidence.json"
V1_EVIDENCE = ROOT / "simplemos_m72_native_gate_reaction_evidence_v1_failed.json"
REPORT = ROOT / "native_gate_reaction/m72_native_gate_reaction_report.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM72EvidenceTest(unittest.TestCase):
    def test_v2_native_charge_and_profile_ledgers(self):
        self.assertEqual(load(FREEZE)["contract_sha256"], digest(CONTRACT))
        report = load(REPORT)
        self.assertTrue(report["acceptance"]["all_checks_pass"])
        self.assertTrue(report["summary"]["native_charge_qualified"])
        self.assertEqual(
            report["classification"],
            "native_charge_qualified_full_profile_material_but_not_dominant")
        for name, count in (
            ("m72_native_gate_reaction_stage_ledger.csv", 32),
            ("m72_interface_profile_ledger.csv", 640),
            ("m72_interface_response_ledger.csv", 320),
            ("m72_native_gate_reaction_pair_ledger.csv", 8),
        ):
            path = ROOT / "native_gate_reaction" / name
            with path.open(newline="", encoding="utf-8") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), count)

    def test_v1_failure_and_v2_evidence_hashes(self):
        v1 = load(V1_EVIDENCE)
        self.assertEqual(v1["status"], "failed")
        for relative, expected in v1["hashes"].items():
            self.assertEqual(digest(REPO / relative), expected)
        evidence = load(EVIDENCE)
        self.assertEqual(evidence["status"], "frozen")
        for relative, expected in evidence["hashes"].items():
            self.assertEqual(digest(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
