import csv, hashlib, json
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m71_gate_electrostatic_partition_contract_v2.json"
FREEZE = ROOT / "simplemos_m71_gate_electrostatic_partition_contract_freeze_v2.json"
EVIDENCE = ROOT / "simplemos_m71_gate_electrostatic_partition_evidence.json"
V1_EVIDENCE = ROOT / "simplemos_m71_gate_electrostatic_partition_evidence_v1_failed.json"
REPORT = ROOT / "gate_electrostatic_partition/m71_gate_electrostatic_partition_report.json"


def load(path): return json.loads(path.read_text(encoding="utf-8-sig"))
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM71EvidenceTest(unittest.TestCase):
    def test_v2_frozen_accepted_ledgers(self):
        self.assertEqual(load(FREEZE)["contract_sha256"], digest(CONTRACT))
        report = load(REPORT)
        self.assertTrue(report["acceptance"]["all_checks_pass"])
        self.assertFalse(report["summary"]["gate_charge_proxy_qualified"])
        for name, count in (("m71_gate_electrostatic_state_ledger.csv", 32),
                            ("m71_gate_electrostatic_pair_ledger.csv", 8)):
            path = ROOT / "gate_electrostatic_partition" / name
            with path.open(newline="", encoding="utf-8") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), count)

    def test_v1_failure_and_v2_artifact_hashes(self):
        v1 = load(V1_EVIDENCE)
        self.assertEqual(v1["status"], "failed")
        for relative, expected in v1["artifacts"].items():
            self.assertEqual(digest(REPO / relative), expected)
        for relative, expected in load(EVIDENCE)["artifacts"].items():
            self.assertEqual(digest(REPO / relative), expected)


if __name__ == "__main__": unittest.main()
