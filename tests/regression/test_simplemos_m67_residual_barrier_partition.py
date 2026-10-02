import csv, hashlib, json
from pathlib import Path
import unittest
REPO=Path(__file__).resolve().parents[2];ROOT=REPO/"reference_tcad/simplemos_sentaurus2022";CONTRACT=ROOT/"simplemos_m67_residual_barrier_partition_contract_v1.json";FREEZE=ROOT/"simplemos_m67_residual_barrier_partition_contract_freeze.json";EVIDENCE=ROOT/"simplemos_m67_residual_barrier_partition_evidence.json";REPORT=ROOT/"residual_barrier_partition/m67_residual_barrier_partition_report.json"
def load(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
class SimpleMosM67EvidenceTest(unittest.TestCase):
    def test_frozen_accepted_ledgers(self):
        self.assertEqual(load(FREEZE)["contract_sha256"],digest(CONTRACT));self.assertTrue(load(REPORT)["acceptance"]["all_checks_pass"])
        for name,count in (("m67_fixed_node_state_ledger.csv",16),("m67_pair_partition_ledger.csv",8)):
            with (ROOT/"residual_barrier_partition"/name).open(newline="",encoding="utf-8") as stream:self.assertEqual(len(list(csv.DictReader(stream))),count)
    def test_artifact_hashes(self):
        for relative,expected in load(EVIDENCE)["artifacts"].items():self.assertEqual(digest(REPO/relative),expected)
if __name__=="__main__":unittest.main()
