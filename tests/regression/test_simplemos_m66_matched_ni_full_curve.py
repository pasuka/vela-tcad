import csv, hashlib, json
from pathlib import Path
import unittest

REPO=Path(__file__).resolve().parents[2];ROOT=REPO/"reference_tcad/simplemos_sentaurus2022"
CONTRACT=ROOT/"simplemos_m66_matched_ni_full_curve_contract_v1.json";FREEZE=ROOT/"simplemos_m66_matched_ni_full_curve_contract_freeze.json";EVIDENCE=ROOT/"simplemos_m66_matched_ni_full_curve_evidence.json";REPORT=ROOT/"matched_ni_full_curve/m66_matched_ni_full_curve_report.json"
def load(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
class SimpleMosM66EvidenceTest(unittest.TestCase):
    def test_frozen_accepted_and_complete(self):
        self.assertEqual(load(FREEZE)["contract_sha256"],digest(CONTRACT));report=load(REPORT);self.assertEqual(report["status"],"accepted");self.assertTrue(report["acceptance"]["all_checks_pass"])
        for relative,count in (("matched_ni_full_curve/m66_matched_ni_full_curve_point_ledger.csv",816),("matched_ni_full_curve/m66_matched_ni_full_curve_pair_ledger.csv",8),("matched_ni_full_curve/m66_execution_case_ledger.csv",16)):
            with (ROOT/relative).open(newline="",encoding="utf-8") as stream:self.assertEqual(len(list(csv.DictReader(stream))),count)
    def test_artifact_hashes(self):
        for relative,expected in load(EVIDENCE)["artifacts"].items():self.assertEqual(digest(REPO/relative),expected)
if __name__=="__main__":unittest.main()
