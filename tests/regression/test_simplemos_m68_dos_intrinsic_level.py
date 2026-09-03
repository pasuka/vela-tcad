import csv, hashlib, json
from pathlib import Path
import unittest
REPO=Path(__file__).resolve().parents[2];ROOT=REPO/"reference_tcad/simplemos_sentaurus2022";CONTRACT=ROOT/"simplemos_m68_dos_intrinsic_level_contract_v1.json";FREEZE=ROOT/"simplemos_m68_dos_intrinsic_level_contract_freeze.json";EVIDENCE=ROOT/"simplemos_m68_dos_intrinsic_level_evidence.json";REPORT=ROOT/"dos_intrinsic_level/m68_dos_intrinsic_level_report.json"
def load(path):return json.loads(path.read_text(encoding="utf-8-sig"))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
class SimpleMosM68EvidenceTest(unittest.TestCase):
    def test_null_intervention_is_exact(self):
        self.assertEqual(load(FREEZE)["contract_sha256"],digest(CONTRACT));report=load(REPORT);self.assertTrue(report["acceptance"]["all_checks_pass"]);self.assertEqual(report["summary"]["absolute_current_difference_A_per_um"],0.0);self.assertTrue(report["summary"]["state_hash_identity"])
        with (ROOT/"dos_intrinsic_level/m68_dos_null_intervention_ledger.csv").open(newline="",encoding="utf-8") as stream:self.assertEqual(len(list(csv.DictReader(stream))),2)
    def test_artifact_hashes(self):
        for relative,expected in load(EVIDENCE)["artifacts"].items():self.assertEqual(digest(REPO/relative),expected)
if __name__=="__main__":unittest.main()
