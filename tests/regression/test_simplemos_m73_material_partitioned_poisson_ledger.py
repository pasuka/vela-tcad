import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_contract_v2.json"
FREEZE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_contract_freeze_v2.json"
EVIDENCE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_evidence.json"
V1_EVIDENCE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_evidence_v1_superseded.json"
REPORT = ROOT / "material_partitioned_poisson_ledger/m73_material_partitioned_poisson_ledger_report.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM73EvidenceTest(unittest.TestCase):
    def test_frozen_accepted_material_partition_ledgers(self):
        self.assertEqual(load(FREEZE)["contract_sha256"], digest(CONTRACT))
        report = load(REPORT)
        self.assertTrue(report["acceptance"]["all_checks_pass"])
        self.assertEqual(
            report["classification"],
            "material_partition_material_but_not_dominant")
        summary = report["summary"]
        self.assertEqual(
            summary["best_material_partition_variant"],
            "region_local_barycentric_si")
        self.assertEqual(summary["best_same_sign_current_pair_count"], 8)
        self.assertLess(summary["best_median_current_closure_fraction"], 0.8)
        self.assertEqual(
            summary["best_electrostatic_partition_variant"],
            "legacy_signed_si")
        self.assertGreater(
            summary["best_median_electrostatic_closure_fraction"], 0.8)
        for name, count in (
            ("m73_poisson_stage_ledger.csv", 160),
            ("m73_poisson_pair_ledger.csv", 40),
            ("m73_variant_summary.csv", 5),
        ):
            path = ROOT / "material_partitioned_poisson_ledger" / name
            with path.open(newline="", encoding="utf-8") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), count)

    def test_v1_supersession_and_v2_evidence_hashes(self):
        self.assertEqual(load(V1_EVIDENCE)["status"], "superseded")
        evidence = load(EVIDENCE)
        self.assertEqual(evidence["status"], "frozen")
        self.assertFalse(evidence["new_sentaurus_execution"])
        self.assertFalse(evidence["new_vela_self_consistent_execution"])
        for relative, expected in evidence["artifacts"].items():
            self.assertEqual(digest(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
