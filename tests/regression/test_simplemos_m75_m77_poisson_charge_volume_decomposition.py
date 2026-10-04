import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"

MILESTONES = {
    "M75": {
        "contract": "simplemos_m75_hole_poisson_charge_volume_contract_v1.json",
        "freeze": "simplemos_m75_hole_poisson_charge_volume_contract_freeze_v1.json",
        "evidence": "simplemos_m75_hole_poisson_charge_volume_evidence.json",
        "report": "hole_poisson_charge_volume/m75_hole_poisson_charge_volume_report.json",
        "classification": "hole_poisson_volume_not_material",
    },
    "M76": {
        "contract": "simplemos_m76_dopant_poisson_charge_volume_contract_v1.json",
        "freeze": "simplemos_m76_dopant_poisson_charge_volume_contract_freeze_v1.json",
        "evidence": "simplemos_m76_dopant_poisson_charge_volume_evidence.json",
        "report": "dopant_poisson_charge_volume/m76_dopant_poisson_charge_volume_report.json",
        "classification": "dopant_poisson_volume_not_material",
    },
    "M77": {
        "contract": "simplemos_m77_combined_poisson_charge_volume_contract_v1.json",
        "freeze": "simplemos_m77_combined_poisson_charge_volume_contract_freeze_v1.json",
        "evidence": "simplemos_m77_combined_poisson_charge_volume_evidence.json",
        "report": "combined_poisson_charge_volume/m77_combined_poisson_charge_volume_report.json",
        "classification": "combined_poisson_volume_not_material",
    },
}


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM75M77EvidenceTest(unittest.TestCase):
    def test_contracts_results_and_artifact_hashes(self):
        reports = {}
        for milestone, item in MILESTONES.items():
            contract = ROOT / item["contract"]
            freeze = load(ROOT / item["freeze"])
            self.assertEqual(freeze["contract_sha256"], digest(contract))
            report = load(ROOT / item["report"])
            reports[milestone] = report
            self.assertTrue(report["acceptance"]["all_checks_pass"])
            self.assertEqual(report["classification"], item["classification"])
            self.assertEqual(report["execution"]["new_sentaurus_solves"], 0)
            self.assertEqual(report["execution"]["new_vela_candidate_workflows"], 16)
            self.assertFalse(report["execution"]["production_default_changed"])
            evidence = load(ROOT / item["evidence"])
            self.assertEqual(evidence["status"], "frozen")
            self.assertFalse(evidence["new_sentaurus_execution"])
            self.assertTrue(evidence["new_vela_self_consistent_execution"])
            for relative, expected in evidence["artifacts"].items():
                self.assertEqual(digest(REPO / relative), expected)

        self.assertLess(
            reports["M75"]["summary"]["median_pair_closure_fraction"],
            1.0e-5)
        self.assertEqual(reports["M76"]["summary"]["improved_pair_count"], 0)
        self.assertLess(
            reports["M76"]["summary"]["median_self_consistent_pair_reduction_dex"],
            0.0)
        self.assertEqual(reports["M77"]["summary"]["improved_pair_count"], 3)
        self.assertLess(
            reports["M77"]["summary"]["maximum_absolute_nonlinear_interaction_dex"],
            2.0e-4)

    def test_portable_ledger_row_counts(self):
        for milestone, slug in (
            ("m75", "hole_poisson_charge_volume"),
            ("m76", "dopant_poisson_charge_volume"),
            ("m77", "combined_poisson_charge_volume"),
        ):
            for suffix, count in (
                ("case_ledger.csv", 16),
                ("pair_ledger.csv", 8),
                ("curve_response_ledger.csv", 274),
            ):
                path = ROOT / slug / f"{milestone}_{suffix}"
                with path.open(newline="", encoding="utf-8") as stream:
                    self.assertEqual(len(list(csv.DictReader(stream))), count)
        decomposition = (ROOT / "combined_poisson_charge_volume" /
                         "m77_charge_volume_decomposition_ledger.csv")
        with decomposition.open(newline="", encoding="utf-8") as stream:
            self.assertEqual(len(list(csv.DictReader(stream))), 8)


if __name__ == "__main__":
    unittest.main()
