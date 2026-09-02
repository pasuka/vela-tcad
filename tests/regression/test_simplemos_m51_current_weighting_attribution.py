import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "current_weighting_attribution"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM51CurrentWeightingAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = ROOT / "simplemos_m51_current_weighting_attribution_contract_v1.json"
        cls.freeze = read_json(
            ROOT / "simplemos_m51_current_weighting_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m51_current_weighting_attribution_evidence.json")
        cls.report = read_json(
            OUTPUT / "m51_current_weighting_attribution_report.json")
        cls.states = rows(OUTPUT / "m51_current_weighting_state_ledger.csv")
        cls.terminals = rows(OUTPUT / "m51_terminal_component_ledger.csv")
        cls.fields = rows(OUTPUT / "m51_field_invariance_ledger.csv")

    def test_contract_was_frozen_before_execution(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])

    def test_six_state_and_terminal_matrix_are_complete(self):
        self.assertEqual(len(self.states), 6)
        self.assertEqual(len(self.terminals), 72)
        self.assertEqual(
            {(row["device"], float(row["gate_voltage_V"])) for row in self.states},
            {(device, gate) for device in ("n23", "n19")
             for gate in (0.0, 0.05, 0.1)})

    def test_state_invariance_and_native_replay(self):
        self.assertGreater(len(self.fields), 0)
        self.assertTrue(self.report["acceptance"]["state_fields_invariant"])
        self.assertTrue(self.report["acceptance"]["native_flux_replays"])
        self.assertEqual(self.report["state_field_summary"]["total_failure_count"], 0)

    def test_classification_is_contract_declared(self):
        self.assertIn(self.report["classification"], {
            "weighted_matches_native_face", "weighted_matches_default",
            "weighted_third_observable", "state_perturbed", "replay_incomplete"})
        self.assertTrue(self.report["acceptance"]["classification_declared"])
        target = self.report["target"]
        self.assertEqual(target["device"], "n23")
        self.assertEqual(float(target["gate_voltage_V"]), 0.05)
        self.assertGreater(
            float(target["absolute_default_native_gap_reduction_fraction"]),
            0.99999)
        if self.report["classification"] == "weighted_third_observable":
            self.assertFalse(target["weighted_matches_native_normal"])

    def test_only_current_weighting_was_intervened(self):
        self.assertEqual(self.evidence["only_intervention"],
                         "Math CurrentWeighting flag")
        self.assertTrue(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["production_defaults_changed"])
        self.assertFalse(self.evidence["closed_topics_reinvestigated"])

    def test_frozen_artifact_and_source_hashes_match(self):
        self.assertEqual(self.evidence["status"], "frozen")
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
