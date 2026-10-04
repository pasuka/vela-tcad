import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "direct_current_attribution"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM52DirectCurrentAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = ROOT / "simplemos_m52_direct_current_attribution_contract_v1.json"
        cls.freeze = read_json(
            ROOT / "simplemos_m52_direct_current_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m52_direct_current_attribution_evidence.json")
        cls.report = read_json(
            OUTPUT / "m52_direct_current_attribution_report.json")
        cls.states = rows(OUTPUT / "m52_direct_current_state_ledger.csv")
        cls.terminals = rows(OUTPUT / "m52_terminal_component_ledger.csv")
        cls.fields = rows(OUTPUT / "m52_field_invariance_ledger.csv")

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
            "direct_matches_native_face", "direct_matches_weighted",
            "direct_matches_default", "direct_third_observable",
            "state_perturbed", "replay_incomplete"})
        self.assertTrue(self.report["acceptance"]["classification_declared"])
        target = self.report["target"]
        self.assertEqual(target["device"], "n23")
        self.assertEqual(float(target["gate_voltage_V"]), 0.05)
        if self.report["classification"] == "direct_matches_weighted":
            self.assertTrue(self.report["direct_weighted_agreement"][
                "all_states_within_contract_tolerance"])
            self.assertLess(
                self.report["direct_weighted_agreement"][
                    "maximum_absolute_difference_A_per_um"], 1e-24)

    def test_only_direct_current_was_intervened(self):
        self.assertEqual(self.evidence["only_intervention"],
                         "Math DirectCurrent flag")
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
