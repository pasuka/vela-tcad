import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "terminal_common_mode_attribution"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM54TerminalCommonModeAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = (
            ROOT / "simplemos_m54_terminal_common_mode_attribution_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m54_terminal_common_mode_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m54_terminal_common_mode_attribution_evidence.json")
        cls.report = read_json(
            OUTPUT / "m54_terminal_common_mode_attribution_report.json")
        cls.states = rows(OUTPUT / "m54_state_decomposition_ledger.csv")
        cls.terminals = rows(OUTPUT / "m54_terminal_component_ledger.csv")
        cls.cases = rows(OUTPUT / "m54_case_summary.csv")
        cls.nwell = rows(OUTPUT / "m54_nwell_pair_summary.csv")

    def test_contract_was_frozen_before_analysis(self):
        self.assertEqual(self.freeze["status"], "frozen_before_analysis")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])

    def test_complete_matrix_and_ledgers(self):
        self.assertEqual(len(self.states), 816)
        self.assertEqual(len(self.terminals), 9792)
        self.assertEqual(len(self.cases), 16)
        self.assertEqual(len(self.nwell), 8)
        self.assertEqual({int(row["state_count"]) for row in self.cases}, {51})

    def test_exact_drain_shift_identities_close(self):
        for row in self.states:
            for component in ("electron", "hole", "total"):
                self.assertEqual(
                    row[f"{component}_drain_identity_within_tolerance"], "True")
        self.assertTrue(self.report["acceptance"]["drain_identities_close"])

    def test_classification_is_declared(self):
        declared = read_json(self.contract_path)["analysis"]["classifications"]
        self.assertIn(self.report["classification"], declared)
        self.assertEqual(self.evidence["classification"],
                         self.report["classification"])
        self.assertEqual(self.report["status"], "complete")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])

    def test_analysis_is_read_only(self):
        self.assertFalse(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["historical_artifacts_rewritten"])
        self.assertFalse(self.evidence["closed_topics_reinvestigated"])
        self.assertEqual(len(self.evidence["raw_terminal_input_hashes"]), 32)

    def test_frozen_artifact_and_source_hashes_match(self):
        self.assertEqual(self.evidence["status"], "frozen")
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
