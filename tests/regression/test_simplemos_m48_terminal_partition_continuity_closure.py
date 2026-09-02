import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM48TerminalPartitionContinuityClosureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = ROOT / "simplemos_m48_terminal_partition_continuity_closure_contract_v1.json"
        cls.freeze = read_json(
            ROOT / "simplemos_m48_terminal_partition_continuity_closure_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m48_terminal_partition_continuity_closure_evidence.json")
        cls.report = read_json(
            ROOT / "terminal_partition_continuity_closure/m48_terminal_partition_continuity_closure_report.json")
        cls.terminals = rows(
            ROOT / "terminal_partition_continuity_closure/m48_four_terminal_component_ledger.csv")
        cls.continuity = rows(
            ROOT / "terminal_partition_continuity_closure/m48_carrier_continuity_closure_ledger.csv")
        cls.contributions = rows(
            ROOT / "terminal_partition_continuity_closure/m48_terminal_identity_contribution_ledger.csv")
        cls.controls = rows(
            ROOT / "terminal_partition_continuity_closure/m48_control_contrast_ledger.csv")

    def test_contract_was_frozen_before_analysis(self):
        self.assertEqual(self.freeze["status"], "frozen_before_analysis")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"], self.freeze["contract_sha256"])

    def test_exact_fixed_state_ledgers_are_complete(self):
        self.assertEqual(len(self.terminals), 144)
        self.assertEqual(len(self.continuity), 12)
        self.assertEqual(len(self.contributions), 84)
        self.assertEqual(len(self.controls), 4)
        self.assertEqual(
            {(row["solver"], row["device"], float(row["gate_voltage_V"]))
             for row in self.continuity},
            {(solver, device, gate)
             for solver in ("sentaurus", "vela")
             for device in ("n23", "n19")
             for gate in (0.0, 0.05, 0.1)})

    def test_terminal_and_component_identities_close(self):
        acceptance = self.report["acceptance"]
        self.assertTrue(acceptance["component_identities_close"])
        self.assertTrue(acceptance["terminal_identities_close"])
        self.assertTrue(acceptance["m47_currents_reproduced"])
        self.assertTrue(acceptance["all_checks_pass"])

    def test_claim_is_partition_not_causality(self):
        finding = self.report["finding"]
        self.assertIn("substrate_partition", finding["classification"])
        self.assertEqual(finding["dominant_bookkeeping_contribution"],
                         "substrate_rebalancing")
        self.assertTrue(finding["target_both_solvers_numerically_resolved"])
        self.assertFalse(finding["target_srh_material_per_contract"])
        self.assertIn("not by themselves identify", finding["causal_guard"])

    def test_closed_topics_and_solver_states_stayed_unchanged(self):
        self.assertTrue(self.evidence["analysis_only"])
        self.assertFalse(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["default_physics_model_changed"])
        self.assertFalse(self.evidence["closed_topics_reinvestigated"])

    def test_frozen_artifact_and_source_hashes_match(self):
        self.assertEqual(self.evidence["status"], "frozen")
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
