import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "substrate_electron_transport_attribution"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM49SubstrateElectronTransportAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = ROOT / "simplemos_m49_substrate_electron_transport_attribution_contract_v1.json"
        cls.freeze = read_json(
            ROOT / "simplemos_m49_substrate_electron_transport_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m49_substrate_electron_transport_attribution_evidence.json")
        cls.report = read_json(
            OUTPUT / "m49_substrate_electron_transport_attribution_report.json")
        cls.boundary = rows(OUTPUT / "m49_substrate_boundary_edge_transport_ledger.csv")
        cls.layers = rows(OUTPUT / "m49_substrate_graph_distance_row_ledger.csv")
        cls.bridge = rows(OUTPUT / "m49_transport_factor_bridge_ledger.csv")
        cls.controls = rows(OUTPUT / "m49_control_localization_ledger.csv")

    def test_contract_was_frozen_before_analysis(self):
        self.assertEqual(self.freeze["status"], "frozen_before_analysis")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"], self.freeze["contract_sha256"])

    def test_exact_state_and_spatial_ledgers_are_complete(self):
        self.assertEqual(len(self.report["states"]), 6)
        self.assertEqual(len(self.layers), 30)
        self.assertEqual(len(self.bridge), 30)
        self.assertEqual(len(self.controls), 6)
        self.assertGreater(len(self.boundary), 0)
        self.assertEqual(
            {(row["device"], float(row["gate_voltage_V"])) for row in self.controls},
            {(device, gate) for device in ("n23", "n19")
             for gate in (0.0, 0.05, 0.1)})

    def test_exact_factor_bridges_close_and_anchor_m48(self):
        acceptance = self.report["acceptance"]
        self.assertTrue(acceptance["m48_terminal_anchors_reproduced"])
        self.assertTrue(acceptance["shapley_identities_close"])
        self.assertTrue(acceptance["exact_bridges_close"])
        self.assertTrue(acceptance["all_checks_pass"])
        self.assertAlmostEqual(
            self.report["target"]["terminal_difference_A_per_um"],
            -4.6947034945330634e-17, delta=1e-30)

    def test_claim_obeys_boundary_observable_gate(self):
        finding = self.report["finding"]
        if finding["boundary_observable_limited"]:
            self.assertEqual(finding["classification"], "boundary_observable_limited")
            self.assertGreater(
                finding["boundary_observable_residual_absolute_fraction"], 0.20)
            self.assertFalse(finding["resolved_mechanism_gate_passed"])
        if finding["resolved_mechanism_gate_passed"]:
            self.assertTrue(finding["dominant_factor_target_localized"])
            self.assertGreaterEqual(
                finding["dominant_absolute_fraction_of_target_difference"], 0.80)
            self.assertFalse(finding["boundary_observable_limited"])

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
