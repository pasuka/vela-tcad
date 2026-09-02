import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "si_oxide_interface_topology"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM56SiOxideInterfaceTopologyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = (
            ROOT / "simplemos_m56_si_oxide_interface_topology_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m56_si_oxide_interface_topology_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m56_si_oxide_interface_topology_evidence.json")
        cls.report = read_json(
            OUTPUT / "m56_si_oxide_interface_topology_report.json")
        cls.topology = rows(OUTPUT / "m56_input_topology_ledger.csv")
        cls.states = rows(OUTPUT / "m56_solved_state_interface_ledger.csv")
        cls.nodes = rows(OUTPUT / "m56_interface_node_potential_ledger.csv")
        cls.support = rows(OUTPUT / "m56_field_support_ledger.csv")

    def test_contract_was_frozen_before_analysis(self):
        self.assertEqual(self.freeze["status"], "frozen_before_analysis")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])

    def test_topology_and_state_matrices_are_complete(self):
        self.assertEqual(len(self.topology), 8)
        self.assertEqual(len(self.states), 6)
        self.assertGreater(len(self.nodes), 0)
        self.assertGreater(len(self.support), 0)
        self.assertTrue(all(int(row["silicon_oxide_shared_edge_count"]) > 0
                            for row in self.topology))

    def test_no_geometric_duplicates_and_contacts_attach_by_region(self):
        self.assertTrue(all(int(row["duplicate_coordinate_group_count"]) == 0
                            for row in self.topology))
        self.assertTrue(all(row["gate_contact_region"] == "Oxide_1"
                            for row in self.topology))
        for row in self.topology:
            self.assertEqual(row["source_contact_region"], "Silicon_1")
            self.assertEqual(row["drain_contact_region"], "Silicon_1")
            self.assertEqual(row["substrate_contact_region"], "Silicon_1")

    def test_potential_is_continuous_and_carrier_support_is_silicon_only(self):
        tolerance = float(read_json(self.contract_path)["acceptance"][
            "maximum_interface_potential_absolute_mismatch_V"])
        self.assertTrue(all(abs(float(row[
            "silicon_minus_oxide_potential_V"])) <= tolerance
                            for row in self.nodes))
        self.assertTrue(all(
            row["carrier_density_current_support_regions"] == "Silicon_1"
            for row in self.states))
        self.assertTrue(self.report["acceptance"][
            "interface_potential_continuous"])
        self.assertTrue(self.report["acceptance"]["carrier_support_silicon_only"])

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
        self.assertGreater(len(self.evidence["raw_input_hashes"]), 0)

    def test_frozen_artifact_and_source_hashes_match(self):
        self.assertEqual(self.evidence["status"], "frozen")
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
