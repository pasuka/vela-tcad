import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "active_species_topology_intervention"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM62ActiveSpeciesTopologyInterventionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.freeze = read_json(
            ROOT / "simplemos_m62_active_species_topology_intervention_contract_freeze.json")
        cls.erratum = read_json(
            ROOT / "simplemos_m62_active_species_topology_intervention_contract_erratum_v1.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m62_active_species_topology_intervention_evidence.json")
        cls.report = read_json(OUTPUT / "m62_active_species_topology_intervention_report.json")
        cls.fields = rows(OUTPUT / "m62_roundtrip_field_invariance_ledger.csv")
        cls.tools = rows(OUTPUT / "m62_writer_capability_ledger.csv")
        cls.candidates = rows(OUTPUT / "m62_candidate_node_ledger.csv")

    def test_contract_erratum_and_upstream_hashes_are_frozen(self):
        contract = ROOT / "simplemos_m62_active_species_topology_intervention_contract_v1.json"
        self.assertEqual(self.freeze["status"], "frozen_before_tool_inventory")
        self.assertEqual(self.freeze["contract_sha256"], sha256(contract))
        self.assertEqual(self.erratum["contract_sha256"], sha256(contract))
        self.assertEqual(self.erratum["corrected_interpretation"], "accepted")
        for relative, expected in self.freeze["upstream_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)

    def test_no_edit_roundtrip_fails_unrelated_field_invariance(self):
        gate = self.report["writer_gate"]
        self.assertFalse(gate["passed"])
        self.assertTrue(gate["topology_files_byte_identical"])
        self.assertTrue(gate["doping_files_byte_identical"])
        self.assertEqual(gate["field_file_count"], 132)
        self.assertEqual(gate["changed_unrelated_field_file_count"], 36)
        self.assertAlmostEqual(gate["representative_displacement_scale_factor"],
                               10000.0)

    def test_stop_rule_prevents_mutation_and_device_solve(self):
        self.assertEqual(self.report["classification"],
                         "e3_stopped_no_identity_preserving_writer")
        self.assertTrue(self.report["writer_gate"]["stop_rule_applied"])
        self.assertFalse(self.report["execution"]["candidate_mutation_attempted"])
        self.assertFalse(self.report["execution"]["new_sentaurus_device_execution"])
        self.assertFalse(self.report["execution"]["mutated_tdr_created"])

    def test_candidates_and_tool_inventory_are_complete(self):
        self.assertEqual(len(self.candidates), 4)
        self.assertEqual({row["device"] for row in self.candidates}, {"n21", "n23"})
        self.assertTrue(all(row["mutation_attempted"] == "False"
                            for row in self.candidates))
        self.assertEqual(len(self.tools), 4)
        tdx = next(row for row in self.tools if row["tool"] == "tdx Tcl interface")
        self.assertEqual(tdx["available"], "True")
        self.assertEqual(tdx["roundtrip_preserves_unrelated_fields"], "False")

    def test_evidence_is_frozen(self):
        self.assertEqual(self.report["status"], "accepted")
        self.assertEqual(self.evidence["status"], "frozen")
        self.assertFalse(self.evidence["new_sentaurus_device_execution"])
        self.assertFalse(self.evidence["mutated_tdr_created"])
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)
        for relative, expected in self.evidence["diagnostic_source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
