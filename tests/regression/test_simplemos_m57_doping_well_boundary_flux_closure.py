import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "doping_well_boundary_flux_closure"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM57DopingWellBoundaryFluxClosureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = (
            ROOT / "simplemos_m57_doping_well_boundary_flux_closure_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m57_doping_well_boundary_flux_closure_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m57_doping_well_boundary_flux_closure_evidence.json")
        cls.report = read_json(
            OUTPUT / "m57_doping_well_boundary_flux_closure_report.json")
        cls.segments = rows(OUTPUT / "m57_boundary_segment_flux_ledger.csv")
        cls.wells = rows(OUTPUT / "m57_well_boundary_summary.csv")
        cls.closures = rows(OUTPUT / "m57_carrier_closure_ledger.csv")
        cls.pairs = rows(OUTPUT / "m57_nwell_pair_summary.csv")

    def test_contract_was_frozen_before_analysis(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])

    def test_read_only_matrix_and_ledgers_are_complete(self):
        self.assertEqual(len(self.wells), 144)
        self.assertEqual(len(self.closures), 288)
        self.assertEqual(len(self.pairs), 48)
        self.assertGreater(len(self.segments), 1000)
        self.assertEqual({row["contact"] for row in self.wells},
                         {"source", "drain", "substrate"})
        self.assertEqual({row["component"] for row in self.closures},
                         {"electron", "hole"})
        self.assertFalse(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["historical_artifacts_rewritten"])

    def test_doping_well_partition_is_unambiguous_and_conservative(self):
        partition = self.report["partition"]
        self.assertEqual(partition["partition_problem_count"], 0)
        self.assertEqual(partition["maximum_labels_per_silicon_triangle"], 2)
        self.assertTrue(partition["topology_invariant_across_bias_states"])
        self.assertEqual(
            partition["maximum_opposite_interface_flux_balance_residual_A_per_um"], 0.0)
        self.assertEqual(len(partition["states_with_unassociated_wells"]), 12)
        self.assertEqual({row["boundary_class"] for row in self.segments},
                         {"physical_contact", "interior_doping_well_interface",
                          "other_semiconductor_exterior"})

    def test_physical_contact_anchor_sets_the_classification(self):
        anchor = self.report["physical_contact_anchor"]
        self.assertFalse(anchor["all_within_tolerance"])
        self.assertEqual(anchor["passing_row_count"], 144)
        self.assertEqual(anchor["row_count"], 288)
        groups = {}
        for contact in ("source", "drain", "substrate"):
            for component in ("electron", "hole"):
                selected = [row for row in self.closures
                            if row["contact"] == contact and row["component"] == component]
                groups[(contact, component)] = sum(
                    row["direct_anchor_within_tolerance"] == "True" for row in selected)
        self.assertEqual(groups[("source", "electron")], 0)
        self.assertEqual(groups[("drain", "electron")], 0)
        self.assertEqual(groups[("drain", "hole")], 0)
        self.assertEqual(groups[("substrate", "electron")], 48)
        self.assertEqual(groups[("substrate", "hole")], 48)
        self.assertEqual(self.report["classification"], "boundary_quadrature_limited")

    def test_target_retains_direct_anchor_and_full_well_residuals(self):
        target = {(row["contact"], row["component"]): row
                  for row in self.report["target"]["closure_rows"]}
        substrate_e = target[("substrate", "electron")]
        self.assertAlmostEqual(
            substrate_e["physical_contact_terminal_oriented_flux_A_per_um"],
            -2.376797223029523e-19, places=30)
        self.assertAlmostEqual(substrate_e["direct_anchor_residual_A_per_um"],
                               -6.835494718731313e-23, places=33)
        self.assertAlmostEqual(substrate_e["default_closure_residual_A_per_um"],
                               -4.1161087816365906e-17, places=28)
        self.assertFalse(substrate_e["default_closure_within_tolerance"])
        self.assertFalse(target[("source", "electron")][
            "direct_anchor_within_tolerance"])
        self.assertFalse(target[("drain", "electron")][
            "direct_anchor_within_tolerance"])

    def test_result_is_accepted_without_overclaiming_closure(self):
        declared = read_json(self.contract_path)["analysis"]["classifications"]
        self.assertIn(self.report["classification"], declared)
        self.assertEqual(self.report["status"], "accepted")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertFalse(self.report["complete_well_closure"]["all_within_tolerance"])
        self.assertEqual(self.report["complete_well_closure"]["passing_row_count"], 39)
        self.assertIn("does not prove a Sentaurus defect", self.report["claim_guard"])

    def test_frozen_artifact_and_source_hashes(self):
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
