import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "compensation_contour_mesh_audit"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM58CompensationContourMeshAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = (
            ROOT / "simplemos_m58_compensation_contour_mesh_audit_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m58_compensation_contour_mesh_audit_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m58_compensation_contour_mesh_audit_evidence.json")
        cls.report = read_json(
            OUTPUT / "m58_compensation_contour_mesh_audit_report.json")
        cls.devices = rows(OUTPUT / "m58_device_compensation_mesh_ledger.csv")
        cls.contours = rows(OUTPUT / "m58_zero_contour_segment_ledger.csv")
        cls.components = rows(OUTPUT / "m58_connectivity_component_ledger.csv")
        cls.critical = rows(OUTPUT / "m58_critical_gate_edge_node_ledger.csv")
        cls.species = rows(OUTPUT / "m58_pair_species_coordinate_ledger.csv")
        cls.pairs = rows(OUTPUT / "m58_nwell_error_attribution_ledger.csv")

    def test_contract_was_frozen_before_analysis(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])
        self.assertIn("before any M58 analysis",
                      self.freeze["pre_execution_correction"])

    def test_read_only_matrix_and_ledgers_are_complete(self):
        self.assertEqual(len(self.devices), 8)
        self.assertEqual(len(self.pairs), 8)
        self.assertEqual(len(self.species), 4)
        self.assertEqual(len(self.contours), 856)
        self.assertEqual(len(self.components), 28)
        self.assertEqual(len(self.critical), 34)
        self.assertFalse(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertTrue(self.evidence["local_immutable_tdr_export_only"])

    def test_active_species_exactly_reconstruct_net_active(self):
        for row in self.devices:
            self.assertEqual(row["active_field_identity_pass"], "True")
            self.assertEqual(
                float(row["active_field_identity_max_abs_residual_cm3"]), 0.0)
        for row in self.species:
            self.assertAlmostEqual(
                float(row["BActive_median_log10_shift_dex"]),
                0.30103, places=5)
            self.assertGreater(int(row["low_only_coordinate_count"]), 100)
            self.assertGreater(int(row["high_only_coordinate_count"]), 100)

    def test_target_has_two_single_node_floating_n_components(self):
        devices = {row["device"]: row for row in self.devices}
        n19, n23 = devices["n19"], devices["n23"]
        self.assertEqual(int(n19["floating_n_type_component_count"]), 0)
        self.assertEqual(int(n23["floating_n_type_component_count"]), 2)
        self.assertEqual(int(n23["single_node_floating_n_type_component_count"]), 2)
        self.assertAlmostEqual(float(n23["gate_edge_min_abs_net_cm3"]),
                               1.59668102306504e15, places=1)
        self.assertAlmostEqual(float(n23["gate_edge_max_compensation_ratio"]),
                               99.90624829734554, places=12)
        floating = [row for row in self.components
                    if row["device"] == "n23" and row["doping_type"] == "n_type"
                    and row["floating"] == "True"]
        self.assertEqual({int(row["node_count"]) for row in floating}, {1})
        critical = [row for row in self.critical
                    if row["device"] == "n23"
                    and int(row["node_id"]) in {298, 1041}]
        self.assertEqual(len(critical), 2)
        self.assertEqual({float(row["y_um"]) for row in critical}, {-0.125, 0.125})

    def test_full_matrix_breaks_single_compensation_explanation(self):
        tests = self.report["systematic_tests"]
        self.assertEqual(tests["error_increased_rows"], 8)
        self.assertEqual(tests["pairs_with_stronger_high_compensation"], 3)
        self.assertEqual(tests["pairs_with_topology_or_closer_contour_change"], 2)
        self.assertEqual(tests["low_ldd_rows_with_topology_or_closer_contour_change"], 4)
        self.assertEqual(tests["high_ldd_rows_with_topology_or_closer_contour_change"], 0)
        n20_n24 = [row for row in self.pairs
                   if row["low_device"] == "n20" and row["high_device"] == "n24"]
        self.assertEqual(len(n20_n24), 2)
        self.assertTrue(all(row["error_increased"] == "True" for row in n20_n24))
        self.assertTrue(all(
            row["high_has_stronger_gate_edge_compensation"] == "False"
            for row in n20_n24))
        self.assertTrue(all(
            row["high_has_topology_or_closer_contour_change"] == "False"
            for row in n20_n24))

    def test_result_is_accepted_without_causal_overclaim(self):
        declared = read_json(self.contract_path)["analysis"]["classifications"]
        self.assertIn(self.report["classification"], declared)
        self.assertEqual(self.report["classification"],
                         "no_systematic_frozen_tdr_relation")
        self.assertEqual(self.report["status"], "accepted")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertIn("does not prove", self.report["claim_guard"])

    def test_frozen_artifact_and_source_hashes(self):
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
