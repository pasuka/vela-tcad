"""Regression checks for frozen SimpleMOS M21 incident-edge evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m21_incident_edge_balance_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM21IncidentEdgeBalanceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"][
            "m21_incident_edge_balance_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))
        edge_entry = cls.evidence["artifacts"]["m21_key_state_incident_edges.csv"]
        with (REPO / edge_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.edges = {int(row["edge_id"]): row for row in csv.DictReader(handle)}

    def test_scope_and_production_drive_are_frozen(self) -> None:
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["sg_probe_count"], 32)
        self.assertTrue(execution["all_states_use_transport_cell_vector"])
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_all_numerical_closures_pass(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_stable_sg_flux_absolute_closure"], 1.0e-12)
        self.assertLess(closure["maximum_edge_factor_absolute_closure"], 1.0e-12)
        self.assertLess(closure[
            "maximum_m20_node_flux_delta_absolute_closure"], 1.0e-12)
        self.assertEqual(closure["maximum_geometry_difference"], 0.0)

    def test_key_and_cross_state_response_is_qf_controlled(self) -> None:
        findings = self.report["findings"]
        self.assertGreater(findings["key_factor_absolute_shares"][
            "qf_log_imbalance"], 0.9999)
        self.assertLess(findings["key_factor_absolute_shares"]["mobility"], 1.0e-4)
        self.assertLess(findings["key_factor_absolute_shares"][
            "sg_secant_conductance"], 1.0e-4)
        self.assertGreater(findings["all_state_qf_absolute_share_range"][0], 0.997)

    def test_leading_edges_exclude_mobility_density_bernoulli_and_geometry(self) -> None:
        for edge_id in (2107, 2094, 2088, 2114):
            row = self.edges[edge_id]
            self.assertGreater(abs(float(row[
                "transport_cell_vector_drive_change_dex"])), 0.05)
            self.assertLess(abs(float(row["mobility_change_dex"])), 3.0e-6)
            self.assertLess(float(row["maximum_endpoint_density_change_dex"]),
                            1.0e-5)
            self.assertLess(float(row["maximum_bernoulli_weight_change_dex"]),
                            2.0e-6)
            self.assertEqual(float(row["length_difference_m"]), 0.0)
            self.assertEqual(float(row["couple_difference_m"]), 0.0)
            self.assertGreater(float(row["qf_log_imbalance_absolute_share"]),
                               0.9999)

    def test_portable_artifacts_and_figures_match_hashes(self) -> None:
        for artifact in self.evidence["artifacts"].values():
            path = REPO / artifact["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), artifact["sha256"], path)
        for figure in self.evidence["figures"]:
            path = REPO / figure["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), figure["sha256"], path)


if __name__ == "__main__":
    unittest.main()
