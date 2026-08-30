"""Regression checks for frozen SimpleMOS M18 QF edge localization."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m18_qf_edge_localization_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM18QfEdgeLocalizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"][
            "m18_qf_edge_localization_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))
        edge_entry = cls.evidence["artifacts"][
            "m18_key_state_edge_contributions.csv"]
        with (REPO / edge_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.edges = list(csv.DictReader(handle))

    def test_scope_and_execution_contract_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"],
                         "complete_with_bias_regime_separation")
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["edge_probe_count"], 16)
        self.assertEqual(execution["probe_variant_count"], 8)
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_edge_and_feedback_ledgers_close(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_edge_to_node_absolute_A_per_um"],
                        1.0e-12)
        self.assertLess(closure["maximum_edge_to_node_relative"], 1.0e-8)
        self.assertLess(closure["maximum_feedback_ledger_absolute_A_per_um"],
                        1.0e-12)
        self.assertLess(closure["maximum_edge_shapley_flux_closure"], 1.0e-12)

    def test_key_state_qf_response_is_drain_cut_localized(self) -> None:
        key = self.report["key_state"]
        self.assertEqual(key["state"], "n21_vd_0p05_vg_0p8")
        self.assertGreater(
            key["combined_contact_cut_absolute_edge_support_fraction"], 0.80)
        self.assertLess(key["internal_channel_absolute_edge_support_fraction"],
                        0.20)
        self.assertGreater(
            key["qf_bucket_absolute_edge_support_fraction"]["drain_contact_cut"],
            0.80)
        self.assertEqual(key["top_edge_count_for_50_percent_absolute_qf_support"],
                         1)

    def test_leading_edges_concentrate_qf_support(self) -> None:
        self.assertEqual(self.edges[0]["bucket"], "drain_contact_cut")
        self.assertEqual(self.edges[1]["bucket"], "drain_contact_cut")
        self.assertGreater(float(self.edges[0]["qf_cumulative_absolute_fraction"]),
                           0.66)
        self.assertGreater(float(self.edges[1]["qf_cumulative_absolute_fraction"]),
                           0.79)

    def test_direct_srh_operator_response_is_negligible(self) -> None:
        ledger = self.report["key_state"][
            "feedback_ledger_relative_to_baseline_current"]
        self.assertGreater(ledger["qf_log_imbalance"], 0.13)
        self.assertLess(abs(ledger["srh"]), 1.0e-9)

    def test_portable_artifacts_and_figures_match_frozen_hashes(self) -> None:
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
