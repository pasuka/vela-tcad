"""Regression checks for frozen SimpleMOS M20 drain-adjacent row evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m20_drain_adjacent_row_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM20DrainAdjacentRowAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"][
            "m20_drain_adjacent_row_audit_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))
        ledger_entry = cls.evidence["artifacts"]["m20_key_state_node_ledger.csv"]
        with (REPO / ledger_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.key = {int(row["node_id"]): row for row in csv.DictReader(handle)}

    def test_scope_and_execution_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"], "complete")
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["probe_count"], 96)
        self.assertTrue(execution["diagnostic_cpp_changed"])
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_flux_term_and_jacobian_partitions_close(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_raw_edge_flux_closure"], 1.0e-10)
        self.assertLess(closure["maximum_raw_term_sum_closure"], 1.0e-10)
        self.assertLess(closure[
            "maximum_exact_jacobian_partition_relative_error"], 1.0e-12)

    def test_key_state_is_flux_balance_not_local_srh(self) -> None:
        key = self.report["findings"]["key_state"]
        self.assertEqual(key["state"], "n21_vd_0p05_vg_0p8")
        self.assertEqual(key["adjacent_nodes"], "990;991;992;995;997")
        self.assertLess(key["absolute_srh_delta_share"], 1.0e-15)
        self.assertGreater(key["absolute_contact_flux_delta_share"], 0.49)
        self.assertLess(key["absolute_contact_flux_delta_share"], 0.51)
        self.assertGreater(key["absolute_internal_flux_delta_share"], 0.49)
        self.assertLess(key["absolute_internal_flux_delta_share"], 0.51)

    def test_nodes_990_992_localization_is_frozen(self) -> None:
        self.assertEqual(set(self.key), {990, 991, 992, 995, 997})
        self.assertEqual(float(self.key[990]["delta_contact_flux"]), 0.0)
        for node, minimum_cancellation in ((991, 1000.0), (992, 900.0)):
            row = self.key[node]
            cancellation = (
                abs(float(row["delta_contact_flux"]))
                + abs(float(row["delta_internal_flux"]))) / abs(
                    float(row["delta_raw_residual"]))
            self.assertGreater(cancellation, minimum_cancellation)
            self.assertGreater(float(row["sentaurus_phin_jacobian_share"]),
                               0.99999999)
            self.assertGreater(float(row["sentaurus_phin_contact_column_share"]),
                               0.44)
            self.assertLess(float(row["sentaurus_phin_contact_column_share"]),
                            0.45)

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
