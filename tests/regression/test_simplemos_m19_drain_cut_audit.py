"""Regression checks for frozen SimpleMOS M19 drain-cut audit."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m19_drain_cut_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM19DrainCutAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"]["m19_drain_cut_audit_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))
        state_entry = cls.evidence["artifacts"]["m19_state_summary.csv"]
        with (REPO / state_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.states = list(csv.DictReader(handle))

    def test_scope_and_execution_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"],
                         "complete_with_bias_regime_separation")
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["sg_probe_count"], 144)
        self.assertEqual(execution["probe_variants_per_state"], 9)
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["cpp_changed"])
        self.assertFalse(execution["default_model_changed"])

    def test_contact_boundary_and_numerical_closures_pass(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_contact_phin_bias_error_V"], 1.0e-12)
        self.assertLess(closure[
            "maximum_strong_state_sg_cut_current_reconstruction_relative_error"],
            1.0e-6)
        self.assertLess(closure["maximum_shapley_closure_A_per_um"], 1.0e-18)

    def test_imported_state_closes_stable_vg08_gap(self) -> None:
        vg08 = [row for row in self.states
                if abs(float(row["gate_voltage_V"]) - 0.8) < 1.0e-12]
        self.assertEqual(len(vg08), 4)
        for row in vg08:
            self.assertLess(abs(float(row[
                "sentaurus_state_vela_sg_error_dex"])), 5.0e-4)
            self.assertGreater(float(row[
                "full_sentaurus_state_gap_closed_fraction"]), 0.97)
            self.assertLess(float(row[
                "full_sentaurus_state_gap_closed_fraction"]), 1.01)

    def test_key_state_is_adjacent_interior_qf_controlled(self) -> None:
        key = self.report["findings"]["key_state"]
        self.assertEqual(key["state"], "n21_vd_0p05_vg_0p8")
        self.assertLess(key["max_contact_phin_difference_V"], 1.0e-12)
        shares = key["frozen_factor_absolute_shares"]
        self.assertGreater(shares["drain_adjacent_interior_phin"], 0.999)
        self.assertLess(shares["drain_contact_endpoint_state"], 1.0e-3)
        self.assertLess(shares["drain_adjacent_interior_other_state"], 1.0e-3)
        frozen = abs(key[
            "frozen_drain_adjacent_interior_phin_relative_to_Id"])
        adjoint = abs(key["m18_drain_cut_qf_response_relative_to_Id"])
        self.assertLess(abs(frozen - adjoint), 2.0e-6)

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
