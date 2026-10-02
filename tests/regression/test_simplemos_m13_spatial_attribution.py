"""Regression checks for frozen SimpleMOS M13 spatial attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m13_spatial_attribution_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM13SpatialAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"]["report"]
        cls.report = json.loads((REPO / report_entry["path"]).read_text(
            encoding="utf-8"))
        summary_entry = cls.evidence["artifacts"]["state_summary"]
        with (REPO / summary_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.rows = list(csv.DictReader(handle))

    def test_scope_and_execution_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"], "complete")
        self.assertEqual(self.report["execution"]["state_count"], 16)
        self.assertEqual(self.report["execution"]["edge_projection_probe_count"], 16)
        self.assertEqual(self.report["execution"]["new_self_consistent_vg0p8_state_count"], 4)
        self.assertFalse(self.report["execution"]["new_sentaurus_execution"])
        self.assertFalse(self.report["execution"]["cpp_changed"])
        self.assertFalse(self.report["execution"]["default_model_changed"])

    def test_drive_discretization_result_is_regime_dependent(self) -> None:
        findings = self.report["findings"]
        self.assertEqual(findings["edge_projection_improves_active_p95_state_count"], 8)
        self.assertEqual(findings["edge_projection_low_drain_improves_state_count"], 8)
        self.assertEqual(findings["edge_projection_high_drain_worsens_state_count"], 8)
        self.assertEqual(findings["sentaurus_drive_improves_active_p95_state_count"], 16)
        means = findings["mean_active_p95_mobility_error_dex"]
        self.assertAlmostEqual(means["transport_cell_vector"], 0.236464370238995,
                               places=12)
        self.assertAlmostEqual(means["sentaurus_drive"], 0.16405366585144,
                               places=12)
        self.assertGreater(means["edge_projection"], means["transport_cell_vector"])

    def test_barrier_proxy_is_directionally_consistent_but_incomplete(self) -> None:
        barrier = self.report["findings"]["high_nwell_low_drain_weak_region_barrier"]
        self.assertEqual(barrier["representative_state_count"], 3)
        self.assertEqual(barrier["sign_consistent_with_terminal_error_count"], 3)
        self.assertGreater(barrier["mean_absolute_prediction_residual_dex"], 0.02)
        self.assertLess(barrier["mean_absolute_prediction_residual_dex"], 0.04)

    def test_old_slotboom_formula_is_not_the_residual_source(self) -> None:
        self.assertLess(
            self.report["findings"]["maximum_old_slotboom_p95_difference_meV"],
            1.0e-10)
        self.assertGreater(
            self.report["findings"]["maximum_surface_psi_p95_difference_mV"],
            15.0)

    def test_portable_artifacts_and_figures_match_frozen_hashes(self) -> None:
        for entry in self.evidence["artifacts"].values():
            path = REPO / entry["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), entry["sha256"], path)
        for entry in self.evidence["figures"]:
            path = REPO / entry["path"]
            self.assertTrue(path.is_file(), path)
            self.assertEqual(sha256(path), entry["sha256"], path)

    def test_state_summary_contains_exact_matrix(self) -> None:
        self.assertEqual(len(self.rows), 16)
        observed = {(row["device"], float(row["drain_voltage_V"]),
                     float(row["gate_voltage_V"])) for row in self.rows}
        expected = {(device, drain, gate)
                    for device in ("n17", "n21")
                    for drain in (0.05, 1.0)
                    for gate in (0.0, 0.05, 0.8, 2.5)}
        self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
