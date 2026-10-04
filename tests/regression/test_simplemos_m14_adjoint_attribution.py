"""Regression checks for frozen SimpleMOS M14 adjoint attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m14_adjoint_attribution_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM14AdjointAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"]["m14_adjoint_attribution_report.json"]
        cls.report = json.loads((REPO / report_entry["path"]).read_text(encoding="utf-8"))
        summary_entry = cls.evidence["artifacts"]["m14_substitution_summary.csv"]
        with (REPO / summary_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.rows = list(csv.DictReader(handle))

    def row(self, state: str, variant: str) -> dict[str, str]:
        return next(row for row in self.rows
                    if row["state"] == state and row["variant"] == variant)

    def test_scope_and_coverage_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"],
                         "complete_with_declared_field_coverage")
        execution = self.report["execution"]
        self.assertEqual(execution["weak_region_terminal_point_count"], 80)
        self.assertEqual(execution["paired_spatial_state_count"], 16)
        self.assertEqual(execution["adjoint_solve_count"], 16)
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_adjoint_linear_system_is_well_resolved(self) -> None:
        self.assertLess(self.report["adjoint_quality"]["maximum_relative_residual"],
                        1.0e-12)

    def test_weak_region_spectrum_is_complete_and_one_sided(self) -> None:
        weak = self.report["weak_region"]
        self.assertTrue(weak["all_vela_above_sentaurus"])
        self.assertAlmostEqual(weak["mean_signed_error_dex"],
                               0.051096304466763776, places=14)
        self.assertAlmostEqual(weak["maximum_absolute_error_dex"],
                               0.10941868092427424, places=14)

    def test_n21_vg08_field_and_region_attribution(self) -> None:
        phin = self.row("n21_vd_0p05_vg_0p8", "field_phin")
        psi = self.row("n21_vd_0p05_vg_0p8", "field_psi")
        channel = self.row("n21_vd_0p05_vg_0p8", "region_channel")
        drain = self.row("n21_vd_0p05_vg_0p8", "region_drain")
        self.assertGreater(abs(float(phin["direct_exact_relative_to_current"])),
                           1.0e4 * abs(float(psi["direct_exact_relative_to_current"])))
        self.assertAlmostEqual(float(channel["two_layer_sum_A_per_um"]),
                               -1.0075693806237023e-08, places=18)
        self.assertLess(abs(float(drain["two_layer_sum_A_per_um"])), 1.0e-13)

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
