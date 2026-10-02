"""Regression checks for frozen SimpleMOS M15 operator-adjoint attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m15_operator_adjoint_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM15OperatorAdjointTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"]["m15_operator_adjoint_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))

        component_entry = cls.evidence["artifacts"]["m15_component_contributions.csv"]
        with (REPO / component_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.components = list(csv.DictReader(handle))
        region_entry = cls.evidence["artifacts"]["m15_region_contributions.csv"]
        with (REPO / region_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.regions = list(csv.DictReader(handle))
        crosscheck_entry = cls.evidence["artifacts"]["m15_m11_mobility_crosscheck.csv"]
        with (REPO / crosscheck_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.crosschecks = list(csv.DictReader(handle))

    @staticmethod
    def select(rows: list[dict[str, str]], state: str,
               field: str, value: str) -> dict[str, str]:
        return next(row for row in rows
                    if row["state"] == state and row[field] == value)

    def test_scope_and_read_only_contract_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"],
                         "complete_with_declared_operator_boundary")
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["carrier_term_probe_count"], 32)
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["cpp_changed"])
        self.assertFalse(execution["default_model_changed"])

    def test_equation_term_sum_closes_to_m14_adjoint_response(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_absolute_A_per_um"], 1.0e-12)
        self.assertLess(closure["maximum_relative"], 1.0e-8)

    def test_key_state_component_attribution(self) -> None:
        state = "n21_vd_0p05_vg_0p8"
        electron = self.select(self.components, state, "component",
                               "electron_transport")
        poisson = self.select(self.components, state, "component", "poisson")
        srh = self.select(self.components, state, "component", "electron_srh")
        self.assertAlmostEqual(float(electron["relative_to_baseline_current"]),
                               0.1269922606019102, places=14)
        self.assertAlmostEqual(float(poisson["relative_to_baseline_current"]),
                               0.0057591010019365715, places=14)
        self.assertLess(abs(float(srh["relative_to_baseline_current"])), 1.0e-8)
        self.assertGreater(float(electron["absolute_contribution_fraction"]), 0.95)

    def test_key_state_spatial_attribution(self) -> None:
        state = "n21_vd_0p05_vg_0p8"
        channel = self.select(self.regions, state, "region", "channel")
        drain = self.select(self.regions, state, "region", "drain")
        self.assertAlmostEqual(float(channel["relative_to_baseline_current"]),
                               0.13525399930103885, places=14)
        self.assertAlmostEqual(float(drain["relative_to_baseline_current"]),
                               -0.002503136537435597, places=14)

    def test_m11_crosscheck_does_not_identify_hfs_uniquely(self) -> None:
        state = "n21_vd_0p05_vg_0p8"
        phumob = self.select(self.crosschecks, state, "term", "phumob")
        hfs = self.select(self.crosschecks, state, "term", "hfs")
        self.assertGreater(abs(float(phumob["effect"])), 1.0)
        self.assertLess(abs(float(hfs["effect"])), 1.0e-6)

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
