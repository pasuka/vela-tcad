"""Regression checks for frozen SimpleMOS M16 transport-factor attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m16_transport_factor_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM16TransportFactorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"]["m16_transport_factor_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))
        grouped_entry = cls.evidence["artifacts"]["m16_grouped_factor_contributions.csv"]
        with (REPO / grouped_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.grouped = list(csv.DictReader(handle))
        factor_entry = cls.evidence["artifacts"]["m16_factor_contributions.csv"]
        with (REPO / factor_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.factors = list(csv.DictReader(handle))

    @staticmethod
    def select(rows: list[dict[str, str]], state: str, factor: str) -> dict[str, str]:
        return next(row for row in rows
                    if row["state"] == state and row["factor"] == factor)

    def test_scope_and_execution_contract_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"],
                         "complete_with_grouped_kernel_interpretation")
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["factor_probe_count"], 16)
        self.assertEqual(execution["probe_variant_count"], 16)
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_factor_sum_closes_to_m15_electron_transport(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_absolute_A_per_um"], 1.0e-12)
        self.assertLess(closure["maximum_relative"], 1.0e-4)
        self.assertLess(closure["maximum_strong_state_relative"], 1.0e-8)
        self.assertLess(closure["maximum_endpoint_absolute_A_per_um"], 1.0e-15)

    def test_key_state_grouped_kernel_is_dominant(self) -> None:
        state = "n21_vd_0p05_vg_0p8"
        kernel = self.select(self.grouped, state, "sg_state_kernel")
        mobility_state = self.select(self.grouped, state, "mobility_state")
        mobility_drive = self.select(self.grouped, state, "mobility_drive")
        self.assertAlmostEqual(float(kernel["relative_to_baseline_current"]),
                               0.13010103701416173, places=14)
        self.assertAlmostEqual(float(mobility_state["relative_to_baseline_current"]),
                               -0.0005210004482502183, places=14)
        self.assertAlmostEqual(float(mobility_drive["relative_to_baseline_current"]),
                               -0.0025877759640011533, places=14)
        self.assertGreater(float(kernel["absolute_contribution_fraction"]), 0.95)
        self.assertLess(float(mobility_drive["absolute_contribution_fraction"]), 0.02)

    def test_independent_bernoulli_population_split_is_cancellation_dominated(self) -> None:
        state = "n21_vd_0p05_vg_0p8"
        bernoulli = self.select(self.factors, state, "bernoulli_weights")
        population = self.select(self.factors, state, "carrier_population")
        self.assertLess(float(bernoulli["relative_to_baseline_current"]), -20.0)
        self.assertGreater(float(population["relative_to_baseline_current"]), 20.0)

    def test_grouped_kernel_dominates_every_strong_state(self) -> None:
        strong_kernel = [float(row["absolute_contribution_fraction"])
                         for row in self.grouped
                         if row["factor"] == "sg_state_kernel"
                         and float(row["gate_voltage_V"]) >= 0.8]
        self.assertEqual(len(strong_kernel), 8)
        self.assertGreater(min(strong_kernel), 0.80)

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
