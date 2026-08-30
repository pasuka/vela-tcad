"""Regression checks for frozen SimpleMOS M17 SG secant attribution."""

from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m17_sg_secant_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM17SgSecantTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        report_entry = cls.evidence["artifacts"]["m17_sg_secant_report.json"]
        cls.report = json.loads((REPO / report_entry["path"])
                                .read_text(encoding="utf-8"))
        factor_entry = cls.evidence["artifacts"]["m17_factor_contributions.csv"]
        with (REPO / factor_entry["path"]).open(newline="", encoding="utf-8") as handle:
            cls.factors = list(csv.DictReader(handle))

    @staticmethod
    def select(rows: list[dict[str, str]], state: str, factor: str) -> dict[str, str]:
        return next(row for row in rows
                    if row["state"] == state and row["factor"] == factor)

    def test_scope_and_execution_contract_are_frozen(self) -> None:
        self.assertEqual(self.evidence["status"],
                         "complete_with_bias_regime_separation")
        execution = self.report["execution"]
        self.assertEqual(execution["state_count"], 16)
        self.assertEqual(execution["factor_probe_count"], 16)
        self.assertEqual(execution["probe_variant_count"], 8)
        self.assertFalse(execution["new_sentaurus_execution"])
        self.assertFalse(execution["default_model_changed"])

    def test_stable_factor_sum_closes_to_m15_transport(self) -> None:
        closure = self.report["closure"]
        self.assertLess(closure["maximum_absolute_A_per_um"], 1.0e-12)
        self.assertLess(closure["maximum_relative"], 1.0e-4)
        self.assertLess(closure["maximum_strong_state_relative"], 1.0e-8)
        self.assertEqual(closure["maximum_endpoint_absolute_A_per_um"], 0.0)

    def test_key_state_qf_imbalance_is_dominant(self) -> None:
        state = "n21_vd_0p05_vg_0p8"
        mobility = self.select(self.factors, state, "mobility")
        conductance = self.select(self.factors, state, "sg_secant_conductance")
        imbalance = self.select(self.factors, state, "qf_log_imbalance")
        self.assertAlmostEqual(float(mobility["relative_to_baseline_current"]),
                               -0.003108508614059112, places=14)
        self.assertAlmostEqual(float(conductance["relative_to_baseline_current"]),
                               -0.007413541243536472, places=14)
        self.assertAlmostEqual(float(imbalance["relative_to_baseline_current"]),
                               0.1375143104595058, places=14)
        self.assertGreater(float(imbalance["absolute_contribution_fraction"]), 0.90)

    def test_stable_parameterization_removes_m16_cancellation(self) -> None:
        cancellation = self.report["key_state"]["cancellation_amplification"]
        self.assertGreater(cancellation["m16_ungrouped_four_factor"], 400.0)
        self.assertLess(cancellation["m17_stable_three_factor"], 1.2)
        self.assertGreater(cancellation["reduction_factor"], 300.0)

    def test_on_state_is_conductance_dominated(self) -> None:
        on_state = [float(row["absolute_contribution_fraction"])
                    for row in self.factors
                    if row["factor"] == "sg_secant_conductance"
                    and float(row["gate_voltage_V"]) == 2.5]
        self.assertEqual(len(on_state), 4)
        self.assertGreater(min(on_state), 0.75)

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
