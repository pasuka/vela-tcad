import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "port_observable_decomposition"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM59PortObservableDecompositionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = read_json(
            ROOT / "simplemos_m59_port_observable_decomposition_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m59_port_observable_decomposition_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m59_port_observable_decomposition_evidence.json")
        cls.report = read_json(
            OUTPUT / "m59_port_observable_decomposition_report.json")
        cls.points = rows(OUTPUT / "m59_pointwise_observable_ledger.csv")
        cls.curves = rows(OUTPUT / "m59_curve_decomposition_ledger.csv")
        cls.pairs = rows(OUTPUT / "m59_pair_decomposition_ledger.csv")
        cls.sensitivity = rows(
            OUTPUT / "m59_threshold_scope_sensitivity_ledger.csv")
        cls.kcl = rows(OUTPUT / "m59_kcl_topology_stratification_ledger.csv")

    def test_contract_and_upstream_hashes_are_frozen(self):
        contract_path = ROOT / "simplemos_m59_port_observable_decomposition_contract_v1.json"
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(self.freeze["contract_sha256"], sha256(contract_path))
        for relative, expected in self.freeze["upstream_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)

    def test_matrix_and_ledgers_are_complete(self):
        self.assertEqual(len(self.points), 816)
        self.assertEqual(len(self.curves), 16)
        self.assertEqual(len(self.pairs), 8)
        self.assertEqual(len(self.kcl), 16)
        self.assertEqual({int(row["primary_burst_flag"]) for row in self.points}, {0, 1})

    def test_direct_observable_is_not_named_as_independent_residual(self):
        headers = set(self.points[0])
        self.assertIn("substrate_default_minus_direct_eCurrent_A_per_um", headers)
        self.assertNotIn("pwell_residual_A_per_um", headers)
        self.assertIn("residual causality remains an M60 hypothesis",
                      self.report["port_observable_component"]["interpretation"])

    def test_burst_free_growth_is_threshold_robust(self):
        selected = [row for row in self.sensitivity
                    if row["scope"] == "all" and float(row["threshold"]) <= 0.03]
        self.assertEqual(len(selected), 5)
        for row in selected:
            self.assertGreaterEqual(
                float(row["minimum_pair_burst_free_growth_dex"]), 0.020)
            self.assertLessEqual(
                float(row["maximum_pair_burst_free_growth_dex"]), 0.027)

    def test_low_and_high_drain_slopes_are_not_conflated(self):
        primary = [row for row in self.sensitivity
                   if abs(float(row["threshold"]) - 0.02) < 1e-12]
        low = next(row for row in primary if row["scope"] == "low_drain")
        high = next(row for row in primary if row["scope"] == "high_drain")
        self.assertAlmostEqual(float(low["zero_intercept_slope"]), 0.370, delta=0.001)
        self.assertGreaterEqual(float(low["pearson"]), 0.99)
        self.assertLess(float(high["zero_intercept_slope"]), 0.0)

    def test_floating_topology_stratifies_default_kcl(self):
        floating = [row for row in self.kcl
                    if int(row["single_node_floating_n_type_component_count"]) > 0]
        nonfloating = [row for row in self.kcl
                       if int(row["single_node_floating_n_type_component_count"]) == 0]
        self.assertEqual({row["device"] for row in floating}, {"n21", "n23"})
        ratio = min(float(row["median_absolute_default_four_terminal_kcl_A_per_um"])
                    for row in floating) / max(
                        float(row["median_absolute_default_four_terminal_kcl_A_per_um"])
                        for row in nonfloating)
        self.assertGreaterEqual(ratio, 100.0)

    def test_evidence_and_acceptance_are_complete(self):
        self.assertEqual(self.report["status"], "accepted")
        self.assertEqual(
            self.report["classification"],
            "port_observable_burst_association_plus_smooth_nwell_growth")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertEqual(self.evidence["status"], "frozen")
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
