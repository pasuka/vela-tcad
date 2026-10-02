import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "tight_convergence_port_burst"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM60TightConvergencePortBurstTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.freeze = read_json(
            ROOT / "simplemos_m60_tight_convergence_port_burst_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m60_tight_convergence_port_burst_evidence.json")
        cls.report = read_json(OUTPUT / "m60_tight_convergence_port_burst_report.json")
        cls.points = rows(OUTPUT / "m60_tight_default_direct_point_ledger.csv")
        cls.terminals = rows(OUTPUT / "m60_tight_terminal_component_ledger.csv")
        cls.cases = rows(OUTPUT / "m60_tight_case_summary.csv")
        cls.logs = rows(OUTPUT / "m60_cnormprint_log_ledger.csv")
        cls.fields = rows(OUTPUT / "m60_default_direct_field_invariance_ledger.csv")

    def test_contract_and_upstream_hashes_are_frozen(self):
        contract = ROOT / "simplemos_m60_tight_convergence_port_burst_contract_v1.json"
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(self.freeze["contract_sha256"], sha256(contract))
        for relative, expected in self.freeze["upstream_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)

    def test_full_paired_matrix_is_complete(self):
        self.assertEqual(len(self.points), 816)
        self.assertEqual(len(self.terminals), 19584)
        self.assertEqual(len(self.cases), 16)
        self.assertEqual(len(self.logs), 32)
        self.assertEqual({int(row["point_count"]) for row in self.cases}, {51})

    def test_tight_convergence_suppresses_all_frozen_bursts(self):
        self.assertEqual(sum(int(row["m59_burst_flag"]) for row in self.points), 15)
        self.assertEqual(sum(int(row["tight_burst_flag"]) for row in self.points), 0)
        self.assertLessEqual(
            self.report["matrix"]["maximum_nonburst_default_vs_m46_log_shift_dex"],
            0.005)

    def test_n23_target_passes_suppression_gates(self):
        target = self.report["target"]
        self.assertAlmostEqual(float(target["tight_default_vela_error_dex"]),
                               0.024398282551395253, delta=1e-15)
        self.assertGreaterEqual(float(target["observable_gap_reduction_fraction"]),
                                0.90)
        self.assertLess(abs(float(target["tight_signed_observable_fraction"])), 0.02)

    def test_default_and_direct_observers_preserve_solved_fields(self):
        self.assertTrue(self.report["state_invariance"]["all_fields_invariant"])
        self.assertEqual(self.report["state_invariance"]["failure_count"], 0)
        self.assertGreater(len(self.fields), 0)
        self.assertTrue(all(row["within_state_invariance"] == "True"
                            for row in self.fields))

    def test_cnormprint_logs_and_secondary_policy_are_retained(self):
        self.assertTrue(all(int(row["console_bytes"]) > 1000 for row in self.logs))
        self.assertTrue(all(int(row["completed_without_exit_failure"]) == 1
                            for row in self.logs))
        qualification = read_json(
            OUTPUT / "m60_sub_fa_secondary_reference_qualification.json")
        self.assertTrue(qualification["eligible"])
        self.assertFalse(qualification["replaces_m8_m46"])

    def test_evidence_is_frozen_and_complete(self):
        self.assertEqual(self.report["status"], "accepted")
        self.assertEqual(self.report["classification"],
                         "tight_convergence_suppresses_port_bursts")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertEqual(self.evidence["status"], "frozen")
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
