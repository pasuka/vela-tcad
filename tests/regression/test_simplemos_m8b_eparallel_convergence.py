from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "simplemos_m8b_eparallel_convergence",
    ROOT / "scripts" / "run_simplemos_m8b_eparallel_convergence.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8b_eparallel_convergence_contract_v1.json"
)


class SimpleMosM8BEparallelConvergenceTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = MODULE.read_json(CONTRACT_PATH)
        MODULE.validate_contract(self.contract)
        self.variant = {
            "id": "eparallel",
            "hfs_option": "HighFieldSaturation(Eparallel)",
            "math_options": [],
        }

    def test_contract_freezes_minimal_two_device_probe(self) -> None:
        self.assertEqual(["n17", "n21"], self.contract["devices"])
        self.assertEqual([0.5, 0.1, 0.01], [
            item["line_search_damping"]
            for item in self.contract["probes"]])
        self.assertEqual(51, self.contract["gate_lattice"]["point_count"])
        self.assertEqual("forbidden", self.contract[
            "gate_lattice"]["interpolation"])

    def test_damping_changes_gate_solver_not_physics_or_output_lattice(self) -> None:
        deck = MODULE.damped_deck(
            "fixture", self.variant, 0.1,
            self.contract["common_numerics"])
        self.assertIn("HighFieldSaturation(Eparallel)", deck)
        self.assertIn(
            "Coupled(Iterations=100 LineSearchDamping=0.1) "
            "{ Poisson Electron Hole }", deck)
        self.assertEqual(1, deck.count("LineSearchDamping=0.1"))
        self.assertIn("CurrentPlot(Time=(Range=(0 1) Intervals=50))", deck)
        self.assertIn("MinStep=1e-08", deck)

    def test_parallel_executor_rejects_zero_workers(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one"):
            MODULE.execute_remote(
                {}, Path("unused"), "host", "ssh", "scp", "remote", 0,
                Path("unused_raw"))


if __name__ == "__main__":
    unittest.main()
