from __future__ import annotations

import importlib.util
import tempfile
import unittest
import csv
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "run_slot_ldmos_noavalanche_continuation.py"
SPEC = importlib.util.spec_from_file_location("slot_adaptive", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SlotAdaptiveContinuationTest(unittest.TestCase):
    def test_attempt_render_keeps_bundle_relative_outputs_and_absolute_restart(self) -> None:
        base = {
            "output_csv": "old.csv",
            "sweep": {
                "bias_points": [60.0],
                "initial_state_file": "old_state.csv",
                "write_state_every_point_prefix": "old/state",
                "write_state_file": "old/final.csv",
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp)
            initial = bundle / "seed.csv"
            document, config, run_dir = MODULE.render_attempt(
                {**base, "solver": {
                    "max_iter": 80,
                    "handoff": {"newton_max_iter": 80},
                    "carrier_row_convergence": {"min_newton_max_iter": 80},
                }}, bundle, 0.2, initial, 3, 8, 1.0e-6,
                "uniform_trust_region", "block_filter", "gummel_newton", 11)
            self.assertEqual(document["sweep"]["bias_points"], [0.2])
            self.assertEqual(document["sweep"]["initial_state_file"], str(initial.resolve()))
            self.assertTrue(document["output_csv"].endswith("iv.csv"))
            self.assertEqual(config.parent, bundle)
            self.assertIn("attempt_003_v_0p2", str(run_dir))
            self.assertEqual(document["solver"]["max_iter"], 8)
            self.assertEqual(document["solver"]["stall_residual_floor"], 1.0e-6)
            self.assertEqual(
                document["solver"]["quasi_fermi_update_limit_mode"],
                "uniform_trust_region")
            self.assertEqual(document["solver"]["line_search_mode"], "block_filter")
            self.assertEqual(document["solver"]["method"], "gummel_newton")
            self.assertEqual(document["solver"]["handoff"]["gummel_max_iter"], 11)
            self.assertEqual(document["solver"]["handoff"]["newton_max_iter"], 8)
            self.assertEqual(
                document["solver"]["carrier_row_convergence"]["min_newton_max_iter"], 8)
            self.assertTrue(
                document["solver"]["carrier_row_convergence"]["diagnostic_csv"]
                .endswith("carrier_row_convergence.csv"))

    def test_poisson_prewarm_preserves_qf_and_writes_separate_state(self) -> None:
        base = {
            "output_csv": "old.csv",
            "solver": {
                "method": "newton",
                "max_iter": 80,
                "handoff": {"newton_max_iter": 80},
            },
            "sweep": {
                "bias_points": [60.0],
                "initial_state_file": "old_state.csv",
                "write_state_every_point_prefix": "old/state",
                "write_state_file": "old/final.csv",
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp)
            initial = bundle / "seed.csv"
            document, config, run_dir = MODULE.render_poisson_prewarm(
                base, bundle, 0.246875, initial, 8, 17)
            self.assertEqual(document["solver"]["method"], "poisson_only")
            self.assertEqual(
                document["solver"]["poisson_only_quasi_fermi_initialization"],
                "preserve_initial")
            self.assertEqual(document["solver"]["max_iter"], 17)
            self.assertEqual(document["sweep"]["bias_points"], [0.246875])
            self.assertEqual(
                document["sweep"]["initial_state_file"], str(initial.resolve()))
            self.assertIn("prewarm_008_v_0p246875", str(run_dir))
            self.assertIn("prewarm_008", config.name)

    def test_batch_predictor_keeps_secant_history_in_one_sweep(self) -> None:
        base = {
            "output_csv": "old.csv",
            "solver": {"max_iter": 80},
            "sweep": {
                "bias_points": [60.0],
                "initial_state_file": "old.csv",
                "write_state_every_point_prefix": "old/state",
                "write_state_file": "old/final.csv",
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp)
            initial = bundle / "seed.csv"
            document, config, run_dir = MODULE.render_batch_predictor(
                base, bundle, [0.2, 0.35, 0.5], initial, 12, 1.0e-6,
                "componentwise_poisson_recorrect", "merit")
            self.assertEqual(document["sweep"]["bias_points"], [0.2, 0.35, 0.5])
            self.assertEqual(
                document["sweep"]["continuation"]["predictor"]["mode"], "secant")
            self.assertEqual(
                document["sweep"]["continuation"]["predictor"]["fields"],
                ["psi", "phin", "phip"])
            self.assertIn("batch_predictor_v_0p5", str(run_dir))
            self.assertEqual(config.name, "adaptive_noavalanche_batch_predictor.json")

    def test_secant_restart_extrapolates_only_coupled_unknowns(self) -> None:
        fields = [
            "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3",
            "electron_qf_increment_V", "hole_qf_increment_V",
            "electron_qf_reference_V", "hole_qf_reference_V"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous = root / "previous.csv"
            current = root / "current.csv"
            output = root / "predicted.csv"
            rows = [
                ["0", "1", "2", "3", "10", "20", "2", "3", "0", "0"],
                ["0", "2", "4", "6", "11", "21", "4", "6", "0", "0"],
            ]
            for path, row in ((previous, rows[0]), (current, rows[1])):
                with path.open("w", newline="", encoding="utf-8") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(fields)
                    writer.writerow(row)
            ratio = MODULE.write_secant_restart(
                previous, current, output, 0.2, 0.35, 0.5)
            self.assertAlmostEqual(ratio, 1.0)
            with output.open(newline="", encoding="utf-8") as stream:
                predicted = next(csv.DictReader(stream))
            self.assertAlmostEqual(float(predicted["psi"]), 3.0)
            self.assertAlmostEqual(float(predicted["phin"]), 6.0)
            self.assertAlmostEqual(float(predicted["phip"]), 9.0)
            self.assertEqual(float(predicted["electrons_m3"]), 11.0)
            self.assertEqual(float(predicted["holes_m3"]), 21.0)
            self.assertAlmostEqual(float(predicted["electron_qf_increment_V"]), 6.0)
            self.assertAlmostEqual(float(predicted["hole_qf_increment_V"]), 9.0)


if __name__ == "__main__":
    unittest.main()
