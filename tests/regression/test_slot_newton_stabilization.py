import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("toggles", ROOT / "scripts/slot_newton_stabilization.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class StabilizationControls(unittest.TestCase):
    def test_factorial_controls_preserve_acceptance_and_physics(self):
        base = {"physics": {"fermi": True}, "contacts": [{"name": "gate", "voltage": 5}],
                "solver": {"carrier_regularization_scale": 1e-8,
                           "carrier_row_convergence": {"mode": "report", "eps_row": .001},
                           "global_continuity_closure": {"mode": "enforce", "relative_tolerance": .01,
                               "linear_system_stabilization": "replace_free_row_sum"}}}
        original = copy.deepcopy(base)
        for reg in (False, True):
            for closure in (False, True):
                with self.subTest(reg=reg, closure=closure):
                    result = module.configure_stabilization(base, regularization=reg, closure_replacement=closure)
                    self.assertEqual(result["solver"]["carrier_regularization_scale"], 1e-8 if reg else 0)
                    self.assertEqual(result["solver"]["global_continuity_closure"]["linear_system_stabilization"],
                                     "replace_free_row_sum" if closure else "off")
                    result["solver"]["carrier_regularization_scale"] = 1e-8
                    result["solver"]["global_continuity_closure"]["linear_system_stabilization"] = "replace_free_row_sum"
                    self.assertEqual(result, original)
                    self.assertEqual(base, original)

    def test_invalid_regularization_is_rejected(self):
        for scale in (0, -1, float("inf"), float("nan")):
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                module.configure_stabilization({}, regularization=True, closure_replacement=False,
                                               regularization_scale=scale)

    def test_cli_preserves_existing_evidence(self):
        scratch = ROOT / "build"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="stabilization_test_", dir=scratch) as directory:
            source = Path(directory) / "input.json"
            output = Path(directory) / "output.json"
            source.write_text(json.dumps({"solver": {"global_continuity_closure": {"mode": "enforce"}}}), encoding="utf8")
            cmd = [sys.executable, "-B", str(ROOT / "scripts/slot_newton_stabilization.py"),
                   str(source), str(output), "--carrier-regularization", "off", "--closure-replacement", "off"]
            subprocess.run(cmd, check=True, capture_output=True)
            before = output.read_bytes()
            self.assertEqual(json.loads(before)["solver"]["global_continuity_closure"]["mode"], "enforce")
            self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
            self.assertEqual(output.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
