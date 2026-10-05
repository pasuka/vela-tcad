from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "schottky_reference_workflow",
    ROOT / "scripts" / "run_schottky_reference_workflow.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

from run_schottky_sentaurus_vm import (  # noqa: E402
    prepare_bundle as prepare_vm_bundle,
    remote_commands as vm_remote_commands,
)


class SchottkyReferenceWorkflowTest(unittest.TestCase):
    def test_checked_in_curve_passes_declared_gate(self) -> None:
        rows = [
            {"bias_V": float(row["bias_V"]),
             "current_total_A_per_um": float(row["current_total_A_per_um"])}
            for row in MODULE.read_csv(
                MODULE.FIXTURE / "vela_schottky_iv_combined.csv")
        ]
        result = MODULE.compare_curves(rows)
        self.assertEqual("pass", result["status"])
        self.assertEqual(24, result["points_compared"])
        self.assertLessEqual(result["maximum_log10_current_error_dex"], 0.5)

    def test_materialized_second_stage_reuses_exact_0p82_state(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first, second = MODULE.materialize(Path(directory))
            for path in (first, second):
                config = json.loads(path.read_text())
                material = Path(config["materials_file"])
                self.assertTrue(material.is_absolute())
                self.assertTrue(material.is_file())
                self.assertEqual("materials_legacy_builtin_300k.json", material.name)
            self.assertTrue(first.is_file())
            config = __import__("json").loads(second.read_text(encoding="utf-8"))
            self.assertTrue(config["sweep"]["continuation"]["arclength"]["enabled"])
            self.assertTrue(config["sweep"]["initial_state_file"].endswith(
                "schottky_bias_0p820000.h5"))
            self.assertEqual(
                "report", config["solver"]["global_continuity_closure"]["mode"])

    def test_sentaurus_vm_bundle_preserves_bounded_physics_source(self) -> None:
        source = ROOT / "reference_tcad" / "schottky_charon_sentaurus2018" / "source"
        with tempfile.TemporaryDirectory() as directory:
            bundle = Path(directory) / "bundle"
            prepare_vm_bundle(source, bundle)
            self.assertEqual(
                {"schottky_n_sde.cmd", "schottky_n_des.cmd"},
                {path.name for path in bundle.iterdir()},
            )
            deck = (bundle / "schottky_n_des.cmd").read_text()
            self.assertIn("Schottky Workfunction=4.75", deck)
            self.assertNotIn("BarrierLowering", deck)
            self.assertNotIn("Tunneling", deck)

        commands = vm_remote_commands("/tmp/schottky")
        self.assertEqual(3, len(commands))
        self.assertIn("sde -e -l", commands[0])
        self.assertIn("sdevice schottky_n_des.cmd", commands[1])


if __name__ == "__main__":
    unittest.main()
