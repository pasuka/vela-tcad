#!/usr/bin/env python3
"""Regression coverage for Sentaurus import contact overrides."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


class SentaurusContactOverridesTest(unittest.TestCase):
    def test_reference_import_applies_vela_contact_overrides(self) -> None:
        module_path = REPO / "scripts" / "sentaurus_import.py"
        spec = importlib.util.spec_from_file_location("sentaurus_import", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_contact_overrides_") as td:
            tmp = Path(td)
            deck_path = tmp / "simulation_idvd.json"
            deck_path.write_text(json.dumps({
                "node_doping_file": "doping.csv",
                "contacts": [
                    {"name": "Source", "bias": 0.0},
                    {"name": "Drain", "bias": 0.0},
                    {"name": "Gate", "bias": 0.0},
                    {"name": "Body", "bias": 0.0},
                ],
                "solver": {"method": "gummel_newton"},
                "sweep": {},
            }), encoding="utf-8")
            cmd_summary = {
                "physics": [],
                "sweeps": [{
                    "contact": "Drain",
                    "stop": 0.5,
                    "step_control": {"MaxStep": 0.02},
                }],
            }
            sim = {
                "name": "idvd",
                "kind": "iv",
                "vela_current_contact": "Drain",
                "vela_node_doping_file": False,
                "vela_contact_overrides": [
                    {"name": "Gate", "type": "metal_gate", "bias": 2.0, "flatband_voltage": 0.0}
                ],
            }

            module.patch_reference_deck(deck_path, cmd_summary, sim, "candidate.csv")
            patched = json.loads(deck_path.read_text(encoding="utf-8"))

        by_name = {item["name"]: item for item in patched["contacts"]}
        self.assertEqual(by_name["Gate"]["type"], "metal_gate")
        self.assertAlmostEqual(by_name["Gate"]["bias"], 2.0)
        self.assertAlmostEqual(by_name["Gate"]["flatband_voltage"], 0.0)
        self.assertNotIn("type", by_name["Drain"])
        self.assertEqual(patched["sweep"]["contact"], "Drain")
        self.assertEqual(patched["sweep"]["current_contact"], "Drain")
        self.assertNotIn("node_doping_file", patched)


if __name__ == "__main__":
    unittest.main()
