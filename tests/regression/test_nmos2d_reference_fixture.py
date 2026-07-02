#!/usr/bin/env python3
"""Regression coverage for the NMOS Sentaurus Id-Vd reference fixture."""

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


class Nmos2dReferenceFixtureTest(unittest.TestCase):
    def test_nmos2d_sentaurus2018_idvd_reference_config(self) -> None:
        path = REPO / "reference_tcad" / "nmos2d_sentaurus2018" / "nmos2d_sentaurus2018_reference.json"
        self.assertTrue(path.is_file(), f"missing reference config: {path}")
        config = json.loads(path.read_text())

        self.assertEqual(config["case"], "nmos2d_sentaurus2018")
        self.assertEqual(config["device"], "nmos2d")
        self.assertEqual(config["mesh_tdr"], "nmos2d_msh.tdr")
        self.assertEqual(config["sde_cmd"], "nmos2d_sde.cmd")

        iv = next(sim for sim in config["simulations"] if sim["name"] == "idvd")
        self.assertEqual(iv["kind"], "iv")
        self.assertEqual(iv["cmd"], "nmos2d_idvd_sdevice.cmd")
        self.assertEqual(iv["plt"], "nmos2d_idvd.plt")
        self.assertEqual(iv["bias_column"], "Drain OuterVoltage")
        self.assertEqual(iv["current_column"], "Drain TotalCurrent")
        self.assertEqual(iv["vela_current_contact"], "Drain")
        self.assertAlmostEqual(iv["vela_stop"], 0.5)
        self.assertAlmostEqual(iv["vela_step"], 0.02)
        self.assertEqual(iv["comparison"]["candidate_column"], "current_total_A_per_um")
        self.assertEqual(iv["comparison"]["interpolation"], "log_current")

        gate_override = next(
            item for item in iv["vela_contact_overrides"] if item["name"] == "Gate"
        )
        self.assertEqual(gate_override["type"], "metal_gate")
        self.assertAlmostEqual(gate_override["bias"], 2.0)
        self.assertAlmostEqual(gate_override["flatband_voltage"], 0.0)
        self.assertIs(iv["vela_node_doping_file"], False)

    def test_nmos2d_sentaurus2018_source_decks_define_rectangular_idvd_case(self) -> None:
        source = REPO / "reference_tcad" / "nmos2d_sentaurus2018" / "source"
        sde = (source / "nmos2d_sde.cmd").read_text()
        sdevice = (source / "nmos2d_idvd_sdevice.cmd").read_text()

        for token in [
            '"R.Si"',
            '"R.Ox"',
            '"Source"',
            '"Drain"',
            '"Gate"',
            '"Body"',
            '"P.Body.Doping"',
            '"N.Source.Doping"',
            '"N.Drain.Doping"',
        ]:
            self.assertIn(token, sde)

        self.assertIn('Grid    = "nmos2d_msh.tdr"', sdevice)
        self.assertIn('Current = "nmos2d_idvd.plt"', sdevice)
        self.assertIn('{ Name="Gate" Voltage=2.0 }', sdevice)
        self.assertIn('Name="Drain"', sdevice)
        self.assertIn("Voltage=0.5", sdevice)
        self.assertIn("Mobility(\n    DopingDependence", sdevice)
        self.assertIn("Recombination(\n    SRH", sdevice)
        self.assertIn("EffectiveIntrinsicDensity(\n    OldSlotboom", sdevice)

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

    def test_nmos2d_sentaurus2018_source_dir_keeps_only_cmd_sources(self) -> None:
        source = REPO / "reference_tcad" / "nmos2d_sentaurus2018" / "source"
        self.assertEqual(
            sorted(item.name for item in source.iterdir()),
            ["nmos2d_idvd_sdevice.cmd", "nmos2d_sde.cmd"],
        )

if __name__ == "__main__":
    unittest.main()
