#!/usr/bin/env python3
"""Gated integration coverage for a local Sentaurus sample project."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SAMPLE_ENV = "VELA_SENTAURUS_SAMPLE_DIR"
BUILD_ENV = "VELA_BUILD_DIR"


class SentaurusSampleIntegrationTest(unittest.TestCase):
    def test_ldmos_n20_sample_inventory_curve_and_cmd_when_enabled(self) -> None:
        sample_root = self._sample_root_or_skip()
        build_root = Path(os.environ.get(BUILD_ENV, REPO / "build"))
        importer = build_root / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")
        if not importer.is_file():
            self.skipTest(f"Sentaurus HDF5 importer is not built: {importer}")

        tdr = sample_root / "n20_des.tdr"
        plt = sample_root / "IdVd_n20_des.plt"
        cmd = sample_root / "pp20_des.cmd"
        for path in [tdr, plt, cmd]:
            self.assertTrue(path.is_file(), f"missing sample artifact: {path}")

        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_sample_") as tmp:
            out = Path(tmp)
            inventory_path = out / "n20_inventory.json"
            curve_path = out / "ldmos_idvd_reference.csv"
            cmd_summary_path = out / "pp20_summary.json"
            runner_path = out / "run_idvd.py"

            subprocess.run(
                [
                    str(importer),
                    "--tdr",
                    str(tdr),
                    "--inventory-json",
                    str(inventory_path),
                ],
                check=True,
                cwd=REPO,
            )
            inventory = json.loads(inventory_path.read_text())
            self.assertEqual(inventory["vertex_count"], 23081)
            self.assertEqual(inventory["region_count"], 7)
            self.assertEqual(inventory["dataset_count"], 63)
            self.assertEqual(
                {region["name"] for region in inventory["regions"] if region["type"] == 1},
                {"gate", "drain", "source", "substrate"},
            )
            self._assert_contact_scalar(inventory, "ContactExternalVoltage", "gate", 10.0)
            self._assert_contact_scalar(inventory, "ContactExternalVoltage", "drain", 30.0)
            self._assert_contact_scalar(inventory, "ContactExternalVoltage", "source", 0.0)
            self._assert_contact_scalar(inventory, "ContactExternalVoltage", "substrate", 0.0)
            self._assert_contact_scalar(inventory, "ContactCurrentFlux", "drain", 0.0021962326764188673)

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "plt",
                    "--input",
                    str(plt),
                    "--output",
                    str(curve_path),
                    "--bias-column",
                    "drain OuterVoltage",
                    "--current-column",
                    "drain TotalCurrent",
                ],
                check=True,
                cwd=REPO,
            )
            with curve_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 81)
            self.assertAlmostEqual(float(rows[-1]["bias_V"]), 30.0, places=12)
            self.assertAlmostEqual(float(rows[-1]["current_total"]), 0.00219623267794, places=14)

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "cmd",
                    "--input",
                    str(cmd),
                    "--summary-json",
                    str(cmd_summary_path),
                    "--python-runner",
                    str(runner_path),
                    "--mesh-json",
                    "vela/mesh.json",
                    "--output-csv",
                    "outputs/ldmos_idvd.csv",
                ],
                check=True,
                cwd=REPO,
            )
            summary = json.loads(cmd_summary_path.read_text())
            self.assertEqual(summary["files"]["grid"], "n13_final_fps.tdr")
            self.assertEqual(summary["files"]["current"], "n20_des.plt")
            self.assertEqual(
                [(sweep["contact"], sweep["stop"]) for sweep in summary["sweeps"]],
                [("gate", 10.0), ("drain", 30.0)],
            )
            self.assertIn("Thermodynamic", summary["unsupported_physics"])
            self.assertIn("from vela.curves import run_iv_curve", runner_path.read_text())

    def _sample_root_or_skip(self) -> Path:
        value = os.environ.get(SAMPLE_ENV)
        if not value:
            self.skipTest(f"{SAMPLE_ENV} is not set")
        root = Path(value)
        self.assertTrue(root.is_dir(), f"{SAMPLE_ENV} is not a directory: {root}")
        return root

    def _assert_contact_scalar(self,
                               inventory: dict[str, object],
                               field_name: str,
                               contact_name: str,
                               expected: float) -> None:
        region_names = {
            int(region["index"]): str(region["name"])
            for region in inventory["regions"]  # type: ignore[index]
        }
        matches = [
            field for field in inventory["fields"]  # type: ignore[index]
            if field["name"] == field_name and region_names[int(field["region"])] == contact_name
        ]
        self.assertEqual(len(matches), 1, f"{field_name} for {contact_name}")
        values = matches[0].get("raw_values", [])
        self.assertEqual(len(values), 1, f"{field_name} for {contact_name}")
        self.assertAlmostEqual(float(values[0]), expected, places=12)


if __name__ == "__main__":
    unittest.main()
