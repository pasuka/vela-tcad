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
PN2D_ENV = "VELA_SENTAURUS_PN2D_DIR"
BUILD_ENV = "VELA_BUILD_DIR"


class SentaurusSampleIntegrationTest(unittest.TestCase):
    def test_pn2d_reference_import_when_enabled(self) -> None:
        pn2d_root = self._pn2d_root_or_skip()
        config_path = pn2d_root / "pn2d_reference.json"
        self.assertTrue(config_path.is_file(), f"missing pn2d reference config: {config_path}")
        build_root = Path(os.environ.get(BUILD_ENV, REPO / "build"))
        importer = build_root / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")
        if not importer.is_file():
            self.skipTest(f"Sentaurus HDF5 importer is not built: {importer}")

        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_pn2d_") as tmp:
            out = Path(tmp)

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "reference",
                    "--config",
                    str(config_path),
                    "--source-dir",
                    str(pn2d_root),
                    "--output-dir",
                    str(out / "reference"),
                    "--tdr-importer",
                    str(importer),
                ],
                check=True,
                cwd=REPO,
            )

            self._assert_doping_csv_has_donors_and_acceptors(out / "reference" / "doping.csv")

            summary = json.loads((out / "reference" / "sde_summary.json").read_text())
            self.assertEqual(summary["defines"]["L"], 2.0)
            self.assertEqual(summary["defines"]["H"], 0.5)
            self.assertEqual(summary["defines"]["Xj"], 1.0)
            self.assertEqual(summary["defines"]["Na"], 1.0e17)
            self.assertEqual(summary["defines"]["Nd"], 1.0e17)

            inventory = json.loads((out / "reference" / "tdr_inventory" / "iv.json").read_text())
            self.assertEqual(
                {region["name"] for region in inventory["regions"] if region["type"] == 1},
                {"Anode", "Cathode"},
            )

            iv_rows = self._read_curve(out / "reference" / "reference_curves" / "pn2d_iv_reference.csv")
            bv_rows = self._read_curve(out / "reference" / "reference_curves" / "pn2d_bv_reference.csv")
            self.assertEqual(len(iv_rows), 62)
            self.assertEqual(len(bv_rows), 39)
            self.assertAlmostEqual(float(iv_rows[-1]["bias_V"]), 1.0, places=12)
            self.assertAlmostEqual(float(iv_rows[-1]["current_total"]), 1.1637177246e-4, places=14)
            self.assertAlmostEqual(float(bv_rows[-1]["bias_V"]), 50.0, places=12)

            iv_deck = json.loads((out / "reference" / "vela" / "simulation_iv.json").read_text())
            bv_deck = json.loads((out / "reference" / "vela" / "simulation_bv.json").read_text())
            self.assertTrue((out / "reference" / "vela" / "simulation_iv.json").is_file())
            self.assertTrue((out / "reference" / "vela" / "simulation_bv.json").is_file())
            self.assertFalse((out / "reference" / "vela" / "simulation_iv_runtime.json").exists())
            self.assertFalse((out / "reference" / "vela" / "simulation_bv_runtime.json").exists())
            self.assertTrue((out / "reference" / "reference_tcad_manifest.json").is_file())
            self.assertEqual(iv_deck["sweep"]["contact"], "Anode")
            self.assertEqual(iv_deck["sweep"]["stop"], 0.3)
            self.assertEqual(iv_deck["node_doping_file"], "doping.csv")
            self.assertEqual(iv_deck["solver"]["method"], "gummel_newton")
            self.assertTrue(iv_deck["solver"]["warm_start"])
            self.assertEqual(iv_deck["solver"]["handoff"]["fallback"], "none")
            self.assertFalse(iv_deck["solver"]["handoff"]["require_gummel_convergence"])
            self.assertEqual(iv_deck["solver"]["handoff"]["gummel_max_iter"], 0)
            self.assertGreater(iv_deck["solver"]["handoff"]["newton_max_iter"], 0)
            self.assertEqual(bv_deck["sweep"]["contact"], "Cathode")
            self.assertEqual(bv_deck["sweep"]["stop"], 0.05)
            self.assertEqual(bv_deck["node_doping_file"], "doping.csv")
            self.assertEqual(bv_deck["solver"]["method"], "gummel_newton")
            self.assertEqual(bv_deck["solver"]["recombination"], ["srh"])
            self.assertAlmostEqual(bv_deck["solver"]["taun"], 1.0e-6)
            self.assertAlmostEqual(bv_deck["solver"]["taup"], 1.0e-6)
            self.assertEqual(bv_deck["solver"]["bandgap_narrowing"], "none")
            self.assertEqual(bv_deck["solver"]["mobility"]["model"], "caughey_thomas")
            self.assertAlmostEqual(
                bv_deck["solver"]["mobility"]["electron_mu_min_m2_V_s"], 46.458)
            self.assertAlmostEqual(
                bv_deck["solver"]["mobility"]["hole_mu_min_m2_V_s"], 39.961)
            self.assertAlmostEqual(
                bv_deck["solver"]["mobility"]["electron_alpha"], 0.6052)
            self.assertAlmostEqual(
                bv_deck["solver"]["mobility"]["hole_alpha"], 0.623)

            faithful_iv = out / "reference" / "vela" / "pn2d_iv.csv"
            faithful_bv = out / "reference" / "vela" / "pn2d_bv.csv"
            self.assertTrue(faithful_iv.is_file())
            self.assertTrue(faithful_bv.is_file())
            self.assertGreaterEqual(len(self._read_curve(faithful_iv)), 2)
            self.assertGreaterEqual(len(self._read_curve(faithful_bv)), 2)
            self._assert_curve_has_finite_currents(faithful_iv)
            self._assert_curve_has_finite_currents(faithful_bv)
            faithful_iv_rows = self._read_curve(faithful_iv)
            faithful_bv_rows = self._read_curve(faithful_bv)
            self.assertLess(abs(float(faithful_iv_rows[0]["current_total"])), 1.0e-9)
            self.assertLess(abs(float(faithful_bv_rows[0]["current_total"])), 1.0e-9)
            for row in faithful_iv_rows:
                self.assertEqual(row["solver_method"], "gummel_newton")
                self.assertEqual(row["handoff_stage"], "newton")
                self.assertGreater(int(row["newton_iterations"]), 0)
            for row in faithful_bv_rows:
                self.assertEqual(row["solver_method"], "gummel_newton")
                self.assertEqual(row["handoff_stage"], "newton")
                self.assertGreater(int(row["newton_iterations"]), 0)

            manifest = json.loads((out / "reference" / "reference_tcad_manifest.json").read_text())
            self.assertFalse(manifest["commit_policy"]["raw_sentaurus_artifacts"])
            self.assertIn("Avalanche", manifest["unsupported_physics"])
            self.assertIn("reports/pn2d_iv_comparison.json", manifest["comparison_reports"])
            self.assertIn("reports/pn2d_bv_comparison.json", manifest["comparison_reports"])
            iv_report = json.loads((out / "reference" / "reports" / "pn2d_iv_comparison.json").read_text())
            bv_report = json.loads((out / "reference" / "reports" / "pn2d_bv_comparison.json").read_text())
            self.assertEqual(iv_report["iv"]["candidate_column"], "current_total_A_per_um")
            self.assertTrue(iv_report["iv"]["trend_match"])
            self.assertLess(iv_report["iv"]["orders_of_magnitude"], 1.0)
            self.assertEqual(bv_report["iv"]["candidate_column"], "current_total_A_per_um")
            self.assertLess(bv_report["iv"]["orders_of_magnitude"], 0.15)

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
        return self._root_or_skip(SAMPLE_ENV)

    def _pn2d_root_or_skip(self) -> Path:
        value = os.environ.get(PN2D_ENV)
        root = Path(value) if value else REPO / "reference_tcad" / "pn2d"
        if not root.is_dir():
            self.skipTest(f"{PN2D_ENV} is not set and bundled pn2d sample is missing")
        return root

    def _root_or_skip(self, env_name: str) -> Path:
        value = os.environ.get(env_name)
        if not value:
            self.skipTest(f"{env_name} is not set")
        root = Path(value)
        self.assertTrue(root.is_dir(), f"{env_name} is not a directory: {root}")
        return root

    def _read_curve(self, path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))

    def _assert_doping_csv_has_donors_and_acceptors(self, path: Path) -> None:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreater(len(rows), 0)
        self.assertTrue(any(float(row["donors_cm3"]) > 0.0 for row in rows))
        self.assertTrue(any(float(row["acceptors_cm3"]) > 0.0 for row in rows))

    def _assert_curve_has_finite_currents(self, path: Path) -> None:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            value = row.get("current_total", "0")
            number = float(value)
            self.assertTrue(number == number and abs(number) != float("inf"), f"{path}: {value}")

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
