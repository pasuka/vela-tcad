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
PN2D_2018_ENV = "VELA_SENTAURUS_PN2D_2018_DIR"
BUILD_ENV = "VELA_BUILD_DIR"


class SentaurusSampleIntegrationTest(unittest.TestCase):
    def test_pn2d_sentaurus2018_fixture_import_and_tdx_parity(self) -> None:
        source_root, config_path = self._pn2d_2018_source_and_config_or_skip()
        build_root = Path(os.environ.get(BUILD_ENV, REPO / "build"))
        importer = build_root / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")
        if not importer.is_file():
            self.skipTest(f"Sentaurus HDF5 importer is not built: {importer}")

        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_pn2d_2018_") as tmp:
            out = Path(tmp) / "reference"
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "reference",
                    "--config",
                    str(config_path),
                    "--source-dir",
                    str(source_root),
                    "--output-dir",
                    str(out),
                    "--tdr-importer",
                    str(importer),
                    "--skip-vela-run",
                ],
                check=True,
                cwd=REPO,
            )
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_sentaurus_tdr_tdx.py"),
                    "--geometry-only",
                    "--coordinate-epsilon",
                    "1e-5",
                    "--tdr-export",
                    str(out),
                    "--tdx-dir",
                    str(source_root),
                    "--output-dir",
                    str(out / "reports"),
                ],
                check=True,
                cwd=REPO,
            )

            self.assertTrue((out / "nodes.csv").is_file())
            self.assertTrue((out / "elements.csv").is_file())
            self.assertTrue((out / "contacts.csv").is_file())

            tdx_report = json.loads((out / "reports" / "tdr_tdx_comparison.json").read_text())
            self.assertEqual(tdx_report["status"], "pass")
            self.assertEqual(sorted(tdx_report), ["boundaries", "elements", "nodes", "status"])
            self.assertEqual(tdx_report["nodes"]["tdr_count"], 1943)
            self.assertTrue(tdx_report["nodes"]["count_match"])
            self.assertEqual(tdx_report["elements"]["tdr_count"], 3680)
            self.assertTrue(tdx_report["elements"]["count_match"])
            self.assertEqual(tdx_report["boundaries"]["contacts"]["Anode"]["grd_element_count"], 16)
            self.assertEqual(tdx_report["boundaries"]["contacts"]["Cathode"]["grd_element_count"], 16)
            self.assertTrue(tdx_report["boundaries"]["contacts"]["Anode"]["node_set_match"])
            self.assertTrue(tdx_report["boundaries"]["contacts"]["Cathode"]["node_set_match"])

    def test_pn2d_sentaurus2018_zero_bias_state_compare_when_enabled(self) -> None:
        source_root, config_path = self._pn2d_2018_source_and_config_or_skip()
        build_root = Path(os.environ.get(BUILD_ENV, REPO / "build"))
        importer = build_root / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")
        if not importer.is_file():
            self.skipTest(f"Sentaurus HDF5 importer is not built: {importer}")
        runner_name = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
        runner_candidates = [
            build_root / runner_name,
            build_root / "src" / "tools" / runner_name,
        ]
        runner = next((candidate for candidate in runner_candidates if candidate.is_file()), None)
        if runner is None:
            self.skipTest(
                "Vela runner is not built: expected one of "
                + ", ".join(str(candidate) for candidate in runner_candidates)
            )

        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_pn2d_2018_0v_") as tmp:
            out = Path(tmp) / "reference"
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "sentaurus_import.py"),
                    "reference",
                    "--config",
                    str(config_path),
                    "--source-dir",
                    str(source_root),
                    "--output-dir",
                    str(out),
                    "--tdr-importer",
                    str(importer),
                    "--skip-vela-run",
                ],
                check=True,
                cwd=REPO,
            )
            balance = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_0v_current_balance.py"),
                    "--reference-root",
                    str(out),
                    "--runner",
                    str(runner),
                    "--output-dir",
                    str(out / "reports" / "0v_current_balance"),
                ],
                cwd=REPO,
            )

            balance_report = json.loads(
                (out / "reports" / "0v_current_balance" / "pn2d_0v_current_balance.json").read_text()
            )
            if balance_report["status"] == "error":
                self.assertNotEqual(balance.returncode, 0)
                self.assertEqual(balance_report["classification"], "input_error")
                self.assertTrue(balance_report["classification_reasons"])
                self.assertIn("converged", balance_report["classification_reasons"][0])
                manifest = json.loads((out / "sim_fields" / "0v" / "field_manifest.json").read_text())
                field_mappings = {
                    item["name"]: item["global_node_mapping"]
                    for item in manifest["fields"]
                    if item["name"] in {"ElectricField", "DopingConcentration"}
                }
                self.assertEqual(field_mappings["ElectricField"], "global_vertex_order")
                self.assertEqual(field_mappings["DopingConcentration"], "global_vertex_order")
                return

            self.assertEqual(balance.returncode, 0)
            compare = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_0v_state.py"),
                    "--reference-root",
                    str(out),
                    "--runner",
                    str(runner),
                    "--output-dir",
                    str(out / "reports" / "0v_state"),
                ],
                cwd=REPO,
            )

            self.assertEqual(balance_report["classification"], "balanced")
            self.assertTrue(balance_report["classification_reasons"])
            self.assertIn("Anode", balance_report["terminal_balance"]["contacts"])
            self.assertIn("Cathode", balance_report["terminal_balance"]["contacts"])
            conservation = balance_report["conservation_summary"]
            self.assertTrue(
                conservation["absolute_floor_gate_pass"] or conservation["relative_gate_pass"]
            )
            self.assertLessEqual(
                conservation["abs_pair_sum_A_per_um"],
                conservation["absolute_floor_gate_A_per_um"],
            )
            self.assertIn("top_edges", balance_report["contact_edges"])
            self.assertEqual(
                balance_report["mesh_reference"]["sentaurus_contact_boundary_elements"]["Anode"],
                16,
            )
            self.assertEqual(
                balance_report["mesh_reference"]["sentaurus_contact_boundary_elements"]["Cathode"],
                16,
            )
            self.assertEqual(balance_report["mesh_reference"]["sentaurus_bulk_triangle_elements"], 3680)
            self.assertEqual(
                balance_report["mesh_reference"]["sentaurus_total_elements_including_contacts"],
                3712,
            )

            report = json.loads((out / "reports" / "0v_state" / "pn2d_0v_state_comparison.json").read_text())
            self.assertEqual(compare.returncode, 0)
            self.assertEqual(report["status"], "pass")
            self.assertEqual(
                report["current_balance_report"],
                str(out / "reports" / "0v_current_balance" / "pn2d_0v_current_balance.json"),
            )
            self.assertEqual(report["terminal_currents"]["status"], "pass")
            self.assertNotIn(
                "terminal_currents_not_equal_and_opposite",
                report["terminal_currents"]["failure_reasons"],
            )
            self.assertLessEqual(report["terminal_currents"]["max_abs_A_per_um"], 1.0e-12)
            self.assertLessEqual(
                abs(report["terminal_currents"]["sum_A_per_um"]),
                report["terminal_currents"]["pair_abs_gate_A_per_um"],
            )
            self.assertGreater(report["field_stats"]["ElectrostaticPotential"]["points_compared"], 1000)
            self.assertIn("ni", report["diagnostic_matrix"]["priorities"])
            self.assertIn("OldSlotboom/BGN", report["diagnostic_matrix"]["priorities"])

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

    def _pn2d_2018_source_and_config_or_skip(self) -> tuple[Path, Path]:
        value = os.environ.get(PN2D_2018_ENV)
        case_root = Path(value) if value else REPO / "reference_tcad" / "pn2d_sentaurus2018"
        if (case_root / "source").is_dir():
            source_root = case_root / "source"
            config_path = case_root / "pn2d_sentaurus2018_reference.json"
        else:
            source_root = case_root
            config_path = case_root / "pn2d_sentaurus2018_reference.json"
            if not config_path.is_file():
                config_path = case_root.parent / "pn2d_sentaurus2018_reference.json"
        if not source_root.is_dir():
            self.skipTest(f"{PN2D_2018_ENV} is not set and bundled pn2d Sentaurus 2018 source is missing")
        self.assertTrue(config_path.is_file(), f"missing pn2d Sentaurus 2018 config: {config_path}")
        return source_root, config_path

    def _root_or_skip(self, env_name: str) -> Path:
        value = os.environ.get(env_name)
        if not value:
            self.skipTest(f"{env_name} is not set")
        root = Path(value)
        self.assertTrue(root.is_dir(), f"{env_name} is not a directory: {root}")
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
