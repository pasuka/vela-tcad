#!/usr/bin/env python3
"""Regression coverage for neutral reference TCAD CSV conversion tools."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


class ReferenceTcadToolsTest(unittest.TestCase):
    def test_reference_tcad_device_configs_exist(self) -> None:
        expected = [
            REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json",
            REPO / "reference_tcad" / "nmos2d" / "nmos2d_reference.json",
            REPO / "reference_tcad" / "pmos2d" / "pmos2d_reference.json",
            REPO / "reference_tcad" / "ldmos2d" / "ldmos2d_reference.json",
            REPO / "reference_tcad" / "igbt2d" / "igbt2d_reference.json",
        ]
        for path in expected:
            self.assertTrue(path.is_file(), f"missing reference config: {path}")

    def test_pn2d_sentaurus2018_zero_bias_uses_stable_newton_cap(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
        config = json.loads(path.read_text())
        zero_bias = next(sim for sim in config["simulations"] if sim["name"] == "0v")

        self.assertEqual(config["vela_solver"]["max_update"], 5.0)
        self.assertEqual(zero_bias["vela_solver"]["max_update"], 2.0)

    def test_pn2d_sentaurus2018_iv_disables_minority_relaxation(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
        config = json.loads(path.read_text())
        iv = next(sim for sim in config["simulations"] if sim["name"] == "iv")

        self.assertIs(
            iv["vela_solver"]["contact_boundary_minority_electron_relaxation"],
            False,
        )

    def test_pn2d_sentaurus2018_iv_uses_models_par_alignment(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
        config = json.loads(path.read_text())
        iv = next(sim for sim in config["simulations"] if sim["name"] == "iv")

        self.assertEqual(iv["vela_materials_file"], "pn2d_sentaurus2018_iv_materials.json")
        self.assertEqual(iv["comparison"]["interpolation"], "log_current")
        self.assertAlmostEqual(iv["vela_solver"]["taun"], 1.0e-5)
        self.assertAlmostEqual(iv["vela_solver"]["taup"], 3.0e-6)

        materials_path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / iv["vela_materials_file"]
        materials = json.loads(materials_path.read_text())
        si = next(item for item in materials["materials"] if item["name"] == "Si")
        self.assertAlmostEqual(si["ni"], 1.4638914958767616e10)

    def test_pn2d_sentaurus2018_iv_cmd_writes_multibias_debug_snapshots(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_iv_sdevice.cmd"
        text = path.read_text()

        self.assertIn("eMobility", text)
        self.assertIn("hMobility", text)
        self.assertIn("MaxStep=0.025", text)
        self.assertRegex(text, r'Goal\s*\{\s*Name="Anode"\s*Voltage=2\.0\s*\}')
        self.assertIn(
            'Plot(FilePrefix="pn2d_iv_multibias" Time=(Range=(0 1) Intervals=40) NoOverWrite)',
            text,
        )

    def test_mobility_scan_resolves_relative_materials_file(self) -> None:
        from scripts.scan_pn2d_iv_mobility_candidates import resolve_input_paths

        with tempfile.TemporaryDirectory(prefix="vela_mobility_scan_paths_") as tmp:
            base = Path(tmp)
            (base / "mesh.json").write_text("{}\n")
            (base / "materials.json").write_text("{}\n")

            resolved = resolve_input_paths(
                {
                    "mesh_file": "mesh.json",
                    "materials_file": "materials.json",
                },
                base,
            )

            self.assertEqual(resolved["mesh_file"], str((base / "mesh.json").resolve()))
            self.assertEqual(
                resolved["materials_file"],
                str((base / "materials.json").resolve()),
            )

    def test_probe_pn2d_0v_qf_drivers_help(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts" / "probe_pn2d_0v_qf_drivers.py"),
                "--help",
            ],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("--reference-root", result.stdout)
        self.assertIn("--runner", result.stdout)
        self.assertIn("--output-dir", result.stdout)

    def test_compare_pn2d_0v_current_related_quantities_help(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts" / "compare_pn2d_0v_current_related_quantities.py"),
                "--help",
            ],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("--reference-root", result.stdout)
        self.assertIn("--current-balance", result.stdout)
        self.assertIn("--edge-csv", result.stdout)

    def test_analyze_pn2d_iv_transport_shape_help(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts" / "analyze_pn2d_iv_transport_shape.py"),
                "--help",
            ],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("--reference-root", result.stdout)
        self.assertIn("--biases", result.stdout)
        self.assertIn("--out-dir", result.stdout)

    def test_analyze_pn2d_iv_transport_shape_writes_contact_edge_transport_proxy(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_iv_transport_proxy_") as tmp:
            root = Path(tmp)
            fields = root / "sim_fields" / "iv" / "fields"
            fixed = root / "reports" / "iv_state" / "fixed"
            vela = root / "vela"
            out_dir = root / "reports" / "transport"
            fields.mkdir(parents=True)
            fixed.mkdir(parents=True)
            vela.mkdir()
            (vela / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                    {"id": 1, "node_ids": [1, 3, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "Si", "material": "Si", "cell_ids": [0, 1]}],
                "contacts": [
                    {"id": 0, "name": "Cathode", "region_id": 0, "node_ids": [3]},
                    {"id": 1, "name": "Anode", "region_id": 0, "node_ids": [0]},
                ],
            }) + "\n")

            for name, values in {
                "eCurrentDensity": [100.0, 100.0, 100.0, 100.0],
                "hCurrentDensity": [40.0, 40.0, 40.0, 40.0],
                "TotalCurrentDensity": [140.0, 140.0, 140.0, 140.0],
                "ElectricField": [1.0, 1.0, 1.0, 1.0],
                "srhRecombination": [0.0, 0.0, 0.0, 0.0],
                "eDensity": [1.0e17, 1.0e17, 1.0e17, 1.0e17],
                "hDensity": [2.0e17, 2.0e17, 2.0e17, 2.0e17],
                "eQuasiFermiPotential": [0.000, 0.010, 0.020, 0.030],
                "hQuasiFermiPotential": [0.000, 0.010, 0.020, 0.030],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            (fixed / "iv_1v_fixed_probe_0001_1V.vtk").write_text(
                "\n".join([
                    "# vtk DataFile Version 3.0",
                    "fixture",
                    "ASCII",
                    "DATASET UNSTRUCTURED_GRID",
                    "POINTS 4 double",
                    "0 0 0",
                    "1e-6 0 0",
                    "0 1e-6 0",
                    "1e-6 1e-6 0",
                    "POINT_DATA 4",
                    "SCALARS Potential float 1",
                    "LOOKUP_TABLE default",
                    "0 0 0 0",
                    "SCALARS ElectronQuasiFermi float 1",
                    "LOOKUP_TABLE default",
                    "0 0.006 0.012 0.018",
                    "SCALARS HoleQuasiFermi float 1",
                    "LOOKUP_TABLE default",
                    "0 0.006 0.012 0.018",
                    "SCALARS Electrons float 1",
                    "LOOKUP_TABLE default",
                    "1e23 1e23 1e23 1e23",
                    "SCALARS Holes float 1",
                    "LOOKUP_TABLE default",
                    "2e23 2e23 2e23 2e23",
                    "SCALARS NetDoping float 1",
                    "LOOKUP_TABLE default",
                    "1e23 1e23 -1e23 -1e23",
                    "",
                ])
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "analyze_pn2d_iv_transport_shape.py"),
                    "--reference-root",
                    str(root),
                    "--biases",
                    "1.0",
                    "--out-dir",
                    str(out_dir),
                ],
                check=True,
                cwd=REPO,
            )

            rows = self._read_csv(out_dir / "contact_edge_transport_proxy_compare_1v.csv")

        cathode_electron_drop = next(
            row for row in rows
            if row["contact"] == "Cathode"
            and row["carrier"] == "electron"
            and row["metric"] == "qf_drop_V"
        )
        self.assertAlmostEqual(float(cathode_electron_drop["vela_to_sentaurus_mean"]), 0.6)
        self.assertEqual(cathode_electron_drop["points_sentaurus"], "2")
        self.assertEqual(cathode_electron_drop["points_vela"], "2")

    def test_analyze_pn2d_iv_transport_shape_writes_multibias_transport_proxy(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_iv_multibias_transport_") as tmp:
            root = Path(tmp)
            fixed = root / "reports" / "iv_state" / "fixed"
            vela = root / "vela"
            export = root / "sentaurus_exports" / "0000" / "fields"
            out_dir = root / "reports" / "transport"
            fixed.mkdir(parents=True)
            vela.mkdir()
            export.mkdir(parents=True)
            (vela / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                    {"id": 1, "node_ids": [1, 3, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "Si", "material": "Si", "cell_ids": [0, 1]}],
                "contacts": [
                    {"id": 0, "name": "Cathode", "region_id": 0, "node_ids": [3]},
                    {"id": 1, "name": "Anode", "region_id": 0, "node_ids": [0]},
                ],
            }) + "\n")
            (fixed / "iv_1v_fixed_probe_0006_0.3V.vtk").write_text(
                "\n".join([
                    "# vtk DataFile Version 3.0",
                    "fixture",
                    "ASCII",
                    "DATASET UNSTRUCTURED_GRID",
                    "POINTS 4 double",
                    "0 0 0",
                    "1e-6 0 0",
                    "0 1e-6 0",
                    "1e-6 1e-6 0",
                    "POINT_DATA 4",
                    "SCALARS Potential float 1",
                    "LOOKUP_TABLE default",
                    "0 0 0 0",
                    "SCALARS ElectronQuasiFermi float 1",
                    "LOOKUP_TABLE default",
                    "0 0.006 0.012 0.018",
                    "SCALARS HoleQuasiFermi float 1",
                    "LOOKUP_TABLE default",
                    "0 0.006 0.012 0.018",
                    "SCALARS Electrons float 1",
                    "LOOKUP_TABLE default",
                    "1e23 1e23 1e23 1e23",
                    "SCALARS Holes float 1",
                    "LOOKUP_TABLE default",
                    "2e23 2e23 2e23 2e23",
                    "SCALARS NetDoping float 1",
                    "LOOKUP_TABLE default",
                    "1e23 1e23 -1e23 -1e23",
                    "",
                ])
            )
            self._write_csv(
                fixed / "iv_1v_fixed_probe_contact_edges.csv",
                ["bias_V", "current_contact", "node0", "node1", "mun", "mup", "current_total_A_per_um"],
                [
                    [0.3, "Cathode", 1, 3, 0.07, 0.03, 1.0e-12],
                    [0.3, "Cathode", 2, 3, 0.07, 0.03, 1.0e-12],
                    [0.3, "Anode", 0, 1, 0.07, 0.03, -1.0e-12],
                    [0.3, "Anode", 0, 2, 0.07, 0.03, -1.0e-12],
                ],
            )
            for name, values in {
                "eCurrentDensity": [100.0, 100.0, 100.0, 100.0],
                "hCurrentDensity": [40.0, 40.0, 40.0, 40.0],
                "eDensity": [1.0e17, 1.0e17, 1.0e17, 1.0e17],
                "hDensity": [2.0e17, 2.0e17, 2.0e17, 2.0e17],
                "eQuasiFermiPotential": [0.000, 0.010, 0.020, 0.030],
                "hQuasiFermiPotential": [0.000, 0.010, 0.020, 0.030],
                "eMobility": [700.0, 700.0, 700.0, 700.0],
                "hMobility": [300.0, 300.0, 300.0, 300.0],
            }.items():
                self._write_csv(
                    export / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            self._write_csv(export / "ContactExternalVoltage_region1.csv", ["node_id", "component0"], [[1, 0.0]])
            self._write_csv(export / "ContactExternalVoltage_region2.csv", ["node_id", "component0"], [[2, 0.3]])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "analyze_pn2d_iv_transport_shape.py"),
                    "--reference-root",
                    str(root),
                    "--biases",
                    "0.3",
                    "--sentaurus-multibias-root",
                    str(root / "sentaurus_exports"),
                    "--out-dir",
                    str(out_dir),
                ],
                check=True,
                cwd=REPO,
            )
            rows = self._read_csv(out_dir / "contact_edge_transport_proxy_compare_multibias.csv")
            decision = json.loads((out_dir / "transport_shape_decision.json").read_text())

        qf_drop = next(
            row for row in rows
            if row["bias_V"] == "0.3"
            and row["contact"] == "Cathode"
            and row["carrier"] == "electron"
            and row["metric"] == "qf_drop_V"
        )
        mobility = next(
            row for row in rows
            if row["bias_V"] == "0.3"
            and row["contact"] == "Cathode"
            and row["carrier"] == "electron"
            and row["metric"] == "mobility_m2_V_s"
        )
        self.assertAlmostEqual(float(qf_drop["vela_to_sentaurus_mean"]), 0.6)
        self.assertAlmostEqual(float(mobility["vela_to_sentaurus_mean"]), 1.0)
        self.assertEqual(decision["recommended_next_branch"], "contact_adjacent_qf_state")

    def test_analyze_pn2d_iv_transport_shape_writes_inferred_ni_eff_report(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_iv_inferred_ni_") as tmp:
            root = Path(tmp)
            vela = root / "vela"
            export = root / "sentaurus_exports" / "0000" / "fields"
            out_dir = root / "reports" / "transport"
            vela.mkdir(parents=True)
            export.mkdir(parents=True)
            (vela / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                    {"id": 1, "node_ids": [1, 3, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "Si", "material": "Si", "cell_ids": [0, 1]}],
                "contacts": [
                    {"id": 0, "name": "Cathode", "region_id": 0, "node_ids": [3]},
                    {"id": 1, "name": "Anode", "region_id": 0, "node_ids": [0]},
                ],
            }) + "\n")
            ni_eff = 1.65e10
            vt = 8.617333262145e-5 * 300.0
            psi = [0.20, 0.25, 0.22, 0.27]
            phin = [0.10, 0.14, 0.11, 0.15]
            phip = [0.30, 0.37, 0.32, 0.39]
            electron_density = [
                ni_eff * math.exp((psi_i - phin_i) / vt)
                for psi_i, phin_i in zip(psi, phin)
            ]
            hole_density = [
                ni_eff * math.exp((phip_i - psi_i) / vt)
                for psi_i, phip_i in zip(psi, phip)
            ]
            for name, values in {
                "ElectrostaticPotential": psi,
                "eQuasiFermiPotential": phin,
                "hQuasiFermiPotential": phip,
                "eDensity": electron_density,
                "hDensity": hole_density,
            }.items():
                self._write_csv(
                    export / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            self._write_csv(export / "ContactExternalVoltage_region1.csv", ["node_id", "component0"], [[1, 0.0]])
            self._write_csv(export / "ContactExternalVoltage_region2.csv", ["node_id", "component0"], [[2, 0.25]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "analyze_pn2d_iv_transport_shape.py"),
                    "--reference-root",
                    str(root),
                    "--biases",
                    "0.25",
                    "--sentaurus-multibias-root",
                    str(root / "sentaurus_exports"),
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "sentaurus_inferred_ni_eff_multibias.csv")
            contact_rows = self._read_csv(
                out_dir / "sentaurus_contact_edge_inferred_ni_eff_multibias.csv"
            )

        electron = next(row for row in rows if row["source"] == "electron_inferred")
        hole = next(row for row in rows if row["source"] == "hole_inferred")
        self.assertAlmostEqual(float(electron["median_ni_eff_cm3"]), ni_eff, delta=ni_eff * 1e-12)
        self.assertAlmostEqual(float(hole["median_ni_eff_cm3"]), ni_eff, delta=ni_eff * 1e-12)
        cathode = next(row for row in contact_rows if row["contact"] == "Cathode")
        anode = next(row for row in contact_rows if row["contact"] == "Anode")
        self.assertAlmostEqual(float(cathode["median_ni_eff_cm3"]), ni_eff, delta=ni_eff * 1e-12)
        self.assertAlmostEqual(float(anode["median_ni_eff_cm3"]), ni_eff, delta=ni_eff * 1e-12)

    def test_runner_writes_newton_residual_probe_for_external_state(self) -> None:
        runner = REPO / "build" / ("vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner")
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_newton_residual_probe_") as tmp:
            root = Path(tmp)
            fields = root / "fields"
            fields.mkdir()
            self._write_csv(root / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 2.0e17, 0.0],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 0.0, 2.0e17],
            ])
            (root / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                    {"id": 1, "node_ids": [1, 3, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "Si", "material": "Si", "cell_ids": [0, 1]}],
                "contacts": [
                    {"id": 0, "name": "Cathode", "region_id": 0, "node_ids": [1, 3]},
                    {"id": 1, "name": "Anode", "region_id": 0, "node_ids": [0, 2]},
                ],
            }) + "\n")
            for name, values in {
                "ElectrostaticPotential": [-0.1, 0.1, -0.08, 0.08],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [0.0, 0.0, 0.0, 0.0],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            output = root / "residual.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_residual_probe",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "output_csv": str(output),
                "state_fields_dir": "fields",
                "scaling": {"mode": "unit_scaling"},
                "doping": [{"region": "Si", "donors": 1.0e17, "acceptors": 1.0e17}],
                "contacts": [
                    {"name": "Cathode", "bias": 0.0},
                    {"name": "Anode", "bias": 0.0},
                ],
                "solver": {
                    "method": "gummel_newton",
                    "bandgap_narrowing": "slotboom",
                    "mobility": {"model": "caughey_thomas_field"},
                },
            }, indent=2) + "\n")

            result = subprocess.run(
                [str(runner), "--config", str(config)],
                cwd=REPO,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(output)
            status = json.loads(result.stdout)

        self.assertEqual(len(rows), 4)
        self.assertIn("block_residuals", status)
        self.assertIn("psi_residual", rows[0])
        self.assertIn("phin_residual", rows[0])
        self.assertIn("phip_residual", rows[0])
        self.assertIn("ni_eff_m3", rows[0])
        self.assertAlmostEqual(float(rows[0]["donors_m3"]), 2.0e23)
        self.assertAlmostEqual(float(rows[3]["acceptors_m3"]), 2.0e23)

    def test_prepare_pn2d_sentaurus_multibias_debug_deck(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_multibias_deck_") as tmp:
            root = Path(tmp)
            input_cmd = root / "pn2d_iv_sdevice.cmd"
            output_cmd = root / "pn2d_iv_multibias_debug_sdevice.cmd"
            input_cmd.write_text(
                """File {
  Grid    = "pn2d_msh.tdr"
  Plot    = "pn2d_iv_des.tdr"
  Current = "pn2d_iv.plt"
  Output  = "pn2d_iv.log"
}

Plot {
  Potential
  eDensity
  hDensity
}

Solve {
  Coupled(Iterations=100) { Poisson }

  Quasistationary(
    Goal {
      Name="Anode"
      Voltage=1.0
    }
  ) {
    Coupled {
      Poisson Electron Hole
    }
  }
}
"""
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "prepare_pn2d_sentaurus_multibias_debug.py"),
                    "--input",
                    str(input_cmd),
                    "--output",
                    str(output_cmd),
                    "--file-prefix",
                    "pn2d_iv_probe",
                    "--time-points",
                    "0,0.3,0.8,1.0",
                ],
                cwd=REPO,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            generated = output_cmd.read_text()

        self.assertIn('Plot(FilePrefix="pn2d_iv_probe" Time=(0;0.3;0.8;1.0) NoOverWrite)', generated)
        self.assertIn("eMobility", generated)
        self.assertIn("hMobility", generated)
        self.assertIn('Current = "pn2d_iv.plt"', generated)

    def test_current_related_terminal_rows_recompute_ratio(self) -> None:
        from scripts.compare_pn2d_0v_current_related_quantities import terminal_rows

        with tempfile.TemporaryDirectory(prefix="vela_current_related_") as tmp:
            terminal_csv = Path(tmp) / "terminal.csv"
            self._write_csv(
                terminal_csv,
                [
                    "contact",
                    "current_electron_A_per_um",
                    "current_hole_A_per_um",
                    "current_total_A_per_um",
                    "current_total",
                ],
                [
                    ["Anode", 0.0, 2.0, -2.0, -2.0e6],
                    ["Cathode", 2.0, 0.0, 2.0, 2.0e6],
                ],
            )
            balance = {
                "sentaurus_current_reference": {
                    "final_coupled": {
                        "contacts": {
                            "Anode": {
                                "electron_current": 0.0,
                                "hole_current": -1.0,
                                "total_current": -1.0,
                            },
                            "Cathode": {
                                "electron_current": 0.0,
                                "hole_current": 1.0,
                                "total_current": 1.0,
                            },
                        }
                    }
                },
                "sentaurus_current_parity": {
                    "by_contact": {
                        "Anode": {"abs_ratio_sentaurus_to_vela": 42.0, "sign_relation": "stale"},
                        "Cathode": {"abs_ratio_sentaurus_to_vela": 42.0, "sign_relation": "stale"},
                    }
                },
            }

            rows = terminal_rows(balance, terminal_csv)

        self.assertEqual(rows[0]["contact"], "Anode")
        self.assertEqual(rows[0]["abs_ratio_sentaurus_to_vela"], 0.5)
        self.assertEqual(rows[0]["sign_relation"], "same_sign")
        self.assertEqual(rows[1]["abs_ratio_sentaurus_to_vela"], 0.5)
        self.assertEqual(rows[1]["sign_relation"], "same_sign")

    def test_current_related_terminal_rows_mark_sub_floor_components_as_numerical_zero(self) -> None:
        from scripts.compare_pn2d_0v_current_related_quantities import terminal_rows

        with tempfile.TemporaryDirectory(prefix="vela_current_related_zero_floor_") as tmp:
            terminal_csv = Path(tmp) / "terminal.csv"
            self._write_csv(
                terminal_csv,
                [
                    "contact",
                    "current_electron_A_per_um",
                    "current_hole_A_per_um",
                    "current_total_A_per_um",
                    "current_total",
                ],
                [
                    ["Anode", -1.6e-27, 0.0, -1.6e-27, -1.6e-21],
                    ["Cathode", 0.0, -3.0e-28, 3.0e-28, 3.0e-22],
                ],
            )
            balance = {
                "sentaurus_current_reference": {
                    "final_coupled": {
                        "contacts": {
                            "Anode": {
                                "electron_current": 8.0e-26,
                                "hole_current": -8.0e-25,
                                "total_current": -7.2e-25,
                            },
                            "Cathode": {
                                "electron_current": -8.0e-26,
                                "hole_current": 8.0e-25,
                                "total_current": 7.2e-25,
                            },
                        }
                    }
                },
            }

            rows = terminal_rows(balance, terminal_csv)

        by_contact = {row["contact"]: row for row in rows}
        self.assertEqual(by_contact["Anode"]["vela_electron_pair_A_per_um"], 0.0)
        self.assertEqual(by_contact["Cathode"]["vela_electron_pair_A_per_um"], 0.0)
        self.assertEqual(by_contact["Anode"]["vela_hole_pair_A_per_um"], 0.0)
        self.assertEqual(by_contact["Cathode"]["vela_hole_pair_A_per_um"], 0.0)
        self.assertEqual(by_contact["Anode"]["vela_electron_pair_status"], "pass_absolute_floor")
        self.assertEqual(by_contact["Cathode"]["vela_hole_pair_status"], "pass_absolute_floor")

    def test_pn_export_converts_to_unit_scaling_deck(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_tcad_") as tmp:
            root = Path(tmp)
            export_dir = root / "export"
            out_dir = root / "vela"
            export_dir.mkdir()

            self._write_csv(
                export_dir / "nodes.csv",
                ["id", "x_um", "y_um"],
                [
                    [0, 0.0, 0.0],
                    [1, 1.0, 0.0],
                    [2, 1.0, 1.0],
                    [3, 0.0, 1.0],
                ],
            )
            self._write_csv(
                export_dir / "elements.csv",
                ["id", "node0", "node1", "node2", "region", "material"],
                [
                    [0, 0, 1, 2, "n_region", "Si"],
                    [1, 0, 2, 3, "p_region", "Si"],
                ],
            )
            self._write_csv(
                export_dir / "contacts.csv",
                ["name", "node_ids", "region"],
                [
                    ["anode", "0;3", "p_region"],
                    ["cathode", "1;2", "n_region"],
                ],
            )
            self._write_csv(
                export_dir / "doping.csv",
                ["node_id", "donors_cm3", "acceptors_cm3"],
                [
                    [0, 0.0, 1.0e17],
                    [1, 1.0e17, 0.0],
                    [2, 1.0e17, 0.0],
                    [3, 0.0, 1.0e17],
                ],
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "convert_tcad_export.py"),
                    "--input-dir",
                    str(export_dir),
                    "--output-dir",
                    str(out_dir),
                    "--device",
                    "pn_diode",
                    "--simulation-types",
                    "iv,cv,bv",
                ],
                check=True,
                cwd=REPO,
            )

            mesh = json.loads((out_dir / "mesh.json").read_text())
            iv = json.loads((out_dir / "simulation_iv.json").read_text())
            cv = json.loads((out_dir / "simulation_cv.json").read_text())
            bv = json.loads((out_dir / "simulation_bv.json").read_text())

            self.assertEqual(mesh["nodes"][1]["x"], 1.0)
            self.assertEqual(mesh["triangles"][0]["node_ids"], [0, 1, 2])
            self.assertEqual(mesh["contacts"][0]["node_ids"], [0, 3])
            self.assertEqual(iv["scaling"], {"mode": "unit_scaling"})
            self.assertEqual(iv["mesh_file"], "mesh.json")
            self.assertEqual(iv["node_doping_file"], "doping.csv")
            self.assertEqual(bv["node_doping_file"], "doping.csv")
            self.assertEqual(iv["doping"][0]["donors"], 1.0e17)
            self.assertEqual(iv["doping"][1]["acceptors"], 1.0e17)
            self.assertEqual(iv["doping"][1]["region"], "p_region")
            self.assertEqual(cv["sweep"]["mode"], "cv_quasistatic")
            self.assertEqual(bv["sweep"]["mode"], "bv_reverse")
            self.assertEqual(
                self._read_csv(out_dir / "doping.csv"),
                self._read_csv(export_dir / "doping.csv"),
            )

    def test_compare_reference_curves_writes_json_and_markdown(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_compare_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            report_json = root / "report.json"
            report_md = root / "report.md"

            header = [
                "bias_V",
                "current_total",
                "capacitance_F_per_m",
                "max_electric_field_V_per_cm",
            ]
            self._write_csv(reference, header, [
                [0.0, 1.0e-12, 2.0e-10, 1.0e4],
                [0.5, 1.0e-9, 2.5e-10, 2.0e4],
                [1.0, 1.0e-6, 3.0e-10, 3.0e4],
            ])
            self._write_csv(candidate, header, [
                [0.0, 1.1e-12, 2.1e-10, 1.1e4],
                [0.5, 1.2e-9, 2.6e-10, 2.1e4],
                [1.0, 1.1e-6, 3.1e-10, 3.1e4],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(candidate),
                    "--output-json",
                    str(report_json),
                    "--output-md",
                    str(report_md),
                ],
                check=True,
                cwd=REPO,
            )

            report = json.loads(report_json.read_text())
            self.assertTrue(report["iv"]["trend_match"])
            self.assertTrue(report["cv"]["trend_match"])
            self.assertTrue(report["bv"]["trend_match"])
            self.assertIn("orders_of_magnitude", report["iv"])
            self.assertIn("Reference TCAD Curve Comparison", report_md.read_text())

    def test_compare_reference_curves_enforces_single_curve_thresholds(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_compare_gate_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            bad_candidate = root / "bad_candidate.csv"
            report_json = root / "report.json"
            report_md = root / "report.md"

            header = ["bias_V", "current_total"]
            self._write_csv(reference, header, [
                [0.0, 1.0e-12],
                [0.5, 1.0e-9],
                [1.0, 1.0e-6],
            ])
            self._write_csv(candidate, header, [
                [0.0, 1.0e-12],
                [0.5, 1.1e-9],
                [1.0, 1.05e-6],
            ])
            self._write_csv(bad_candidate, header, [
                [0.0, 1.0e-6],
                [0.5, 1.0e-9],
                [1.0, 1.0e-12],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(candidate),
                    "--output-json",
                    str(root / "default_report.json"),
                    "--output-md",
                    str(root / "default_report.md"),
                ],
                check=True,
                cwd=REPO,
            )
            default_report = json.loads((root / "default_report.json").read_text())
            self.assertEqual(default_report["status"], "pass")
            self.assertFalse(default_report["cv"]["available"])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(candidate),
                    "--output-json",
                    str(report_json),
                    "--output-md",
                    str(report_md),
                    "--kind",
                    "iv",
                    "--require-trend-match",
                    "--min-points",
                    "3",
                    "--max-relative-error",
                    "0.2",
                    "--max-orders-of-magnitude",
                    "0.1",
                ],
                check=True,
                cwd=REPO,
            )
            report = json.loads(report_json.read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["iv"]["points_compared"], 3)
            self.assertEqual(report["checked_kinds"], ["iv"])

            failed = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_reference_curves.py"),
                    "--reference",
                    str(reference),
                    "--candidate",
                    str(bad_candidate),
                    "--output-json",
                    str(root / "bad_report.json"),
                    "--output-md",
                    str(root / "bad_report.md"),
                    "--kind",
                    "iv",
                    "--require-trend-match",
                    "--min-points",
                    "3",
                    "--max-orders-of-magnitude",
                    "0.1",
                ],
                cwd=REPO,
            )
            self.assertNotEqual(failed.returncode, 0)
            bad_report = json.loads((root / "bad_report.json").read_text())
            self.assertEqual(bad_report["status"], "fail")
            self.assertTrue(any("trend" in failure for failure in bad_report["failures"]))

    def test_compare_reference_curves_interpolates_by_bias(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_bias_match_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            out_json = root / "report.json"
            out_md = root / "report.md"
            self._write_csv(reference, ["bias_V", "current_total"], [
                [0.0, 1.0e-12],
                [0.5, 1.0e-9],
                [1.0, 1.0e-6],
            ])
            self._write_csv(candidate, ["bias_V", "current_total"], [
                [0.0, -1.0e-12],
                [0.25, -3.0e-11],
                [0.75, -2.0e-9],
                [1.0, -1.0e-6],
            ])

            subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "compare_reference_curves.py"),
                "--reference", str(reference),
                "--candidate", str(candidate),
                "--output-json", str(out_json),
                "--output-md", str(out_md),
                "--kind", "iv",
                "--candidate-scale", "-1.0",
                "--bias-min", "0.5",
                "--bias-max", "1.0",
                "--max-orders-of-magnitude", "0.25",
                "--require-trend-match",
            ], check=True, cwd=REPO)

            report = json.loads(out_json.read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["iv"]["points_compared"], 2)
            self.assertEqual(report["iv"]["reference_bias_range"], [0.5, 1.0])
            self.assertEqual(report["iv"]["candidate_scale"], -1.0)

    def test_compare_reference_curves_log_current_interpolation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_log_interp_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            out_json = root / "report.json"
            out_md = root / "report.md"
            self._write_csv(reference, ["bias_V", "current_total"], [
                [0.0, 1.0e-12],
                [0.5, 1.0e-9],
                [1.0, 1.0e-6],
            ])
            self._write_csv(candidate, ["bias_V", "current_total"], [
                [0.0, -1.0e-12],
                [1.0, -1.0e-6],
            ])

            subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "compare_reference_curves.py"),
                "--reference", str(reference),
                "--candidate", str(candidate),
                "--output-json", str(out_json),
                "--output-md", str(out_md),
                "--kind", "iv",
                "--candidate-scale", "-1.0",
                "--interpolation", "log_current",
                "--max-orders-of-magnitude", "1e-12",
                "--max-relative-error", "1e-12",
                "--require-trend-match",
            ], check=True, cwd=REPO)

            report = json.loads(out_json.read_text())
            self.assertEqual(report["iv"]["interpolation"], "log_current")
            self.assertEqual(report["iv"]["points_compared"], 3)
            self.assertAlmostEqual(report["iv"]["max_relative_error"], 0.0, delta=1.0e-12)

    def test_compare_reference_curves_can_select_candidate_column(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_reference_column_match_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            out_json = root / "report.json"
            out_md = root / "report.md"
            self._write_csv(reference, ["bias_V", "current_total"], [
                [0.0, 0.0],
                [0.5, 1.0e-15],
                [1.0, 2.0e-15],
            ])
            self._write_csv(candidate, ["bias_V", "current_total", "current_total_A_per_um"], [
                [0.0, 0.0, 0.0],
                [0.5, 1.0e-9, 1.0e-15],
                [1.0, 2.0e-9, 2.0e-15],
            ])

            subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "compare_reference_curves.py"),
                "--reference", str(reference),
                "--candidate", str(candidate),
                "--output-json", str(out_json),
                "--output-md", str(out_md),
                "--kind", "iv",
                "--candidate-column", "current_total_A_per_um",
                "--max-orders-of-magnitude", "0.01",
                "--require-trend-match",
            ], check=True, cwd=REPO)

            report = json.loads(out_json.read_text())
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["iv"]["candidate_column"], "current_total_A_per_um")
            self.assertEqual(report["iv"]["reference_trend"], "increasing")
            self.assertEqual(report["iv"]["candidate_trend"], "increasing")

    def test_checked_in_pn_validation_assets_are_complete(self) -> None:
        pn_dir = REPO / "reference_tcad" / "pn_diode"
        vela_dir = pn_dir / "vela"
        reports_dir = pn_dir / "reports"
        validation_doc = REPO / "docs" / "validation" / "pn_diode_unit_scaling_validation.md"

        readme = (pn_dir / "README.md").read_text()
        for required in [
            "2D silicon",
            "p_region",
            "n_region",
            "anode",
            "cathode",
            "abrupt junction",
            "forward IV",
            "reverse quasi-static CV",
            "reverse BV",
            "unit_scaling",
        ]:
            self.assertIn(required, readme)

        for csv_name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
            self.assertTrue((pn_dir / csv_name).is_file(), csv_name)

        mesh = json.loads((vela_dir / "mesh.json").read_text())
        self.assertEqual({region["name"] for region in mesh["regions"]}, {"n_region", "p_region"})
        self.assertEqual({contact["name"] for contact in mesh["contacts"]}, {"anode", "cathode"})

        for deck_name in ["simulation_iv.json", "simulation_cv.json", "simulation_bv.json"]:
            deck = json.loads((vela_dir / deck_name).read_text())
            self.assertEqual(deck["scaling"], {"mode": "unit_scaling"})
            self.assertEqual(deck["mesh_file"], "mesh.json")

        cv_deck = json.loads((vela_dir / "simulation_cv.json").read_text())
        bv_deck = json.loads((vela_dir / "simulation_bv.json").read_text())
        self.assertLess(cv_deck["sweep"]["stop"], 0.0)
        self.assertLess(bv_deck["sweep"]["stop"], 0.0)

        report = json.loads((reports_dir / "pn_diode_comparison.json").read_text())
        for section in ["iv", "cv", "bv"]:
            self.assertTrue(report[section]["available"], section)
            self.assertTrue(report[section]["trend_match"], section)

        doc = validation_doc.read_text()
        for required in [
            "IV monotonic",
            "finite capacitance",
            "max field non-decreasing",
            "trend and order-of-magnitude",
            "no calibration claim",
        ]:
            self.assertIn(required, doc)

    def test_checked_in_reference_configs_cover_device_fixtures(self) -> None:
        expected = {
            "nmos2d": ["idvd", "idvg", "idvg_surface", "cv", "bv"],
            "pmos2d": ["idvd", "idvg", "idvg_surface", "cv", "bv"],
            "ldmos2d": ["iv", "bv", "fieldplate"],
            "igbt2d": ["iv", "high_injection_iv", "charge_cv", "bv", "bv_ii"],
        }
        for device, simulations in expected.items():
            with self.subTest(device=device):
                path = REPO / "reference_tcad" / device / f"{device}_reference.json"
                self.assertTrue(path.is_file(), path)
                config = json.loads(path.read_text())
                self.assertEqual(config["case"], device)
                self.assertEqual(config["schema"], "vela.reference_tcad.checked_in.v1")
                self.assertEqual(
                    [sim["name"] for sim in config["simulations"]],
                    simulations,
                )
                for sim in config["simulations"]:
                    vela = REPO / "reference_tcad" / device / "vela" / sim["deck"]
                    reference = (
                        REPO
                        / "reference_tcad"
                        / device
                        / "reference_curves"
                        / sim["reference_curve"]
                    )
                    report = REPO / "reference_tcad" / device / "reports" / sim["report_json"]
                    self.assertTrue(vela.is_file(), vela)
                    self.assertTrue(reference.is_file(), reference)
                    self.assertTrue(report.is_file(), report)

    def test_checked_in_mos_validation_assets_are_complete(self) -> None:
        cases = [
            {
                "device": "nmos2d",
                "regions": {"p_body", "n_source", "n_drain", "gate_oxide"},
                "idvg_stop_sign": 1.0,
                "idvd_current_sign": 1.0,
                "validation_doc": "nmos2d_unit_scaling_validation.md",
            },
            {
                "device": "pmos2d",
                "regions": {"n_body", "p_source", "p_drain", "gate_oxide"},
                "idvg_stop_sign": -1.0,
                "idvd_current_sign": 1.0,
                "validation_doc": "pmos2d_unit_scaling_validation.md",
            },
        ]

        for case in cases:
            with self.subTest(device=case["device"]):
                device_dir = REPO / "reference_tcad" / case["device"]
                vela_dir = device_dir / "vela"
                reports_dir = device_dir / "reports"

                readme = (device_dir / "README.md").read_text()
                for required in [
                    "mixed Si/SiO2",
                    "metal_gate",
                    "interface charge",
                    "surface mobility",
                    "multi-terminal quasi-static CV",
                    "all-Si MOS baseline",
                    "unit_scaling",
                ]:
                    self.assertIn(required, readme)

                for csv_name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
                    self.assertTrue((device_dir / csv_name).is_file(), csv_name)

                mesh = json.loads((vela_dir / "mesh.json").read_text())
                self.assertEqual({region["name"] for region in mesh["regions"]}, case["regions"])
                self.assertEqual(
                    {contact["name"] for contact in mesh["contacts"]},
                    {"source", "drain", "body", "gate"},
                )
                self.assertIn("SiO2", {region["material"] for region in mesh["regions"]})

                for deck_name in [
                    "simulation_idvd.json",
                    "simulation_idvg.json",
                    "simulation_idvg_surface.json",
                    "simulation_cv.json",
                    "simulation_bv.json",
                ]:
                    deck = json.loads((vela_dir / deck_name).read_text())
                    self.assertEqual(deck["scaling"], {"mode": "unit_scaling"})
                    self.assertEqual(deck["mesh_file"], "mesh.json")

                surface_deck = json.loads((vela_dir / "simulation_idvg_surface.json").read_text())
                self.assertEqual(surface_deck["solver"]["mobility"]["model"], "caughey_thomas_field_surface")
                self.assertEqual(surface_deck["contacts"][2]["type"], "metal_gate")
                self.assertIn("interfaces", surface_deck)

                idvg_rows = self._read_csv(vela_dir / f"{case['device']}_idvg.csv")
                idvg_currents = [abs(float(row["current_total"])) for row in idvg_rows]
                self.assertGreater(idvg_currents[-1], idvg_currents[0])
                self.assertEqual(
                    self._sign(float(idvg_rows[-1]["bias_V"])),
                    case["idvg_stop_sign"],
                )

                idvd_rows = self._read_csv(vela_dir / f"{case['device']}_idvd.csv")
                nonzero_idvd = [float(row["current_total"]) for row in idvd_rows if abs(float(row["bias_V"])) > 0.0]
                self.assertTrue(nonzero_idvd)
                self.assertEqual(self._sign(nonzero_idvd[-1]), case["idvd_current_sign"])

                cv_rows = self._read_csv(vela_dir / f"{case['device']}_cv.csv")
                cv_columns = set(cv_rows[0])
                for column in [
                    "charge_gate_C_per_m",
                    "capacitance_Cgate_gate_F_per_m",
                    "charge_drain_C_per_m",
                    "capacitance_Cgate_drain_F_per_m",
                    "charge_source_C_per_m",
                    "capacitance_Cgate_source_F_per_m",
                    "charge_body_C_per_m",
                    "capacitance_Cgate_body_F_per_m",
                ]:
                    self.assertIn(column, cv_columns)
                    values = [float(row[column]) for row in cv_rows]
                    self.assertTrue(all(math.isfinite(value) for value in values), column)

                bv_rows = self._read_csv(vela_dir / f"{case['device']}_bv.csv")
                fields = [float(row["max_electric_field_V_per_cm"]) for row in bv_rows]
                self.assertTrue(all(b >= a for a, b in zip(fields, fields[1:])))

                for report_name in ["idvd", "idvg", "cv", "bv"]:
                    report = json.loads((reports_dir / f"{case['device']}_{report_name}_comparison.json").read_text())
                    section = "iv" if report_name in {"idvd", "idvg"} else report_name
                    self.assertTrue(report[section]["available"], report_name)
                    self.assertTrue(report[section]["trend_match"], report_name)

                doc = (REPO / "docs" / "validation" / case["validation_doc"]).read_text()
                for required in [
                    "Id-Vg",
                    "Id-Vd",
                    "multi-terminal CV",
                    "BV max field",
                    "trend and order-of-magnitude",
                    "no calibration claim",
                ]:
                    self.assertIn(required, doc)

    def test_checked_in_power_device_validation_assets_are_complete(self) -> None:
        cases = [
            {
                "device": "ldmos2d",
                "regions": {"p_body", "n_source", "n_drift", "n_drain", "gate_oxide"},
                "contacts": {"source", "body", "drain", "gate"},
                "decks": [
                    "simulation_iv.json",
                    "simulation_bv.json",
                    "simulation_bv_fieldplate.json",
                ],
                "reports": ["iv", "bv", "fieldplate"],
                "validation_doc": "ldmos2d_unit_scaling_validation.md",
                "required_doc": [
                    "low-bias DD-IV",
                    "field-plate",
                    "max field",
                    "engineering trend validation",
                    "no calibration claim",
                ],
            },
            {
                "device": "igbt2d",
                "regions": {"p_collector", "n_buffer", "n_drift", "p_base", "n_emitter"},
                "contacts": {"collector", "gate", "emitter"},
                "decks": [
                    "simulation_iv.json",
                    "simulation_high_injection_iv.json",
                    "simulation_charge_cv.json",
                    "simulation_bv.json",
                    "simulation_bv_ii.json",
                ],
                "reports": ["iv", "high_injection_iv", "charge_cv", "bv", "bv_ii"],
                "validation_doc": "igbt2d_unit_scaling_validation.md",
                "required_doc": [
                    "high-injection",
                    "stored charge",
                    "breakdown diagnostic",
                    "engineering trend validation",
                    "no calibration claim",
                ],
            },
        ]

        for case in cases:
            with self.subTest(device=case["device"]):
                device_dir = REPO / "reference_tcad" / case["device"]
                vela_dir = device_dir / "vela"
                reports_dir = device_dir / "reports"
                reference_dir = device_dir / "reference_curves"

                readme = (device_dir / "README.md").read_text()
                for required in [
                    "unit_scaling",
                    "engineering trend validation",
                    "no calibration claim",
                ]:
                    self.assertIn(required, readme)

                for csv_name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
                    self.assertTrue((device_dir / csv_name).is_file(), csv_name)

                mesh = json.loads((vela_dir / "mesh.json").read_text())
                self.assertEqual({region["name"] for region in mesh["regions"]}, case["regions"])
                self.assertEqual({contact["name"] for contact in mesh["contacts"]}, case["contacts"])

                for deck_name in case["decks"]:
                    deck = json.loads((vela_dir / deck_name).read_text())
                    self.assertEqual(deck["scaling"], {"mode": "unit_scaling"})
                    self.assertEqual(deck["mesh_file"], "mesh.json")

                for report_name in case["reports"]:
                    self.assertTrue((reference_dir / f"{case['device']}_{report_name}_reference.csv").is_file())
                    report = json.loads(
                        (reports_dir / f"{case['device']}_{report_name}_comparison.json").read_text()
                    )
                    available = [
                        key for key in ("iv", "cv", "bv")
                        if report[key]["available"] and report[key]["trend_match"]
                    ]
                    self.assertTrue(available, report_name)

                doc = (REPO / "docs" / "validation" / case["validation_doc"]).read_text()
                for required in case["required_doc"]:
                    self.assertIn(required, doc)

    @staticmethod
    def _write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    @staticmethod
    def _read_csv(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))

    @staticmethod
    def _sign(value: float) -> float:
        if value > 0.0:
            return 1.0
        if value < 0.0:
            return -1.0
        return 0.0


if __name__ == "__main__":
    unittest.main()
