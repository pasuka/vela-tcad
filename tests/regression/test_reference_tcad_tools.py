#!/usr/bin/env python3
"""Regression coverage for neutral reference TCAD CSV conversion tools."""

from __future__ import annotations

import csv
import importlib.util
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
    def test_pn2d_bv_curve_loader_prefers_per_um_current(self) -> None:
        module_path = REPO / "scripts" / "compare_pn2d_bv_multibias_fields.py"
        spec = importlib.util.spec_from_file_location("compare_pn2d_bv_multibias_fields", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_bv_curve_loader_") as tmp:
            curve = Path(tmp) / "vela_curve.csv"
            self._write_csv(curve, [
                "bias_V",
                "current_total",
                "current_total_A_per_um",
                "converged",
            ], [[-0.2, -1.0e-12, -1.0e-18, 1]])

            self.assertEqual(module.load_curve_points(curve), [(-0.2, -1.0e-18)])

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

    def test_pn2d_sentaurus2018_preserves_reported_compensated_doping(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
        config = json.loads(path.read_text())

        self.assertEqual(config["tdr_doping"]["compensated_node_policy"], "reported")

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

    def test_pn2d_sentaurus2018_bv_uses_full_sentaurus_physics(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
        config = json.loads(path.read_text())
        bv = next(sim for sim in config["simulations"] if sim["name"] == "bv")
        solver = bv["vela_solver"]

        self.assertAlmostEqual(solver["abstol"], 1.0e-9)
        self.assertEqual(solver["mobility"]["model"], "masetti_field")
        self.assertEqual(
            solver["mobility"]["high_field_driving_force"],
            "quasi_fermi_gradient",
        )
        self.assertIs(
            solver["contact_boundary_minority_electron_relaxation"],
            False,
        )
        self.assertEqual(bv["vela_materials_file"], "pn2d_sentaurus2018_iv_materials.json")
        self.assertEqual(solver["bandgap_narrowing"], "old_slotboom")
        self.assertEqual(solver["recombination"], ["srh"])
        self.assertEqual(solver["impact_ionization"]["model"], "van_overstraeten")
        self.assertEqual(
            solver["impact_ionization"]["driving_force"],
            "quasi_fermi_gradient",
        )
        self.assertEqual(solver["impact_ionization"]["generation"], "current_density")
        self.assertEqual(
            solver["impact_ionization"]["current_approximation"],
            "density_gradient",
        )
        self.assertEqual(bv["vela_stop"], -0.05)
        self.assertEqual(bv["vela_step"], -0.05)
        self.assertNotIn("candidate_bias_scale", bv["comparison"])
        self.assertEqual(bv["runtime_diagnostic"]["step"], 0.05)
        self.assertEqual(bv["runtime_diagnostic"]["bias_points"], [-0.5, -2.0, -5.0, -10.0, -20.0])

    def test_pn2d_bv_jacobian_block_audit_writes_contract_csv(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_jacobian_block_audit.py"
        spec = importlib.util.spec_from_file_location("diagnose_pn2d_bv_jacobian_block_audit", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_jacobian_block_audit_") as tmp:
            output = Path(tmp) / "audit.csv"
            module.write_placeholder_report(output, [-13.2])
            rows = self._read_csv(output)

        self.assertEqual(
            list(rows[0].keys()),
            [
                "bias_V",
                "block",
                "analytic_norm",
                "fd_norm",
                "diff_norm",
                "rel_diff",
            ],
        )
        self.assertEqual(rows[0]["bias_V"], "-13.2")
        self.assertEqual(
            {row["block"] for row in rows},
            {"srh_auger", "sg_avalanche", "transport", "poisson", "dirichlet_or_gauge"},
        )

    def test_pn2d_bv_jacobian_block_audit_uses_probe_for_finite_rows(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_jacobian_block_audit.py"
        spec = importlib.util.spec_from_file_location("diagnose_pn2d_bv_jacobian_block_audit", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_jacobian_block_probe_") as tmp:
            root = Path(tmp)
            probe = root / "fake_probe.py"
            output = root / "audit.csv"
            probe.write_text(
                "\n".join([
                    "import argparse, csv",
                    "parser = argparse.ArgumentParser()",
                    "parser.add_argument('--output', required=True)",
                    "parser.add_argument('--bias', action='append', required=True)",
                    "args = parser.parse_args()",
                    "with open(args.output, 'w', newline='') as handle:",
                    "    writer = csv.writer(handle)",
                    "    writer.writerow(['bias_V','block','analytic_norm','fd_norm','diff_norm','rel_diff'])",
                    "    writer.writerow([args.bias[0], 'sg_avalanche', '2.0', '2.0', '1.0e-8', '5.0e-9'])",
                ]),
                encoding="utf-8",
            )

            module.write_probe_report(output, [-13.2], [sys.executable, str(probe)])
            rows = self._read_csv(output)

        self.assertEqual(rows[0]["block"], "sg_avalanche")
        self.assertTrue(math.isfinite(float(rows[0]["rel_diff"])))
        self.assertLess(float(rows[0]["rel_diff"]), 1.0e-6)

    def test_pn2d_bv_knee_shape_computes_thresholds_and_log_error(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_knee_shape.py"
        spec = importlib.util.spec_from_file_location("diagnose_pn2d_bv_knee_shape", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_bv_knee_shape_") as tmp:
            root = Path(tmp)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            self._write_csv(reference, ["bias_V", "current_total"], [
                [-10.0, -1.0e-16],
                [-11.0, -1.2e-16],
                [-12.0, -2.5e-16],
            ])
            self._write_csv(candidate, ["bias_V", "current_total_A_per_um"], [
                [-10.0, -1.0e-16],
                [-11.0, -2.1e-16],
                [-12.0, -4.5e-16],
            ])

            ref_points = module.load_curve(reference)
            cand_points = module.load_curve(candidate)
            ref_summary = module.summarize_curve(ref_points, -10.0, -12.0)
            cand_summary = module.summarize_curve(cand_points, -10.0, -12.0)

        self.assertEqual(ref_summary.first_growth_over_1p5, -12.0)
        self.assertEqual(cand_summary.first_growth_over_2p0, -11.0)
        self.assertAlmostEqual(
            module.max_abs_log10_error(cand_points, ref_points, -10.0, -12.0),
            math.log10(4.5 / 2.5),
        )

    def test_runner_real_state_jacobian_block_probe_config_contract(self) -> None:
        config = {
            "simulation_type": "newton_jacobian_block_probe",
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "contacts": [
                {"name": "Anode", "bias": -13.2},
                {"name": "Cathode", "bias": 0.0},
            ],
            "solver": {
                "method": "gummel_newton",
                "recombination": ["srh"],
                "impact_ionization": {
                    "model": "van_overstraeten",
                    "driving_force": "quasi_fermi_gradient",
                    "generation": "current_density",
                    "current_approximation": "density_gradient",
                },
            },
            "state_file": "states/bv_state_bias_m13p200000.csv",
            "output_csv": "reports/jacobian_blocks_m13p2.csv",
            "finite_difference_step": 1.0e-7,
        }

        self.assertEqual(config["simulation_type"], "newton_jacobian_block_probe")
        self.assertEqual(config["state_file"], "states/bv_state_bias_m13p200000.csv")
        self.assertEqual(config["finite_difference_step"], 1.0e-7)

    def test_pn2d_bv_field_compare_script_help(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts" / "compare_pn2d_bv_multibias_fields.py"),
                "--help",
            ],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("--sentaurus-root", result.stdout)
        self.assertIn("--vela-vtk-root", result.stdout)
        self.assertIn("--curve-reference", result.stdout)

    def test_pn2d_bv_contact_current_extraction_filters_edges_and_classifies_parity(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_contact_current_") as tmp:
            root = Path(tmp)
            contact_edges = root / "contact_edges.csv"
            vela_curve = root / "vela_iv.csv"
            reference_curve = root / "sentaurus_iv.csv"
            out = root / "out"

            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "node0",
                "node1",
                "psi0",
                "psi1",
                "phin0",
                "phin1",
                "phip0",
                "phip1",
                "n0",
                "n1",
                "p0",
                "p1",
                "mun",
                "mup",
                "electron_continuity_flux",
                "hole_continuity_flux",
                "current_electron",
                "current_hole",
                "current_total_A_per_um",
            ], [
                [-0.5, "Anode", 10, 1, 2, -0.4, -0.4, 0.0, 0.1, 0.0, -0.1, 1.0e9, 2.0e9, 3.0e9, 4.0e9, 0.1, 0.05, 1.0, -2.0, -1.0e-12, 2.0e-12, -3.0e-18],
                [-0.5, "Anode", 11, 2, 3, -0.3, -0.3, 0.1, 0.2, -0.1, -0.2, 2.0e9, 3.0e9, 4.0e9, 5.0e9, 0.1, 0.05, 2.0, -3.0, -2.0e-12, 4.0e-12, -6.0e-18],
                [-0.5, "Cathode", 12, 4, 5, -0.2, -0.2, 0.0, 0.0, 0.0, 0.0, 1.0e9, 1.0e9, 1.0e9, 1.0e9, 0.1, 0.05, 1.0, 1.0, None, None, None],
            ])
            self._write_csv(vela_curve, [
                "bias_V",
                "current_total",
                "current_total_A_per_um",
            ], [[-0.5, -9.0e-12, -9.0e-18]])
            self._write_csv(reference_curve, [
                "bias_V",
                "current_total",
            ], [
                [-0.5, -1.8e-17],
                [-0.5, -9.0e-17],
            ])
            sentaurus = root / "sentaurus" / "sentaurus_-0.5v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [2, 0.0, 0.0],
            ])
            (sentaurus / "field_manifest.json").write_text(json.dumps({
                "fields": [
                    {
                        "name": "ContactCurrentFlux",
                        "region": 2,
                        "region_name": "Anode",
                    },
                ],
            }), encoding="utf-8")
            self._write_csv(fields / "ContactCurrentFlux_region2.csv", ["node_id", "component0"], [
                [2, -1.8e-17],
            ])
            plt = root / "pn2d_bv.plt"
            plt.write_text(
                """
DF-ISE text
Info { datasets = ["time" "Anode OuterVoltage" "Anode TotalCurrent"] }
Data {
0 -0.5 -1.8e-17
}
""".lstrip(),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_contact_current_extraction.py"),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--curve-candidate",
                    str(vela_curve),
                    "--curve-reference",
                    str(reference_curve),
                    "--out-dir",
                    str(out),
                    "--biases",
                    "-0.5",
                    "--contact",
                    "Anode",
                    "--parity-rtol",
                    "0.01",
                    "--sentaurus-root",
                    str(root / "sentaurus"),
                    "--sentaurus-plt",
                    str(plt),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            summary_rows = self._read_csv(out / "edge_summary.csv")
            self.assertEqual(len(summary_rows), 1)
            row = summary_rows[0]
            self.assertEqual(row["classification"], "terminal_extraction_consistent")
            self.assertEqual(row["sentaurus_current_field_status"], "sentaurus_current_field_unavailable")
            self.assertEqual(int(row["edge_rows"]), 2)
            self.assertEqual(int(row["skipped_incomplete_rows"]), 1)
            self.assertAlmostEqual(float(row["edge_current_A_per_um"]), -9.0e-18)
            self.assertAlmostEqual(float(row["vela_iv_current_A_per_um"]), -9.0e-18)
            self.assertAlmostEqual(float(row["log10_vela_over_sentaurus_current"]), math.log10(0.5))
            self.assertAlmostEqual(float(row["sentaurus_plt_current_A"]), -1.8e-17)
            self.assertAlmostEqual(float(row["sentaurus_contact_current_flux_A"]), -1.8e-17)
            self.assertEqual(row["sentaurus_contact_current_flux_region"], "2")
            self.assertEqual(row["sentaurus_terminal_crosscheck_status"], "sentaurus_plt_contact_flux_consistent")
            self.assertEqual(row["current_gap_branch"], "transport_or_field_mismatch")

            filtered_rows = self._read_csv(out / "contact_edges_filtered.csv")
            self.assertEqual([row["edge_id"] for row in filtered_rows], ["10", "11"])

            report = json.loads((out / "summary.json").read_text())
            self.assertEqual(report["biases"], [-0.5])
            self.assertEqual(report["rows"][0]["classification"], "terminal_extraction_consistent")

    def test_pn2d_bv_contact_current_extraction_aligns_sentaurus_nearest_fields(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_contact_current_sentaurus_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            contact_edges = root / "contact_edges.csv"
            vela_curve = root / "vela_iv.csv"
            reference_curve = root / "sentaurus_iv.csv"
            sentaurus = root / "sentaurus" / "sentaurus_-0.5v"
            fields = sentaurus / "fields"
            out = root / "out"
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 1, "x": 0.0, "y": 0.0},
                    {"id": 2, "x": 2.0, "y": 0.0},
                ],
                "triangles": [],
                "contacts": [],
                "regions": [],
            }), encoding="utf-8")
            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "node0",
                "node1",
                "psi0",
                "psi1",
                "phin0",
                "phin1",
                "phip0",
                "phip1",
                "n0",
                "n1",
                "p0",
                "p1",
                "mun",
                "mup",
                "current_electron",
                "current_hole",
                "current_total_A_per_um",
            ], [[
                -0.5,
                "Anode",
                10,
                1,
                2,
                -0.4,
                -0.2,
                0.0,
                0.2,
                -0.1,
                -0.3,
                1.0e15,
                3.0e15,
                2.0e15,
                4.0e15,
                0.12,
                0.04,
                -1.0e-12,
                2.0e-12,
                -3.0e-18,
            ]])
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [42, 1.0, 0.0],
                [43, 10.0, 0.0],
            ])
            for name, value in {
                "TotalCurrentDensity": 7.0,
                "eCurrentDensity": 2.0,
                "hCurrentDensity": 5.0,
                "eDensity": 2.5e9,
                "hDensity": 3.5e9,
                "eMobility": 1200.0,
                "hMobility": 450.0,
                "eQuasiFermiPotential": 0.15,
                "hQuasiFermiPotential": -0.2,
                "ElectrostaticPotential": -0.25,
            }.items():
                rows = [
                    [42, value],
                    [43, value * 2.0],
                ]
                if name == "hQuasiFermiPotential":
                    rows.extend([[1, -0.10], [2, -0.30]])
                if name == "eQuasiFermiPotential":
                    rows.extend([[1, 0.00], [2, 0.10]])
                if name == "ElectrostaticPotential":
                    rows.extend([[1, -0.40], [2, -0.20]])
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], rows)
            self._write_csv(vela_curve, ["bias_V", "current_total_A_per_um"], [[-0.5, -3.0e-18]])
            self._write_csv(reference_curve, ["bias_V", "current_total"], [[-0.5, -6.0e-18]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_contact_current_extraction.py"),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--curve-candidate",
                    str(vela_curve),
                    "--curve-reference",
                    str(reference_curve),
                    "--mesh",
                    str(mesh),
                    "--sentaurus-root",
                    str(root / "sentaurus"),
                    "--out-dir",
                    str(out),
                    "--biases",
                    "-0.5",
                    "--contact",
                    "Anode",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            summary = self._read_csv(out / "edge_summary.csv")[0]
            self.assertEqual(summary["sentaurus_current_field_status"], "sentaurus_current_field_available")
            self.assertEqual(summary["classification"], "terminal_extraction_consistent")

            edge = self._read_csv(out / "contact_edges_filtered.csv")[0]
            self.assertEqual(edge["sentaurus_nearest_node_id"], "42")
            self.assertAlmostEqual(float(edge["sentaurus_nearest_distance_um"]), 0.0)
            self.assertAlmostEqual(float(edge["sentaurus_total_current_density"]), 7.0)
            self.assertAlmostEqual(float(edge["vela_hole_qf_drop_V"]), -0.2)
            self.assertAlmostEqual(float(edge["sentaurus_hole_qf_drop_V"]), -0.2)
            self.assertAlmostEqual(float(edge["sentaurus_over_vela_hole_qf_drop"]), 1.0)
            self.assertAlmostEqual(float(edge["sentaurus_electron_density_cm3"]), 2.5e9)
            self.assertAlmostEqual(float(edge["vela_electron_density_avg_cm3"]), 2.0e9)
            self.assertAlmostEqual(float(edge["sentaurus_electron_mobility_cm2_V_s"]), 1200.0)
            self.assertAlmostEqual(float(edge["vela_electron_mobility_cm2_V_s"]), 1200.0)
            self.assertAlmostEqual(float(edge["sentaurus_electron_qf_V"]), 0.15)
            self.assertAlmostEqual(float(edge["sentaurus_potential_V"]), -0.25)
            self.assertAlmostEqual(float(edge["vela_potential_avg_V"]), -0.3)
            self.assertAlmostEqual(float(edge["vela_electron_qf_avg_V"]), 0.1)
            self.assertAlmostEqual(float(edge["vela_hole_qf_avg_V"]), -0.2)

    def test_pn2d_bv_anode_hole_qf_sensitivity_recomputes_edge_current(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_hole_qf_sensitivity_") as tmp:
            root = Path(tmp)
            contact_edges = root / "contact_edges_filtered.csv"
            summary_csv = root / "edge_summary.csv"
            out = root / "out"

            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "node0",
                "node1",
                "edge_length_m",
                "edge_couple_m",
                "outward_sign",
                "psi0",
                "psi1",
                "phip0",
                "phip1",
                "ni0",
                "ni1",
                "mup",
                "current_electron",
                "current_hole",
                "current_total_A_per_um",
                "sentaurus_hole_qf_drop_V",
            ], [[
                -12.8,
                "Anode",
                7,
                0,
                1,
                1.0e-6,
                2.0e-6,
                -1.0,
                -13.20,
                -13.20,
                -12.80,
                -12.80,
                1.6e16,
                1.6e16,
                0.03,
                -5.0e-13,
                0.0,
                -5.0e-19,
                -2.5e-14,
            ]])
            self._write_csv(summary_csv, [
                "bias_V",
                "sentaurus_current_A",
                "sentaurus_plt_current_A",
            ], [[-12.8, -2.0e-17, -2.0e-17]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_anode_hole_qf_sensitivity.py"),
                    "--contact-edges",
                    str(contact_edges),
                    "--edge-summary",
                    str(summary_csv),
                    "--out-dir",
                    str(out),
                    "--bias",
                    "-12.8",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "hole_qf_sensitivity_edges.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertLess(float(row["baseline_hole_current_A_per_um"]), 1.0e-30)
            self.assertGreater(float(row["sentaurus_drop_hole_current_A_per_um"]), 0.0)
            self.assertLess(
                float(row["sentaurus_drop_total_current_A_per_um"]),
                float(row["baseline_total_current_A_per_um"]),
            )

            summary = json.loads((out / "hole_qf_sensitivity_summary.json").read_text())
            self.assertEqual(summary["bias_V"], -12.8)
            self.assertEqual(summary["edge_count"], 1)
            self.assertGreater(summary["sentaurus_drop_hole_current_A_per_um"], 0.0)
            self.assertGreater(summary["abs_log10_error_baseline"], summary["abs_log10_error_sentaurus_drop"])

    def test_pn2d_bv_contact_boundary_projection_classifies_erased_qf_drop(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_contact_projection_") as tmp:
            root = Path(tmp)
            restart = root / "sentaurus_state_restart.csv"
            final = root / "last_state.csv"
            contact_edges = root / "contact_edges.csv"
            out = root / "out"

            self._write_csv(restart, [
                "node_id",
                "psi",
                "phin",
                "phip",
                "electrons_m3",
                "holes_m3",
            ], [
                [0, -13.0, -12.599, -12.6, 1.0e10, 1.0e23],
                [1, -13.0, -12.6, -12.600000000000009, 1.0e10, 1.0e23],
            ])
            self._write_csv(final, [
                "node_id",
                "psi",
                "phin",
                "phip",
                "electrons_m3",
                "holes_m3",
            ], [
                [0, -13.0, -12.599, -12.6, 1.0e10, 1.0e23],
                [1, -13.0, -12.6, -12.6, 1.0e10, 1.0e23],
            ])
            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "node0",
                "node1",
                "outward_sign",
                "phip0",
                "phip1",
            ], [[
                -12.6,
                "Anode",
                7,
                0,
                1,
                -1.0,
                -12.6,
                -12.6,
            ]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_contact_boundary_projection.py"),
                    "--restart-state",
                    str(restart),
                    "--final-state",
                    str(final),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--out-dir",
                    str(out),
                    "--bias",
                    "-12.6",
                    "--contact",
                    "Anode",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "contact_boundary_projection_edges.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["drop_source"], "contact_projection_erased_restart_drop")
            self.assertEqual(row["contact_node"], "1")
            self.assertEqual(row["interior_node"], "0")
            self.assertGreater(abs(float(row["restart_interior_minus_contact_phip_V"])), 1.0e-15)
            self.assertLess(abs(float(row["projected_interior_minus_contact_phip_V"])), 1.0e-30)
            self.assertLess(abs(float(row["final_interior_minus_contact_phip_V"])), 1.0e-30)

            summary = json.loads((out / "contact_boundary_projection_summary.json").read_text())
            self.assertEqual(summary["edge_count"], 1)
            self.assertEqual(summary["projection_erased_restart_drop_edges"], 1)
            self.assertEqual(summary["dominant_drop_source"], "contact_projection_erased_restart_drop")

    def test_pn2d_bv_contact_qf_floor_reporting_compares_restart_and_ulp_policies(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_contact_qf_floor_") as tmp:
            root = Path(tmp)
            contact_edges = root / "contact_edges_filtered.csv"
            summary_csv = root / "edge_summary.csv"
            restart = root / "sentaurus_state_restart.csv"
            out = root / "out"
            ulp5 = math.ulp(-12.8) * 5.0

            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "node0",
                "node1",
                "edge_length_m",
                "edge_couple_m",
                "outward_sign",
                "psi0",
                "psi1",
                "phip0",
                "phip1",
                "ni0",
                "ni1",
                "mup",
                "current_electron",
                "current_hole",
                "current_total_A_per_um",
                "sentaurus_hole_qf_drop_V",
            ], [[
                -12.8,
                "Anode",
                7,
                0,
                1,
                1.0e-6,
                2.0e-6,
                -1.0,
                -13.20,
                -13.20,
                -12.80,
                -12.80,
                1.6e16,
                1.6e16,
                0.03,
                -5.0e-13,
                0.0,
                -5.0e-19,
                -ulp5,
            ]])
            self._write_csv(summary_csv, [
                "bias_V",
                "sentaurus_current_A",
                "sentaurus_plt_current_A",
            ], [[-12.8, -8.0e-18, -8.0e-18]])
            self._write_csv(restart, [
                "node_id",
                "psi",
                "phin",
                "phip",
                "electrons_m3",
                "holes_m3",
            ], [
                [0, -13.20, -12.8, -12.8, 1.0e10, 1.0e23],
                [1, -13.20, -12.8, -12.8 - ulp5, 1.0e10, 1.0e23],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_contact_qf_floor_reporting.py"),
                    "--contact-edges",
                    str(contact_edges),
                    "--edge-summary",
                    str(summary_csv),
                    "--restart-state",
                    str(restart),
                    "--out-dir",
                    str(out),
                    "--bias",
                    "-12.8",
                    "--ulp-floor-multiplier",
                    "5",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "contact_qf_floor_reporting_edges.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertLess(float(row["restart_hole_qf_drop_V"]), 0.0)
            self.assertAlmostEqual(float(row["restart_hole_qf_drop_V"]), -ulp5)
            self.assertAlmostEqual(float(row["ulp_floor_hole_qf_drop_V"]), -ulp5)
            self.assertGreater(float(row["restart_drop_hole_current_A_per_um"]), 0.0)
            self.assertGreater(float(row["ulp_floor_hole_current_A_per_um"]), 0.0)

            summary = json.loads((out / "contact_qf_floor_reporting_summary.json").read_text())
            self.assertEqual(summary["edge_count"], 1)
            self.assertEqual(summary["restart_drop_edges"], 1)
            self.assertGreater(
                summary["abs_log10_error_baseline"],
                summary["abs_log10_error_restart_drop"],
            )
            self.assertGreater(
                summary["abs_log10_error_baseline"],
                summary["abs_log10_error_ulp_floor"],
            )

    def test_pn2d_bv_qf_floor_restart_point_verifier_checks_enabled_current(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_qf_floor_verify_") as tmp:
            root = Path(tmp)
            iv_csv = root / "iv.csv"
            contact_edges = root / "contact_edges.csv"
            restart_summary = root / "restart_summary.json"
            out_json = root / "verification.json"

            self._write_csv(iv_csv, [
                "bias_V",
                "converged",
                "current_total_A_per_um",
            ], [
                [-12.7, 1, -7.288777764227923e-17],
            ])
            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "current_hole",
                "current_total_A_per_um",
                "hole_qf_drop_override_applied",
            ], [
                [-12.7, "Anode", 1, 3.0e-17, -3.1e-17, 1],
                [-12.7, "Anode", 2, 4.0e-17, -4.2e-17, 1],
                [-12.7, "Cathode", 9, 0.0, 1.0e-20, 0],
            ])
            restart_summary.write_text(json.dumps({
                "bias_V": -12.7,
                "contact": "Anode",
                "reference_current_A": -7.90491101526414e-17,
                "restart_drop_total_current_A_per_um": -7.288777764227923e-17,
                "abs_log10_error_restart_drop": 0.03524227675351084,
            }))

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "verify_pn2d_bv_qf_floor_restart_point.py"),
                    "--iv-csv",
                    str(iv_csv),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--restart-summary",
                    str(restart_summary),
                    "--bias",
                    "-12.7",
                    "--contact",
                    "Anode",
                    "--expected-override-edges",
                    "2",
                    "--out-json",
                    str(out_json),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads(out_json.read_text())
            self.assertTrue(summary["passed"])
            self.assertEqual(summary["override_edge_count"], 2)
            self.assertLess(summary["abs_log10_error_vs_reference"], 0.07)
            self.assertLess(summary["abs_delta_vs_restart_current_A_per_um"], 1.0e-28)

            failed = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "verify_pn2d_bv_qf_floor_restart_point.py"),
                    "--iv-csv",
                    str(iv_csv),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--restart-summary",
                    str(restart_summary),
                    "--bias",
                    "-12.7",
                    "--contact",
                    "Anode",
                    "--expected-override-edges",
                    "3",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("override", failed.stderr)

    def test_pn2d_bv_qf_floor_reporting_stability_compares_default_and_policy(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_qf_floor_stability_") as tmp:
            root = Path(tmp)
            restart_summary = root / "restart_summary.json"
            default_iv = root / "default_iv.csv"
            default_edges = root / "default_edges.csv"
            qf_iv = root / "qf_iv.csv"
            qf_edges = root / "qf_edges.csv"
            out_json = root / "stability.json"

            restart_summary.write_text(json.dumps({
                "bias_V": -13.2,
                "contact": "Anode",
                "reference_current_A": -8.384720888068e-17,
                "baseline_total_current_A_per_um": -5.858870183437891e-17,
                "restart_drop_total_current_A_per_um": -7.324143619438298e-17,
                "abs_log10_error_baseline": 0.155674734983936,
                "abs_log10_error_restart_drop": 0.058731758506640475,
            }))
            self._write_csv(default_iv, [
                "bias_V",
                "converged",
                "current_total_A_per_um",
            ], [[-13.2, 1, -5.858870183437891e-17]])
            self._write_csv(default_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "current_total_A_per_um",
                "hole_qf_drop_override_applied",
            ], [
                [-13.2, "Anode", 1, -2.0e-17, 0],
                [-13.2, "Anode", 2, -3.858870183437891e-17, 0],
            ])
            self._write_csv(qf_iv, [
                "bias_V",
                "converged",
                "current_total_A_per_um",
            ], [[-13.2, 1, -7.324143619438298e-17]])
            self._write_csv(qf_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "current_total_A_per_um",
                "hole_qf_drop_override_applied",
            ], [
                [-13.2, "Anode", 1, -3.0e-17, 1],
                [-13.2, "Anode", 2, -4.324143619438298e-17, 1],
                [-13.2, "Cathode", 9, 1.0e-20, 0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "verify_pn2d_bv_qf_floor_reporting_stability.py"),
                    "--default-iv-csv",
                    str(default_iv),
                    "--default-contact-edge-csv",
                    str(default_edges),
                    "--qf-floor-iv-csv",
                    str(qf_iv),
                    "--qf-floor-contact-edge-csv",
                    str(qf_edges),
                    "--restart-summary",
                    str(restart_summary),
                    "--bias",
                    "-13.2",
                    "--contact",
                    "Anode",
                    "--expected-override-edges",
                    "2",
                    "--out-json",
                    str(out_json),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads(out_json.read_text())
            self.assertTrue(summary["passed"])
            self.assertEqual(summary["default_override_edge_count"], 0)
            self.assertEqual(summary["qf_floor_override_edge_count"], 2)
            self.assertLess(summary["qf_floor_abs_log10_error"], summary["default_abs_log10_error"])
            self.assertLess(summary["default_delta_vs_baseline_current_A_per_um"], 1.0e-28)
            self.assertLess(summary["qf_floor_delta_vs_restart_current_A_per_um"], 1.0e-28)

    def test_pn2d_bv_sentaurus_terminal_current_crosscheck_reports_flux_plt_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_sentaurus_terminal_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus" / "sentaurus_-13p2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            out_json = root / "terminal_crosscheck.json"

            (sentaurus / "field_manifest.json").write_text(json.dumps({
                "fields": [
                    {
                        "name": "ContactCurrentFlux",
                        "region": 2,
                        "region_name": "Anode",
                    },
                ],
            }), encoding="utf-8")
            self._write_csv(fields / "ContactCurrentFlux_region2.csv", ["node_id", "component0"], [
                [2, -7.269471996838386e-17],
            ])
            plt = root / "pn2d_bv.plt"
            plt.write_text(
                """
DF-ISE text
Info { datasets = ["time" "Anode OuterVoltage" "Anode TotalCurrent"] }
Data {
0 -13.2 -8.384720888068e-17
}
""".lstrip(),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "verify_pn2d_sentaurus_terminal_current_crosscheck.py"),
                    "--sentaurus-root",
                    str(root / "sentaurus"),
                    "--sentaurus-plt",
                    str(plt),
                    "--bias",
                    "-13.2",
                    "--contact",
                    "Anode",
                    "--max-relative-difference",
                    "0.01",
                    "--out-json",
                    str(out_json),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads(out_json.read_text())
            self.assertEqual(summary["status"], "sentaurus_plt_contact_flux_mismatch")
            self.assertEqual(summary["contact_current_flux_region"], 2)
            self.assertAlmostEqual(summary["contact_current_flux_A"], -7.269471996838386e-17)
            self.assertAlmostEqual(summary["plt_current_A"], -8.384720888068e-17)
            self.assertGreater(summary["relative_difference"], 0.13)
            self.assertAlmostEqual(
                summary["log10_contact_flux_over_plt"],
                math.log10(7.269471996838386e-17 / 8.384720888068e-17),
            )

    def test_pn2d_bv_debug_direction_summarizes_terminal_and_field_evidence(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_debug_direction_") as tmp:
            root = Path(tmp)
            ranking = root / "debug_ranking.json"
            terminal = root / "terminal.json"
            qf_stability = root / "qf_stability.json"
            out_json = root / "direction.json"

            ranking.write_text(json.dumps({
                "curve": {
                    "max_abs_log10_error": 0.05874,
                },
                "field_rankings": [
                    {
                        "quantity": "electric_field",
                        "field_error": 0.775,
                        "junction_error": 0.108,
                        "p_contact_error": 0.007,
                        "n_contact_error": 0.263,
                    },
                    {
                        "quantity": "avalanche_generation",
                        "avalanche_status": "thresholded_peak",
                        "field_error": 12.85,
                        "junction_error": 0.474,
                        "sentaurus_p99": 3.19365234469398e15,
                        "vela_p99": 1.2178566e15,
                    },
                    {
                        "quantity": "electron_density",
                        "field_error": 0.140,
                        "junction_error": 0.159,
                    },
                ],
            }), encoding="utf-8")
            terminal.write_text(json.dumps({
                "status": "sentaurus_plt_contact_flux_mismatch",
                "relative_difference": 0.1330096619932436,
                "log10_contact_flux_over_plt": -0.061985742401219186,
            }), encoding="utf-8")
            qf_stability.write_text(json.dumps({
                "qf_floor_abs_log10_error": 0.058731758506640475,
                "default_abs_log10_error": 0.155674734983936,
            }), encoding="utf-8")

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "summarize_pn2d_bv_debug_direction.py"),
                    "--debug-ranking",
                    str(ranking),
                    "--terminal-crosscheck",
                    str(terminal),
                    "--qf-floor-stability",
                    str(qf_stability),
                    "--out-json",
                    str(out_json),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads(out_json.read_text())
            self.assertEqual(summary["terminal_assessment"], "plt_gap_explained_by_sentaurus_terminal_convention")
            self.assertEqual(summary["next_actions"][0]["id"], "thresholded_avalanche_support")
            self.assertEqual(summary["next_actions"][1]["id"], "junction_electric_field_reconstruction")
            self.assertIn("electron_density", summary["secondary_quantities"])

    def test_pn2d_bv_thresholded_avalanche_support_reports_mask_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_avalanche_support_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            vtk = root / "vela_0000_-13.2V.vtk"
            out_dir = root / "out"

            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
                [3, 1.0, 1.0],
            ])
            self._write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 10.0],
                [2, 20.0],
                [3, 100.0],
            ])
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 double
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
POINT_DATA 4
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
0 5e6 80e6 30e6
""".lstrip(),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_thresholded_avalanche_support.py"),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--vela-vtk",
                    str(vtk),
                    "--out-dir",
                    str(out_dir),
                    "--percentile",
                    "75",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((out_dir / "thresholded_avalanche_support_summary.json").read_text())
            self.assertEqual(summary["matched_points"], 4)
            self.assertEqual(summary["sentaurus_active_count"], 1)
            self.assertEqual(summary["vela_active_count"], 1)
            self.assertEqual(summary["overlap_count"], 0)
            self.assertEqual(summary["false_positive_count"], 1)
            self.assertEqual(summary["false_negative_count"], 1)
            self.assertAlmostEqual(summary["jaccard"], 0.0)
            self.assertAlmostEqual(summary["peak_separation_um"], 1.0)

    def test_pn2d_bv_support_local_factors_summarizes_source_ratios(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_support_factors_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            support = root / "support.csv"
            vtk = root / "vela_0000_-13.2V.vtk"
            out_dir = root / "out"

            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
            ])
            for name, values in {
                "ImpactIonization": [100.0, 200.0],
                "ElectricField": [4.0, 8.0],
                "eMobility": [1000.0, 500.0],
                "hMobility": [400.0, 200.0],
                "eDensity": [1.0e10, 2.0e10],
                "hDensity": [3.0e10, 6.0e10],
                "eCurrentDensity": [1.0, 4.0],
                "hCurrentDensity": [2.0, 8.0],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [
                [0, 0.0, 0.0, 100.0, 50.0, 1, 1, "overlap"],
                [1, 1.0, 0.0, 200.0, 20.0, 1, 0, "false_negative"],
            ])
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 double
0 0 0
1e-6 0 0
POINT_DATA 2
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
50e6 20e6
SCALARS ElectricField float 1
LOOKUP_TABLE default
2 2
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.05 0.10
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.02 0.04
SCALARS Electrons float 1
LOOKUP_TABLE default
5e15 10e15
SCALARS Holes float 1
LOOKUP_TABLE default
15e15 30e15
SCALARS ElectronHighFieldDrive float 1
LOOKUP_TABLE default
1 3
SCALARS HoleHighFieldDrive float 1
LOOKUP_TABLE default
2 4
""".lstrip(),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_support_local_factors.py"),
                    "--support-csv",
                    str(support),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--vela-vtk",
                    str(vtk),
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "support_local_factors.csv")
            self.assertEqual(len(rows), 2)
            overlap = next(row for row in rows if row["support_class"] == "overlap")
            self.assertAlmostEqual(float(overlap["source_ratio_vela_over_sentaurus"]), 0.5)
            self.assertAlmostEqual(float(overlap["electric_field_ratio_vela_over_sentaurus"]), 0.5)
            self.assertAlmostEqual(float(overlap["electron_mobility_ratio_vela_over_sentaurus"]), 0.5)
            summary = json.loads((out_dir / "support_local_factors_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)
            self.assertAlmostEqual(
                summary["support_classes"]["false_negative"]["source_ratio_median"],
                0.1,
            )

    def test_pn2d_bv_sg_source_ownership_backprojects_edges_to_support_nodes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_sg_ownership_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            support = root / "support.csv"
            sg_edges = root / "sg_edges.csv"
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 3.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 2.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                ],
                "contacts": [],
                "regions": [],
            }), encoding="utf-8")
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [
                [0, 0.0, 0.0, 4.0e6, 3.0e6, 1, 1, "overlap"],
                [1, 3.0, 0.0, 8.0e6, 5.0e6, 1, 0, "false_negative"],
            ])
            self._write_csv(sg_edges, [
                "edge_id",
                "node0",
                "node1",
                "source_integral",
                "electron_source_integral",
                "hole_source_integral",
            ], [
                [10, 0, 1, 6.0, 5.0, 1.0],
                [11, 1, 2, 4.0, 3.0, 1.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_source_ownership.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(sg_edges),
                    "--mesh",
                    str(mesh),
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "sg_source_ownership_nodes.csv")
            self.assertEqual(len(rows), 2)
            node0 = next(row for row in rows if row["node_id"] == "0")
            node1 = next(row for row in rows if row["node_id"] == "1")
            self.assertAlmostEqual(float(node0["reconstructed_node_source_integral"]), 3.0)
            self.assertAlmostEqual(float(node1["reconstructed_node_source_integral"]), 5.0)
            self.assertAlmostEqual(float(node0["vtk_node_source_integral"]), 3.0)
            self.assertAlmostEqual(float(node1["vtk_node_source_integral"]), 5.0)
            summary = json.loads((out_dir / "sg_source_ownership_summary.json").read_text())
            self.assertAlmostEqual(summary["support_classes"]["overlap"]["reconstructed_over_vtk"], 1.0)
            self.assertAlmostEqual(summary["support_classes"]["false_negative"]["reconstructed_over_vtk"], 1.0)

    def test_pn2d_bv_sg_source_ownership_filters_cxx_edge_dump_bias(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_sg_ownership_cxx_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            support = root / "support.csv"
            sg_edges = root / "sg_edges.csv"
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 3.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 2.0},
                ],
                "triangles": [{"id": 0, "node_ids": [0, 1, 2]}],
            }), encoding="utf-8")
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 4.0e6, 3.0e6, 1, 1, "overlap"]])
            self._write_csv(sg_edges, [
                "point_index",
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "edge_source_integral",
                "electron_source_integral",
                "hole_source_integral",
            ], [
                [0, -12.85, 10, 0, 1, 100.0, 80.0, 20.0],
                [1, -13.2, 10, 0, 1, 6.0, 5.0, 1.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_source_ownership.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(sg_edges),
                    "--mesh",
                    str(mesh),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "sg_source_ownership_nodes.csv")
            self.assertAlmostEqual(float(rows[0]["reconstructed_node_source_integral"]), 3.0)

    def test_pn2d_bv_cxx_edge_flux_current_comparator_infers_weighted_alpha(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_cxx_flux_current_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"

            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0e6, 5.0e5, 1, 1, "overlap"]])
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0]])
            for name, value in {
                "ImpactIonization": 1.0e6,
                "eCurrentDensity": 1.602176634e-19 * 3.0e7,
                "hCurrentDensity": 1.602176634e-19 * 2.0e7,
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [[0, value]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
                "electron_source_integral",
                "hole_source_integral",
                "edge_source_integral",
            ], [
                [-13.2, 10, 0, 1, 3.0e11, 2.0e11, 4.0, 2.0, 12.0, 4.0, 16.0],
                [-13.2, 11, 0, 2, 1.0e11, 1.0e11, 6.0, 1.0, 6.0, 1.0, 7.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_cxx_edge_flux_current.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "cxx_edge_flux_current_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["sentaurus_weighted_alpha_m_inv"]), 2.0)
            self.assertAlmostEqual(float(row["cxx_weighted_alpha_m_inv"]), 23.0 / 7.0)
            self.assertAlmostEqual(float(row["cxx_particle_flux_sum"]), 7.0e11)
            summary = json.loads((out_dir / "cxx_edge_flux_current_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_source_geometry_decomposes_area_volume_factor(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_source_geometry_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            mesh = root / "mesh.json"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 3.0e6, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 2.0e6},
                ],
                "triangles": [
                    {"node_ids": [0, 1, 2]},
                ],
            }), encoding="utf-8")
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 6.0e-6, 0.75e-6, 1, 1, "overlap"]])
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0]])
            self._write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [[0, 6.0e-6]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
                "electron_source_integral",
                "hole_source_integral",
                "edge_source_integral",
                "node0_source_integral",
                "node1_source_integral",
            ], [
                [-13.2, 10, 0, 1, 0.5, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.5, 0.5],
                [-13.2, 11, 0, 2, 0.5, 4.0, 0.0, 1.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0],
                [-13.2, 12, 0, 3, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_source_geometry.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--mesh",
                    str(mesh),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "source_geometry_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["node_volume_m2"]), 1.0)
            self.assertAlmostEqual(float(row["cxx_endpoint_area_sum_m2"]), 0.75)
            self.assertAlmostEqual(float(row["cxx_endpoint_area_over_node_volume"]), 0.75)
            self.assertAlmostEqual(float(row["cxx_active_endpoint_area_fraction"]), 2.0 / 3.0)
            self.assertAlmostEqual(float(row["cxx_source_density_proxy_m3_s"]), 2.0)
            self.assertAlmostEqual(float(row["cxx_source_density_proxy_over_sentaurus_generation"]), 1.0 / 3.0)
            self.assertAlmostEqual(float(row["cxx_active_source_density_proxy_m3_s"]), 3.0)
            self.assertAlmostEqual(float(row["cxx_active_source_density_proxy_over_sentaurus_generation"]), 0.5)
            self.assertAlmostEqual(float(row["cxx_source_over_sentaurus_source"]), 0.25)
            summary = json.loads((out_dir / "source_geometry_summary.json").read_text())
            self.assertAlmostEqual(
                summary["support_classes"]["overlap"]["source_ratio_factor_product_median"],
                0.25,
            )

    def test_pn2d_bv_edge_direction_policy_compares_active_reconstructions(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_edge_direction_policy_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            mesh = root / "mesh.json"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 3.0e6, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 2.0e6},
                ],
                "triangles": [
                    {"node_ids": [0, 1, 2]},
                ],
            }), encoding="utf-8")
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 4.0e-6, 2.0e-6, 1, 1, "overlap"]])
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0]])
            self._write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [[0, 4.0e-6]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
                "edge_source_integral",
                "node0_source_integral",
                "node1_source_integral",
            ], [
                [-13.2, 10, 0, 1, 0.0, 0.0, 1.0, 0.0, 0.5, 4.0, 0.0, 1.0, 0.0, 2.0, 1.0, 1.0],
                [-13.2, 11, 0, 1, 0.0, 0.0, 1.0, 0.0, 0.5, 4.0, 0.0, 1.0, 0.0, 2.0, 1.0, 1.0],
                [-13.2, 12, 0, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_edge_direction_source_policy.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--mesh",
                    str(mesh),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "edge_direction_source_policy_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["dominant_active_axis"], "x")
            self.assertAlmostEqual(float(row["junction_normal_active_area_fraction"]), 0.5)
            self.assertAlmostEqual(float(row["junction_tangent_total_area_fraction"]), 0.5)
            self.assertAlmostEqual(float(row["full_face_average_ratio_to_sentaurus"]), 0.5)
            self.assertAlmostEqual(float(row["active_edge_average_ratio_to_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(row["active_edge_density_sum_ratio_to_sentaurus"]), 2.0)
            summary = json.loads((out_dir / "edge_direction_source_policy_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)
            self.assertEqual(summary["support_classes"]["overlap"]["closest_policy_by_median"], "active_edge_average")

    def test_pn2d_bv_active_edge_current_density_matches_sentaurus_currents(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_active_edge_current_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            fields.mkdir(parents=True)
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            charge = 1.602176634e-19

            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 8.0e-6, 4.0e-6, 1, 1, "overlap"]])
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [[0, 0.0, 0.0]])
            for name, value in {
                "ImpactIonization": 8.0e-6,
                "eCurrentDensity": charge * 3.0 / 1.0e4,
                "hCurrentDensity": charge * 1.0 / 1.0e4,
                "TotalCurrentDensity": charge * 4.0 / 1.0e4,
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [[0, value]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
            ], [
                [-13.2, 10, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 4.0],
                [-13.2, 11, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0, 4.0, 1.0, 1.0, 4.0],
                [-13.2, 12, 0, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_edge_current_density.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "active_edge_current_density_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["active_electron_flux_avg_m2_s"]), 3.0)
            self.assertAlmostEqual(float(row["active_hole_flux_avg_m2_s"]), 1.0)
            self.assertAlmostEqual(float(row["active_particle_flux_avg_m2_s"]), 4.0)
            self.assertAlmostEqual(float(row["active_generation_avg_m3_s"]), 8.0)
            self.assertAlmostEqual(float(row["active_weighted_alpha_m_inv"]), 2.0)
            self.assertAlmostEqual(float(row["active_particle_flux_over_sentaurus_eh_flux"]), 1.0)
            self.assertAlmostEqual(float(row["active_generation_over_sentaurus_generation"]), 1.0)
            self.assertAlmostEqual(float(row["active_particle_flux_over_sentaurus_total_current_equiv"]), 1.0)
            summary = json.loads((out_dir / "active_edge_current_density_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_active_edge_flux_factors_reconstructs_matching_states(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_active_edge_flux_factors_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            charge = 1.602176634e-19
            vt = 8.617333262145e-5 * 300.0
            expected_electron_flux = vt * 2.0 / 1.0e-6
            expected_hole_flux = vt * 1.0 / 1.0e-6

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                ],
            }), encoding="utf-8")
            self._write_csv(doping, ["node_id", "donors_m3", "acceptors_m3"], [
                [0, 0.0, 0.0],
                [1, 0.0, 0.0],
                [2, 0.0, 0.0],
            ])
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0, expected_electron_flux, expected_hole_flux, 1.0, 1.0],
                [-13.2, 2, 0, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
0 1 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
4 2 4
SCALARS Holes float 1
LOOKUP_TABLE default
3 2 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [0.0, 0.0, 0.0],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [3.0e-6, 2.0e-6, 3.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])
            for name, value in {
                "ImpactIonization": 1.0,
                "eCurrentDensity": charge * expected_electron_flux / 1.0e4,
                "hCurrentDensity": charge * expected_hole_flux / 1.0e4,
                "TotalCurrentDensity": charge * (expected_electron_flux + expected_hole_flux) / 1.0e4,
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [[0, value], [1, value], [2, value]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_edge_flux_factors.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--mesh",
                    str(mesh),
                    "--doping-csv",
                    str(doping),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--vela-vtk-root",
                    str(vtk_root),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                    "--bandgap-narrowing",
                    "none",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "active_edge_flux_factors_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["vela_density_particle_flux_over_sentaurus_state"]), 1.0)
            self.assertAlmostEqual(float(row["sentaurus_density_particle_flux_over_sentaurus_current"]), 1.0)
            self.assertAlmostEqual(float(row["cxx_particle_flux_over_vela_density_flux"]), 1.0)
            self.assertAlmostEqual(float(row["electron_density_geommean_vela_over_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(row["electron_mobility_vela_over_sentaurus"]), 1.0)
            summary = json.loads((out_dir / "active_edge_flux_factors_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_edge_source_current_consistency_writes_active_edge_rows(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_edge_source_current_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            charge = 1.602176634e-19

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 8.0e-6, 4.0e-6, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um",
                "x1_um", "y1_um", "edge_area_proxy_m2", "electron_flux_proxy",
                "hole_flux_proxy", "electron_alpha_m_inv", "hole_alpha_m_inv",
                "edge_source_integral",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 1.0, 2.0, 4.0, 16.0],
                [-13.2, 2, 0, 2, 0.0, 0.0, 0.0, 1.0, 2.0, 9.0, 9.0, 1.0, 1.0, 36.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
edge source current
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
4 2 4
SCALARS Holes float 1
LOOKUP_TABLE default
3 2 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [0.0, 0.0, 0.0],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [3.0e-6, 2.0e-6, 3.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
                "ImpactIonization": [8.0e-6, 8.0e-6, 1.0e-6],
                "eCurrentDensity": [charge * 2.0 / 1.0e4, charge * 2.0 / 1.0e4, 0.0],
                "hCurrentDensity": [charge * 1.0 / 1.0e4, charge * 1.0 / 1.0e4, 0.0],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_edge_source_current_consistency.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(edge_csv),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vtk_root),
                    "--bias", "-13.2",
                    "--out-dir", str(out_dir),
                    "--bandgap-narrowing", "none",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "edge_source_current_consistency_edges.csv")
            summary = json.loads((out_dir / "edge_source_current_consistency_summary.json").read_text())

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["edge_id"], "0")
        self.assertAlmostEqual(float(row["cxx_particle_flux_m2_s"]), 3.0)
        self.assertAlmostEqual(float(row["sentaurus_support_particle_flux_m2_s"]), 3.0)
        self.assertAlmostEqual(float(row["cxx_particle_flux_over_sentaurus_support_current"]), 1.0)
        self.assertIn("vela_density_particle_flux_m2_s", row)
        self.assertIn("sentaurus_qf_old_slotboom_particle_flux_m2_s", row)
        self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_edge_state_substitution_current_writes_variant_rows(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_edge_state_substitution_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            charge = 1.602176634e-19

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 0.0],
                [1, 0.0, 0.0],
                [2, 0.0, 0.0],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um",
                "x1_um", "y1_um", "edge_area_proxy_m2", "electron_flux_proxy",
                "hole_flux_proxy", "electron_alpha_m_inv", "hole_alpha_m_inv",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 1.0, 1.0, 1.0],
                [-13.2, 2, 0, 2, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 1.0, 1.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
edge state substitution
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
2 2 2
SCALARS Holes float 1
LOOKUP_TABLE default
1 1 1
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [0.0, 0.0, 0.0],
                "eDensity": [2.0e-6, 2.0e-6, 2.0e-6],
                "hDensity": [1.0e-6, 1.0e-6, 1.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
                "ImpactIonization": [1.0, 1.0, 1.0],
                "eCurrentDensity": [charge * 2.0 / 1.0e4, charge * 2.0 / 1.0e4, 0.0],
                "hCurrentDensity": [charge * 1.0 / 1.0e4, charge * 1.0 / 1.0e4, 0.0],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_edge_state_substitution_current.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(edge_csv),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vtk_root),
                    "--bias", "-13.2",
                    "--out-dir", str(out_dir),
                    "--bandgap-narrowing", "none",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "edge_state_substitution_current_edges.csv")
            summary = json.loads((out_dir / "edge_state_substitution_current_summary.json").read_text())

        variants = {row["variant"]: row for row in rows}
        self.assertEqual(set(variants), {
            "sentaurus_baseline",
            "sentaurus_psi_vela_qf",
            "vela_baseline",
            "vela_psi_sentaurus_qf",
            "vela_qf_shift",
            "vela_qf_shift_sentaurus_mobility",
        })
        self.assertAlmostEqual(
            float(variants["sentaurus_baseline"]["sentaurus_support_particle_flux_m2_s"]),
            3.0,
        )
        self.assertIn(
            "density_particle_flux_over_sentaurus_support_current",
            variants["sentaurus_baseline"],
        )
        self.assertIn(
            "qf_old_slotboom_particle_flux_over_sentaurus_support_current",
            variants["vela_qf_shift"],
        )
        self.assertEqual(summary["support_classes"]["overlap"]["variants"]["sentaurus_baseline"]["count"], 1)

    def test_pn2d_bv_qf_shift_continuity_balance_writes_focused_variants(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_qf_shift_balance_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            charge = 1.602176634e-19

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 0.0],
                [1, 0.0, 0.0],
                [2, 0.0, 0.0],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um",
                "x1_um", "y1_um", "edge_area_proxy_m2", "electron_flux_proxy",
                "hole_flux_proxy", "electron_alpha_m_inv", "hole_alpha_m_inv",
                "edge_source_integral",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 1.0, 1.0, 1.0, 6.0],
                [-13.2, 2, 0, 2, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 1.0, 1.0, 0.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
qf shift balance
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
2 2 2
SCALARS Holes float 1
LOOKUP_TABLE default
1 1 1
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
0 0 0
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [0.0, 0.0, 0.0],
                "eDensity": [2.0e-6, 2.0e-6, 2.0e-6],
                "hDensity": [1.0e-6, 1.0e-6, 1.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
                "ImpactIonization": [5.0e-6, 5.0e-6, 0.0],
                "srhRecombination": [0.0, 0.0, 0.0],
                "eCurrentDensity": [charge * 2.0 / 1.0e4, charge * 2.0 / 1.0e4, 0.0],
                "hCurrentDensity": [charge * 1.0 / 1.0e4, charge * 1.0 / 1.0e4, 0.0],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_qf_shift_continuity_balance.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(edge_csv),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vtk_root),
                    "--bias", "-13.2",
                    "--out-dir", str(out_dir),
                    "--bandgap-narrowing", "none",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "qf_shift_continuity_balance_nodes.csv")
            summary = json.loads((out_dir / "qf_shift_continuity_balance_summary.json").read_text())

        variants = {row["variant"]: row for row in rows}
        self.assertEqual(set(variants), {"vela_baseline", "vela_qf_shift", "vela_psi_sentaurus_qf"})
        self.assertAlmostEqual(float(variants["vela_baseline"]["sentaurus_active_current_integral_s_inv"]), 3.0)
        self.assertAlmostEqual(float(variants["vela_baseline"]["cxx_edge_source_s_inv"]), 3.0)
        self.assertIn("electron_residual_over_impact", variants["vela_qf_shift"])
        self.assertIn("particle_flux_over_sentaurus_support_current", variants["vela_psi_sentaurus_qf"])
        self.assertEqual(summary["support_classes"]["overlap"]["variants"]["vela_baseline"]["count"], 1)

    def test_pn2d_bv_exact_carrier_term_states_prepares_variant_configs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_exact_carrier_terms_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            base_config = root / "base.json"
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Cathode", "region_id": 0, "node_ids": [1]},
                    {"id": 1, "name": "Anode", "region_id": 0, "node_ids": [0]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 0.0],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            base_config.write_text(json.dumps({
                "mesh_file": str(mesh),
                "node_doping_file": str(doping),
                "scaling": {"mode": "unit_scaling"},
                "doping": [{"region": "Si", "donors": 1.0e17, "acceptors": 1.0e17}],
                "contacts": [
                    {"name": "Cathode", "bias": 0.0},
                    {"name": "Anode", "bias": -13.2},
                ],
                "solver": {
                    "method": "gummel_newton",
                    "bandgap_narrowing": "slotboom",
                    "mobility": {"model": "caughey_thomas_field"},
                    "impact_ionization": {
                        "model": "van_overstraeten",
                        "generation": "current_density",
                        "current_approximation": "density_gradient",
                    },
                },
            }, indent=2) + "\n")
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
exact carrier terms
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
2 2 2
SCALARS Holes float 1
LOOKUP_TABLE default
1 1 1
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.2, 0.0],
                "eQuasiFermiPotential": [-0.01, -0.01, -0.01],
                "hQuasiFermiPotential": [0.02, 0.02, 0.02],
                "eDensity": [2.0e-6, 2.0e-6, 2.0e-6],
                "hDensity": [1.0e-6, 1.0e-6, 1.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_exact_carrier_term_states.py"),
                    "--base-config", str(base_config),
                    "--support-csv", str(support),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vtk_root),
                    "--bias", "-13.2",
                    "--electron-qf-shift-v", "-0.0065",
                    "--hole-qf-shift-v", "0.0096",
                    "--electron-impact-scale", "1.24",
                    "--hole-impact-scale", "1.51",
                    "--preserve-contact-qf-on-shift",
                    "--out-dir", str(out_dir),
                    "--prepare-only",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((out_dir / "exact_carrier_term_state_summary.json").read_text())
            config = json.loads(
                (out_dir / "configs" / "vela_qf_shift_newton_carrier_term_probe.json").read_text()
            )
            shifted_phin = self._read_csv(
                out_dir / "states" / "vela_qf_shift" / "fields" / "eQuasiFermiPotential_region0.csv"
            )
            hybrid_h = self._read_csv(
                out_dir / "states" / "vela_psi_sentaurus_qf" / "fields" / "hQuasiFermiPotential_region0.csv"
            )

        self.assertEqual(summary["variants"], ["vela_baseline", "vela_psi_sentaurus_qf", "vela_qf_shift"])
        self.assertEqual(summary["selected_node_count"], 1)
        self.assertEqual(config["simulation_type"], "newton_carrier_term_probe")
        self.assertEqual(config["carrier_term_probe"]["electron_impact_scale"], 1.24)
        self.assertEqual(config["carrier_term_probe"]["hole_impact_scale"], 1.51)
        self.assertAlmostEqual(float(shifted_phin[0]["component0"]), 0.0)
        self.assertAlmostEqual(float(shifted_phin[2]["component0"]), -0.0065)
        self.assertAlmostEqual(float(hybrid_h[0]["component0"]), 0.02)

    def test_pn2d_bv_exact_carrier_term_states_exports_impact_source_components(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_exact_carrier_term_states.py"
        scripts_dir = str(REPO / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        spec = importlib.util.spec_from_file_location("diagnose_exact_carrier_terms", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        self.assertIn("impact_electron_source", module.NODE_FIELDS)
        self.assertIn("impact_hole_source", module.NODE_FIELDS)
        self.assertIn("impact_combined_source", module.NODE_FIELDS)

    def test_pn2d_bv_source_policy_matrix_ranks_policy_pairs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_source_policy_matrix_") as tmp:
            root = Path(tmp)
            exact = root / "exact_nodes.csv"
            full = root / "full_terms.csv"
            out_dir = root / "out"
            headers = [
                "variant",
                "bias_V",
                "node_id",
                "support_class",
                "electron_flux",
                "electron_recombination",
                "electron_impact",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_impact",
                "hole_gauge",
                "hole_boundary",
            ]
            self._write_csv(exact, headers, [
                ["vela_psi_sentaurus_qf", -13.2, 0, "overlap", 20.0, 0.0, -10.0, 0.0, 0.0, 13.0, 0.0, -10.0, 0.0, 0.0],
                ["vela_psi_sentaurus_qf", -13.2, 1, "overlap", 4.0, 0.0, -2.0, 0.0, 0.0, 2.6, 0.0, -2.0, 0.0, 0.0],
                ["other", -13.2, 2, "overlap", 999.0, 0.0, -10.0, 0.0, 0.0, 999.0, 0.0, -10.0, 0.0, 0.0],
            ])
            self._write_csv(full, [
                "node_id",
                "electron_flux",
                "electron_recombination",
                "electron_impact",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_impact",
                "hole_gauge",
                "hole_boundary",
                "impact_electron_source",
                "impact_hole_source",
                "impact_combined_source",
            ], [
                [0, 20.0, 0.0, -10.0, 0.0, 0.0, 13.0, 0.0, -10.0, 0.0, 0.0, 7.0, 3.0, 10.0],
                [1, 4.0, 0.0, -2.0, 0.0, 0.0, 2.6, 0.0, -2.0, 0.0, 0.0, 1.4, 0.6, 2.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_source_policy_matrix.py"),
                    "--exact-node-csv",
                    str(exact),
                    "--carrier-term-csv",
                    str(full),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((out_dir / "source_policy_matrix_summary.json").read_text())
            best_overlap = summary["best_by_scope"]["support:overlap"]
            self.assertEqual(best_overlap["electron_policy"], "double_combined")
            self.assertEqual(best_overlap["hole_policy"], "combined_plus_hole")
            self.assertAlmostEqual(best_overlap["combined_carrier_l2"], 0.0)
            best_global = summary["best_by_scope"]["global"]
            self.assertEqual(best_global["electron_policy"], "double_combined")
            self.assertEqual(best_global["hole_policy"], "combined_plus_hole")
            rows = self._read_csv(out_dir / "source_policy_matrix_summary.csv")
            self.assertTrue(any(row["scope"] == "support:overlap" for row in rows))
            self.assertTrue(any(row["scope"] == "global" for row in rows))

    def test_pn2d_bv_source_component_coefficients_merges_edge_features(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_source_coefficients_") as tmp:
            root = Path(tmp)
            exact = root / "exact_nodes.csv"
            full = root / "full_terms.csv"
            edge_features = root / "edge_features.csv"
            out_dir = root / "out"
            self._write_csv(exact, [
                "variant",
                "bias_V",
                "node_id",
                "x",
                "y",
                "support_class",
                "electron_flux",
                "electron_recombination",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_gauge",
                "hole_boundary",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 0, 1.0e-6, 2.0e-8, "overlap", 15.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0],
                ["vela_psi_sentaurus_qf", -13.2, 1, 1.0e-6, 4.0e-8, "overlap", 7.5, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0],
                ["other", -13.2, 2, 1.0e-6, 6.0e-8, "overlap", 999.0, 0.0, 0.0, 0.0, 999.0, 0.0, 0.0, 0.0],
            ])
            self._write_csv(full, [
                "node_id",
                "impact_electron_source",
                "impact_hole_source",
                "impact_combined_source",
            ], [
                [0, 8.0, 2.0, 10.0],
                [1, 4.0, 1.0, 5.0],
            ])
            self._write_csv(edge_features, [
                "node_id",
                "dominant_active_axis",
                "active_edge_count",
                "active_endpoint_area_fraction",
                "junction_normal_active_area_fraction",
                "junction_tangent_active_area_fraction",
            ], [
                [0, "x", 2, 0.5, 0.5, 0.0],
                [1, "y", 1, 0.25, 0.0, 0.25],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_source_component_coefficients.py"),
                    "--exact-node-csv",
                    str(exact),
                    "--carrier-term-csv",
                    str(full),
                    "--edge-feature-csv",
                    str(edge_features),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "source_component_coefficients_nodes.csv")
            self.assertEqual(len(rows), 2)
            node0 = next(row for row in rows if row["node_id"] == "0")
            self.assertEqual(node0["dominant_active_axis"], "x")
            self.assertAlmostEqual(float(node0["electron_required_over_combined"]), 1.5)
            self.assertAlmostEqual(float(node0["hole_required_over_combined"]), 1.2)
            self.assertAlmostEqual(float(node0["electron_min_norm_electron_coeff"]), 15.0 * 8.0 / 68.0)
            self.assertAlmostEqual(float(node0["electron_min_norm_hole_coeff"]), 15.0 * 2.0 / 68.0)
            summary = json.loads((out_dir / "source_component_coefficients_summary.json").read_text())
            overlap = summary["support_classes"]["overlap"]
            self.assertEqual(overlap["count"], 2)
            self.assertAlmostEqual(overlap["electron_required_over_combined_median"], 1.5)
            self.assertIn("axis:x", summary["feature_groups"])

    def test_pn2d_bv_edge_local_source_replay_groups_incident_edges(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_edge_local_replay_") as tmp:
            root = Path(tmp)
            exact = root / "exact_nodes.csv"
            edges = root / "sg_edges.csv"
            edge_current = root / "edge_source_current_consistency_edges.csv"
            out_dir = root / "out"
            self._write_csv(exact, [
                "variant",
                "bias_V",
                "node_id",
                "x",
                "y",
                "support_class",
                "electron_flux",
                "electron_recombination",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_gauge",
                "hole_boundary",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 10, 1.0e-6, 1.0e-6, "overlap", 30.0, 0.0, 0.0, 0.0, 18.0, 0.0, 0.0, 0.0],
            ])
            self._write_csv(edges, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "electron_source_integral",
                "hole_source_integral",
                "edge_source_integral",
            ], [
                [-13.2, 1, 1, 10, 0.0, 1.0, 1.0, 1.0, 8.0, 2.0, 10.0],
                [-13.2, 2, 10, 2, 1.0, 1.0, 2.0, 1.0, 4.0, 6.0, 10.0],
            ])
            self._write_csv(edge_current, [
                "bias_V",
                "support_node_id",
                "edge_id",
                "endpoint_area_m2",
                "cxx_particle_flux_m2_s",
                "cxx_endpoint_source_integral_s_inv",
                "sentaurus_support_particle_flux_m2_s",
                "sentaurus_edgeavg_particle_flux_m2_s",
                "sentaurus_edgeavg_electron_current_flux_m2_s",
                "sentaurus_edgeavg_hole_current_flux_m2_s",
                "sentaurus_support_generation_m3_s",
                "sentaurus_edgeavg_generation_m3_s",
            ], [
                [-13.2, 10, 1, 1.0, 10.0, 5.0, 20.0, 30.0, 6.0, 24.0, 15.0, 10.0],
                [-13.2, 10, 2, 1.0, 10.0, 5.0, 20.0, 30.0, 6.0, 24.0, 15.0, 10.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_edge_local_source_replay.py"),
                    "--exact-node-csv",
                    str(exact),
                    "--sg-edge-csv",
                    str(edges),
                    "--edge-current-csv",
                    str(edge_current),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "edge_local_source_replay_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["left_combined_source"]), 5.0)
            self.assertAlmostEqual(float(row["right_combined_source"]), 5.0)
            self.assertAlmostEqual(float(row["all_combined_source"]), 10.0)
            self.assertAlmostEqual(float(row["current_cxx_endpoint_source"]), 10.0)
            self.assertAlmostEqual(float(row["sentaurus_support_generation_source"]), 30.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_generation_source"]), 20.0)
            self.assertAlmostEqual(float(row["sentaurus_support_current_scaled_source"]), 20.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_current_scaled_source"]), 30.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_electron_current_scaled_source"]), 6.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_hole_current_scaled_source"]), 24.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_abs_electron_current_scaled_source"]), 6.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_abs_hole_current_scaled_source"]), 24.0)
            self.assertAlmostEqual(float(row["electron_required_over_sentaurus_support_generation_source"]), 1.0)
            self.assertAlmostEqual(float(row["hole_required_over_sentaurus_edgeavg_generation_source"]), 0.9)
            self.assertAlmostEqual(float(row["hybrid_e2cxx_h_edgeavg_current_electron_residual"]), 10.0)
            self.assertAlmostEqual(float(row["hybrid_e2cxx_h_edgeavg_current_hole_residual"]), -12.0)
            self.assertAlmostEqual(float(row["hybrid_e2cxx_h_edgeavg_current_residual_l2"]), math.sqrt(244.0))
            self.assertAlmostEqual(float(row["electron_double_all_residual"]), 10.0)
            self.assertAlmostEqual(float(row["hole_combined_plus_hole_residual"]), 4.0)
            summary = json.loads((out_dir / "edge_local_source_replay_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)
            self.assertAlmostEqual(
                summary["support_classes"]["overlap"]["hybrid_e2cxx_h_edgeavg_current_combined_l2"],
                math.sqrt(244.0),
            )
            self.assertAlmostEqual(
                summary["support_classes"]["overlap"]["hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined_median"],
                1.0,
            )
            self.assertAlmostEqual(
                summary["support_classes"]["overlap"]["hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined_median"],
                -1.2,
            )
            self.assertAlmostEqual(
                summary["support_classes"]["overlap"]["electron_double_all_residual_over_combined_median"],
                1.0,
            )

    def test_pn2d_bv_exact_carrier_term_states_decomposes_srh(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_exact_carrier_term_states.py"
        scripts_dir = str(REPO / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        spec = importlib.util.spec_from_file_location("diagnose_exact_carrier_terms", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        result = module.srh_decomposition(
            psi=0.12,
            phin=0.05,
            phip=0.16,
            ni_eff=2.0e16,
            vt=0.025,
            node_volume=3.0e-14,
            taun=1.0e-5,
            taup=3.0e-6,
            continuity_scale=2.0,
            exact_recombination=0.0,
            exact_impact=-4.0,
        )

        n = 2.0e16 * math.exp((0.12 - 0.05) / 0.025)
        p = 2.0e16 * math.exp((0.16 - 0.12) / 0.025)
        excess = 2.0e16 * 2.0e16 * math.expm1((0.16 - 0.05) / 0.025)
        denominator = 3.0e-6 * (n + 2.0e16) + 1.0e-5 * (p + 2.0e16)
        rate = excess / denominator
        integral = rate * 3.0e-14 / 2.0
        self.assertAlmostEqual(result["electron_density_m3"], n)
        self.assertAlmostEqual(result["hole_density_m3"], p)
        self.assertAlmostEqual(result["electron_density_over_ni"], n / 2.0e16)
        self.assertAlmostEqual(result["hole_density_over_ni"], p / 2.0e16)
        self.assertAlmostEqual(result["srh_excess_product_m6"], excess)
        self.assertAlmostEqual(result["srh_excess_over_ni2"], excess / (2.0e16 * 2.0e16))
        self.assertAlmostEqual(result["srh_taup_n_plus_ni_m3_s"], 3.0e-6 * (n + 2.0e16))
        self.assertAlmostEqual(result["srh_taun_p_plus_ni_m3_s"], 1.0e-5 * (p + 2.0e16))
        self.assertAlmostEqual(result["srh_denominator_m3_s"], denominator)
        self.assertAlmostEqual(result["srh_rate_m3_s"], rate)
        self.assertAlmostEqual(result["srh_recombination_integral_scaled_s_inv"], integral)
        self.assertAlmostEqual(result["exact_recombination_over_abs_impact"], 0.0)

    def test_pn2d_bv_active_overlap_outlier_table_merges_state_flux_and_residuals(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_overlap_outliers_") as tmp:
            root = Path(tmp)
            exact = root / "exact_nodes.csv"
            replay = root / "edge_replay.csv"
            consistency = root / "edge_current.csv"
            mobility = root / "mobility.csv"
            mixed = root / "mixed.csv"
            residual = root / "residual_proxy.csv"
            out = root / "out"
            self._write_csv(exact, [
                "variant",
                "bias_V",
                "node_id",
                "x",
                "y",
                "support_class",
                "electron_flux",
                "electron_recombination",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_gauge",
                "hole_boundary",
                "psi_V",
                "phin_V",
                "phip_V",
                "electron_density_m3",
                "hole_density_m3",
                "ni_eff_m3",
                "node_volume_m2",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 10, 1.0e-6, 2.0e-8, "overlap", 30.0, 1.0, 2.0, 3.0, 18.0, 0.5, 1.5, 2.5, -1.0, -0.9, -1.1, 1.0e20, 2.0e18, 1.0e16, 4.0e-16],
                ["other", -13.2, 11, 1.0e-6, 1.0e-8, "overlap", 999.0, 0.0, 0.0, 0.0, 999.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            ])
            self._write_csv(replay, [
                "variant",
                "bias_V",
                "node_id",
                "x",
                "y",
                "support_class",
                "all_combined_source",
                "sentaurus_edgeavg_current_scaled_source",
                "hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined",
                "hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined",
                "hybrid_e2cxx_h_edgeavg_current_residual_l2",
                "sentaurus_edgeavg_electron_current_scaled_source",
                "sentaurus_edgeavg_hole_current_scaled_source",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 10, 1.0e-6, 2.0e-8, "overlap", 10.0, 22.0, 0.25, -0.5, 5.0, 9.0, 13.0],
            ])
            self._write_csv(consistency, [
                "bias_V",
                "support_node_id",
                "edge_id",
                "cxx_particle_flux_m2_s",
                "sentaurus_edgeavg_particle_flux_m2_s",
                "sentaurus_support_particle_flux_m2_s",
                "sentaurus_edgeavg_electron_current_flux_m2_s",
                "sentaurus_edgeavg_hole_current_flux_m2_s",
                "cxx_generation_density_m3_s",
                "sentaurus_edgeavg_generation_m3_s",
                "sentaurus_support_generation_m3_s",
                "cxx_generation_over_sentaurus_edgeavg_generation",
                "cxx_particle_flux_over_sentaurus_edgeavg_current",
            ], [
                [-13.2, 10, 1, 20.0, 40.0, 50.0, 16.0, 24.0, 100.0, 120.0, 140.0, 0.8, 0.5],
                [-13.2, 10, 2, 30.0, 60.0, 70.0, 24.0, 36.0, 200.0, 240.0, 280.0, 0.9, 0.6],
            ])
            self._write_csv(mobility, [
                "node_id",
                "electron_final_over_sentaurus_mobility",
                "hole_final_over_sentaurus_mobility",
                "electric_field_edge_abs_vela_over_sentaurus",
                "electron_drive_vela_over_sentaurus_qf",
                "hole_drive_vela_over_sentaurus_qf",
            ], [
                [10, 0.7, 0.8, 1.1, 1.2, 1.3],
            ])
            self._write_csv(mixed, [
                "variant",
                "node_id",
                "electron_density_geommean_over_sentaurus",
                "hole_density_geommean_over_sentaurus",
                "particle_flux_over_sentaurus",
                "generation_over_sentaurus",
            ], [
                ["vela_psi_sentaurus_qf", 10, 0.9, 1.1, 1.2, 1.3],
            ])
            self._write_csv(residual, [
                "variant",
                "node_id",
                "electron_residual_over_impact",
                "hole_residual_over_impact",
            ], [
                ["sentaurus_state_sentaurus_node_source", 10, 0.14, 0.17],
                ["shifted_vela_qf_sentaurus_node_source", 10, -0.09, 0.03],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_overlap_outlier_table.py"),
                    "--exact-node-csv",
                    str(exact),
                    "--edge-replay-csv",
                    str(replay),
                    "--edge-current-csv",
                    str(consistency),
                    "--mobility-csv",
                    str(mobility),
                    "--mixed-state-csv",
                    str(mixed),
                    "--residual-proxy-csv",
                    str(residual),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "active_overlap_outlier_table.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["node_id"], "10")
            self.assertAlmostEqual(float(row["electron_required_source_s_inv"]), 36.0)
            self.assertAlmostEqual(float(row["hole_required_source_s_inv"]), 22.5)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_particle_flux_m2_s_mean"]), 50.0)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_electron_current_fraction_mean"]), 0.4)
            self.assertAlmostEqual(float(row["sentaurus_edgeavg_hole_current_fraction_mean"]), 0.6)
            self.assertAlmostEqual(float(row["electron_final_over_sentaurus_mobility"]), 0.7)
            self.assertAlmostEqual(float(row["mixed_particle_flux_over_sentaurus"]), 1.2)
            self.assertAlmostEqual(float(row["sentaurus_state_electron_residual_over_impact"]), 0.14)
            summary = json.loads((out / "active_overlap_outlier_summary.json").read_text())
            self.assertEqual(summary["row_count"], 1)
            self.assertEqual(summary["max_hybrid_residual_node_id"], "10")

    def test_pn2d_bv_continuity_scaling_decomposition_compares_proxy_to_exact_units(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_continuity_scaling_") as tmp:
            root = Path(tmp)
            exact = root / "exact_nodes.csv"
            replay = root / "edge_replay.csv"
            residual = root / "residual_proxy.csv"
            out = root / "out"
            self._write_csv(exact, [
                "variant",
                "bias_V",
                "node_id",
                "x",
                "y",
                "support_class",
                "electron_flux",
                "electron_recombination",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_gauge",
                "hole_boundary",
                "node_volume_m2",
                "inferred_continuity_scale_from_srh",
            ], [
                [
                    "vela_psi_sentaurus_qf",
                    -13.2,
                    10,
                    1.0e-6,
                    2.0e-8,
                    "overlap",
                    30.0,
                    1.0,
                    2.0,
                    3.0,
                    18.0,
                    0.5,
                    1.5,
                    2.5,
                    4.0e-16,
                    10.0,
                ],
                [
                    "other",
                    -13.2,
                    11,
                    1.0e-6,
                    4.0e-8,
                    "overlap",
                    999.0,
                    0.0,
                    0.0,
                    0.0,
                    999.0,
                    0.0,
                    0.0,
                    0.0,
                    4.0e-16,
                    10.0,
                ],
            ])
            self._write_csv(replay, [
                "variant",
                "bias_V",
                "node_id",
                "support_class",
                "all_combined_source",
                "sentaurus_edgeavg_current_scaled_source",
                "electron_required_source",
                "hole_required_source",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 10, "overlap", 20.0, 30.0, 36.0, 22.5],
            ])
            self._write_csv(residual, [
                "variant",
                "node_id",
                "support_class",
                "transport_model",
                "impact_source_s_inv",
                "electron_transport_s_inv",
                "hole_transport_s_inv",
                "srh_source_s_inv",
                "electron_residual_proxy_s_inv",
                "hole_residual_proxy_s_inv",
            ], [
                [
                    "probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source",
                    10,
                    "overlap",
                    "qf_old_slotboom_ni",
                    3000.0,
                    3600.0,
                    3300.0,
                    300.0,
                    900.0,
                    600.0,
                ],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_continuity_scaling_decomposition.py"),
                    "--exact-node-csv",
                    str(exact),
                    "--edge-replay-csv",
                    str(replay),
                    "--residual-proxy-csv",
                    str(residual),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--proxy-variant",
                    "probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source",
                    "--transport-model",
                    "qf_old_slotboom_ni",
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "continuity_scaling_decomposition_nodes.csv")
            summary = json.loads((out / "continuity_scaling_decomposition_summary.json").read_text())

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["node_id"], "10")
        self.assertAlmostEqual(float(row["exact_edgeavg_current_over_proxy_impact_scale"]), 0.01)
        self.assertAlmostEqual(float(row["exact_electron_required_over_proxy_required_scale"]), 36.0 / 3900.0)
        self.assertAlmostEqual(float(row["exact_hole_required_over_proxy_required_scale"]), 22.5 / 3600.0)
        self.assertAlmostEqual(
            float(row["electron_required_scale_over_edgeavg_scale"]),
            (36.0 / 3900.0) / 0.01,
        )
        self.assertAlmostEqual(
            float(row["hole_required_scale_over_edgeavg_scale"]),
            (22.5 / 3600.0) / 0.01,
        )
        self.assertAlmostEqual(float(row["inverse_inferred_continuity_scale_from_srh"]), 0.1)
        self.assertEqual(summary["row_count"], 1)
        self.assertAlmostEqual(
            summary["metrics"]["electron_required_scale_over_edgeavg_scale"]["median"],
            (36.0 / 3900.0) / 0.01,
        )

    def test_pn2d_bv_source_ownership_replay_ranks_sentaurus_branch_sources(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_source_ownership_replay_") as tmp:
            root = Path(tmp)
            edge_replay = root / "edge_local_source_replay.csv"
            out_dir = root / "out"
            self._write_csv(edge_replay, [
                "variant",
                "bias_V",
                "node_id",
                "y",
                "support_class",
                "all_electron_source",
                "all_hole_source",
                "all_combined_source",
                "current_cxx_endpoint_source",
                "sentaurus_edgeavg_current_scaled_source",
                "sentaurus_edgeavg_electron_current_scaled_source",
                "sentaurus_edgeavg_hole_current_scaled_source",
                "sentaurus_edgeavg_abs_electron_current_scaled_source",
                "sentaurus_edgeavg_abs_hole_current_scaled_source",
                "electron_required_source",
                "hole_required_source",
            ], [
                [
                    "vela_psi_sentaurus_qf", -13.2, 10, 0.0, "overlap",
                    4.0, 6.0, 10.0, 10.0, 12.0, 5.0, 7.0, 5.0, 7.0,
                    5.0, 7.0,
                ],
                [
                    "vela_psi_sentaurus_qf", -13.2, 11, 1.0, "overlap",
                    8.0, 12.0, 20.0, 20.0, 24.0, 10.0, 14.0, 10.0, 14.0,
                    10.0, 14.0,
                ],
                [
                    "other", -13.2, 12, 2.0, "overlap",
                    8.0, 12.0, 20.0, 20.0, 24.0, 999.0, 999.0, 999.0, 999.0,
                    10.0, 14.0,
                ],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_source_ownership_replay.py"),
                    "--edge-replay-csv",
                    str(edge_replay),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--support-class",
                    "overlap",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((out_dir / "source_ownership_replay_summary.json").read_text())
            best = summary["best_by_combined_l2"]
            self.assertEqual(best["electron_policy"], "sentaurus_edgeavg_electron_current")
            self.assertEqual(best["hole_policy"], "sentaurus_edgeavg_hole_current")
            self.assertAlmostEqual(best["combined_carrier_l2"], 0.0)
            self.assertEqual(summary["row_count"], 2)
            rows = self._read_csv(out_dir / "source_ownership_replay_summary.csv")
            best_row = min(rows, key=lambda row: float(row["combined_carrier_l2"]))
            self.assertEqual(best_row["electron_policy"], "sentaurus_edgeavg_electron_current")
            self.assertEqual(best_row["hole_policy"], "sentaurus_edgeavg_hole_current")

    def test_pn2d_bv_row_term_decomposition_identifies_flux_y_trend(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_row_term_decomp_") as tmp:
            root = Path(tmp)
            exact = root / "exact_nodes.csv"
            edge_replay = root / "edge_replay.csv"
            out = root / "out"
            self._write_csv(exact, [
                "variant",
                "bias_V",
                "node_id",
                "x",
                "y",
                "support_class",
                "electron_flux",
                "electron_recombination",
                "electron_impact",
                "electron_gauge",
                "electron_boundary",
                "hole_flux",
                "hole_recombination",
                "hole_impact",
                "hole_gauge",
                "hole_boundary",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 1, 0.0, 0.0, "overlap", 8.0, -1.0, -5.0, 0.0, 0.0, 6.0, -1.0, -5.0, 0.0, 0.0],
                ["vela_psi_sentaurus_qf", -13.2, 2, 0.0, 1.0, "overlap", 8.0, -1.0, -5.0, 0.0, 0.0, 9.0, -1.0, -5.0, 0.0, 0.0],
                ["vela_psi_sentaurus_qf", -13.2, 3, 0.0, 2.0, "overlap", 8.0, -1.0, -5.0, 0.0, 0.0, 12.0, -1.0, -5.0, 0.0, 0.0],
                ["other", -13.2, 4, 0.0, 3.0, "overlap", 99.0, -1.0, -5.0, 0.0, 0.0, 99.0, -1.0, -5.0, 0.0, 0.0],
            ])
            self._write_csv(edge_replay, [
                "variant",
                "bias_V",
                "node_id",
                "support_class",
                "all_combined_source",
                "sentaurus_edgeavg_current_scaled_source",
            ], [
                ["vela_psi_sentaurus_qf", -13.2, 1, "overlap", 5.0, 5.5],
                ["vela_psi_sentaurus_qf", -13.2, 2, "overlap", 5.0, 5.5],
                ["vela_psi_sentaurus_qf", -13.2, 3, "overlap", 5.0, 5.5],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_row_term_decomposition.py"),
                    "--exact-node-csv",
                    str(exact),
                    "--edge-replay-csv",
                    str(edge_replay),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--support-class",
                    "overlap",
                    "--source-policy",
                    "all_combined_source",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "row_term_decomposition_nodes.csv")
            self.assertEqual(len(rows), 3)
            self.assertAlmostEqual(float(rows[-1]["hole_required_over_source"]), 11.0 / 5.0)
            summary = json.loads((out / "row_term_decomposition_summary.json").read_text())
            self.assertEqual(summary["row_count"], 3)
            self.assertEqual(summary["dominant_y_correlation"]["hole"], "hole_flux")
            self.assertAlmostEqual(summary["correlations"]["hole_flux"]["corr_y"], 1.0)
            self.assertAlmostEqual(summary["correlations"]["hole_residual_under_policy"]["corr_y"], 1.0)
            self.assertIsNone(summary["correlations"]["hole_recombination"]["corr_y"])

    def test_pn2d_bv_edge_flux_y_trend_localizes_hole_flux_ratio(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_edge_flux_y_trend_") as tmp:
            root = Path(tmp)
            edges = root / "edge_state_substitution_current_edges.csv"
            out = root / "out"
            self._write_csv(edges, [
                "bias_V",
                "variant",
                "support_node_id",
                "support_class",
                "y_um",
                "edge_id",
                "edge_axis",
                "sentaurus_edgeavg_electron_current_flux_m2_s",
                "sentaurus_edgeavg_hole_current_flux_m2_s",
                "sentaurus_edgeavg_particle_flux_m2_s",
                "qf_old_slotboom_electron_flux_m2_s",
                "qf_old_slotboom_hole_flux_m2_s",
                "qf_old_slotboom_particle_flux_m2_s",
            ], [
                [-13.2, "vela_psi_sentaurus_qf", 10, "overlap", 0.0, 100, "x", 10.0, 20.0, 30.0, 10.0, 10.0, 20.0],
                [-13.2, "vela_psi_sentaurus_qf", 10, "overlap", 0.0, 101, "x", 10.0, 20.0, 30.0, 10.0, 10.0, 20.0],
                [-13.2, "vela_psi_sentaurus_qf", 11, "overlap", 1.0, 102, "x", 10.0, 20.0, 30.0, 10.0, 20.0, 30.0],
                [-13.2, "vela_psi_sentaurus_qf", 11, "overlap", 1.0, 103, "x", 10.0, 20.0, 30.0, 10.0, 20.0, 30.0],
                [-13.2, "other", 12, "overlap", 2.0, 104, "y", 10.0, 20.0, 30.0, 99.0, 99.0, 99.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_edge_flux_y_trend.py"),
                    "--edge-state-csv",
                    str(edges),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--support-class",
                    "overlap",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "edge_flux_y_trend_nodes.csv")
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["edge_axes"], "x")
            self.assertAlmostEqual(float(rows[0]["hole_qf_over_sentaurus_edgeavg_median"]), 0.5)
            self.assertAlmostEqual(float(rows[1]["hole_qf_over_sentaurus_edgeavg_median"]), 1.0)
            summary = json.loads((out / "edge_flux_y_trend_summary.json").read_text())
            self.assertEqual(summary["row_count"], 2)
            self.assertEqual(summary["axis_counts"], {"x": 4})
            self.assertAlmostEqual(summary["correlations"]["hole_qf_over_sentaurus_edgeavg_median"]["corr_y"], 1.0)

    def test_pn2d_bv_sg_carrier_split_decomposition_merges_edge_state_and_flux_forms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_sg_split_decomp_") as tmp:
            root = Path(tmp)
            edge_state = root / "edge_state.csv"
            flux_forms = root / "sg_flux_forms.csv"
            out = root / "out"
            self._write_csv(edge_state, [
                "bias_V",
                "variant",
                "support_node_id",
                "support_class",
                "y_um",
                "edge_id",
                "edge_axis",
                "sentaurus_edgeavg_electron_current_flux_m2_s",
                "sentaurus_edgeavg_hole_current_flux_m2_s",
                "sentaurus_edgeavg_particle_flux_m2_s",
                "qf_old_slotboom_electron_flux_m2_s",
                "qf_old_slotboom_hole_flux_m2_s",
                "qf_old_slotboom_particle_flux_m2_s",
            ], [
                [-13.2, "vela_psi_sentaurus_qf", 10, "overlap", 0.0, 100, "x", 8.0, 2.0, 10.0, 6.0, 4.0, 10.0],
                [-13.2, "vela_psi_sentaurus_qf", 11, "overlap", 1.0, 101, "x", 2.0, 8.0, 10.0, 4.0, 6.0, 10.0],
            ])
            self._write_csv(flux_forms, [
                "bias_V",
                "edge_id",
                "source",
                "node0",
                "node1",
                "psi0_V",
                "psi1_V",
                "phin0_V",
                "phin1_V",
                "phip0_V",
                "phip1_V",
                "ni_model_0_m3",
                "ni_model_1_m3",
                "electron_density_0_m3",
                "electron_density_1_m3",
                "hole_density_0_m3",
                "hole_density_1_m3",
                "electron_flux_qf_model_abs",
                "hole_flux_qf_model_abs",
            ], [
                [-13.2, 100, "vela", 0, 1, 0.0, 0.1, -0.2, -0.1, 0.2, 0.3, 1.0e16, 2.0e16, 1.0, 2.0, 3.0, 4.0, 7.0, 3.0],
                [-13.2, 100, "sentaurus", 0, 1, 0.01, 0.11, -0.21, -0.12, 0.21, 0.31, 1.0e16, 2.0e16, 1.5, 2.5, 3.5, 4.5, 8.0, 2.0],
                [-13.2, 101, "vela", 0, 1, 0.0, 0.1, -0.2, -0.1, 0.2, 0.3, 1.0e16, 2.0e16, 1.0, 2.0, 3.0, 4.0, 3.0, 7.0],
                [-13.2, 101, "sentaurus", 0, 1, 0.01, 0.11, -0.21, -0.12, 0.21, 0.31, 1.0e16, 2.0e16, 1.5, 2.5, 3.5, 4.5, 2.0, 8.0],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_carrier_split_decomposition.py"),
                    "--edge-state-csv",
                    str(edge_state),
                    "--sg-flux-form-csv",
                    str(flux_forms),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--support-class",
                    "overlap",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "sg_carrier_split_decomposition_edges.csv")
            self.assertEqual(len(rows), 2)
            self.assertAlmostEqual(float(rows[0]["mixed_electron_fraction"]), 0.6)
            self.assertAlmostEqual(float(rows[0]["sentaurus_electron_fraction"]), 0.8)
            self.assertAlmostEqual(float(rows[0]["electron_fraction_delta"]), -0.2)
            self.assertIn("mixed_electron_eta", rows[0])
            summary = json.loads((out / "sg_carrier_split_decomposition_summary.json").read_text())
            self.assertEqual(summary["row_count"], 2)
            self.assertAlmostEqual(summary["correlations"]["electron_fraction_delta"]["corr_y"], 1.0)
            self.assertAlmostEqual(summary["correlations"]["hole_fraction_delta"]["corr_y"], -1.0)

    def test_pn2d_bv_potential_qf_alignment_replay_fits_constant_psi_shift(self) -> None:
        scripts_dir = str(REPO / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag

        with tempfile.TemporaryDirectory(prefix="vela_bv_potential_qf_replay_") as tmp:
            root = Path(tmp)
            edge_state = root / "edge_state.csv"
            flux_forms = root / "sg_flux_forms.csv"
            out = root / "out"
            vt = sgdiag.K_B_OVER_Q * 300.0
            ni = 1.0e16
            mobility = 1.0e-3
            length = 1.0e-6
            coef = mobility * vt / length

            def component_flux(psi0: float, psi1: float, phin0: float, phin1: float,
                               phip0: float, phip1: float) -> tuple[float, float, float]:
                electron = abs(sgdiag.sg_electron_flux_qf_variable_ni(
                    ni, ni, psi0, psi1, phin0, phin1, vt, coef))
                hole = abs(sgdiag.sg_hole_flux_qf_variable_ni(
                    ni, ni, psi0, psi1, phip0, phip1, vt, coef))
                return electron, hole, electron + hole

            sent_edge_100 = component_flux(0.0, 0.02, -0.20, -0.18, 0.20, 0.22)
            sent_edge_101 = component_flux(0.0, 0.02, -0.19, -0.17, 0.21, 0.23)

            self._write_csv(edge_state, [
                "bias_V",
                "variant",
                "support_node_id",
                "support_class",
                "y_um",
                "edge_id",
                "edge_axis",
                "sentaurus_edgeavg_electron_current_flux_m2_s",
                "sentaurus_edgeavg_hole_current_flux_m2_s",
                "sentaurus_edgeavg_particle_flux_m2_s",
            ], [
                [-13.2, "vela_psi_sentaurus_qf", 10, "overlap", 0.0, 100, "x", *sent_edge_100],
                [-13.2, "vela_psi_sentaurus_qf", 11, "overlap", 1.0, 101, "x", *sent_edge_101],
            ])
            self._write_csv(flux_forms, [
                "bias_V",
                "edge_id",
                "source",
                "node0",
                "node1",
                "length_m",
                "psi0_V",
                "psi1_V",
                "phin0_V",
                "phin1_V",
                "phip0_V",
                "phip1_V",
                "ni_model_0_m3",
                "ni_model_1_m3",
                "electron_mobility_m2_V_s",
                "hole_mobility_m2_V_s",
            ], [
                [-13.2, 100, "vela", 0, 1, length, 0.05, 0.07, -0.20, -0.18, 0.20, 0.22, ni, ni, mobility, mobility],
                [-13.2, 100, "sentaurus", 0, 1, length, 0.00, 0.02, -0.20, -0.18, 0.20, 0.22, ni, ni, mobility, mobility],
                [-13.2, 101, "vela", 0, 1, length, 0.05, 0.07, -0.19, -0.17, 0.21, 0.23, ni, ni, mobility, mobility],
                [-13.2, 101, "sentaurus", 0, 1, length, 0.00, 0.02, -0.19, -0.17, 0.21, 0.23, ni, ni, mobility, mobility],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_potential_qf_alignment_replay.py"),
                    "--edge-state-csv",
                    str(edge_state),
                    "--sg-flux-form-csv",
                    str(flux_forms),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--support-class",
                    "overlap",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "potential_qf_alignment_replay_edges.csv")
            self.assertEqual({row["replay_variant"] for row in rows}, {
                "sentaurus_psi_sentaurus_qf_sentaurus_mobility",
                "sentaurus_psi_sentaurus_qf_vela_mobility",
                "vela_psi_sentaurus_qf_sentaurus_mobility",
                "vela_psi_sentaurus_qf_vela_mobility",
                "vela_psi_shifted_const_sentaurus_qf_vela_mobility",
            })
            summary = json.loads((out / "potential_qf_alignment_replay_summary.json").read_text())
            baseline = summary["variants"]["vela_psi_sentaurus_qf_vela_mobility"]
            pure = summary["variants"]["sentaurus_psi_sentaurus_qf_sentaurus_mobility"]
            shifted = summary["variants"]["vela_psi_shifted_const_sentaurus_qf_vela_mobility"]
            self.assertLess(pure["combined_fraction_l2"], baseline["combined_fraction_l2"])
            self.assertLess(shifted["combined_fraction_l2"], baseline["combined_fraction_l2"])
            self.assertAlmostEqual(shifted["psi_shift_V"], -0.05, places=3)
            node_shift_rows = self._read_csv(out / "potential_qf_alignment_node_shift_scan.csv")
            self.assertEqual(len(node_shift_rows), 2)
            self.assertAlmostEqual(float(node_shift_rows[0]["best_psi_shift_V"]), -0.05, places=3)
            self.assertAlmostEqual(
                float(node_shift_rows[0]["direct_sentaurus_minus_vela_psi_shift_median_V"]),
                -0.05,
                places=6,
            )
            self.assertEqual(summary["node_shift_scan"]["row_count"], 2)

    def test_pn2d_bv_psi_field_alignment_reports_neighbor_residuals(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_psi_field_alignment_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            edge_state = root / "edge_state.csv"
            vela_vtk = root / "vela_-13.2V.vtk"
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            out = root / "out"
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                    {"id": 1, "node_ids": [1, 3, 2]},
                ],
                "contacts": [{"name": "anode", "node_ids": [0]}],
            }))
            self._write_csv(doping, [
                "node_id", "donors_m3", "acceptors_m3",
            ], [
                [0, 1.0e21, 0.0],
                [1, 1.0e21, 0.0],
                [2, 0.0, 2.0e21],
                [3, 0.0, 2.0e21],
            ])
            self._write_csv(edge_state, [
                "bias_V", "variant", "support_node_id", "support_class", "y_um", "edge_id",
            ], [
                [-13.2, "vela_psi_sentaurus_qf", 0, "overlap", 0.0, ""],
            ])

            def vtk_scalar(name: str, values: list[float]) -> str:
                body = "\n".join(str(value) for value in values)
                return f"SCALARS {name} double 1\nLOOKUP_TABLE default\n{body}\n"

            vela_vtk.write_text(
                "# vtk DataFile Version 3.0\npsi alignment fixture\nASCII\nDATASET POLYDATA\n"
                "POINT_DATA 4\n"
                + vtk_scalar("Potential", [0.14, 0.15, 0.12, 0.13])
                + vtk_scalar("ElectronQuasiFermi", [0.0, 0.0, 0.0, 0.0])
                + vtk_scalar("HoleQuasiFermi", [0.0, 0.0, 0.0, 0.0])
                + vtk_scalar("Electrons", [1.0e16, 1.0e16, 2.0e16, 2.0e16])
                + vtk_scalar("Holes", [2.0e16, 2.0e16, 1.0e16, 1.0e16])
                + vtk_scalar("ElectronMobility", [0.1, 0.1, 0.1, 0.1])
                + vtk_scalar("HoleMobility", [0.05, 0.05, 0.05, 0.05])
            )

            def write_field(name: str, values: list[float]) -> None:
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", name], [
                    [index, value] for index, value in enumerate(values)
                ])

            write_field("ElectrostaticPotential", [0.0, 0.01, -0.01, 0.0])
            write_field("eQuasiFermiPotential", [0.0, 0.0, 0.0, 0.0])
            write_field("hQuasiFermiPotential", [0.0, 0.0, 0.0, 0.0])
            write_field("eDensity", [1.0e10, 1.0e10, 2.0e10, 2.0e10])
            write_field("hDensity", [2.0e10, 2.0e10, 1.0e10, 1.0e10])
            write_field("eMobility", [1000.0, 1000.0, 1000.0, 1000.0])
            write_field("hMobility", [500.0, 500.0, 500.0, 500.0])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_psi_field_alignment.py"),
                    "--mesh",
                    str(mesh),
                    "--doping-csv",
                    str(doping),
                    "--edge-state-csv",
                    str(edge_state),
                    "--vela-vtk",
                    str(vela_vtk),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--variant",
                    "vela_psi_sentaurus_qf",
                    "--bias",
                    "-13.2",
                    "--support-class",
                    "overlap",
                    "--neighbor-depth",
                    "1",
                    "--out-dir",
                    str(out),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "psi_field_alignment_nodes.csv")
            by_node = {row["node_id"]: row for row in rows}
            self.assertEqual(by_node["0"]["node_class"], "active_support")
            self.assertEqual(by_node["0"]["contact_names"], "anode")
            self.assertAlmostEqual(float(by_node["0"]["psi_delta_V"]), -0.14)
            self.assertAlmostEqual(float(by_node["0"]["psi_delta_residual_after_active_median_V"]), 0.0)
            self.assertEqual(by_node["2"]["node_class"], "neighbor_depth_1")
            self.assertAlmostEqual(float(by_node["2"]["psi_delta_residual_after_active_median_V"]), 0.01)
            self.assertIn("ni_model_m3", by_node["0"])
            self.assertIn("vela_inferred_ni_hole_m3", by_node["0"])
            self.assertIn("sentaurus_inferred_ni_hole_m3", by_node["0"])
            summary = json.loads((out / "psi_field_alignment_summary.json").read_text())
            self.assertEqual(summary["row_count"], 3)
            self.assertAlmostEqual(summary["active_seed_psi_delta_median_V"], -0.14)

    def test_pn2d_bv_exact_carrier_term_states_computes_impact_closure(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_exact_carrier_term_states.py"
        scripts_dir = str(REPO / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        spec = importlib.util.spec_from_file_location("diagnose_exact_carrier_terms", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        result = module.impact_closure(
            flux=10.0,
            recombination=-3.0,
            impact=-5.0,
            residual=2.0,
        )

        self.assertAlmostEqual(result["required_impact_multiplier"], 1.4)
        self.assertAlmostEqual(result["impact_multiplier_delta"], 0.4)
        self.assertAlmostEqual(result["residual_over_abs_impact"], 0.4)
        self.assertAlmostEqual(result["closed_residual"], 0.0)

    def test_pn2d_bv_active_edge_mobility_inputs_decompose_limiter(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_active_edge_mobility_inputs_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                ],
            }), encoding="utf-8")
            self._write_csv(doping, ["node_id", "donors_m3", "acceptors_m3"], [
                [0, 0.0, 0.0],
                [1, 0.0, 0.0],
                [2, 0.0, 0.0],
            ])
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 2.0, 1.0e5, 1.0e5, 1.0, 1.0],
                [-13.2, 1, 0, 2, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 1.0, 1.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
0 1 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.2 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 -0.1 0
SCALARS Electrons float 1
LOOKUP_TABLE default
4 2 4
SCALARS Holes float 1
LOOKUP_TABLE default
3 2 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.25 0.25 0.25
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.12 0.12 0.12
SCALARS ElectronLowFieldMobility float 1
LOOKUP_TABLE default
0.5 0.5 0.5
SCALARS HoleLowFieldMobility float 1
LOOKUP_TABLE default
0.4 0.4 0.4
SCALARS ElectronMobilityLimiter float 1
LOOKUP_TABLE default
0.5 0.5 0.5
SCALARS HoleMobilityLimiter float 1
LOOKUP_TABLE default
0.3 0.3 0.3
SCALARS ElectricField float 1
LOOKUP_TABLE default
10 20 30
SCALARS ElectronHighFieldDrive float 1
LOOKUP_TABLE default
100 200 300
SCALARS HoleHighFieldDrive float 1
LOOKUP_TABLE default
400 500 600
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.05, 0.0],
                "eQuasiFermiPotential": [0.0, 0.10, 0.0],
                "hQuasiFermiPotential": [0.0, -0.05, 0.0],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [3.0e-6, 2.0e-6, 3.0e-6],
                "eMobility": [2500.0, 2500.0, 2500.0],
                "hMobility": [600.0, 600.0, 600.0],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_edge_mobility_inputs.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(edge_csv),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vtk_root),
                    "--bias", "-13.2",
                    "--out-dir", str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "active_edge_mobility_inputs_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["active_edge_ids"], "0")
            self.assertAlmostEqual(float(row["electron_limited_mobility_product_m2_V_s"]), 0.25)
            self.assertAlmostEqual(float(row["electron_product_over_final_mobility"]), 1.0)
            self.assertAlmostEqual(float(row["electron_final_over_sentaurus_mobility"]), 1.0)
            self.assertAlmostEqual(float(row["hole_limited_mobility_product_m2_V_s"]), 0.12)
            self.assertAlmostEqual(float(row["hole_final_over_sentaurus_mobility"]), 2.0)
            self.assertAlmostEqual(float(row["vela_electric_field_scalar_avg_V_cm"]), 15.0)
            self.assertAlmostEqual(float(row["vela_electric_field_scalar_avg_V_m"]), 1500.0)
            self.assertAlmostEqual(float(row["vela_electron_high_field_drive_avg_V_m"]), 15000.0)
            summary = json.loads((out_dir / "active_edge_mobility_inputs_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_active_edge_density_factors_split_ni_and_exponent(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_active_edge_density_factors_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            vt = 8.617333262145e-5 * 300.0
            sentaurus_hole_qf = vt * math.log(2.0)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                ],
            }), encoding="utf-8")
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [-13.2, 2, 0, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
0 1 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
2 1 2
SCALARS Holes float 1
LOOKUP_TABLE default
3 3 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [sentaurus_hole_qf, sentaurus_hole_qf, sentaurus_hole_qf],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [6.0e-6, 6.0e-6, 6.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_edge_density_factors.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--mesh",
                    str(mesh),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--vela-vtk-root",
                    str(vtk_root),
                    "--bias",
                    "-13.2",
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "active_edge_density_factors_nodes.csv")
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertAlmostEqual(float(row["electron_density_geommean_vela_over_sentaurus"]), 0.5)
            self.assertAlmostEqual(float(row["electron_inferred_ni_geommean_vela_over_sentaurus"]), 0.5)
            self.assertAlmostEqual(float(row["electron_boltzmann_factor_geommean_vela_over_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(row["hole_density_geommean_vela_over_sentaurus"]), 0.5)
            self.assertAlmostEqual(float(row["hole_inferred_ni_geommean_vela_over_sentaurus"]), 1.0)
            self.assertAlmostEqual(float(row["hole_boltzmann_factor_geommean_vela_over_sentaurus"]), 0.5)
            summary = json.loads((out_dir / "active_edge_density_factors_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_active_edge_mixed_state_replay_applies_qf_shift(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_active_edge_mixed_state_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            vtk_root = root / "vtk"
            fields.mkdir(parents=True)
            vtk_root.mkdir()
            mesh = root / "mesh.json"
            support = root / "support.csv"
            edge_csv = root / "sg_edges.csv"
            out_dir = root / "out"
            vt = 8.617333262145e-5 * 300.0
            hole_shift = vt * math.log(2.0)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2]},
                ],
            }), encoding="utf-8")
            self._write_csv(support, [
                "node_id",
                "x_um",
                "y_um",
                "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s",
                "sentaurus_active",
                "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0, 1.0, 1, 1, "overlap"]])
            self._write_csv(edge_csv, [
                "bias_V",
                "edge_id",
                "node0",
                "node1",
                "x0_um",
                "y0_um",
                "x1_um",
                "y1_um",
                "edge_area_proxy_m2",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [-13.2, 2, 0, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
            ])
            (vtk_root / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
0 1 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
4 2 4
SCALARS Holes float 1
LOOKUP_TABLE default
3 2 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
""".strip() + "\n",
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [hole_shift, hole_shift, hole_shift],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [6.0e-6, 4.0e-6, 6.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_edge_mixed_state_replay.py"),
                    "--support-csv",
                    str(support),
                    "--sg-edge-csv",
                    str(edge_csv),
                    "--mesh",
                    str(mesh),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--vela-vtk-root",
                    str(vtk_root),
                    "--bias",
                    "-13.2",
                    "--hole-qf-shift-v",
                    str(hole_shift),
                    "--out-dir",
                    str(out_dir),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out_dir / "active_edge_mixed_state_replay_nodes.csv")
            by_variant = {row["variant"]: row for row in rows}
            self.assertAlmostEqual(
                float(by_variant["vela_baseline"]["hole_density_geommean_over_sentaurus"]),
                0.5,
            )
            self.assertAlmostEqual(
                float(by_variant["vela_baseline"]["particle_flux_over_sentaurus"]),
                0.75,
            )
            self.assertAlmostEqual(
                float(by_variant["vela_qf_shift"]["hole_density_geommean_over_sentaurus"]),
                1.0,
            )
            self.assertAlmostEqual(
                float(by_variant["vela_qf_shift"]["particle_flux_over_sentaurus"]),
                1.0,
            )
            self.assertAlmostEqual(
                float(by_variant["vela_qf_shift_sentaurus_mobility"]["particle_flux_over_sentaurus"]),
                1.0,
            )
            summary = json.loads((out_dir / "active_edge_mixed_state_replay_summary.json").read_text())
            self.assertEqual(summary["support_classes"]["overlap"]["variants"]["vela_qf_shift"]["count"], 1)

    def test_pn2d_bv_field_source_summary_reports_peak_and_sum_proxy(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_field_source_summary_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus" / "sentaurus_-0.5v"
            fields = sentaurus / "fields"
            vela = root / "vela"
            out = root / "out"
            fields.mkdir(parents=True)
            vela.mkdir()
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 2.0, 0.0],
            ])
            self._write_csv(fields / "ElectricField_region0.csv", ["node_id", "component0"], [
                [0, 2.0],
                [1, 4.0],
                [2, 8.0],
            ])
            self._write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [
                [0, 1.0e3],
                [1, 2.0e3],
                [2, 4.0e3],
            ])
            (vela / "mini_0000_-0.5V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1 0 0
2 0 0
POINT_DATA 3
SCALARS ElectricField float 1
LOOKUP_TABLE default
1
2
4
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
500
1000
2000
""".strip(),
                encoding="utf-8",
            )
            field_compare = root / "field_compare.csv"
            self._write_csv(field_compare, [
                "bias_V",
                "quantity",
                "avalanche_status",
            ], [
                [-0.5, "avalanche_generation", "thresholded_peak"],
                [-0.5, "electric_field", ""],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_field_source_summary.py"),
                    "--sentaurus-root",
                    str(root / "sentaurus"),
                    "--vela-vtk-root",
                    str(vela),
                    "--field-compare",
                    str(field_compare),
                    "--out-dir",
                    str(out),
                    "--biases",
                    "-0.5",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            rows = self._read_csv(out / "field_source_summary.csv")
            by_quantity = {row["quantity"]: row for row in rows}
            self.assertAlmostEqual(float(by_quantity["electric_field"]["sentaurus_max"]), 8.0)
            self.assertAlmostEqual(float(by_quantity["electric_field"]["vela_max"]), 4.0)
            self.assertAlmostEqual(float(by_quantity["electric_field"]["log10_vela_over_sentaurus_sum_proxy"]), math.log10(0.5))
            self.assertEqual(by_quantity["avalanche_generation"]["field_compare_status"], "thresholded_peak")
            self.assertAlmostEqual(float(by_quantity["avalanche_generation"]["sentaurus_sum_proxy"]), 7000.0)
            self.assertAlmostEqual(float(by_quantity["avalanche_generation"]["vela_sum_proxy"]), 3500.0e-6)
            self.assertAlmostEqual(
                float(by_quantity["avalanche_generation"]["log10_vela_over_sentaurus_sum_proxy"]),
                math.log10(0.5e-6),
            )

    def test_pn2d_bv_transport_driver_summary_aggregates_density_mobility_and_flux(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_transport_driver_summary_") as tmp:
            root = Path(tmp)
            contact_edges = root / "contact_edges_filtered.csv"
            feedback = root / "continuity_feedback_edges.csv"
            out = root / "out"
            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "vela_electron_density_avg_cm3",
                "sentaurus_electron_density_cm3",
                "vela_hole_density_avg_cm3",
                "sentaurus_hole_density_cm3",
                "vela_electron_mobility_cm2_V_s",
                "sentaurus_electron_mobility_cm2_V_s",
                "vela_hole_mobility_cm2_V_s",
                "sentaurus_hole_mobility_cm2_V_s",
            ], [
                [-0.5, "Anode", 1, 50.0, 100.0, 25.0, 100.0, 80.0, 100.0, 50.0, 100.0],
                [-0.5, "Anode", 2, 25.0, 100.0, 50.0, 100.0, 40.0, 100.0, 25.0, 100.0],
            ])
            self._write_csv(feedback, [
                "bias_V",
                "relation",
                "vela_electron_flux_abs_m2_s",
                "sentaurus_electron_flux_abs_m2_s",
                "vela_hole_flux_abs_m2_s",
                "sentaurus_hole_flux_abs_m2_s",
            ], [[-0.5, "focus", 30.0, 100.0, 10.0, 100.0]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_transport_driver_summary.py"),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--continuity-feedback-csv",
                    str(feedback),
                    "--out-dir",
                    str(out),
                    "--biases",
                    "-0.5",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            rows = self._read_csv(out / "transport_driver_summary.csv")
            by_key = {(row["driver"], row["carrier"]): row for row in rows}
            self.assertAlmostEqual(
                float(by_key[("contact_edge_density", "electron")]["log10_p50_vela_over_sentaurus"]),
                0.5 * (math.log10(0.5) + math.log10(0.25)),
            )
            self.assertAlmostEqual(
                float(by_key[("contact_edge_mobility", "hole")]["log10_p50_vela_over_sentaurus"]),
                0.5 * (math.log10(0.5) + math.log10(0.25)),
            )
            self.assertAlmostEqual(
                float(by_key[("interior_sg_flux", "electron")]["log10_p50_vela_over_sentaurus"]),
                math.log10(0.3),
            )
            self.assertEqual(by_key[("interior_sg_flux", "hole")]["classification"], "driver_mismatch")

    def test_pn2d_bv_contact_bgn_density_diagnostic_separates_ni_and_qf_drive(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_contact_bgn_density_") as tmp:
            root = Path(tmp)
            contact_edges = root / "contact_edges_filtered.csv"
            out = root / "out"
            self._write_csv(contact_edges, [
                "bias_V",
                "current_contact",
                "edge_id",
                "vela_electron_density_avg_cm3",
                "sentaurus_electron_density_cm3",
                "vela_hole_density_avg_cm3",
                "sentaurus_hole_density_cm3",
                "vela_potential_avg_V",
                "sentaurus_potential_V",
                "vela_electron_qf_avg_V",
                "sentaurus_electron_qf_V",
                "vela_hole_qf_avg_V",
                "sentaurus_hole_qf_V",
            ], [[
                -0.5,
                "Anode",
                1,
                1.0e3,
                2.5e3,
                1.0e17,
                1.0e17,
                -0.3,
                -0.25,
                0.1,
                0.15,
                -0.2,
                -0.15,
            ]])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_contact_bgn_density.py"),
                    "--contact-edge-csv",
                    str(contact_edges),
                    "--out-dir",
                    str(out),
                    "--biases",
                    "-0.5",
                    "--contact",
                    "Anode",
                    "--thermal-voltage",
                    "0.1",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            summary_rows = self._read_csv(out / "contact_bgn_density_summary.csv")
            by_carrier = {row["carrier"]: row for row in summary_rows}
            self.assertAlmostEqual(
                float(by_carrier["electron"]["log10_density_p50_vela_over_sentaurus"]),
                math.log10(0.4),
            )
            self.assertAlmostEqual(
                float(by_carrier["electron"]["log10_inferred_ni_eff_p50_vela_over_sentaurus"]),
                math.log10(0.4),
            )
            self.assertAlmostEqual(float(by_carrier["electron"]["drive_delta_V_p50"]), 0.0)
            self.assertEqual(by_carrier["electron"]["classification"], "ni_eff_mismatch")
            self.assertEqual(by_carrier["hole"]["classification"], "density_close")

            edge_rows = self._read_csv(out / "contact_bgn_density_edges.csv")
            self.assertAlmostEqual(
                float(edge_rows[0]["vela_electron_inferred_ni_eff_cm3"]),
                1.0e3 / math.exp((-0.3 - 0.1) / 0.1),
            )
            report = json.loads((out / "summary.json").read_text())
            self.assertEqual(report["contact"], "Anode")

    def test_pn2d_solver_sensitivity_matrix_writes_configs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_solver_sensitivity_") as tmp:
            root = Path(tmp)
            out = root / "out"
            base_config = root / "base_config.json"
            base_config.write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "materials_file": "materials.json",
                "node_doping_file": "doping.csv",
                "output_csv": "original.csv",
                "solver": {
                    "method": "gummel_newton",
                    "preserved_knob": 123,
                },
                "sweep": {
                    "mode": "original",
                    "start": 1.0,
                    "stop": 2.0,
                    "step": 0.25,
                    "write_vtk": False,
                },
            }), encoding="utf-8")
            cmd = [
                sys.executable,
                "scripts/run_pn2d_solver_sensitivity_matrix.py",
                "--base-config",
                str(base_config),
                "--out-dir",
                str(out),
                "--dry-run",
                "--bias-window",
                "reverse_low",
            ]
            result = subprocess.run(cmd, cwd=REPO, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads((out / "manifest.json").read_text())
            self.assertEqual(manifest["base_config"], str(base_config.resolve()))
            cases_by_name = {case["name"]: case for case in manifest["cases"]}
            names = set(cases_by_name)
            self.assertIn("baseline", names)
            self.assertIn("bank_rose_like_damped", names)
            self.assertIn("qf_hard_limit_0p0259", names)
            self.assertIn("bgn_contact_ni_match_offset_0p0197", names)
            self.assertEqual(cases_by_name["baseline"]["status"], "generated")
            self.assertEqual(cases_by_name["qf_hard_limit_0p0259"]["status"], "generated")
            self.assertEqual(
                cases_by_name["qf_hard_limit_0p0259"]["reason"],
                "dry-run config generation only",
            )
            self.assertEqual(cases_by_name["baseline"]["config"], str((out / "baseline" / "simulation.json").resolve()))
            self.assertEqual(cases_by_name["baseline"]["out_dir"], str((out / "baseline").resolve()))
            baseline = json.loads((out / "baseline" / "simulation.json").read_text())
            self.assertEqual(baseline["simulation_type"], "dc_sweep")
            self.assertEqual(baseline["mesh_file"], str((root / "mesh.json").resolve()))
            self.assertEqual(baseline["materials_file"], str((root / "materials.json").resolve()))
            self.assertEqual(baseline["node_doping_file"], str((root / "doping.csv").resolve()))
            self.assertEqual(baseline["solver"]["preserved_knob"], 123)
            self.assertEqual(baseline["sweep"]["mode"], "original")
            self.assertEqual(baseline["sweep"]["start"], 0)
            self.assertEqual(baseline["sweep"]["stop"], -5)
            self.assertEqual(baseline["sweep"]["step"], -0.05)
            self.assertIs(baseline["sweep"]["write_vtk"], True)
            self.assertEqual(baseline["output_csv"], str((out / "baseline" / "iv.csv").resolve()))
            self.assertEqual(
                baseline["sweep"]["vtk_prefix"],
                str((out / "baseline" / "vtk" / "baseline").resolve()),
            )
            branch_guard = json.loads((out / "branch_guard_0p05" / "simulation.json").read_text())
            self.assertNotIn("continuation", branch_guard)
            self.assertTrue(branch_guard["sweep"]["continuation"]["branch_acceptance"]["psi_phin_jump"])
            self.assertEqual(
                branch_guard["sweep"]["continuation"]["branch_acceptance"]["max_psi_phin_jump_V"],
                0.05,
            )
            qf_limited = json.loads((out / "qf_hard_limit_0p0259" / "simulation.json").read_text())
            self.assertEqual(qf_limited["solver"]["quasi_fermi_update_limit_V"], 0.0259)
            bgn_contact = json.loads((out / "bgn_contact_ni_match_offset_0p0197" / "simulation.json").read_text())
            self.assertEqual(bgn_contact["solver"]["bandgap_narrowing"]["model"], "old_slotboom")
            self.assertAlmostEqual(
                bgn_contact["solver"]["bandgap_narrowing"]["offset_eV"],
                0.0197,
            )
            split = json.loads((out / "sentaurus_split_ni_slotboom" / "simulation.json").read_text())
            self.assertEqual(split["solver"]["bandgap_narrowing"], "slotboom")
            self.assertEqual(
                split["materials_file"],
                str(
                    REPO
                    / "reference_tcad"
                    / "pn2d_sentaurus2018"
                    / "source"
                    / "pn2d_sentaurus2018_iv_materials.json"
                ),
            )

    def test_pn2d_bv_mobility_debug_script_writes_decomposition_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_mobility_debug_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus" / "sentaurus_0v"
            fields = sentaurus / "fields"
            vela = root / "vela"
            out = root / "out"
            fields.mkdir(parents=True)
            vela.mkdir()
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            self._write_csv(sentaurus / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                [0, 0, 1, 2, "R.Si", "Si"],
            ])
            for name, values in {
                "eMobility": [700.0, 300.0, 700.0],
                "hMobility": [320.0, 180.0, 320.0],
                "eDensity": [1.0e3, 1.0e8, 1.0e3],
                "hDensity": [1.0e8, 1.0e3, 1.0e8],
                "ElectrostaticPotential": [0.0, 0.1, 0.0],
                "ElectricField": [0.0, 1.0e3, 0.0],
                "eQuasiFermiPotential": [0.0, 0.2, 0.0],
                "hQuasiFermiPotential": [0.0, -0.1, 0.0],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])
            (vela / "mini_0000_0V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.07 0.06 0.07
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.032 0.03 0.032
SCALARS ElectronLowFieldMobility float 1
LOOKUP_TABLE default
0.07 0.07 0.07
SCALARS HoleLowFieldMobility float 1
LOOKUP_TABLE default
0.032 0.032 0.032
SCALARS ElectronHighFieldDrive float 1
LOOKUP_TABLE default
0 10000 0
SCALARS HoleHighFieldDrive float 1
LOOKUP_TABLE default
0 5000 0
SCALARS ElectronMobilityLimiter float 1
LOOKUP_TABLE default
1 0.857142857 1
SCALARS HoleMobilityLimiter float 1
LOOKUP_TABLE default
1 0.9375 1
SCALARS Electrons float 1
LOOKUP_TABLE default
1e9 1e15 1e9
SCALARS Holes float 1
LOOKUP_TABLE default
1e15 1e9 1e15
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectricField float 1
LOOKUP_TABLE default
0 1000 0
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_mobility.py"),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "0",
                    "--electron-saturation-velocity-m-s", "1.07e5",
                    "--hole-saturation-velocity-m-s", "8.37e4",
                    "--electron-field-beta", "1.109",
                    "--hole-field-beta", "1.213",
                ],
                check=True,
                cwd=REPO,
            )

            summary_rows = self._read_csv(out / "mobility_decomposition_summary.csv")
            self.assertEqual(summary_rows[0]["bias_V"], "0.0")
            self.assertIn("sentaurus_eMobility_min", summary_rows[0])
            self.assertIn("vela_electron_limiter_min", summary_rows[0])
            self.assertIn("electron_density_to_mobility_top_error_distance_um", summary_rows[0])
            self.assertIn("electron_limiter_to_mobility_top_error_distance_um", summary_rows[0])
            self.assertIn("sentaurus_grad_eQuasiFermiPotential_V_per_cm_p95", summary_rows[0])
            self.assertIn("sentaurus_electric_field_V_per_cm_p95", summary_rows[0])
            self.assertIn("sentaurus_electron_implied_high_field_drive_V_per_cm_p95", summary_rows[0])
            self.assertIn("electron_implied_to_vela_drive_ratio_p95", summary_rows[0])
            self.assertIn("electron_implied_drive_p95_over_vela_drive_p95", summary_rows[0])
            self.assertEqual(summary_rows[0]["mobility_model_electron_field_beta"], "1.109")
            self.assertEqual(summary_rows[0]["mobility_model_hole_saturation_velocity_m_s"], "83700.0")
            self.assertGreater(
                float(summary_rows[0]["sentaurus_electron_implied_high_field_drive_V_per_cm_p95"]),
                0.0,
            )
            neighborhood_rows = self._read_csv(out / "top_error_neighborhood.csv")
            self.assertTrue(neighborhood_rows)
            anchors = {row["anchor"] for row in neighborhood_rows}
            self.assertIn("electron_mobility_top_error", anchors)
            self.assertIn("electron_density_top_error", anchors)
            electron_mobility = next(
                row for row in neighborhood_rows
                if row["anchor"] == "electron_mobility_top_error"
            )
            self.assertIn("sentaurus_eMobility", electron_mobility)
            self.assertIn("vela_electron_mobility", electron_mobility)
            self.assertIn("vela_electron_high_field_drive_V_per_cm", electron_mobility)
            self.assertIn("sentaurus_electron_implied_high_field_drive_V_per_cm", electron_mobility)
            self.assertIn("electron_density_log10_error", electron_mobility)
            self.assertIn("nearest_vela_node_distance_um", electron_mobility)
            self.assertIn("sentaurus_eQuasiFermiPotential", electron_mobility)
            self.assertIn("vela_electron_quasi_fermi", electron_mobility)
            self.assertIn("sentaurus_grad_eQuasiFermiPotential_V_per_cm", electron_mobility)
            self.assertIn("vela_grad_phin_V_per_cm", electron_mobility)
            profile_rows = self._read_csv(out / "quasi_fermi_anchor_profiles.csv")
            self.assertTrue(profile_rows)
            self.assertIn("profile_axis", profile_rows[0])
            self.assertIn("sentaurus_eQuasiFermiPotential", profile_rows[0])
            self.assertIn("vela_electron_quasi_fermi", profile_rows[0])

    def test_pn2d_bv_qf_anchor_script_reports_internal_phin_shift(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_qf_anchor_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus" / "sentaurus_0v"
            fields = sentaurus / "fields"
            vela = root / "vela"
            out = root / "out"
            fields.mkdir(parents=True)
            vela.mkdir()
            mesh = {
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 2.0, "y": 0.0},
                ],
                "triangles": [],
                "contacts": [
                    {"name": "Anode", "node_ids": [0]},
                    {"name": "Cathode", "node_ids": [2]},
                ],
            }
            (root / "mesh.json").write_text(json.dumps(mesh), encoding="utf-8")
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 2.0, 0.0],
            ])
            self._write_csv(sentaurus / "contacts.csv", ["name", "node_ids", "region"], [
                ["Anode", "0", "R.Si"],
                ["Cathode", "2", "R.Si"],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.0, 0.0],
                "eQuasiFermiPotential": [0.0, 0.0, 0.0],
                "hQuasiFermiPotential": [0.0, 0.0, 0.0],
                "eDensity": [1.0e3, 1.0e3, 1.0e3],
                "hDensity": [1.0e3, 1.0e3, 1.0e3],
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, values[0]],
                    [1, values[1]],
                    [2, values[2]],
                ])
            ratio = math.exp(0.18 / (8.617333262145e-5 * 300.0))
            (vela / "mini_0000_0V.vtk").write_text(
                f"""
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
2e-6 0 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.18 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 0 0
SCALARS Electrons float 1
LOOKUP_TABLE default
1e9 {1.0e9 * ratio} 1e9
SCALARS Holes float 1
LOOKUP_TABLE default
1e9 1e9 1e9
""".lstrip(),
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_qf_anchor.py"),
                    "--mesh", str(root / "mesh.json"),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "0",
                    "--bands", "middle:0.5:1.5",
                ],
                cwd=REPO,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            band = self._read_csv(out / "qf_anchor_band_summary.csv")[0]
            self.assertEqual(band["band"], "middle")
            self.assertAlmostEqual(float(band["delta_psi_minus_phin_median_V"]), 0.18)
            self.assertAlmostEqual(float(band["uniform_phin_shift_to_match_electron_median_V"]), 0.18)
            self.assertGreater(float(band["contact_phin_violation_if_uniform_shift_V"]), 0.17)
            self.assertAlmostEqual(
                float(band["electron_log10_ratio_median"]),
                float(band["electron_log10_ratio_from_exponent_median"]),
                places=6,
            )

            contacts = self._read_csv(out / "qf_anchor_contact_summary.csv")
            self.assertEqual({row["contact"] for row in contacts}, {"Anode", "Cathode"})
            self.assertTrue(all(abs(float(row["delta_phin_median_V"])) < 1.0e-12 for row in contacts))

    def test_pn2d_bv_heatmap_accepts_sentaurus_impact_ionization_field(self) -> None:
        text = (REPO / "scripts" / "plot_pn2d_multibias_heatmaps.py").read_text()

        self.assertIn('"sentaurus_field": ["AvalancheGeneration", "ImpactIonization"]', text)
        self.assertIn('"resolved_sentaurus_fields"', text)

    def test_pn2d_bv_field_compare_writes_validation_summary_artifacts(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_summary_compare_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus" / "sentaurus_0v"
            fields = sentaurus / "fields"
            vela = root / "vela"
            out = root / "reports"
            fields.mkdir(parents=True)
            vela.mkdir()

            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
                [3, 1.0, 1.0],
            ])
            self._write_csv(sentaurus / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                [0, 0, 1, 2, "R.Si", "Si"],
                [1, 1, 3, 2, "R.Si", "Si"],
            ])
            self._write_csv(sentaurus / "contacts.csv", ["name", "node_ids", "region"], [
                ["Anode", "0;2", "R.Si"],
                ["Cathode", "1;3", "R.Si"],
            ])
            self._write_csv(sentaurus / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 1.0e17, 0.0],
            ])
            self._write_csv(fields / "ElectrostaticPotential_region0.csv", ["node_id", "component0"], [
                [0, -0.4],
                [1, 0.4],
                [2, -0.4],
                [3, 0.4],
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [0, 1.0e3],
                [1, 1.0e10],
                [2, 1.0e3],
                [3, 1.0e10],
            ])
            self._write_csv(fields / "ImpactIonization_region0.csv", ["node_id", "component0"], [
                [0, 0.0],
                [1, 1.0e2],
                [2, 0.0],
                [3, 1.0e4],
            ])
            (vela / "mini_0000_0V.vtk").write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
-0.39 0.41 -0.42 0.38
SCALARS Electrons float 1
LOOKUP_TABLE default
1e12 1e16 1e12 1e16
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
0 1e8 0 1e10
""".lstrip()
            )
            reference_curve = root / "sentaurus_curve.csv"
            candidate_curve = root / "vela_curve.csv"
            self._write_csv(reference_curve, ["bias_V", "current_total"], [[0.0, 1.0e-26]])
            self._write_csv(candidate_curve, [
                "bias_V", "current_total_A_per_um", "converged",
            ], [[0.0, 1.0e-20, 1]])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_bv_multibias_fields.py"),
                    "--sentaurus-root",
                    str(root / "sentaurus"),
                    "--vela-vtk-root",
                    str(vela),
                    "--curve-reference",
                    str(reference_curve),
                    "--curve-candidate",
                    str(candidate_curve),
                    "--out-dir",
                    str(out),
                    "--biases",
                    "0",
                    "--quantities",
                    "potential,electron_density,avalanche_generation",
                ],
                check=True,
                cwd=REPO,
            )

            self.assertTrue((out / "field_compare.csv").is_file())
            self.assertTrue((out / "curve_compare.csv").is_file())
            self.assertTrue((out / "debug_ranking.json").is_file())
            self.assertTrue((out / "README.md").is_file())

            curve_rows = self._read_csv(out / "curve_compare.csv")
            self.assertEqual(curve_rows[0]["status"], "below_current_floor")
            self.assertEqual(curve_rows[0]["log10_abs_ratio"], "")

            field_rows = self._read_csv(out / "field_compare.csv")
            electron = next(row for row in field_rows if row["quantity"] == "electron_density")
            self.assertEqual(electron["status"], "ok")
            self.assertIn("top_error_x_um", electron)
            self.assertIn("p_contact_error", electron)
            self.assertIn("n_contact_error", electron)
            self.assertIn("junction_error", electron)
            self.assertIn("bulk_error", electron)

            ranking = json.loads((out / "debug_ranking.json").read_text())
            self.assertEqual(ranking["field_rankings"][0]["quantity"], "electron_density")

    def test_pn2d_bv_signed_partition_diagnostic_reports_direction(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_signed_partition_") as tmp:
            root = Path(tmp)
            sentaurus = root / "sentaurus_0v"
            fields = sentaurus / "fields"
            out = root / "out"
            fields.mkdir(parents=True)

            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
                [3, 1.0, 1.0],
            ])
            self._write_csv(sentaurus / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                [0, 0, 1, 2, "R.Si", "Si"],
                [1, 1, 3, 2, "R.Si", "Si"],
            ])
            self._write_csv(sentaurus / "contacts.csv", ["name", "node_ids", "region"], [
                ["Anode", "0;2", "R.Si"],
                ["Cathode", "1;3", "R.Si"],
            ])
            self._write_csv(sentaurus / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 1.0e17, 0.0],
            ])
            self._write_csv(fields / "eDensity_region0.csv", ["node_id", "component0"], [
                [0, 100.0],
                [1, 100.0],
                [2, 100.0],
                [3, 100.0],
            ])
            self._write_csv(fields / "hDensity_region0.csv", ["node_id", "component0"], [
                [0, 100.0],
                [1, 100.0],
                [2, 100.0],
                [3, 100.0],
            ])
            vtk = root / "vela_0000_0V.vtk"
            vtk.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS Electrons float 1
LOOKUP_TABLE default
5e7 5e7 5e7 5e7
SCALARS Holes float 1
LOOKUP_TABLE default
1e8 1e7 1e8 1e7
""".lstrip()
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_signed_field_partitions.py"),
                    "--sentaurus-dir",
                    str(sentaurus),
                    "--vela-vtk",
                    str(vtk),
                    "--out-dir",
                    str(out),
                    "--bias",
                    "0",
                    "--quantities",
                    "electron_density,hole_density",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "signed_field_partition_summary.csv")
            by_key = {(row["quantity"], row["region"]): row for row in rows}
            electron_all = by_key[("electron_density", "all")]
            self.assertAlmostEqual(
                float(electron_all["signed_log10_median_vela_over_sentaurus"]),
                math.log10(0.5),
            )
            self.assertEqual(electron_all["direction"], "vela_lower")
            hole_n = by_key[("hole_density", "n_contact")]
            self.assertAlmostEqual(
                float(hole_n["signed_log10_median_vela_over_sentaurus"]),
                -1.0,
            )
            self.assertEqual(hole_n["direction"], "vela_lower")

            top = self._read_csv(out / "signed_field_partition_top_errors.csv")
            self.assertTrue(any(row["quantity"] == "hole_density" for row in top))
            report = json.loads((out / "signed_field_partition_summary.json").read_text())
            self.assertEqual(report["bias_V"], 0.0)
            self.assertEqual(report["summary_row_count"], len(rows))

    def test_pn2d_bv_branch_transition_jump_report_finds_density_step(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_bv_branch_transition_") as tmp:
            root = Path(tmp)
            signed_root = root / "signed"
            s0 = root / "sentaurus_0v"
            s1 = root / "sentaurus_m0p1v"
            vtk0 = root / "state_0V.vtk"
            vtk1 = root / "state_-0.1V.vtk"
            out = root / "out"
            signed_root.mkdir()

            for sentaurus, current in [(s0, -1.0e-17), (s1, -2.0e-17)]:
                fields = sentaurus / "fields"
                fields.mkdir(parents=True)
                self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                    [0, 0.0, 0.0],
                    [1, 1.0, 0.0],
                    [2, 0.0, 1.0],
                    [3, 1.0, 1.0],
                ])
                self._write_csv(sentaurus / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                    [0, 0, 1, 2, "R.Si", "Si"],
                    [1, 1, 3, 2, "R.Si", "Si"],
                ])
                self._write_csv(sentaurus / "contacts.csv", ["name", "node_ids", "region"], [
                    ["Anode", "0;2", "R.Si"],
                    ["Cathode", "1;3", "R.Si"],
                ])
                self._write_csv(sentaurus / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                    [0, 0.0, 1.0e17],
                    [1, 1.0e17, 0.0],
                    [2, 0.0, 1.0e17],
                    [3, 1.0e17, 0.0],
                ])
                self._write_csv(fields / "ContactCurrentFlux_region2.csv", ["node_id", "component0"], [
                    [2, current],
                ])

            for label, bias, sentaurus, vtk, electron_median in [
                ("b0", 0.0, s0, vtk0, 0.0),
                ("b1", -0.1, s1, vtk1, 0.25),
            ]:
                signed = signed_root / label
                signed.mkdir()
                self._write_csv(signed / "signed_field_partition_summary.csv", [
                    "bias_V",
                    "quantity",
                    "region",
                    "status",
                    "signed_log10_median_vela_over_sentaurus",
                    "signed_log10_p95_abs",
                    "signed_delta_median",
                    "signed_delta_p95_abs",
                    "direction",
                ], [
                    [bias, "electron_density", "all", "ok", electron_median, 0.5, "", "", "vela_higher"],
                    [bias, "electron_density", "junction", "ok", electron_median + 0.1, 0.6, "", "", "vela_higher"],
                    [bias, "hole_density", "all", "ok", -0.1, 0.2, "", "", "vela_lower"],
                    [bias, "hole_density", "junction", "ok", -0.2, 0.3, "", "", "vela_lower"],
                    [bias, "potential", "all", "ok", "", "", 0.01, 0.02, "vela_higher"],
                ])
                (signed / "signed_field_partition_summary.json").write_text(
                    json.dumps({
                        "bias_V": bias,
                        "sentaurus_dir": str(sentaurus),
                        "vela_vtk": str(vtk),
                    }),
                    encoding="utf-8",
                )

            vtk0.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0 0.1
SCALARS ElectricField float 1
LOOKUP_TABLE default
10 20 30 40
SCALARS Electrons float 1
LOOKUP_TABLE default
1e8 1e8 1e8 1e8
SCALARS Holes float 1
LOOKUP_TABLE default
1e8 1e8 1e8 1e8
""".lstrip(),
                encoding="utf-8",
            )
            vtk1.write_text(
                """
# vtk DataFile Version 2.0
mini
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0 0.2
SCALARS ElectricField float 1
LOOKUP_TABLE default
10 20 35 50
SCALARS Electrons float 1
LOOKUP_TABLE default
1e8 1e9 1e10 1e8
SCALARS Holes float 1
LOOKUP_TABLE default
1e8 1e8 1e8 1e9
""".lstrip(),
                encoding="utf-8",
            )
            candidate_curve = root / "vela_iv.csv"
            self._write_csv(candidate_curve, ["bias_V", "current_total_A_per_um", "converged"], [
                [0.0, -1.0e-17, 1],
                [-0.1, -4.0e-17, 1],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_branch_transition_jumps.py"),
                    "--signed-root",
                    str(signed_root),
                    "--out-dir",
                    str(out),
                    "--points",
                    "0:b0,-0.1:b1",
                    "--curve-candidate",
                    str(candidate_curve),
                    "--quantities",
                    "electron_density,hole_density",
                    "--regions",
                    "all,junction",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = self._read_csv(out / "branch_transition_summary.csv")
            electron_b1 = next(
                row for row in summary
                if row["bias_V"] == "-0.1" and row["quantity"] == "electron_density" and row["region"] == "all"
            )
            self.assertAlmostEqual(float(electron_b1["signed_log10_median_vela_over_sentaurus"]), 0.25)
            self.assertAlmostEqual(float(electron_b1["current_log10_abs_ratio"]), math.log10(2.0))

            jumps = self._read_csv(out / "branch_transition_adjacent_jumps.csv")
            electron_jump = next(row for row in jumps if row["quantity"] == "electron_density" and row["region"] == "all")
            self.assertEqual(int(electron_jump["top_node_id"]), 2)
            self.assertAlmostEqual(float(electron_jump["top_signed_log10_jump"]), 2.0)
            hole_jump = next(row for row in jumps if row["quantity"] == "hole_density" and row["region"] == "all")
            self.assertEqual(int(hole_jump["top_node_id"]), 3)
            self.assertAlmostEqual(float(hole_jump["top_signed_log10_jump"]), 1.0)

            report = json.loads((out / "branch_transition_summary.json").read_text())
            self.assertEqual(report["point_count"], 2)
            self.assertEqual(report["adjacent_jump_row_count"], len(jumps))

    def test_pn2d_bv_avalanche_hotspot_diagnostic_reports_geometry(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_avalanche_hotspots_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            vtk = root / "state.vtk"
            out = root / "reports"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.2, "y": 0.1},
                    {"id": 3, "x": -0.1, "y": 0.4},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                    {"id": 1, "region_id": 0, "node_ids": [0, 2, 3]},
                ],
                "regions": [],
                "contacts": [],
            }))
            vtk.write_text(
                """
# vtk DataFile Version 2.0
hotspots
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
2e-7 1e-7 0
-1e-7 4e-7 0
CELLS 2 8
3 0 1 2
3 0 2 3
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
10 100 1000 0
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20 4e20
SCALARS Holes float 1
LOOKUP_TABLE default
4e20 3e20 2e20 1e20
SCALARS ElectricField float 1
LOOKUP_TABLE default
1e5 2e5 3e5 4e5
SCALARS ElectronHighFieldDrive float 1
LOOKUP_TABLE default
1e6 2e6 3e6 4e6
SCALARS HoleHighFieldDrive float 1
LOOKUP_TABLE default
4e6 3e6 2e6 1e6
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_avalanche_hotspots.py"),
                    "--vtk",
                    str(vtk),
                    "--mesh",
                    str(mesh),
                    "--out-dir",
                    str(out),
                    "--top",
                    "4",
                ],
                check=True,
                cwd=REPO,
            )

            rows = self._read_csv(out / "avalanche_hotspots.csv")
            self.assertEqual(rows[0]["node_id"], "2")
            self.assertAlmostEqual(float(rows[0]["avalanche_generation_m3_s"]), 1000.0)
            self.assertEqual(rows[0]["adjacent_element_count"], "2")
            self.assertEqual(rows[0]["obtuse_adjacent_element_count"], "1")
            self.assertGreater(float(rows[0]["node_control_volume_m2"]), 0.0)
            self.assertAlmostEqual(float(rows[0]["electric_field_V_m"]), 3.0e7)
            self.assertAlmostEqual(float(rows[0]["electron_high_field_drive_V_m"]), 3.0e8)
            self.assertAlmostEqual(float(rows[0]["hole_high_field_drive_V_m"]), 2.0e8)

            summary = json.loads((out / "avalanche_hotspots_summary.json").read_text())
            self.assertEqual(summary["node_count"], 4)
            self.assertEqual(summary["top_node"]["node_id"], 2)
            self.assertAlmostEqual(summary["concentration"]["top1_fraction"], 1000.0 / 1110.0)
            self.assertEqual(summary["geometry"]["top1_obtuse_adjacent_element_count"], 1)

    def test_pn2d_bv_sg_avalanche_edge_diagnostic_decomposes_sources(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_sg_avalanche_edges_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            vtk = root / "state.vtk"
            out = root / "reports"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.2, "y": 1.0},
                    {"id": 3, "x": 1.2, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                    {"id": 1, "region_id": 0, "node_ids": [1, 3, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "region_id": 0, "node_ids": [0, 1]},
                ],
            }))
            vtk.write_text(
                """
# vtk DataFile Version 2.0
sg avalanche edges
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
2e-7 1e-6 0
1.2e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 8 1 4
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 1 8 4
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20 4e20
SCALARS Holes float 1
LOOKUP_TABLE default
4e20 3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21 1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e30 2e30 3e30 4e30
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_avalanche_edges.py"),
                    "--vtk",
                    str(vtk),
                    "--mesh",
                    str(mesh),
                    "--out-dir",
                    str(out),
                    "--bias",
                    "-13.2",
                    "--top",
                    "8",
                ],
                check=True,
                cwd=REPO,
            )

            rows = self._read_csv(out / "sg_avalanche_edges.csv")
            self.assertGreaterEqual(len(rows), 5)
            self.assertEqual(rows[0]["rank"], "1")
            self.assertGreater(float(rows[0]["source_integral"]), 0.0)
            classes = {row["edge_class"] for row in rows}
            self.assertIn("contact_boundary", classes)
            self.assertIn("boundary_noncontact", classes)

            summary = json.loads((out / "sg_avalanche_edge_summary.json").read_text())
            self.assertEqual(summary["bias_V"], -13.2)
            self.assertEqual(summary["edge_count"], 5)
            self.assertGreater(summary["total_source_integral_reconstructed"], 0.0)
            self.assertEqual(summary["edge_area_policy"], "0.5 * edge.length * edge.couple")
            by_class = {item["edge_class"]: item for item in summary["by_edge_class"]}
            self.assertIn("contact_boundary", by_class)
            self.assertGreater(by_class["contact_boundary"]["source_integral"], 0.0)
            by_node_class = {item["node_class"]: item for item in summary["by_node_class"]}
            self.assertIn("contact", by_node_class)
            self.assertGreater(by_node_class["contact"]["vtk_source_integral"], 0.0)

    def test_pn2d_bv_sg_edge_source_dump_compare_summarizes_cpp_and_python(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_sg_edge_compare_") as tmp:
            root = Path(tmp)
            cpp_csv = root / "cpp.csv"
            python_csv = root / "python.csv"
            out_json = root / "summary.json"

            cpp_csv.write_text(
                "\n".join([
                    "point_index,bias_V,edge_id,node0,node1,edge_source_integral,node0_source_integral,node1_source_integral,edge_class",
                    "7,-13.2,10,351,986,8.0,4.0,4.0,interior_bulk",
                    "7,-13.2,11,1,351,2.0,1.0,1.0,contact_edge",
                    "6,-10,99,1,2,100.0,50.0,50.0,interior_bulk",
                    "",
                ]),
                encoding="utf-8",
            )
            python_csv.write_text(
                "\n".join([
                    "rank,edge_id,node0,node1,source_integral,edge_class",
                    "1,10,351,986,7.5,interior_bulk",
                    "2,11,1,351,2.5,contact_boundary",
                    "",
                ]),
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "compare_pn2d_bv_sg_edge_source_dump.py"),
                    "--cpp-csv",
                    str(cpp_csv),
                    "--python-csv",
                    str(python_csv),
                    "--out-json",
                    str(out_json),
                    "--bias",
                    "-13.2",
                    "--node",
                    "351",
                    "--node",
                    "986",
                ],
                check=True,
                cwd=REPO,
            )

            summary = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(summary["bias_V"], -13.2)
            self.assertEqual(summary["cpp"]["top_edge_ids"], [10, 11])
            self.assertEqual(summary["python"]["top_edge_ids"], [10, 11])
            self.assertAlmostEqual(summary["cpp"]["total_source_integral"], 10.0)
            self.assertAlmostEqual(summary["python"]["total_source_integral"], 10.0)
            self.assertAlmostEqual(summary["cpp"]["interior_bulk_source_fraction"], 0.8)
            self.assertAlmostEqual(summary["python"]["contact_edge_source_fraction"], 0.25)
            self.assertAlmostEqual(summary["comparison"]["total_source_log10_cpp_over_python"], 0.0)
            self.assertAlmostEqual(summary["cpp"]["node_source_integrals"]["351"], 5.0)
            self.assertAlmostEqual(summary["python"]["node_source_integrals"]["986"], 3.75)

    def test_pn2d_bv_ionization_integral_diagnostic_writes_required_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_ionization_integral_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            vtk = root / "state_0000_-5V.vtk"
            out = root / "reports"

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }))
            vtk.write_text(
                """
# vtk DataFile Version 2.0
ionization integral
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS ElectricField float 1
LOOKUP_TABLE default
2e6 1e6 5e5
SCALARS ElectronHighFieldDrive float 1
LOOKUP_TABLE default
2e6 1e6 5e5
SCALARS HoleHighFieldDrive float 1
LOOKUP_TABLE default
1.5e6 7.5e5 5e5
""".lstrip(),
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_ionization_integral.py"),
                    "--vtk",
                    str(vtk),
                    "--mesh",
                    str(mesh),
                    "--out-dir",
                    str(out),
                    "--bias",
                    "-5",
                ],
                check=True,
                cwd=REPO,
            )

            rows = self._read_csv(out / "ionization_integral_paths.csv")
            self.assertEqual({row["carrier"] for row in rows}, {"electron", "hole"})
            for row in rows:
                self.assertEqual(row["bias_V"], "-5.0")
                self.assertGreater(float(row["max_field_V_per_m"]), 0.0)
                self.assertGreater(float(row["integral_alpha_dx"]), 0.0)
                self.assertIn("dominant_edge_or_cell_id", row)

            summary = json.loads((out / "ionization_integral_summary.json").read_text())
            self.assertEqual(summary["bias_V"], -5.0)
            self.assertEqual(summary["path_count"], 2)
            self.assertGreater(summary["max_integral_alpha_dx"], 0.1)

    def test_pn2d_bv_local_avalanche_factor_diagnostic_writes_edge_summary(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_local_avalanche_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_0v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.2, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }))
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            (vela / "mini_0000_0V.vtk").write_text(
                """
# vtk DataFile Version 2.0
local factors
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
2e-7 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 8 1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 1 8
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e30 2e30 3e30
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.5, 0.0],
            ])
            for name, value in {
                "ElectricField": 2.0e5,
                "ImpactIonization": 1.0e12,
                "eDensity": 2.0e14,
                "hDensity": 3.0e14,
                "eMobility": 1000.0,
                "hMobility": 500.0,
                "eCurrentDensity": 1.0e-6,
                "hCurrentDensity": 2.0e-6,
            }.items():
                self._write_csv(fields / f"{name}_region0.csv", ["node_id", "component0"], [
                    [0, value],
                ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_local_avalanche_factors.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "0",
                    "--edge-id", "0",
                ],
                check=True,
                cwd=REPO,
            )

            rows = self._read_csv(out / "local_avalanche_factors.csv")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["edge_id"], "0")
            self.assertGreater(float(rows[0]["vela_source_density_m3_s"]), 0.0)
            self.assertGreater(float(rows[0]["sentaurus_weighted_alpha_m_inv"]), 0.0)
            summary = json.loads((out / "local_avalanche_factors.json").read_text())
            self.assertEqual(summary["edge_id"], 0)
            self.assertEqual(summary["rows"][0]["sentaurus_node_id"], 0)

    def test_pn2d_bv_continuity_feedback_diagnostic_writes_local_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_continuity_feedback_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_0v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.2, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }))
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            (vela / "mini_0000_0V.vtk").write_text(
                """
# vtk DataFile Version 2.0
continuity feedback
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
2e-7 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0.0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.2 0.1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 -0.1 -0.2
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e30 2e30 3e30
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
1e25 1e25 1e25
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.2, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.05, 0.0],
                "eQuasiFermiPotential": [0.0, 0.1, 0.05],
                "hQuasiFermiPotential": [0.0, -0.05, -0.1],
                "eDensity": [1.0e14, 2.0e14, 3.0e14],
                "hDensity": [3.0e14, 2.0e14, 1.0e14],
                "ImpactIonization": [1.0e12, 2.0e12, 3.0e12],
                "srhRecombination": [1.0e10, 1.0e10, 1.0e10],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            for name, values in {
                "eCurrentDensity": [(1.0e-6, 0.0), (2.0e-6, 0.0), (3.0e-6, 0.0)],
                "hCurrentDensity": [(4.0e-6, 0.0), (5.0e-6, 0.0), (6.0e-6, 0.0)],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0", "component1"],
                    [[idx, x, y] for idx, (x, y) in enumerate(values)],
                )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_continuity_feedback.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "0",
                    "--edge-id", "0",
                ],
                check=True,
                cwd=REPO,
            )

            edge_rows = self._read_csv(out / "continuity_feedback_edges.csv")
            node_rows = self._read_csv(out / "continuity_feedback_nodes.csv")
            summary = json.loads((out / "continuity_feedback_summary.json").read_text())

        focus_edge = next(row for row in edge_rows if row["relation"] == "focus")
        self.assertEqual(focus_edge["edge_id"], "0")
        self.assertAlmostEqual(float(focus_edge["vela_dphin_V"]), 0.2)
        self.assertAlmostEqual(float(focus_edge["sentaurus_dphin_V"]), 0.1)
        self.assertGreater(float(focus_edge["vela_electron_flux_abs_m2_s"]), 0.0)
        focus_node = next(row for row in node_rows if row["node_id"] == "0")
        self.assertIn("vela_electron_residual_estimate_s_inv", focus_node)
        self.assertGreater(float(focus_node["sentaurus_generation_node_integral_s_inv"]), 0.0)
        self.assertEqual(summary["edge_id"], 0)
        self.assertGreaterEqual(summary["edge_rows"], 3)

    def test_pn2d_bv_active_support_continuity_balance_summarizes_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_support_continuity_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            sg_edges = root / "sg_edges.csv"
            edge_local_source = root / "edge_local_source.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [
                [0, 0.0, 0.0, 1.0e12, 8.0e11, 1, 1, "overlap"],
                [2, 0.0, 1.0, 0.0, 0.0, 0, 0, "inactive"],
            ])
            self._write_csv(sg_edges, [
                "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um",
                "x1_um", "y1_um", "edge_area_proxy_m2", "electron_alpha_m_inv",
                "hole_alpha_m_inv", "electron_flux_proxy", "hole_flux_proxy",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 1.0e-15, 1.0e6, 2.0e6, 3.0e14, 4.0e14],
                [-13.2, 1, 0, 2, 0.0, 0.0, 0.0, 1.0, 1.0e-15, 1.0e2, 2.0e2, 3.0e10, 4.0e10],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
support continuity
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0.0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.2 0.1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 -0.1 -0.2
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
8e17 2e17 1e17
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
1e15 1e15 1e15
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.05, 0.0],
                "eQuasiFermiPotential": [0.0, 0.1, 0.05],
                "hQuasiFermiPotential": [0.0, -0.05, -0.1],
                "eDensity": [1.0e14, 2.0e14, 3.0e14],
                "hDensity": [3.0e14, 2.0e14, 1.0e14],
                "ImpactIonization": [1.0e12, 2.0e12, 3.0e12],
                "srhRecombination": [1.0e10, 1.0e10, 1.0e10],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            for name, values in {
                "eCurrentDensity": [(1.0e-6, 0.0), (2.0e-6, 0.0), (3.0e-6, 0.0)],
                "hCurrentDensity": [(4.0e-6, 0.0), (5.0e-6, 0.0), (6.0e-6, 0.0)],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0", "component1"],
                    [[idx, x, y] for idx, (x, y) in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_support_continuity_balance.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(sg_edges),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "active_support_continuity_balance_nodes.csv")
            summary = json.loads((out / "active_support_continuity_balance_summary.json").read_text())

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["support_class"], "overlap")
        self.assertEqual(rows[0]["active_edge_ids"], "0")
        self.assertIn("vela_electron_residual_estimate_s_inv", rows[0])
        self.assertGreater(float(rows[0]["sentaurus_generation_node_integral_s_inv"]), 0.0)
        self.assertEqual(summary["support_classes"]["overlap"]["count"], 1)

    def test_pn2d_bv_active_support_residual_proxy_compares_state_and_source_policy(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_residual_proxy_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            support = root / "support.csv"
            sg_edges = root / "sg_edges.csv"
            edge_local_source = root / "edge_local_source.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            vt = 8.617333262145e-5 * 300.0
            hole_shift = vt * math.log(2.0)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [
                [0, 0.0, 0.0, 1.0e12, 8.0e11, 1, 1, "overlap"],
                [2, 0.0, 1.0, 0.0, 0.0, 0, 0, "inactive"],
            ])
            self._write_csv(sg_edges, [
                "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um",
                "x1_um", "y1_um", "edge_area_proxy_m2", "electron_alpha_m_inv",
                "hole_alpha_m_inv", "electron_flux_proxy", "hole_flux_proxy",
                "edge_source_integral",
            ], [
                [-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 2.0e-15, 1.0e6, 2.0e6, 3.0e14, 4.0e14, 1.1e6],
                [-13.2, 1, 0, 2, 0.0, 0.0, 0.0, 1.0, 2.0e-15, 1.0e2, 2.0e2, 3.0e10, 4.0e10, 1.1e2],
            ])
            self._write_csv(edge_local_source, [
                "node_id", "current_cxx_endpoint_source", "sentaurus_edgeavg_current_scaled_source",
            ], [[0, 2.5e5, 7.5e5]])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
residual proxy
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0.0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.2 0.1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 -0.1 -0.2
SCALARS Electrons float 1
LOOKUP_TABLE default
4 2 4
SCALARS Holes float 1
LOOKUP_TABLE default
3 2 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
1e9 1e9 1e9
""".lstrip(),
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.05, 0.0],
                "eQuasiFermiPotential": [0.0, 0.1, 0.05],
                "hQuasiFermiPotential": [hole_shift, -0.05 + hole_shift, -0.1 + hole_shift],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [6.0e-6, 4.0e-6, 6.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
                "ImpactIonization": [4.0e12, 2.0e12, 3.0e12],
                "srhRecombination": [1.0e10, 1.0e10, 1.0e10],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_support_residual_proxy.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(sg_edges),
                    "--mesh", str(mesh),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--edge-local-source-csv", str(edge_local_source),
                    "--hole-qf-shift-v", str(hole_shift),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "active_support_residual_proxy_nodes.csv")
            summary = json.loads((out / "active_support_residual_proxy_summary.json").read_text())

        by_variant = {row["variant"]: row for row in rows}
        self.assertIn("vela_state_vela_edge_source", by_variant)
        self.assertIn("vela_state_sentaurus_node_source", by_variant)
        self.assertIn("sentaurus_state_sentaurus_node_source", by_variant)
        self.assertIn("shifted_vela_qf_sentaurus_node_source", by_variant)
        self.assertIn("probe_vela_density_vela_mobility_vela_source", by_variant)
        self.assertIn("probe_sentaurus_density_vela_mobility_vela_source", by_variant)
        self.assertIn("probe_sentaurus_density_sentaurus_mobility_vela_source", by_variant)
        self.assertIn("probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source", by_variant)
        self.assertAlmostEqual(
            float(by_variant["vela_state_vela_edge_source"]["impact_source_s_inv"]),
            5.5e5,
        )
        self.assertAlmostEqual(
            float(by_variant["probe_sentaurus_density_sentaurus_mobility_vela_source"]["impact_source_s_inv"]),
            5.5e5,
        )
        self.assertAlmostEqual(
            float(
                by_variant[
                    "probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source"
                ]["impact_source_s_inv"]
            ),
            1.65e6,
        )
        self.assertGreater(
            float(by_variant["sentaurus_state_sentaurus_node_source"]["impact_source_s_inv"]),
            float(by_variant["vela_state_vela_edge_source"]["impact_source_s_inv"]),
        )
        self.assertIn("electron_residual_over_impact", by_variant["sentaurus_state_sentaurus_node_source"])
        self.assertEqual(summary["support_classes"]["overlap"]["variants"]["vela_state_vela_edge_source"]["count"], 1)

    def test_pn2d_bv_active_support_residual_proxy_splits_transport_modes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_residual_transport_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            sg_edges = root / "sg_edges.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus_-13.2v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 5.0e18],
                [1, 5.0e18, 0.0],
                [2, 0.0, 5.0e18],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [[0, 0.0, 0.0, 1.0e12, 8.0e11, 1, 1, "overlap"]])
            self._write_csv(sg_edges, [
                "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um",
                "x1_um", "y1_um", "edge_area_proxy_m2", "electron_alpha_m_inv",
                "hole_alpha_m_inv", "electron_flux_proxy", "hole_flux_proxy",
                "edge_source_integral",
            ], [[-13.2, 0, 0, 1, 0.0, 0.0, 1.0, 0.0, 2.0e-15, 1.0e6, 2.0e6, 3.0e14, 4.0e14, 1.1e6]])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
residual transport proxy
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0.0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.2 0.1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 -0.1 -0.2
SCALARS Electrons float 1
LOOKUP_TABLE default
4 2 4
SCALARS Holes float 1
LOOKUP_TABLE default
3 2 3
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
1 1 1
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
1e9 1e9 1e9
""".lstrip(),
                encoding="utf-8",
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.05, 0.0],
                "eQuasiFermiPotential": [0.0, 0.1, 0.05],
                "hQuasiFermiPotential": [0.0, -0.05, -0.1],
                "eDensity": [4.0e-6, 2.0e-6, 4.0e-6],
                "hDensity": [6.0e-6, 4.0e-6, 6.0e-6],
                "eMobility": [1.0e4, 1.0e4, 1.0e4],
                "hMobility": [1.0e4, 1.0e4, 1.0e4],
                "ImpactIonization": [4.0e12, 2.0e12, 3.0e12],
                "srhRecombination": [1.0e10, 1.0e10, 1.0e10],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_support_residual_proxy.py"),
                    "--support-csv", str(support),
                    "--sg-edge-csv", str(sg_edges),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--transport-modes", "density_sg,qf_inferred_ni,qf_old_slotboom_ni",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "active_support_residual_proxy_nodes.csv")
            summary = json.loads((out / "active_support_residual_proxy_summary.json").read_text())

        modes = {row["transport_model"] for row in rows if row["variant"] == "vela_state_sentaurus_node_source"}
        self.assertEqual(modes, {"density_sg", "qf_inferred_ni", "qf_old_slotboom_ni"})
        by_mode = {
            row["transport_model"]: row
            for row in rows
            if row["variant"] == "vela_state_sentaurus_node_source"
        }
        self.assertNotEqual(
            by_mode["qf_inferred_ni"]["electron_transport_s_inv"],
            by_mode["qf_old_slotboom_ni"]["electron_transport_s_inv"],
        )
        variant_summary = summary["support_classes"]["overlap"]["variants"]["vela_state_sentaurus_node_source"]
        self.assertEqual(variant_summary["transport_models"]["density_sg"]["count"], 1)
        self.assertEqual(variant_summary["transport_models"]["qf_old_slotboom_ni"]["count"], 1)

    def test_pn2d_bv_active_support_sensitivity_reports_term_derivatives(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_support_sensitivity_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            vela = root / "vela"
            out = root / "reports"
            vela.mkdir()

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0e-6, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0e-6},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [
                [0, 0.0, 0.0, 1.0e12, 8.0e11, 1, 1, "overlap"],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
support sensitivity
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.1 0.0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.2 0.1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0 -0.1 -0.2
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
8e17 2e17 1e17
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
1e15 1e15 1e15
""".lstrip()
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_active_support_sensitivity.py"),
                    "--support-csv", str(support),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--nodes", "0",
                    "--delta-v", "1e-3",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "active_support_sensitivity_rows.csv")
            summary = json.loads((out / "active_support_sensitivity_summary.json").read_text())

        self.assertEqual({row["carrier"] for row in rows}, {"electron", "hole"})
        electron = next(row for row in rows if row["carrier"] == "electron")
        hole = next(row for row in rows if row["carrier"] == "hole")
        self.assertEqual(electron["node_id"], "0")
        self.assertNotEqual(float(electron["d_flux_dqf_s_inv_per_V"]), 0.0)
        self.assertNotEqual(float(hole["d_flux_dqf_s_inv_per_V"]), 0.0)
        self.assertIn("electron", summary["by_carrier"])
        self.assertIn("max_abs_d_residual_dqf_s_inv_per_V", summary["by_carrier"]["electron"])

    def test_pn2d_bv_sg_coupling_paths_reports_contact_path_drops(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_coupling_paths_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                    {"id": 1, "region_id": 0, "node_ids": [1, 3, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "region_id": 0, "node_ids": [0]},
                    {"id": 1, "name": "Cathode", "region_id": 0, "node_ids": [3]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 1.0e17, 0.0],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [
                [1, 1.0, 0.0, 1.0e12, 8.0e11, 1, 1, "overlap"],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
coupling paths
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
0.0 0.1 0.0 0.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0.0 0.05 0.0 0.10
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0.0 -0.04 0.0 -0.08
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 1e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21 1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e17 8e17 1e17 2e17
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
                [3, 1.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.11, 0.0, 0.21],
                "eQuasiFermiPotential": [0.0, 0.04, 0.0, 0.12],
                "hQuasiFermiPotential": [0.0, -0.03, 0.0, -0.09],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_coupling_paths.py"),
                    "--support-csv", str(support),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--nodes", "1",
                    "--target-contact", "Cathode",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            edge_rows = self._read_csv(out / "sg_coupling_path_edges.csv")
            summary_rows = self._read_csv(out / "sg_coupling_path_summary.csv")
            summary = json.loads((out / "sg_coupling_path_summary.json").read_text())

        self.assertGreaterEqual(len(edge_rows), 2)
        self.assertEqual({row["carrier"] for row in summary_rows}, {"electron", "hole"})
        electron_path = next(row for row in summary_rows if row["carrier"] == "electron")
        self.assertEqual(electron_path["start_node"], "1")
        self.assertEqual(electron_path["target_contact"], "Cathode")
        self.assertGreater(float(electron_path["min_coupling_abs_s_inv_per_V"]), 0.0)
        self.assertIn("sentaurus_total_qf_drop_V", electron_path)
        self.assertEqual(summary["path_count"], 2)

    def test_pn2d_bv_sg_coupling_paths_replays_sentaurus_state_for_coupling(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_coupling_replay_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            support = root / "support.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                    {"id": 3, "x": 1.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                    {"id": 1, "region_id": 0, "node_ids": [1, 3, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "region_id": 0, "node_ids": [0]},
                    {"id": 1, "name": "Cathode", "region_id": 0, "node_ids": [3]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 1.0e17, 0.0],
            ])
            self._write_csv(support, [
                "node_id", "x_um", "y_um", "sentaurus_avalanche_cm3_s",
                "vela_avalanche_cm3_s", "sentaurus_active", "vela_active",
                "support_class",
            ], [
                [1, 1.0, 0.0, 1.0e12, 8.0e11, 1, 1, "overlap"],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
coupling replay
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
CELLS 2 8
3 0 1 2
3 1 3 2
CELL_TYPES 2
5
5
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
0.0 0.1 0.0 0.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0.0 0.05 0.0 0.10
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0.0 -0.04 0.0 -0.08
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 1e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
-1e21 1e21 -1e21 1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e17 8e17 1e17 2e17
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.0, 0.0],
                [1, 1.0, 0.0],
                [2, 0.0, 1.0],
                [3, 1.0, 1.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.12, 0.0, 0.22],
                "eQuasiFermiPotential": [0.0, 0.06, 0.0, 0.13],
                "hQuasiFermiPotential": [0.0, -0.02, 0.0, -0.07],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_coupling_paths.py"),
                    "--support-csv", str(support),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--nodes", "1",
                    "--target-contact", "Cathode",
                    "--coupling-state", "sentaurus",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            edge_rows = self._read_csv(out / "sg_coupling_path_edges.csv")
            summary = json.loads((out / "sg_coupling_path_summary.json").read_text())

        self.assertEqual(summary["coupling_state"], "sentaurus")
        self.assertTrue(edge_rows)
        self.assertTrue(all(row["coupling_state"] == "sentaurus" for row in edge_rows))
        self.assertGreater(float(edge_rows[0]["coupling_abs_s_inv_per_V"]), 0.0)
        self.assertIn("coupling_density_exp_from", edge_rows[0])
        self.assertIn("coupling_eta", edge_rows[0])
        self.assertIn("coupling_derivative_from_s_inv_per_V", edge_rows[0])
        self.assertIn("coupling_forward_flux_from_m2_s", edge_rows[0])

    def test_python_sg_variable_ni_qf_flux_matches_density_form_at_large_bias(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_sg_avalanche_edges.py"
        spec = importlib.util.spec_from_file_location("diagnose_pn2d_bv_sg_avalanche_edges", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        sgdiag = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sgdiag)

        vt = 8.617333262145e-5 * 300.0
        ni = 1.6556319846864456e16
        coef = 1.0

        psi0 = -13.555289531002014
        psi1 = -13.596737204970617
        phin0 = -13.060389803084155
        phin1 = -13.1147070328538
        phip0 = -12.973676068344172
        phip1 = -13.199999109476362

        n0 = ni * math.exp((psi0 - phin0) / vt)
        n1 = ni * math.exp((psi1 - phin1) / vt)
        p0 = ni * math.exp((phip0 - psi0) / vt)
        p1 = ni * math.exp((phip1 - psi1) / vt)

        electron_density = sgdiag.sg_electron_flux(n0, n1, psi1 - psi0, vt, coef)
        electron_qf = sgdiag.sg_electron_flux_qf_variable_ni(
            ni, ni, psi0, psi1, phin0, phin1, vt, coef)
        hole_density = sgdiag.sg_hole_flux(p0, p1, psi1 - psi0, vt, coef)
        hole_qf = sgdiag.sg_hole_flux_qf_variable_ni(
            ni, ni, psi0, psi1, phip0, phip1, vt, coef)

        self.assertNotEqual(electron_density, 0.0)
        self.assertNotEqual(hole_density, 0.0)
        self.assertAlmostEqual(electron_qf, electron_density, delta=abs(electron_density) * 1.0e-12)
        self.assertAlmostEqual(hole_qf, hole_density, delta=abs(hole_density) * 1.0e-12)

    def test_pn2d_bv_transition_edge_state_comparator_ranks_state_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_transition_edge_state_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            path_edges = root / "path_edges.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.70, "y": 0.0},
                    {"id": 1, "x": 0.75, "y": 0.02},
                    {"id": 2, "x": 0.80, "y": 0.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(path_edges, [
                "bias_V", "coupling_state", "carrier", "start_node", "target_contact",
                "path_found", "step", "edge_id", "node_from", "node_to",
                "x_from_um", "y_from_um", "x_to_um", "y_to_um",
                "delta_qf_drop_V",
            ], [
                [-13.2, "vela", "electron", 2, "Anode", "True", 1, 10, 2, 1, 0.80, 0.0, 0.75, 0.0, -0.02],
                [-13.2, "vela", "electron", 2, "Anode", "True", 2, 11, 1, 0, 0.75, 0.0, 0.70, 0.0, 0.01],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
transition state
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0.70e-6 0 0
0.75e-6 0.02e-6 0
0.80e-6 0 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
-1.0 -1.2 -1.5
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.8 -0.95 -1.10
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-1.1 -1.25 -1.45
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e17 2e17 3e17
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.70, 0.0],
                [1, 0.75, 0.0],
                [2, 0.80, 0.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [-1.02, -1.18, -1.44],
                "eQuasiFermiPotential": [-0.82, -0.94, -1.05],
                "hQuasiFermiPotential": [-1.12, -1.22, -1.40],
                "eDensity": [1.1e14, 2.1e14, 3.1e14],
                "hDensity": [3.1e14, 2.1e14, 1.1e14],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_transition_edge_state.py"),
                    "--path-edge-csv", str(path_edges),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--carrier", "electron",
                    "--x-min-um", "0.70",
                    "--x-max-um", "0.80",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "transition_edge_state_compare.csv")
            summary = json.loads((out / "transition_edge_state_summary.json").read_text())

        self.assertEqual(len(rows), 2)
        self.assertIn("delta_electron_qf_field_V_m", rows[0])
        self.assertIn("delta_electric_field_V_m", rows[0])
        self.assertIn("vela_electron_density_avg_cm3", rows[0])
        self.assertIn("sentaurus_electron_density_avg_cm3", rows[0])
        self.assertIn("hybrid_vela_psi_sentaurus_phin_delta_exp_avg_V", rows[0])
        self.assertIn("hybrid_sentaurus_psi_vela_phin_delta_exp_avg_V", rows[0])
        self.assertIn("dominant_drop_mismatch_metric", summary)
        self.assertGreater(summary["edge_count"], 0)

    def test_pn2d_bv_transition_flux_replay_compares_hybrid_states(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_transition_flux_replay_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            path_edges = root / "path_edges.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.70, "y": 0.0},
                    {"id": 1, "x": 0.75, "y": 0.02},
                    {"id": 2, "x": 0.80, "y": 0.0},
                ],
                "triangles": [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]}],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(path_edges, [
                "bias_V", "coupling_state", "carrier", "start_node", "target_contact",
                "path_found", "step", "edge_id", "node_from", "node_to",
                "x_from_um", "y_from_um", "x_to_um", "y_to_um",
                "length_m", "couple_m", "delta_qf_drop_V",
            ], [
                [-13.2, "vela", "electron", 2, "Anode", "True", 1, 10, 2, 1, 0.80, 0.0, 0.75, 0.0, 5.0e-8, 2.0e-8, -0.02],
                [-13.2, "vela", "electron", 2, "Anode", "True", 2, 11, 1, 0, 0.75, 0.0, 0.70, 0.0, 5.0e-8, 2.0e-8, 0.01],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
transition flux replay
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0.70e-6 0 0
0.75e-6 0 0
0.80e-6 0 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
-1.0 -1.2 -1.5
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.8 -0.95 -1.10
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-1.1 -1.25 -1.45
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e17 2e17 3e17
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
2e16 4e16 6e16
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.70, 0.0],
                [1, 0.75, 0.02],
                [2, 0.80, 0.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [-1.02, -1.18, -1.44],
                "eQuasiFermiPotential": [-0.82, -0.94, -1.05],
                "ImpactIonization": [1.0e11, 2.0e11, 3.0e11],
                "srhRecombination": [1.0e10, 2.0e10, 3.0e10],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_transition_flux_replay.py"),
                    "--path-edge-csv", str(path_edges),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--start-node", "2",
                    "--target-contact", "Anode",
                    "--x-min-um", "0.70",
                    "--x-max-um", "0.80",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "transition_flux_replay_edges.csv")
            summary = json.loads((out / "transition_flux_replay_summary.json").read_text())

        self.assertEqual(
            {row["state"] for row in rows},
            {"vela", "sentaurus", "vela_psi_sentaurus_phin", "sentaurus_psi_vela_phin"},
        )
        self.assertIn("electron_flux_integral_s_inv", rows[0])
        self.assertIn("electron_flux_integral_delta_vs_sentaurus_s_inv", rows[0])
        self.assertIn("best_flux_match_state", summary)
        self.assertGreater(summary["edge_count"], 0)

    def test_pn2d_bv_transition_row_residual_decomposition_groups_node_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_transition_row_residual_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            path_edges = root / "path_edges.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.70, "y": 0.0},
                    {"id": 1, "x": 0.75, "y": 0.02},
                    {"id": 2, "x": 0.80, "y": 0.0},
                ],
                "triangles": [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]}],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            self._write_csv(path_edges, [
                "bias_V", "coupling_state", "carrier", "start_node", "target_contact",
                "path_found", "step", "edge_id", "node_from", "node_to",
                "x_from_um", "y_from_um", "x_to_um", "y_to_um",
                "length_m", "couple_m", "delta_qf_drop_V",
            ], [
                [-13.2, "vela", "electron", 2, "Anode", "True", 1, 10, 2, 1, 0.80, 0.0, 0.75, 0.0, 5.0e-8, 2.0e-8, -0.02],
                [-13.2, "vela", "electron", 2, "Anode", "True", 2, 11, 1, 0, 0.75, 0.0, 0.70, 0.0, 5.0e-8, 2.0e-8, 0.01],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
transition row residual
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0.70e-6 0 0
0.75e-6 0.02e-6 0
0.80e-6 0 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
-1.0 -1.2 -1.5
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.8 -0.95 -1.10
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-1.1 -1.25 -1.45
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20 3e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20 1e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1 0.1
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
SCALARS NetDoping float 1
LOOKUP_TABLE default
1e21 1e21 -1e21
SCALARS AvalancheGeneration float 1
LOOKUP_TABLE default
1e17 2e17 3e17
SCALARS SRHRecombination float 1
LOOKUP_TABLE default
2e16 4e16 6e16
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.70, 0.0],
                [1, 0.75, 0.02],
                [2, 0.80, 0.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [-1.02, -1.18, -1.44],
                "eQuasiFermiPotential": [-0.82, -0.94, -1.05],
                "ImpactIonization": [1.0e11, 2.0e11, 3.0e11],
                "srhRecombination": [1.0e10, 2.0e10, 3.0e10],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_transition_row_residual.py"),
                    "--path-edge-csv", str(path_edges),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--bias", "-13.2",
                    "--start-node", "2",
                    "--target-contact", "Anode",
                    "--x-min-um", "0.70",
                    "--x-max-um", "0.80",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "transition_row_residual_nodes.csv")
            summary = json.loads((out / "transition_row_residual_summary.json").read_text())

        self.assertEqual(
            {row["state"] for row in rows},
            {"vela", "sentaurus", "vela_psi_sentaurus_phin", "sentaurus_psi_vela_phin"},
        )
        self.assertIn("electron_residual_contribution_s_inv", rows[0])
        self.assertIn("electron_residual_delta_vs_sentaurus_s_inv", rows[0])
        self.assertIn("d_residual_d_self_phin_s_inv_per_V", rows[0])
        self.assertIn("local_newton_step_phin_V", rows[0])
        self.assertIn("vela_to_sentaurus_phin_delta_V", rows[0])
        self.assertIn("vela_avalanche_node_integral_s_inv", rows[0])
        self.assertIn("vela_srh_node_integral_s_inv", rows[0])
        self.assertIn("electron_full_residual_s_inv", rows[0])
        self.assertIn("sentaurus_electron_full_residual_s_inv", rows[0])
        self.assertIn("electron_full_residual_delta_vs_sentaurus_s_inv", rows[0])
        self.assertIn("best_residual_match_state", summary)
        self.assertIn("best_full_residual_match_state", summary)
        self.assertIn("median_abs_source_delta_vs_sentaurus_s_inv", summary)
        self.assertGreater(summary["node_count"], 0)
        vela_node = next(row for row in rows if row["state"] == "vela" and row["node_id"] == "1")
        self.assertGreater(float(vela_node["vela_avalanche_node_integral_s_inv"]), 0.0)
        expected_full = (
            float(vela_node["electron_residual_contribution_s_inv"])
            + float(vela_node["vela_srh_node_integral_s_inv"])
            - float(vela_node["vela_avalanche_node_integral_s_inv"])
        )
        self.assertAlmostEqual(float(vela_node["electron_full_residual_s_inv"]), expected_full)

    def test_pn2d_bv_transition_diagnostics_reject_empty_edge_selection(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_transition_empty_edges_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            path_edges = root / "path_edges.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            fields = sentaurus / "fields"
            vela.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.70, "y": 0.0},
                    {"id": 1, "x": 0.75, "y": 0.0},
                ],
                "triangles": [],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 0.0, 1.0e17],
            ])
            self._write_csv(path_edges, [
                "bias_V", "coupling_state", "carrier", "start_node", "target_contact",
                "path_found", "step", "edge_id", "node_from", "node_to",
                "x_from_um", "y_from_um", "x_to_um", "y_to_um",
                "length_m", "couple_m", "delta_qf_drop_V",
            ], [
                [-13.2, "vela", "electron", 1, "Cathode", "True", 1, 10, 1, 0,
                 0.75, 0.0, 0.70, 0.0, 5.0e-8, 2.0e-8, -0.02],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
transition empty selection
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0.70e-6 0 0
0.75e-6 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
-1.0 -1.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.8 -0.95
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-1.1 -1.25
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 2e20
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.1 0.1
""".lstrip()
            )
            self._write_csv(sentaurus / "nodes.csv", ["id", "x_um", "y_um"], [
                [0, 0.70, 0.0],
                [1, 0.75, 0.0],
            ])
            for name, values in {
                "ElectrostaticPotential": [-1.02, -1.18],
                "eQuasiFermiPotential": [-0.82, -0.94],
                "hQuasiFermiPotential": [-1.12, -1.22],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            common_args = [
                "--path-edge-csv", str(path_edges),
                "--mesh", str(mesh),
                "--doping-csv", str(doping),
                "--sentaurus-dir", str(sentaurus),
                "--vela-vtk-root", str(vela),
                "--bias", "-13.2",
                "--target-contact", "Anode",
            ]
            for script_name in (
                "diagnose_pn2d_bv_transition_edge_state.py",
                "diagnose_pn2d_bv_transition_flux_replay.py",
            ):
                result = subprocess.run(
                    [
                        sys.executable,
                        str(REPO / "scripts" / script_name),
                        *common_args,
                        "--out-dir", str(root / script_name),
                    ],
                    cwd=REPO,
                    text=True,
                    capture_output=True,
                )
                self.assertNotEqual(result.returncode, 0, script_name)
                self.assertIn("no transition path edges selected", result.stderr)

    def test_pn2d_bv_transition_row_residual_requires_vtk_mobility_scalar(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_transition_missing_scalar_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            path_edges = root / "path_edges.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus"
            vela.mkdir()
            sentaurus.mkdir()
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.70, "y": 0.0},
                    {"id": 1, "x": 0.75, "y": 0.0},
                ],
                "triangles": [],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 0.0, 1.0e17],
            ])
            self._write_csv(path_edges, [
                "bias_V", "coupling_state", "carrier", "start_node", "target_contact",
                "path_found", "step", "edge_id", "node_from", "node_to",
                "x_from_um", "y_from_um", "x_to_um", "y_to_um",
                "length_m", "couple_m", "delta_qf_drop_V",
            ], [
                [-13.2, "vela", "electron", 1, "Anode", "True", 1, 10, 1, 0,
                 0.75, 0.0, 0.70, 0.0, 5.0e-8, 2.0e-8, -0.02],
            ])
            (vela / "mini_0000_-13.2V.vtk").write_text(
                """
# vtk DataFile Version 2.0
transition missing scalar
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0.70e-6 0 0
0.75e-6 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
-1.0 -1.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.8 -0.95
""".lstrip()
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_transition_row_residual.py"),
                    "--path-edge-csv", str(path_edges),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-dir", str(sentaurus),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(root / "out"),
                    "--bias", "-13.2",
                    "--target-contact", "Anode",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("VTK is missing required scalars ElectronMobility", result.stderr)

    def test_pn2d_bv_transition_continuation_overshoot_flags_sparse_jump(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_transition_overshoot_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            vela = root / "vela"
            sentaurus_root = root / "sentaurus"
            out = root / "reports"
            vela.mkdir()
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.65, "y": 0.0},
                    {"id": 1, "x": 0.70, "y": 0.0},
                    {"id": 2, "x": 0.75, "y": 0.0},
                ],
                "triangles": [],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")

            for step, bias, phin_values, psi_values, electron_values in [
                (0, -12.9, [-1.001, -0.999, -1.001], [-1.2, -1.3, -1.4], [1.0e18, 2.0e18, 3.0e18]),
                (1, -13.2, [-1.020, -0.980, -1.020], [-1.25, -1.35, -1.45], [1.0e20, 2.0e20, 3.0e20]),
            ]:
                (vela / f"mini_{step:04d}_{bias:g}V.vtk").write_text(
                    f"""
# vtk DataFile Version 2.0
transition overshoot
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0.65e-6 0 0
0.70e-6 0 0
0.75e-6 0 0
CELLS 0 0
CELL_TYPES 0
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
{psi_values[0]} {psi_values[1]} {psi_values[2]}
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
{phin_values[0]} {phin_values[1]} {phin_values[2]}
SCALARS Electrons float 1
LOOKUP_TABLE default
{electron_values[0]} {electron_values[1]} {electron_values[2]}
""".lstrip()
                )

            for bias, psi_values, phin_values, density_values in [
                (-12.9, [-1.2, -1.3, -1.4], [-1.0, -1.0, -1.0], [1.0e12, 2.0e12, 3.0e12]),
                (-13.2, [-1.25, -1.35, -1.45], [-1.0, -1.0, -1.0], [1.0e14, 2.0e14, 3.0e14]),
            ]:
                sent = sentaurus_root / f"sentaurus_{bias:g}v"
                fields = sent / "fields"
                fields.mkdir(parents=True)
                self._write_csv(sent / "nodes.csv", ["id", "x_um", "y_um"], [
                    [0, 0.65, 0.0],
                    [1, 0.70, 0.0],
                    [2, 0.75, 0.0],
                ])
                for name, values in {
                    "ElectrostaticPotential": psi_values,
                    "eQuasiFermiPotential": phin_values,
                    "eDensity": density_values,
                }.items():
                    self._write_csv(
                        fields / f"{name}_region0.csv",
                        ["node_id", "component0"],
                        [[idx, value] for idx, value in enumerate(values)],
                    )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_transition_continuation_overshoot.py"),
                    "--mesh", str(mesh),
                    "--vela-vtk-root", str(vela),
                    "--sentaurus-root", str(sentaurus_root),
                    "--out-dir", str(out),
                    "--bias-min", "-13.2",
                    "--bias-max", "-12.9",
                    "--nodes", "0,1,2",
                    "--phin-threshold-v", "0.01",
                    "--max-accepted-gap-v", "0.1",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "transition_continuation_nodes.csv")
            summary = json.loads((out / "transition_continuation_summary.json").read_text())

        self.assertEqual(len(rows), 6)
        self.assertIn("delta_phin_V", rows[0])
        self.assertIn("delta_psi_minus_phin_V", rows[0])
        self.assertIn("log10_vela_over_sentaurus_electron_density", rows[0])
        self.assertAlmostEqual(summary["first_large_phin_mismatch_bias_V"], -13.2)
        self.assertAlmostEqual(summary["first_alternating_large_phin_mismatch_bias_V"], -13.2)
        self.assertTrue(summary["needs_intermediate_restart"])
        self.assertAlmostEqual(summary["largest_accepted_gap_V"], 0.3)

    def test_prepare_pn2d_bv_focused_restart_writes_restart_and_config(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_focused_restart_") as tmp:
            root = Path(tmp)
            base = root / "base.json"
            vtk = root / "restart.vtk"
            out = root / "focused"
            base.write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "contacts": [{"name": "Anode", "bias": 0.0}],
                "solver": {"method": "gummel_newton"},
                "sweep": {
                    "mode": "bv_reverse",
                    "contact": "Anode",
                    "current_contact": "Anode",
                    "start": -12.0,
                    "stop": -13.2,
                    "step": -0.2,
                    "write_vtk": True,
                    "vtk_prefix": "old/vtk",
                    "csv_file": "old.csv",
                    "initial_state_file": "old_restart.csv",
                    "write_state_file": "old_latest.csv",
                    "diagnostics": {
                        "newton_history": {"enabled": True, "csv_file": "old_newton.csv"},
                        "sg_avalanche_edges": {"enabled": True, "csv_file": "old_sg.csv"},
                    },
                },
            }) + "\n")
            vtk.write_text(
                """
# vtk DataFile Version 2.0
restart
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1e-6 0 0
CELLS 0 0
CELL_TYPES 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
0.1 0.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0.0 0.01
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-0.1 -0.2
SCALARS Electrons float 1
LOOKUP_TABLE default
1e20 2e20
SCALARS Holes float 1
LOOKUP_TABLE default
3e20 4e20
""".lstrip()
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "prepare_pn2d_bv_focused_restart.py"),
                    "--base-config", str(base),
                    "--restart-vtk", str(vtk),
                    "--out-dir", str(out),
                    "--bias-points", "-12.9078,-13.0,-13.1,-13.2",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            restart_rows = self._read_csv(out / "restart_from_vtk.csv")
            config = json.loads((out / "simulation.json").read_text())

        self.assertEqual(restart_rows[0]["node_id"], "0")
        self.assertAlmostEqual(float(restart_rows[1]["phin"]), 0.01)
        self.assertEqual(config["sweep"]["bias_points"], [-12.9078, -13.0, -13.1, -13.2])
        self.assertEqual(config["sweep"]["start"], -12.9078)
        self.assertEqual(config["sweep"]["stop"], -13.2)
        self.assertAlmostEqual(config["sweep"]["step"], -0.0922)
        self.assertTrue(config["sweep"]["initial_state_file"].endswith("restart_from_vtk.csv"))
        self.assertEqual(Path(config["sweep"]["vtk_prefix"]).parts[-2:], ("vtk", "focused_restart"))
        self.assertTrue(config["sweep"]["diagnostics"]["newton_history"]["csv_file"].endswith("newton_history.csv"))

    def test_pn2d_bv_newton_residual_state_diagnostic_prepares_probe_inputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_residual_states_") as tmp:
            root = Path(tmp)
            base = root / "base.json"
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_-1v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                ],
                "triangles": [],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Cathode", "region_id": 0, "node_ids": [0]},
                    {"id": 1, "name": "Anode", "region_id": 0, "node_ids": [1]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 0.0, 1.0e17],
            ])
            base.write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "output_csv": "unused.csv",
                "scaling": {"mode": "unit_scaling"},
                "mesh_geometry": {"node_volume_policy": "mixed_voronoi"},
                "contacts": [
                    {"name": "Cathode", "bias": 0.0},
                    {"name": "Anode", "bias": 0.0},
                ],
                "solver": {
                    "method": "gummel_newton",
                    "bandgap_narrowing": "old_slotboom",
                    "mobility": {"model": "masetti_field"},
                },
                "sweep": {"mode": "bv_reverse", "start": 0.0, "stop": -1.0, "step": -1.0},
            }, indent=2) + "\n")
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
residual states
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 2 float
0 0 0
1e-6 0 0
POINT_DATA 2
SCALARS Potential float 1
LOOKUP_TABLE default
-0.2 -0.1
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.3 -0.2
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-0.4 -0.3
""".lstrip()
            )
            for name, values in {
                "ElectrostaticPotential": [-0.25, -0.15],
                "eQuasiFermiPotential": [-0.35, -0.25],
                "hQuasiFermiPotential": [-0.45, -0.35],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_newton_residual_states.py"),
                    "--base-config", str(base),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--states", "vela:-1,sentaurus:-1,hybrid_vpsi_sqf:-1,hybrid_spsi_shift_vqf:-1",
                    "--nodes", "0,1",
                    "--prepare-only",
                ],
                check=True,
                cwd=REPO,
            )

            vela_state = self._read_csv(
                out / "states" / "vela_m1v" / "fields" / "ElectrostaticPotential_region0.csv"
            )
            config = json.loads(
                (out / "configs" / "sentaurus_m1v_newton_residual_probe.json").read_text()
            )
            hybrid_vpsi = self._read_csv(
                out / "states" / "hybrid_vpsi_sqf_m1v" / "fields" / "eQuasiFermiPotential_region0.csv"
            )
            hybrid_shift = self._read_csv(
                out / "states" / "hybrid_spsi_shift_vqf_m1v" / "fields" / "ElectrostaticPotential_region0.csv"
            )
            summary = json.loads((out / "newton_residual_state_summary.json").read_text())

        self.assertEqual(vela_state[0]["component0"], "-0.20000000000000001")
        self.assertEqual(hybrid_vpsi[0]["component0"], "-0.34999999999999998")
        self.assertEqual(hybrid_shift[0]["component0"], "-0.20000000000000001")
        self.assertEqual(config["simulation_type"], "newton_residual_probe")
        self.assertEqual(config["mesh_geometry"]["node_volume_policy"], "mixed_voronoi")
        self.assertNotIn("sweep", config)
        anode = next(contact for contact in config["contacts"] if contact["name"] == "Anode")
        self.assertEqual(anode["bias"], -1.0)
        self.assertEqual(summary["node_rows"], 0)
        self.assertEqual(len(summary["states"]), 4)
        shifted = next(item for item in summary["states"] if item["source"] == "hybrid_spsi_shift_vqf")
        self.assertAlmostEqual(shifted["psi_shift_V"], 0.05)

    def test_pn2d_bv_newton_residual_state_diagnostic_allows_vela_without_sentaurus(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_residual_vela_only_") as tmp:
            root = Path(tmp)
            base = root / "base.json"
            mesh = root / "mesh.json"
            vela = root / "vela"
            out = root / "reports"
            vela.mkdir()
            mesh.write_text(json.dumps({
                "nodes": [{"id": 0, "x": 0.0, "y": 0.0}],
                "triangles": [],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
            }) + "\n")
            base.write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "output_csv": "unused.csv",
                "contacts": [
                    {"name": "Cathode", "bias": 0.0},
                    {"name": "Anode", "bias": 0.0},
                ],
                "sweep": {"mode": "bv_reverse", "start": 0.0, "stop": -1.0, "step": -1.0},
            }, indent=2) + "\n")
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
residual states
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 1 float
0 0 0
POINT_DATA 1
SCALARS Potential float 1
LOOKUP_TABLE default
-0.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.3
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-0.4
""".lstrip()
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_newton_residual_states.py"),
                    "--base-config", str(base),
                    "--sentaurus-root", str(root / "missing_sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--states", "vela:-1",
                    "--nodes", "0",
                    "--prepare-only",
                ],
                check=True,
                cwd=REPO,
            )

            summary = json.loads((out / "newton_residual_state_summary.json").read_text())
            self.assertEqual(summary["states"][0]["source"], "vela")
            self.assertTrue(
                (out / "states" / "vela_m1v" / "fields" / "ElectrostaticPotential_region0.csv").exists()
            )

    def test_pn2d_bv_newton_residual_state_finds_m_token_sentaurus_exports(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_newton_residual_states.py"
        spec = importlib.util.spec_from_file_location("diagnose_pn2d_bv_newton_residual_states", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        scripts_dir = str(REPO / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_residual_sentaurus_token_") as tmp:
            root = Path(tmp)
            expected = root / "sentaurus_m12p7v"
            (expected / "fields").mkdir(parents=True)

            self.assertEqual(module.sentaurus_dir(root, -12.7), expected)

    def test_pn2d_bv_continuity_feedback_finds_m_token_sentaurus_exports(self) -> None:
        module_path = REPO / "scripts" / "diagnose_pn2d_bv_continuity_feedback.py"
        spec = importlib.util.spec_from_file_location("diagnose_pn2d_bv_continuity_feedback", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        scripts_dir = str(REPO / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_feedback_sentaurus_token_") as tmp:
            root = Path(tmp)
            expected = root / "sentaurus_m12p7v"
            expected.mkdir(parents=True)
            (expected / "nodes.csv").write_text("id,x_um,y_um\n0,0,0\n", encoding="utf-8")

            self.assertEqual(module.sentaurus_dir(root, -12.7), expected)

    def test_pn2d_bv_sg_flux_form_diagnostic_compares_state_sources(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_sg_flux_forms_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_m1v"
            fields = sentaurus / "fields"
            out = root / "out"
            vela.mkdir()
            fields.mkdir(parents=True)

            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n", encoding="utf-8")
            self._write_csv(doping, ["node_id", "donors_m3", "acceptors_m3"], [
                [0, 1.0e23, 0.0],
                [1, 1.0e23, 0.0],
                [2, 1.0e23, 0.0],
            ])
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
sg flux forms
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.025 0
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.10 -0.08 -0.10
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0.10 0.11 0.10
SCALARS Electrons float 1
LOOKUP_TABLE default
1e16 2e16 1e16
SCALARS Holes float 1
LOOKUP_TABLE default
1e10 1e10 1e10
SCALARS ElectronMobility float 1
LOOKUP_TABLE default
0.12 0.10 0.11
SCALARS HoleMobility float 1
LOOKUP_TABLE default
0.05 0.05 0.05
""".lstrip(),
                encoding="utf-8",
            )
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.025, 0.0],
                "eQuasiFermiPotential": [-0.1, -0.09, -0.1],
                "hQuasiFermiPotential": [0.1, 0.11, 0.1],
                "eDensity": [1.0e10, 1.5e10, 1.0e10],
                "hDensity": [1.0e4, 1.0e4, 1.0e4],
                "eMobility": [1200.0, 1000.0, 1100.0],
                "hMobility": [500.0, 500.0, 500.0],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_sg_flux_forms.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "-1",
                    "--edge-ids", "0",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "sg_flux_form_edges.csv")
            by_source = {row["source"]: row for row in rows}
            self.assertEqual(set(by_source), {"vela", "sentaurus"})
            self.assertIn("electron_flux_qf_model_abs", by_source["vela"])
            self.assertIn("electron_flux_qf_inferred_abs", by_source["sentaurus"])
            self.assertGreater(float(by_source["vela"]["electron_qf_field_V_m"]), 0.0)
            self.assertAlmostEqual(float(by_source["sentaurus"]["electron_mobility_m2_V_s"]), 0.11)
            summary = json.loads((out / "sg_flux_form_summary.json").read_text())
            self.assertEqual(summary["edge_ids"], [0])
            self.assertEqual(summary["row_count"], 2)

    def test_pn2d_bv_accepted_state_evolution_reports_density_jump(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_accepted_evolution_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            vtk_root = root / "vtk"
            out = root / "out"
            iv = root / "iv.csv"
            vtk_root.mkdir()
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n", encoding="utf-8")
            (vtk_root / "state_0001_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
evolution before
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Electrons float 1
LOOKUP_TABLE default
1e8 1e8 1e8
SCALARS Holes float 1
LOOKUP_TABLE default
1e8 1e8 1e8
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.1 0
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.2 0
""".lstrip(),
                encoding="utf-8",
            )
            (vtk_root / "state_0002_-1.1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
evolution after
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Electrons float 1
LOOKUP_TABLE default
1e9 1e10 1e8
SCALARS Holes float 1
LOOKUP_TABLE default
1e8 1e8 1e8
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
0 0.3 0
SCALARS Potential float 1
LOOKUP_TABLE default
0 0.25 0
""".lstrip(),
                encoding="utf-8",
            )
            self._write_csv(iv, [
                "bias_V",
                "current_total_A_per_um",
                "iterations",
                "step_diagnostics",
                "converged",
            ], [
                [-1.0, -1.0e-18, 3, "attempted_step=-0.1;accepted_step=-0.1;retry_count=0", 1],
                [-1.10004, -2.0e-18, 5, "attempted_step=-0.1;accepted_step=-0.1;retry_count=1", 1],
            ])

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_accepted_state_evolution.py"),
                    "--mesh", str(mesh),
                    "--vtk-root", str(vtk_root),
                    "--iv-csv", str(iv),
                    "--out-dir", str(out),
                    "--bias-min", "-1.1",
                    "--bias-max", "-1.0",
                    "--edge-ids", "0",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            rows = self._read_csv(out / "accepted_state_evolution.csv")
            self.assertEqual(len(rows), 2)
            after = rows[1]
            self.assertAlmostEqual(float(after["electron_density_jump_median_dex"]), 1.0)
            self.assertAlmostEqual(float(after["electron_density_jump_max_dex"]), 2.0)
            self.assertEqual(after["top_electron_jump_node_id"], "1")
            self.assertEqual(after["iterations"], "5")
            edge_rows = self._read_csv(out / "accepted_state_edge_evolution.csv")
            after_edge = edge_rows[1]
            self.assertAlmostEqual(float(after_edge["electron_qf_drop_delta_V"]), 0.2)
            self.assertAlmostEqual(float(after_edge["electron_density_avg_jump_dex"]), math.log10(55.0 / 1.0))
            summary = json.loads((out / "accepted_state_evolution_summary.json").read_text())
            self.assertEqual(summary["state_count"], 2)
            self.assertEqual(summary["max_electron_density_jump"]["node_id"], 1)

    def test_pn2d_bv_potential_profile_diagnostic_writes_plateaus(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_potential_profiles_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_-1v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
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
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
potential profiles
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
0.0 0.2 0.0 0.2
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-0.1 0.1 -0.1 0.1
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
0.1 0.3 0.1 0.3
""".lstrip()
            )
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.1, 0.0, 0.1],
                "eQuasiFermiPotential": [-0.2, 0.0, -0.2, 0.0],
                "hQuasiFermiPotential": [0.2, 0.2, 0.2, 0.2],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_potential_profiles.py"),
                    "--mesh", str(mesh),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "-1",
                    "--horizontal-y-um", "0",
                    "--vertical-x-um", "1",
                    "--bands", "left:0:0,right:1:1",
                ],
                check=True,
                cwd=REPO,
            )

            profile = self._read_csv(out / "potential_profile_samples.csv")
            plateaus = self._read_csv(out / "potential_plateau_offsets.csv")
            summary = json.loads((out / "potential_profile_summary.json").read_text())

        self.assertTrue(profile)
        self.assertIn("delta_psi_minus_phin_V", profile[0])
        right = next(row for row in plateaus if row["band"] == "right")
        self.assertAlmostEqual(float(right["delta_psi_median_V"]), 0.1)
        self.assertAlmostEqual(float(right["delta_psi_minus_phin_median_V"]), 0.0)
        self.assertEqual(summary["plateau_rows"], 2)

    def test_pn2d_bv_branch_profile_diagnostic_compares_three_states(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_branch_profiles_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            avalanche = root / "avalanche"
            noimpact = root / "noimpact"
            sentaurus = root / "sentaurus" / "sentaurus_-1v"
            fields = sentaurus / "fields"
            out = root / "reports"
            avalanche.mkdir()
            noimpact.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
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
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [],
            }) + "\n")
            vtk_template = """
# vtk DataFile Version 2.0
branch profiles
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
1e-6 0 0
0 1e-6 0
1e-6 1e-6 0
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
{psi}
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
{phin}
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
{phip}
SCALARS Electrons float 1
LOOKUP_TABLE default
{electrons}
SCALARS Holes float 1
LOOKUP_TABLE default
{holes}
""".lstrip()
            (avalanche / "mini_0000_-1V.vtk").write_text(vtk_template.format(
                psi="0.0 0.4 0.0 0.4",
                phin="-0.1 0.2 -0.1 0.2",
                phip="0.1 0.5 0.1 0.5",
                electrons="1e18 1e20 1e18 1e20",
                holes="1e17 1e19 1e17 1e19",
            ))
            (noimpact / "mini_0000_-1V.vtk").write_text(vtk_template.format(
                psi="0.0 0.2 0.0 0.2",
                phin="-0.1 0.0 -0.1 0.0",
                phip="0.1 0.3 0.1 0.3",
                electrons="1e18 1e18 1e18 1e18",
                holes="1e17 1e17 1e17 1e17",
            ))
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.1, 0.0, 0.1],
                "eQuasiFermiPotential": [-0.1, -0.1, -0.1, -0.1],
                "hQuasiFermiPotential": [0.1, 0.2, 0.1, 0.2],
                "eDensity": [1.0e12, 1.0e12, 1.0e12, 1.0e12],
                "hDensity": [1.0e11, 1.0e11, 1.0e11, 1.0e11],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_branch_profiles.py"),
                    "--mesh", str(mesh),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-avalanche-vtk-root", str(avalanche),
                    "--vela-no-impact-vtk-root", str(noimpact),
                    "--out-dir", str(out),
                    "--biases", "-1",
                    "--bands", "right:1:1",
                    "--focus-nodes", "1",
                ],
                check=True,
                cwd=REPO,
            )

            medians = self._read_csv(out / "branch_profile_medians.csv")
            comparison = self._read_csv(out / "branch_profile_comparison.csv")
            focus = self._read_csv(out / "branch_profile_focus_nodes.csv")
            summary = json.loads((out / "branch_profile_summary.json").read_text())

        self.assertEqual(len(medians), 3)
        self.assertEqual(len(focus), 3)
        row = comparison[0]
        self.assertAlmostEqual(float(row["avalanche_minus_sentaurus_psi_median_V"]), 0.3)
        self.assertAlmostEqual(float(row["noimpact_minus_sentaurus_psi_median_V"]), 0.1)
        self.assertAlmostEqual(float(row["avalanche_log10_electron_over_noimpact_median"]), 2.0)
        self.assertEqual(summary["comparison_rows"], 1)

    def test_pn2d_bv_noimpact_variant_scan_prepares_isolation_configs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_noimpact_variant_scan_") as tmp:
            root = Path(tmp)
            base = root / "simulation.json"
            out = root / "scan"
            base.write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "output_csv": "base.csv",
                "solver": {
                    "method": "gummel_newton",
                    "mobility": {
                        "model": "masetti_field",
                        "high_field_driving_force": "quasi_fermi_gradient",
                    },
                    "recombination": ["srh"],
                    "impact_ionization": {"model": "none"},
                    "contact_boundary_minority_electron_relaxation": False,
                },
                "sweep": {
                    "mode": "bv_reverse",
                    "contact": "Anode",
                    "current_contact": "Anode",
                    "start": 0.0,
                    "stop": -20.0,
                    "step": -0.1,
                    "write_vtk": False,
                },
            }) + "\n")

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_noimpact_variant_scan.py"),
                    "--base-config", str(base),
                    "--runner", str(REPO / "build-release" / "vela_example_runner.exe"),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--out-dir", str(out),
                    "--biases", "-1,-2",
                    "--variants",
                    "baseline,srh_tau_equal_1e_m6,no_srh,low_field_masetti,contact_relax_n_0p10",
                    "--prepare-only",
                ],
                check=True,
                cwd=REPO,
            )

            baseline = json.loads(
                (out / "cases" / "baseline" / "baseline.json").read_text()
            )
            no_srh = json.loads((out / "cases" / "no_srh" / "no_srh.json").read_text())
            tau = json.loads(
                (out / "cases" / "srh_tau_equal_1e_m6" /
                 "srh_tau_equal_1e_m6.json").read_text()
            )
            low_field = json.loads(
                (out / "cases" / "low_field_masetti" / "low_field_masetti.json").read_text()
            )
            relax = json.loads(
                (out / "cases" / "contact_relax_n_0p10" /
                 "contact_relax_n_0p10.json").read_text()
            )
            summary = json.loads((out / "noimpact_variant_scan_summary.json").read_text())

        self.assertEqual(baseline["sweep"]["bias_points"], [0.0, -1.0, -2.0])
        self.assertTrue(baseline["sweep"]["write_vtk"])
        self.assertEqual(baseline["solver"]["impact_ionization"]["model"], "none")
        self.assertEqual(tau["solver"]["taun"], 1.0e-6)
        self.assertEqual(tau["solver"]["taup"], 1.0e-6)
        self.assertEqual(no_srh["solver"]["recombination"], ["none"])
        self.assertEqual(low_field["solver"]["mobility"]["model"], "masetti")
        self.assertTrue(relax["solver"]["contact_boundary_minority_electron_relaxation"])
        self.assertEqual(
            relax["solver"]["contact_boundary_minority_electron_relaxation_contact_side"],
            "n_contact_only",
        )
        self.assertEqual(len(summary["prepared"]), 5)

    def test_pn2d_bv_poisson_boundary_charge_diagnostic_writes_doping_contact_charge(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_poisson_charge_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_-1v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 0.5, "y": 0.0},
                    {"id": 2, "x": 1.0, "y": 0.0},
                    {"id": 3, "x": 1.5, "y": 0.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                    {"id": 1, "node_ids": [1, 3, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "node_ids": [0]},
                    {"id": 1, "name": "Cathode", "node_ids": [3]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 1.0e17],
                [2, 1.0e17, 0.0],
                [3, 1.0e17, 0.0],
            ])
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
poisson charge
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 4 float
0 0 0
5e-7 0 0
1e-6 0 0
1.5e-6 0 0
POINT_DATA 4
SCALARS Potential float 1
LOOKUP_TABLE default
-1.4 -0.4 0.1 0.4
SCALARS ElectronQuasiFermi float 1
LOOKUP_TABLE default
-1.0 -0.5 0.0 0.0
SCALARS HoleQuasiFermi float 1
LOOKUP_TABLE default
-1.0 -0.5 0.0 0.0
SCALARS Electrons float 1
LOOKUP_TABLE default
1e9 2e12 3e12 4e12
SCALARS Holes float 1
LOOKUP_TABLE default
4e12 3e12 2e12 1e9
""".lstrip()
            )
            for name, values in {
                "ElectrostaticPotential": [-1.39, -0.4, 0.1, 0.39],
                "eQuasiFermiPotential": [-1.0, -0.5, 0.0, 0.0],
                "hQuasiFermiPotential": [-1.0, -0.5, 0.0, 0.0],
                "eDensity": [1.0e3, 2.0e6, 3.0e6, 4.0e6],
                "hDensity": [4.0e6, 3.0e6, 2.0e6, 1.0e3],
                "DonorConcentration": [0.0, 1.0e17, 1.0e17, 1.0e17],
                "AcceptorConcentration": [1.0e17, 1.0e17, 0.0, 0.0],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            self._write_csv(
                fields / "ContactExternalVoltage_region1.csv",
                ["node_id", "component0"],
                [[0, -1.0]],
            )
            self._write_csv(
                fields / "ContactExternalVoltage_region2.csv",
                ["node_id", "component0"],
                [[3, 0.0]],
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_poisson_boundary_charge.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "-1",
                    "--bands", "left:0:0,junction:0.5:0.5,right:1.5:1.5",
                ],
                check=True,
                cwd=REPO,
            )

            doping_rows = self._read_csv(out / "doping_distribution_summary.csv")
            contact_rows = self._read_csv(out / "contact_boundary_summary.csv")
            charge_rows = self._read_csv(out / "charge_balance_bands.csv")
            summary = json.loads((out / "poisson_boundary_charge_summary.json").read_text())

        vela_groups = [row for row in doping_rows if row["source"] == "vela"]
        self.assertEqual(len(vela_groups), 3)
        self.assertIn("expected_psi_builtin_V", contact_rows[0])
        anode = next(row for row in contact_rows if row["contact"] == "Anode")
        self.assertLess(float(anode["expected_psi_builtin_V"]), 0.0)
        right = next(row for row in charge_rows if row["band"] == "right")
        self.assertAlmostEqual(float(right["vela_electron_median_cm3"]), 4.0e6)
        self.assertAlmostEqual(float(right["sentaurus_electron_median_cm3"]), 4.0e6)
        self.assertAlmostEqual(float(right["delta_net_doping_median_cm3"]), 0.0)
        self.assertEqual(summary["charge_rows"], 3)

    def test_pn2d_bv_poisson_flux_balance_diagnostic_writes_top_nodes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_poisson_flux_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_-1v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "node_ids": [0]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
            ])
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
poisson flux
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
0.0 0.1 0.0
SCALARS Electrons float 1
LOOKUP_TABLE default
1e10 2e10 3e10
SCALARS Holes float 1
LOOKUP_TABLE default
3e10 2e10 1e10
""".lstrip()
            )
            for name, values in {
                "ElectrostaticPotential": [0.0, 0.05, 0.0],
                "eDensity": [1.0e4, 2.0e4, 3.0e4],
                "hDensity": [3.0e4, 2.0e4, 1.0e4],
                "DonorConcentration": [0.0, 1.0e17, 0.0],
                "AcceptorConcentration": [1.0e17, 0.0, 1.0e17],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_poisson_flux_balance.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "-1",
                    "--bands", "left:0:0,right:1:1",
                    "--top", "4",
                ],
                check=True,
                cwd=REPO,
            )

            top_rows = self._read_csv(out / "poisson_flux_balance_top_nodes.csv")
            band_rows = self._read_csv(out / "poisson_flux_balance_bands.csv")
            compare_rows = self._read_csv(out / "poisson_flux_balance_band_compare.csv")
            summary = json.loads((out / "poisson_flux_balance_summary.json").read_text())

        self.assertEqual(len(top_rows), 4)
        self.assertIn("flux_term_C_per_m", top_rows[0])
        self.assertIn("charge_term_C_per_m", top_rows[0])
        self.assertFalse(any(row["contact_names"] == "Anode" for row in top_rows))
        vela_right = next(
            row for row in band_rows
            if row["state"] == "vela" and row["band"] == "right"
        )
        self.assertEqual(vela_right["node_count"], "1")
        self.assertTrue(compare_rows)
        self.assertEqual(summary["compare_rows"], 2)

    def test_pn2d_bv_fixed_charge_rhs_diagnostic_writes_integrated_doping(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_fixed_charge_rhs_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            out = root / "reports"
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 2.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 2.0},
                    {"id": 3, "x": 2.0, "y": 2.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                    {"id": 1, "node_ids": [1, 3, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "node_ids": [0]},
                    {"id": 1, "name": "Cathode", "node_ids": [1]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 1.0e17, 0.0],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_fixed_charge_rhs.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--out-dir", str(out),
                    "--bands", "left:0:0.5,right:1.5:2.0",
                ],
                check=True,
                cwd=REPO,
            )

            node_rows = self._read_csv(out / "fixed_charge_rhs_nodes.csv")
            band_rows = self._read_csv(out / "fixed_charge_rhs_bands.csv")
            profile_rows = self._read_csv(out / "fixed_charge_rhs_profile.csv")
            integrated = self._read_csv(out / "integrated_doping.csv")
            mixed_volume = self._read_csv(out / "mixed_volume_doping.csv")
            summary = json.loads((out / "fixed_charge_rhs_summary.json").read_text())

        self.assertEqual(len(node_rows), 4)
        self.assertEqual(len(integrated), 4)
        self.assertEqual(len(mixed_volume), 4)
        self.assertGreaterEqual(len(profile_rows), 2)
        right_node = next(row for row in node_rows if row["node_id"] == "1")
        self.assertLess(float(right_node["integrated_to_nodal_ratio"]), 1.0)
        self.assertIn("mixed_to_barycentric_volume_ratio", right_node)
        self.assertGreater(float(right_node["abs_mixed_volume_delta_fixed_charge_C_per_m"]), 0.0)
        self.assertGreater(float(right_node["abs_delta_fixed_charge_C_per_m"]), 0.0)
        right_band = next(row for row in band_rows if row["band"] == "right")
        self.assertEqual(right_band["node_count"], "2")
        self.assertIn("median_mixed_volume_to_nodal_ratio", right_band)
        self.assertEqual(summary["node_rows"], 4)

    def test_pn2d_bv_junction_geometry_diagnostic_writes_box_terms(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_junction_geometry_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            inventory = root / "inventory_mesh.json"
            doping = root / "doping.csv"
            out = root / "reports"
            mesh_data = {
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
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "node_ids": [0]},
                ],
            }
            mesh.write_text(json.dumps(mesh_data) + "\n")
            inventory.write_text(json.dumps({
                "geometry": {
                    "vertices": [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
                    "regions": [{
                        "index": 0,
                        "name": "R.Si",
                        "material": "Silicon",
                        "type": 0,
                        "triangles": [[0, 1, 2], [1, 3, 2]],
                        "edges": [],
                        "points": [],
                    }],
                }
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 0.0, 1.0e17],
                [3, 1.0e17, 0.0],
            ])

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_junction_geometry.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--tdr-inventory-mesh", str(inventory),
                    "--out-dir", str(out),
                    "--focus-nodes", "1",
                    "--x-window", "0.75:1.0",
                    "--radius-um", "0.1",
                    "--bands", "right:0.75:1.0",
                ],
                check=True,
                cwd=REPO,
            )

            node_rows = self._read_csv(out / "junction_geometry_nodes.csv")
            edge_rows = self._read_csv(out / "junction_geometry_edges.csv")
            cell_rows = self._read_csv(out / "junction_geometry_cells.csv")
            band_rows = self._read_csv(out / "junction_geometry_bands.csv")
            summary = json.loads((out / "junction_geometry_summary.json").read_text())

        focus = next(row for row in node_rows if row["node_id"] == "1")
        self.assertIn("focus", focus["selected_reason"])
        self.assertGreater(float(focus["control_volume_m2"]), 0.0)
        self.assertTrue(edge_rows)
        self.assertIn("couple_over_length", edge_rows[0])
        self.assertTrue(cell_rows)
        self.assertEqual(band_rows[0]["band"], "right")
        self.assertTrue(summary["inventory_comparison"]["triangle_sets_match"])
        self.assertAlmostEqual(summary["inventory_comparison"]["max_coordinate_delta_um"], 0.0)

    def test_pn2d_bv_poisson_reconstruction_diagnostic_solves_frozen_states(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_bv_poisson_recon_") as tmp:
            root = Path(tmp)
            mesh = root / "mesh.json"
            doping = root / "doping.csv"
            vela = root / "vela"
            sentaurus = root / "sentaurus" / "sentaurus_-1v"
            fields = sentaurus / "fields"
            out = root / "reports"
            vela.mkdir()
            fields.mkdir(parents=True)
            mesh.write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "R.Si", "material": "Silicon"}],
                "contacts": [
                    {"id": 0, "name": "Anode", "node_ids": [0]},
                    {"id": 1, "name": "Cathode", "node_ids": [1]},
                ],
            }) + "\n")
            self._write_csv(doping, ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 0.0, 1.0e17],
                [1, 1.0e17, 0.0],
                [2, 1.0e17, 1.0e17],
            ])
            (vela / "mini_0000_-1V.vtk").write_text(
                """
# vtk DataFile Version 2.0
poisson reconstruction
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 float
0 0 0
1e-6 0 0
0 1e-6 0
CELLS 1 4
3 0 1 2
CELL_TYPES 1
5
POINT_DATA 3
SCALARS Potential float 1
LOOKUP_TABLE default
-1.4 0.4 -0.5
SCALARS Electrons float 1
LOOKUP_TABLE default
1e10 2e10 3e10
SCALARS Holes float 1
LOOKUP_TABLE default
3e10 2e10 1e10
""".lstrip()
            )
            for name, values in {
                "ElectrostaticPotential": [-1.39, 0.39, -0.45],
                "eDensity": [1.0e4, 2.0e4, 3.0e4],
                "hDensity": [3.0e4, 2.0e4, 1.0e4],
                "DonorConcentration": [0.0, 1.0e17, 1.0e17],
                "AcceptorConcentration": [1.0e17, 0.0, 1.0e17],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "diagnose_pn2d_bv_poisson_reconstruction.py"),
                    "--mesh", str(mesh),
                    "--doping-csv", str(doping),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--vela-vtk-root", str(vela),
                    "--out-dir", str(out),
                    "--biases", "-1",
                    "--focus-nodes", "0,1,2",
                    "--bc-sources", "vela_expected,sentaurus_state",
                    "--charge-sources", "depletion,vela_frozen,sentaurus_frozen",
                    "--bands", "all:0:1",
                ],
                check=True,
                cwd=REPO,
            )

            band_rows = self._read_csv(out / "poisson_reconstruction_bands.csv")
            node_rows = self._read_csv(out / "poisson_reconstruction_focus_nodes.csv")
            contact_rows = self._read_csv(out / "poisson_reconstruction_contacts.csv")
            summary = json.loads((out / "poisson_reconstruction_summary.json").read_text())

        self.assertEqual(len(band_rows), 6)
        self.assertEqual(len(node_rows), 18)
        self.assertEqual(len(contact_rows), 4)
        self.assertIn("median_reconstructed_minus_vela_V", band_rows[0])
        self.assertIn("reconstructed_minus_sentaurus_V", node_rows[0])
        self.assertEqual(summary["band_rows"], 6)

    def test_pn2d_sentaurus2018_iv_cmd_writes_multibias_debug_snapshots(self) -> None:
        path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_iv_sdevice.cmd"
        text = path.read_text()

        self.assertIn("eMobility", text)
        self.assertIn("hMobility", text)
        self.assertIn("MaxStep=0.025", text)
        self.assertRegex(text, r'Goal\s*\{\s*Name="Anode"\s*Voltage=10\.0\s*\}')
        self.assertIn("0.05 V spacing over the 0-10 V normalized sweep", text)
        self.assertIn(
            'Plot(FilePrefix="pn2d_iv_multibias" Time=(Range=(0 1) Intervals=200) NoOverWrite)',
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
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
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

    def test_runner_writes_newton_step_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_newton_step_probe_") as tmp:
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
            output = root / "step.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_step_probe",
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
        self.assertIn("step_norm", status)
        self.assertIn("trial_block_residuals", status)
        self.assertIn("delta_psi_V", rows[0])
        self.assertIn("delta_phin_V", rows[0])
        self.assertIn("delta_psi_minus_phin_V", rows[0])
        self.assertIn("trial_electron_density_m3", rows[0])

    def test_runner_writes_newton_jvp_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_newton_jvp_probe_") as tmp:
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
            output = root / "jvp.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_jvp_probe",
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
                "directions": [
                    {
                        "name": "middle_psi_minus_phin",
                        "mode": "psi_minus_phin",
                        "amplitude_V": 1.0e-6,
                        "exclude_contacts": False,
                    }
                ],
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

        self.assertEqual(status["direction_count"], 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["direction"], "middle_psi_minus_phin")
        self.assertIn("analytic_norm", rows[0])
        self.assertIn("finite_difference_norm", rows[0])
        self.assertLess(float(rows[0]["relative_error"]), 1.0e-4)
        self.assertGreater(int(rows[0]["selected_nodes"]), 0)

    def test_runner_writes_newton_block_step_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_newton_block_step_probe_") as tmp:
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
            output = root / "block_step.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_block_step_probe",
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
                "block_modes": ["poisson_only", "carrier_only"],
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

        self.assertEqual(status["block_step_count"], 2)
        self.assertEqual(len(rows), 8)
        self.assertEqual({row["mode"] for row in rows}, {"poisson_only", "carrier_only"})
        poisson = next(row for row in rows if row["mode"] == "poisson_only")
        carrier = next(row for row in rows if row["mode"] == "carrier_only")
        self.assertIn("delta_psi_minus_phin_V", poisson)
        self.assertAlmostEqual(float(poisson["delta_phin_V"]), 0.0)
        self.assertAlmostEqual(float(poisson["delta_phip_V"]), 0.0)
        self.assertAlmostEqual(float(carrier["delta_psi_V"]), 0.0)

    def test_runner_writes_newton_carrier_row_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_newton_carrier_row_probe_") as tmp:
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
            output = root / "carrier_rows.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_carrier_row_probe",
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
                    "max_update": 1.0,
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

        self.assertEqual(status["row_count"], 4)
        self.assertEqual(len(rows), 4)
        self.assertIn("electron_diagonal", rows[0])
        self.assertIn("electron_row_abs_sum", rows[0])
        self.assertIn("raw_delta_phin_V", rows[0])
        self.assertIn("capped_delta_phin_V", rows[0])
        self.assertLessEqual(abs(float(rows[0]["capped_delta_phin_V"])), 0.026)

    def test_runner_writes_newton_carrier_term_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_newton_carrier_term_probe_") as tmp:
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
            output = root / "carrier_terms.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_carrier_term_probe",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "output_csv": str(output),
                "state_fields_dir": "fields",
                "carrier_term_probe": {
                    "electron_impact_scale": 2.0,
                    "hole_impact_scale": 3.0,
                },
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

        self.assertEqual(status["row_count"], 4)
        self.assertEqual(len(rows), 4)
        self.assertIn("electron_flux", rows[0])
        self.assertIn("electron_recombination", rows[0])
        self.assertIn("electron_impact", rows[0])
        self.assertIn("electron_adjusted_impact", rows[0])
        self.assertIn("electron_adjusted_residual", rows[0])
        self.assertIn("hole_adjusted_impact", rows[0])
        self.assertIn("hole_adjusted_residual", rows[0])
        self.assertIn("impact_electron_source", rows[0])
        self.assertIn("impact_hole_source", rows[0])
        self.assertIn("impact_combined_source", rows[0])
        self.assertIn("electron_residual", rows[0])
        self.assertIn("hole_flux", rows[0])
        self.assertEqual(status["carrier_term_probe"]["electron_impact_scale"], 2.0)
        self.assertEqual(status["carrier_term_probe"]["hole_impact_scale"], 3.0)
        self.assertIn("adjusted_block_residuals", status)

    def test_runner_writes_edge_mobility_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_edge_mobility_probe_") as tmp:
            root = Path(tmp)
            fields = root / "fields"
            fields.mkdir()
            self._write_csv(root / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], [
                [0, 1.0e17, 0.0],
                [1, 1.0e17, 0.0],
                [2, 1.0e17, 0.0],
            ])
            (root / "mesh.json").write_text(json.dumps({
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 2], "region_id": 0},
                ],
                "regions": [{"id": 0, "name": "Si", "material": "Si", "cell_ids": [0]}],
                "contacts": [],
            }) + "\n")
            for name, values in {
                "ElectrostaticPotential": [0.0, 1.0, 0.0],
                "eQuasiFermiPotential": [0.0, 2.0, 0.0],
                "hQuasiFermiPotential": [0.0, -1.0, 0.0],
            }.items():
                self._write_csv(
                    fields / f"{name}_region0.csv",
                    ["node_id", "component0"],
                    [[idx, value] for idx, value in enumerate(values)],
                )
            output = root / "edge_mobility.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "edge_mobility_probe",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "output_csv": str(output),
                "state_fields_dir": "fields",
                "scaling": {"mode": "unit_scaling"},
                "doping": [{"region": "Si", "donors": 1.0e17, "acceptors": 0.0}],
                "contacts": [],
                "solver": {
                    "mobility": {
                        "model": "caughey_thomas_field",
                        "high_field_driving_force": "quasi_fermi_gradient",
                        "electron_saturation_velocity_m_s": 1.0e5,
                        "hole_saturation_velocity_m_s": 1.0e5,
                    },
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

        self.assertEqual(status["edge_count"], 3)
        self.assertEqual(len(rows), 3)
        focus = next(row for row in rows if row["edge_id"] == "0")
        self.assertAlmostEqual(float(focus["electric_field_V_m"]), 1.0e6)
        self.assertAlmostEqual(float(focus["electron_mobility_field_V_m"]), 2.0e6)
        self.assertAlmostEqual(float(focus["hole_mobility_field_V_m"]), 1.0e6)
        self.assertGreater(float(focus["electron_low_field_mobility_m2_V_s"]), 0.0)
        self.assertLess(
            float(focus["electron_final_mobility_m2_V_s"]),
            float(focus["electron_low_field_mobility_m2_V_s"]),
        )
        self.assertIn("electron_mobility_limiter", focus)

    def test_runner_writes_newton_regularized_carrier_step_probe_for_external_state(self) -> None:
        exe_name = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
        runner = REPO / "build-release" / exe_name
        if not runner.is_file():
            runner = REPO / "build" / exe_name
        if not runner.is_file():
            self.skipTest(f"missing built runner: {runner}")

        with tempfile.TemporaryDirectory(prefix="vela_regularized_carrier_step_probe_") as tmp:
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
            output = root / "regularized_carrier_step.csv"
            config = root / "probe.json"
            config.write_text(json.dumps({
                "simulation_type": "newton_regularized_carrier_step_probe",
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
                "regularization_scales": [0.0, 10.0],
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

        self.assertEqual(status["regularized_step_count"], 2)
        self.assertEqual(len(rows), 8)
        self.assertEqual({row["regularization_scale"] for row in rows}, {"0", "10"})
        self.assertIn("delta_psi_minus_phin_V", rows[0])
        self.assertIn("raw_step_norm", rows[0])

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
