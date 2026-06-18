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
        self.assertEqual(solver["bandgap_narrowing"], "old_slotboom")
        self.assertEqual(solver["recombination"], ["srh"])
        self.assertEqual(solver["impact_ionization"]["model"], "van_overstraeten")
        self.assertEqual(
            solver["impact_ionization"]["driving_force"],
            "quasi_fermi_gradient",
        )
        self.assertEqual(solver["impact_ionization"]["generation"], "current_density")
        self.assertEqual(bv["vela_stop"], -0.05)
        self.assertEqual(bv["vela_step"], -0.05)
        self.assertNotIn("candidate_bias_scale", bv["comparison"])
        self.assertEqual(bv["runtime_diagnostic"]["step"], 0.05)
        self.assertEqual(bv["runtime_diagnostic"]["bias_points"], [-0.5, -2.0, -5.0, -10.0, -20.0])

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
        self.assertNotIn("sweep", config)
        anode = next(contact for contact in config["contacts"] if contact["name"] == "Anode")
        self.assertEqual(anode["bias"], -1.0)
        self.assertEqual(summary["node_rows"], 0)
        self.assertEqual(len(summary["states"]), 4)
        shifted = next(item for item in summary["states"] if item["source"] == "hybrid_spsi_shift_vqf")
        self.assertAlmostEqual(shifted["psi_shift_V"], 0.05)

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
