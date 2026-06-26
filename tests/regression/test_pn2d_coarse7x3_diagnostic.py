#!/usr/bin/env python3
"""Regression coverage for the PN2D coarse Sentaurus/Vela BV diagnostic."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
COARSE_ROOT = REPO / "reference_tcad" / "pn2d_sentaurus2018_coarse7x3"
COARSE_SOURCE = COARSE_ROOT / "source"
REQUIRED_PLOT_FIELDS = {
    "Potential",
    "ElectricField",
    "eDensity",
    "hDensity",
    "eQuasiFermi",
    "hQuasiFermi",
    "eCurrent",
    "hCurrent",
    "TotalCurrent",
    "Doping",
    "DonorConcentration",
    "AcceptorConcentration",
    "SRHRecombination",
    "AvalancheGeneration",
    "eMobility",
    "hMobility",
    "eVelocity",
    "hVelocity",
    "eAlphaAvalanche",
    "hAlphaAvalanche",
    "eIonIntegral",
    "hIonIntegral",
    "MeanIonIntegral",
}


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Pn2dCoarse7x3DiagnosticTest(unittest.TestCase):
    def test_coarse_case_is_independent_and_exports_required_sentaurus_fields(self) -> None:
        self.assertTrue((COARSE_SOURCE / "pn2d_sde.cmd").is_file())
        self.assertTrue((COARSE_SOURCE / "pn2d_bv_sdevice.cmd").is_file())
        self.assertTrue((COARSE_SOURCE / "models.par").is_file())
        self.assertFalse(str(COARSE_ROOT).endswith("pn2d_sentaurus2018"))

        sde = (COARSE_SOURCE / "pn2d_sde.cmd").read_text(encoding="utf-8")
        self.assertIn("(define L 2.0)", sde)
        self.assertIn("(define H 0.5)", sde)
        self.assertIn("(define XJ 1.0)", sde)
        self.assertIn("0.3333333333333333 0.25", sde)
        self.assertNotIn("Junction.Mesh", sde)

        deck = (COARSE_SOURCE / "pn2d_bv_sdevice.cmd").read_text(encoding="utf-8")
        missing = sorted(field for field in REQUIRED_PLOT_FIELDS if field not in deck)
        self.assertEqual(missing, [])

    def test_coarse_reference_config_does_not_overwrite_fine_reference(self) -> None:
        config_path = COARSE_ROOT / "pn2d_sentaurus2018_coarse7x3_reference.json"
        self.assertTrue(config_path.is_file())
        config = json.loads(config_path.read_text(encoding="utf-8"))
        self.assertEqual(config["case"], "pn2d_sentaurus2018_coarse7x3")
        self.assertEqual(config["mesh_tdr"], "pn2d_msh.tdr")
        self.assertEqual(config["sde_cmd"], "pn2d_sde.cmd")
        self.assertEqual(len(config["simulations"]), 1)
        sim = config["simulations"][0]
        self.assertEqual(sim["name"], "bv")
        self.assertEqual(sim["vela_materials_file"], "pn2d_sentaurus2018_iv_materials.json")
        self.assertFalse(sim.get("execute", True))

    def test_aligned_fixed_bias_config_uses_exact_bias_points_and_distinct_outputs(self) -> None:
        module = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        with tempfile.TemporaryDirectory(prefix="vela_coarse7x3_aligned_cfg_") as tmp:
            root = Path(tmp)
            base = root / "simulation_bv.json"
            materials = root / "pn2d_sentaurus2018_iv_materials.json"
            (root / "mesh.json").write_text(
                json.dumps({"nodes": [{"id": 0, "x": 0.0, "y": 0.0}], "triangles": [], "contacts": []}),
                encoding="utf-8",
            )
            (root / "doping.csv").write_text("node_id,donors_cm3,acceptors_cm3\n0,0,1e17\n", encoding="utf-8")
            base.write_text(
                json.dumps({"mesh_file": "mesh.json", "node_doping_file": "doping.csv", "solver": {}, "sweep": {}}),
                encoding="utf-8",
            )
            materials.write_text(json.dumps({"materials": [{"name": "Si", "ni": 1.0e10}]}), encoding="utf-8")

            patched = module.write_previous_full20_config(
                base_config=base,
                out_dir=root,
                output_csv_name="coarse_previous_full20_aligned.csv",
                bias_points=[0.0, -1.0, -5.0, -10.0, -16.0, -18.0, -20.0],
                config_name="simulation_coarse_previous_full20_aligned.json",
                vtk_subdir="vtk_aligned",
                diagnostics_suffix="_aligned",
            )

            data = json.loads(patched.read_text(encoding="utf-8"))
            self.assertEqual(patched.name, "simulation_coarse_previous_full20_aligned.json")
            self.assertEqual(data["output_csv"], "coarse_previous_full20_aligned.csv")
            self.assertEqual(data["sweep"]["bias_points"], [0.0, -1.0, -5.0, -10.0, -16.0, -18.0, -20.0])
            self.assertEqual(data["sweep"]["vtk_prefix"], "vtk_aligned/dc_sweep")
            self.assertTrue(data["sweep"]["diagnostics"]["sg_avalanche_edges"]["csv_file"].endswith("sg_avalanche_edges_aligned.csv"))
    def test_previous_full20_patch_requires_materials_and_grad_qf_cell_gradient(self) -> None:
        module = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        with tempfile.TemporaryDirectory(prefix="vela_coarse7x3_patch_") as tmp:
            root = Path(tmp)
            base = root / "simulation_bv.json"
            materials = root / "pn2d_sentaurus2018_iv_materials.json"
            (root / "mesh.json").write_text(
                json.dumps({
                    "nodes": [{"id": 0, "x": 0.0, "y": 0.0}],
                    "triangles": [],
                    "contacts": [],
                }),
                encoding="utf-8",
            )
            (root / "doping.csv").write_text(
                "node_id,donors_cm3,acceptors_cm3\n0,0,1e17\n",
                encoding="utf-8",
            )
            base.write_text(
                json.dumps(
                    {
                        "mesh_file": "mesh.json",
                        "node_doping_file": "doping.csv",
                        "output_csv": "old.csv",
                        "solver": {
                            "impact_ionization": {
                                "model": "van_overstraeten",
                                "current_approximation": "density_gradient",
                            }
                        },
                        "sweep": {"write_vtk": False},
                    }
                ),
                encoding="utf-8",
            )
            materials.write_text(json.dumps({"materials": [{"name": "Si", "ni": 1.0e10}]}), encoding="utf-8")

            patched = module.write_previous_full20_config(
                base_config=base,
                out_dir=root,
                output_csv_name="coarse_previous_full20.csv",
            )

            data = json.loads(patched.read_text(encoding="utf-8"))
            self.assertEqual(data["output_csv"], "coarse_previous_full20.csv")
            self.assertEqual(Path(data["materials_file"]).name, materials.name)
            self.assertTrue(data["sweep"]["write_vtk"])
            impact = data["solver"]["impact_ionization"]
            self.assertEqual(impact["current_approximation"], "grad_qf")
            self.assertEqual(impact["quasi_fermi_gradient_discretization"], "cell_gradient")
            self.assertEqual(data["solver"]["quasi_fermi_update_limit_V"], 0.1)
            self.assertEqual(data["solver"]["max_update"], 0)

    def test_multibias_compare_knows_sentaurus_alpha_velocity_and_ion_integrals(self) -> None:
        module = load_module(REPO / "scripts" / "compare_pn2d_bv_multibias_fields.py")
        for quantity in [
            "electron_velocity",
            "hole_velocity",
            "electron_alpha_avalanche",
            "hole_alpha_avalanche",
            "electron_ion_integral",
            "hole_ion_integral",
            "mean_ion_integral",
        ]:
            self.assertIn(quantity, module.FIELD_SPECS)
            self.assertTrue(module.FIELD_SPECS[quantity]["sentaurus"])
        self.assertEqual(module.FIELD_SPECS["electron_alpha_avalanche"]["sentaurus_scale"], 100.0)
        self.assertEqual(module.FIELD_SPECS["hole_alpha_avalanche"]["sentaurus_scale"], 100.0)

        coarse = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        self.assertEqual(coarse.NODE_FIELD_MAP["electron_alpha_avalanche"]["sentaurus_scale"], 100.0)
        self.assertEqual(coarse.NODE_FIELD_MAP["hole_alpha_avalanche"]["sentaurus_scale"], 100.0)

    def test_curve_plot_writes_requested_sentaurus_reference_alias(self) -> None:
        module = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        with tempfile.TemporaryDirectory(prefix="vela_coarse7x3_curve_") as tmp:
            root = Path(tmp)
            ref = root / "importer_reference.csv"
            cand = root / "coarse_previous_full20.csv"
            ref.write_text("bias_V,current_total\n0,1e-12\n-1,2e-12\n", encoding="utf-8")
            cand.write_text("bias_V,current_total_A_per_um\n0,1.1e-12\n-1,2.1e-12\n", encoding="utf-8")

            os.environ["MPLCONFIGDIR"] = str(root / "mplconfig")
            result = module.write_curve_plot(ref, cand, root)

            alias = root / "sentaurus_coarse_bv_reference.csv"
            self.assertTrue(alias.is_file())
            self.assertEqual(alias.read_text(encoding="utf-8"), ref.read_text(encoding="utf-8"))
            self.assertEqual(Path(result["sentaurus_coarse_bv_reference_csv"]), alias)
    def test_node_field_compare_uses_nearest_vela_vtk_snapshot(self) -> None:
        module = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        with tempfile.TemporaryDirectory(prefix="vela_coarse7x3_nearest_vtk_") as tmp:
            root = Path(tmp)
            sent = root / "sentaurus_-10v"
            (sent / "fields").mkdir(parents=True)
            (sent / "nodes.csv").write_text("id,x_um,y_um\n0,0,0\n", encoding="utf-8")
            (sent / "fields" / "Potential_region0.csv").write_text(
                "node_id,value\n0,10\n", encoding="utf-8"
            )
            vtk_far = root / "dc_sweep_0V.vtk"
            vtk_near = root / "dc_sweep_-10.25V.vtk"
            vtk_far.write_text("SCALARS Potential double 1\nLOOKUP_TABLE default\n1\n", encoding="utf-8")
            vtk_near.write_text("SCALARS Potential double 1\nLOOKUP_TABLE default\n11\n", encoding="utf-8")

            rows = module.node_field_compare(
                sentaurus_exports={-10.0: sent},
                vela_vtks={0.0: vtk_far, -10.25: vtk_near},
                biases=[-10.0],
            )
            potential = next(row for row in rows if row["quantity"] == "potential")
            self.assertEqual(potential["vela_bias_V"], -10.25)
            self.assertEqual(potential["vela_value_scaled_to_sentaurus_units"], 11.0)
    def test_node_field_compare_can_require_exact_vela_bias(self) -> None:
        module = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        with tempfile.TemporaryDirectory(prefix="vela_coarse7x3_exact_vtk_") as tmp:
            root = Path(tmp)
            sent = root / "sentaurus_-10v"
            (sent / "fields").mkdir(parents=True)
            (sent / "nodes.csv").write_text("id,x_um,y_um\n0,0,0\n", encoding="utf-8")
            (sent / "fields" / "Potential_region0.csv").write_text("node_id,value\n0,10\n", encoding="utf-8")
            vtk_near = root / "dc_sweep_-10.25V.vtk"
            vtk_near.write_text("SCALARS Potential double 1\nLOOKUP_TABLE default\n11\n", encoding="utf-8")

            rows = module.node_field_compare(
                sentaurus_exports={-10.0: sent},
                vela_vtks={-10.25: vtk_near},
                biases=[-10.0],
                exact_vela_bias=True,
            )
            potential = next(row for row in rows if row["quantity"] == "potential")
            self.assertIsNone(potential["vela_bias_V"])
            self.assertEqual(potential["vela_field"], "")
            self.assertIsNone(potential["vela_value_scaled_to_sentaurus_units"])
    def test_active_edge_top_table_ranks_alpha_flux_and_source_integral(self) -> None:
        module = load_module(REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py")
        rows = [
            {
                "bias_V": "-10.2432", "nearest_sentaurus_bias_V": "-10.0", "edge_id": "1",
                "node0": "0", "node1": "1", "edge_class": "interior",
                "electron_flux_abs": "2.0", "hole_flux_abs": "5.0",
                "electron_alpha_m_inv": "3.0", "hole_alpha_m_inv": "7.0",
                "electron_source_integral": "6.0", "hole_source_integral": "35.0",
                "source_integral": "41.0",
                "sent_e_current": "1.0e-18", "sent_h_current": "2.0e-18",
                "sent_e_alpha": "0.3", "sent_h_alpha": "0.7",
            },
            {
                "bias_V": "-10.2432", "nearest_sentaurus_bias_V": "-10.0", "edge_id": "2",
                "node0": "1", "node1": "2", "edge_class": "interior",
                "electron_flux_abs": "10.0", "hole_flux_abs": "1.0",
                "electron_alpha_m_inv": "11.0", "hole_alpha_m_inv": "13.0",
                "electron_source_integral": "110.0", "hole_source_integral": "13.0",
                "source_integral": "123.0",
                "sent_e_current": "1.0e-18", "sent_h_current": "2.0e-18",
                "sent_e_alpha": "1.1", "sent_h_alpha": "1.3",
            },
        ]

        ranked = module.build_active_edge_top_rows(rows, target_bias=-10.0, nearest_vela_bias=-10.2432, limit=2)

        self.assertEqual([row["edge_id"] for row in ranked], ["2", "1"])
        self.assertAlmostEqual(ranked[0]["electron_alpha_flux"], 110.0)
        self.assertAlmostEqual(ranked[0]["hole_alpha_flux"], 13.0)
        self.assertAlmostEqual(ranked[0]["combined_alpha_flux"], 123.0)
        self.assertAlmostEqual(ranked[0]["source_integral"], 123.0)
        self.assertIn("sent_e_flux_from_J_over_q", ranked[0])
        self.assertIn("electron_flux_over_sentaurus", ranked[0])
    def test_vm_dry_run_accepts_coarse_source_dir(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_coarse7x3_vm_dry_") as tmp:
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "run_sentaurus_vm_reference.py"),
                    "pn2d",
                    "--source-dir",
                    str(COARSE_SOURCE),
                    "--local-output-dir",
                    tmp,
                    "--run-id",
                    "coarse_dry",
                    "--stages",
                    "bv",
                    "--dry-run",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads((Path(tmp) / "coarse_dry" / "sentaurus_vm_run_manifest.json").read_text())
            self.assertIn("sde -e -l pn2d_sde.cmd", manifest["commands"][0])
            self.assertIn("sdevice pn2d_bv_sdevice.cmd", manifest["commands"][1])


if __name__ == "__main__":
    unittest.main()
