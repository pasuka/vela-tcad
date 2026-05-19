import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import vela


class PythonApiTest(unittest.TestCase):
    def setUp(self):
        source_dir_env = os.environ.get("VELA_SOURCE_DIR")
        if source_dir_env is None:
            raise unittest.SkipTest("VELA_SOURCE_DIR must point to the repository root")
        source_dir = Path(source_dir_env)
        self.tmpdir = Path(tempfile.mkdtemp(prefix="vela_python_api_"))
        self.mesh_file = self.tmpdir / "mesh.json"
        shutil.copyfile(source_dir / "examples" / "pn_diode" / "mesh.json", self.mesh_file)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def write_poisson_config(self):
        cfg = {
            "mesh_file": str(self.mesh_file),
            "output_vtk": str(self.tmpdir / "poisson.vtk"),
            "doping": [
                {"region": "n_region", "donors": 1.0e23, "acceptors": 0.0},
                {"region": "p_region", "donors": 0.0, "acceptors": 1.0e23},
            ],
            "contacts": [
                {"name": "anode", "bias": 0.0},
                {"name": "cathode", "bias": 0.0},
            ],
        }
        path = self.tmpdir / "poisson.json"
        path.write_text(json.dumps(cfg), encoding="utf-8")
        return path

    def write_sweep_config(self, mode="iv", output_name="iv.csv"):
        cfg = {
            "mesh_file": str(self.mesh_file),
            "output_csv": str(self.tmpdir / output_name),
            "doping": [
                {"region": "n_region", "donors": 1.0e23, "acceptors": 0.0},
                {"region": "p_region", "donors": 0.0, "acceptors": 1.0e23},
            ],
            "contacts": [
                {"name": "anode", "bias": 0.0},
                {"name": "cathode", "bias": 0.0},
            ],
            "solver": {
                "max_iter": 80,
                "reltol": 1.0e-5,
                "damping_psi": 0.5,
            },
            "sweep": {
                "mode": mode,
                "contact": "anode",
                "start": 0.0,
                "stop": 0.5,
                "step": 0.25,
                "current_contact": "anode",
                "write_vtk": True,
                "vtk_prefix": str(self.tmpdir / f"pn_{mode}_sweep"),
                "terminal_charge": {
                    "contact": "anode",
                    "regions": ["p_region", "n_region"],
                    "per_meter": True,
                },
                "breakdown": {
                    "max_electric_field_V_per_m": 1.0,
                    "current_jump_ratio": 100.0,
                    "non_convergence": True,
                },
            },
        }
        path = self.tmpdir / f"sweep_{mode}.json"
        path.write_text(json.dumps(cfg), encoding="utf-8")
        return path

    def relative_sweep_config(self):
        (self.tmpdir / "outputs").mkdir(exist_ok=True)
        return {
            "mesh_file": "mesh.json",
            "output_csv": "outputs/mapping_iv.csv",
            "doping": [
                {"region": "n_region", "donors": 1.0e23, "acceptors": 0.0},
                {"region": "p_region", "donors": 0.0, "acceptors": 1.0e23},
            ],
            "contacts": [
                {"name": "anode", "bias": 0.0},
                {"name": "cathode", "bias": 0.0},
            ],
            "solver": {
                "max_iter": 80,
                "reltol": 1.0e-5,
                "damping_psi": 0.5,
            },
            "sweep": {
                "mode": "iv",
                "contact": "anode",
                "start": 0.0,
                "stop": 0.25,
                "step": 0.25,
                "current_contact": "anode",
                "write_vtk": True,
                "vtk_prefix": "outputs/mapping_iv_sweep",
            },
        }

    def test_import_and_run_api(self):
        mesh = vela.load_mesh(str(self.mesh_file))
        self.assertEqual(mesh.num_nodes(), 4)

        poisson_cfg = self.write_poisson_config()
        psi = vela.run_poisson(str(poisson_cfg))
        self.assertEqual(len(psi), 4)
        self.assertTrue((self.tmpdir / "poisson.vtk").exists())

        exported_vtk = self.tmpdir / "python_export.vtk"
        vela.write_vtk(str(exported_vtk))
        self.assertTrue(exported_vtk.exists())

        sweep_cfg = self.write_sweep_config()
        points = vela.run_iv_curve(str(sweep_cfg))
        self.assertEqual(len(points), 3)
        self.assertTrue(all(point["converged"] for point in points))
        self.assertEqual(points[0]["curve_type"], "iv")
        self.assertEqual(points[0]["bias_contact"], "anode")
        self.assertEqual(points[0]["current_contact"], "anode")
        self.assertEqual(points[0]["attempted_step"], 0.0)
        self.assertEqual(points[0]["accepted_step"], 0.0)
        self.assertEqual(points[0]["retry_count"], 0)
        self.assertEqual(points[0]["convergence_diagnostics"]["retry_count"], 0)
        self.assertIn("validation_diagnostics", points[0])
        self.assertEqual(
            points[0]["convergence_diagnostics"]["validation_diagnostics"],
            points[0]["validation_diagnostics"],
        )
        self.assertIn("total_current", points[1])
        self.assertAlmostEqual(points[1]["attempted_step"], 0.25)
        self.assertAlmostEqual(points[1]["accepted_step"], 0.25)
        self.assertEqual(points[1]["retry_count"], 0)
        self.assertEqual(points[0]["output_csv"], str(self.tmpdir / "iv.csv"))
        self.assertEqual(points[0]["output_vtk"], str(self.tmpdir / "pn_iv_sweep_0000_0V.vtk"))
        self.assertTrue(Path(points[0]["output_csv"]).exists())
        self.assertIn(points[0]["output_csv"], points[0]["output_files"])
        self.assertIn(points[0]["output_vtk"], points[0]["output_files"])
        self.assertTrue(Path(points[0]["output_vtk"]).exists())

        validation_failure_cfg = self.write_sweep_config(output_name="validation_failure.csv")
        validation_failure_data = json.loads(validation_failure_cfg.read_text(encoding="utf-8"))
        validation_failure_data["validation"] = {
            "enforce_minimum_carrier_density": True,
            "minimum_carrier_density": 1.0e30,
        }
        validation_failure_cfg.write_text(json.dumps(validation_failure_data), encoding="utf-8")
        validation_failure_points = vela.run_iv_curve(str(validation_failure_cfg))
        self.assertEqual(validation_failure_points[0]["failure_reason"], "validation_failed")
        self.assertIn("below minimum", validation_failure_points[0]["validation_diagnostics"])
        self.assertEqual(
            validation_failure_points[0]["convergence_diagnostics"]["validation_diagnostics"],
            validation_failure_points[0]["validation_diagnostics"],
        )

        cv_cfg = self.write_sweep_config(mode="cv_quasistatic", output_name="cv.csv")
        cv_points = vela.run_cv_curve(str(cv_cfg))
        self.assertEqual(cv_points[0]["curve_type"], "cv_quasistatic")
        self.assertIn("capacitance", cv_points[1])
        self.assertIn("terminal_charge", cv_points[1])
        self.assertTrue(Path(cv_points[0]["output_csv"]).exists())

        multi_cv_cfg = self.write_sweep_config(mode="cv_quasistatic", output_name="cv_multi.csv")
        multi_cv_data = json.loads(multi_cv_cfg.read_text(encoding="utf-8"))
        multi_cv_data["sweep"]["terminal_charges"] = [
            {"name": "gate", "contact": "anode", "regions": ["p_region"], "per_meter": True},
            {"name": "drain", "contact": "cathode", "regions": ["n_region"], "per_meter": True},
        ]
        multi_cv_cfg.write_text(json.dumps(multi_cv_data), encoding="utf-8")
        multi_cv_points = vela.run_cv_curve(str(multi_cv_cfg))
        self.assertIn("charge_gate_C_per_m", multi_cv_points[1])
        self.assertIn("capacitance_Canode_gate_F_per_m", multi_cv_points[1])
        self.assertIn("charge_drain_C_per_m", multi_cv_points[1])
        self.assertIn("capacitance_Canode_drain_F_per_m", multi_cv_points[1])
        self.assertIn("gate", multi_cv_points[1]["terminal_charges"])
        self.assertIn("drain", multi_cv_points[1]["terminal_capacitances"])

        stored_cv_data = json.loads(multi_cv_cfg.read_text(encoding="utf-8"))
        stored_cv_data["sweep"]["stored_charge"] = {"regions": ["p_region", "n_region"], "per_meter": True}
        multi_cv_cfg.write_text(json.dumps(stored_cv_data), encoding="utf-8")
        stored_cv_points = vela.run_cv_curve(str(multi_cv_cfg))
        self.assertIn("stored_charge_C_per_m", stored_cv_points[1])

        bv_cfg = self.write_sweep_config(mode="bv_reverse", output_name="bv.csv")
        bv_points = vela.run_bv_curve(str(bv_cfg))
        self.assertEqual(bv_points[0]["curve_type"], "bv_reverse")
        self.assertIn("breakdown_voltage", bv_points[0])
        self.assertTrue(any(point["breakdown_detected"] for point in bv_points))
        self.assertTrue(Path(bv_points[0]["output_csv"]).exists())
        self.assertTrue(Path(bv_points[0]["output_vtk"]).exists())

        old_cwd = Path.cwd()
        try:
            os.chdir(self.tmpdir)
            mapping_points = vela.run_iv_curve(self.relative_sweep_config())
        finally:
            os.chdir(old_cwd)
        self.assertEqual(len(mapping_points), 2)
        self.assertEqual(mapping_points[0]["output_csv"], str(self.tmpdir / "outputs" / "mapping_iv.csv"))
        self.assertEqual(mapping_points[0]["output_vtk"], str(self.tmpdir / "outputs" / "mapping_iv_sweep_0000_0V.vtk"))
        self.assertTrue(Path(mapping_points[0]["output_csv"]).exists())
        self.assertTrue(Path(mapping_points[0]["output_vtk"]).exists())


if __name__ == "__main__":
    unittest.main()
