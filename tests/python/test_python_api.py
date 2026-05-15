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

        cv_cfg = self.write_sweep_config(mode="cv_quasistatic", output_name="cv.csv")
        cv_points = vela.run_cv_curve(str(cv_cfg))
        self.assertEqual(cv_points[0]["curve_type"], "cv_quasistatic")
        self.assertIn("capacitance", cv_points[1])
        self.assertIn("terminal_charge", cv_points[1])
        self.assertTrue(Path(cv_points[0]["output_csv"]).exists())

        bv_cfg = self.write_sweep_config(mode="bv_reverse", output_name="bv.csv")
        bv_points = vela.run_bv_curve(str(bv_cfg))
        self.assertEqual(bv_points[0]["curve_type"], "bv_reverse")
        self.assertIn("breakdown_voltage", bv_points[0])
        self.assertTrue(any(point["breakdown_detected"] for point in bv_points))
        self.assertTrue(Path(bv_points[0]["output_csv"]).exists())
        self.assertTrue(Path(bv_points[0]["output_vtk"]).exists())


if __name__ == "__main__":
    unittest.main()
