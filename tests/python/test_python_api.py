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

    def write_sweep_config(self):
        cfg = {
            "mesh_file": str(self.mesh_file),
            "output_csv": str(self.tmpdir / "iv.csv"),
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
                "contact": "anode",
                "start": 0.0,
                "stop": 0.5,
                "step": 0.25,
                "current_contact": "anode",
                "write_vtk": True,
                "vtk_prefix": str(self.tmpdir / "pn_sweep"),
            },
        }
        path = self.tmpdir / "sweep.json"
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
        points = vela.run_dc_sweep(str(sweep_cfg))
        self.assertEqual(len(points), 3)
        self.assertTrue(all(point["converged"] for point in points))
        self.assertEqual(points[0]["attempted_step"], 0.0)
        self.assertEqual(points[0]["accepted_step"], 0.0)
        self.assertEqual(points[0]["retry_count"], 0)
        self.assertAlmostEqual(points[1]["attempted_step"], 0.25)
        self.assertAlmostEqual(points[1]["accepted_step"], 0.25)
        self.assertEqual(points[1]["retry_count"], 0)
        self.assertTrue((self.tmpdir / "iv.csv").exists())
        self.assertTrue((self.tmpdir / "pn_sweep_0000_0V.vtk").exists())


if __name__ == "__main__":
    unittest.main()
