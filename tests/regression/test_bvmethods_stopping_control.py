"""The stopping experiment must preserve physics and declare changed gates."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("stopping_control", ROOT/"scripts/run_bvmethods_stopping_control.py")
control = importlib.util.module_from_spec(spec)
spec.loader.exec_module(control)


class StoppingControlTests(unittest.TestCase):
    def test_physics_thresholds_and_original_are_preserved(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            for file in ("mesh.json", "doping.csv", "materials.json"):
                (root/file).write_text("{}")
            original = {
                "mesh_file":"mesh.json", "node_doping_file":"doping.csv", "materials_file":"materials.json",
                "solver":{"abstol":1e-10,"mobility":{"model":"constant"},
                    "impact_ionization":{"coupling_mode":"self_consistent"},
                    "carrier_row_convergence":{"mode":"enforce","eps_row":.001,"recovery":{"mode":"off"}},
                    "global_continuity_closure":{"mode":"enforce","tolerance":.01,"source_floor":1e-14}},
                "sweep":{"start":.1,"stop":1.,"bias_points":[.1,1.],"initial_state_file":"native.h5",
                    "breakdown":{"current_abs_A_per_um":.0013},
                    "external_circuit":{"solver":"coupled_newton","resistance_ohm_um":1e7,
                        "coupled_equation_tolerance":1e-8,"residual_tolerance_V":1e-6}}}
            source = root/"baseline.json"; source.write_text(json.dumps(original))
            arms = control.prepare(source, root/"experiment", 100000.)
            for arm, path in arms.items():
                result=json.loads((path/"config.json").read_text())
                solver=copy.deepcopy(result["solver"])
                solver["carrier_row_convergence"]["mode"]="enforce"
                solver["global_continuity_closure"]["mode"]="enforce"
                self.assertEqual(solver,original["solver"])
                self.assertNotIn("initial_state_file",result["sweep"])
                self.assertEqual(result["sweep"]["initialization"]["mode"],"poisson_block")
                self.assertEqual(result["sweep"]["breakdown"],original["sweep"]["breakdown"])
                circuit=copy.deepcopy(result["sweep"]["external_circuit"])
                update=circuit.pop("coupled_update_convergence")
                self.assertEqual(circuit,original["sweep"]["external_circuit"])
                self.assertEqual(update["mode"],"enforce" if arm=="density_update" else "report")
                self.assertEqual(result["solver"]["carrier_row_convergence"]["mode"],
                                 "enforce" if arm=="closure_report" else "report")
            self.assertEqual(json.loads(source.read_text()),original)
            with self.assertRaises(FileExistsError):
                control.prepare(source,root/"experiment",1.)


if __name__ == "__main__":
    unittest.main()
