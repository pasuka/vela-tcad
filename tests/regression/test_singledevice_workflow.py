from __future__ import annotations

import importlib.util
import csv
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "singledevice_workflow", ROOT / "scripts" / "run_singledevice_workflow.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SingleDeviceWorkflowTest(unittest.TestCase):
    def test_exact_bias_comparison_csv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.csv"
            curve = root / "curve.csv"
            output = root / "comparison.csv"
            reference.write_text("bias_V,current_total\n-0.5,2e-14\n0.31,3e-7\n")
            curve.write_text(
                "bias_V,current_total_A_per_um,converged\n"
                "-0.5,-1.9e-14,1\n0.31,3.1e-7,1\n")
            MODULE.comparison_csv(reference, curve, output)
            with output.open(newline="") as handle:
                result = list(csv.DictReader(handle))
            self.assertEqual(2, len(result))
            self.assertEqual(1.9e-14, float(result[0]["vela_current_A_per_um"]))

    def write_reference(self, path: Path) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["bias_V", "current"])
            writer.writeheader()
            for bias in (-0.5, -0.365, -0.23):
                writer.writerow({"bias_V": bias, "current": 1.0})

    def base(self, drain: float) -> dict:
        return {
            "simulation_type": "dc_sweep",
            "mesh_file": "mesh.json", "node_doping_file": "doping.csv",
            "materials_file": "materials.json", "output_csv": "curve.csv",
            "contacts": [
                {"name": "source", "bias": 0.0},
                {"name": "drain", "bias": drain},
                {"name": "gate", "bias": -0.5},
                {"name": "substrate", "bias": 0.0},
            ],
            "solver": {
                "method": "newton", "quasi_fermi_update_limit_V": 0.025},
            "sweep": {"mode": "iv", "contact": "gate", "current_contact": "drain",
                      "start": -0.5, "stop": -0.23, "step": 0.135},
        }

    def test_save_load_dependency_and_reference_bias_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("mesh.json", "materials.json"):
                (root / name).write_text("{}")
            (root / "doping.csv").write_text("node_id,donors,acceptors\n")
            lin_base = root / "lin.json"
            sat_base = root / "sat.json"
            lin_base.write_text(json.dumps(self.base(0.1)))
            sat_base.write_text(json.dumps(self.base(1.1)))
            lin_ref = root / "lin.csv"
            sat_ref = root / "sat.csv"
            self.write_reference(lin_ref)
            self.write_reference(sat_ref)
            manifest = MODULE.materialize(
                lin_base, sat_base, lin_ref, sat_ref, root / "run")

            by_name = {stage["name"]: stage for stage in manifest["stages"]}
            common = manifest["common_saved_state"]
            self.assertEqual(Path(common).suffix, ".h5")
            self.assertTrue(all(Path(stage["final_state_file"]).suffix == ".h5"
                                for stage in manifest["stages"]))
            self.assertEqual(common, by_name["linear_drain_ramp"]["initial_state_file"])
            self.assertEqual(common, by_name["saturation_drain_ramp"]["initial_state_file"])
            self.assertEqual(
                by_name["linear_drain_ramp"]["final_state_file"],
                by_name["linear_idvg"]["initial_state_file"])
            self.assertEqual(
                by_name["saturation_drain_ramp"]["final_state_file"],
                by_name["saturation_idvg"]["initial_state_file"])
            self.assertEqual([-0.5, -0.365, -0.23],
                             manifest["reference_gate_biases_V"])
            self.assertEqual(5, len(by_name["saturation_idvg"]["bias_points"]))
            self.assertEqual([], by_name["equilibrium"]["depends_on"])
            self.assertEqual(["equilibrium"],
                             by_name["saturation_drain_ramp"]["depends_on"])
            equilibrium_cfg = json.loads(Path(by_name["equilibrium"]["config"]).read_text())
            linear_ramp_cfg = json.loads(
                Path(by_name["linear_drain_ramp"]["config"]).read_text())
            linear_curve_cfg = json.loads(
                Path(by_name["linear_idvg"]["config"]).read_text())
            saturation_curve_cfg = json.loads(
                Path(by_name["saturation_idvg"]["config"]).read_text())
            self.assertEqual(100, equilibrium_cfg["solver"]["max_iter"])
            self.assertEqual(
                0.1, equilibrium_cfg["solver"]["quasi_fermi_update_limit_V"])
            self.assertEqual(100, linear_ramp_cfg["solver"]["max_iter"])
            self.assertEqual(
                0.1, linear_curve_cfg["solver"]["quasi_fermi_update_limit_V"])
            self.assertEqual(
                0.025,
                saturation_curve_cfg["solver"]["quasi_fermi_update_limit_V"])


if __name__ == "__main__":
    unittest.main()
