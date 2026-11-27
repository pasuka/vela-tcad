"""Unit provenance and production-CLI coverage for IALMob point comparisons."""
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.prepare_slot_ldmos_ialmob_probe import extract_parameters, nearest_segment, prepare
from scripts.compare_slot_ldmos_ialmob_probe import compare

REPO = Path(__file__).resolve().parents[2]


class IALMobToolsTest(unittest.TestCase):
    def test_source_parameters_preserve_cm_coefficients_and_carrier_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "material.par"
            path.write_text('IALMob:\n{\n n_ref = 9.68e16, 2.23e17 # cm^-3\n alpha_sr = 5e-22, 0\n me_over_m0 = 1\n}\nIALMob "100":\n{\n n_ref=123\n}\n')
            values = extract_parameters(path)
            self.assertEqual(values["electron"]["n_ref"], 9.68e16)
            self.assertEqual(values["hole"]["n_ref"], 2.23e17)
            self.assertEqual(values["electron"]["alpha_sr"], 5e-22)
            self.assertEqual(values["me_over_m0"], 1)
            for body in ["n_ref=1\nn_ref=2", "n_ref=nan", "n_ref=1,2,3"]:
                path.write_text("IALMob:\n{\n" + body + "\n}\n")
                with self.assertRaises(ValueError):
                    extract_parameters(path)

    def test_geometric_distance_gradient_includes_segment_end(self):
        self.assertEqual(nearest_segment((.5, 2), [(0, 0, 1, 0)]), (2, 0, 1))
        distance, nx, ny = nearest_segment((2, 1), [(0, 0, 1, 0)])
        self.assertAlmostEqual(distance, math.sqrt(2))
        self.assertAlmostEqual(nx, 1 / math.sqrt(2))
        self.assertAlmostEqual(ny, nx)
        with self.assertRaises(ValueError):
            nearest_segment((0, 0), [])

    def fixture(self, root):
        export = root / "export"
        (export / "fields").mkdir(parents=True)
        fields = {"eDensity": ("cm^-3", [2e10]), "hDensity": ("cm^-3", [3e10]),
            "DonorConcentration": ("cm^-3", [1e17]), "AcceptorConcentration": ("cm^-3", [4e16]),
            "ElectricField": ("V*cm^-1", [3, 4]), "eGradQuasiFermi": ("V*cm^-1", [6, 8]),
            "hGradQuasiFermi": ("V*cm^-1", [5, 12]), "eMobility": ("cm^2*V^-1*s^-1", [100]),
            "hMobility": ("cm^2*V^-1*s^-1", [50]), "eEnormal": ("V*cm^-1", [7]),
            "hEnormal": ("V*cm^-1", [9]), "LatticeTemperature": ("K", [300])}
        metadata = {"fields": [{"name": key, "unit": unit, "region": 0} for key, (unit, _) in fields.items()],
                    "regions": [{"index": 0, "name": "Silicon_1", "material": "Silicon"}]}
        (export / "metadata.json").write_text(json.dumps(metadata))
        for name, (_, values) in fields.items():
            with (export / "fields" / (name + "_region0.csv")).open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["node_id"] + [f"component{i}" for i in range(len(values))])
                writer.writerow([0] + values)
        (export / "nodes.csv").write_text("id,x_um,y_um\n0,0,1\n1,0,0\n2,1,0\n3,0,-1\n")
        (export / "elements.csv").write_text("id,node0,node1,node2,region,material\n0,0,1,2,Silicon_1,Si\n1,1,3,2,Oxide_1,SiO2\n")
        parameters = root / "source.par"
        parameters.write_text("IALMob:\n{\n mumin=56,46\n}\n")
        return export, parameters

    def test_native_fields_are_case_sensitive_and_units_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export, parameters = self.fixture(root)
            output = root / "probe"
            prepare(export, parameters, output)
            point = json.loads((output / "points.json").read_text())[0]
            self.assertEqual(point["donors_m3"], 1e23)
            self.assertEqual(point["electrons_m3"], 2e16)
            self.assertEqual(point["electron_driving_field_V_per_m"], 1000)
            self.assertEqual(point["electron_normal_field_V_per_m"], 700)
            self.assertEqual(point["hole_normal_field_V_per_m"], 900)
            self.assertEqual(point["interface_distance_m"], 1e-6)
            temperature = export / "fields/LatticeTemperature_region0.csv"
            temperature.write_text("node_id,component0\n0,350\n")
            with self.assertRaisesRegex(ValueError, "uniform 300 K"):
                prepare(export, parameters, output)
            temperature.write_text("node_id,component0\n0,300\n")
            metadata = json.loads((export / "metadata.json").read_text())
            metadata["fields"][-1]["unit"] = "V/m"
            (export / "metadata.json").write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, "Unexpected unit"):
                prepare(export, parameters, output)

    def test_production_cli_and_comparison_require_both_carriers(self):
        configured = os.environ.get("VELA_EXAMPLE_RUNNER")
        runner = Path(configured) if configured else REPO / "build" / ("vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner")
        if not runner.is_file() and not configured:
            self.skipTest("Build the example runner to exercise the production probe")
        self.assertTrue(runner.is_file())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export, parameters = self.fixture(root)
            output = root / "probe"
            prepare(export, parameters, output)
            result = subprocess.run([str(runner), "--config", str(output / "config.json"), "--log", "off"], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(compare(output)["electron"]["all"]["count"], 1)
            lines = (output / "mobility.csv").read_text().splitlines()
            (output / "mobility.csv").write_text("\n".join(lines[:-1]) + "\n")
            with self.assertRaisesRegex(ValueError, "Missing mobility"):
                compare(output)
            config = json.loads((output / "config.json").read_text())
            config["scaling"] = {"mode": "tcad"}
            (output / "config.json").write_text(json.dumps(config))
            rejected = subprocess.run([str(runner), "--config", str(output / "config.json"), "--log", "off"], capture_output=True, text=True, timeout=30)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("explicitly SI", rejected.stdout + rejected.stderr)


if __name__ == "__main__":
    unittest.main()
