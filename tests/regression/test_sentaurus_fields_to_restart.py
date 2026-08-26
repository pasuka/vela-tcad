import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class SentaurusFieldsToRestartTest(unittest.TestCase):
    def test_merges_regions_and_converts_density_to_m3(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = root / "fields"
            fields.mkdir()
            (root / "nodes.csv").write_text("node_id,x,y\n0,0,0\n1,1,0\n")
            entries = []
            for name in ("ElectrostaticPotential", "eQuasiFermiPotential",
                         "hQuasiFermiPotential", "eQuantumPotential",
                         "ConductionBandEnergy", "ElectronAffinity",
                         "BandgapNarrowing"):
                for region, node in enumerate((0, 1)):
                    filename = f"{name}_region{region}.csv"
                    value = (node + 0.1)
                    if name == "ConductionBandEnergy":
                        value = 4.6 - (node + 0.1) - 0.9
                    elif name == "ElectronAffinity":
                        value = 0.9
                    elif name == "BandgapNarrowing":
                        value = 0.0
                    (fields / filename).write_text(
                        f"node_id,component0\n{node},{value}\n")
                    entries.append({"name": name, "csv_file": filename,
                                    "unit": "V", "mapping_status": "complete"})
            for name in ("eDensity", "hDensity"):
                filename = f"{name}_region0.csv"
                (fields / filename).write_text("node_id,component0\n0,2e10\n")
                entries.append({"name": name, "csv_file": filename,
                                "unit": "cm^-3", "mapping_status": "complete"})
            (root / "field_manifest.json").write_text(json.dumps({"fields": entries}))
            (root / "elements.csv").write_text(
                "id,node0,node1,node2,region,material\n"
                "0,0,0,0,R.Silicon,Silicon\n"
                "1,1,1,1,R.Oxide,SiO2\n")
            output = root / "state.csv"
            subprocess.run([sys.executable, str(REPO / "scripts" / "sentaurus_fields_to_restart.py"),
                            "--export-dir", str(root), "--output", str(output)], check=True)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(float(rows[0]["electrons_m3"]), 2e16)
            self.assertEqual(float(rows[1]["electrons_m3"]), 0.0)
            self.assertEqual(float(rows[1]["phin"]), 0.0)
            self.assertEqual(float(rows[1]["phip"]), 0.0)
            self.assertEqual(float(rows[1]["electron_quantum_potential_V"]), 0.0)
            self.assertAlmostEqual(
                float(rows[0]["electron_quantum_potential_like_V"]),
                -0.9 - 1.5 * 0.025851999786 * math.log(1.0618016171622988))
            expected_insulator = (
                -0.9 - 1.5 * 0.025851999786 * math.log(0.42))
            self.assertAlmostEqual(
                float(rows[1]["electron_quantum_potential_like_V"]),
                expected_insulator)

    def test_classical_dd_explicitly_allows_missing_quantum_potential(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = root / "fields"
            fields.mkdir()
            (root / "nodes.csv").write_text("node_id,x,y\n0,0,0\n")
            entries = []
            values = {
                "ElectrostaticPotential": 0.1,
                "eQuasiFermiPotential": 0.0,
                "hQuasiFermiPotential": 0.0,
                "ConductionBandEnergy": 3.6,
                "ElectronAffinity": 0.9,
                "eDensity": 2.0e10,
                "hDensity": 3.0e10,
            }
            for name, value in values.items():
                filename = f"{name}_region0.csv"
                (fields / filename).write_text(
                    f"node_id,component0\n0,{value}\n")
                entries.append({
                    "name": name,
                    "csv_file": filename,
                    "unit": "cm^-3" if name.endswith("Density") else "V",
                    "mapping_status": "complete",
                })
            (root / "field_manifest.json").write_text(
                json.dumps({"fields": entries}))
            (root / "elements.csv").write_text(
                "id,node0,node1,node2,region,material\n"
                "0,0,0,0,R.Silicon,Silicon\n")
            output = root / "state.csv"
            subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "sentaurus_fields_to_restart.py"),
                "--export-dir", str(root),
                "--output", str(output),
                "--allow-missing-electron-quantum-potential",
            ], check=True)
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(float(row["electron_quantum_potential_V"]), 0.0)
            self.assertEqual(float(row["electrons_m3"]), 2.0e16)
            self.assertEqual(float(row["holes_m3"]), 3.0e16)


if __name__ == "__main__":
    unittest.main()
