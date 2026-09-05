from __future__ import annotations

import importlib.util
import csv
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts/run_pn2d_same_mesh_sdevice_vector_export.py"
)
SPEC = importlib.util.spec_from_file_location("pn2d_same_mesh_export", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Pn2dSameMeshSdeviceVectorExportTest(unittest.TestCase):
    def test_dfise_grid_preserves_nodes_triangles_and_contact_segments(self) -> None:
        nodes = {
            0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0), 3: (1.0, 1.0)
        }
        triangles = [(0, 1, 3), (0, 3, 2)]
        contacts = {"anode": [0, 2], "cathode": [1, 3]}
        text = MODULE.render_dfise_grid(nodes, triangles, contacts)
        self.assertIn("nb_vertices = 4", text)
        self.assertIn("nb_elements = 4", text)
        self.assertIn("Elements (2)", text)
        self.assertIn('Region ("Anode")', text)
        self.assertIn('Region ("Cathode")', text)

    def test_dfise_doping_preserves_signed_and_species_concentrations(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            doping = Path(directory) / "doping.csv"
            doping.write_text(
                "node_id,donors_cm3,acceptors_cm3\n"
                "0,0,1e17\n1,1e17,0\n",
                encoding="utf-8",
            )
            nodes = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}
            doping.write_text(
                "node_id,donors_cm3,acceptors_cm3\n"
                "0,0,1e17\n1,1e17,0\n2,0,1e17\n",
                encoding="utf-8",
            )
            text = MODULE.render_dfise_doping(
                nodes,
                [(0, 1, 2)],
                {"anode": [0, 2], "cathode": [1, 2]},
                doping,
            )
            self.assertIn("-1.000000000000000e+17", text)
            self.assertIn("1.000000000000000e+17", text)
            self.assertIn("nb_edges = 3", text)
            self.assertIn("nb_elements = 3", text)
            self.assertIn("PhosphorusActiveConcentration", text)
            self.assertIn("BoronActiveConcentration", text)

    def test_coordinate_mapping_accepts_reordered_exact_nodes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            export = Path(directory)
            (export / "nodes.csv").write_text(
                "id,x_um,y_um\n0,1,0\n1,0,0\n",
                encoding="utf-8",
            )
            mapping = MODULE.coordinate_mapping({0: (0.0, 0.0), 1: (1.0, 0.0)}, export)
            self.assertEqual(mapping, {0: 1, 1: 0})

    def test_deck_requests_two_component_current_density_vectors(self) -> None:
        deck = MODULE.render_sdevice_deck(-20.0)
        self.assertIn('Grid = "pn2d_same_mesh.tdr"', deck)
        self.assertIn('Doping = "pn2d_same_mesh.tdr"', deck)
        self.assertIn("eCurrentDensity/Vector", deck)
        self.assertIn("hCurrentDensity/Vector", deck)
        self.assertIn('Name="Anode" Voltage=-20', deck)

    def test_classical_dd_state_converts_density_to_per_cubic_metre(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = root / "fields"
            fields.mkdir()
            values = {
                "ElectrostaticPotential": 1.0,
                "eQuasiFermiPotential": 2.0,
                "hQuasiFermiPotential": 3.0,
                "eDensity": 4.0,
                "hDensity": 5.0,
            }
            for name, value in values.items():
                (fields / f"{name}_region0.csv").write_text(
                    f"node_id,component0\n0,{value}\n", encoding="utf-8"
                )
            output = root / "state.csv"
            MODULE.write_classical_dd_state(root, output)
            with output.open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(float(row["electrons_m3"]), 4.0e6)
            self.assertEqual(float(row["holes_m3"]), 5.0e6)


if __name__ == "__main__":
    unittest.main()
