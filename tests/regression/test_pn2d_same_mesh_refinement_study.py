from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/run_pn2d_same_mesh_refinement_study.py"
SPEC = importlib.util.spec_from_file_location("pn_refinement", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Pn2dSameMeshRefinementStudyTest(unittest.TestCase):
    def test_nested_levels_preserve_expected_counts_and_old_vertices(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = MODULE.prepare(Path(directory))
            by = {level["label"]: level for level in manifest["levels"]}
            self.assertEqual((by["7x3"]["nodes"], by["7x3"]["triangles"]), (21, 24))
            self.assertEqual((by["13x5"]["nodes"], by["13x5"]["triangles"]), (65, 96))
            self.assertEqual((by["25x9"]["nodes"], by["25x9"]["triangles"]), (225, 384))
            coarse = MODULE.read_json(Path(by["7x3"]["mesh"]))
            fine = MODULE.read_json(Path(by["25x9"]["mesh"]))
            coarse_xy = {(row["x"], row["y"]) for row in coarse["nodes"]}
            fine_xy = {(row["x"], row["y"]) for row in fine["nodes"]}
            self.assertTrue(coarse_xy <= fine_xy)

    def test_refined_doping_preserves_coarse_node_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = MODULE.prepare(Path(directory))
            by = {level["label"]: level for level in manifest["levels"]}
            coarse_mesh = MODULE.read_json(Path(by["7x3"]["mesh"]))
            fine_mesh = MODULE.read_json(Path(by["25x9"]["mesh"]))
            coarse_doping = MODULE.read_csv(Path(by["7x3"]["doping"]))
            fine_doping = MODULE.read_csv(Path(by["25x9"]["doping"]))
            coarse = {(row["x"], row["y"]): coarse_doping[row["id"]] for row in coarse_mesh["nodes"]}
            fine = {(row["x"], row["y"]): fine_doping[row["id"]] for row in fine_mesh["nodes"]}
            for coordinate, expected in coarse.items():
                self.assertEqual(float(fine[coordinate]["donors_cm3"]), float(expected["donors_cm3"]))
                self.assertEqual(float(fine[coordinate]["acceptors_cm3"]), float(expected["acceptors_cm3"]))

    def test_final_plt_currents_uses_last_native_row(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "current.plt"
            path.write_text(
                'Info { datasets = [ "time" "Cathode TotalCurrent" "Anode TotalCurrent" ] }\n'
                'Data { 0 1 -1  1 2 -2 }\n',
                encoding="utf-8",
            )
            currents = MODULE.final_plt_currents(path)
            self.assertEqual(currents["cathode_A_per_um"], 2.0)
            self.assertEqual(currents["anode_A_per_um"], -2.0)


if __name__ == "__main__":
    unittest.main()
