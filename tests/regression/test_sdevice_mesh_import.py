"""WP2 regression tests: Sentaurus grid import -> Vela mesh bundle + C2 qualification report."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from vela_sdevice import schema_validator as sv  # noqa: E402
from vela_sdevice.mesh_import import MeshImportError, check_electrodes, find_importer, import_grid  # noqa: E402
import run_pn2d_same_mesh_sdevice_vector_export as dfise_writers  # noqa: E402

SCHEMA = sv.load_schema(REPO / "schemas" / "vela.mesh_qualification.v1.schema.json")
TMP_ROOT = REPO / "build" if (REPO / "build").is_dir() else None
RUNNER_CANDIDATES = [REPO / "build" / "vela_example_runner.exe", Path(r"D:\code-repo\vela-tcad\build\vela_example_runner.exe"),
                     REPO / "build" / "vela_example_runner"]

# 3x2 node strip on [0,2]x[0,1] um: four right triangles, cathode on x=0, anode on x=2.
NODES = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (2.0, 0.0), 3: (0.0, 1.0), 4: (1.0, 1.0), 5: (2.0, 1.0)}
TRIANGLES = [(0, 1, 4), (0, 4, 3), (1, 2, 5), (1, 5, 4)]
CONTACTS = {"cathode": [0, 3], "anode": [2, 5]}
DOPING = {0: (1e17, 0.0), 1: (1e17, 0.0), 2: (0.0, 2e17), 3: (1e17, 0.0), 4: (1e17, 0.0), 5: (0.0, 2e17)}

# Hand-written two-material grid with one obtuse triangle (node 2 sits close to the base 0-1).
OBTUSE_GRD = """DF-ISE text

Info {
  version = 1.1
  type = grid
  dimension = 2
  nb_vertices = 5
  nb_edges = 8
  nb_faces = 0
  nb_elements = 5
  nb_regions = 3
  regions = [ "R.Si" "R.Ox" "gate" ]
  materials = [ Silicon Oxide Contact ]
}

Data {
  CoordSystem {
    translate = [ 0 0 0 ]
    transform = [ 1 0 0 0 1 0 0 0 1 ]
  }

  Vertices (5) {
 0.0 0.0
 1.0 0.0
 0.5 0.2
 0.0 1.0
 1.0 1.0
  }

  Edges (8) {
 0 1
 1 2
 2 0
 2 3
 3 0
 1 4
 4 2
 4 3
  }

  Locations (8) {
eiiieiei
  }

  Elements (5) {
 2 0 1 2
 2 -3 3 4
 2 -2 5 6
 2 -7 7 -4
 1 3 4
  }

  Region ("R.Si") {
    material = Silicon
    Elements (1) {
 0
    }
  }

  Region ("R.Ox") {
    material = Oxide
    Elements (3) {
 1 2 3
    }
  }

  Region ("gate") {
    material = Contact
    Elements (1) {
 4
    }
  }
}
"""

OBTUSE_DAT = """DF-ISE text

Info {
  version = 1.0
  type = dataset
  dimension = 2
  nb_vertices = 5
  nb_edges = 8
  nb_faces = 0
  nb_elements = 5
  nb_regions = 3
  datasets = [ "DopingConcentration" ]
  functions = [ DopingConcentration ]
}

Data {
  Dataset ("DopingConcentration") {
    function = DopingConcentration
    type = scalar
    dimension = 1
    location = vertex
    validity = [ "R.Si" ]
    Values (5) {
 1e16
 1e16
 -5e15
 0
 0
    }
  }
}
"""


def write_strip_fixture(directory: Path) -> Path:
    doping_csv = directory / "vela_doping.csv"
    with doping_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "donors_cm3", "acceptors_cm3"])
        for node_id, (donors, acceptors) in DOPING.items():
            writer.writerow([node_id, repr(donors), repr(acceptors)])
    grd = directory / "strip_msh.grd"
    grd.write_text(dfise_writers.render_dfise_grid(NODES, TRIANGLES, CONTACTS), encoding="utf-8")
    grd.with_suffix(".dat").write_text(dfise_writers.render_dfise_doping(NODES, TRIANGLES, CONTACTS, doping_csv),
                                       encoding="utf-8")
    return grd


def read_doping(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(r["node_id"]): (float(r["donors_cm3"]), float(r["acceptors_cm3"])) for r in csv.DictReader(handle)}


def runner() -> Path | None:
    for candidate in RUNNER_CANDIDATES:
        if candidate.is_file():
            return candidate
    return None


class DfiseImportTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="wp2_mesh_", dir=TMP_ROOT)
        self.tmp = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_strip_fixture_round_trips_and_validates(self) -> None:
        grd = write_strip_fixture(self.tmp)
        out = self.tmp / "out"
        report = import_grid(grd, out)
        self.assertEqual([], sv.collect_errors(SCHEMA, report))
        self.assertEqual([], sv.collect_errors(SCHEMA, json.loads((out / "mesh_qualification.json").read_text())))
        mesh = json.loads((out / "mesh.json").read_text(encoding="utf-8"))
        self.assertEqual(len(NODES), len(mesh["nodes"]))
        self.assertEqual([[NODES[n["id"]][0], NODES[n["id"]][1]] for n in mesh["nodes"]],
                         [[n["x"], n["y"]] for n in mesh["nodes"]])
        self.assertEqual(sorted(sorted(t) for t in TRIANGLES), sorted(sorted(t["node_ids"]) for t in mesh["triangles"]))
        for tri in mesh["triangles"]:
            p = [mesh["nodes"][i] for i in tri["node_ids"]]
            area2 = (p[1]["x"] - p[0]["x"]) * (p[2]["y"] - p[0]["y"]) - (p[2]["x"] - p[0]["x"]) * (p[1]["y"] - p[0]["y"])
            self.assertGreater(area2, 0.0, "triangles must be counter-clockwise")
        self.assertEqual([("R.Si", "Si")], [(r["name"], r["material"]) for r in mesh["regions"]])
        contacts = {c["name"]: c for c in mesh["contacts"]}
        self.assertEqual({"Cathode": [0, 3], "Anode": [2, 5]}, {k: v["node_ids"] for k, v in contacts.items()})
        self.assertEqual(1, len(contacts["Cathode"]["edge_node_ids"]))
        self.assertEqual(DOPING, read_doping(out / "doping.csv"))
        # report content
        self.assertEqual("dfise", report["source"]["format"])
        self.assertEqual("strip_msh.dat", report["source"]["companion_file"])
        self.assertEqual({"min": [0.0, 0.0], "max": [2.0, 1.0]}, report["bounding_box_um"])
        self.assertEqual(6, report["nodes"]["count"])
        self.assertEqual(4, report["elements"]["count"])
        self.assertEqual(4, report["elements"]["by_kind"]["triangle"])
        self.assertEqual(2, report["elements"]["by_kind"]["line"])
        self.assertEqual(["Si"], report["materials"])
        self.assertEqual([("Cathode", 2, 1), ("Anode", 2, 1)],
                         [(c["name"], c["node_count"], c["edge_count"]) for c in report["contacts"]])
        self.assertAlmostEqual(90.0, report["angles"]["max_deg"], places=9)
        self.assertEqual(0, report["angles"]["obtuse_element_count"])
        self.assertTrue(report["qualification"]["non_obtuse"])
        self.assertTrue(report["qualification"]["coordinates_plausible"])
        self.assertEqual("barycentric", report["qualification"]["recommended_node_volume_policy"])
        self.assertEqual(["PhosphorusActiveConcentration", "BoronActiveConcentration"], report["doping"]["source_fields"])
        self.assertFalse(report["doping"]["net_only"])

    def test_electrode_check(self) -> None:
        grd = write_strip_fixture(self.tmp)
        report = import_grid(grd, self.tmp / "out")
        result = check_electrodes(report, ["Anode", "Cathode", "gate"])
        self.assertEqual([True, True, False], [r["matched"] for r in result])
        self.assertEqual(result, report["qualification"]["electrodes_checked"])
        self.assertEqual([], sv.collect_errors(SCHEMA, report))

    def test_obtuse_two_material_grid_is_reported_not_rejected(self) -> None:
        grd = self.tmp / "obtuse_msh.grd"
        grd.write_text(OBTUSE_GRD, encoding="utf-8")
        grd.with_suffix(".dat").write_text(OBTUSE_DAT, encoding="utf-8")
        report = import_grid(grd, self.tmp / "out")
        self.assertEqual([], sv.collect_errors(SCHEMA, report))
        self.assertEqual(1, report["angles"]["obtuse_element_count"])
        self.assertFalse(report["qualification"]["non_obtuse"])
        self.assertTrue(report["qualification"]["allow_obtuse"])
        self.assertEqual(["Si", "SiO2"], report["materials"])
        self.assertEqual([("R.Si", "Si", 1), ("R.Ox", "SiO2", 3)],
                         [(r["name"], r["material"], r["element_count"]) for r in report["regions"]])
        self.assertEqual("Oxide", report["regions"][1]["source_material"])
        self.assertEqual([{"name": "gate", "node_count": 2, "edge_count": 1, "adjacent_materials": ["SiO2"]}],
                         report["contacts"])
        self.assertTrue(report["doping"]["net_only"])
        doping = read_doping(self.tmp / "out" / "doping.csv")
        self.assertEqual((1e16, 0.0), doping[0])
        self.assertEqual((0.0, 5e15), doping[2])
        mesh = json.loads((self.tmp / "out" / "mesh.json").read_text(encoding="utf-8"))
        self.assertEqual(1, mesh["contacts"][0]["region_id"])

    def test_missing_dat_and_bad_extension(self) -> None:
        grd = write_strip_fixture(self.tmp)
        grd.with_suffix(".dat").unlink()
        with self.assertRaises(MeshImportError) as ctx:
            import_grid(grd, self.tmp / "out")
        self.assertEqual("E304", ctx.exception.message_id)
        other = self.tmp / "x.xyz"
        other.write_text("nothing", encoding="utf-8")
        with self.assertRaises(MeshImportError) as ctx:
            import_grid(other, self.tmp / "out")
        self.assertEqual("E114", ctx.exception.message_id)
        with self.assertRaises(MeshImportError) as ctx:
            import_grid(self.tmp / "absent.grd", self.tmp / "out")
        self.assertEqual("E302", ctx.exception.message_id)

    def test_quadrilateral_grid_is_rejected(self) -> None:
        grd = self.tmp / "quad_msh.grd"
        text = OBTUSE_GRD.replace(" 2 0 1 2\n", " 3 0 1 2 -5\n")
        grd.write_text(text, encoding="utf-8")
        grd.with_suffix(".dat").write_text(OBTUSE_DAT, encoding="utf-8")
        with self.assertRaises(MeshImportError) as ctx:
            import_grid(grd, self.tmp / "out")
        self.assertEqual("E402", ctx.exception.message_id)
        self.assertIn("quadrilaterals=1", ctx.exception.message)

    def test_tdr_without_importer_fails_closed(self) -> None:
        tdr = self.tmp / "x_msh.tdr"
        tdr.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\x00" * 16)
        if find_importer(None) is not None:
            self.skipTest("a sentaurus_import executable is available; the no-importer path is not exercised here")
        with self.assertRaises(MeshImportError) as ctx:
            import_grid(tdr, self.tmp / "out")
        self.assertEqual("E402", ctx.exception.message_id)

    def test_cli_exit_codes(self) -> None:
        grd = write_strip_fixture(self.tmp)
        env = {**os.environ, "PYTHONPATH": str(REPO / "scripts")}
        proc = subprocess.run([sys.executable, "-m", "vela_sdevice.mesh_import", str(grd), str(self.tmp / "cli_out")],
                              capture_output=True, text=True, env=env, cwd=str(REPO))
        self.assertEqual(0, proc.returncode, proc.stderr)
        self.assertEqual(SCHEMA["$id"], json.loads(proc.stdout)["schema"])
        proc = subprocess.run([sys.executable, "-m", "vela_sdevice.mesh_import", str(self.tmp / "nope.grd"), str(self.tmp / "o")],
                              capture_output=True, text=True, env=env, cwd=str(REPO))
        self.assertEqual(3, proc.returncode)
        self.assertIn("E302", proc.stderr)

    def test_runner_loads_generated_mesh(self) -> None:
        exe = runner()
        if exe is None:
            self.skipTest("vela_example_runner not built; mesh.json load smoke test skipped")
        grd = write_strip_fixture(self.tmp)
        out = self.tmp / "out"
        import_grid(grd, out)
        deck = {
            "simulation_type": "poisson",
            "mesh_file": "mesh.json",
            "output_vtk": "poisson.vtk",
            "scaling": {"mode": "unit_scaling"},
            "doping": [{"region": "R.Si", "donors": 1e17, "acceptors": 0.0}],
            "contacts": [{"name": "Cathode", "bias": 0.0}, {"name": "Anode", "bias": 0.0}],
            "runtime_log": {"enabled": False},
        }
        (out / "smoke.json").write_text(json.dumps(deck, indent=2), encoding="utf-8")
        env = {**os.environ, "Path": r"D:\msys64\ucrt64\bin;D:\msys64\usr\bin;" + os.environ.get("Path", "")}
        proc = subprocess.run([str(exe), "--config", "smoke.json"], cwd=str(out), capture_output=True, text=True, env=env)
        self.assertEqual(0, proc.returncode, proc.stdout[-2000:] + proc.stderr[-2000:])
        self.assertTrue((out / "poisson.vtk").is_file())


if __name__ == "__main__":
    unittest.main()
