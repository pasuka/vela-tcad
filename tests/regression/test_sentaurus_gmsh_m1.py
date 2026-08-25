"""M1 real Gmsh OCC, bundle, exporter, and fail-closed tests."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402
from sentaurus_gmsh_backend import GmshBackendError, validate_m1_device_ir  # noqa: E402
from sentaurus_ir_v2 import migrate_device_ir_v1  # noqa: E402

FIXTURE = REPO / "tests" / "fixtures" / "sentaurus_gmsh_m1" / "two_regions_sde.cmd"


class M1CapabilityTest(unittest.TestCase):
    def _ir(self, text: str) -> dict:
        with tempfile.TemporaryDirectory(prefix="vela_m1_capability_") as tmp:
            path = Path(tmp) / "case.cmd"
            path.write_text(text, encoding="utf-8")
            return migrate_device_ir_v1(parse_sde_device_ir(path))

    def test_checked_in_fixture_is_inside_m1_subset(self) -> None:
        validate_m1_device_ir(migrate_device_ir_v1(parse_sde_device_ir(FIXTURE)))

    def test_overlap_fails_closed_as_m2_boolean_semantics(self) -> None:
        ir = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.A")
(sdegeo:create-rectangle (position 0.5 0 0) (position 1.5 1 0) "Silicon" "R.B")
""")
        with self.assertRaises(GmshBackendError) as context:
            validate_m1_device_ir(ir)
        self.assertIn("ABA/BAB", str(context.exception))
        self.assertEqual("boolean-ownership", context.exception.details["capability"])

    def test_touching_interface_fails_closed(self) -> None:
        touching = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.A")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "R.B")
""")
        with self.assertRaises(GmshBackendError) as context:
            validate_m1_device_ir(touching)
        self.assertEqual("conformal-interface", context.exception.details["capability"])


class M1GmshIntegrationTest(unittest.TestCase):
    def test_two_region_pipeline_roundtrips_every_m1_format(self) -> None:
        prefix = Path(sys.prefix)
        wrapper = prefix / "lib" / "gmsh.py"
        library = prefix / "bin" / "libgmsh.dll"
        executable = prefix / "bin" / "gmsh.exe"
        if not all(path.is_file() for path in (wrapper, library, executable)):
            self.skipTest("Gmsh 4.15.2 UCRT64 runtime is not installed")
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m1_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable, str(REPO / "scripts" / "sentaurus_gmsh_pipeline.py"),
                    "--sde", str(FIXTURE), "--output-dir", str(output),
                    "--runtime-dir", str(root / "runtime"),
                    "--gmsh-python", str(wrapper), "--gmsh-library", str(library),
                    "--gmsh-executable", str(executable), "--mesh-size", "0.25",
                    "--capability-level", "M1",
                ],
                cwd=REPO, text=True, capture_output=True, check=False,
            )
            self.assertEqual(0, completed.returncode, completed.stderr)
            expected = {
                "device.msh", "device.vtu", "device.vtk", "device.grd", "device.dat",
                "mesh.json", "doping.csv", "mesh_field_bundle.json", "id_mapping.json",
                "device_ir.json", "mesh_qualification.json", "cross_format_report.json",
                "gmsh_runtime_manifest.json", "run_manifest.json",
            }
            self.assertEqual(expected, {path.name for path in output.iterdir() if path.is_file()})
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("generated", manifest["status"])
            self.assertEqual("passed", manifest["cross_format"]["status"])
            self.assertEqual("4.15.2", manifest["gmsh_version"])
            bundle = json.loads((output / "mesh_field_bundle.json").read_text(encoding="utf-8"))
            self.assertEqual("vela.mesh_field_bundle.v1", bundle["schema"])
            self.assertEqual([0, 1], [item["id"] for item in bundle["regions"]])
            self.assertEqual(list(range(len(bundle["nodes"]))),
                             [item["vela_id"] for item in bundle["nodes"]])
            self.assertEqual(list(range(len(bundle["triangles"]))),
                             [item["vela_id"] for item in bundle["triangles"]])
            self.assertTrue(all(item["gmsh_tag"] > 0 for item in bundle["nodes"]))
            self.assertTrue(all(item["gmsh_tag"] > 0 for item in bundle["triangles"]))
            self.assertTrue(all(item["gmsh_tag"] > 0 for item in bundle["lines"]))
            self.assertTrue(manifest["mesh_qualification"]["all_triangles_ccw"])
            self.assertEqual(2, manifest["mesh_qualification"]["connected_components"])
            with (output / "doping.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            left = [row for row in rows if float(row["x_um"]) <= 1.0]
            right = [row for row in rows if float(row["x_um"]) >= 1.5]
            self.assertTrue(left and right)
            self.assertEqual({(0.0, 1.0e16)},
                             {(float(row["donors_cm3"]), float(row["acceptors_cm3"])) for row in left})
            self.assertEqual({(2.0e16, 0.0)},
                             {(float(row["donors_cm3"]), float(row["acceptors_cm3"])) for row in right})
            msh_text = (output / "device.msh").read_text(encoding="utf-8", errors="replace")
            self.assertEqual(1, msh_text.count("$MeshFormat"))
            self.assertEqual(1, msh_text.count("$Nodes"))
            self.assertEqual(1, msh_text.count("$Elements"))
            self.assertEqual(5, msh_text.count("$NodeData"))
            mesh = json.loads((output / "mesh.json").read_text(encoding="utf-8"))
            self.assertEqual({"Si"}, {region["material"] for region in mesh["regions"]})
            vtu_text = (output / "device.vtu").read_text(encoding="utf-8")
            self.assertIn("ContactId", vtu_text)
            self.assertIn("PhysicalCurveTag", vtu_text)

            repeated = root / "output-repeat"
            repeat_command = list(completed.args)
            repeat_command[repeat_command.index(str(output))] = str(repeated)
            repeat_runtime = root / "runtime-repeat"
            repeat_command[repeat_command.index(str(root / "runtime"))] = str(repeat_runtime)
            repeat = subprocess.run(
                repeat_command, cwd=REPO, text=True, capture_output=True, check=False,
            )
            self.assertEqual(0, repeat.returncode, repeat.stderr)
            deterministic = expected - {"gmsh_runtime_manifest.json", "run_manifest.json"}
            for name in sorted(deterministic):
                first = hashlib.sha256((output / name).read_bytes()).hexdigest()
                second = hashlib.sha256((repeated / name).read_bytes()).hexdigest()
                self.assertEqual(first, second, f"non-deterministic M1 artifact: {name}")

            overlap = root / "overlap.cmd"
            overlap.write_text(
                '(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.A")\n'
                '(sdegeo:create-rectangle (position 0.5 0 0) (position 1.5 1 0) "Silicon" "R.B")\n',
                encoding="utf-8",
            )
            failed_output = root / "failed-output"
            failed_output.mkdir()
            (failed_output / "run_manifest.json").write_text('{"status":"generated"}\n')
            fail_command = list(completed.args)
            fail_command[fail_command.index(str(FIXTURE))] = str(overlap)
            fail_command[fail_command.index(str(output))] = str(failed_output)
            fail_runtime = root / "runtime-fail"
            fail_command[fail_command.index(str(root / "runtime"))] = str(fail_runtime)
            failed = subprocess.run(
                fail_command, cwd=REPO, text=True, capture_output=True, check=False,
            )
            self.assertNotEqual(0, failed.returncode)
            self.assertFalse((failed_output / "run_manifest.json").exists())
            report = json.loads((failed_output / "fail_report.json").read_text(encoding="utf-8"))
            self.assertEqual("failed", report["status"])
            self.assertEqual("gmsh_backend", report["stage"])
            self.assertEqual("boolean-ownership", report["details"]["capability"])


if __name__ == "__main__":
    unittest.main()
