"""M2 ABA/BAB, contact, ordered-doping, Gaussian, and exporter tests."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402
from sentaurus_doping import evaluate_m2_doping_at  # noqa: E402
from sentaurus_gmsh_backend import (  # noqa: E402
    GmshBackendError,
    validate_m2_device_ir,
)
from sentaurus_ir_v2 import migrate_device_ir_v1  # noqa: E402
from sentaurus_mesh_exporters import bundle_to_vela_mesh  # noqa: E402

FIXTURE = REPO / "tests" / "fixtures" / "sentaurus_gmsh_m2" / "m2_core_sde.cmd"
ORACLE = FIXTURE.with_name("sentaurus2022_oracle_summary.json")


def parse_ir(path: Path) -> dict:
    return migrate_device_ir_v1(parse_sde_device_ir(path))


class M2DopingUnitTest(unittest.TestCase):
    def _ir(self, text: str) -> dict:
        with tempfile.TemporaryDirectory(prefix="vela_m2_doping_") as tmp:
            path = Path(tmp) / "case.cmd"
            path.write_text(text, encoding="utf-8")
            return parse_ir(path)

    def test_fixture_has_complete_provenance_and_analytic_gaussian_values(self) -> None:
        ir = parse_ir(FIXTURE)
        oracle = json.loads(ORACLE.read_text(encoding="utf-8"))
        self.assertEqual(oracle["source_sha256"], hashlib.sha256(FIXTURE.read_bytes()).hexdigest())
        self.assertEqual("Sentaurus Structure Editor T-2022.03-SP2", oracle["tool"])
        self.assertEqual(0, oracle["exit_code"])
        self.assertTrue(ir["capability_manifest"]["provenance_complete"])
        validate_m2_device_ir(ir)

        donor, acceptor, species = evaluate_m2_doping_at(ir, 2.0, 1.0)
        self.assertAlmostEqual(1.02e18, donor, delta=1.0e5)
        self.assertEqual(0.0, acceptor)
        self.assertAlmostEqual(1.0e18, species["ArsenicActiveConcentration"], delta=1.0e5)
        self.assertTrue(math.isclose(
            evaluate_m2_doping_at(ir, 1.0, 1.0)[2]["ArsenicActiveConcentration"],
            oracle["region_local_fields_cm3"]["R.N"][
                "ArsenicActiveConcentration_at_x1_y1"],
            rel_tol=1.0e-12,
        ))

        # East-going baseline: the T-2022.03-SP2 Positive side is its left
        # normal (positive y); the lower side has no Gaussian contribution.
        self.assertEqual(2.0e16, evaluate_m2_doping_at(ir, 2.0, 0.6)[0])
        # NoEval excludes the baseline itself.
        self.assertEqual(2.0e16, evaluate_m2_doping_at(ir, 2.0, 0.8)[0])

        placement = ir["doping"]["placements"][-1]
        placement["direction"] = "Negative"
        placement["evaluation"] = "Eval"
        self.assertAlmostEqual(1.02e18, evaluate_m2_doping_at(ir, 2.0, 0.6)[0], delta=1.0e5)
        self.assertEqual(2.0e16, evaluate_m2_doping_at(ir, 2.0, 1.0)[0])
        placement["direction"] = "Both"
        self.assertAlmostEqual(1.02e18, evaluate_m2_doping_at(ir, 2.0, 0.6)[0], delta=1.0e5)
        self.assertAlmostEqual(1.02e18, evaluate_m2_doping_at(ir, 2.0, 1.0)[0], delta=1.0e5)

    def test_ordered_noreplace_localreplace_and_replace(self) -> None:
        prefix = """
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdedr:define-constant-profile "P.old" "PhosphorusActiveConcentration" 10)
(sdedr:define-constant-profile-region "P.old.place" "P.old" "R.Si")
(sdedr:define-constant-profile "As" "ArsenicActiveConcentration" 20)
(sdedr:define-constant-profile-region "As.place" "As" "R.Si")
(sdedr:define-constant-profile "B" "BoronActiveConcentration" 30)
(sdedr:define-constant-profile-region "B.place" "B" "R.Si")
(sdedr:define-constant-profile "P.new" "PhosphorusActiveConcentration" 5)
(sdedr:define-constant-profile-region "P.new.place" "P.new" "R.Si" 0 "LocalReplace")
"""
        ir = self._ir(prefix)
        donor, acceptor, species = evaluate_m2_doping_at(ir, 0.5, 0.5)
        self.assertEqual(25.0, donor)
        self.assertEqual(30.0, acceptor)
        self.assertEqual(5.0, species["PhosphorusActiveConcentration"])
        self.assertEqual(20.0, species["ArsenicActiveConcentration"])

        replace = self._ir(prefix + """
(sdedr:define-constant-profile "B.replace" "BoronActiveConcentration" 7)
(sdedr:define-constant-profile-region "B.replace.place" "B.replace" "R.Si" 0 "Replace")
""")
        donor, acceptor, species = evaluate_m2_doping_at(replace, 0.5, 0.5)
        self.assertEqual(0.0, donor)
        self.assertEqual(7.0, acceptor)
        self.assertEqual({"BoronActiveConcentration": 7.0}, species)

    def test_constant_window_target_and_gaussian_depth_dose_erf_restriction(self) -> None:
        window = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 2 1 0) "Silicon" "R.Si")
(sdedr:define-refeval-window "W" "Rectangle" (position 0 0 0) (position 1 1 0))
(sdedr:define-constant-profile "D" "PhosphorusActiveConcentration" 4e15)
(sdedr:define-constant-profile-placement "D.place" "D" "W")
""")
        self.assertEqual(4e15, evaluate_m2_doping_at(window, 0.5, 0.5)[0])
        self.assertEqual(0.0, evaluate_m2_doping_at(window, 1.5, 0.5)[0])

        depth = self._ir("""
(sdegeo:create-rectangle (position -2 -1 0) (position 2 1 0) "Silicon" "R.Si")
(sdedr:define-refeval-window "B" "Line" (position -1 0 0) (position 1 0 0))
(sdedr:define-refeval-window "Limit" "Rectangle" (position -0.5 0 0) (position 0.5 1 0))
(sdedr:define-gaussian-profile "G" "PhosphorusActiveConcentration"
 "PeakPos" 0.2 "PeakVal" 1e18 "ValueAtDepth" 1e16 "Depth" 0.4
 "Erf" "Factor" 1)
(sdedr:define-analytical-profile-placement
 "G.place" "G" "B" "Positive" "NoReplace" "Eval" "Limit" 0 "evalwin")
""")
        self.assertTrue(math.isclose(
            evaluate_m2_doping_at(depth, 0.0, 0.4)[0], 1e16, rel_tol=1.0e-12))
        self.assertEqual(0.0, evaluate_m2_doping_at(depth, 0.8, 0.2)[0])

        dose = self._ir("""
(sdegeo:create-rectangle (position -2 -1 0) (position 2 1 0) "Silicon" "R.Si")
(sdedr:define-refeval-window "B" "Line" (position -1 0 0) (position 1 0 0))
(sdedr:define-gaussian-profile "Dose.G" "ArsenicActiveConcentration"
 "PeakPos" 0 "Dose" 1e12 "StdDev" 0.1 "Gauss" "Factor" 1)
(sdedr:define-analytical-profile-placement
 "Dose.place" "Dose.G" "B" "Both" "NoReplace" "Eval")
""")
        expected_peak = 1e12 / (math.sqrt(2.0 * math.pi) * 0.1e-4)
        self.assertTrue(math.isclose(
            evaluate_m2_doping_at(dose, 0.0, 0.0)[0], expected_peak,
            rel_tol=1.0e-12))

    def test_boolean_boundary_policy_matches_aba_and_bab(self) -> None:
        common = """
(sdedr:define-constant-profile "Old" "BoronActiveConcentration" 10)
(sdedr:define-constant-profile-region "Old.place" "Old" "R.Old")
(sdedr:define-constant-profile "New" "PhosphorusActiveConcentration" 20)
(sdedr:define-constant-profile-region "New.place" "New" "R.New")
"""
        aba = self._ir("""
(sdegeo:set-default-boolean "ABA")
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Old")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "R.New")
""" + common)
        self.assertEqual((20.0, 0.0), evaluate_m2_doping_at(aba, 1.0, 0.5)[:2])

        bab = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Old")
(sdegeo:set-default-boolean "BAB")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "R.New")
""" + common)
        self.assertEqual((0.0, 10.0), evaluate_m2_doping_at(bab, 1.0, 0.5)[:2])

    def test_constant_decay_and_refinement_fail_closed(self) -> None:
        decay = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdedr:define-constant-profile "D" "PhosphorusActiveConcentration" 1e15)
(sdedr:define-constant-profile-region "D.place" "D" "R.Si" 0.1 "Gauss")
""")
        with self.assertRaises(GmshBackendError) as context:
            validate_m2_device_ir(decay)
        self.assertEqual("constant-profile-decay", context.exception.details["capability"])

        refined = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdedr:define-refinement-size "Ref" 0.2 0.2 0.1 0.1)
(sdedr:define-refinement-region "Ref.place" "Ref" "R.Si")
""")
        with self.assertRaises(GmshBackendError) as context:
            validate_m2_device_ir(refined)
        self.assertEqual("refinement", context.exception.details["capability"])

    def test_vela_contact_owner_uses_curve_adjacency_not_endpoint_triangles(self) -> None:
        # Node 1 belongs to triangles on both sides of an interface.  The
        # contact curve itself is adjacent only to R.Left, which is the owner
        # already audited by the Gmsh backend.
        bundle = {
            "nodes": [
                {"vela_id": 0, "x": 0.0, "y": 1.0},
                {"vela_id": 1, "x": 1.0, "y": 1.0},
                {"vela_id": 2, "x": 0.0, "y": 0.0},
                {"vela_id": 3, "x": 2.0, "y": 1.0},
            ],
            "triangles": [
                {"vela_id": 0, "region_id": 0, "material_id": 0,
                 "node_ids": [0, 2, 1]},
                {"vela_id": 1, "region_id": 1, "material_id": 0,
                 "node_ids": [1, 2, 3]},
            ],
            "regions": [
                {"id": 0, "name": "R.Left"}, {"id": 1, "name": "R.Right"},
            ],
            "materials": [{"id": 0, "name": "Silicon"}],
            "contacts": [{"id": 0, "name": "TopLeft", "physical_tag": 3001}],
            "lines": [{"node_ids": [0, 1], "physical_curve_tags": [3001]}],
            "provenance": {
                "material_aliases": {"Silicon": "Si"},
                "contact_owners": {"TopLeft": ["R.Left"]},
            },
        }
        mesh = bundle_to_vela_mesh(bundle)
        self.assertEqual(0, mesh["contacts"][0]["region_id"])


class M2GmshIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        prefix = Path(sys.prefix)
        cls.wrapper = prefix / "lib" / "gmsh.py"
        cls.library = prefix / "bin" / "libgmsh.dll"
        cls.executable = prefix / "bin" / "gmsh.exe"
        if not all(path.is_file() for path in (cls.wrapper, cls.library, cls.executable)):
            raise unittest.SkipTest("Gmsh 4.15.2 UCRT64 runtime is not installed")

    def _run(self, source: Path, output: Path, runtime: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable, str(REPO / "scripts" / "sentaurus_gmsh_pipeline.py"),
                "--sde", str(source), "--output-dir", str(output),
                "--runtime-dir", str(runtime), "--gmsh-python", str(self.wrapper),
                "--gmsh-library", str(self.library),
                "--gmsh-executable", str(self.executable), "--mesh-size", "0.25",
                "--capability-level", "M2",
            ],
            cwd=REPO, text=True, capture_output=True, check=False,
        )

    def test_core_fixture_roundtrips_m2_semantics(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m2_core_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = self._run(FIXTURE, output, root / "runtime")
            self.assertEqual(0, completed.returncode, completed.stderr)
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("generated", manifest["status"])
            self.assertEqual("M2", manifest["milestone"])
            self.assertEqual("passed", manifest["cross_format"]["status"])
            bundle = json.loads((output / "mesh_field_bundle.json").read_text(encoding="utf-8"))
            audit = bundle["provenance"]["geometry_audit"]
            self.assertEqual({"R.N": 2.0, "R.P": 1.0}, audit["region_areas_um2"])
            self.assertEqual({"interface:R.P+R.N"},
                             {item["physical_name"] for item in bundle["interfaces"]})
            self.assertEqual({"Anode": ["R.P"], "Cathode": ["R.N"]},
                             bundle["provenance"]["contact_owners"])
            self.assertTrue(any(item["kind"] == "interface" for item in bundle["lines"]))
            interface_nodes = [item["vela_id"] for item in bundle["nodes"]
                               if math.isclose(item["x"], 1.0, abs_tol=1.0e-12)]
            self.assertTrue(interface_nodes)
            donors = bundle["fields"]["DonorConcentration"]["values"]
            acceptors = bundle["fields"]["AcceptorConcentration"]["values"]
            for node in interface_nodes:
                self.assertTrue(math.isclose(donors[node], 2.01e16, rel_tol=1.0e-9))
                self.assertTrue(math.isclose(acceptors[node], 1.0e16, rel_tol=1.0e-12))
            anode_tag = next(item["physical_tag"] for item in bundle["contacts"]
                             if item["name"] == "Anode")
            anode_nodes = {node for line in bundle["lines"]
                           if anode_tag in line["physical_curve_tags"]
                           for node in line["node_ids"]}
            anode_y = [bundle["nodes"][node]["y"] for node in anode_nodes]
            self.assertAlmostEqual(0.0, min(anode_y))
            self.assertAlmostEqual(0.5, max(anode_y))
            self.assertTrue((output / "device.msh").is_file())
            self.assertTrue((output / "device.vtu").is_file())
            self.assertTrue((output / "device.vtk").is_file())
            self.assertTrue((output / "device.grd").is_file())
            self.assertTrue((output / "device.dat").is_file())

    def test_boolean_fragment_area_and_interface_invariants(self) -> None:
        cases = {
            "partial_bab": ("""
(sdegeo:create-rectangle (position 0 0 0) (position 2 1 0) "Silicon" "Old")
(sdegeo:set-default-boolean "BAB")
(sdegeo:create-rectangle (position 1 0 0) (position 3 1 0) "Silicon" "New")
""", {"Old": 2.0, "New": 1.0}, 1),
            "contained_aba": ("""
(sdegeo:create-rectangle (position 0 0 0) (position 3 1 0) "Silicon" "Old")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "New")
""", {"Old": 2.0, "New": 1.0}, 1),
            "contained_bab": ("""
(sdegeo:create-rectangle (position 0 0 0) (position 3 1 0) "Silicon" "Old")
(sdegeo:set-default-boolean "BAB")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "New")
""", {"Old": 3.0}, 0),
            "multiple_old_aba": ("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "Left")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "Right")
(sdegeo:create-rectangle (position 0.5 0 0) (position 1.5 1 0) "Silicon" "New")
""", {"Left": 0.5, "New": 1.0, "Right": 0.5}, 2),
            "edge_touch": ("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "Left")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "Right")
""", {"Left": 1.0, "Right": 1.0}, 1),
            "point_touch": ("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "Left")
(sdegeo:create-rectangle (position 1 1 0) (position 2 2 0) "Silicon" "Right")
""", {"Left": 1.0, "Right": 1.0}, 0),
        }
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m2_bool_") as tmp:
            root = Path(tmp)
            for name, (text, expected_areas, interface_count) in cases.items():
                with self.subTest(name=name):
                    source = root / f"{name}.cmd"
                    source.write_text(text, encoding="utf-8")
                    output = root / f"{name}-output"
                    completed = self._run(source, output, root / f"{name}-runtime")
                    self.assertEqual(0, completed.returncode, completed.stderr)
                    bundle = json.loads((output / "mesh_field_bundle.json").read_text(
                        encoding="utf-8"))
                    actual = bundle["provenance"]["geometry_audit"]["region_areas_um2"]
                    self.assertEqual(set(expected_areas), set(actual))
                    for region, area in expected_areas.items():
                        self.assertTrue(math.isclose(area, actual[region], abs_tol=1.0e-10))
                    self.assertEqual(interface_count, len(bundle["interfaces"]))

    def test_endpoint_contact_fails_closed_with_source_location(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m2_contact_fail_") as tmp:
            root = Path(tmp)
            source = root / "endpoint.cmd"
            source.write_text("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdegeo:define-contact-set "Bad")
(sdegeo:set-contact (find-edge-id (position 0 0 0)) "Bad")
""", encoding="utf-8")
            output = root / "output"
            completed = self._run(source, output, root / "runtime")
            self.assertNotEqual(0, completed.returncode)
            report = json.loads((output / "fail_report.json").read_text(encoding="utf-8"))
            self.assertEqual("gmsh_backend", report["stage"])
            self.assertIn("does not select one curve interior", report["message"])
            self.assertRegex(report["details"]["candidates"].__repr__(), r"True")
            self.assertTrue(report["details"]["source_location"].endswith(":4:1"))


if __name__ == "__main__":
    unittest.main()
