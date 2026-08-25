"""M3 refinement field DAG and complete PN2D preprocessing tests."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import SdeParseError, parse_sde_device_ir  # noqa: E402
from sentaurus_gmsh_backend import GmshBackendError, validate_m3_device_ir  # noqa: E402
from sentaurus_ir_v2 import migrate_device_ir_v1  # noqa: E402

FIXTURE = REPO / "tests" / "fixtures" / "sentaurus_gmsh_m3" / "m3_refinement_sde.cmd"
PN2D = REPO / "reference_tcad" / "pn2d_sentaurus2022_oracle" / "input" / "pn2d_sde.cmd"


def parse_ir(path: Path) -> dict:
    return migrate_device_ir_v1(parse_sde_device_ir(path))


class M3RefinementUnitTest(unittest.TestCase):
    def _ir(self, text: str) -> dict:
        with tempfile.TemporaryDirectory(prefix="vela_m3_ir_") as tmp:
            path = Path(tmp) / "case.cmd"
            path.write_text(text, encoding="utf-8")
            return parse_ir(path)

    def test_fixture_freezes_polygon_and_maxlenint_semantics(self) -> None:
        ir = parse_ir(FIXTURE)
        validate_m3_device_ir(ir)
        self.assertTrue(ir["capability_manifest"]["provenance_complete"])
        windows = {item["name"]: item for item in ir["mesh_control"]["windows"]}
        self.assertEqual("polygon", windows["Polygon.Window"]["shape"])
        self.assertEqual(4, len(windows["Polygon.Window"]["points"]))
        function = ir["mesh_control"]["refinement_functions"][0]
        self.assertEqual("MaxLenInt", function["function"])
        self.assertEqual(["R.Si", "R.Ox"], function["interface"])
        self.assertEqual(0.03, function["value"])
        self.assertEqual(1.5, function["factor"])
        self.assertTrue(function["double_side"])
        self.assertTrue(function["use_region_names"])

    def test_maxlenint_default_factor_and_material_names(self) -> None:
        ir = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Oxide" "R.Ox")
(sdedr:define-refeval-window "W" "Rectangle" (position 0 0 0) (position 2 1 0))
(sdedr:define-refinement-size "Ref" 0.2 0.2 0.05 0.05)
(sdedr:define-refinement-placement "Ref.Place" "Ref" "W")
(sdedr:define-refinement-function "Ref" "MaxLenInt" "Silicon" "Oxide" 0.02)
""")
        validate_m3_device_ir(ir)
        function = ir["mesh_control"]["refinement_functions"][0]
        self.assertIsNone(function["factor"])
        self.assertFalse(function["double_side"])
        self.assertFalse(function["use_region_names"])

    def test_phase4b_function_is_parsed_then_fails_closed(self) -> None:
        ir = self._ir("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdedr:define-refinement-size "Ref" 0.2 0.2 0.05 0.05)
(sdedr:define-refinement-region "Ref.Place" "Ref" "R.Si")
(sdedr:define-refinement-function
  "Ref" "DopingConcentration" "MaxTransDiff" 1)
""")
        function = ir["mesh_control"]["refinement_functions"][0]
        self.assertEqual("DopingConcentration", function["variable"])
        self.assertEqual("MaxTransDiff", function["function"])
        with self.assertRaises(GmshBackendError) as context:
            validate_m3_device_ir(ir)
        self.assertEqual("refinement-4b", context.exception.details["capability"])

    def test_polygon_requires_closed_counterclockwise_vertices(self) -> None:
        for points in (
            "(list (position 0 0 0) (position 1 0 0) (position 0 1 0))",
            "(list (position 0 0 0) (position 0 1 0) (position 1 0 0) (position 0 0 0))",
        ):
            with self.subTest(points=points):
                with tempfile.TemporaryDirectory(prefix="vela_m3_polygon_fail_") as tmp:
                    path = Path(tmp) / "bad.cmd"
                    path.write_text(
                        f'(sdedr:define-refeval-window "P" "Polygon" {points})\n',
                        encoding="utf-8",
                    )
                    with self.assertRaises(SdeParseError):
                        parse_sde_device_ir(path)


class M3GmshIntegrationTest(unittest.TestCase):
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
                "--capability-level", "M3",
            ],
            cwd=REPO, text=True, capture_output=True, check=False,
        )

    def test_refinement_fixture_roundtrips_field_dag_and_maxlenint(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m3_refinement_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = self._run(FIXTURE, output, root / "runtime")
            self.assertEqual(0, completed.returncode, completed.stderr)
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            qualification = manifest["mesh_qualification"]["refinement"]
            self.assertEqual("M3", manifest["milestone"])
            self.assertEqual("passed", manifest["cross_format"]["status"])
            self.assertEqual("passed", qualification["status"])
            self.assertEqual(4, len(qualification["placements"]))
            self.assertTrue(all(item["passed"] for item in qualification["placements"]))
            self.assertEqual(1, len(qualification["max_len_int"]))
            self.assertTrue(qualification["max_len_int"][0]["passed"])
            reports = {item["name"]: item for item in qualification["placements"]}
            self.assertLess(
                reports["Line.Ref.Place"]["p95_edge_component"][0],
                reports["Polygon.Ref.Place"]["p95_edge_component"][0],
            )
            dag = json.loads((output / "refinement_field_dag.json").read_text(
                encoding="utf-8"))
            self.assertEqual("MinAniso", dag["nodes"][-1]["kind"])
            self.assertEqual("Gmsh MathEvalAniso/Distance/MinAniso with BAMG", dag["engine"])
            self.assertTrue(any(item["id"].startswith("maxlenint:") and
                                item["kind"] == "MathEvalAniso" for item in dag["nodes"]))

    def test_complete_pn2d_preprocessing_meets_cross_format_and_quality_contracts(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m3_pn2d_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = self._run(PN2D, output, root / "runtime")
            self.assertEqual(0, completed.returncode, completed.stderr)
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            quality = manifest["mesh_qualification"]
            self.assertEqual("passed", manifest["cross_format"]["status"])
            self.assertEqual("passed", quality["refinement"]["status"])
            self.assertGreater(quality["node_count"], 900)
            self.assertGreater(quality["cell_count"], 1800)
            self.assertEqual(1, quality["connected_components"])
            placements = {item["name"]: item for item in quality["refinement"]["placements"]}
            self.assertEqual([0.01, 0.02],
                             placements["Junction.Mesh.Place"]["requested_max_spacing"])
            self.assertTrue(placements["Junction.Mesh.Place"]["passed"])
            bundle = json.loads((output / "mesh_field_bundle.json").read_text(encoding="utf-8"))
            self.assertEqual(7, bundle["gmsh"]["options"]["Mesh.Algorithm"])
            self.assertEqual(1000, bundle["gmsh"]["options"]["Mesh.AnisoMax"])
            self.assertEqual({"Anode", "Cathode"},
                             {item["name"] for item in bundle["contacts"]})


if __name__ == "__main__":
    unittest.main()
