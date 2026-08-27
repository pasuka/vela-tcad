"""M5 selectable Transfinite Rectangle mesh-policy regression tests."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import sentaurus_gmsh_backend as backend  # noqa: E402
import sentaurus_gmsh_pipeline as pipeline  # noqa: E402
from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402
from sentaurus_ir_v2 import migrate_device_ir_v1  # noqa: E402
from tests.regression.gmsh_test_runtime import require_gmsh_runtime  # noqa: E402


PN2D = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sde.cmd"
M2_FIXTURE = REPO / "tests" / "fixtures" / "sentaurus_gmsh_m2" / "m2_core_sde.cmd"


class M5PolicyUnitTest(unittest.TestCase):
    def test_pipeline_defaults_to_auto_alternate_right(self) -> None:
        args = pipeline.build_parser().parse_args([
            "--sde", str(PN2D), "--output-dir", "out",
        ])
        self.assertEqual("M5", args.capability_level)
        self.assertEqual("auto", args.mesh_strategy)
        self.assertEqual("AlternateRight", args.transfinite_diagonal)

    def test_eligibility_is_fail_closed_and_auditable(self) -> None:
        device_ir = migrate_device_ir_v1(parse_sde_device_ir(PN2D))
        eligible = backend._m5_transfinite_eligibility(device_ir)
        self.assertTrue(eligible["eligible"])
        modified = json.loads(json.dumps(device_ir))
        modified["geometry"]["bodies"].append(
            {**modified["geometry"]["bodies"][0], "id": "body:1",
             "creation_order": 1})
        rejected = backend._m5_transfinite_eligibility(modified)
        self.assertFalse(rejected["eligible"])
        self.assertIn("body_count", {item["code"] for item in rejected["reasons"]})


class M5GmshIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        sources = require_gmsh_runtime()
        cls.wrapper = sources.wrapper
        cls.library = sources.library
        cls.executable = sources.executable

    def test_pn2d_auto_generates_exact_right_tri3_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m5_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = subprocess.run([
                sys.executable, str(REPO / "scripts" / "sentaurus_gmsh_pipeline.py"),
                "--sde", str(PN2D), "--output-dir", str(output),
                "--runtime-dir", str(root / "runtime"),
                "--gmsh-python", str(self.wrapper),
                "--gmsh-library", str(self.library),
                "--gmsh-executable", str(self.executable),
                "--mesh-size", "0.2", "--capability-level", "M5",
                "--mesh-strategy", "auto",
                "--transfinite-diagonal", "AlternateRight",
            ], cwd=REPO, text=True, capture_output=True, check=False)
            self.assertEqual(0, completed.returncode, completed.stderr)
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            qualification = manifest["mesh_qualification"]
            contract = json.loads((output / "transfinite_mesh_contract.json").read_text(
                encoding="utf-8"))
            self.assertEqual("M5", manifest["milestone"])
            self.assertEqual("transfinite-rect", qualification["mesh_strategy"])
            self.assertEqual("mixed_voronoi", qualification["node_volume_policy"])
            self.assertEqual(1482, qualification["node_count"])
            self.assertEqual(2800, qualification["cell_count"])
            self.assertTrue(qualification["transfinite"]["all_triangles_right"])
            self.assertTrue(qualification["transfinite"]["non_obtuse"])
            self.assertEqual([18, 20, 18], contract["x_intervals"])
            self.assertEqual([25], contract["y_intervals"])
            self.assertEqual("passed", manifest["cross_format"]["status"])

    def test_explicit_transfinite_rejects_ineligible_multibody_deck(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m5_fail_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = subprocess.run([
                sys.executable, str(REPO / "scripts" / "sentaurus_gmsh_pipeline.py"),
                "--sde", str(M2_FIXTURE), "--output-dir", str(output),
                "--runtime-dir", str(root / "runtime"),
                "--gmsh-python", str(self.wrapper),
                "--gmsh-library", str(self.library),
                "--gmsh-executable", str(self.executable),
                "--capability-level", "M5", "--mesh-strategy", "transfinite-rect",
            ], cwd=REPO, text=True, capture_output=True, check=False)
            self.assertEqual(2, completed.returncode)
            report = json.loads((output / "fail_report.json").read_text(encoding="utf-8"))
            self.assertEqual("failed", report["status"])
            self.assertEqual("mesh_strategy", report["details"]["stage"])
            self.assertFalse(report["details"]["eligibility"]["eligible"])
            self.assertFalse((output / "run_manifest.json").exists())

    def test_partial_window_refines_axis_projection_and_preserves_full_contact(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m5_partial_") as tmp:
            root = Path(tmp)
            output = root / "output"
            source = root / "partial_window.cmd"
            source.write_text(
                PN2D.read_text(encoding="utf-8").replace(
                    "(position 0.9 0.0 0.0)", "(position 0.9 0.25 0.0)", 1),
                encoding="utf-8",
            )
            completed = subprocess.run([
                sys.executable, str(REPO / "scripts" / "sentaurus_gmsh_pipeline.py"),
                "--sde", str(source), "--output-dir", str(output),
                "--runtime-dir", str(root / "runtime"),
                "--gmsh-python", str(self.wrapper),
                "--gmsh-library", str(self.library),
                "--gmsh-executable", str(self.executable),
                "--mesh-size", "0.2", "--capability-level", "M5",
                "--mesh-strategy", "auto",
            ], cwd=REPO, text=True, capture_output=True, check=False)
            self.assertEqual(0, completed.returncode, completed.stderr)
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            contract = json.loads((output / "transfinite_mesh_contract.json").read_text(
                encoding="utf-8"))
            contacts = {item["name"]: item["node_count"]
                        for item in manifest["mesh_qualification"]["contacts"]}
            self.assertEqual([5, 13], contract["y_intervals"])
            self.assertEqual(19, contacts["Anode"])
            self.assertEqual(19, contacts["Cathode"])

    def test_all_supported_diagonal_arrangements_remain_right_tri3(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_m5_diagonals_") as tmp:
            root = Path(tmp)
            for diagonal in ("Left", "Right", "AlternateLeft", "AlternateRight"):
                with self.subTest(diagonal=diagonal):
                    output = root / diagonal
                    completed = subprocess.run([
                        sys.executable,
                        str(REPO / "scripts" / "sentaurus_gmsh_pipeline.py"),
                        "--sde", str(PN2D), "--output-dir", str(output),
                        "--runtime-dir", str(root / f"runtime-{diagonal}"),
                        "--gmsh-python", str(self.wrapper),
                        "--gmsh-library", str(self.library),
                        "--gmsh-executable", str(self.executable),
                        "--mesh-size", "0.2", "--capability-level", "M5",
                        "--mesh-strategy", "transfinite-rect",
                        "--transfinite-diagonal", diagonal,
                    ], cwd=REPO, text=True, capture_output=True, check=False)
                    self.assertEqual(0, completed.returncode, completed.stderr)
                    qualification = json.loads(
                        (output / "mesh_qualification.json").read_text(encoding="utf-8"))
                    self.assertEqual(diagonal, qualification["transfinite"]["diagonal"])
                    self.assertTrue(qualification["transfinite"]["all_triangles_right"])
                    self.assertTrue(qualification["transfinite"]["non_obtuse"])


if __name__ == "__main__":
    unittest.main()
