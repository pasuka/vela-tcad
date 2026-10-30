"""P1 corpus contracts, native identity, reproducibility and safe failed reruns."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import sentaurus_sde_cases as cases
import sentaurus_gmsh_pipeline as pipeline
from jsonschema import ValidationError, validate
from tests.regression.gmsh_test_runtime import require_gmsh_runtime
from vela_sdevice.mesh_import import import_grid


class CaseManifestTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="vela P1 contract ")
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.data = json.loads(cases.DEFAULT_MANIFEST.read_text())

    def altered(self):
        p = self.root / "cases.json"
        p.write_text(json.dumps(self.data))
        return p

    def test_manifest_matches_all_current_primary_sources(self):
        data = cases.load_manifest(cases.DEFAULT_MANIFEST)
        actual = {p.relative_to(REPO).as_posix() for p in REPO.glob("reference_tcad/*/source/*")
                  if p.suffix in {".cmd", ".scm"} and "(sdegeo:create-" in p.read_text(encoding="utf-8")}
        self.assertEqual(actual, {c["source"] for c in data["cases"]})
        self.assertEqual(7, sum(c["expected"] == "generated" for c in data["cases"]))
        self.assertEqual(1, sum(c["expected"] == "rejected" for c in data["cases"]))

    def test_source_edit_cannot_silently_change_case(self):
        self.data["cases"][0]["source_sha256_lf"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "source hash mismatch"):
            cases.load_manifest(self.altered())

    def test_nested_contract_schema_is_enforced(self):
        del self.data["cases"][0]["contract"]["contacts"]["Anode"]["length_um"]
        with self.assertRaises(ValidationError):
            cases.load_manifest(self.altered())

    def test_material_mapping_has_its_own_identity(self):
        self.data["cases"][0]["contract"]["materials_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "material contract hash mismatch"):
            cases.load_manifest(self.altered())

    def test_duplicate_ids_are_rejected(self):
        self.data["cases"][1]["id"] = self.data["cases"][0]["id"]
        with self.assertRaisesRegex(ValueError, "duplicate case"):
            cases.load_manifest(self.altered())

    def test_invalid_rerun_invalidates_old_batch_success(self):
        (self.root / "case_summary.json").write_text(json.dumps({"schema": "vela.reference_sde_run.v1", "status": "passed"}))
        self.data["cases"][0]["source_sha256_lf"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "source hash mismatch"):
            cases.run_suite(self.altered(), self.root, self.root / "runtime")
        self.assertFalse((self.root / "case_summary.json").exists())


class OutputOwnershipTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="vela P1 ownership ")
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.source = self.root / "input.cmd"
        self.source.write_text("owned source")
        self.out = self.root / "output"
        self.out.mkdir()

    def test_unrelated_reserved_filename_is_preserved(self):
        original = self.out / "device.grd"
        original.write_bytes(b"user's grid")
        with self.assertRaisesRegex(pipeline.PipelineError, "unrelated file"):
            pipeline._prepare_output(self.out, self.source)
        self.assertEqual(b"user's grid", original.read_bytes())

    def test_input_output_collision_preserves_source(self):
        p = self.out / "device.grd"
        p.write_bytes(self.source.read_bytes())
        with self.assertRaisesRegex(pipeline.PipelineError, "input source collides"):
            pipeline._prepare_output(self.out, p)
        self.assertEqual(self.source.read_bytes(), p.read_bytes())

    def test_unrelated_temporary_filename_is_preserved(self):
        note = self.out / "fail_report.json.tmp"
        note.write_text("user data")
        pipeline._strict_json(self.out / "fail_report.json", {"status": "failed"})
        self.assertEqual("user data", note.read_text())

    def test_runtime_staging_cannot_overwrite_input(self):
        args = pipeline.build_parser().parse_args([
            "--sde", str(self.source), "--output-dir", str(self.out),
            "--runtime-dir", str(self.root)])
        with self.assertRaisesRegex(pipeline.PipelineError, "runtime staging"):
            pipeline.run_pipeline(args)
        self.assertEqual("owned source", self.source.read_text())

    def test_manifest_cannot_delete_outside_output(self):
        manifest = self.out / "run_manifest.json"
        manifest.write_text(json.dumps({"schema": pipeline.RUN_SCHEMA,
            "generated": ["../input.cmd"], "sha256": {"../input.cmd": cases.digest(self.source)}}))
        with self.assertRaisesRegex(pipeline.PipelineError, "ownership manifest"):
            pipeline._prepare_output(self.out, self.source)
        self.assertEqual("owned source", self.source.read_text())
        self.assertTrue(manifest.is_file())


class CorpusIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require_gmsh_runtime()
        cls.temp = tempfile.TemporaryDirectory(prefix="vela P1 full corpus ")
        cls.addClassCleanup(cls.temp.cleanup)
        cls.root = Path(cls.temp.name)
        cls.runtime = cls.root / "runtime"
        cls.output = Path(os.environ.get("VELA_REFERENCE_SDE_TEST_OUTPUT", cls.root / "case outputs with spaces"))
        cls.data = cases.load_manifest(cases.DEFAULT_MANIFEST)
        cls.summary = cases.run_suite(cases.DEFAULT_MANIFEST, cls.output, cls.runtime, repeat=2)

    def cli(self, source, output):
        return subprocess.run([sys.executable, str(REPO / "scripts/sentaurus_gmsh_pipeline.py"),
            "--sde", str(source), "--output-dir", str(output), "--runtime-dir", str(self.runtime), "--parameter", "node=0"],
            cwd=REPO, text=True, capture_output=True)

    def test_seven_exports_native_identity_and_one_rejection(self):
        self.assertEqual("passed", self.summary["status"], self.summary)
        self.assertEqual({"generated": 7, "expected_rejection": 1}, self.summary["counts"])
        for row in self.summary["cases"]:
            with self.subTest(case=row["id"]):
                self.assertEqual("not_run", row["stages"]["reference"])
                if row["status"] == "generated":
                    self.assertEqual("passed", row["stages"]["native_import"])
                    self.assertEqual(2, row["repeats"])
                else:
                    self.assertEqual("not_run", row["stages"]["export"])

    def test_bjt_local_variant_refines_authored_local_windows(self):
        coarse = json.loads((self.output / "bjt/mesh_qualification.json").read_text())
        local = json.loads((self.output / "bjt_local/mesh_qualification.json").read_text())
        bounds = {"RefPlace.BaseCollector": [.75,.42,5.25,.88],
                  "RefPlace.EmitterBaseLeft": [2.52,0,2.98,.55],
                  "RefPlace.EmitterBaseRight": [4.02,0,4.48,.55]}
        # Compare both meshes in the same authored windows, independently of field tags.
        for name, (x0,y0,x1,y1) in bounds.items():
            spacings=[]
            for case in ("bjt", "bjt_local"):
                b=json.loads((self.output / case / "mesh_field_bundle.json").read_text())
                nodes={n["vela_id"]:n for n in b["nodes"]}
                edges={tuple(sorted((t["node_ids"][i],t["node_ids"][(i+1)%3]))) for t in b["triangles"] for i in range(3)}
                lengths=sorted(math.hypot(nodes[a]["x"]-nodes[c]["x"],nodes[a]["y"]-nodes[c]["y"]) for a,c in edges
                    if x0 <= (nodes[a]["x"]+nodes[c]["x"])/2 <= x1 and y0 <= (nodes[a]["y"]+nodes[c]["y"])/2 <= y1)
                self.assertTrue(lengths)
                spacings.append(lengths[math.ceil(.95*len(lengths))-1])
            self.assertLess(spacings[1],spacings[0],(name,spacings))
        for q in (coarse,local):
            validate(q["maxtransdiff"], json.loads((REPO/"schemas/vela.gmsh_maxtransdiff.v1.schema.json").read_text()))
            self.assertEqual("passed",q["maxtransdiff"]["status"])
            self.assertGreater(q["maxtransdiff"]["remesh_iterations"],0)

    def test_bjt_integral_and_junction_corruption_are_rejected(self):
        from sentaurus_bjt_contract import check
        b=json.loads((self.output/"bjt/mesh_field_bundle.json").read_text())
        b["fields"]["AcceptorConcentration"]["values"]=[v*1.05 for v in b["fields"]["AcceptorConcentration"]["values"]]
        self.assertIn("integrated concentration",check(b)["failures"])
        b=json.loads((self.output/"bjt/mesh_field_bundle.json").read_text())
        b["fields"]["DopingConcentration"]["values"]=[1e17]*len(b["nodes"])
        self.assertTrue(any("junction count" in f for f in check(b)["failures"]))

    def test_analytic_contract_detects_doping_corruption(self):
        out = self.root / "corrupted field"
        shutil.copytree(self.output / "pn2d", out)
        path = out / "mesh_field_bundle.json"
        bundle = json.loads(path.read_text())
        bundle["fields"]["DonorConcentration"]["values"][0] += 1e16
        path.write_text(json.dumps(bundle))
        with self.assertRaisesRegex(ValueError, "analytic doping mismatch"):
            cases.check_contract(out, self.data["cases"][0]["contract"])

    def test_mesh_size_contract_detects_under_refinement(self):
        contract = json.loads(json.dumps(self.data["cases"][0]["contract"]))
        contract["mesh_quality"]["region_max_spacing_um"]["R.Si"] = [0.0001, 0.0001]
        with self.assertRaisesRegex(ValueError, "P95 spacing"):
            cases.check_contract(self.output / "pn2d", contract)

    def test_consumer_reversed_triangle_is_detected(self):
        out = self.root / "native reversed"
        import_grid(self.output / "gate_reference/device.grd", out)
        path = out / "mesh.json"
        mesh = json.loads(path.read_text())
        mesh["triangles"][0]["node_ids"].reverse()
        path.write_text(json.dumps(mesh))
        with self.assertRaisesRegex(ValueError, "triangle orientation"):
            cases.compare_native(self.output / "gate_reference", out)

    def test_consumer_contact_corruption_is_detected(self):
        out = self.root / "native mutated"
        import_grid(self.output / "gate_reference/device.grd", out)
        path = out / "mesh.json"
        mesh = json.loads(path.read_text())
        mesh["contacts"][0]["edge_node_ids"].pop()
        path.write_text(json.dumps(mesh))
        with self.assertRaisesRegex(ValueError, "contact edges/nodes"):
            cases.compare_native(self.output / "gate_reference", out)

    def test_success_failure_success_keeps_inputs_and_unrelated_files(self):
        output = self.root / "retry output"
        source = self.root / "source with spaces.cmd"
        original = (REPO / self.data["cases"][0]["source"]).read_bytes()
        source.write_bytes(original)
        self.assertEqual(0, self.cli(source, output).returncode)
        expected_hashes = {name: cases.digest(output / name) for name in cases.PRIMARY}
        note = output / "user-notes.txt"
        note.write_text("keep this")
        failed_source = REPO / self.data["cases"][7]["source"]
        result = self.cli(failed_source, output)
        self.assertEqual(2, result.returncode, result.stderr)
        self.assertFalse((output / "run_manifest.json").exists())
        self.assertFalse(any((output / name).exists() for name in pipeline.GENERATED_FILES))
        self.assertEqual("sde_parse", json.loads((output / "fail_report.json").read_text())["stage"])
        self.assertEqual(0, self.cli(source, output).returncode)
        self.assertEqual(expected_hashes, {name: cases.digest(output / name) for name in cases.PRIMARY})
        self.assertFalse((output / "fail_report.json").exists())
        self.assertEqual(original, source.read_bytes())
        self.assertEqual("keep this", note.read_text())

    def test_user_modified_previous_grid_is_not_overwritten(self):
        output = self.root / "modified previous output"
        source = REPO / self.data["cases"][2]["source"]
        self.assertEqual(0, self.cli(source, output).returncode)
        before = (output / "run_manifest.json").read_bytes()
        grid = output / "device.grd"
        grid.write_bytes(b"user modifications")
        result = self.cli(source, output)
        self.assertEqual(2, result.returncode)
        self.assertIn("output_safety", result.stderr)
        self.assertEqual(b"user modifications", grid.read_bytes())
        self.assertEqual(before, (output / "run_manifest.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
