"""Contract tests for the ``vela device`` (SDevice cmd direct-run) work packages.

WP0 freezes the interface contracts C1, C2, C4, C7 and C10 as JSON schemas, the
message catalogue in ``docs/sentaurus_sdevice_workflow.md`` and the golden IR
fixtures under ``tests/fixtures/sdevice``.  These tests pin that layer without
any parser or solver code:

* every schema is well formed, has a minimal valid sample and rejects an
  invalid one (checked with the stdlib validator and, when installed, with
  ``jsonschema`` for agreement);
* every golden IR validates against C1 and is internally consistent (hashes,
  line counts, verdicts, output names, message ids);
* the message catalogue has unique ids whose severity and exit-code column
  agree with the id itself.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from vela_sdevice import schema_validator as sv  # noqa: E402

SCHEMAS = REPO / "schemas"
FIXTURES = REPO / "tests" / "fixtures" / "sdevice"
WORKFLOW_DOC = REPO / "docs" / "sentaurus_sdevice_workflow.md"

C1 = SCHEMAS / "vela.sentaurus_execution_ir.v2.schema.json"
C2 = SCHEMAS / "vela.mesh_qualification.v1.schema.json"
C4 = SCHEMAS / "vela.stage_plan.v1.schema.json"
C7 = SCHEMAS / "vela.reference_artifact_manifest.v1.schema.json"
C10 = SCHEMAS / "vela.sdevice_sidecar.v1.schema.json"
ALL_SCHEMAS = (C1, C2, C4, C7, C10)

SHA256 = "0" * 64
CATALOGUE_ROW = re.compile(r"^\| (?P<id>[WE][1-8]\d{2}) \| (?P<severity>warning|fatal) \| (?P<exit>[0-8]|\u2014) \|")


def normalised_text(path: Path) -> str:
    return path.read_bytes().decode("utf-8-sig").replace("\r\n", "\n")


def normalised_sha256(path: Path) -> str:
    return hashlib.sha256(normalised_text(path).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def catalogue() -> dict:
    entries = {}
    for line in WORKFLOW_DOC.read_text(encoding="utf-8").splitlines():
        match = CATALOGUE_ROW.match(line)
        if match:
            entries[match["id"]] = {"severity": match["severity"], "exit": match["exit"]}
    return entries


def try_jsonschema():
    try:
        import jsonschema  # type: ignore[import-not-found]
    except ImportError:
        return None
    return jsonschema


# --- minimal samples -----------------------------------------------------------------

def minimal_ir() -> dict:
    return {
        "schema": "vela.sentaurus_execution_ir.v2",
        "source": {"path": "x_des.cmd", "sha256": SHA256, "line_count": 1, "syntax_version": "T-2022.03"},
        "verdict": {"accepted": True, "exit_code": 0, "fatal_ids": []},
        "preprocessor_findings": [],
        "errors": [],
        "warnings": [],
        "files": {
            "grid": {"raw": "x_msh.tdr", "resolved": "x_msh.tdr", "kind": "tdr", "line": 1},
            "parameter": {"raw": "m.par", "resolved": "m.par", "keyword": "Parameter", "line": 1},
            "current": {"raw": None, "resolved": "x_des.plt", "line": None},
            "plot": {"raw": None, "resolved": "x_des.tdr", "vela": "x_des.vtu", "line": None},
            "output": {"raw": None, "resolved": "output_des.log", "line": None},
            "ignored": [],
        },
        "electrodes": [{"name": "a", "voltage_V": 0.0, "attributes": {}, "line": 1}],
        "physics": {"temperature_K": 300.0, "global": [], "scoped": [], "effective": [], "fixed_charges": []},
        "plot": {"explicit": False, "fields": []},
        "current_plot": {"entries": []},
        "math": {
            "iterations": 50, "digits": 5, "rel_err_control": True,
            "err_ref_electron_cm3": None, "err_ref_hole_cm3": None,
            "extrapolate": False, "derivatives": True, "aval_derivatives": True,
            "compute_ionization_integrals": False, "break_at_ion_integral": None,
            "break_criteria": [], "exit_on_failure": False, "number_of_threads": 1, "ignored": [],
        },
        "solve": {"statements": [{"kind": "coupled", "line": 1, "equations": ["Poisson"], "iterations": None}]},
        "stages": [{
            "index": 0, "kind": "initial", "statement_index": 0, "equations": ["Poisson"],
            "iterations": None, "current_prefix": "",
            "outputs": {"current_file": "x_des.plt", "plot_files": [], "current_plot_times": None},
            "depends_on": [], "line": 1,
        }],
        "analysis": {"kind": "equilibrium", "reason": "R1: no Quasistationary statement",
                     "sweep_contacts": [], "final_contact": None, "final_voltage": None},
    }


def minimal_mesh_qualification() -> dict:
    return {
        "schema": "vela.mesh_qualification.v1",
        "source": {"grid_file": "x_msh.tdr", "format": "tdr", "sha256": SHA256},
        "dimension": 2,
        "coordinate_unit": {"declared": "um", "applied_scale_to_um": 1.0},
        "bounding_box_um": {"min": [0.0, 0.0], "max": [2.0, 0.5]},
        "nodes": {"count": 4},
        "elements": {"count": 2, "by_kind": {"triangle": 2}, "converted_quadrilaterals": 0},
        "regions": [{"name": "R.Si", "material": "Si", "element_count": 2}],
        "materials": ["Si"],
        "contacts": [{"name": "Anode", "node_count": 2, "edge_count": 1}],
        "angles": {"max_deg": 90.0, "obtuse_element_count": 0, "obtuse_fraction": 0.0},
        "doping": {"file": "doping.csv", "columns": ["node_id", "donors_cm3", "acceptors_cm3"],
                   "source_fields": ["DonorConcentration", "AcceptorConcentration"]},
        "outputs": {"mesh_json": "mesh.json", "doping_csv": "doping.csv"},
        "qualification": {"two_dimensional": True, "triangles_only": True, "coordinates_plausible": True,
                          "non_obtuse": True, "recommended_node_volume_policy": "barycentric",
                          "allow_obtuse": True},
    }


def minimal_stage_plan() -> dict:
    return {
        "schema": "vela.stage_plan.v1",
        "cmd": {"path": "x_des.cmd", "sha256": SHA256, "execution_ir": "execution_ir.json"},
        "work_dir": "x_des_vela",
        "exit_on_failure": False,
        "threads": 1,
        "inputs": {"mesh_json": "mesh.json", "doping_csv": "doping.csv", "materials_json": "materials.json",
                   "mesh_qualification": "mesh_qualification.json"},
        "outputs": {"log": "x_des.log", "final_vtu": "x_des.vtu", "pvd": "x_des.pvd",
                    "plot_final_only": False, "vtu_binary": False, "keep_intermediate": True},
        "electrodes": [{"name": "a", "contact": "a", "initial_voltage_V": 0.0}],
        "stages": [{
            "index": 0, "ir_stage_index": 0, "kind": "equilibrium", "deck": "simulation_stage00.json",
            "csv": "stage00.csv", "input_state": None, "output_state": "stage00_state.h5",
            "electrode_voltages_start_V": {"a": 0.0}, "electrode_voltages_goal_V": {"a": 0.0},
            "current_file": "x_des.plt", "plt_write": "all_accepted_steps",
            "vtu_write": {"per_accepted_step": True, "series": []}, "required_t": [],
        }],
    }


def minimal_manifest() -> dict:
    return {
        "schema": "vela.reference_artifact_manifest.v1",
        "case": "pn2d_sentaurus2018",
        "version": "v1",
        "created_utc": "2026-09-26T00:00:00Z",
        "store": {"kind": "ssh", "host": "sentaurus", "root": "/home/tcad/vela_reference_artifacts"},
        "sentaurus_version": "T-2022.03-SP2 (0.7745337)",
        "deck_matches_repository": True,
        "artifacts": [{
            "path": "plt/pn2d_0v.plt", "category": "B", "sha256": SHA256, "bytes": 1541,
            "producer": {"tool": "sdevice", "cmd": "pn2d_0v_sdevice.cmd", "cmd_sha256": SHA256},
            "verify": {"status": "registered"},
        }],
    }


def minimal_sidecar() -> dict:
    return {"schema": "vela.sdevice_sidecar.v1", "mesh_geometry": {"node_volume_policy": "mixed_voronoi"}}


class SchemaFileTest(unittest.TestCase):
    def test_schema_headers(self) -> None:
        for path in ALL_SCHEMAS:
            with self.subTest(schema=path.name):
                schema = sv.load_schema(path)
                self.assertEqual("http://json-schema.org/draft-07/schema#", schema["$schema"])
                self.assertEqual(path.name[: -len(".schema.json")], schema["$id"])
                self.assertEqual(schema["$id"], schema["properties"]["schema"]["const"])
                self.assertIs(False, schema["additionalProperties"])

    def test_minimal_samples_are_valid(self) -> None:
        samples = {C1: minimal_ir(), C2: minimal_mesh_qualification(), C4: minimal_stage_plan(),
                   C7: minimal_manifest(), C10: minimal_sidecar()}
        for path, sample in samples.items():
            with self.subTest(schema=path.name):
                self.assertEqual([], sv.collect_errors(sv.load_schema(path), sample))

    def test_invalid_samples_are_rejected(self) -> None:
        ir = minimal_ir()
        ir["verdict"]["exit_code"] = 4                      # not in the parse-level enum
        mesh = minimal_mesh_qualification()
        mesh["dimension"] = 4
        plan = minimal_stage_plan()
        plan["work_dir"] = "\u5de5\u4f5c\u76ee\u5f55_vela"        # ASCII only (CJK escaped)
        manifest = minimal_manifest()
        manifest["artifacts"][0]["path"] = "raw/pn2d_0v.plt"   # category directory not allowed
        sidecar = minimal_sidecar()
        sidecar["solver"] = {"mobility": {"model": "masetti"}}  # physics choice belongs to the cmd
        unknown_key = minimal_sidecar()
        unknown_key["electrodes"] = []
        cases = [(C1, ir), (C2, mesh), (C4, plan), (C7, manifest), (C10, sidecar), (C10, unknown_key)]
        for path, sample in cases:
            with self.subTest(schema=path.name):
                self.assertNotEqual([], sv.collect_errors(sv.load_schema(path), sample))

    def test_agreement_with_jsonschema_when_available(self) -> None:
        jsonschema = try_jsonschema()
        if jsonschema is None:
            self.skipTest("jsonschema not installed; stdlib validator verdicts not cross-checked")
        samples = [(C1, minimal_ir(), True), (C2, minimal_mesh_qualification(), True),
                   (C4, minimal_stage_plan(), True), (C7, minimal_manifest(), True),
                   (C10, minimal_sidecar(), True)]
        broken = minimal_sidecar()
        broken["solver"] = {"recombination": ["srh"]}
        samples.append((C10, broken, False))
        for path, sample, expected in samples:
            schema = sv.load_schema(path)
            with self.subTest(schema=path.name, expected=expected):
                self.assertEqual(expected, sv.is_valid(schema, sample))
                self.assertEqual(expected, jsonschema.Draft7Validator(schema).is_valid(sample))


class GoldenFixtureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.index = load_json(FIXTURES / "fixtures.json")
        cls.schema = sv.load_schema(C1)
        cls.catalogue = catalogue()

    def golden_for(self, entry: dict) -> dict:
        return load_json(FIXTURES / entry["golden"])

    def test_index_lists_existing_files_with_matching_hashes(self) -> None:
        for entry in self.index["fixtures"]:
            with self.subTest(deck=entry["deck"]):
                deck = FIXTURES / entry["deck"]
                origin = REPO / entry["origin"]
                self.assertTrue(deck.is_file(), deck)
                self.assertTrue((FIXTURES / entry["golden"]).is_file(), entry["golden"])
                self.assertEqual(entry["sha256"], normalised_sha256(deck))
                if origin.is_file():
                    self.assertEqual(normalised_text(origin), normalised_text(deck), "fixture drifted from origin")
                self.assertIn(entry["compare"], self.index["compare_modes"])
        for entry in self.index["parameter_files"]:
            with self.subTest(par=entry["path"]):
                self.assertTrue((FIXTURES / entry["path"]).is_file())

    def test_full_goldens_validate_and_are_consistent(self) -> None:
        for entry in self.index["fixtures"]:
            if entry["compare"] != "full":
                continue
            with self.subTest(deck=entry["deck"]):
                golden = self.golden_for(entry)
                deck = FIXTURES / entry["deck"]
                self.assertEqual([], sv.collect_errors(self.schema, golden))
                self.assertEqual(deck.name, golden["source"]["path"])
                self.assertEqual(entry["sha256"], golden["source"]["sha256"])
                self.assertEqual(len(normalised_text(deck).splitlines()), golden["source"]["line_count"])
                fatal = golden["preprocessor_findings"] + golden["errors"]
                self.assertEqual(not fatal, golden["verdict"]["accepted"])
                self.assertEqual(entry["expected_exit_code"], golden["verdict"]["exit_code"])
                ordered_ids = []
                for finding in fatal:
                    if finding["id"] not in ordered_ids:
                        ordered_ids.append(finding["id"])
                self.assertEqual(ordered_ids, golden["verdict"]["fatal_ids"])
                for finding in fatal:
                    self.assertEqual("fatal", finding["severity"])
                for finding in golden["warnings"]:
                    self.assertEqual("warning", finding["severity"])
                self.check_message_ids(golden)
                self.check_stages(golden)

    def check_message_ids(self, golden: dict) -> None:
        seen = []
        for finding in golden["preprocessor_findings"] + golden["errors"] + golden["warnings"]:
            seen.append((finding["id"], finding["severity"]))
        for item in golden["math"]["ignored"] + golden["files"]["ignored"]:
            seen.append((item["finding_id"], "warning"))
        for model in golden["physics"]["effective"]:
            if "finding_id" in model:
                seen.append((model["finding_id"], "warning" if model["status"] == "metadata_only" else "fatal"))
        for message_id, severity in seen:
            self.assertIn(message_id, self.catalogue, f"{message_id} missing from the message catalogue")
            self.assertEqual(self.catalogue[message_id]["severity"], severity, message_id)
        # every ignored Math keyword must also be reported as a warning finding
        warning_lines = {(w["id"], w["line"]) for w in golden["warnings"]}
        for item in golden["math"]["ignored"]:
            self.assertIn((item["finding_id"], item["line"]), warning_lines)

    def check_stages(self, golden: dict) -> None:
        statements = golden["solve"]["statements"]
        current_resolved = golden["files"]["current"]["resolved"]
        previous = None
        for stage in golden["stages"]:
            self.assertEqual(stage["index"], golden["stages"].index(stage))
            self.assertEqual([] if previous is None else [previous], stage["depends_on"])
            previous = stage["index"]
            statement = statements[stage["statement_index"]]
            self.assertIn(statement["kind"], {"poisson", "coupled", "quasistationary", "save", "load"})
            self.assertEqual(stage["current_prefix"] + current_resolved, stage["outputs"]["current_file"])
            if stage["kind"] == "sweep":
                self.assertEqual("quasistationary", statement["kind"])
                self.assertEqual(statement["goal"], stage["goal"])
                self.assertEqual(statement["t_control"], stage["t_control"])
            for plot in stage["outputs"]["plot_files"]:
                suffix = "_NNNN_des.tdr" if plot["time"] is not None else "_des.tdr"
                self.assertEqual(plot["file_prefix"] + suffix, plot["sentaurus_pattern"])
                self.assertEqual(plot["sentaurus_pattern"][: -len(".tdr")] + ".vtu", plot["vela_pattern"])
        sweeps = [s for s in golden["stages"] if s["kind"] == "sweep"]
        analysis = golden["analysis"]
        self.assertEqual([s["goal"]["contact"] for s in sweeps], analysis["sweep_contacts"])
        if sweeps:
            self.assertEqual(sweeps[-1]["goal"]["contact"], analysis["final_contact"])
            self.assertEqual(sweeps[-1]["goal"].get("voltage_V"), analysis["final_voltage"])
        else:
            self.assertEqual("equilibrium", analysis["kind"])
        if len(sweeps) > 1:
            self.assertEqual("multi_stage", analysis["kind"])

    def test_findings_goldens(self) -> None:
        for entry in self.index["fixtures"]:
            if entry["compare"] != "findings":
                continue
            with self.subTest(deck=entry["deck"]):
                golden = self.golden_for(entry)
                deck = FIXTURES / entry["deck"]
                self.assertEqual("vela.sdevice_fixture_findings.v1", golden["schema"])
                self.assertEqual(entry["sha256"], golden["source"]["sha256"])
                self.assertEqual(len(normalised_text(deck).splitlines()), golden["source"]["line_count"])
                self.assertFalse(golden["verdict"]["accepted"])
                self.assertEqual(entry["expected_exit_code"], golden["verdict"]["exit_code"])
                findings = golden.get("preprocessor_findings", []) + golden.get("errors", [])
                self.assertTrue(findings)
                ordered = []
                for finding in findings:
                    self.assertIn(finding["id"], self.catalogue)
                    self.assertEqual("fatal", self.catalogue[finding["id"]]["severity"])
                    self.assertEqual(str(entry["expected_exit_code"]), self.catalogue[finding["id"]]["exit"])
                    if finding["id"] not in ordered:
                        ordered.append(finding["id"])
                    lines = normalised_text(deck).splitlines()
                    self.assertIn(finding["text"], lines[finding["line"] - 1],
                                  f"{finding['id']} text not found on line {finding['line']}")
                self.assertEqual(ordered, golden["verdict"]["fatal_ids"])


class SchemaValidatorTest(unittest.TestCase):
    """Behavioural checks of the stdlib validator on the keywords the contracts rely on."""

    SCHEMA = {
        "type": "object",
        "additionalProperties": False,
        "required": ["id"],
        "properties": {
            "id": {"type": "string", "pattern": "^[WE][1-8][0-9]{2}$"},
            "count": {"type": "integer", "minimum": 1},
            "ratio": {"type": "number", "exclusiveMinimum": 0, "maximum": 1},
            "time": {"$ref": "#/definitions/time"},
            "flag": {"const": True},
            "kind": {"enum": ["a", "b", 1]},
            "solver": {"type": "object", "not": {"anyOf": [{"required": ["mobility"]}]}},
        },
        "definitions": {
            "time": {"oneOf": [
                {"type": "null"},
                {"type": "object", "required": ["list"], "additionalProperties": False,
                 "properties": {"list": {"type": "array", "minItems": 1, "items": {"type": "number"}}}},
            ]},
        },
    }

    def errors(self, document: dict) -> list:
        return sv.collect_errors(self.SCHEMA, document)

    def test_accepts_conforming_document(self) -> None:
        self.assertEqual([], self.errors({"id": "W101", "count": 2, "ratio": 0.5, "time": None,
                                          "flag": True, "kind": 1, "solver": {"max_iter": 3}}))
        self.assertEqual([], self.errors({"id": "E301", "time": {"list": [0.0, 1.0]}}))

    def test_rejects_each_violation(self) -> None:
        self.assertTrue(self.errors({}))                                   # required
        self.assertTrue(self.errors({"id": "X101"}))                        # pattern
        self.assertTrue(self.errors({"id": "W101", "extra": 1}))            # additionalProperties
        self.assertTrue(self.errors({"id": "W101", "count": 1.5}))          # integer
        self.assertTrue(self.errors({"id": "W101", "count": True}))         # bool is not integer
        self.assertTrue(self.errors({"id": "W101", "ratio": 0}))            # exclusiveMinimum
        self.assertTrue(self.errors({"id": "W101", "time": {"list": []}}))  # minItems via $ref/oneOf
        self.assertTrue(self.errors({"id": "W101", "time": {}}))            # oneOf: matches none
        self.assertTrue(self.errors({"id": "W101", "flag": 1}))             # const True != 1
        self.assertTrue(self.errors({"id": "W101", "kind": True}))          # enum 1 != True
        self.assertTrue(self.errors({"id": "W101", "solver": {"mobility": {}}}))  # not/anyOf

    def test_integer_valued_float_counts_as_integer(self) -> None:
        self.assertEqual([], self.errors({"id": "W101", "count": 2.0}))

    def test_unsupported_keyword_is_reported_not_ignored(self) -> None:
        with self.assertRaises(sv.UnsupportedSchemaKeyword):
            sv.collect_errors({"type": "object", "patternProperties": {}}, {})


class ArtifactManifestTest(unittest.TestCase):
    """reference_tcad/<case>/artifact_manifest.json files registered against the VM archive."""

    CATEGORY_BY_DIR = {"source": "A", "plt": "B", "mesh": "C", "tdr": "D", "log": "E"}

    @classmethod
    def setUpClass(cls) -> None:
        cls.schema = sv.load_schema(C7)
        cls.manifests = sorted((REPO / "reference_tcad").glob("*/artifact_manifest.json"))

    def test_first_batch_is_registered(self) -> None:
        cases = {path.parent.name for path in self.manifests}
        self.assertTrue({"pn2d_sentaurus2018", "schottky_charon_sentaurus2018",
                         "genius_bjt_sentaurus2022", "transportmodels_sentaurus2022"} <= cases, cases)

    def test_manifests_validate_and_are_consistent(self) -> None:
        for path in self.manifests:
            manifest = load_json(path)
            case = path.parent.name
            with self.subTest(case=case):
                self.assertEqual([], sv.collect_errors(self.schema, manifest))
                self.assertEqual(case, manifest["case"])
                paths = [a["path"] for a in manifest["artifacts"]]
                self.assertEqual(len(paths), len(set(paths)), "duplicate artifact paths")
                archived_sources = {a["path"][len("source/"):] for a in manifest["artifacts"]
                                    if a["path"].startswith("source/")}
                repo_source = REPO / "reference_tcad" / case / "source"
                for artifact in manifest["artifacts"]:
                    top = artifact["path"].split("/", 1)[0]
                    self.assertEqual(self.CATEGORY_BY_DIR[top], artifact["category"], artifact["path"])
                    producer = artifact["producer"]
                    if producer["tool"] == "repository":
                        name = artifact["path"].split("/")[-1]
                        repo_file = repo_source / name
                        self.assertTrue(repo_file.is_file(), f"{artifact['path']} claims a repository origin")
                        # C7 records archived bytes, whereas cmd_sha256 records LF
                        # text. Preserve exact archives first; allow checkout EOL
                        # expansion only when the archived byte count agrees.
                        content = repo_file.read_bytes()
                        if len(content) != artifact["bytes"]:
                            content = content.decode("utf-8-sig").replace("\r\n", "\n").encode("utf-8")
                        self.assertEqual(len(content), artifact["bytes"], artifact["path"])
                        self.assertEqual(hashlib.sha256(content).hexdigest(), artifact["sha256"], artifact["path"])
                    if producer["tool"] in {"sdevice", "sde", "snmesh", "sprocess"}:
                        self.assertIn("cmd", producer, artifact["path"])
                        variant = artifact["path"].split("/")[1] if artifact["path"].count("/") == 2 else None
                        expected = f"{variant}/{producer['cmd']}" if variant else producer["cmd"]
                        self.assertIn(expected, archived_sources, f"{artifact['path']}: producer cmd not archived")
                    if producer["tool"] == "sdevice":
                        self.assertIn("cmd_sha256", producer, artifact["path"])
                        if manifest["deck_matches_repository"]:
                            self.assertTrue((repo_source / producer["cmd"]).is_file(),
                                            f"{case}: {producer['cmd']} missing from the repository source")


class MessageCatalogueTest(unittest.TestCase):
    def test_catalogue_ids_are_unique_and_consistent(self) -> None:
        rows = []
        for line in WORKFLOW_DOC.read_text(encoding="utf-8").splitlines():
            match = CATALOGUE_ROW.match(line)
            if match:
                rows.append(match)
        ids = [m["id"] for m in rows]
        self.assertGreaterEqual(len(ids), 40)
        self.assertEqual(len(ids), len(set(ids)), "duplicate message ids")
        for match in rows:
            message_id, severity, exit_code = match["id"], match["severity"], match["exit"]
            if message_id.startswith("W"):
                self.assertEqual("warning", severity, message_id)
                self.assertEqual("\u2014", exit_code, message_id)
            else:
                self.assertEqual("fatal", severity, message_id)
                self.assertEqual(message_id[1], exit_code, message_id)

    def test_golden_message_ids_use_catalogue_templates(self) -> None:
        entries = catalogue()
        for path in FIXTURES.rglob("*.ir.json"):
            golden = load_json(path)
            for finding in golden["errors"] + golden["warnings"] + golden["preprocessor_findings"]:
                self.assertIn(finding["id"], entries, f"{path.name}: {finding['id']}")


if __name__ == "__main__":
    unittest.main()
