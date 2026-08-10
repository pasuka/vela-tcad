"""Regression tests for the fail-closed Sentaurus SDE/SDevice frontend.

Phase 0 and Phase 1 of the Sentaurus import plan replaced the historical regex
based SDE reader with an S-expression frontend that emits a versioned device IR,
and inverted the SDevice physics blacklist into a whitelist that emits a
versioned execution IR.  These tests pin both contracts:

* every command that the frozen fixtures use must parse into golden IR values;
* anything outside the whitelist must fail closed with file, line and command
  text; and
* both IR documents must validate against their JSON schema.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import (  # noqa: E402
    SdeParseError,
    legacy_summary,
    parse_sde_device_ir,
)
from sentaurus_execution_ir import (  # noqa: E402
    ExecutionIrError,
    build_execution_ir,
    classify_models,
    derive_analysis,
)
from sentaurus_import import parse_cmd, sentaurus_models  # noqa: E402
from sentaurus_sexp import SexpError, parse_file  # noqa: E402

PN2D_SDE = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sde.cmd"
PN2D_IV_CMD = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_iv_sdevice.cmd"
PN2D_BV_CMD = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_bv_sdevice.cmd"
DEVICE_IR_SCHEMA = REPO / "schemas" / "vela.sentaurus_device_ir.v1.schema.json"
EXECUTION_IR_SCHEMA = REPO / "schemas" / "vela.sentaurus_execution_ir.v1.schema.json"


def write_sde(root: Path, body: str) -> Path:
    path = root / "case_sde.cmd"
    path.write_text(body.strip() + "\n", encoding="utf-8")
    return path


def validate_against_schema(testcase: unittest.TestCase,
                            document: dict,
                            schema_path: Path) -> None:
    """Validate ``document`` against ``schema_path``.

    ``jsonschema`` is not a declared dependency of this repository, so the check
    degrades to a structural assertion when the module is unavailable.  The
    schema's ``required`` and ``additionalProperties`` rules are still enforced
    for the top-level object in that case.
    """
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    try:
        import jsonschema  # type: ignore[import-not-found]
    except ImportError:
        for key in schema.get("required", []):
            testcase.assertIn(key, document, f"missing required key {key!r}")
        if schema.get("additionalProperties") is False:
            unexpected = set(document) - set(schema.get("properties", {}))
            testcase.assertEqual(set(), unexpected,
                                 f"unexpected top-level keys {sorted(unexpected)}")
        return
    jsonschema.validate(document, schema)


class SdeSexpFrontendTest(unittest.TestCase):
    def test_tokenizer_reports_location_for_unbalanced_input(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_sexp_", dir=REPO / "build") as tmp:
            path = write_sde(Path(tmp), "(define L 2.0")
            with self.assertRaises(SexpError) as ctx:
                parse_file(path)
            self.assertIn("case_sde.cmd", str(ctx.exception))

    def test_tokenizer_skips_line_and_block_comments(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_sexp_", dir=REPO / "build") as tmp:
            path = write_sde(Path(tmp), """
; leading line comment
#| block
   comment |#
(define L 2.0) ; trailing comment
""")
            forms = parse_file(path)
            self.assertEqual(1, len(forms))
            self.assertEqual("define", forms[0].head())

    def test_nested_arithmetic_is_evaluated(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_sexp_", dir=REPO / "build") as tmp:
            path = write_sde(Path(tmp), """
(define H 0.5)
(define XJ 1.0)
(sdegeo:create-rectangle
  (position (- XJ 0.15) (/ H 2.0) 0.0)
  (position (+ XJ (* 0.15 2.0)) H 0.0)
  "Silicon" "R.Si")
""")
            ir = parse_sde_device_ir(path)
            region = ir["geometry"]["regions"][0]
            self.assertAlmostEqual(0.85, region["lower_left"][0])
            self.assertAlmostEqual(0.25, region["lower_left"][1])
            self.assertAlmostEqual(1.3, region["upper_right"][0])

    def test_undefined_symbol_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_sexp_", dir=REPO / "build") as tmp:
            path = write_sde(Path(tmp), """
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position Lmissing 0.5 0.0)
  "Silicon" "R.Si")
""")
            with self.assertRaises((SdeParseError, SexpError)) as ctx:
                parse_sde_device_ir(path)
            message = str(ctx.exception)
            self.assertIn("Lmissing", message)
            self.assertIn("case_sde.cmd:", message)


class SdeDeviceIrGoldenTest(unittest.TestCase):
    def setUp(self) -> None:
        self.ir = parse_sde_device_ir(PN2D_SDE)

    def test_schema_and_source(self) -> None:
        self.assertEqual("vela.sentaurus_device_ir.v1", self.ir["schema"])
        self.assertTrue(self.ir["source"].endswith("pn2d_sde.cmd"))
        validate_against_schema(self, self.ir, DEVICE_IR_SCHEMA)

    def test_defines_are_captured(self) -> None:
        self.assertEqual(
            {"L": 2.0, "H": 0.5, "XJ": 1.0},
            {key: self.ir["defines"][key] for key in ("L", "H", "XJ")},
        )

    def test_single_region_with_material(self) -> None:
        regions = self.ir["geometry"]["regions"]
        self.assertEqual(1, len(regions))
        region = regions[0]
        self.assertEqual("R.Si", region["name"])
        self.assertEqual("rectangle", region["shape"])
        self.assertEqual(0, region["order"])
        self.assertEqual([0.0, 0.0], region["lower_left"])
        self.assertEqual([2.0, 0.5], region["upper_right"])
        self.assertEqual(
            [{"region": "R.Si", "sde_material": "Silicon", "vela_material": "Si"}],
            self.ir["materials"])

    def test_two_contacts_carry_pick_points(self) -> None:
        contacts = {contact["name"]: contact for contact in self.ir["contacts"]}
        self.assertEqual({"Anode", "Cathode"}, set(contacts))
        anode = contacts["Anode"]["pick_points"]
        cathode = contacts["Cathode"]["pick_points"]
        self.assertEqual(1, len(anode))
        self.assertEqual(1, len(cathode))
        # The pick points were previously dropped entirely by the regex reader.
        self.assertAlmostEqual(0.0, anode[0]["x"])
        self.assertAlmostEqual(0.25, anode[0]["y"])
        self.assertAlmostEqual(2.0, cathode[0]["x"])
        self.assertAlmostEqual(0.25, cathode[0]["y"])

    def test_two_constant_profiles_and_placements(self) -> None:
        profiles = {entry["name"]: entry for entry in self.ir["doping"]["profiles"]}
        self.assertEqual({"P.Doping", "N.Doping"}, set(profiles))
        self.assertEqual("constant", profiles["P.Doping"]["type"])
        self.assertEqual("acceptor", profiles["P.Doping"]["carrier"])
        self.assertEqual("BoronActiveConcentration", profiles["P.Doping"]["species"])
        self.assertEqual(1e17, profiles["P.Doping"]["value_cm3"])
        self.assertEqual("donor", profiles["N.Doping"]["carrier"])
        self.assertEqual("PhosphorusActiveConcentration", profiles["N.Doping"]["species"])
        self.assertEqual(1e17, profiles["N.Doping"]["value_cm3"])

        placements = {entry["name"]: entry for entry in self.ir["doping"]["placements"]}
        self.assertEqual({"P.Place", "N.Place"}, set(placements))
        self.assertEqual("P.Doping", placements["P.Place"]["profile"])
        self.assertEqual("window", placements["P.Place"]["target_kind"])
        self.assertEqual("P.Window", placements["P.Place"]["target"])
        self.assertEqual("N.Window", placements["N.Place"]["target"])
        # Later placements win where windows overlap, so priority must be ordered.
        self.assertLess(placements["P.Place"]["priority"], placements["N.Place"]["priority"])

    def test_refeval_windows_and_refinement_sizes(self) -> None:
        windows = {entry["name"]: entry for entry in self.ir["doping"]["windows"]}
        self.assertEqual({"P.Window", "N.Window", "Junction.Window"}, set(windows))
        p_window = windows["P.Window"]
        self.assertEqual([0.0, 0.0], p_window["lower_left"])
        self.assertEqual([1.0, 0.5], p_window["upper_right"])
        self.assertEqual(windows, {entry["name"]: entry
                                   for entry in self.ir["mesh_control"]["windows"]})

        sizes = {entry["name"]: entry
                 for entry in self.ir["mesh_control"]["refinement_sizes"]}
        self.assertEqual({"Global.Mesh", "Junction.Mesh"}, set(sizes))
        self.assertEqual(0.05, sizes["Global.Mesh"]["max_x"])
        self.assertEqual(0.05, sizes["Global.Mesh"]["max_y"])
        self.assertEqual(0.01, sizes["Junction.Mesh"]["max_x"])
        self.assertEqual(0.02, sizes["Junction.Mesh"]["max_y"])

        placements = {entry["name"]: entry
                      for entry in self.ir["mesh_control"]["refinement_placements"]}
        self.assertEqual({"Global.Mesh.Place", "Junction.Mesh.Place"}, set(placements))
        self.assertEqual("Global.Mesh", placements["Global.Mesh.Place"]["refinement"])
        self.assertEqual("region", placements["Global.Mesh.Place"]["target_kind"])
        self.assertEqual("R.Si", placements["Global.Mesh.Place"]["target"])
        self.assertEqual("window", placements["Junction.Mesh.Place"]["target_kind"])
        self.assertEqual("Junction.Window", placements["Junction.Mesh.Place"]["target"])

    def test_build_mesh_directive_is_recorded(self) -> None:
        build = self.ir["mesh_control"]["build"]
        self.assertEqual("snmesh", build["engine"])
        self.assertEqual("pn2d", build["prefix"])

    def test_legacy_summary_projection_is_stable(self) -> None:
        summary = legacy_summary(self.ir)
        self.assertEqual({"Anode", "Cathode"},
                         {entry["name"] for entry in summary["contacts"]})
        self.assertEqual(2, len(summary["doping_profiles"]))
        self.assertEqual(2, len(summary["doping_placements"]))
        self.assertEqual(1, len(summary["geometry"]["rectangles"]))
        self.assertIn("device_ir", summary)


class SdeFailClosedTest(unittest.TestCase):
    HEADER = """
(define L 2.0)
(define H 0.5)
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L H 0.0)
  "Silicon" "R.Si")
"""

    def parse_body(self, body: str) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_fail_", dir=REPO / "build") as tmp:
            parse_sde_device_ir(write_sde(Path(tmp), self.HEADER + body))

    FAIL_CLOSED = (SdeParseError, SexpError)

    def test_unknown_command_fails_closed_with_location(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            self.parse_body('(sdegeo:invent-a-command "R.Si" 1.0)\n')
        message = str(ctx.exception)
        self.assertIn("case_sde.cmd:", message)
        self.assertIn("sdegeo:invent-a-command", message)

    def test_unknown_top_level_namespace_fails_closed(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            self.parse_body('(sdeplugin:do-something "x")\n')
        self.assertIn("sdeplugin:do-something", str(ctx.exception))

    def test_known_but_unimplemented_command_reports_a_reason(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            self.parse_body(
                '(sdedr:define-gaussian-profile "G" "BoronActiveConcentration"\n'
                '  "PeakPos" 0.0 "PeakVal" 1e18 "ValueAtDepth" 1e15 "Depth" 0.2\n'
                '  "Gauss" "Factor" 0.8)\n')
        message = str(ctx.exception)
        self.assertIn("sdedr:define-gaussian-profile", message)
        self.assertIn("case_sde.cmd:", message)

    def test_unknown_dopant_species_fails_closed(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            self.parse_body('(sdedr:define-constant-profile "X" "UnobtainiumConcentration" 1e17)\n')
        self.assertIn("Unobtainium", str(ctx.exception))

    def test_unknown_material_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sde_fail_", dir=REPO / "build") as tmp:
            path = write_sde(Path(tmp), """
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position 1.0 1.0 0.0)
  "Unobtainium" "R.X")
""")
            with self.assertRaises(SdeParseError) as ctx:
                parse_sde_device_ir(path)
            self.assertIn("Unobtainium", str(ctx.exception))

    def test_placement_referencing_unknown_profile_fails_closed(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            self.parse_body('(sdedr:define-constant-profile-region "P" "Missing" "R.Si")\n')
        self.assertIn("Missing", str(ctx.exception))

    def test_duplicate_region_name_fails_closed(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            self.parse_body("""
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L H 0.0)
  "Silicon" "R.Si")
""")
        self.assertIn("R.Si", str(ctx.exception))

    def test_every_frozen_fixture_either_parses_or_fails_closed(self) -> None:
        fixtures = sorted((REPO / "reference_tcad").glob("*/source/*sde.cmd"))
        self.assertGreaterEqual(len(fixtures), 3)
        for fixture in fixtures:
            with self.subTest(fixture=fixture.name):
                try:
                    ir = parse_sde_device_ir(fixture)
                except (SdeParseError, SexpError) as error:
                    # Fail-closed is an acceptable outcome, silence is not.
                    self.assertIn(fixture.name, str(error))
                    continue
                self.assertEqual("vela.sentaurus_device_ir.v1", ir["schema"])
                validate_against_schema(self, ir, DEVICE_IR_SCHEMA)


class ExecutionIrTest(unittest.TestCase):
    def build(self, path: Path, **kwargs):
        summary = parse_cmd(path, {})
        return build_execution_ir(summary, str(path), sentaurus_models(summary), **kwargs)

    def test_iv_fixture_produces_valid_execution_ir(self) -> None:
        ir = self.build(PN2D_IV_CMD)
        self.assertEqual("vela.sentaurus_execution_ir.v1", ir["schema"])
        validate_against_schema(self, ir, EXECUTION_IR_SCHEMA)
        self.assertEqual({"Anode", "Cathode"},
                         {entry["name"] for entry in ir["electrodes"]})
        self.assertEqual("iv", ir["analysis"]["kind"])
        self.assertEqual("Anode", ir["analysis"]["final_contact"])
        self.assertEqual([], ir["unsupported"])

    def test_bv_fixture_is_classified_from_the_execution_ir(self) -> None:
        ir = self.build(PN2D_BV_CMD)
        validate_against_schema(self, ir, EXECUTION_IR_SCHEMA)
        self.assertEqual("bv", ir["analysis"]["kind"])
        self.assertIn("avalanche", ir["analysis"]["reason"].lower())

    def test_stage_order_is_preserved_with_dependency_edges(self) -> None:
        ir = self.build(PN2D_IV_CMD)
        stages = ir["stages"]
        self.assertGreaterEqual(len(stages), 1)
        self.assertEqual(list(range(len(stages))), [stage["index"] for stage in stages])
        self.assertEqual([], stages[0]["depends_on"])
        for previous, stage in zip(stages, stages[1:]):
            self.assertEqual([previous["index"]], stage["depends_on"])

    def test_analysis_kind_is_derived_not_supplied(self) -> None:
        self.assertEqual(
            "equilibrium",
            derive_analysis({"sweeps": []}, classify_models(set()))["kind"])

        self.assertEqual(
            "iv",
            derive_analysis(
                {"sweeps": [{"contact": "Anode", "stop": 1.0}]},
                classify_models({"Mobility", "SRH"}))["kind"])

        bv = derive_analysis(
            {"sweeps": [{"contact": "Cathode", "stop": -30.0}]},
            classify_models({"Avalanche"}))
        self.assertEqual("bv", bv["kind"])
        self.assertEqual(-30.0, bv["final_voltage"])

        self.assertEqual(
            "cv",
            derive_analysis(
                {"sweeps": [{"contact": "Anode", "stop": 1.0, "equations": ["ACCoupled"]}]},
                classify_models({"Mobility"}))["kind"])

    def test_ambiguous_analysis_fails_closed(self) -> None:
        with self.assertRaises(ExecutionIrError):
            derive_analysis(
                {"sweeps": [{"contact": "Cathode", "stop": -30.0,
                             "equations": ["ACCoupled"]}]},
                classify_models({"Avalanche"}))

    def test_unexpanded_placeholder_goal_fails_closed(self) -> None:
        with self.assertRaises(ExecutionIrError) as ctx:
            derive_analysis(
                {"sweeps": [{"contact": "Gate", "stop": "@Vg@"}]},
                classify_models({"Mobility"}))
        self.assertIn("--template-var", str(ctx.exception))

    def test_physics_whitelist_is_inverted(self) -> None:
        classified = classify_models({"Mobility", "SRH", "Trap", "Hydrodynamic",
                                      "TotallyMadeUpModel"})
        supported = {entry["model"] for entry in classified["supported"]}
        self.assertIn("Mobility", supported)
        self.assertIn("SRH", supported)
        self.assertIn("Trap", {entry["model"] for entry in classified["metadata_only"]})
        unsupported = {entry["model"] for entry in classified["unsupported"]}
        self.assertIn("Hydrodynamic", unsupported)
        # An unknown token must not be silently accepted.
        self.assertIn("TotallyMadeUpModel", unsupported)

    def test_unsupported_physics_fails_closed_unless_explicitly_allowed(self) -> None:
        summary = parse_cmd(PN2D_IV_CMD, {})
        models = sentaurus_models(summary) | {"Hydrodynamic"}
        with self.assertRaises(ExecutionIrError):
            build_execution_ir(summary, str(PN2D_IV_CMD), models)
        report = build_execution_ir(summary, str(PN2D_IV_CMD), models,
                                    allow_unsupported=True)
        self.assertIn("Hydrodynamic",
                      {entry["model"] for entry in report["unsupported"]})


if __name__ == "__main__":
    unittest.main()
