"""WP1 regression tests: SDevice cmd parser -> execution IR v2.

* every fixture in tests/fixtures/sdevice/fixtures.json reproduces its golden
  (full JSON equality or findings-only comparison, see fixtures.json);
* one injected deck per fail-closed class of plan section 4.7 with the expected
  message id and line;
* property tests (comment stripping, keyword case, BOM/CRLF invariance);
* every produced IR validates against the C1 schema; the CLI exit code equals
  ``verdict.exit_code``.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from vela_sdevice import schema_validator as sv  # noqa: E402
from vela_sdevice.cmd_parser import catalogue, parse_cmd_file, parse_cmd_text  # noqa: E402
from vela_sdevice.cmd_parser.lexer import strip_comments  # noqa: E402

FIXTURES = REPO / "tests" / "fixtures" / "sdevice"
SCHEMA = sv.load_schema(REPO / "schemas" / "vela.sentaurus_execution_ir.v2.schema.json")
WORKFLOW_DOC = REPO / "docs" / "sentaurus_sdevice_workflow.md"
CATALOGUE_ROW = re.compile(r"^\| (?P<id>[WE][1-8]\d{2}) \| (?P<severity>warning|fatal) \| (?P<exit>[0-8]|\u2014) \| `(?P<template>[^`]*)` \|")

BASE_DECK = """File {{
  Grid = "x_msh.tdr"
  Parameter = "models.par"
}}
Electrode {{
  {{ Name="a" Voltage=0.0 }}
  {{ Name="c" Voltage=0.0 }}
}}
Physics {{
  {physics}
}}
Math {{ Iterations=20 }}
Solve {{
  Coupled {{ Poisson }}
  {solve}
}}
"""


def deck(physics: str = "Recombination(SRH)", solve: str = "", extra: str = "") -> str:
    return BASE_DECK.format(physics=physics, solve=solve) + extra


def json_diff(actual, expected, path="$", out=None):
    out = [] if out is None else out
    if isinstance(actual, dict) and isinstance(expected, dict):
        for key in sorted(set(actual) | set(expected)):
            if key not in actual:
                out.append(f"{path}.{key}: missing (golden={json.dumps(expected[key])[:100]})")
            elif key not in expected:
                out.append(f"{path}.{key}: unexpected (parser={json.dumps(actual[key])[:100]})")
            else:
                json_diff(actual[key], expected[key], f"{path}.{key}", out)
    elif isinstance(actual, list) and isinstance(expected, list):
        if len(actual) != len(expected):
            out.append(f"{path}: length parser={len(actual)} golden={len(expected)}")
        for index, (a, e) in enumerate(zip(actual, expected)):
            json_diff(a, e, f"{path}[{index}]", out)
    elif actual != expected or type(actual) is bool and type(expected) is not bool:
        out.append(f"{path}: parser={json.dumps(actual)[:120]} golden={json.dumps(expected)[:120]}")
    return out


class GoldenFixtureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.index = json.loads((FIXTURES / "fixtures.json").read_text(encoding="utf-8"))

    def test_every_fixture_matches_its_golden(self) -> None:
        for entry in self.index["fixtures"]:
            with self.subTest(deck=entry["deck"]):
                ir = parse_cmd_file(FIXTURES / entry["deck"])
                golden = json.loads((FIXTURES / entry["golden"]).read_text(encoding="utf-8"))
                self.assertEqual([], sv.collect_errors(SCHEMA, ir), "IR does not validate against C1")
                self.assertEqual(entry["expected_exit_code"], ir["verdict"]["exit_code"])
                if entry["compare"] == "full":
                    diff = json_diff(ir, golden)
                else:
                    diff = json_diff(ir["verdict"], golden["verdict"], "$.verdict")
                    for key in ("preprocessor_findings", "errors"):
                        if key in golden:
                            mine = [{k: f[k] for k in ("id", "line", "text")} for f in ir[key]]
                            diff += json_diff(mine, golden[key], f"$.{key}")
                self.assertEqual([], diff, "\n".join(diff))

    def test_findings_use_catalogue_templates(self) -> None:
        rows = {}
        for line in WORKFLOW_DOC.read_text(encoding="utf-8").splitlines():
            match = CATALOGUE_ROW.match(line)
            if match:
                rows[match["id"]] = match
        for message_id, (severity, template) in catalogue.TEMPLATES.items():
            self.assertIn(message_id, rows, f"{message_id} missing from docs/sentaurus_sdevice_workflow.md")
            self.assertEqual(rows[message_id]["severity"], severity, message_id)
            self.assertEqual(rows[message_id]["template"], template, message_id)


class FailClosedTest(unittest.TestCase):
    def findings(self, text: str, path_name: str = "case_des.cmd") -> dict:
        return parse_cmd_text(text, path_name=path_name)

    def assert_fatal(self, ir: dict, message_id: str, line: int, text_fragment: str = "") -> None:
        pool = ir["preprocessor_findings"] + ir["errors"]
        matches = [f for f in pool if f["id"] == message_id]
        self.assertTrue(matches, f"{message_id} not reported; got {[f['id'] for f in pool]}")
        self.assertIn(line, [f["line"] for f in matches], f"{message_id} lines {[f['line'] for f in matches]}")
        if text_fragment:
            self.assertTrue(any(text_fragment in f["text"] for f in matches), [f["text"] for f in matches])
        self.assertFalse(ir["verdict"]["accepted"])
        self.assertIn(message_id, ir["verdict"]["fatal_ids"])
        self.assertEqual([], sv.collect_errors(SCHEMA, ir))

    def test_preprocessor_directive_E101(self) -> None:
        ir = self.findings("#if @Vg@ > 0\n" + deck())
        self.assert_fatal(ir, "E101", 1, "#if")
        self.assert_fatal(ir, "E102", 1, "@Vg@")
        self.assertEqual(1, ir["verdict"]["exit_code"])
        self.assertEqual("rejected", ir["analysis"]["kind"])

    def test_tcl_block_E103(self) -> None:
        ir = self.findings(deck(extra='!( set x 1 )!\n'))
        self.assert_fatal(ir, "E103", 17)

    def test_syntax_error_E104(self) -> None:
        ir = self.findings('File { Grid = "x.tdr" Parameter="m.par" }\nElectrode { { Name="a" Voltage=0 }\nSolve { Coupled { Poisson } }\n')
        self.assertTrue(any(f["id"] == "E104" for f in ir["errors"]), ir["errors"])
        self.assertEqual(1, ir["verdict"]["exit_code"])

    def test_unsupported_section_E105(self) -> None:
        ir = self.findings(deck(extra='Thermode {\n  { Name="a" Temperature=300 }\n}\nSystem {\n  Vsource_pset v1 (a 0) { dc = 0 }\n}\n'))
        self.assert_fatal(ir, "E105", 17, "Thermode")
        self.assert_fatal(ir, "E105", 20, "System")

    def test_unsupported_physics_E106(self) -> None:
        ir = self.findings(deck(physics="Mobility(DopingDep IALMob)\n  Thermodynamic\n  Traps((Donor Level EnergyMid=0 Conc=1e15))"))
        self.assert_fatal(ir, "E106", 10, "IALMob")
        self.assert_fatal(ir, "E106", 11, "Thermodynamic")
        self.assert_fatal(ir, "E106", 12)
        keys = {m["key"]: m for m in ir["physics"]["effective"]}
        self.assertEqual("unsupported", keys["Mobility.IALMob"]["status"])
        self.assertEqual("E106", keys["Mobility.IALMob"]["finding_id"])

    def test_unknown_physics_token_E107(self) -> None:
        ir = self.findings(deck(physics="Mobility(DopingDep ENORMAL_MODEL)\n  FooBar"))
        self.assert_fatal(ir, "E107", 10, "ENORMAL_MODEL")
        self.assert_fatal(ir, "E107", 11, "FooBar")

    def test_unsupported_electrode_attribute_E108(self) -> None:
        text = deck().replace('{ Name="a" Voltage=0.0 }', '{ Name="a" Voltage=0.0 Charge=1e-12 Poisson=Neumann AreaFactor=2 }')
        ir = self.findings(text)
        self.assert_fatal(ir, "E108", 6, "Charge")
        self.assert_fatal(ir, "E108", 6, "Poisson=Neumann")
        self.assert_fatal(ir, "E108", 6, "AreaFactor=2")

    def test_unsupported_solve_statement_E109(self) -> None:
        ir = self.findings(deck(solve='System("rm x.plt")\n  Transient(InitialTime=0 FinalTime=1) { Coupled { Poisson Electron Hole } }\n  Unset("a")'))
        self.assert_fatal(ir, "E109", 15, 'System("rm x.plt")')
        self.assert_fatal(ir, "E109", 16, "Transient")
        self.assert_fatal(ir, "E109", 17, "Unset")
        self.assertEqual(["unsupported"] * 3, [s["kind"] for s in ir["solve"]["statements"][1:]])

    def test_goal_forms_E110(self) -> None:
        two_goals = 'Quasistationary(InitialStep=0.1 Goal{Name="a" Voltage=1} Goal{Name="c" Voltage=1}) { Coupled { Poisson Electron Hole } }'
        ir = self.findings(deck(solve=two_goals))
        self.assert_fatal(ir, "E110", 15, 'Goal{Name="c"')
        ir = self.findings(deck(solve='Quasistationary(Goal{Model="DeviceTemperature" Parameter="Temperature" Value=400}) { Coupled { Poisson } }'))
        self.assert_fatal(ir, "E110", 15, "Model=")

    def test_unknown_electrode_attribute_E111(self) -> None:
        ir = self.findings(deck().replace('{ Name="c" Voltage=0.0 }', '{ Name="c" Voltage=0.0 Wibble=3 }'))
        self.assert_fatal(ir, "E111", 7, "Wibble")

    def test_file_keyword_E112_and_W107(self) -> None:
        text = deck().replace('Parameter = "models.par"', 'Parameter = "models.par"\n  MobilityDoping = "m.tdr"\n  NewtonPlot = "np"')
        ir = self.findings(text)
        self.assert_fatal(ir, "E112", 4, "MobilityDoping")
        self.assertEqual([("W107", 5)], [(w["id"], w["line"]) for w in ir["warnings"] if w["id"] == "W107"])
        self.assertEqual([{"keyword": "NewtonPlot", "value": "np", "line": 5, "finding_id": "W107"}], ir["files"]["ignored"])

    def test_inconsistent_scoped_physics_E113(self) -> None:
        text = deck() + 'Physics(Material="Silicon") { Mobility(DopingDep) }\nPhysics(Material="Germanium") { Mobility(DopingDep Enormal) }\n'
        ir = self.findings(text)
        self.assert_fatal(ir, "E113", 18, 'Physics(Material="Germanium")')
        self.assertEqual(2, len(ir["physics"]["scoped"]))

    def test_grid_extension_E114(self) -> None:
        ir = self.findings(deck().replace('Grid = "x_msh.tdr"', 'Grid = "x.grd2"'))
        self.assert_fatal(ir, "E114", 2, "x.grd2")
        self.assertEqual("unknown", ir["files"]["grid"]["kind"])

    def test_undecodable_E115(self) -> None:
        with tempfile.TemporaryDirectory(dir=REPO / "build" if (REPO / "build").is_dir() else None) as tmp:
            path = Path(tmp) / "bad_des.cmd"
            path.write_bytes(b"\xff\xfe\x00\x81 File {")
            ir = parse_cmd_file(path)
        self.assertEqual(["E115"], ir["verdict"]["fatal_ids"])
        self.assertEqual([], sv.collect_errors(SCHEMA, ir))

    def test_unsupported_equation_E116(self) -> None:
        ir = self.findings(deck(solve="Coupled { Poisson Electron Hole Temperature }"))
        self.assert_fatal(ir, "E116", 15, "Temperature")

    def test_regexp_electrode_E117(self) -> None:
        ir = self.findings(deck().replace('{ Name="c" Voltage=0.0 }', '{ Name=Regexp("c.*") Voltage=0.0 }'))
        self.assert_fatal(ir, "E117", 7, "Regexp")

    def test_missing_parameter_E301_exit_3(self) -> None:
        ir = self.findings(deck().replace('  Parameter = "models.par"\n', ""))
        self.assert_fatal(ir, "E301", 1, "File {")
        self.assertEqual(3, ir["verdict"]["exit_code"])
        self.assertIsNone(ir["files"]["parameter"])

    def test_exit_code_prefers_input_fatal_over_missing_declaration(self) -> None:
        ir = self.findings(deck(physics="Thermodynamic").replace('  Parameter = "models.par"\n', ""))
        self.assertEqual(["E301", "E106"], ir["verdict"]["fatal_ids"])  # first-occurrence order
        self.assertEqual(1, ir["verdict"]["exit_code"])


class NamingAndMappingTest(unittest.TestCase):
    def test_file_name_completion_rules(self) -> None:
        ir = parse_cmd_text(deck().replace('Grid = "x_msh.tdr"', 'Grid = "x_msh"\n  Current = "iv"\n  Plot = "dev"\n  Output = "run.log"'),
                            path_name="case_des.cmd")
        files = ir["files"]
        self.assertEqual(("x_msh.tdr", "tdr"), (files["grid"]["resolved"], files["grid"]["kind"]))
        self.assertEqual("iv_des.plt", files["current"]["resolved"])
        self.assertEqual(("dev_des.tdr", "dev_des.vtu"), (files["plot"]["resolved"], files["plot"]["vela"]))
        self.assertEqual("run.log_des.log", files["output"]["resolved"])
        defaults = parse_cmd_text(deck(), path_name="case_des.cmd")["files"]
        self.assertEqual("case_des.plt", defaults["current"]["resolved"])
        self.assertEqual("case_des.tdr", defaults["plot"]["resolved"])
        self.assertEqual("output_des.log", defaults["output"]["resolved"])
        self.assertEqual("x.grd", parse_cmd_text(deck().replace("x_msh.tdr", "x.grd"), path_name="c.cmd")["files"]["grid"]["resolved"])
        self.assertEqual("dfise", parse_cmd_text(deck().replace("x_msh.tdr", "x.grd"), path_name="c.cmd")["files"]["grid"]["kind"])

    def test_set_mode_current_and_current_goal(self) -> None:
        solve = 'Quasistationary(InitialStep=0.1 Goal{Name="a" Voltage=1}) { Coupled { Poisson Electron Hole } }\n' \
                '  Set("a" mode current)\n' \
                '  Quasistationary(InitialStep=0.1 Goal{Name="a" Current=1e-4}) { Coupled { Poisson Electron Hole } }'
        ir = parse_cmd_text(deck(solve=solve), path_name="c_des.cmd")
        self.assertTrue(ir["verdict"]["accepted"], ir["errors"])
        stages = ir["stages"]
        self.assertEqual(["initial", "sweep", "sweep"], [s["kind"] for s in stages])
        self.assertNotIn("current_mode_contacts", stages[1])
        self.assertEqual(["a"], stages[2]["current_mode_contacts"])
        self.assertEqual({"contact": "a", "current_A_per_um": 0.0001}, stages[2]["goal"])
        self.assertEqual("multi_stage", ir["analysis"]["kind"])

    def test_scoped_physics_union_and_fixed_charge(self) -> None:
        text = deck(physics="Recombination(SRH)") + \
            'Physics(Material="Silicon") { Mobility(DopingDep) }\n' \
            'Physics(MaterialInterface="Silicon/Oxide") { Traps(Conc=1e11 FixedCharge) }\n' \
            'Physics(Material="Oxide") { Traps((FixedCharge SpatialShape=Uniform Conc=1e16)) }\n'
        ir = parse_cmd_text(text, path_name="c_des.cmd")
        self.assertTrue(ir["verdict"]["accepted"], ir["errors"])
        self.assertEqual(["Mobility.DopingDependence", "Recombination.SRH"], [m["key"] for m in ir["physics"]["effective"]])
        self.assertEqual([{"scope": "MaterialInterface", "target": "Silicon/Oxide", "line": 18, "conc_cm2": 1e11},
                          {"scope": "Material", "target": "Oxide", "line": 19, "conc_cm3": 1e16}],
                         ir["physics"]["fixed_charges"])


class PropertyTest(unittest.TestCase):
    def test_comment_stripping_protects_quoted_text(self) -> None:
        text = 'File { Grid = "a*b#c.tdr" } * trailing\n# whole line\nPlot { eCurrent/Vector }\n'
        stripped = strip_comments(text)
        self.assertIn('"a*b#c.tdr"', stripped)
        self.assertNotIn("trailing", stripped)
        self.assertNotIn("whole line", stripped)
        self.assertIn("eCurrent/Vector", stripped)
        self.assertEqual(len(text), len(stripped))

    def test_keyword_case_variants_give_same_effective_set(self) -> None:
        a = parse_cmd_text(deck(physics="Mobility(DopingDep HighFieldsaturation(GradQuasiFermi))\n  Recombination(SRH(DopingDep))"), path_name="c_des.cmd")
        b = parse_cmd_text(deck(physics="Mobility(DopingDependence HighFieldSaturation(GradQuasiFermi))\n  Recombination(SRH(DopingDependence))"), path_name="c_des.cmd")
        self.assertEqual(a["physics"]["effective"], b["physics"]["effective"])
        lower = parse_cmd_text(deck().replace("Iterations=20", "iterations=20 notdamped=5"), path_name="c_des.cmd")
        self.assertEqual(20, lower["math"]["iterations"])
        self.assertEqual("notdamped", lower["math"]["ignored"][0]["keyword"])

    def test_bom_and_crlf_variants_are_invariant(self) -> None:
        path = FIXTURES / "pn2d_sentaurus2018" / "pn2d_iv_sdevice.cmd"
        text = path.read_bytes().decode("utf-8-sig").replace("\r\n", "\n")
        reference = json.dumps(parse_cmd_text(text, path_name=path.name), sort_keys=True)
        crlf = json.dumps(parse_cmd_text(text.replace("\n", "\r\n"), path_name=path.name), sort_keys=True)
        bom = json.dumps(parse_cmd_text("\ufeff" + text, path_name=path.name), sort_keys=True)
        self.assertEqual(reference, crlf)
        self.assertEqual(reference, bom)

    def test_section_order_does_not_matter(self) -> None:
        text = deck()
        blocks = ['File {\n  Grid = "x_msh.tdr"\n  Parameter = "models.par"\n}',
                  'Electrode {\n  { Name="a" Voltage=0.0 }\n  { Name="c" Voltage=0.0 }\n}',
                  'Physics {\n  Recombination(SRH)\n}', 'Math { Iterations=20 }',
                  'Solve {\n  Coupled { Poisson }\n  \n}']
        self.assertEqual(text, "\n".join(blocks) + "\n")
        reordered = "\n".join(reversed(blocks)) + "\n"
        a = parse_cmd_text(text, path_name="c_des.cmd")
        b = parse_cmd_text(reordered, path_name="c_des.cmd")
        for key in ("files", "electrodes", "physics", "math", "analysis"):
            strip = json.dumps(a[key], sort_keys=True)
            strip_b = json.dumps(b[key], sort_keys=True)
            self.assertEqual(re.sub(r'"line": \d+', "", strip), re.sub(r'"line": \d+', "", strip_b), key)


class CliTest(unittest.TestCase):
    def test_module_exit_code_equals_verdict(self) -> None:
        cases = [("pn2d_sentaurus2018/pn2d_0v_sdevice.cmd", 3), ("bvmethods_sentaurus2018/aba_poisson_sdevice.cmd", 0),
                 ("singledevice_sentaurus2018/singledevice_sdevice.cmd", 1)]
        env = {"PYTHONPATH": str(REPO / "scripts"), "PYTHONIOENCODING": "utf-8"}
        import os
        env = {**os.environ, **env}
        for rel, expected in cases:
            with self.subTest(deck=rel):
                proc = subprocess.run([sys.executable, "-m", "vela_sdevice.cmd_parser", str(FIXTURES / rel)],
                                      capture_output=True, text=True, encoding="utf-8", env=env, cwd=str(REPO))
                self.assertEqual(expected, proc.returncode, proc.stderr)
                ir = json.loads(proc.stdout)
                self.assertEqual(expected, ir["verdict"]["exit_code"])

    def test_missing_deck_exits_3(self) -> None:
        import os
        env = {**os.environ, "PYTHONPATH": str(REPO / "scripts")}
        proc = subprocess.run([sys.executable, "-m", "vela_sdevice.cmd_parser", str(FIXTURES / "nope_des.cmd")],
                              capture_output=True, text=True, env=env, cwd=str(REPO))
        self.assertEqual(3, proc.returncode)
        self.assertIn("E305", proc.stderr)


if __name__ == "__main__":
    unittest.main()
