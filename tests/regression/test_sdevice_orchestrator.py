"""WP5 gate, lifecycle, publication and hermetic archive acceptance."""
from __future__ import annotations
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import tempfile
import time
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from build_vela_sdevice_pyz import build_archive
from vela_sdevice.log_levels import WORKBENCH_COMMAND_ENV
from vela_sdevice.materials import default_models_par_path
from vela_sdevice.plt import read_plt_document
from vela_sdevice.orchestrator.outputs import vtu_metadata
import run_pn2d_same_mesh_sdevice_vector_export as dfise_writers

FAKE = REPO / "tests/fixtures/sdevice_fake_runner/fake_runner.py"


def inputs(root, *, multistage=False, fixture=False):
    nodes = {0: (0., 0.), 1: (1., 0.), 2: (1., 1.), 3: (0., 1.)}
    triangles, contacts = [(0, 1, 2), (0, 2, 3)], {"anode": [0, 3], "cathode": [1, 2]}
    grid = dfise_writers.render_dfise_grid(nodes, triangles, contacts)
    grid = grid.replace('"Anode"', '"anode"').replace('"Cathode"', '"cathode"')
    (root / "mesh.grd").write_text(grid)
    doping = root / "doping_source.csv"
    doping.write_text("node_id,donors_cm3,acceptors_cm3\n" + "".join(f"{i},1e15,0\n" for i in range(4)))
    (root / "mesh.dat").write_text(dfise_writers.render_dfise_doping(nodes, triangles, contacts, doping))
    (root / "models.par").write_bytes(default_models_par_path().read_bytes())
    cmd = root / "device_des.cmd"
    if fixture:
        text = (REPO / "tests/fixtures/sdevice/pn2d_sentaurus2018/pn2d_iv_sdevice.cmd").read_text()
        import re
        text = re.sub(r'Grid\s*=\s*"[^"]+"', 'Grid="mesh.grd" Parameter="models.par"', text)
        text = text.replace('"Anode"', '"anode"').replace('"Cathode"', '"cathode"')
        cmd.write_text(text)
    else:
        second = 'Quasistationary(Goal {Name="anode" Voltage=2}) {Coupled {Poisson Electron Hole}}' if multistage else ""
        cmd.write_text('''File {Grid="mesh.grd" Parameter="models.par" Output="run" Current="current" Plot="final"}
            Electrode { {Name="anode" Voltage=0} {Name="cathode" Voltage=0} }
            Solve { Coupled {Poisson} Coupled {Poisson Electron Hole}
              Quasistationary(Goal {Name="anode" Voltage=1}) {Coupled {Poisson Electron Hole}}
        ''' + second + "}")
    return cmd


class OrchestratorTests(unittest.TestCase):
    def test_workbench_main_log_keeps_progress_at_info_and_streams_it_at_error(self):
        inputs(self.root, multistage=True)
        invocation = 'sd-sim.exe device_des.cmd'
        log_path = self.root / "run_des.log"
        for level in ("info", "error"):
            with self.subTest(level=level):
                log_path.unlink(missing_ok=True)
                done = self.call(str(self.cmd), "--log-level", level, env={WORKBENCH_COMMAND_ENV: invocation})
                self.assert_exit(0, done)
                points = [line for line in done.stdout.splitlines() if line.endswith(" V accepted")]
                self.assertGreater(len(points), 1)
                self.assertEqual(len(points), len(set(points)))
                self.assertGreater(len({line.split()[0] for line in points}), 1)
                lines = log_path.read_text(encoding="utf-8").splitlines()
                self.assertEqual(lines.count(f"[INVOCATION] {invocation}"), 1)
                logged_points = [line.removeprefix("[PROGRESS] ") for line in lines if line.startswith("[PROGRESS] ")]
                self.assertEqual(logged_points, points if level == "info" else [])
                self.assertTrue(any(line.startswith("[SUMMARY] ") for line in lines))

    def test_workbench_launcher_error_writes_main_log(self):
        from vela_workbench_cases import Case, run_case

        log_name = "logs/test_des.log"
        case = Case(name="test", source=self.root, arguments=[self.cmd.name],
                    hashes={}, work=self.work.name, log=log_name)
        with self.assertRaises(FileNotFoundError):
            run_case(self.root / "missing-sd-sim", case, self.root, log_level="error")
        lines = (self.root / log_name).read_text(encoding="utf-8").splitlines()
        self.assertTrue(lines[0].startswith("[INVOCATION] "))
        self.assertTrue(lines[1].startswith("[FATAL] "))
        self.assertEqual(list(self.root.rglob("*.log")), [self.root / log_name])

    def test_workbench_frontend_error_preserves_console_diagnostic(self):
        from vela_workbench_cases import Case, run_case

        case = Case(name="test", source=self.root, hashes={}, work=self.work.name, log="run_des.log",
                    arguments=["-c", "import sys; sys.stderr.write('invalid input\\n'); sys.exit(4)"])
        self.assertEqual(run_case(sys.executable, case, self.root), ("Solver exited (code {})", 4))
        lines = (self.root / case.log).read_text(encoding="utf-8").splitlines()
        self.assertTrue(lines[0].startswith("[INVOCATION] "))
        self.assertEqual(lines[1:], ["[FATAL] sd-sim exited with code 4", "[CONSOLE] invalid input"])

    def test_current_stop_completes_stage_before_voltage_goal(self):
        text = self.cmd.read_text().replace('Solve {',
            'Math {BreakCriteria {Current(Contact="anode" AbsVal=4e-8)}} Solve {')
        self.cmd.write_text(text)
        done = self.call(str(self.cmd), '--exit-on-failure', mode='partial')
        self.assert_exit(0, done)
        runtime = json.loads((self.work/'run_manifest.json').read_text())['stages'][1]
        self.assertFalse(runtime['goal_reached'])
        self.assertEqual(runtime['stop_criterion_reached'], 'current')
        self.assertAlmostEqual(runtime['reached_voltages_V']['anode'], .5)

    def test_unmet_current_stop_does_not_hide_incomplete_stage(self):
        text = self.cmd.read_text().replace('Solve {',
            'Math {BreakCriteria {Current(Contact="anode" AbsVal=9e-8)}} Solve {')
        self.cmd.write_text(text)
        done = self.call(str(self.cmd), '--exit-on-failure', mode='partial')
        self.assert_exit(5, done)
        runtime = json.loads((self.work/'run_manifest.json').read_text())['stages'][1]
        self.assertFalse(runtime['goal_reached'])
        self.assertNotIn('stop_criterion_reached', runtime)

    def test_slow_background_publisher_finishes_before_final_publication(self):
        from vela_sdevice.__main__ import main
        from vela_sdevice.orchestrator.outputs import OutputPublisher
        runner = self.root / "delayed_runner.py"
        runner.write_text("import time,runpy\ntime.sleep(.3)\nrunpy.run_path(" + repr(str(FAKE)) +
                          ",run_name='__main__')\n", encoding="utf-8")
        started, active = threading.Event(), threading.Event()
        original = OutputPublisher.update

        def slow_update(publisher, index, *, final=False):
            if threading.current_thread() is not threading.main_thread() and not started.is_set():
                started.set()
                active.set()
                try:
                    time.sleep(6)
                finally:
                    active.clear()
            if final:
                self.assertFalse(active.is_set(), "Final publication overlapped an unfinished background update")
            return original(publisher, index, final=final)

        environment = dict(self.env, VELA_EXAMPLE_RUNNER=str(runner))
        with patch.dict(os.environ, environment), patch.object(OutputPublisher, "update", slow_update):
            self.assertEqual(main(["device", str(self.cmd), "-q"]), 0)
        self.assertTrue(started.is_set(), "Slow background publication was not exercised")
        self.assertTrue((self.root / "final_des.vtu").is_file())
        self.assertFalse((self.work / ".lock").exists())

    def test_unchanged_csv_is_published_once_until_final_read(self):
        from vela_sdevice.orchestrator.outputs import OutputPublisher

        stage = {
            "csv": "stage00.csv", "current_file": "current_des.plt",
            "electrode_voltages_start_V": {"anode": 0.},
            "electrode_voltages_goal_V": {"anode": 1.},
            "plt_write": "all_accepted_steps", "vtu_write": {"series": []},
        }
        result = SimpleNamespace(stage_plan={"stages": [stage], "electrodes": [{"name": "anode"}]},
                                 decks=[{"sweep": {"contact": "anode"}}])
        publisher = OutputPublisher(self.root, self.root, result)
        path = self.root / stage["csv"]
        path.write_text("bias_V,current_total_A_per_um\n0,1e-12\n")
        with patch.object(publisher, "_publish_current", wraps=publisher._publish_current) as publish:
            self.assertEqual(len(publisher.update(0)), 1)
            self.assertEqual(len(publisher.update(0)), 1)
            self.assertEqual(publish.call_count, 1)
            with path.open("a") as stream:
                stream.write("1,2e-12")
            self.assertEqual(len(publisher.update(0)), 1)
            with path.open("a") as stream:
                stream.write("\n")
            self.assertEqual(len(publisher.update(0)), 2)
            self.assertEqual(len(publisher.update(0, final=True)), 2)
            self.assertEqual(publish.call_count, 4)
        self.assertEqual(len(publisher.completed_rows[stage["current_file"]]), 2)
        _, rows = read_plt_document((self.root / stage["current_file"]).read_text())
        self.assertEqual([row["time"] for row in rows], [0., 1.])

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=REPO / "build", prefix="wp5_test_")
        self.root = Path(self.temp.name)
        self.cmd = inputs(self.root)
        self.work = self.root / "device_des_vela"
        self.count = self.root / "count.jsonl"
        self.env = dict(os.environ, PYTHONPATH=str(REPO / "scripts"), VELA_EXAMPLE_RUNNER=str(FAKE),
                        FAKE_RUNNER_COUNT_FILE=str(self.count), FAKE_RUNNER_MODE="ok")

    def tearDown(self):
        self.temp.cleanup()

    def call(self, *args, mode="ok", target=None, env=None):
        environment = dict(self.env, FAKE_RUNNER_MODE=mode)
        environment.update(env or {})
        profile=[]
        runner=environment.get("VELA_EXAMPLE_RUNNER")
        if runner and Path(runner).resolve()!=FAKE.resolve() and "--physical-constants-profile" not in args:
            # Exercise explicit historical builds as well as the fresh default;
            # never feed compatible materials to a CODATA-only binary.
            profile=self.runner_constants_args(runner)
        return subprocess.run([sys.executable, "-m", "vela_sdevice", "device", *(args or [str(target or self.cmd)]), *profile],
                              cwd=self.root, env=environment, capture_output=True, text=True, encoding="utf-8", timeout=30)

    def runner_constants_args(self, runner):
        done=subprocess.run([runner,"--physical-constants"],capture_output=True,text=True,check=True,timeout=10)
        return ["--physical-constants-profile",json.loads(done.stdout)["profile"]]

    def assert_exit(self, code, done):
        self.assertEqual(done.returncode, code, done.stdout + done.stderr)

    def calls(self):
        return [json.loads(line) for line in self.count.read_text().splitlines()] if self.count.exists() else []

    def test_arithmetic_precision_dispatches_before_first_solve(self):
        source = self.cmd.read_text()
        missing = str(self.root / "missing-runner")
        for setting in ("", "-ExtendedPrecision"):
            with self.subTest(setting=setting):
                self.cmd.write_text(source.replace("Solve {", f"Math {{ {setting} }} Solve {{", 1))
                done = self.call(str(self.cmd), "-i", env={"VELA_EXAMPLE_RUNNER_EXTENDED80": missing})
                self.assert_exit(0, done)
                manifest = json.loads((self.work / "run_manifest.json").read_text())
                self.assertEqual(manifest["stages"][0]["requested_arithmetic_precision"], "double")
                self.assertFalse(any(s["require_end_to_end"] for s in manifest["arithmetic"]["planned_stages"]))

        self.count.unlink()
        extended = self.root / "extended.py"
        mixed = {"arithmetic_precision": "extended80", "radix": 2, "significand_bits": 64,
                 "max_exponent": 16384, "scope": "mixed", "end_to_end_linear_backends": []}
        for setting, capabilities, backend, expected in (
            ("ExtendedPrecision", None, "", "missing executable"),
            ("ExtendedPrecision", mixed | {"significand_bits": 53}, "", "80-bit arithmetic format"),
            ("ExtendedPrecision(80)", mixed, "pardiso_mkl", "VELA_LINEAR_SOLVER=pardiso_mkl"),
            ("ExtendedPrecision(80)", mixed, "", "reports mixed precision"),
            ("ExtendedPrecision(80)", mixed | {"scope": "end_to_end"}, "", "Frontend input conversion is binary64"),
        ):
            with self.subTest(setting=setting, expected=expected):
                self.cmd.write_text(source.replace("Solve {", f"Math {{ {setting} }} Solve {{", 1))
                extended.write_text(
                    "import sys\nassert sys.argv[1:] == ['--arithmetic-capabilities']\n"
                    f"print({json.dumps(capabilities)!r})\n", encoding="utf-8")
                done = self.call(str(self.cmd), env={
                    "VELA_EXAMPLE_RUNNER_EXTENDED80": str(extended) if capabilities else missing,
                    "VELA_LINEAR_SOLVER": backend,
                })
                self.assertNotEqual(done.returncode, 0)
                self.assertIn(expected, done.stderr)
                manifest = json.loads((self.work / "run_manifest.json").read_text())
                self.assertEqual(manifest["arithmetic"]["capability_check"], "failed")
                self.assertEqual(manifest["stages"], [])
                self.assertEqual(self.calls(), [])
                if capabilities:
                    self.assertEqual(manifest["arithmetic"]["reported_capabilities"], capabilities)
                    stages = manifest["arithmetic"]["planned_stages"]
                    self.assertGreater(len(stages), 1)
                    self.assertTrue(all(s["requested"] == "extended80" and s["require_end_to_end"]
                                        for s in stages))
                    input_check = manifest["arithmetic"]["frontend_input_check"]
                    self.assertEqual(input_check["status"], "blocked")
                    self.assertTrue(all(s["binary64_consumers"] for s in input_check["stages"]))

    def test_strict_arithmetic_requires_coverage_of_each_selected_stage(self):
        from vela_sdevice.orchestrator import DeviceRun, RunError

        required = ["mesh", "doping", "materials", "temperature", "contacts", "sweep", "ialmob", "state_archive"]
        complete = {"arithmetic_precision": "extended80", "radix": 2, "significand_bits": 64,
                    "max_exponent": 16384, "scope": "end_to_end",
                    "end_to_end_linear_backends": ["sparselu_extended80"],
                    "end_to_end_input_consumers": required}

        def verify(capabilities, missing):
            stages = [{"index": 0, "required_consumers": required[:-1], "binary64_consumers": [],
                       "missing_inputs": []},
                      {"index": 1, "required_consumers": required, "binary64_consumers": [],
                       "missing_inputs": missing}]
            coverage = {"status": "blocked", "stages": stages}
            audit = {"linear_solver_override": None, "frontend_input_check": coverage}
            device = SimpleNamespace(manifest={"arithmetic": audit}, work=self.work)
            with patch("vela_sdevice.orchestrator.arithmetic_capabilities", return_value=capabilities):
                DeviceRun._require_extended80(device, ["runner"])
            return audit

        audit = verify(complete, [])
        self.assertEqual(audit["capability_check"], "passed")
        self.assertEqual(audit["frontend_input_check"]["status"], "passed")
        with self.assertRaisesRegex(RunError, "Frontend input conversion is binary64"):
            verify(complete | {"end_to_end_input_consumers": required[:-1]}, [])
        with self.assertRaisesRegex(RunError, "Frontend input conversion is binary64"):
            verify(complete, ["decimal_stage_controls"])
        with self.assertRaisesRegex(RunError, "no end-to-end extended80 linear backend"):
            verify(complete | {"end_to_end_linear_backends": []}, [])

    def test_scoped_extended80_capability_does_not_enable_unrelated_stages(self):
        from vela_sdevice.orchestrator import DeviceRun, RunError

        profile = "slot_population_bvds_v1"
        capabilities = {"arithmetic_precision": "extended80", "radix": 2, "significand_bits": 64,
                        "max_exponent": 16384, "scope": "mixed", "end_to_end_linear_backends": [],
                        "scoped_extended80_profiles": {profile: {
                            "linear_backends": ["sparse_lu"],
                            "input_consumers": ["mesh", "state_archive"]}}}

        def verify(profiles, backend=None):
            coverage = {"status": "blocked", "stages": [
                {"profile": value, "required_consumers": ["mesh", "state_archive"],
                 "binary64_consumers": [], "missing_inputs": []} for value in profiles]}
            audit = {"linear_solver_override": backend, "frontend_input_check": coverage}
            device = SimpleNamespace(manifest={"arithmetic": audit}, work=self.work)
            with patch("vela_sdevice.orchestrator.arithmetic_capabilities", return_value=capabilities):
                DeviceRun._require_extended80(device, ["runner"])
            return audit["capability_check"]

        self.assertEqual(verify([profile, profile], "sparse_lu"), "passed")
        with self.assertRaisesRegex(RunError, "reports mixed precision"):
            verify([profile, None])
        with self.assertRaisesRegex(RunError, "unsupported for end-to-end extended80"):
            verify([profile], "pardiso_mkl")

    def test_extended80_stage_handoff_requires_archive_identity_before_publication(self):
        literal = "1.000000000000000000867361737988403547205962240695953369140625"
        text = self.cmd.read_text().replace("Voltage=0", "Voltage=1", 1)
        self.cmd.write_text(text.replace("Solve {", "Math { ExtendedPrecision(80) } Solve {", 1))
        runner = self.root / "future_extended.py"
        caps = {"arithmetic_precision": "extended80", "radix": 2, "significand_bits": 64,
                "max_exponent": 16384, "scope": "end_to_end",
                "end_to_end_linear_backends": ["sparselu_extended80"],
                "end_to_end_input_consumers":
                    ["mesh", "doping", "materials", "temperature", "contacts", "sweep", "state_archive"]}
        runner.write_text(f'''
import csv, hashlib, json, os, runpy, sys
from pathlib import Path
if sys.argv[1:] == ["--arithmetic-capabilities"]:
    print({json.dumps(caps)!r})
    raise SystemExit(0)
if sys.argv[1:] == ["--verify-extended80-stage-state"]:
    request = json.load(sys.stdin)
    assert request["device_voltage_V"] == {literal!r}
    assert request["supply_voltage_V"] == {literal!r}
    if os.environ.get("EXTENDED_STATE_MODE") == "diagnostic_in_file":
        print("Archive has a diagnostic role", file=sys.stderr)
        raise SystemExit(2)
    print(json.dumps({{"verified": True, "arithmetic_precision": "extended80"}}))
    raise SystemExit(0)
code = runpy.run_path({str(FAKE)!r})["main"]()
deck = Path(sys.argv[sys.argv.index("--config") + 1])
config = json.loads(deck.read_text())
archive = Path(config["sweep"]["write_state_file"])
csv_name = Path(config["output_csv"])
with csv_name.open() as stream:
    count = len(list(csv.DictReader(stream)))
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
mode = os.environ.get("EXTENDED_STATE_MODE")
record = {{"schema": "vela.extended80-stage-state.v1",
           "role": "diagnostic" if mode == "diagnostic" else "accepted",
           "archive": archive.name, "accepted_points": count, "deck_sha256": sha(deck),
           "node_count": 4, "mesh_sha256": "a"*64, "input_sha256": "b"*64,
           "csv_sha256": sha(csv_name),
           "archive_sha256": "0" * 64 if mode == "wrong_archive" else sha(archive),
           "voltages_V": {{"anode": {literal!r}, "cathode": "0"}},
           "outer_voltage_V": {literal!r},
           "currents_A_per_um": {{"anode": "1e-7", "cathode": "-1e-7"}}}}
print(json.dumps({{"config": deck.name, "arithmetic_precision": "extended80",
                  "accepted_extended80_state": record}}))
raise SystemExit(code)
''', encoding="utf-8")
        environment = {"VELA_EXAMPLE_RUNNER_EXTENDED80": str(runner), "VELA_LINEAR_SOLVER": ""}
        for invalid in ("diagnostic", "wrong_archive", "diagnostic_in_file"):
            with self.subTest(invalid=invalid):
                done = self.call(str(self.cmd), "-i", env=environment | {"EXTENDED_STATE_MODE": invalid})
                self.assertNotEqual(done.returncode, 0)
                self.assertIn("Cannot verify accepted extended80 stage state", done.stderr)
                self.assertFalse((self.root / "current_des.plt").exists())
        done = self.call(str(self.cmd), "-i", env=environment)
        self.assert_exit(0, done)
        audit = json.loads((self.work / "run_manifest.json").read_text())
        self.assertEqual(audit["arithmetic"]["capability_check"], "passed")
        stage = audit["stages"][0]
        self.assertTrue(stage["goal_reached"])
        self.assertEqual(stage["accepted_extended80_state"]["voltages_V"]["anode"], literal)
        self.assertTrue((self.root / "current_des.plt").exists())

    def test_strict_load_checks_the_named_parent_before_starting_runner(self):
        from hashlib import sha256
        from vela_sdevice.orchestrator import DeviceRun
        from vela_sdevice.orchestrator.extended_stage import AcceptedExtendedStage
        from vela_sdevice.orchestrator.support import RunError

        device = DeviceRun.__new__(DeviceRun)
        device.work = self.root
        saved, latest = "saved.h5", "latest.h5"
        checkpoint = self.root / saved
        checkpoint.write_bytes(b"the accepted state")
        device.result = SimpleNamespace(
            stage_plan={"stages": [{"index": 0, "output_state": saved},
                                   {"index": 1, "output_state": latest},
                                   {"index": 2, "input_state": saved}]},
            decks=[{}, {}, {"solver": {"require_end_to_end_precision": True}}],
        )
        parent = AcceptedExtendedStage(
            {"archive": saved, "archive_sha256": sha256(checkpoint.read_bytes()).hexdigest()},
            {}, None, {},
        )
        device.extended_states = {1: parent}
        with self.assertRaisesRegex(RunError, "missing verified accepted parent"):
            device._execute(2, ["runner"])
        device.extended_states[0] = parent
        checkpoint.write_bytes(b"replaced contents")
        with self.assertRaisesRegex(RunError, "parent archive changed after acceptance"):
            device._execute(2, ["runner"])

    def test_arithmetic_input_keeps_original_numeric_tokens(self):
        from vela_sdevice.cmd_parser import parse_cmd_file
        from vela_sdevice.orchestrator.arithmetic_input import collect_input

        literal = "1.000000000000000000867361737988403547205962240695953369140625"
        temperature = "300.0000000000000002220446049250313080847263336181640625"
        source = self.cmd.read_text().replace("Voltage=0", f"Voltage={literal}", 1)
        self.cmd.write_text(source.replace("Solve {", f"Physics {{Temperature={temperature}}} Solve {{", 1))
        par = self.root / "models.par"
        par.write_text(par.read_text().replace("4.05", literal, 1))
        grd = self.root / "mesh.grd"
        head, vertices, tail = grd.read_text().partition("Vertices (")
        self.assertTrue(vertices)
        grd.write_text(head + vertices + tail.replace("1.000000000000000e+00", literal, 1))
        dat = self.root / "mesh.dat"
        dat.write_text(dat.read_text().replace("1.000000000000000e+15", literal, 1))
        sidecar = self.cmd.with_suffix(".vela.json")
        sidecar.write_text('{"~field/":' + literal + ', "description":"' + literal + '"}')
        captures = {}
        ir = parse_cmd_file(self.cmd, arithmetic_input=captures)
        result = collect_input(self.cmd, captures, par, grd, sidecar, "sentaurus_u2022_12")
        self.assertTrue(ir["verdict"]["accepted"])
        self.assertEqual(result["temperature_K"], temperature)
        self.assertEqual(ir["electrodes"][0]["voltage_V"], 1.0)
        self.assertEqual(result["sources"]["sidecar"]["numbers"], {"/~0field~1": literal})
        for role in ("cmd", "par", "grd", "dat"):
            with self.subTest(role=role):
                self.assertIn(literal, result["sources"][role]["numbers"].values())

    def test_tdr_arithmetic_input_labels_binary64_csv(self):
        from vela_sdevice.orchestrator.arithmetic_input import collect_tdr_export
        tdr = self.root / "sample.tdr"
        tdr.write_bytes(b"bound native input")
        (self.root / "nodes.csv").write_text("id,x_um,y_um\n0,0.10000000000000001,0\n")
        (self.root / "doping.csv").write_text("node_id,donors_cm3,acceptors_cm3\n0,1e15,0\n")
        sources = collect_tdr_export(tdr, self.root)
        self.assertEqual(sources["tdr_coordinates"]["numbers"]["0/x_um"], "0.10000000000000001")
        self.assertEqual(sources["tdr_coordinates"]["number_encoding"], "binary64_roundtrip")
        self.assertEqual(sources["tdr_doping"]["original_tdr_sha256"],
                         sources["tdr_coordinates"]["original_tdr_sha256"])

    def test_doping_declaration_and_explicit_net_recompute_reach_mesh_import(self):
        original = self.cmd.read_text()
        source = original.replace('Grid="mesh.grd"', 'Grid="mesh.grd" Doping="mesh.dat"')
        self.cmd.write_text(source.replace("Solve {", "Math {ComputeDopingConcentration} Solve {", 1))
        self.assert_exit(0, self.call(str(self.cmd), "--check"))
        manifest = json.loads((self.work / "run_manifest.json").read_text())
        self.assertTrue(manifest["doping"]["compute_net"])
        self.assertEqual(manifest["doping"]["regions"]["R.Si"]["mode"], "species")
        self.cmd.write_text(source.replace('Doping="mesh.dat"', 'Doping="different.tdr"'))
        self.assert_exit(4, self.call(str(self.cmd), "--check"))
        failure = json.loads((self.work / "fail_report.json").read_text())
        self.assertIn("separate mesh input", failure["message"])

    def test_slot_arithmetic_input_reads_bound_files_and_stage_lexemes(self):
        from hashlib import sha256
        from vela_sdevice.orchestrator.arithmetic_input import collect_slot_input
        from vela_sdevice.orchestrator.support import RunError

        literal = "1.000000000000000000867361737988403547205962240695953369140625"
        files = {
            "mesh": self.root / "mesh.json",
            "doping": self.root / "doping.csv",
            "deck": self.root / "reference.json",
        }
        files["mesh"].write_text('{"nodes":[{"id":0,"x":' + literal + ',"y":0}]}')
        files["doping"].write_text(f"node_id,donors_cm3,acceptors_cm3\n0,{literal},0\n")
        files["deck"].write_text('{"solver":{"mobility":{"electron_field_beta":' + literal + '}}}')
        material = self.root / "materials.json"
        material.write_text('{"materials":[{"name":"Si","ni":1,"decimal_semiconductor":{"temperature_K":"300"}},'
                            '{"name":"PolySilicon","ni":2}]}')
        stage = self.root / "simulation_stage00.json"
        stage.write_text(json.dumps({"decimal_stage_controls": {"sweep": {"vtk_bias_points": ["1", literal]}}}))
        binding = {"inputs": {key: {"sha256": sha256(path.read_bytes()).hexdigest()} for key, path in files.items()},
                   "decimal_material": {"staged_sha256": sha256(material.read_bytes()).hexdigest()}}
        sources, materials = collect_slot_input(self.root, binding, files["deck"], [stage])
        self.assertEqual(sources["slot_bound_coordinates"]["numbers"]["0/x"], literal)
        self.assertEqual(sources["slot_bound_doping"]["numbers"]["0/donors_cm3"], literal)
        self.assertEqual(sources["slot_reference_deck"]["numbers"]["/solver/mobility/electron_field_beta"], literal)
        self.assertEqual(sources[stage.stem]["numbers"]["/sweep/vtk_bias_points/1"], literal)
        self.assertEqual(materials, {"Si": {"temperature_K": "300"}})
        self.assertIn("/materials/1/ni", sources["slot_staged_materials"]["numbers"])
        files["doping"].write_text("node_id,donors_cm3,acceptors_cm3\n0,1,0\n")
        with self.assertRaisesRegex(RunError, "SLOT doping differs"):
            collect_slot_input(self.root, binding, files["deck"], [stage])

    def test_ok_outputs_log_and_manifest(self):
        done = self.call()
        self.assert_exit(0, done)
        self.assertEqual(done.stdout.count(" V accepted"), 4)
        log = (self.root / "run_des.log").read_text()
        for section in ("[HEADER]", "[CMD]", "[DECISION]", "[STAGE]", "[SUMMARY]"):
            self.assertIn(section, log)
        manifest = json.loads((self.work / "run_manifest.json").read_text())
        self.assertTrue(all(s["goal_reached"] for s in manifest["stages"]))
        names, rows = read_plt_document((self.root / "current_des.plt").read_text())
        self.assertEqual(names, ["anode", "cathode"])
        self.assertEqual(rows[-1]["cathode"]["total_current_A_per_um"], -1e-7)
        self.assertTrue((self.root / "final_des.vtu").is_file())
        pvd = ET.parse(self.root / "final_des.pvd")
        for dataset in pvd.findall(".//DataSet"):
            self.assertTrue((self.root / dataset.attrib["file"]).is_file())
        self.assertFalse((self.work / ".lock").exists())

    def test_check_and_initial_only(self):
        self.assert_exit(0, self.call(str(self.cmd), "--check"))
        self.assertEqual(self.calls(), [])
        self.assertTrue(json.loads((self.work / "run_manifest.json").read_text())["check_only"])
        checked = json.loads((self.work / "mesh_qualification.json").read_text())["qualification"]["electrodes_checked"]
        self.assertEqual([entry["matched"] for entry in checked], [True, True])
        self.assert_exit(0, self.call(str(self.cmd), "-i"))
        self.assertEqual(len(self.calls()), 1)
        self.assertTrue((self.root / "final_des.vtu").exists())

    def test_bvmethods_cmd_default_constants_and_explicit_override(self):
        text = self.cmd.read_text().replace('{Name="anode" Voltage=0}', '{Name="anode" Voltage=0 Resistor=1e7}')
        text = text.replace("Solve {", """
            Physics { EffectiveIntrinsicDensity(OldSlotboom) Fermi
              Mobility(DopingDep HighFieldSaturation(GradQuasiFermi) Enormal)
              Recombination(SRH(DopingDep) Avalanche(GradQuasiFermi)) }
            Plot { Potential eDensity hDensity SpaceCharge Doping }
            Math { Extrapolate AvalDerivatives Iterations=20 }
            Solve {""")
        self.cmd.write_text(text)
        self.cmd.with_suffix(".vela.json").write_text(json.dumps({
            "schema": "vela.sdevice_sidecar.v1", "bvmethods_operator_profile": "bvmethods_tri3_qf"
        }))
        for args, expected in (
            ([], "sentaurus_u2022_12"),
            (["--physical-constants-profile", "codata2018"], "codata2018"),
        ):
            with self.subTest(constants=expected):
                self.assert_exit(0, self.call(str(self.cmd), "--check", "--allow-lossy", *args))
                report = json.loads((self.work / "materials_report.json").read_text())
                deck = json.loads((self.work / "simulation_stage01.json").read_text())
                self.assertEqual(report.get("physical_constants_profile", "codata2018"), expected)
                self.assertEqual(deck.get("physical_constants_profile", "codata2018"), expected)
                self.assertEqual(deck["sweep"]["external_circuit"]["max_iterations"], 20)
        self.assertEqual(self.calls(), [])

    def test_explicit_plot_names_cannot_be_overwritten_by_default_series(self):
        self.cmd.write_text(self.cmd.read_text().replace('Coupled {Poisson Electron Hole}}',
            'Coupled {Poisson Electron Hole} Plot(FilePrefix="final" Time=(0;0.5;1)) }', 1))
        self.assert_exit(0, self.call())
        for index, t in enumerate((0., .5, 1.)):
            metadata = vtu_metadata(self.root / f"final_{index:04d}_des.vtu")
            self.assertEqual(metadata["Stage"], 1)
            self.assertEqual(metadata["t"], t)
        self.assertEqual(len(list(self.root.glob("final_????_des.vtu"))), 7)

    def test_cli_identity_usage_parameters_and_filename_resolution(self):
        for option in ("-h", "-v", "-P"):
            self.assert_exit(0, self.call(option))
        self.assertEqual((self.root / "models.par").read_bytes(), default_models_par_path().read_bytes())
        self.assert_exit(0, self.call(str(self.cmd.with_suffix("")), "--check"))
        self.assert_exit(0, self.call(str(self.root / "device"), "--check"))
        for args in (("-d",), (str(self.cmd), "--threads", "0"), ()):
            command = [sys.executable, "-m", "vela_sdevice", "device", *args]
            done = subprocess.run(command, cwd=self.root, env=self.env, capture_output=True, text=True)
            self.assert_exit(2, done)

    def test_input_failures_exit_1_and_3(self):
        self.cmd.write_text(self.cmd.read_text() + "\n#include \"other.cmd\"\n")
        self.assert_exit(1, self.call())
        self.assertEqual(json.loads((self.work / "fail_report.json").read_text())["exit_code"], 1)
        inputs(self.root)
        (self.root / "mesh.grd").unlink()
        self.assert_exit(3, self.call())
        self.assert_exit(3, self.call(str(self.root / "missing.cmd")))

    def test_reserved_debug_export_fails_before_any_side_effects(self):
        before = {p.name: p.read_bytes() for p in self.root.iterdir() if p.is_file()}
        for args in (("-d",), (str(self.cmd), "-d"), ("missing.cmd", "-d"),
                     ("-d", "--check"), ("-d", "-P"), ("-d", "-v")):
            with self.subTest(args=args):
                done = self.call(*args)
                self.assert_exit(2, done)
                self.assertIn("E202: -d", done.stderr)
                self.assertIn("not implemented yet", done.stderr)
                self.assertEqual(done.stdout, "")
                self.assertEqual(self.calls(), [])
                self.assertFalse(self.work.exists())
                self.assertEqual(set(p.name for p in self.root.iterdir()), set(before))
                self.assertEqual({p.name: p.read_bytes() for p in self.root.iterdir()}, before)
        help_output = self.call("-h")
        self.assert_exit(0, help_output)
        self.assertIn("not implemented", help_output.stdout)

    def test_all_reserved_options_report_without_starting_or_writing(self):
        from vela_sdevice.__main__ import RESERVED_OPTIONS
        before = {p.name: p.read_bytes() for p in self.root.iterdir() if p.is_file()}
        cases = [(flag,) for flag in RESERVED_OPTIONS]
        cases += [("-P:Silicon",), ("-L:Silicon/Oxide",), ("-P:AlGaAs:0.3",),
                  ("--max_threads=4",), ("-ver", "T-2022.03"),
                  ("-rel", "T-2022.03-SP2"), ("-P", str(self.cmd))]
        for args in cases:
            with self.subTest(args=args):
                done = self.call(*args)
                self.assert_exit(2, done)
                self.assertIn("not implemented yet", done.stderr)
                self.assertIn("No simulation was started", done.stderr)
                self.assertEqual(done.stdout, "")
                self.assertEqual(set(p.name for p in self.root.iterdir()), set(before))
                self.assertEqual({p.name: p.read_bytes() for p in self.root.iterdir()}, before)

    def test_sd_sim_direct_zipapp_interface(self):
        archive = build_archive(self.root / "direct.pyz")
        for args, code, text in ((("-h",), 0, "usage: sd-sim"),
                                 (("-v",), 0, "Vela"),
                                 (("-d",), 2, "not implemented yet"),
                                 (("-P:Silicon",), 2, "not implemented yet"),
                                 (("--not-a-real-option",), 2, "E202"),
                                 ((str(self.cmd), "--check"), 0, ""),
                                 ((str(self.cmd),), 0, "")):
            with self.subTest(args=args):
                done = subprocess.run([sys.executable, "-I", "-S", str(archive), "sd-sim", *args],
                                      cwd=self.root, env=self.env, capture_output=True,
                                      text=True, encoding="utf-8", timeout=30)
                self.assert_exit(code, done)
                self.assertIn(text, done.stdout + done.stderr)

    def test_sidecar_and_runner_rejection_exit_4(self):
        self.cmd.with_suffix(".vela.json").write_text('{"schema":"vela.sdevice_sidecar.v1","contacts":[]}')
        self.assert_exit(4, self.call())
        self.assertEqual(self.calls(), [])
        self.cmd.with_suffix(".vela.json").unlink()
        self.assert_exit(4, self.call(mode="reject"))
        report = json.loads((self.work / "fail_report.json").read_text())
        self.assertIn("deck rejected", report["details"]["stderr"])

    def test_stage_log_profile_and_verbose_precedence(self):
        for profile, verbose in (("minimal", False), ("default", False), ("debug", False), ("minimal", True)):
            with self.subTest(profile=profile, verbose=verbose):
                self.count.unlink(missing_ok=True)
                self.cmd.with_suffix(".vela.json").write_text(json.dumps({
                    "schema": "vela.sdevice_sidecar.v1", "apply_to": [1],
                    "runtime_log": {"profile": profile}}))
                self.assert_exit(0, self.call(str(self.cmd), *(["--verbose"] if verbose else [])))
                calls = self.calls()
                self.assertEqual([c["log_profile"] for c in calls],
                                 ["debug", "debug"] if verbose else ["default", profile])
                self.assertTrue(all(c["deck"]["runtime_log"]["append"] for c in calls))
                self.assertEqual(calls[1]["deck"]["runtime_log"]["profile"], profile)

    def test_real_runner_stage_log_profile(self):
        runner = os.environ.get("VELA_EXAMPLE_RUNNER")
        if not runner or not Path(runner).is_file():
            self.skipTest("VELA_EXAMPLE_RUNNER unavailable")
        for profile, verbose in (("minimal", False), ("debug", False), ("minimal", True)):
            with self.subTest(profile=profile, verbose=verbose):
                self.cmd.with_suffix(".vela.json").write_text(json.dumps({
                    "schema": "vela.sdevice_sidecar.v1", "apply_to": [1],
                    "runtime_log": {"profile": profile}}))
                self.assert_exit(0, self.call(str(self.cmd), *(["--verbose"] if verbose else []),
                                            env={"VELA_EXAMPLE_RUNNER": runner}))
                expected = ["debug", "debug"] if verbose else ["default", profile]
                for index, actual in enumerate(expected):
                    log = (self.work / f"stage{index:02d}_runner.log").read_text()
                    self.assertIn("profile: " + actual, log)
                manifest = json.loads((self.work / "run_manifest.json").read_text())
                self.assertTrue(all(stage["goal_reached"] for stage in manifest["stages"]))

    def test_solver_failure_exit_5(self):
        self.assert_exit(5, self.call(mode="fail"))
        self.assertEqual(json.loads((self.work / "fail_report.json").read_text())["message_id"], "E501")
        self.assert_exit(5, self.call(mode="crash"))

    def test_output_failure_exit_6(self):
        self.assert_exit(6, self.call(mode="bad_output"))
        self.assertFalse((self.work / ".lock").exists())

    def test_manifest_write_failure_releases_lock(self):
        self.work.mkdir()
        (self.work / "run_manifest.json").mkdir()
        done = self.call(str(self.cmd), "--check")
        self.assert_exit(6, done)
        self.assertNotIn("Traceback", done.stderr)
        self.assertEqual(json.loads((self.work / "fail_report.json").read_text())["exit_code"], 6)
        self.assertFalse((self.work / ".lock").exists())

    def test_live_lock_conflict_and_stale_lock_reclamation(self):
        self.work.mkdir()
        lock = self.work / ".lock"
        owner = {"pid": os.getpid(), "host": socket.gethostname(), "time": time.time(), "token": "owned"}
        lock.write_text(json.dumps(owner))
        self.assert_exit(7, self.call())
        self.assertEqual(json.loads(lock.read_text()), owner)
        owner["pid"] = 2147483647
        lock.write_text(json.dumps(owner))
        self.assert_exit(0, self.call())
        self.assertIn("W701", (self.root / "run_des.log").read_text())
        self.assertFalse(lock.exists())

    def test_partial_rebase_and_exit_on_failure(self):
        self.cmd = inputs(self.root, multistage=True)
        self.assert_exit(0, self.call(mode="partial"))
        self.assertEqual(self.calls()[-1]["deck"]["sweep"]["start"], .5)
        runtime = json.loads((self.work / "run_manifest.json").read_text())["stages"]
        self.assertFalse(runtime[1]["goal_reached"])
        self.assertTrue(runtime[2]["goal_reached"])
        self.count.unlink()
        self.assert_exit(5, self.call(str(self.cmd), "--exit-on-failure", mode="partial"))
        self.assertEqual(len(self.calls()), 2)

    def test_threads_quiet_binary_and_final_only(self):
        done = self.call(str(self.cmd), "-q", "--threads", "2", "--vtu-binary", "--plot-final-only", "--no-keep-intermediate")
        self.assert_exit(0, done)
        self.assertEqual(done.stdout, "")
        self.assertTrue(all(c["threads"] == "2" and c["hdf5_locking"] == "FALSE" for c in self.calls()))
        self.assertIn(b"\x00", (self.root / "final_des.vtu").read_bytes())
        self.assertEqual(len(list(self.root.glob("*.vtu"))), 1)
        self.assertEqual(vtu_metadata(self.root / "final_des.vtu")["Stage"], 1)
        self.assertFalse(list(self.work.glob("stage*.vtu")))

    def test_signal_retains_partial_plt_and_removes_lock_exit_8(self):
        env = dict(self.env, FAKE_RUNNER_MODE="hang")
        proc = subprocess.Popen([sys.executable, "-m", "vela_sdevice", "device", str(self.cmd)],
            cwd=self.root, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0)
        try:
            deadline = time.monotonic() + 10
            while not (self.root / "current_des.plt").exists() and time.monotonic() < deadline:
                time.sleep(.05)
            self.assertTrue((self.root / "current_des.plt").exists())
            try:
                os.kill(proc.pid, signal.CTRL_BREAK_EVENT if os.name == "nt" else signal.SIGTERM)
            except OSError as error:
                self.skipTest(f"CTRL_BREAK_EVENT unavailable in this console: {error}")
            stdout, stderr = proc.communicate(timeout=10)
            self.assertEqual(proc.returncode, 8, stdout + stderr)
            self.assertTrue(read_plt_document((self.root / "current_des.plt").read_text())[1])
            self.assertFalse((self.work / ".lock").exists())
        finally:
            if proc.poll() is None and self.calls():
                # A console without CTRL_BREAK support must not leave the hang
                # fixture running after the test is explicitly skipped.
                if os.name == "nt":
                    import ctypes
                    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
                    kernel.OpenProcess.restype = ctypes.c_void_p
                    kernel.TerminateProcess.argtypes = [ctypes.c_void_p, ctypes.c_uint]
                    kernel.CloseHandle.argtypes = [ctypes.c_void_p]
                    child = kernel.OpenProcess(1, False, self.calls()[-1]["pid"])
                    if child:
                        kernel.TerminateProcess(child, 1); kernel.CloseHandle(child)
                else:
                    os.kill(self.calls()[-1]["pid"], signal.SIGTERM)
            if proc.poll() is None:
                proc.terminate(); proc.wait(timeout=5)
            proc.stdout.close(); proc.stderr.close()

    def test_hermetic_zipapp_runs_version_defaults_and_check(self):
        archive = build_archive(self.root / "bin/vela_sdevice.pyz")
        env = dict(self.env, PYTHONPATH="")
        for args in (("-v",), ("-P",), (str(self.cmd), "--check"), (str(self.cmd),)):
            done = subprocess.run([sys.executable, "-I", "-S", str(archive), "device", *args], cwd=self.root,
                                  env=env, capture_output=True, text=True, encoding="utf-8", timeout=30)
            self.assert_exit(0, done)

    def test_zipapp_preserves_numbered_nonzero_exit_codes(self):
        archive = build_archive(self.root / "exit_codes.pyz")
        self.cmd.write_text(self.cmd.read_text() + '\nThermode {{Name="anode" Temperature=300}}\n')
        for arguments, expected, message in (([str(self.cmd),"--check"],1,"E105"),
                                            (["-d"],2,"not implemented yet"),
                                            (["--unsupported-option"],2,"E202"),
                                            (["missing_des.cmd"],3,"E305")):
            done = subprocess.run([sys.executable,"-I","-S",str(archive),"device",*arguments],cwd=self.root,
                env=self.env,capture_output=True,text=True,encoding="utf-8",timeout=30)
            self.assert_exit(expected,done)
            self.assertIn(message,done.stderr)

    def test_native_launcher_preserves_rejection_exit_code(self):
        launcher=os.environ.get("VELA_LAUNCHER")
        if not launcher or not Path(launcher).is_file():self.skipTest("VELA_LAUNCHER unavailable")
        self.cmd.write_text(self.cmd.read_text() + '\nThermode {{Name="anode" Temperature=300}}\n')
        env=dict(self.env,VELA_PYTHON=sys.executable);env.pop("VELA_SDEVICE_PYZ",None)
        done=subprocess.run([launcher,"device",str(self.cmd),"--check"],cwd=self.root,
            env=env,capture_output=True,text=True,encoding="utf-8",timeout=30)
        self.assert_exit(1,done);self.assertIn("E105",done.stderr)

    def test_real_runner_pn2d_fixture_in_unicode_directory(self):
        runner = os.environ.get("VELA_EXAMPLE_RUNNER")
        if not runner or not Path(runner).is_file():
            self.skipTest("VELA_EXAMPLE_RUNNER unavailable")
        root = self.root / ("space " + chr(0x4e2d) + chr(0x6587))
        root.mkdir()
        cmd = inputs(root, fixture=True)
        # Full fixture voltages and Plot Time controls, on the synthetic mesh.
        launcher = os.environ.get("VELA_LAUNCHER")
        if launcher and Path(launcher).is_file():
            env = dict(self.env, VELA_EXAMPLE_RUNNER=runner, PYTHONPATH="", VELA_PYTHON=sys.executable)
            env.pop("VELA_SDEVICE_PYZ", None)
            done = subprocess.run([launcher, "device", str(cmd), "--allow-lossy", *self.runner_constants_args(runner)], cwd=root,
                env=env, capture_output=True, text=True, encoding="utf-8", timeout=30)
        else:
            done = self.call(str(cmd), "--allow-lossy", env={"VELA_EXAMPLE_RUNNER": runner})
        self.assert_exit(0, done)
        self.assertTrue(list(root.glob("*.log")))
        self.assertTrue(list(root.glob("*.plt")))
        self.assertTrue(list(root.glob("*.vtu")))
        self.assertTrue(list(root.glob("*.pvd")))
        self.assertFalse((root / "device_des_vela/.lock").exists())

    def test_current_time_mapping_and_next_voltage_rebase(self):
        text = self.cmd.read_text().rstrip()
        self.cmd.write_text(text[:-1] + '''
            Set("anode" mode current)
            Quasistationary(Goal {Name="anode" Current=2e-7}) {
              Coupled {Poisson Electron Hole} CurrentPlot(Time=(0;0.5;1)) }
            Set("anode" mode voltage)
            Quasistationary(Goal {Name="anode" Voltage=3}) {Coupled {Poisson Electron Hole}} }''')
        self.assert_exit(0, self.call(str(self.cmd), "--check"))
        self.assertEqual(self.calls(), [])
        self.assert_exit(0, self.call())
        calls = self.calls()
        control = calls[2]["deck"]["sweep"]["voltage_to_current"]
        self.assertEqual(control["current_points_A_per_um"], [1e-7, 1.5e-7, 2e-7])
        self.assertEqual(calls[3]["deck"]["sweep"]["start"], 2)

    def test_saved_partial_voltage_is_restored_on_load(self):
        text = self.cmd.read_text().rstrip()
        self.cmd.write_text(text[:-1] + ' Save(FilePrefix="partial") Load(FilePrefix="partial") }')
        self.assert_exit(0, self.call(mode="partial"))
        self.assertEqual(self.calls()[2]["deck"]["sweep"]["start"], .5)
        self.assertEqual(self.calls()[2]["deck"]["sweep"]["initial_state_file"], "partial_des.h5")

    def test_native_launcher_uses_bundled_archive(self):
        launcher = os.environ.get("VELA_LAUNCHER")
        if not launcher or not Path(launcher).is_file():
            self.skipTest("VELA_LAUNCHER unavailable")
        archive = Path(launcher).parent / "bin/vela_sdevice.pyz"
        self.assertTrue(archive.is_file())
        env = dict(self.env, PYTHONPATH="", VELA_PYTHON=sys.executable)
        env.pop("VELA_SDEVICE_PYZ", None)
        done = subprocess.run([launcher, "device", str(self.cmd)], cwd=self.root, env=env,
            capture_output=True, text=True, encoding="utf-8", timeout=30)
        self.assert_exit(0, done)
        self.assertTrue((self.root / "final_des.vtu").is_file())

    def test_real_runner_current_save_load_chain(self):
        runner = os.environ.get("VELA_EXAMPLE_RUNNER")
        if not runner or not Path(runner).is_file():
            self.skipTest("VELA_EXAMPLE_RUNNER unavailable")
        text = self.cmd.read_text().replace('Voltage=1}', 'Voltage=0.001}')
        text = text.replace('Coupled {Poisson Electron Hole}\n',
                            'Coupled {Poisson Electron Hole} Save(FilePrefix="eq")\n', 1)
        self.cmd.write_text(text.rstrip()[:-1] + '''
          Set("anode" mode current)
          Quasistationary(Goal {Name="anode" Current=1e-7}) {
            Coupled {Poisson Electron Hole} CurrentPlot(Time=(0;0.5;1)) }
          Load(FilePrefix="eq") }''')
        self.cmd.with_suffix(".vela.json").write_text(json.dumps({
            "schema": "vela.sdevice_sidecar.v1", "apply_to": [0, 1, 2], "solver": {"method": "newton"}}))
        done = self.call(str(self.cmd), env={"VELA_EXAMPLE_RUNNER": runner})
        self.assert_exit(0, done)
        manifest = json.loads((self.work / "run_manifest.json").read_text())
        self.assertEqual(len(manifest["stages"]), 4)
        self.assertTrue(all(s["goal_reached"] for s in manifest["stages"]))
        self.assertEqual(manifest["stages"][-1]["reached_voltages_V"]["anode"], 0)
        self.assertTrue((self.work / "eq_des.h5").is_file())


if __name__ == "__main__":
    unittest.main()
