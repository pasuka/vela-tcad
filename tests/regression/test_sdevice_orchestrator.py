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
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from build_vela_sdevice_pyz import build_archive
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
        return subprocess.run([sys.executable, "-m", "vela_sdevice", "device", *(args or [str(target or self.cmd)])],
                              cwd=self.root, env=environment, capture_output=True, text=True, encoding="utf-8", timeout=30)

    def assert_exit(self, code, done):
        self.assertEqual(done.returncode, code, done.stdout + done.stderr)

    def calls(self):
        return [json.loads(line) for line in self.count.read_text().splitlines()] if self.count.exists() else []

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
            done = subprocess.run([launcher, "device", str(cmd), "--allow-lossy"], cwd=root,
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
