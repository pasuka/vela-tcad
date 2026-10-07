"""WP11 real-runner determinism, thread controls, and duplicate lock tests."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

from tests.regression.sdevice_package_support import (synthetic_cmd, isolated_environment,
    command, run, numeric_outputs)
REPO = Path(__file__).resolve().parents[2]
SUFFIX = ".exe" if os.name == "nt" else ""
LAUNCHER = Path(os.environ.get("VELA_LAUNCHER", REPO / "build-release" / ("vela" + SUFFIX)))
LINUX_REASON = "WP11 Windows scope is lock conflict only; ps/thread and numerical concurrency matrix run on Linux"

class ConcurrentDeviceTests(unittest.TestCase):
    def setUp(self):
        if not LAUNCHER.is_file():
            self.fail(f"WP11 launcher not built: {LAUNCHER}")
        self.temp = tempfile.TemporaryDirectory(prefix="wp11_concurrent_", dir=REPO / "build-release" if os.name == "nt" else None)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = isolated_environment()
        self.env["VELA_PYTHON"] = sys.executable
        # The explicit interpreter override also covers a relocated bundled
        # Python. Protect its shipped caches in launcher and protocol children.
        self.env["PYTHONDONTWRITEBYTECODE"] = "1"
        self.env.update(COLUMNS="80", LINES="24")

    def assert_ok(self, done):
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def baseline(self, cmd, stem, *args):
        self.assert_ok(run(LAUNCHER, cmd, *args, env=self.env))
        result = (cmd.parent / (stem + ".plt")).read_bytes(), numeric_outputs(cmd.parent, stem)
        # A missing parallel publication must not be masked by serial files.
        (cmd.parent / (stem + ".plt")).unlink()
        for pattern in (stem + "*.vtu", stem + "*.pvd"):
            for path in cmd.parent.glob(pattern):
                path.unlink()
        return result

    def parallel(self, commands, cap=1):
        processes = []
        samples = []
        logs = []
        try:
            for cmd, args in commands:
                log = tempfile.TemporaryFile()
                logs.append(log)
                processes.append(subprocess.Popen(command(LAUNCHER, cmd, *args), cwd=cmd.parent,
                    env=self.env, stdout=log, stderr=subprocess.STDOUT))
            limits = {p.pid: cap[index] if isinstance(cap, tuple) else cap for index, p in enumerate(processes)}
            # Sample only the solver children, not the Python publication threads.
            deadline = time.monotonic() + 90
            while any(p.poll() is None for p in processes):
                self.assertLess(time.monotonic(), deadline, "Concurrent solver timed out")
                if os.name != "nt":
                    # Read the executable by PID: whitespace/Unicode in a
                    # relocated package path cannot be parsed from ps argv.
                    table = subprocess.check_output(["/bin/ps", "-eo", "pid=,ppid=,nlwp="], text=True, env=self.env)
                    parents = {p.pid for p in processes}
                    for line in table.splitlines():
                        fields = line.split()
                        if len(fields) != 3:
                            continue
                        pid, parent, threads = fields
                        if int(parent) not in parents:
                            continue
                        try:
                            executable = Path(f"/proc/{pid}/exe").resolve(strict=True).name
                        except OSError:
                            continue  # Child exited between ps and /proc reads.
                        if executable == "vela_example_runner":
                            samples.append(int(threads))
                            self.assertLessEqual(int(threads), limits[int(parent)], line)
                time.sleep(.01)
            for proc, log in zip(processes, logs):
                proc.wait(timeout=10)
                log.seek(0)
                self.assertEqual(proc.returncode, 0, log.read().decode("utf-8", errors="replace"))
            if os.name != "nt":
                self.assertTrue(samples, "ps did not observe any real solver child")
            return samples
        finally:
            for proc in processes:
                if proc.poll() is None:
                    proc.kill()
                proc.wait(timeout=10)
            for log in logs:
                log.close()

    @unittest.skipIf(os.name == "nt", LINUX_REASON)
    def test_four_cmds_same_directory_match_serial_exactly(self):
        cmds = [synthetic_cmd(self.root, f"iv_{index}") for index in range(4)]
        baselines = [self.baseline(cmd, f"iv_{index}") for index, cmd in enumerate(cmds)]
        samples = self.parallel([(cmd, []) for cmd in cmds])
        for index, cmd in enumerate(cmds):
            self.assertEqual((cmd.parent / f"iv_{index}.plt").read_bytes(), baselines[index][0])
            self.assertEqual(numeric_outputs(cmd.parent, f"iv_{index}"), baselines[index][1])
            self.assertFalse((cmd.parent / (cmd.stem + "_vela/.lock")).exists())
        print(f"WP11 same-directory: four cmd, exact PLT/VTU, solver nlwp max={max(samples)}")

    @unittest.skipIf(os.name == "nt", LINUX_REASON)
    def test_four_directories_match_serial_exactly(self):
        cmds = [synthetic_cmd(self.root / str(index), "pn2d_iv") for index in range(4)]
        baselines = [self.baseline(cmd, "pn2d_iv") for cmd in cmds]
        samples = self.parallel([(cmd, []) for cmd in cmds])
        for cmd, baseline in zip(cmds, baselines):
            self.assertEqual((cmd.parent / "pn2d_iv.plt").read_bytes(), baseline[0])
            self.assertEqual(numeric_outputs(cmd.parent, "pn2d_iv"), baseline[1])
        print(f"WP11 four-directory: exact PLT/VTU, solver nlwp max={max(samples)}")

    @unittest.skipIf(os.name == "nt", LINUX_REASON)
    def test_math_and_cli_thread_precedence_and_process_caps(self):
        commands, limits, baselines = [], [], []
        for math_threads, cli_threads, expected in ((None, None, 1), (2, None, 2), (2, 1, 1), (1, 2, 2)):
            with self.subTest(math=math_threads, cli=cli_threads):
                cmd = synthetic_cmd(self.root / f"{math_threads}_{cli_threads}", "pn2d_iv")
                if math_threads:
                    cmd.write_text(cmd.read_text().replace("Math {", f"Math {{ NumberOfThreads={math_threads}"))
                args = [] if cli_threads is None else ["--threads", str(cli_threads)]
                baselines.append(self.baseline(cmd, "pn2d_iv", *args))
                commands.append((cmd, args))
                limits.append(expected)
        samples = self.parallel(commands, cap=tuple(limits))
        for (cmd, _), expected, baseline in zip(commands, limits, baselines):
            plan = json.loads((cmd.parent / (cmd.stem + "_vela/stage_plan.json")).read_text())
            self.assertEqual(plan["threads"], expected)
            self.assertEqual((cmd.parent / "pn2d_iv.plt").read_bytes(), baseline[0])
            self.assertEqual(numeric_outputs(cmd.parent, "pn2d_iv"), baseline[1])
        print(f"WP11 simultaneous thread combinations: limits={limits}, exact PLT/VTU, nlwp max={max(samples)}")

    def test_duplicate_cmd_second_exit_7_preserves_owner(self):
        # A controlled slow protocol runner makes lock timing deterministic.
        # The numerical matrix above uses the actual C++ runner.
        cmd = synthetic_cmd(self.root)
        slow = self.root / "slow_runner.py"
        fake = REPO / "tests/fixtures/sdevice_fake_runner/fake_runner.py"
        slow.write_text("import time,runpy\ntime.sleep(2)\nrunpy.run_path(" + repr(str(fake)) +
                        ",run_name='__main__')\n", encoding="utf-8")
        env = dict(self.env, VELA_EXAMPLE_RUNNER=str(slow))
        proc = subprocess.Popen(command(LAUNCHER, cmd), cwd=cmd.parent, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8")
        lock = self.root / (cmd.stem + "_vela/.lock")
        try:
            deadline = time.monotonic() + 10
            owner = None
            while time.monotonic() < deadline:
                if proc.poll() is not None:
                    stdout, stderr = proc.communicate()
                    self.fail(stdout + stderr)
                try:
                    candidate = lock.read_bytes()
                    if json.loads(candidate).get("token"):
                        owner = candidate
                        break
                except (FileNotFoundError, ValueError):
                    pass
                time.sleep(.02)
            self.assertIsNotNone(owner, "First process did not publish lock ownership")
            second = run(LAUNCHER, cmd, env=env)
            self.assertEqual(second.returncode, 7, second.stdout + second.stderr)
            self.assertIn("E701", second.stderr)
            self.assertEqual(lock.read_bytes(), owner)
            stdout, stderr = proc.communicate(timeout=30)
            self.assertEqual(proc.returncode, 0, stdout + stderr)
            self.assertFalse(lock.exists())
        finally:
            if proc.poll() is None:
                proc.kill()
            proc.communicate(timeout=10)

if __name__ == "__main__":
    unittest.main()
