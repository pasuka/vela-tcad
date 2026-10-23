"""WP7 regression tests: native `vela` launcher (contract C9).

Runs the built launcher against a fake ``vela_sdevice.pyz`` that echoes its
arguments, working directory and environment as JSON and exits with the code
given by ``FAKE_EXIT``.  Skips when the launcher has not been built.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipapp
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EXE_SUFFIX = ".exe" if os.name == "nt" else ""
LAUNCHER = Path(os.environ.get("VELA_LAUNCHER") or (REPO / "build" / f"vela{EXE_SUFFIX}"))

FAKE_MAIN = """\
import json
import os
import sys

print(json.dumps({
    "argv": sys.argv[1:],
    "cwd": os.getcwd(),
    "executable": sys.executable,
    "path": os.environ.get("PATH", ""),
    "nousersite": os.environ.get("PYTHONNOUSERSITE"),
}))
sys.stdout.flush()
sys.exit(int(os.environ.get("FAKE_EXIT", "0")))
"""

# Spaces, double quotes, a backslash before a quote, a trailing backslash, an
# empty argument, a tab and Chinese text (escaped: tests/ stays ASCII-only).
TRICKY_ARGS = [
    "plain",
    "with space",
    'dq"inside',
    'back\\"slash-quote',
    "ends with backslash\\",
    "",
    "tab\there",
    "\u4e2d\u6587\u8def\u5f84 \u53c2\u6570",
    "-P",
    "--threads",
    "4",
    "\u6d4b\u8bd5_des.cmd",
]


def same_path(left: str | os.PathLike[str], right: str | os.PathLike[str]) -> bool:
    return os.path.normcase(os.path.realpath(left)) == os.path.normcase(os.path.realpath(right))


class VelaLauncherTest(unittest.TestCase):
    tmp: Path
    pyz: Path
    root_layout: Path  # <root>/vela + <root>/bin/vela_sdevice.pyz
    flat_layout: Path  # <dir>/vela + <dir>/vela_sdevice.pyz
    bare_layout: Path  # <dir>/vela only
    empty_dir: Path

    @classmethod
    def setUpClass(cls) -> None:
        if not LAUNCHER.is_file():
            raise unittest.SkipTest(f"launcher not built: {LAUNCHER}")
        cls.tmp = Path(tempfile.mkdtemp(prefix="vela_launcher_"))
        source = cls.tmp / "fake_src"
        source.mkdir()
        (source / "__main__.py").write_text(FAKE_MAIN, encoding="ascii")
        cls.pyz = cls.tmp / "fake_vela_sdevice.pyz"
        zipapp.create_archive(source, cls.pyz)

        cls.root_layout = cls.tmp / "root"
        (cls.root_layout / "bin").mkdir(parents=True)
        shutil.copy2(LAUNCHER, cls.root_layout / LAUNCHER.name)
        shutil.copy2(cls.pyz, cls.root_layout / "bin" / "vela_sdevice.pyz")

        cls.flat_layout = cls.tmp / "flat"
        cls.flat_layout.mkdir()
        shutil.copy2(LAUNCHER, cls.flat_layout / LAUNCHER.name)
        shutil.copy2(cls.pyz, cls.flat_layout / "vela_sdevice.pyz")

        cls.bare_layout = cls.tmp / "bare"
        cls.bare_layout.mkdir()
        shutil.copy2(LAUNCHER, cls.bare_layout / LAUNCHER.name)

        cls.empty_dir = cls.tmp / "empty"
        cls.empty_dir.mkdir()

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # -- helpers -----------------------------------------------------------

    def base_env(self, *, fake_exit: int = 0, use_override: bool = True) -> dict[str, str]:
        env = os.environ.copy()
        env.pop("VELA_SDEVICE_PYZ", None)
        env.pop("VELA_PYTHON", None)
        if use_override:
            env["VELA_SDEVICE_PYZ"] = str(self.pyz)
        env["VELA_PYTHON"] = sys.executable
        env["PYTHONIOENCODING"] = "utf-8"
        env["PYTHONUTF8"] = "1"
        env["FAKE_EXIT"] = str(fake_exit)
        return env

    def run_launcher(self, args: list[str], *, launcher: Path = LAUNCHER, env: dict[str, str] | None = None,
                     cwd: Path | None = None) -> subprocess.CompletedProcess[bytes]:
        return subprocess.run([str(launcher), *args], env=env if env is not None else self.base_env(),
                              cwd=str(cwd) if cwd is not None else None, capture_output=True, check=False)

    @staticmethod
    def echo(result: subprocess.CompletedProcess[bytes]) -> dict:
        return json.loads(result.stdout.decode("utf-8"))

    # -- (a) argument pass-through, VELA_PYTHON override ---------------------

    def test_arguments_pass_through_unchanged(self) -> None:
        result = self.run_launcher(["device", *TRICKY_ARGS])
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", "replace"))
        payload = self.echo(result)
        self.assertEqual(payload["argv"], ["device", *TRICKY_ARGS])
        self.assertTrue(same_path(payload["executable"], sys.executable), payload["executable"])
        self.assertEqual(payload["nousersite"], "1")
        self.assertNotIn(b"using Python from PATH", result.stderr)

    def test_sd_sim_direct_arguments_and_exit_code(self) -> None:
        direct = Path(os.environ.get("SD_SIM_LAUNCHER") or LAUNCHER.with_name(f"sd-sim{EXE_SUFFIX}"))
        if not direct.is_file():
            self.skipTest("sd-sim launcher not built")
        for layout in (self.root_layout, self.flat_layout):
            target = layout / direct.name
            shutil.copy2(direct, target)
            result = self.run_launcher(TRICKY_ARGS, launcher=target,
                                       env=self.base_env(fake_exit=2, use_override=False))
            self.assertEqual(result.returncode, 2, result.stderr.decode("utf-8", "replace"))
            self.assertEqual(self.echo(result)["argv"], ["sd-sim", *TRICKY_ARGS])

    # -- (b) exit code propagation --------------------------------------------

    def test_exit_codes_are_propagated(self) -> None:
        for code in range(9):
            with self.subTest(code=code):
                result = self.run_launcher(["device", "x"], env=self.base_env(fake_exit=code))
                self.assertEqual(result.returncode, code, result.stderr.decode("utf-8", "replace"))
                self.assertEqual(self.echo(result)["argv"], ["device", "x"])

    # -- (c) no interpreter ----------------------------------------------------

    def test_no_python_interpreter_exits_2_with_e202(self) -> None:
        env = self.base_env()
        env["VELA_PYTHON"] = str(self.tmp / "missing" / f"python{EXE_SUFFIX}")
        env["PATH"] = str(self.empty_dir)
        result = self.run_launcher(["device", "x"], launcher=self.bare_layout / LAUNCHER.name, env=env)
        self.assertEqual(result.returncode, 2)
        stderr = result.stderr.decode("utf-8", "replace")
        self.assertIn("E202", stderr)
        self.assertIn("no Python interpreter available", stderr)
        self.assertEqual(result.stdout, b"")

    # -- (d) usage errors -------------------------------------------------------

    def test_no_subcommand_prints_usage_and_exits_2(self) -> None:
        result = self.run_launcher([])
        self.assertEqual(result.returncode, 2)
        stderr = result.stderr.decode("utf-8", "replace")
        self.assertIn("E201", stderr)
        self.assertIn("Usage: vela", stderr)

    def test_unknown_subcommand_exits_2(self) -> None:
        result = self.run_launcher(["frobnicate", "x"])
        self.assertEqual(result.returncode, 2)
        stderr = result.stderr.decode("utf-8", "replace")
        self.assertIn("frobnicate", stderr)
        self.assertIn("Usage: vela", stderr)

    # -- (e) version -------------------------------------------------------------

    def test_version_flags(self) -> None:
        for flag in ("--version", "-v"):
            with self.subTest(flag=flag):
                result = self.run_launcher([flag])
                self.assertEqual(result.returncode, 0)
                lines = result.stdout.decode("utf-8").splitlines()
                self.assertEqual(len(lines), 1, lines)
                self.assertTrue(lines[0].startswith("vela "), lines[0])
                self.assertIn("sdevice-syntax T-2022.03", lines[0])

    # -- (f) working directory ---------------------------------------------------

    def test_child_inherits_working_directory(self) -> None:
        workdir = self.tmp / "work dir \u5de5\u4f5c"
        workdir.mkdir()
        result = self.run_launcher(["device", "case_des.cmd"], cwd=workdir)
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", "replace"))
        self.assertTrue(same_path(self.echo(result)["cwd"], workdir))

    # -- bin/ discovery and child environment ------------------------------------

    def test_root_bin_layout_is_discovered_and_bin_prepended_to_path(self) -> None:
        env = self.base_env(use_override=False)
        launcher = self.root_layout / LAUNCHER.name
        result = self.run_launcher(["device", "a", "b c"], launcher=launcher, env=env)
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", "replace"))
        payload = self.echo(result)
        self.assertEqual(payload["argv"], ["device", "a", "b c"])
        self.assertTrue(same_path(payload["path"].split(os.pathsep)[0], self.root_layout / "bin"))
        self.assertEqual(payload["nousersite"], "1")

    def test_flat_layout_is_discovered(self) -> None:
        env = self.base_env(use_override=False)
        launcher = self.flat_layout / LAUNCHER.name
        result = self.run_launcher(["device"], launcher=launcher, env=env)
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", "replace"))
        payload = self.echo(result)
        self.assertEqual(payload["argv"], ["device"])
        self.assertTrue(same_path(payload["path"].split(os.pathsep)[0], self.flat_layout))

    def test_missing_pyz_exits_2(self) -> None:
        env = self.base_env(use_override=False)
        result = self.run_launcher(["device", "x"], launcher=self.bare_layout / LAUNCHER.name, env=env)
        self.assertEqual(result.returncode, 2)
        stderr = result.stderr.decode("utf-8", "replace")
        self.assertIn("E202", stderr)
        self.assertIn("vela_sdevice.pyz", stderr)

    def test_path_fallback_warns(self) -> None:
        python_dir = Path(sys.executable).resolve().parent
        names = [f"python3{EXE_SUFFIX}", f"python{EXE_SUFFIX}"]
        if not any((python_dir / name).is_file() for name in names):
            self.skipTest(f"no python3/python next to {sys.executable}")
        env = self.base_env()
        del env["VELA_PYTHON"]
        env["PATH"] = str(python_dir)
        result = self.run_launcher(["device", "x"], launcher=self.bare_layout / LAUNCHER.name, env=env)
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", "replace"))
        self.assertIn(b"vela: using Python from PATH", result.stderr)
        self.assertEqual(self.echo(result)["argv"], ["device", "x"])


if __name__ == "__main__":
    unittest.main()
