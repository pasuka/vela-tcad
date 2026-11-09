"""Contract and integration tests for the pinned Gmsh Python API runtime."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gmsh_runtime import (  # noqa: E402
    EXPECTED_GMSH_VERSION,
    RuntimeSources,
    expected_library_name,
    main,
    stage_runtime,
    wrapper_version,
)
from tests.regression.gmsh_test_runtime import require_gmsh_runtime  # noqa: E402


class GmshRuntimeUnitTest(unittest.TestCase):
    def test_wrapper_version_and_windows_library_contract(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_wrapper_") as tmp:
            wrapper = Path(tmp) / "gmsh.py"
            wrapper.write_text(
                'GMSH_API_VERSION = "4.15.2"\n', encoding="utf-8")
            self.assertEqual(EXPECTED_GMSH_VERSION, wrapper_version(wrapper))
        self.assertEqual(
            "gmsh-4.15.dll",
            expected_library_name(EXPECTED_GMSH_VERSION, "Windows"),
        )

    def test_runtime_staging_uses_wrapper_expected_library_name(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_stage_") as tmp:
            root = Path(tmp)
            wrapper = root / "source_gmsh.py"
            library = root / "libgmsh.dll"
            executable = root / "gmsh.exe"
            wrapper.write_text(
                'GMSH_API_VERSION = "4.15.2"\n', encoding="utf-8")
            library.write_bytes(b"test-library")
            executable.write_bytes(b"test-executable")
            staged = stage_runtime(
                RuntimeSources(wrapper, library, executable, "4.15.2"),
                root / "runtime",
            )
            self.assertEqual(
                expected_library_name(EXPECTED_GMSH_VERSION),
                Path(staged["library"]).name,
            )
            self.assertTrue(Path(staged["wrapper"]).is_file())
            self.assertTrue(Path(staged["library"]).is_file())
            self.assertEqual(64, len(staged["wrapper_sha256"]))
            self.assertEqual(64, len(staged["library_sha256"]))

    def test_missing_explicit_dependency_fails_closed_with_strict_json(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_missing_") as tmp:
            root = Path(tmp)
            output = root / "output"
            output.mkdir()
            (output / "gmsh_probe_manifest.json").write_text(
                "stale\n", encoding="utf-8")
            status = main([
                "--output-dir", str(output),
                "--runtime-dir", str(root / "runtime"),
                "--gmsh-python", str(root / "missing-gmsh.py"),
                "--gmsh-library", str(root / "missing-gmsh.dll"),
                "--gmsh-executable", str(root / "missing-gmsh.exe"),
            ])
            self.assertEqual(2, status)
            report = json.loads(
                (output / "fail_report.json").read_text(encoding="utf-8"),
                parse_constant=lambda value: self.fail(
                    f"non-standard JSON constant {value}"),
            )
            self.assertEqual("failed", report["status"])
            self.assertEqual("gmsh_runtime", report["stage"])
            self.assertIn("does not exist", report["message"])
            self.assertIn("install Gmsh 4.15.2", report["recommendation"])
            self.assertFalse((output / "gmsh_probe_manifest.json").exists())


class GmshRuntimeIntegrationTest(unittest.TestCase):
    def test_packaged_runtime_generates_tri3(self) -> None:
        sources = require_gmsh_runtime()
        wrapper = sources.wrapper
        library = sources.library
        executable = sources.executable

        # Run the public entry point in a child process so a platform can unload
        # the staged shared library before TemporaryDirectory removes it.
        with tempfile.TemporaryDirectory(prefix="vela_gmsh_probe_") as tmp:
            root = Path(tmp)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "probe_gmsh_python.py"),
                    "--output-dir", str(root / "output"),
                    "--runtime-dir", str(root / "runtime"),
                    "--gmsh-python", str(wrapper),
                    "--gmsh-library", str(library),
                    "--gmsh-executable", str(executable),
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(0, completed.returncode, completed.stderr)
            payload = json.loads(
                (root / "output" / "gmsh_probe_manifest.json").read_text(
                    encoding="utf-8"))
            self.assertEqual("passed", payload["status"])
            self.assertEqual(
                {EXPECTED_GMSH_VERSION}, set(payload["versions"].values()))
            self.assertEqual(
                EXPECTED_GMSH_VERSION,
                payload["smoke_test"]["shared_library_version"],
            )
            self.assertGreater(payload["smoke_test"]["node_count"], 0)
            self.assertGreater(payload["smoke_test"]["triangle_count"], 0)
            self.assertEqual([2], payload["smoke_test"]["element_types"])
            self.assertEqual(
                Path(payload["staging"]["library"]).resolve(),
                Path(payload["staging"]["loaded_library"]).resolve(),
            )
            manifest = json.loads(
                (root / "output" / "gmsh_probe_manifest.json").read_text(
                    encoding="utf-8"),
                parse_constant=lambda value: self.fail(
                    f"non-standard JSON constant {value}"),
            )
            self.assertEqual(payload, manifest)
            self.assertFalse((root / "output" / "fail_report.json").exists())

            mismatch = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "probe_gmsh_python.py"),
                    "--output-dir", str(root / "mismatch-output"),
                    "--runtime-dir", str(root / "unused-runtime"),
                    "--gmsh-python", str(wrapper),
                    "--gmsh-library", str(library),
                    "--gmsh-executable", str(executable),
                    "--expected-version", "4.15.1",
                ],
                cwd=REPO,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(2, mismatch.returncode)
            mismatch_report = json.loads(
                (root / "mismatch-output" / "fail_report.json").read_text(
                    encoding="utf-8"))
            self.assertIn("versions do not match", mismatch_report["message"])
            self.assertEqual(
                "4.15.1", mismatch_report["details"]["versions"]["expected"])


if __name__ == "__main__":
    unittest.main()
