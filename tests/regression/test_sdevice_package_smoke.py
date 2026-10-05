"""WP11 offline packager contracts and extracted-package real-runner smoke."""
from __future__ import annotations
import hashlib
import json
import os
import re
from pathlib import Path
from types import SimpleNamespace
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import zipfile

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from package_vela import safe_member, unpack_python, audit, package, tar_permissions, distribution_platform
from build_vela_sdevice_pyz import build_archive
from tests.regression.sdevice_package_support import synthetic_cmd, isolated_environment, run
from run_sdevice_package_e2e import verify_package, runtime_environment
from vela_sdevice.plt import read_plt_document

class OfflinePackageTests(unittest.TestCase):
    def test_native_package_admission_rejects_corruption_and_escape(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            binary = root / 'runner'
            binary.write_bytes(b'qualified runtime')
            item = {'path': 'runner', 'size': binary.stat().st_size,
                    'sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}
            manifest = root / 'PACKAGE_MANIFEST.json'
            manifest.write_text(json.dumps({'files': [item]}))
            verify_package(root)
            binary.write_bytes(b'corrupted runtime')
            with self.assertRaisesRegex(ValueError, 'fingerprint'):
                verify_package(root)
            item['path'] = '../outside'
            manifest.write_text(json.dumps({'files': [item]}))
            with self.assertRaisesRegex(ValueError, 'escapes'):
                verify_package(root)

    def test_native_package_environment_ignores_host_overrides(self):
        with patch.dict(os.environ, {'VELA_PYTHON': 'host-python', 'VELA_EXAMPLE_RUNNER': 'host-runner',
                                     'PYTHONPATH': 'host-modules', 'FAKE_RUNNER_MODE': 'fake'}):
            env = runtime_environment()
        for key in ('VELA_PYTHON', 'VELA_EXAMPLE_RUNNER', 'PYTHONPATH', 'FAKE_RUNNER_MODE'):
            self.assertNotIn(key, env)
        self.assertEqual(env['PYTHONDONTWRITEBYTECODE'], '1')
        self.assertEqual(env['PATH'], str(Path(os.environ.get('SystemRoot', 'C:/Windows')) / 'System32')
                         if os.name == 'nt' else '/nonexistent')

    def test_linux_target_requires_ubuntu_26_04_x86_64(self):
        with patch("package_vela.os.name", "posix"), patch("package_vela.platform.machine", return_value="x86_64"):
            for distro, version, accepted in (("ubuntu", "26.04", True), ("ubuntu", "24.04", False),
                                             ("debian", "26.04", False)):
                with self.subTest(distro=distro, version=version), patch(
                    "package_vela.platform.freedesktop_os_release", return_value={"ID": distro, "VERSION_ID": version}
                ):
                    if accepted:
                        self.assertEqual(distribution_platform(), "ubuntu26.04-x86_64")
                        with patch("package_vela.platform.machine", return_value="aarch64"):
                            with self.assertRaisesRegex(ValueError, "x86_64"):
                                distribution_platform()
                    else:
                        with self.assertRaisesRegex(ValueError, "26.04"):
                            distribution_platform()

    def test_cross_host_frontend_records_explicit_source_commit(self):
        with tempfile.TemporaryDirectory() as temp:
            revision = "0123456789abcdef" * 2 + "01234567"
            archive = build_archive(Path(temp) / "frontend.pyz", revision)
            with zipfile.ZipFile(archive) as zipped:
                info = json.loads(zipped.read("vela_sdevice/_assets/version.json"))
            self.assertEqual(info["commit"], revision)
            with self.assertRaises(ValueError):
                build_archive(Path(temp) / "bad.pyz", "unknown")

    def test_wrong_python_hash_rejected_before_install(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            archive = root / "python.zip"
            archive.write_bytes(b"offline input")
            args = SimpleNamespace(build_dir=root / "missing-build", python_archive=archive,
                                   output_dir=root / "output", python_sha256="0" * 64)
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                package(args)
            self.assertFalse(args.output_dir.exists())

    def test_linux_archive_modes_are_safe_and_executable(self):
        for name, expected in (("vela/bin/vela", 0o755), ("vela/bin/python/bin/python3.12", 0o755),
                               ("vela/share/vela/materials/models.par", 0o644)):
            entry = tarfile.TarInfo(name)
            entry.mode, entry.uid = 0o777, 1000
            normalized = tar_permissions(entry)
            self.assertEqual((normalized.mode, normalized.uid, normalized.gid), (expected, 0, 0))

    def test_paths_reject_archive_escape(self):
        for name in ("/absolute", "../escape", "python/../../escape", "C:/escape", "python\\escape"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                safe_member(name)
        self.assertEqual(safe_member("python/bin/python3").as_posix(), "python/bin/python3")

    def test_python_archive_escape_rejected_before_extract(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for format in ("zip", "tar"):
                archive = root / format
                if format == "zip":
                    with zipfile.ZipFile(archive, "w") as zipped:
                        zipped.writestr("../escape", b"wrong")
                else:
                    with tarfile.open(archive, "w") as tar:
                        info = tarfile.TarInfo("python/bin/python3")
                        info.type, info.linkname = tarfile.SYMTYPE, "../../../escape"
                        tar.addfile(info)
                with self.assertRaises(ValueError):
                    unpack_python(archive, root / ("dest_" + format),
                                  "win64" if format == "zip" else "ubuntu26.04-x86_64")
            self.assertFalse((root / "escape").exists())

    def test_reference_artifacts_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "reference_staging").mkdir()
            (root / "reference_staging/private.tdr").write_bytes(b"reference")
            with self.assertRaises(ValueError):
                audit(root)

class ExtractedPackageSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        archive = os.environ.get("VELA_PACKAGE_ARCHIVE")
        if not archive:
            raise unittest.SkipTest("VELA_PACKAGE_ARCHIVE not set; supply a WP11 offline runtime package")
        cls.archive = Path(archive).resolve()
        if not cls.archive.is_file():
            raise AssertionError(f"Requested package does not exist: {cls.archive}")
        cls.temp = tempfile.TemporaryDirectory(prefix="wp11_", dir=REPO / "build-release" if os.name == "nt" else None)
        cls.root = Path(cls.temp.name) / "space \u4e2d\u6587"
        cls.root.mkdir()
        if cls.archive.suffix == ".zip":
            with zipfile.ZipFile(cls.archive) as archive:
                archive.extractall(cls.root)
        else:
            with tarfile.open(cls.archive) as archive:
                archive.extractall(cls.root, filter="data")
        cls.package = cls.root / "vela"
        cls.launcher = cls.package / "bin" / ("vela.exe" if os.name == "nt" else "vela")
        cls.env = isolated_environment()
        cls.addClassCleanup(cls.temp.cleanup)

    def test_archive_allowlist_hashes_and_notices(self):
        info = json.loads((self.package / "PACKAGE_MANIFEST.json").read_text(encoding="utf-8"))
        self.assertEqual(info["platform"], "win64" if os.name == "nt" else "ubuntu26.04-x86_64")
        self.assertIn("python_sha256", info)
        recorded = {entry["path"]: entry for entry in info["files"]}
        present = {p.relative_to(self.package).as_posix() for p in self.package.rglob("*") if p.is_file()}
        self.assertEqual(present, set(recorded) | {"PACKAGE_MANIFEST.json"})
        for name, entry in recorded.items():
            content = (self.package / name).read_bytes()
            self.assertEqual(hashlib.sha256(content).hexdigest(), entry["sha256"], name)
            self.assertEqual(len(content), entry["size"], name)
            self.assertFalse(set(Path(name).parts) & {"reference_staging", "reference_tcad", ".git"})
            self.assertNotIn(Path(name).suffix, (".tdr", ".tdx", ".plt", ".cmd"))
        with zipfile.ZipFile(self.package / "bin/vela_sdevice.pyz") as zipped:
            self.assertFalse(any("reference_staging" in name or "reference_tcad" in name for name in zipped.namelist()))
        notices = (self.package / "THIRD_PARTY_NOTICES").read_text(encoding="utf-8")
        self.assertIn("Python Software Foundation", notices)
        self.assertIn("HDF", notices)
        self.assertIn("LICENSE.openssl-3.txt", notices)
        if os.name == "nt":
            self.assertEqual((self.package / "bin/python/python312._pth").read_text().splitlines(),
                             ["python312.zip", ".", ".."])

    @unittest.skipIf(os.name == "nt", "ELF/ldd check is Linux only; Windows DLL closure is checked by isolated PATH smoke")
    def test_linux_elf_dependencies_resolve_inside_package(self):
        allowed_system = re.compile(r"^(ld-linux.*|lib(c|m|pthread|dl|rt|resolv)\.so.*)$")
        env = dict(self.env)
        env.pop("LD_LIBRARY_PATH", None)
        for name in ("vela_example_runner", "sentaurus_import"):
            table = subprocess.check_output(["/usr/bin/ldd", str(self.package / "bin" / name)],
                                             env=env, text=True)
            self.assertNotIn("not found", table)
            resolved = re.findall(r"=> (/.+?) \(0x[0-9a-fA-F]+\)", table)
            self.assertTrue(resolved)
            for library in resolved:
                if not allowed_system.fullmatch(Path(library).name):
                    self.assertTrue(Path(library).resolve().is_relative_to(self.package), library)
        info = json.loads((self.package / "PACKAGE_MANIFEST.json").read_text())
        self.assertIn("libstdc++6", info["licenses"])

    def test_no_system_python_unicode_iv_and_parameter_export(self):
        self.assertIsNone(shutil.which("python", path=self.env["PATH"]))
        self.assertIsNone(shutil.which("python3", path=self.env["PATH"]))
        original_files = {p.relative_to(self.package) for p in self.package.rglob("*") if p.is_file()}
        root = self.root / "simulation \u6d4b\u8bd5"
        cmd = synthetic_cmd(root, "iv \u4e2d\u6587")
        # -P must work from the extracted root launcher, then be consumed by File Parameter.
        (root / "models.par").unlink()
        done = subprocess.run([str(self.package / self.launcher.name), "device", "-P"], cwd=root,
            env=self.env, capture_output=True, text=True, encoding="utf-8", timeout=30)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual({p.relative_to(self.package) for p in self.package.rglob("*") if p.is_file()}, original_files,
                         "Running bundled Python must not write shared package caches")
        self.assertEqual((root / "models.par").read_bytes(),
                         (self.package / "share/vela/materials/models.par").read_bytes())
        done = run(self.launcher, cmd, env=self.env)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertNotIn("using Python from PATH", done.stderr)
        names, rows = read_plt_document(next(root.glob("*.plt")).read_text(encoding="utf-8"))
        self.assertEqual(names, ["anode", "cathode"])
        self.assertEqual(rows[-1]["anode"]["outer_voltage_V"], 10.0)
        self.assertTrue(list(root.glob("*.vtu")))
        self.assertTrue(list(root.glob("*.pvd")))
        # Unicode stems have the documented content-hash fallback workspace.
        reports = list(root.glob("*_vela/run_manifest.json"))
        self.assertEqual(len(reports), 1)
        info = json.loads(reports[0].read_text(encoding="utf-8"))
        self.assertTrue(all(stage["goal_reached"] for stage in info["stages"]))
        self.assertFalse(list(root.glob("*_vela/.lock")))
        importer = self.package / "bin" / ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import")
        done = subprocess.run([str(importer), "--help"], env=self.env, capture_output=True,
                              text=True, encoding="utf-8", timeout=30)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual({p.relative_to(self.package) for p in self.package.rglob("*") if p.is_file()}, original_files)

if __name__ == "__main__":
    unittest.main()
