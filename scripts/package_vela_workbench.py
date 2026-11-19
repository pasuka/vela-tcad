"""Build a portable Windows SDS package from an installed SDevice runtime."""
from __future__ import annotations

import argparse
from collections import deque
from hashlib import file_digest
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import zipfile

from vela_workbench_cases import prepare_cases


SCRIPTS = Path(__file__).resolve().parent
SOURCE_FILES = {".py", ".pyc", ".pyo", ".pyw", ".pyz"}
ARCHIVE_FILES = {".zip", ".whl", ".egg"}
MKL_COMPONENTS = ("rt", "sequential", "core", "def", "avx2", "avx512", "avx10", "mc3")


def compile_frontends(python, source_archive, package, scratch):
    compile_dir = scratch / "nuitka"
    compile_dir.mkdir()
    entry = compile_dir / "vela_sdevice_entry.py"
    entry.write_text("from vela_sdevice.__main__ import main\nraise SystemExit(main())\n", encoding="utf-8")
    with zipfile.ZipFile(source_archive) as archive:
        for name in archive.namelist():
            if name.startswith("vela_sdevice/_assets/"):
                archive.extract(name, compile_dir)

    env = dict(os.environ, PYTHONPATH=str(SCRIPTS))
    env.pop("PYTHONHOME", None)
    common = [python, "-m", "nuitka", "--mode=standalone", "--zig",
              "--assume-yes-for-downloads", f"--jobs={min(os.cpu_count() or 1, 4)}",
              f"--output-dir={compile_dir}"]
    builds = (
        (SCRIPTS / "vela_workbench.py", "SDS-GUI.exe", package / "app", [
            "--enable-plugin=tk-inter",
            "--nofollow-import-to=vela_workbench_catalog", "--windows-console-mode=disable",
        ]),
        (entry, "vela-sdevice.exe", package / "vela/bin/frontend", [
            f"--include-data-dir={compile_dir / 'vela_sdevice/_assets'}=vela_sdevice/_assets",
        ]),
    )
    for source, name, target, options in builds:
        subprocess.run([*common, f"--output-filename={name}", *options, source],
                       cwd=compile_dir, env=env, check=True)
        distribution = compile_dir / (source.stem + ".dist")
        if not (distribution / name).is_file():
            raise FileNotFoundError(distribution / name)
        shutil.copytree(distribution, target)


def copy_python_licenses(python, destination):
    paths = subprocess.check_output([
        python, "-c", "import json, sys, sysconfig; "
        "print(json.dumps([sys.base_prefix, sysconfig.get_path('purelib')]))",
    ], text=True)
    prefix, site = map(Path, json.loads(paths))
    shutil.copy2(prefix / "LICENSE.txt", destination / "CPython.txt")
    for metadata in site.glob("*.dist-info"):
        target = destination / metadata.stem
        for path in metadata.rglob("*"):
            name = path.name.lower()
            if path.is_file() and (name.startswith(("license", "copying", "notice"))
                                   or "licenses" in path.relative_to(metadata).parts):
                output = target / path.relative_to(metadata)
                output.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, output)


def check_archive(archive, prefix=""):
    for entry in archive.infolist():
        if entry.is_dir():
            continue
        suffix = Path(entry.filename).suffix.lower()
        if suffix in SOURCE_FILES:
            raise ValueError(f"Python code in release archive: {prefix}{entry.filename}")
        with archive.open(entry) as stream:
            signature = stream.read(4)
        if suffix in ARCHIVE_FILES or signature in (b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08"):
            with zipfile.ZipFile(io.BytesIO(archive.read(entry))) as nested:
                check_archive(nested, f"{prefix}{entry.filename}!")


def copy_pe_imports(paths, destination, search, objdump):
    win = Path(os.environ["WINDIR"]) / "System32"
    binaries = deque(paths)
    bundled = {path.name.lower() for path in destination.iterdir()}
    while binaries:
        batch = [binaries.popleft() for _ in range(min(24, len(binaries)))]
        dump = subprocess.check_output([objdump, "-p", *batch], text=True, errors="replace")
        for name in set(re.findall(r"^\s*DLL Name:\s*(\S+)", dump, re.MULTILINE)):
            if name.lower() in bundled or name.lower().startswith(("api-ms-", "ext-ms-")):
                continue
            source = next((folder / name for folder in search if (folder / name).is_file()), None)
            if source is None and (win / name).is_file():
                if not name.lower().startswith(("vcruntime", "msvcp", "concrt", "vcomp")):
                    continue
                source = win / name
            if source is None:
                raise FileNotFoundError(f"Missing dependent DLL: {name}")
            output = destination / name
            shutil.copy2(source, output)
            bundled.add(name.lower())
            binaries.append(output)


def copy_dll_dependencies(runtime, prefix):
    binaries = runtime / "bin"
    engines = (path for path in binaries.iterdir() if path.suffix.lower() in {".exe", ".dll"})
    copy_pe_imports(engines, binaries, (prefix / "bin",), prefix / "bin/objdump.exe")


def stage_build_runtime(build, runtime):
    """Stage the SDeviceRuntime files; the packager resolves their PE imports."""
    repo = SCRIPTS.parent
    binaries = runtime / "bin"
    binaries.mkdir(parents=True)
    for name in ("sd-sim", "vela", "vela_example_runner",
                 "vela_example_runner_extended80", "sentaurus_import"):
        shutil.copy2(build / f"{name}.exe", binaries / f"{name}.exe")
    shutil.copy2(binaries / "sd-sim.exe", runtime / "sd-sim.exe")
    share = runtime / "share/vela"
    for source in (repo / "schemas", repo / "configs/schema"):
        target = share / "schemas"
        target.mkdir(parents=True, exist_ok=True)
        for path in source.glob("*.json"):
            shutil.copy2(path, target)
    shutil.copytree(repo / "share/vela/materials", share / "materials")
    shutil.copytree(repo / "configs/templates", share / "templates")
    shutil.copy2(repo / "LICENSE", runtime)


def build_package(args, package, scratch):
    prefix = Path(sys.base_prefix).resolve()
    runtime = package / "vela"
    if args.build_dir:
        stage_build_runtime(args.build_dir, runtime)
        source_archive = args.build_dir / "bin/vela_sdevice.pyz"
    else:
        shutil.copytree(args.runtime, runtime, ignore=shutil.ignore_patterns(
            "python", "__pycache__", "vela_sdevice.pyz",
        ))
        source_archive = args.runtime / "bin/vela_sdevice.pyz"
    for path in (runtime / "sd-sim.exe", runtime / "bin/vela_example_runner.exe",
                 runtime / "bin/sentaurus_import.exe", source_archive):
        if not path.is_file():
            raise FileNotFoundError(f"Missing SDevice runtime file: {path}")
    for component in MKL_COMPONENTS:
        name = f"mkl_{component}.3.dll"
        shutil.copy2(args.mkl / "bin" / name, runtime / "bin" / name)
    mkl_licenses = args.mkl / "share/doc/mkl/licensing"
    shutil.copytree(mkl_licenses, package / "licenses/oneMKL")

    compile_frontends(args.nuitka_python, source_archive, package, scratch)
    prepare_cases(package / "cases")
    subprocess.run([prefix / "bin/g++.exe", "-std=c++20", "-O2", "-static", "-mwindows",
                    "-municode", SCRIPTS / "vela_workbench_launcher.cpp", "-o", package / "SDS.exe"],
                   check=True)
    copy_dll_dependencies(runtime, prefix)
    copy_python_licenses(args.nuitka_python, package / "licenses")
    env = dict(os.environ)
    # Verify that a runtime supplied via --runtime uses the compiled frontend.
    env.pop("VELA_SDEVICE_PYZ", None)
    subprocess.run([runtime / "sd-sim.exe", "--version"], cwd=package, env=env, check=True)
    (package / "README.txt").write_text(
        "Semiconductor Device Simulator\n\n"
        "Extract the entire ZIP, then double-click SDS.exe. The application uses the bundled\n"
        "reference inputs and writes runs under %LOCALAPPDATA%\\VelaTCAD\\runs.\n"
        "Click the Chinese-language switch in the toolbar to switch languages.\n",
        encoding="utf-8",
    )
    with (runtime / "sd-sim.exe").open("rb") as stream:
        solver_digest = file_digest(stream, "sha256").hexdigest()
    (package / "PACKAGE.json").write_text(json.dumps(
        {"frontend": "Nuitka standalone", "mkl": "2026.1",
         "sd_sim_sha256": solver_digest,
         "case_count": len(json.loads((package / "cases/cases.json").read_text(encoding="utf-8")))},
        indent=2,
    ) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--runtime", type=Path, help="SDeviceRuntime install directory")
    source.add_argument("--build-dir", type=Path, help="Completed MSYS2 UCRT64 Release build directory")
    parser.add_argument("--mkl", type=Path, required=True, help="oneMKL 2026.1 installation root")
    parser.add_argument("--nuitka-python", type=Path, required=True,
                        help="Windows x64 CPython 3.12 with Nuitka and Matplotlib installed")
    parser.add_argument("--output", type=Path, required=True, help="New directory for SDS and its ZIP")
    args = parser.parse_args()
    for name in ("build_dir", "runtime", "mkl", "nuitka_python"):
        value = getattr(args, name)
        if value:
            setattr(args, name, value.resolve())
    output = args.output.resolve()
    if output.exists():
        parser.error(f"Output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sds-", dir=output.parent) as directory:
        staging = Path(directory) / "release"
        package = staging / "SDS"
        package.mkdir(parents=True)
        build_package(args, package, Path(directory))
        with zipfile.ZipFile(staging / "SDS-windows-x64.zip", "w", zipfile.ZIP_DEFLATED,
                             compresslevel=3) as archive:
            for path in sorted(package.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(staging))
        with zipfile.ZipFile(staging / "SDS-windows-x64.zip") as archive:
            check_archive(archive)
        staging.rename(output)
    print(output / "SDS-windows-x64.zip")


if __name__ == "__main__":
    main()
