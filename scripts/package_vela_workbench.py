"""Build a portable Windows SDS package from an installed SDevice runtime."""
from __future__ import annotations

import argparse
from collections import deque
from hashlib import file_digest
from importlib import metadata
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

from packaging.requirements import Requirement
from vela_workbench_cases import prepare_cases


SCRIPTS = Path(__file__).resolve().parent
EXCLUDED = shutil.ignore_patterns("__pycache__", "*.pyc")
MKL_COMPONENTS = ("rt", "sequential", "core", "def", "avx2", "avx512", "avx10", "mc3")


def copy_python(target, prefix):
    stdlib = Path(sysconfig.get_path("stdlib"))
    site = Path(sysconfig.get_path("purelib"))
    versioned_lib = target / stdlib.relative_to(prefix)
    shutil.copytree(stdlib, versioned_lib, ignore=shutil.ignore_patterns(
        "__pycache__", "*.pyc", "site-packages", "test", "tests", "ensurepip", "idlelib", "config-*",
    ))
    bundled_site = target / site.relative_to(prefix)
    bundled_site.mkdir()
    packages = set()
    pending = ["matplotlib"]
    while pending:
        package = metadata.distribution(pending.pop())
        name = package.metadata["Name"]
        if name in packages:
            continue
        packages.add(name)
        for entry in package.files:
            path = Path(package.locate_file(entry)).resolve()
            if not path.is_relative_to(site) or not path.is_file() or path.suffix == ".pyc":
                continue
            to = bundled_site / path.relative_to(site)
            to.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, to)
        for declaration in package.requires or ():
            requirement = Requirement(declaration)
            if not requirement.marker or requirement.marker.evaluate({"extra": ""}):
                pending.append(requirement.name)

    binaries = target / "bin"
    binaries.mkdir()
    for name in ("python.exe", "pythonw.exe"):
        shutil.copy2(prefix / "bin" / name, binaries / name)
    for module in sorted((prefix / "lib").iterdir()):
        if module.is_dir() and re.fullmatch(r"(tcl|tk)\d+(\.\d+)?", module.name):
            shutil.copytree(module, target / "lib" / module.name, ignore=EXCLUDED)
    return sorted(packages)


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


def copy_dll_dependencies(package, prefix):
    # Python extension imports need their DLLs beside the interpreter,
    # regardless of the Windows DLL search policy for the solver.
    objdump = prefix / "bin" / "objdump.exe"
    python = package / "python"
    python_bin = python / "bin"
    vela_bin = package / "vela" / "bin"
    extensions = (path for path in python.rglob("*") if path.suffix.lower() in {".exe", ".dll", ".pyd"})
    copy_pe_imports(extensions, python_bin, (prefix / "bin", vela_bin), objdump)
    engines = (path for path in vela_bin.iterdir() if path.suffix.lower() in {".exe", ".dll"})
    copy_pe_imports(engines, vela_bin, (prefix / "bin", python_bin), objdump)


def stage_build_runtime(build, runtime):
    """Stage the SDeviceRuntime files; the packager resolves their PE imports."""
    repo = SCRIPTS.parent
    binaries = runtime / "bin"
    binaries.mkdir(parents=True)
    for name in ("sd-sim", "vela", "vela_example_runner",
                 "vela_example_runner_extended80", "sentaurus_import"):
        shutil.copy2(build / f"{name}.exe", binaries / f"{name}.exe")
    shutil.copy2(binaries / "sd-sim.exe", runtime / "sd-sim.exe")
    shutil.copy2(build / "bin/vela_sdevice.pyz", binaries)
    share = runtime / "share/vela"
    for source in (repo / "schemas", repo / "configs/schema"):
        target = share / "schemas"
        target.mkdir(parents=True, exist_ok=True)
        for path in source.glob("*.json"):
            shutil.copy2(path, target)
    shutil.copytree(repo / "share/vela/materials", share / "materials")
    shutil.copytree(repo / "configs/templates", share / "templates")
    shutil.copy2(repo / "LICENSE", runtime)


def build_package(args, package):
    prefix = Path(sys.base_prefix).resolve()
    runtime = package / "vela"
    if args.build_dir:
        stage_build_runtime(args.build_dir, runtime)
    else:
        shutil.copytree(args.runtime, runtime, ignore=shutil.ignore_patterns("python", "__pycache__"))
    required = [runtime / "sd-sim.exe", runtime / "bin/vela_example_runner.exe",
                runtime / "bin/sentaurus_import.exe", runtime / "bin/vela_sdevice.pyz"]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f"Missing SDevice runtime file: {path}")
    for component in MKL_COMPONENTS:
        name = f"mkl_{component}.3.dll"
        shutil.copy2(args.mkl / "bin" / name, runtime / "bin" / name)
    mkl_licenses = args.mkl / "share/doc/mkl/licensing"
    shutil.copytree(mkl_licenses, package / "licenses/oneMKL")

    packages = copy_python(package / "python", prefix)
    app = package / "app"
    app.mkdir()
    for name in ("vela_workbench.py", "vela_workbench_cases.py", "vela_workbench_results.py",
                 "vela_workbench_text.py"):
        shutil.copy2(SCRIPTS / name, app / name)
    with zipfile.ZipFile(runtime / "bin/vela_sdevice.pyz") as frontend:
        frontend.extractall(app)
    prepare_cases(package / "cases")
    subprocess.run([prefix / "bin/g++.exe", "-std=c++20", "-O2", "-static", "-mwindows",
                    "-municode", SCRIPTS / "vela_workbench_launcher.cpp", "-o", package / "SDS.exe"],
                   check=True)
    copy_dll_dependencies(package, prefix)
    licensed = ["python", "tcl", "tk", *[f"python-{name.lower()}" for name in packages]]
    for name in licensed:
        directory = prefix / "share/licenses" / name
        if directory.is_dir():
            shutil.copytree(directory, package / "licenses" / name)
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
        {"python": sys.version.split()[0], "packages": packages, "mkl": "2026.1",
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
    parser.add_argument("--output", type=Path, required=True, help="New directory for SDS and its ZIP")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"Output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sds-", dir=output.parent) as directory:
        staging = Path(directory) / "release"
        package = staging / "SDS"
        package.mkdir(parents=True)
        build_package(args, package)
        with zipfile.ZipFile(staging / "SDS-windows-x64.zip", "w", zipfile.ZIP_DEFLATED,
                             compresslevel=3) as archive:
            for path in sorted(package.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(staging))
        staging.rename(output)
    print(output / "SDS-windows-x64.zip")


if __name__ == "__main__":
    main()
