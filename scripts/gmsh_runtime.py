#!/usr/bin/env python3
"""Discover, stage, and smoke-test the pinned Gmsh Python API runtime."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


EXPECTED_GMSH_VERSION = "4.15.2"
PROBE_SCHEMA = "vela.gmsh_runtime_probe.v1"
FAIL_SCHEMA = "vela.gmsh_runtime_probe.fail_report.v1"
REPO = Path(__file__).resolve().parents[1]


class GmshRuntimeError(RuntimeError):
    """A fail-closed runtime discovery or probe error."""

    def __init__(self, message: str, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.details = details or {}


@dataclass(frozen=True)
class RuntimeSources:
    wrapper: Path
    library: Path
    executable: Path
    version: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def wrapper_version(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        r'^GMSH_API_VERSION\s*=\s*["\'](?P<version>\d+\.\d+\.\d+)["\']',
        text,
        re.MULTILINE,
    )
    if match is None:
        raise GmshRuntimeError(
            f"cannot read GMSH_API_VERSION from Python wrapper: {path}",
            {"wrapper": str(path)},
        )
    return match.group("version")


def expected_library_name(version: str, system: str | None = None) -> str:
    major, minor, _patch = version.split(".")
    current = system or platform.system()
    if current == "Windows":
        return f"gmsh-{major}.{minor}.dll"
    if current == "Darwin":
        return f"libgmsh.{major}.{minor}.dylib"
    return f"libgmsh.so.{major}.{minor}"


def _explicit_file(path: Path | None, label: str) -> Path | None:
    if path is None:
        return None
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise GmshRuntimeError(
            f"explicit {label} does not exist: {resolved}",
            {label: str(resolved)},
        )
    return resolved


def _first_file(candidates: Sequence[Path]) -> Path | None:
    seen: set[Path] = set()
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_file():
            return resolved
    return None


def discover_runtime_sources(
    gmsh_python: Path | None = None,
    gmsh_library: Path | None = None,
    gmsh_executable: Path | None = None,
    environ: dict[str, str] | None = None,
) -> RuntimeSources:
    """Resolve a wrapper, shared library, and CLI without changing the system."""
    env = os.environ if environ is None else environ
    wrapper = _explicit_file(
        gmsh_python or (Path(env["VELA_GMSH_PYTHON"])
                        if env.get("VELA_GMSH_PYTHON") else None),
        "gmsh_python",
    )
    executable = _explicit_file(
        gmsh_executable or (Path(env["VELA_GMSH_EXECUTABLE"])
                            if env.get("VELA_GMSH_EXECUTABLE") else None),
        "gmsh_executable",
    )
    if executable is None:
        found = shutil.which("gmsh") or shutil.which("gmsh.exe")
        if found:
            executable = Path(found).resolve()
        else:
            python_prefixes = [
                Path(sys.prefix), Path(sys.executable).resolve().parent.parent,
            ]
            executable = _first_file([
                prefix / "bin" / name
                for prefix in python_prefixes
                for name in ("gmsh.exe", "gmsh")
            ])

    prefixes = [Path(sys.prefix), Path(sys.executable).resolve().parent.parent]
    if executable is not None:
        prefixes.append(executable.parent.parent)
    if wrapper is None:
        wrapper = _first_file([
            *(prefix / "lib" / "gmsh.py" for prefix in prefixes),
            *(prefix / "Lib" / "site-packages" / "gmsh.py" for prefix in prefixes),
        ])

    library = _explicit_file(
        gmsh_library or (Path(env["VELA_GMSH_LIBRARY"])
                         if env.get("VELA_GMSH_LIBRARY") else None),
        "gmsh_library",
    )
    version = wrapper_version(wrapper) if wrapper is not None else "unknown"
    if library is None:
        names = [
            expected_library_name(version) if version != "unknown" else "",
            "libgmsh.dll",
            "gmsh.dll",
            "libgmsh.so",
            "libgmsh.dylib",
        ]
        directories = [
            *(prefix / "bin" for prefix in prefixes),
            *(prefix / "lib" for prefix in prefixes),
        ]
        library = _first_file([
            directory / name
            for directory in directories
            for name in names
            if name
        ])

    missing = [
        name for name, value in (
            ("gmsh Python wrapper", wrapper),
            ("Gmsh shared library", library),
            ("gmsh executable", executable),
        ) if value is None
    ]
    if missing:
        raise GmshRuntimeError(
            "missing Gmsh runtime component(s): " + ", ".join(missing),
            {
                "missing": missing,
                "environment_variables": [
                    "VELA_GMSH_PYTHON",
                    "VELA_GMSH_LIBRARY",
                    "VELA_GMSH_EXECUTABLE",
                ],
                "recommendation": (
                    "install Gmsh 4.15.2 with its Python wrapper and shared "
                    "library, or pass all three explicit paths"
                ),
            },
        )
    assert wrapper is not None and library is not None and executable is not None
    return RuntimeSources(wrapper, library, executable, version)


def executable_version(path: Path) -> str:
    completed = subprocess.run(
        [str(path), "-version"],
        text=True,
        capture_output=True,
        check=False,
    )
    combined = "\n".join([completed.stdout, completed.stderr]).strip()
    match = re.search(r"\b\d+\.\d+\.\d+\b", combined)
    if completed.returncode != 0 or match is None:
        raise GmshRuntimeError(
            f"unable to query gmsh executable version: {path}",
            {
                "executable": str(path),
                "returncode": completed.returncode,
                "stdout": completed.stdout,
                "stderr": completed.stderr,
            },
        )
    return match.group(0)


def stage_runtime(sources: RuntimeSources, runtime_dir: Path) -> dict[str, Any]:
    """Stage an importable wrapper/DLL layout under an ignored build path."""
    runtime = runtime_dir.resolve()
    runtime.mkdir(parents=True, exist_ok=True)
    wrapper_target = runtime / "gmsh.py"
    library_target = runtime / expected_library_name(sources.version)
    shutil.copy2(sources.wrapper, wrapper_target)
    shutil.copy2(sources.library, library_target)
    return {
        "directory": str(runtime),
        "wrapper": str(wrapper_target),
        "wrapper_sha256": sha256(wrapper_target),
        "library": str(library_target),
        "library_sha256": sha256(library_target),
    }


def _load_staged_module(wrapper: Path):
    importlib.invalidate_caches()
    sys.modules.pop("gmsh", None)
    spec = importlib.util.spec_from_file_location("gmsh", wrapper)
    if spec is None or spec.loader is None:
        raise GmshRuntimeError(
            f"unable to create an import spec for staged wrapper: {wrapper}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["gmsh"] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop("gmsh", None)
        raise
    return module


def _mesh_smoke_test(gmsh) -> dict[str, Any]:
    initialized = False
    try:
        gmsh.initialize(["vela-gmsh-probe", "-nopopup"], False)
        initialized = True
        library_version = gmsh.option.getString("General.Version")
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.RandomFactor", 0)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.option.setNumber("Mesh.Optimize", 0)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
        gmsh.model.add("vela_gmsh_runtime_probe")
        points = [
            gmsh.model.occ.addPoint(0.0, 0.0, 0.0, 2.0),
            gmsh.model.occ.addPoint(1.0, 0.0, 0.0, 2.0),
            gmsh.model.occ.addPoint(0.0, 1.0, 0.0, 2.0),
        ]
        lines = [
            gmsh.model.occ.addLine(points[0], points[1]),
            gmsh.model.occ.addLine(points[1], points[2]),
            gmsh.model.occ.addLine(points[2], points[0]),
        ]
        loop = gmsh.model.occ.addCurveLoop(lines)
        surface = gmsh.model.occ.addPlaneSurface([loop])
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(2, [surface], 1001, "region:probe")
        gmsh.model.mesh.generate(2)
        node_tags, coordinates, _parameters = gmsh.model.mesh.getNodes()
        element_types, element_tags, element_nodes = gmsh.model.mesh.getElements(2)
        if len(node_tags) == 0 or len(element_types) == 0:
            raise GmshRuntimeError("Gmsh generated an empty 2-D probe mesh")
        if any(int(kind) != 2 for kind in element_types):
            raise GmshRuntimeError(
                "Gmsh probe generated a non-Tri3 2-D element",
                {"element_types": [int(kind) for kind in element_types]},
            )
        count = sum(len(tags) for tags in element_tags)
        connectivity_count = sum(len(nodes) for nodes in element_nodes)
        if connectivity_count != 3 * count:
            raise GmshRuntimeError("Gmsh Tri3 connectivity length is inconsistent")
        return {
            "model": "vela_gmsh_runtime_probe",
            "shared_library_version": library_version,
            "node_count": len(node_tags),
            "coordinate_count": len(coordinates),
            "triangle_count": count,
            "element_types": [int(kind) for kind in element_types],
            "node_tags_positive": all(int(tag) > 0 for tag in node_tags),
        }
    finally:
        if initialized:
            gmsh.finalize()


def run_probe(
    output_dir: Path,
    runtime_dir: Path,
    gmsh_python: Path | None = None,
    gmsh_library: Path | None = None,
    gmsh_executable: Path | None = None,
    expected_version: str = EXPECTED_GMSH_VERSION,
    environ: dict[str, str] | None = None,
) -> dict[str, Any]:
    output = output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    for stale in ("fail_report.json", "gmsh_probe_manifest.json"):
        path = output / stale
        if path.exists():
            path.unlink()
    sources = discover_runtime_sources(
        gmsh_python, gmsh_library, gmsh_executable, environ)
    cli_version = executable_version(sources.executable)
    versions = {
        "expected": expected_version,
        "wrapper": sources.version,
        "executable": cli_version,
    }
    if len(set(versions.values())) != 1:
        raise GmshRuntimeError(
            "Gmsh wrapper, executable, and pinned versions do not match",
            {"versions": versions},
        )
    staged = stage_runtime(sources, runtime_dir)
    try:
        gmsh = _load_staged_module(Path(staged["wrapper"]))
    except Exception as error:  # noqa: BLE001 - normalized below
        raise GmshRuntimeError(
            f"failed to import staged Gmsh Python API: {error}",
            {"staging": staged, "exception_type": type(error).__name__},
        ) from error
    api_version = str(getattr(gmsh, "__version__", "unknown"))
    loaded_library = Path(str(getattr(gmsh, "libpath", ""))).resolve()
    if api_version != expected_version:
        raise GmshRuntimeError(
            "loaded Gmsh Python API version does not match the pin",
            {"versions": {**versions, "api": api_version}},
        )
    if loaded_library != Path(staged["library"]).resolve():
        raise GmshRuntimeError(
            "Gmsh wrapper loaded an undeclared shared library",
            {
                "loaded_library": str(loaded_library),
                "expected_staged_library": staged["library"],
            },
        )
    smoke = _mesh_smoke_test(gmsh)
    library_version = smoke["shared_library_version"]
    if library_version != expected_version:
        raise GmshRuntimeError(
            "loaded Gmsh shared library version does not match the pin",
            {
                "versions": {
                    **versions,
                    "api": api_version,
                    "library": library_version,
                }
            },
        )
    payload = {
        "schema": PROBE_SCHEMA,
        "status": "passed",
        "versions": {
            **versions,
            "api": api_version,
            "library": library_version,
        },
        "python": {
            "executable": str(Path(sys.executable).resolve()),
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "sources": {
            "wrapper": str(sources.wrapper),
            "wrapper_sha256": sha256(sources.wrapper),
            "library": str(sources.library),
            "library_sha256": sha256(sources.library),
            "executable": str(sources.executable),
            "executable_sha256": sha256(sources.executable),
        },
        "staging": {**staged, "loaded_library": str(loaded_library)},
        "options": {
            "General.NumThreads": 1,
            "Mesh.Algorithm": 6,
            "Mesh.RandomFactor": 0,
            "Mesh.RecombineAll": 0,
            "Mesh.ElementOrder": 1,
            "Mesh.MshFileVersion": 4.1,
            "Mesh.Binary": 0,
            "Mesh.Optimize": 0,
            "Mesh.OptimizeNetgen": 0,
        },
        "smoke_test": smoke,
    }
    _atomic_json(output / "gmsh_probe_manifest.json", payload)
    return payload


def write_fail_report(
    output_dir: Path,
    error: BaseException,
    expected_version: str = EXPECTED_GMSH_VERSION,
) -> Path:
    details = error.details if isinstance(error, GmshRuntimeError) else {}
    payload = {
        "schema": FAIL_SCHEMA,
        "status": "failed",
        "stage": "gmsh_runtime",
        "message": str(error),
        "expected_version": expected_version,
        "exception_type": type(error).__name__,
        "recommendation": (
            "install Gmsh 4.15.2 with its Python wrapper, shared library and "
            "CLI, or configure VELA_GMSH_PYTHON, VELA_GMSH_LIBRARY and "
            "VELA_GMSH_EXECUTABLE"
        ),
        "details": details,
    }
    path = output_dir.resolve() / "fail_report.json"
    _atomic_json(path, payload)
    manifest = output_dir.resolve() / "gmsh_probe_manifest.json"
    if manifest.exists():
        manifest.unlink()
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--runtime-dir", type=Path,
        default=REPO / "build" / "gmsh-python-runtime",
        help="Ignored directory used to stage gmsh.py and its shared library",
    )
    parser.add_argument("--gmsh-python", type=Path)
    parser.add_argument("--gmsh-library", type=Path)
    parser.add_argument("--gmsh-executable", type=Path)
    parser.add_argument("--expected-version", default=EXPECTED_GMSH_VERSION)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = run_probe(
            args.output_dir,
            args.runtime_dir,
            args.gmsh_python,
            args.gmsh_library,
            args.gmsh_executable,
            args.expected_version,
        )
    except Exception as error:  # noqa: BLE001 - CLI must always emit a report
        report = write_fail_report(args.output_dir, error, args.expected_version)
        print(f"Gmsh runtime probe failed; inspect {report}", file=sys.stderr)
        return 2
    print(json.dumps(payload, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
