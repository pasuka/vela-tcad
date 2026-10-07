"""Collect license texts for the installed runtime and compiled dependencies."""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import shutil
import subprocess

# MSYS2 splits the GCC license texts among runtime packages.
DLL_LICENSES = {
    "libgcc": "libgcc", "libstdc++": "libstdc++", "libwinpthread": "libwinpthread",
    "libgfortran": "libgfortran", "libgomp": "libgomp", "libquadmath": "libquadmath",
    "libspdlog": "spdlog", "libfmt": "fmt", "libhdf5": "hdf5", "libaec": "libaec",
    "zlib": "zlib", "libbz2": "bzip2", "liblzma": "xz",
    "libopenblas": "OpenBLAS", "libmetis": "metis", "libscotch": "scotch",
    "libesmumps": "scotch", "libmumps": "mumps", "libstrumpack": "strumpack",
    "libsuperlu": "superlu_mt",
}
BASE_LICENSES = ("boost", "eigen", "json", "spdlog", "highfive", "hdf5", "zlib", "libaec")
DEBIAN_BASE = ("libboost-dev", "libeigen3-dev", "nlohmann-json3-dev", "libspdlog-dev",
               "libhdf5-dev", "zlib1g", "libaec0")

def generate(stage, build, license_roots=(), *, extra_linux_libraries=()):
    stage, build = Path(stage), Path(build)
    target = stage / "share/vela/licenses"
    target.mkdir(parents=True, exist_ok=True)
    entries = []

    def collect(name, source):
        source = Path(source)
        if not source.exists():
            raise ValueError(f"Missing license text for {name}: {source}")
        dest = target / name
        if not dest.exists():
            if source.is_dir():
                shutil.copytree(source, dest)
            else:
                dest.mkdir()
                shutil.copy2(source, dest / source.name)
            entries.append(name)

    if (stage / "bin/vela.exe").exists():
        roots = [Path(p) for p in license_roots]
        cache = (build / "CMakeCache.txt").read_text(encoding="utf-8")
        match = re.search(r"^CMAKE_CXX_COMPILER:(?:FILEPATH|STRING)=(.+)$", cache, re.M)
        if match is None:
            raise ValueError("CMake cache has no configured C++ compiler path")
        compiler = match.group(1)
        roots.append(Path(compiler).parent.parent / "share/licenses")
        required = set(BASE_LICENSES)
        for dll in (stage / "bin").glob("*.dll"):
            name = dll.name.lower()
            if any(name.startswith("lib" + p) for p in ("amd", "camd", "ccolamd", "cholmod", "colamd", "spqr", "umfpack", "suitesparseconfig")):
                required.add("suitesparse")
                continue
            matches = [license for prefix, license in DLL_LICENSES.items() if name.startswith(prefix)]
            if not matches:
                raise ValueError(f"No license mapping for runtime dependency {dll.name}")
            required.update(matches)
        for name in sorted(required):
            source = next((r / name for r in roots if (r / name).is_dir()), None)
            if source is None:
                raise ValueError(f"Missing MSYS2 license directory: {name}")
            collect(name, source)
    else:
        # Ask dpkg about the actual ldd closure, including dependencies with
        # SONAME symlinks (the installed package contains copies of those files).
        libraries = {str(path) for path in extra_linux_libraries}
        for name in ("vela_example_runner", "sentaurus_import"):
            result = subprocess.run(["ldd", str(build / name)], check=True, capture_output=True, text=True)
            libraries.update(re.findall(r"=> (/.+?) \(0x[0-9a-fA-F]+\)", result.stdout))
        packages = set(DEBIAN_BASE)
        for library in sorted(libraries):
            owner = subprocess.run(["dpkg-query", "-S", library], capture_output=True, text=True)
            if owner.returncode:
                owner = subprocess.run(["dpkg-query", "-S", str(Path(library).resolve())], check=True, capture_output=True, text=True)
            packages.add(owner.stdout.split(": ")[0])
        for package in sorted(packages):
            name = package.split(":")[0]
            collect(name, Path("/usr/share/doc") / name / "copyright")
        # Debian copyright records refer to these full shared license texts.
        collect("debian-common-licenses", Path("/usr/share/common-licenses"))
        highfive = next((p for p in [build / "_deps/highfive-src/LICENSE", *[Path(r) / "highfive/LICENSE" for r in license_roots]] if p.is_file()), None)
        if highfive is None:
            raise ValueError("Provide --license-root containing highfive/LICENSE")
        collect("highfive", highfive)

    # Both official runtimes contain their complete CPython license file.
    python = stage / "bin/python"
    texts = [p for p in python.rglob("*") if p.is_file() and ("license" in p.name.lower() or "copyright" in p.name.lower())]
    if not texts:
        raise ValueError("Python archive is missing license texts")
    for index, source in enumerate(sorted(texts)):
        collect(f"python-{index:03d}", source)
    lines = ["THIRD_PARTY_NOTICES", "",
        "Vela includes the third-party components whose unmodified license texts are listed below.",
        "CPython and its bundled-library license texts (including OpenSSL) are included below.",
        "GCC runtime copies include the GCC Runtime Library Exception.",
        "These notices do not replace the terms in the individual license texts.", ""]
    for name in entries:
        lines += [f"Component: {name}"]
        for source in sorted((target / name).rglob("*")):
            if source.is_file():
                lines += [f"--- {source.relative_to(stage).as_posix()} ---",
                          source.read_text(encoding="utf-8", errors="replace"), ""]
    (stage / "THIRD_PARTY_NOTICES").write_text("\n".join(lines), encoding="utf-8")
    return entries

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--license-root", type=Path, action="append", default=[])
    args = parser.parse_args()
    generate(args.stage, args.build_dir, args.license_root)

if __name__ == "__main__":
    main()
