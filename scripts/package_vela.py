"""Build an offline, relocatable Vela runtime distribution from CMake install."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import posixpath
import re
import shutil
import subprocess
import tarfile
import tempfile
import zipfile

from generate_third_party_notices import generate

def sha256(path):
    with Path(path).open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        return digest.hexdigest()

def tar_permissions(item):
    """Normalize archive ownership/modes even when staging on a Windows mount."""
    item.uid = item.gid = 0
    item.uname = item.gname = "root"
    if item.isdir():
        item.mode = 0o755
    elif item.isfile():
        name = PurePosixPath(item.name)
        executable = (item.name in ("vela/vela", "vela/bin/vela", "vela/bin/vela_example_runner", "vela/bin/sentaurus_import")
                      or item.name.startswith("vela/bin/python/bin/") or ".so" in name.name)
        item.mode = 0o755 if executable else 0o644
    return item

def safe_member(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise ValueError(f"Unsafe archive member: {name}")
    return path

def unpack_python(archive, destination, system):
    destination.mkdir(parents=True)
    if system == "win64":
        with zipfile.ZipFile(archive) as zipped:
            for item in zipped.infolist():
                safe_member(item.filename)
                if (item.external_attr >> 16) & 0o170000 == 0o120000:
                    raise ValueError(f"Python zip symlink is unsupported: {item.filename}")
            zipped.extractall(destination)
        required = ("python.exe", "python312.dll", "python312.zip", "python312._pth")
        if any(not (destination / name).is_file() for name in required):
            raise ValueError("Expected the official CPython 3.12 amd64 embeddable zip")
        # Extension modules (_socket, _ctypes, etc.) live beside python.exe.
        # All paths stay in the package; no import site or environment paths.
        (destination / "python312._pth").write_text("python312.zip\n.\n..\n", encoding="ascii")
    else:
        with tarfile.open(archive) as tar:
            members = tar.getmembers()
            for item in members:
                path = safe_member(item.name)
                if not path.parts or path.parts[0] != "python":
                    raise ValueError(f"Expected install_only archive rooted at python/: {item.name}")
                if item.isdev() or item.isfifo() or item.islnk():
                    raise ValueError(f"Unsupported Python archive member: {item.name}")
                if item.issym():
                    # Standalone terminfo links use ../ within python/share/.
                    link = PurePosixPath(item.linkname)
                    resolved = posixpath.normpath(posixpath.join(str(path.parent), item.linkname))
                    if link.is_absolute() or "\\" in item.linkname or not resolved.startswith("python/"):
                        raise ValueError(f"Escaping Python symlink: {item.name}")
            # Keep the archive's python/ layout, extracting relative to bin/.
            tar.extractall(destination.parent, filter="data")
        if not (destination / "bin/python3").is_file():
            raise ValueError("Expected python-build-standalone CPython 3.12 install_only archive")
    interpreter = destination / ("python.exe" if system == "win64" else "bin/python3")
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
    done = subprocess.run([str(interpreter), "-B", "-I", "-S", "-c",
        "import json,sys,platform; print(json.dumps([list(sys.version_info[:2]),platform.machine()]))"],
        env=env, check=True, capture_output=True, text=True)
    version, machine = json.loads(done.stdout)
    if version != [3, 12] or machine.lower() not in ("amd64", "x86_64"):
        raise ValueError(f"Expected CPython 3.12 x86_64, found {version} {machine}")

def distribution_platform():
    if os.name == "nt":
        return "win64"
    info = platform.freedesktop_os_release()
    if info.get("ID") != "ubuntu" or info.get("VERSION_ID") != "26.04":
        raise ValueError("Ubuntu package must be built on Ubuntu 26.04")
    if platform.machine() != "x86_64":
        raise ValueError("Ubuntu package requires x86_64")
    return "ubuntu26.04-x86_64"

def audit(stage):
    forbidden = {"reference_staging", "reference_tcad", ".git"}
    names = []
    for path in sorted(stage.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(stage)
        if "__pycache__" in relative.parts and relative.parts[:2] != ("bin", "python"):
            raise ValueError(f"Generated source cache in distribution: {relative}")
        if forbidden.intersection(relative.parts) or path.suffix.lower() in (".tdr", ".tdx", ".plt", ".grd", ".dat", ".cmd"):
            raise ValueError(f"Reference/proprietary artifact in distribution: {relative}")
        names.append({"path": relative.as_posix(), "size": path.stat().st_size, "sha256": sha256(path)})
    with zipfile.ZipFile(stage / "bin/vela_sdevice.pyz") as archive:
        for name in archive.namelist():
            if (forbidden | {"__pycache__"}).intersection(PurePosixPath(name).parts) or name.lower().endswith((".tdr", ".tdx", ".plt", ".cmd")):
                raise ValueError(f"Reference artifact inside frontend archive: {name}")
    return names

def package(args):
    build, python_archive, output = args.build_dir.resolve(), args.python_archive.resolve(), args.output_dir.resolve()
    actual = sha256(python_archive)
    if actual != args.python_sha256.lower():
        raise ValueError(f"Python archive SHA-256 mismatch: {actual}")
    system = distribution_platform()
    output.mkdir(parents=True, exist_ok=True)
    # POSIX Python distributions contain case-sensitive paths and symlinks.
    # Use native /tmp even when the output directory is a Windows-mounted tree.
    with tempfile.TemporaryDirectory(prefix="vela_package_", dir=output if os.name == "nt" else None) as temp:
        stage = Path(temp) / "vela"
        subprocess.run([args.cmake, "--install", str(build), "--prefix", str(stage),
                        "--config", args.config, "--component", "SDeviceRuntime"], check=True)
        suffix = ".exe" if system == "win64" else ""
        for name in ("vela", "vela_example_runner", "sentaurus_import"):
            if not (stage / "bin" / (name + suffix)).is_file():
                raise ValueError(f"Missing required executable {name}; build with detected HDF5/TDR support")
        unpack_python(python_archive, stage / "bin/python", system)
        # install_only / embeddable archives omit some static-library notices.
        # Supply the matching upstream license bundle offline as well.
        if not (args.python_license_dir / "LICENSE.openssl-3.txt").is_file():
            raise ValueError("Python license directory must include LICENSE.openssl-3.txt")
        shutil.copytree(args.python_license_dir, stage / "bin/python/licenses")
        licenses = generate(stage, build, args.license_root)
        with zipfile.ZipFile(stage / "bin/vela_sdevice.pyz") as zipped:
            version = json.loads(zipped.read("vela_sdevice/_assets/version.json"))
        if not re.fullmatch(r"[0-9a-fA-F]{40}", version.get("commit", "")):
            raise ValueError("Frontend source commit is unavailable; configure VELA_SOURCE_REVISION from verified Git HEAD")
        (stage / "VERSION").write_text(json.dumps(version, indent=2) + "\n", encoding="utf-8")
        manifest = {**version, "platform": system, "python_archive": python_archive.name,
                    "python_sha256": actual, "licenses": licenses, "files": audit(stage)}
        (stage / "PACKAGE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        basename = f"vela-{version['version']}-{system}"
        archive = output / (basename + (".zip" if system == "win64" else ".tar.gz"))
        temporary = archive.with_name(archive.name + ".tmp")
        if system == "win64":
            with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as zipped:
                for path in sorted(stage.rglob("*")):
                    if path.is_file():
                        zipped.write(path, "vela/" + path.relative_to(stage).as_posix())
        else:
            with tarfile.open(temporary, "w:gz", dereference=False) as tar:
                tar.add(stage, arcname="vela", filter=tar_permissions)
        os.replace(temporary, archive)
        digest = sha256(archive)
        archive.with_name(archive.name + ".sha256").write_text(f"{digest}  {archive.name}\n", encoding="ascii")
        # Keep an inspectable manifest alongside the compressed package.
        shutil.copy2(stage / "PACKAGE_MANIFEST.json", output / (basename + ".manifest.json"))
        return archive

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--python-archive", type=Path, required=True)
    parser.add_argument("--python-sha256", required=True, help="Trusted SHA-256 of the offline CPython archive")
    parser.add_argument("--python-license-dir", type=Path, required=True,
                        help="Offline upstream CPython/bundled-library license texts")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--config", default="Release")
    parser.add_argument("--license-root", type=Path, action="append", default=[])
    args = parser.parse_args()
    try:
        print(package(args))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"package_vela: {error}\n")

if __name__ == "__main__":
    main()
