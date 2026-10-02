"""Build a stdlib-only SDevice frontend archive with its resource closure."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipapp

REPO = Path(__file__).resolve().parents[1]
# The emitters' complete non-stdlib import closure. convert_tcad_export is a
# deferred import on the TDR path; the DF-ISE parser needs compare_sentaurus_tdr_tdx.
MODULES = ("sentaurus_parameter_ir.py", "sentaurus_parameter_map.py",
           "compare_sentaurus_tdr_tdx.py", "convert_tcad_export.py")


def build_archive(target: Path) -> Path:
    target = Path(target).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=target.parent, prefix="sdevice_pack_") as temp:
        root = Path(temp)
        shutil.copytree(REPO / "scripts/vela_sdevice", root / "vela_sdevice",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for name in MODULES:
            shutil.copy2(REPO / "scripts" / name, root / name)
        assets = root / "vela_sdevice/_assets"
        for name in ("schemas", "configs/schema", "share/vela/materials"):
            shutil.copytree(REPO / name, assets / name)
        version = re.search(r"project\([^\s)]+\s+VERSION\s+(\S+)",
                            (REPO / "CMakeLists.txt").read_text()).group(1)
        revision = subprocess.run(["git", "-c", "core.fsmonitor=false", "rev-parse", "HEAD"],
                                  cwd=REPO, capture_output=True, text=True)
        (assets / "version.json").write_text(json.dumps({"version": version,
            "commit": revision.stdout.strip() if revision.returncode == 0 else "unknown"}), encoding="utf-8")
        zipapp.create_archive(root, target=target, main="vela_sdevice.__main__:main", interpreter=None)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=REPO / "build/bin/vela_sdevice.pyz")
    args = parser.parse_args()
    print(build_archive(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
