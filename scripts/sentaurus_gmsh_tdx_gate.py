#!/usr/bin/env python3
"""Internal licensed TDX round-trip gate for a Gmsh preprocessing bundle.

The public workflow never invokes this module.  A live run stages only the
generated DF-ISE ``device.grd/device.dat`` pair in an isolated VM directory,
freezes ``tdx -h``, converts the pair to TDR, imports the returned TDR with the
checked-in Vela importer and compares topology, contacts and node datasets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
SAFE_RUN_ID = re.compile(r"[A-Za-z0-9_.-]+")
SCHEMA = "vela.sentaurus_gmsh_tdx_gate.v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def default_windows_openssh(name: str) -> str:
    if os.name == "nt":
        candidate = Path(os.environ.get("SystemRoot", r"C:\Windows")) / (
            f"System32/OpenSSH/{name}.exe")
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run_checked(argv: Sequence[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return subprocess.run(list(argv), text=True, check=True, **kwargs)


def validate_preprocess(preprocess: Path) -> dict[str, Any]:
    required = [
        preprocess / "device.grd", preprocess / "device.dat",
        preprocess / "mesh.json", preprocess / "doping.csv",
        preprocess / "gmsh_preprocess_manifest.json",
    ]
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"preprocessing bundle is missing: {', '.join(missing)}")
    manifest = json.loads(required[-1].read_text(encoding="utf-8"))
    if manifest.get("mesh_backend") != "gmsh":
        raise ValueError("TDX gate accepts only a Gmsh preprocessing manifest")
    if manifest.get("cross_format", {}).get("status") != "passed":
        raise ValueError("Gmsh preprocessing cross-format gate did not pass")
    for name in ("device.grd", "device.dat"):
        expected = manifest.get("sha256", {}).get(name)
        actual = sha256_file(preprocess / name)
        if expected != actual:
            raise ValueError(f"{name} hash does not match the preprocessing manifest")
    return manifest


def remote_commands(remote_dir: str) -> list[str]:
    return [
        "tdx -h",
        f"mkdir -p {remote_dir}",
        (
            f"cd {remote_dir} && tdx -d device.grd device.dat device.tdr "
            "> tdx.stdout.txt 2> tdx.stderr.txt; code=$?; "
            "printf '%s\\n' \"$code\" > tdx.exit_code.txt; exit $code"
        ),
    ]


def run_gate(args: argparse.Namespace) -> dict[str, Any]:
    preprocess = args.preprocess_dir.resolve()
    preprocess_manifest = validate_preprocess(preprocess)
    if not SAFE_RUN_ID.fullmatch(args.run_id):
        raise ValueError("run ID contains unsupported characters")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    remote_dir = f"{args.remote_root.rstrip('/')}/{args.run_id}"
    commands = remote_commands(remote_dir)
    manifest: dict[str, Any] = {
        "schema": SCHEMA,
        "status": "prepared" if args.dry_run else "running",
        "internal_only": True,
        "dry_run": bool(args.dry_run),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": args.run_id,
        "ssh_target": args.ssh_target,
        "remote_dir": remote_dir,
        "commands": commands,
        "preprocess_manifest_sha256": sha256_file(
            preprocess / "gmsh_preprocess_manifest.json"),
        "input_sha256": {
            name: sha256_file(preprocess / name)
            for name in ("device.grd", "device.dat")
        },
        "gmsh_version": preprocess_manifest.get("gmsh_version"),
    }
    manifest_path = output / "tdx_gate_manifest.json"
    write_json(manifest_path, manifest)
    if args.dry_run:
        return manifest

    ssh = args.ssh_bin or default_windows_openssh("ssh")
    scp = args.scp_bin or default_windows_openssh("scp")
    try:
        help_run = subprocess.run(
            [ssh, args.ssh_target, commands[0]], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        help_path = output / "tdx_help.txt"
        help_path.write_text(help_run.stdout or "", encoding="utf-8")
        if help_run.returncode != 0 or not help_path.stat().st_size:
            raise RuntimeError(
                f"tdx -h failed with exit code {help_run.returncode}")
        run_checked([ssh, args.ssh_target, commands[1]])
        for name in ("device.grd", "device.dat"):
            run_checked([scp, str(preprocess / name),
                         f"{args.ssh_target}:{remote_dir}/{name}"])
        converted = subprocess.run(
            [ssh, args.ssh_target, commands[2]], text=True, check=False)
        for name in ("tdx.stdout.txt", "tdx.stderr.txt", "tdx.exit_code.txt"):
            run_checked([scp, f"{args.ssh_target}:{remote_dir}/{name}", str(output) + "/"])
        if converted.returncode != 0:
            raise RuntimeError(
                f"tdx conversion failed with exit code {converted.returncode}")
        run_checked([
            scp, f"{args.ssh_target}:{remote_dir}/device.tdr", str(output) + "/"])

        export = output / "tdr_export"
        importer = args.tdr_importer.resolve()
        if not importer.is_file():
            raise FileNotFoundError(f"TDR importer does not exist: {importer}")
        run_checked([
            str(importer), "--tdr", str(output / "device.tdr"),
            "--inventory-json", str(export / "inventory.json"),
            "--export-dir", str(export),
        ])
        comparison = output / "comparison"
        compare_run = subprocess.run([
            sys.executable, str(REPO / "scripts" / "compare_sentaurus_tdr_tdx.py"),
            "--tdr-export", str(export), "--tdx-dir", str(preprocess),
            "--output-dir", str(comparison), "--coordinate-epsilon", "1e-12",
        ], text=True, check=False)
        report_path = comparison / "tdr_tdx_comparison.json"
        if not report_path.is_file():
            raise RuntimeError("TDX comparison did not produce its JSON report")
        report = json.loads(report_path.read_text(encoding="utf-8"))
        manifest.update({
            "status": "passed" if compare_run.returncode == 0 else "failed",
            "tdx_help_sha256": sha256_file(help_path),
            "tdx_exit_code": int((output / "tdx.exit_code.txt").read_text().strip()),
            "returned_tdr_sha256": sha256_file(output / "device.tdr"),
            "comparison": report,
        })
        if compare_run.returncode != 0:
            raise RuntimeError("TDX/TDR round-trip comparison failed")
        write_json(manifest_path, manifest)
        return manifest
    except Exception as error:
        manifest["status"] = "failed"
        manifest["error"] = str(error)
        write_json(manifest_path, manifest)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preprocess-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", default="pn2d_gmsh_m4")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default="~/sentaurus_runs/vela_gmsh_tdx")
    parser.add_argument(
        "--tdr-importer", type=Path,
        default=REPO / "build-release" /
        ("sentaurus_import.exe" if os.name == "nt" else "sentaurus_import"),
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = run_gate(args)
    except Exception as error:  # noqa: BLE001 - failed manifest is the contract
        print(str(error), file=sys.stderr)
        return 1
    print(json.dumps(manifest, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
