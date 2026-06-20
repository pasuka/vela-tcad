#!/usr/bin/env python3
"""Run opt-in Sentaurus reference decks on the SSH-accessible VM.

The default path is intentionally conservative: dry-run planning and local
manifest generation work without a VM or license. Live runs stage artifacts in a
build directory and never overwrite committed reference inputs.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REMOTE_ROOT = "~/sentaurus_runs/vela_oracle"
DEFAULT_OUTPUT_DIR = (
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "sentaurus_vm_runs"
)
DEFAULT_SOURCE_DIR = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source"
ARTIFACT_GLOBS = ["*.tdr", "*.plt", "*.log", "*.grd", "*.dat", "run_pn2d_*.out"]
COMMON_FILES = ["pn2d_sde.cmd", "models.par"]
STAGE_FILES = {
    "0v": "pn2d_0v_sdevice.cmd",
    "iv": "pn2d_iv_sdevice.cmd",
    "bv": "pn2d_bv_sdevice.cmd",
}


def default_windows_openssh(name: str) -> str:
    if os.name == "nt":
        candidate = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "OpenSSH" / f"{name}.exe"
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def split_stages(raw: str) -> list[str]:
    stages = [part.strip() for part in raw.split(",") if part.strip()]
    if not stages:
        raise ValueError("--stages must contain at least one stage")
    invalid = [stage for stage in stages if stage not in STAGE_FILES]
    if invalid:
        raise ValueError(f"unsupported stage(s): {', '.join(invalid)}")
    return stages


def normalize_remote_root(remote_root: str) -> str:
    return remote_root.rstrip("/")


def remote_source_dir(remote_root: str, run_id: str) -> str:
    return f"{normalize_remote_root(remote_root)}/{run_id}/source"


def stage_commands(remote_dir: str, stages: Sequence[str]) -> list[str]:
    commands = [f"cd {remote_dir} && sde -e -l pn2d_sde.cmd"]
    for stage in stages:
        commands.append(
            f"cd {remote_dir} && sdevice {STAGE_FILES[stage]} > run_pn2d_{stage}.out 2>&1"
        )
    return commands


def required_files(stages: Sequence[str]) -> list[str]:
    files = list(COMMON_FILES)
    for stage in stages:
        files.append(STAGE_FILES[stage])
    return files


def validate_required_files(source_dir: Path, stages: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    for name in required_files(stages):
        path = source_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"missing required source file: {path}")
        paths.append(path)
    return paths


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def run_checked(argv: Sequence[str]) -> None:
    subprocess.run(list(argv), check=True)


def run_live(args: argparse.Namespace,
             source_files: Sequence[Path],
             remote_dir: str,
             local_source_copy_dir: Path) -> list[str]:
    ssh_bin = args.ssh_bin or default_windows_openssh("ssh")
    scp_bin = args.scp_bin or default_windows_openssh("scp")
    warnings: list[str] = []

    local_source_copy_dir.mkdir(parents=True, exist_ok=True)
    run_checked([ssh_bin, args.ssh_target, f"mkdir -p {remote_dir}"])
    for path in source_files:
        run_checked([scp_bin, str(path), f"{args.ssh_target}:{remote_dir}/"])

    for command in stage_commands(remote_dir, args.stages):
        run_checked([ssh_bin, args.ssh_target, command])

    for pattern in ARTIFACT_GLOBS:
        completed = subprocess.run(
            [scp_bin, f"{args.ssh_target}:{remote_dir}/{pattern}", str(local_source_copy_dir) + "/"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if completed.returncode != 0:
            warnings.append(
                f"artifact copy warning for {pattern}: {completed.stderr.strip() or completed.stdout.strip()}"
            )
    return warnings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("device", choices=["pn2d"])
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None, help="SSH executable; defaults to Windows OpenSSH when available")
    parser.add_argument("--scp-bin", default=None, help="SCP executable; defaults to Windows OpenSSH when available")
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--local-output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--stages", default="bv", help="Comma-separated stages: 0v,iv,bv")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.stages = split_stages(args.stages)
        source_dir = args.source_dir.resolve()
        source_files = validate_required_files(source_dir, args.stages)
        run_id = args.run_id or datetime.now().strftime("pn2d_vm_%Y%m%d_%H%M%S")
        remote_dir = remote_source_dir(args.remote_root, run_id)
        local_run_dir = args.local_output_dir.resolve() / run_id
        local_source_copy_dir = local_run_dir / "source"
        commands = stage_commands(remote_dir, args.stages)

        warnings: list[str] = []
        if not args.dry_run:
            warnings = run_live(args, source_files, remote_dir, local_source_copy_dir)

        manifest = {
            "schema": "vela.sentaurus_vm_run.v1",
            "device": args.device,
            "dry_run": args.dry_run,
            "ssh_target": args.ssh_target,
            "remote_root": args.remote_root,
            "remote_source_dir": remote_dir,
            "run_id": run_id,
            "stages": args.stages,
            "source_dir": str(source_dir),
            "local_run_dir": str(local_run_dir),
            "local_source_copy_dir": str(local_source_copy_dir),
            "required_files": [path.name for path in source_files],
            "commands": commands,
            "warnings": warnings,
        }
        write_manifest(local_run_dir / "sentaurus_vm_run_manifest.json", manifest)
        print(json.dumps(manifest, indent=2))
        return 0
    except Exception as exc:  # noqa: BLE001 - CLI should report cleanly.
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
