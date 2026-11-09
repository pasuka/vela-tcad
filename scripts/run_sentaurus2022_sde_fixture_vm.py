#!/usr/bin/env python3
"""Generate a provenance-sealed Sentaurus 2022 SDE oracle over SSH.

The runner is intentionally independent of a local Sentaurus installation.  It
stages one owned SDE command file in an isolated VM directory, enables the
official evaluated command log, executes SDE, and retrieves the geometry and
mesh artifacts without modifying the committed source deck.
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
from typing import Sequence

REPO = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402


DEFAULT_SOURCE = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sde.cmd"
DEFAULT_OUTPUT = REPO / "build-release" / "sentaurus2022-fixtures"
DEFAULT_REMOTE_ROOT = "~/sentaurus_runs/vela_sde_fixtures"
ARTIFACT_GLOBS = (
    "*.log.eval",
    "*_bnd.tdr",
    "*_msh.cmd",
    "*_msh.tdr",
    "sde.stdout.txt",
    "sde.stderr.txt",
    "sde.version.txt",
    "sde.exit_code.txt",
)


def default_windows_openssh(name: str) -> str:
    if os.name == "nt":
        candidate = (
            Path(os.environ.get("SystemRoot", r"C:\Windows"))
            / "System32" / "OpenSSH" / f"{name}.exe"
        )
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def remote_run_dir(remote_root: str, run_id: str) -> str:
    return f"{remote_root.rstrip('/')}/{run_id}"


def portable_source_path(source: Path) -> str:
    try:
        return source.relative_to(REPO).as_posix()
    except ValueError:
        return str(source)


def mesh_summary(stdout_path: Path) -> dict[str, int] | None:
    if not stdout_path.is_file():
        return None
    text = stdout_path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(
        r"The final mesh has:\s*([0-9]+) points\.\s*([0-9]+) elements\.\s*([0-9]+) regions\.",
        text,
    )
    if not matches:
        return None
    points, elements, regions = matches[-1]
    return {
        "points": int(points),
        "elements": int(elements),
        "regions": int(regions),
    }


def portable_ir_provenance(value: object, evaluated_log: Path) -> object:
    """Replace host-specific evaluated-log paths throughout a Device IR."""
    portable = "artifacts/" + evaluated_log.name
    original = str(evaluated_log)
    if isinstance(value, dict):
        return {
            key: portable_ir_provenance(item, evaluated_log)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [portable_ir_provenance(item, evaluated_log) for item in value]
    if isinstance(value, str):
        return value.replace(original, portable)
    return value


def evaluated_wrapper(source_name: str, source_text: str) -> str:
    """Return one evaluated unit because T-2022.03 does not flatten ``load``.

    On the licensed T-2022.03-SP2 VM, an evaluated wrapper containing ``load``
    records only the load form.  Embedding the byte-for-byte source text after
    enabling the log makes SDE emit the required variable-expanded command
    stream while the separately copied owned source remains hash-sealed.
    """
    return (
        "(set! evaluate-log-file #t)\n"
        f"; Vela owned source follows: {source_name}\n"
        f"{source_text.rstrip()}\n"
    )


def remote_commands(remote_dir: str) -> list[str]:
    return [
        f"mkdir -p {remote_dir}",
        f"cd {remote_dir} && sde -v > sde.version.txt 2>&1",
        (
            f"cd {remote_dir} && "
            "sde -e -l vela_sde_evaluated_wrapper.scm "
            "> sde.stdout.txt 2> sde.stderr.txt; "
            "code=$?; printf '%s\\n' \"$code\" > sde.exit_code.txt; exit $code"
        ),
    ]


def run_checked(argv: Sequence[str]) -> None:
    subprocess.run(list(argv), check=True)


def copy_remote_artifacts(
    scp_bin: str, ssh_target: str, remote_dir: str, artifact_dir: Path
) -> list[str]:
    warnings: list[str] = []
    for pattern in ARTIFACT_GLOBS:
        completed = subprocess.run(
            [scp_bin, f"{ssh_target}:{remote_dir}/{pattern}", str(artifact_dir) + "/"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if completed.returncode != 0:
            warnings.append(
                f"artifact copy warning for {pattern}: "
                f"{completed.stderr.strip() or completed.stdout.strip()}"
            )
    return warnings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--fixture-name", default="pn2d_sentaurus2022")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--local-output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        source = args.source.resolve()
        if not source.is_file():
            raise FileNotFoundError(f"SDE source does not exist: {source}")
        run_id = args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        local_run_dir = args.local_output_dir.resolve() / args.fixture_name / run_id
        input_dir = local_run_dir / "input"
        artifact_dir = local_run_dir / "artifacts"
        input_dir.mkdir(parents=True, exist_ok=True)
        artifact_dir.mkdir(parents=True, exist_ok=True)
        source_copy = input_dir / source.name
        shutil.copy2(source, source_copy)
        wrapper = input_dir / "vela_sde_evaluated_wrapper.scm"
        wrapper.write_text(
            evaluated_wrapper(source.name, source.read_text(encoding="utf-8")),
            encoding="utf-8",
        )

        remote_dir = remote_run_dir(args.remote_root, f"{args.fixture_name}_{run_id}")
        commands = remote_commands(remote_dir)
        warnings: list[str] = []
        if not args.dry_run:
            ssh_bin = args.ssh_bin or default_windows_openssh("ssh")
            scp_bin = args.scp_bin or default_windows_openssh("scp")
            run_checked([ssh_bin, args.ssh_target, commands[0]])
            for staged in (source_copy, wrapper):
                run_checked([scp_bin, str(staged), f"{args.ssh_target}:{remote_dir}/"])
            run_checked([ssh_bin, args.ssh_target, commands[1]])
            sde_run = subprocess.run([ssh_bin, args.ssh_target, commands[2]], check=False)
            warnings.extend(
                copy_remote_artifacts(
                    scp_bin, args.ssh_target, remote_dir, artifact_dir
                )
            )
            if sde_run.returncode != 0:
                raise RuntimeError(
                    f"Sentaurus SDE exited with {sde_run.returncode}; "
                    f"retrieved diagnostics are in {artifact_dir}"
                )

        files = sorted(
            path for path in local_run_dir.rglob("*")
            if path.is_file() and path.name != "manifest.json"
        )
        version_path = artifact_dir / "sde.version.txt"
        exit_path = artifact_dir / "sde.exit_code.txt"
        evaluated_logs = sorted(artifact_dir.glob("*.log.eval"))
        boundary_tdrs = sorted(artifact_dir.glob("*_bnd.tdr"))
        mesh_commands = sorted(artifact_dir.glob("*_msh.cmd"))
        mesh_tdrs = sorted(artifact_dir.glob("*_msh.tdr"))
        expected_ir_path = local_run_dir / "expected_device_ir.json"
        if evaluated_logs:
            expected_ir = parse_sde_device_ir(
                evaluated_logs[-1], frontend="sentaurus2022-evaluated"
            )
            expected_ir = portable_ir_provenance(
                expected_ir, evaluated_logs[-1]
            )
            expected_ir_path.write_text(
                json.dumps(expected_ir, indent=2) + "\n", encoding="utf-8"
            )
            files = sorted(
                path for path in local_run_dir.rglob("*")
                if path.is_file() and path.name != "manifest.json"
            )
        manifest = {
            "schema": "vela.sentaurus2022_sde_oracle.v1",
            "fixture": args.fixture_name,
            "run_id": run_id,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "dry_run": args.dry_run,
            "ssh_target": args.ssh_target,
            "remote_run_dir": remote_dir,
            "source_original": portable_source_path(source),
            "source_sha256": sha256_file(source),
            "wrapper_sha256": sha256_file(wrapper),
            "commands": commands,
            "sde_version": (
                version_path.read_text(encoding="utf-8", errors="replace").strip()
                if version_path.is_file() else None
            ),
            "exit_code": (
                int(exit_path.read_text(encoding="utf-8").strip())
                if exit_path.is_file() else None
            ),
            "evaluated_logs": [path.name for path in evaluated_logs],
            "boundary_tdrs": [path.name for path in boundary_tdrs],
            "mesh_commands": [path.name for path in mesh_commands],
            "mesh_tdrs": [path.name for path in mesh_tdrs],
            "mesh_summary": mesh_summary(artifact_dir / "sde.stdout.txt"),
            "expected_device_ir": (
                expected_ir_path.relative_to(local_run_dir).as_posix()
                if expected_ir_path.is_file() else None
            ),
            "files": {
                path.relative_to(local_run_dir).as_posix(): {
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for path in files
            },
            "warnings": warnings,
        }
        (local_run_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(manifest, indent=2))
        if not args.dry_run and (
            not evaluated_logs
            or evaluated_logs[-1].stat().st_size <= 50
            or not boundary_tdrs
            or not mesh_commands
            or not mesh_tdrs
        ):
            print(
                "SDE succeeded but the evaluated log or required geometry/mesh oracle "
                "artifact is missing",
                file=sys.stderr,
            )
            return 1
        return 0
    except Exception as exc:  # noqa: BLE001 - CLI must emit a concise failure.
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
