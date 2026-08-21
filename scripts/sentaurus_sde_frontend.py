#!/usr/bin/env python3
"""Sentaurus 2022 evaluated-log and offline-static frontend orchestration."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from sentaurus_device_ir import parse_sde_device_ir


MANIFEST_SCHEMA = "vela.sentaurus_sde_frontend_manifest.v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def offline_static_frontend(source: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    source = source.resolve()
    ir = parse_sde_device_ir(source)
    return ir, {
        "schema": MANIFEST_SCHEMA,
        "frontend": "offline-static",
        "sde_version": None,
        "source": str(source),
        "source_sha256": sha256_file(source),
        "evaluated_log": None,
        "evaluated_log_sha256": None,
        "capability": "declarative Scheme subset; unsupported forms fail closed",
    }


def sentaurus2022_evaluated_frontend(source: Path, output_dir: Path,
                                     executable: str = "sde",
                                     build_mesh: bool = False) \
        -> tuple[dict[str, Any], dict[str, Any]]:
    """Execute SDE and parse its variable-expanded ``.log.eval`` command stream."""
    source = source.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    wrapper = output_dir / "vela_sde_evaluated_wrapper.scm"
    escaped = str(source).replace("\\", "/").replace('"', '\\"')
    forms = ["(set! evaluate-log-file #t)", f'(load "{escaped}")']
    if build_mesh:
        prefix = (output_dir / "sentaurus_official").as_posix().replace('"', '\\"')
        forms.append(f'(sde:build-mesh "{prefix}")')
    wrapper.write_text("\n".join(forms) + "\n", encoding="utf-8")

    version_run = subprocess.run(
        [executable, "-version"], cwd=output_dir, capture_output=True, text=True,
        check=False)
    completed = subprocess.run(
        [executable, "-e", "-l", str(wrapper)], cwd=output_dir,
        capture_output=True, text=True, check=False)
    (output_dir / "sde.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "sde.stderr.txt").write_text(completed.stderr, encoding="utf-8")

    logs = sorted(output_dir.glob("*.log.eval"), key=lambda path: path.stat().st_mtime)
    manifest: dict[str, Any] = {
        "schema": MANIFEST_SCHEMA,
        "frontend": "sentaurus2022-evaluated",
        "sde_version": (version_run.stdout or version_run.stderr).strip(),
        "source": str(source),
        "source_sha256": sha256_file(source),
        "wrapper": str(wrapper),
        "command": [executable, "-e", "-l", str(wrapper)],
        "exit_code": completed.returncode,
        "stdout": "sde.stdout.txt",
        "stderr": "sde.stderr.txt",
        "evaluated_log": str(logs[-1]) if logs else None,
        "evaluated_log_sha256": sha256_file(logs[-1]) if logs else None,
    }
    (output_dir / "sde_frontend_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Sentaurus SDE exited with {completed.returncode}; see "
            f"{output_dir / 'sde.stderr.txt'}")
    if not logs:
        raise RuntimeError(
            "Sentaurus SDE succeeded but produced no .log.eval; verify that "
            "evaluate-log-file is supported by the installed T-2022.03 build")
    return parse_sde_device_ir(logs[-1], frontend="sentaurus2022-evaluated"), manifest
