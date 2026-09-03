#!/usr/bin/env python3
"""Run a prepared Templates/LDMOS Id-Vd ablation chain on the Sentaurus VM."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any, Sequence

from run_templates_ldmos_sentaurus_vm import executable, records, sha256_file, write_json


def run(argv: Sequence[str], *, check: bool = True,
        network_retries: int = 3) -> subprocess.CompletedProcess[str]:
    completed: subprocess.CompletedProcess[str] | None = None
    for attempt in range(network_retries + 1):
        completed = subprocess.run(list(argv), check=False, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        transient = completed.returncode == 255 and any(marker in completed.stdout for marker in (
            "Could not resolve hostname", "Connection timed out", "Connection reset",
            "Connection refused", "No route to host",
        ))
        if not transient or attempt == network_retries:
            break
        time.sleep(2.0 * (attempt + 1))
    assert completed is not None
    if check and completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, completed.args, output=completed.stdout)
    return completed


def execute(prepared: Path, parent_run: Path, output: Path, remote_name: str,
            ssh_target: str, ssh_bin: str, scp_bin: str,
            variants: list[str]) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite sealed run: {output}")
    prepared_manifest = json.loads((prepared / "manifest.json").read_text(encoding="utf-8"))
    known = {item["id"] for item in prepared_manifest["single_factor_chain"]}
    invalid = [item for item in variants if item not in known]
    if invalid or not variants:
        raise ValueError(f"invalid variants: {invalid or variants}")
    parent_manifest_path = parent_run / "manifest" / "run_manifest.json"
    parent_manifest = json.loads(parent_manifest_path.read_text(encoding="utf-8"))
    parent_remote = parent_manifest["metadata"]["remote_dir"]
    remote = f"{parent_remote}/{remote_name}"
    quoted_remote = shlex.quote(remote)
    quoted_parent = shlex.quote(parent_remote)
    create = (
        f"umask 077; if test -e {quoted_remote}; then exit 17; fi; "
        f"mkdir -p {quoted_remote}; cp {quoted_parent}/work/n1_fps.tdr {quoted_remote}/"
    )
    run([ssh_bin, ssh_target, create])
    for variant in variants:
        remote_variant = f"{remote}/{variant}"
        run([ssh_bin, ssh_target, f"mkdir -p {shlex.quote(remote_variant)}"])
        run([scp_bin, str(prepared / variant / "IdVd.cmd"),
             str(prepared / variant / "sdevice.par"),
             f"{ssh_target}:{remote_variant}/"])
    results: list[dict[str, Any]] = []
    for variant in variants:
        remote_variant = f"{remote}/{variant}"
        command = (
            f"cd {shlex.quote(remote_variant)}; ln -s ../n1_fps.tdr n1_fps.tdr; "
            "/usr/bin/time -v -o timing.txt sdevice IdVd.cmd > run.out 2>&1; "
            "rc=$?; printf '%s\n' \"$rc\" > exitcode; exit $rc"
        )
        started = time.perf_counter()
        completed = run([ssh_bin, ssh_target, command], check=False)
        results.append({
            "id": variant,
            "status": "pass" if completed.returncode == 0 else "fail",
            "exit_code": completed.returncode,
            "wall_clock_seconds_host": time.perf_counter() - started,
            "ssh_output": completed.stdout.strip(),
        })
    archive_command = (
        f"cd {quoted_remote}; find . -mindepth 2 -maxdepth 2 -type f -print0 | "
        "tar --null -czf idvd_ablation_results.tgz --files-from=-"
    )
    run([ssh_bin, ssh_target, archive_command])
    output.mkdir(parents=True)
    archive = output / "idvd_ablation_results.tgz"
    run([scp_bin, f"{ssh_target}:{remote}/idvd_ablation_results.tgz", str(archive)])
    subprocess.run(["tar", "-xzf", str(archive), "-C", str(output)], check=True)
    manifest = {
        "schema": "vela.templates_ldmos.sentaurus_idvd_ablation_run.v1",
        "classification": "derived_single_factor_control_not_official_oracle",
        "prepared_manifest_sha256": sha256_file(prepared / "manifest.json"),
        "parent_run_manifest_sha256": sha256_file(parent_manifest_path),
        "remote_directory": remote,
        "variants": results,
        "records": records(output),
    }
    write_json(output / "run_manifest.json", manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-dir", type=Path, required=True)
    parser.add_argument("--parent-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--remote-name", required=True)
    parser.add_argument("--variants", default="D1-isothermal,D2-no-hRecVelocity,D5-no-IALMob")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    args = parser.parse_args(argv)
    variants = [item.strip() for item in args.variants.split(",") if item.strip()]
    result = execute(args.prepared_dir.resolve(), args.parent_run.resolve(),
                     args.output_dir.resolve(), args.remote_name, args.ssh_target,
                     args.ssh_bin, args.scp_bin, variants)
    print(json.dumps(result, indent=2))
    return 0 if all(item["status"] == "pass" for item in result["variants"]) else 2


if __name__ == "__main__":
    raise SystemExit(main())
