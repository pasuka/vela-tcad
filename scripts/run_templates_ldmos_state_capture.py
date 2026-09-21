#!/usr/bin/env python3
"""Execute output-only Templates/LDMOS representative-state decks on the VM."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any, Sequence

from run_templates_ldmos_sentaurus_vm import executable, records, safe_extract, sha256_file, write_json


STAGES = (("idvg", "IdVg.cmd"), ("idvd", "IdVd.cmd"), ("bv", "BVdss.cmd"))


def run(argv: Sequence[str], check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv), check=check, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--state-decks", type=Path, required=True)
    parser.add_argument("--stages", default="idvg,idvd,bv")
    parser.add_argument("--output-name", default="representative_states")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    state_decks = args.state_decks.resolve()
    requested = [item.strip().lower() for item in args.stages.split(",") if item.strip()]
    invalid = [item for item in requested if item not in {stage for stage, _ in STAGES}]
    if invalid or not requested:
        raise ValueError(f"invalid state-capture stages: {invalid or requested}")
    stages = [item for item in STAGES if item[0] in requested]
    if not args.output_name or Path(args.output_name).name != args.output_name:
        raise ValueError("--output-name must be one directory name")
    output = run_dir / args.output_name
    sealed_outputs = (
        output / "state_capture_manifest.json",
        output / "state_capture_results.tgz",
        output / "raw",
    )
    if any(path.exists() for path in sealed_outputs):
        raise FileExistsError(f"refusing to overwrite sealed representative states: {output}")
    run_manifest = json.loads((run_dir / "manifest" / "run_manifest.json").read_text(encoding="utf-8"))
    remote_parent = run_manifest["metadata"]["remote_dir"]
    remote = f"{remote_parent}/{args.output_name}"
    quoted_remote = shlex.quote(remote)
    quoted_parent = shlex.quote(remote_parent)
    create = (
        f"umask 077; if test -e {quoted_remote}; then exit 17; fi; "
        f"mkdir -p {quoted_remote}; cp {quoted_parent}/work/n1_fps.tdr {quoted_remote}/"
    )
    run([args.ssh_bin, args.ssh_target, create])
    upload = [state_decks / name for _, name in stages] + [state_decks / "sdevice.par"]
    run([args.scp_bin, *(str(path) for path in upload), f"{args.ssh_target}:{remote}/"])
    results: list[dict[str, Any]] = []
    for stage, deck in stages:
        command = (
            f"cd {quoted_remote}; /usr/bin/time -v -o timing_{stage}.txt "
            f"sdevice {shlex.quote(deck)} > run_{stage}.out 2>&1; "
            f"rc=$?; printf '%s\\n' \"$rc\" > {stage}.exitcode; exit $rc"
        )
        started = time.perf_counter()
        completed = run([args.ssh_bin, args.ssh_target, command], check=False)
        results.append({
            "stage": stage,
            "return_code": completed.returncode,
            "wall_clock_seconds_host": time.perf_counter() - started,
            "status": "pass" if completed.returncode == 0 else "fail",
            "ssh_output": completed.stdout,
        })
    archive_command = (
        f"cd {quoted_remote}; find . -maxdepth 1 -type f "
        "! -name state_capture_results.tgz -print0 | "
        "tar --null -czf state_capture_results.tgz --files-from=-"
    )
    run([args.ssh_bin, args.ssh_target, archive_command])
    output.mkdir(parents=True, exist_ok=True)
    archive = output / "state_capture_results.tgz"
    run([args.scp_bin, f"{args.ssh_target}:{remote}/state_capture_results.tgz", str(archive)])
    raw = output / "raw"
    safe_extract(archive, raw)
    manifest = {
        "schema": "vela.templates_ldmos.state_capture_manifest.v1",
        "classification": "derived_output_only_control_not_official_oracle",
        "parent_run_id": run_dir.name,
        "parent_run_manifest_sha256": sha256_file(run_dir / "manifest" / "run_manifest.json"),
        "state_deck_manifest_sha256": sha256_file(state_decks / "state_deck_manifest.json"),
        "remote_directory": remote,
        "stage_results": results,
        "records": records(raw),
    }
    write_json(output / "state_capture_manifest.json", manifest)
    print(json.dumps({"results": results, "tdr_count": len(list(raw.glob("*.tdr")))}, indent=2))
    return 0 if all(item["status"] == "pass" for item in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
