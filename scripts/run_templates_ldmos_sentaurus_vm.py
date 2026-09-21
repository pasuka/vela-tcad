#!/usr/bin/env python3
"""Run and seal the T-2022.03-SP2 Applications Library Templates/LDMOS oracle.

Without ``--live`` the command materializes an isolated bundle from
``--source-dir`` and writes reviewable manifests. A live run copies the
Applications Library source from the configured SSH target, never edits it in
place, executes the selected stages in a unique remote directory, and stores
all proprietary artifacts below ignored ``reference_staging/``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from templates_ldmos_contracts import (
    BENCHMARK,
    draft_governance_contracts,
    render_summary,
    validate_document,
)


REPO = Path(__file__).resolve().parents[1]
DEFAULT_STAGING = REPO / "reference_staging" / "templates_ldmos_sentaurus2022"
DEFAULT_APPLICATIONS_PATH = (
    "/atctools/Synopsys/tcad/T-2022.03/tcad/T-2022.03-SP2/"
    "Applications_Library/Templates/LDMOS"
)
DEFAULT_REMOTE_ROOT = "/root/sentaurus_runs/vela_oracle_2022/templates_ldmos"
EXPECTED_VERSION = "T-2022.03-SP2"
REQUIRED_SOURCE_FILES = (
    "sprocess_fps.cmd",
    "IdVg_des.cmd",
    "IdVd_des.cmd",
    "BVdss_des.cmd",
    "sdevice.par",
    "gtree.dat",
)
STAGE_ORDER = ("sprocess", "idvg", "idvd", "bv")
STAGE_COMMANDS = {
    "sprocess": ("sprocess", "sprocess.cmd"),
    "idvg": ("sdevice", "IdVg.cmd"),
    "idvd": ("sdevice", "IdVd.cmd"),
    "bv": ("sdevice", "BVdss.cmd"),
}
NODE_VARIABLES = {
    "IdVg_des.cmd": {
        "node": "2", "tdr": "n1_fps.tdr", "parameter": "sdevice.par",
        "log": "n2_des.log", "plot": "n2_des.plt", "tdrdat": "n2_des.tdr",
    },
    "IdVd_des.cmd": {
        "node": "4", "tdr": "n1_fps.tdr", "parameter": "sdevice.par",
        "log": "n4_des.log", "plot": "n4_des.plt", "tdrdat": "n4_des.tdr",
    },
    "BVdss_des.cmd": {
        "node": "6", "tdr": "n1_fps.tdr", "parameter": "sdevice.par",
        "log": "n6_des.log", "plot": "n6_des.plt", "tdrdat": "n6_des.tdr",
    },
}
BUNDLE_NAMES = {
    "sprocess_fps.cmd": "sprocess.cmd",
    "IdVg_des.cmd": "IdVg.cmd",
    "IdVd_des.cmd": "IdVd.cmd",
    "BVdss_des.cmd": "BVdss.cmd",
    "sdevice.par": "sdevice.par",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def executable(name: str) -> str:
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
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def records(root: Path) -> list[dict[str, Any]]:
    if not root.exists():
        return []
    result = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        result.append({
            "path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    return result


def write_json(path: Path, payload: dict[str, Any], *, validate: bool = False) -> None:
    if validate:
        validate_document(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def manifest(manifest_type: str,
             run_id: str,
             file_records: list[dict[str, Any]],
             metadata: dict[str, Any]) -> dict[str, Any]:
    payload = {
        "schema": "vela.templates_ldmos.phase01_manifest.v1",
        "manifest_type": manifest_type,
        "benchmark": BENCHMARK,
        "run_id": run_id,
        "generated_at": utc_now(),
        "metadata": metadata,
        "records": file_records,
    }
    validate_document(payload)
    return payload


def preprocess(text: str, variables: dict[str, str]) -> str:
    for name, value in variables.items():
        text = text.replace(f"@{name}@", value)
    return text


def prepare_bundle(source: Path, bundle: Path) -> dict[str, Any]:
    missing = [name for name in REQUIRED_SOURCE_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing Templates/LDMOS source files: {', '.join(missing)}")
    bundle.mkdir(parents=True, exist_ok=True)
    materialization: list[dict[str, Any]] = []
    for source_name, bundle_name in BUNDLE_NAMES.items():
        source_path = source / source_name
        text = source_path.read_text(encoding="utf-8", errors="strict")
        variables = {"node": "1"} if source_name == "sprocess_fps.cmd" else NODE_VARIABLES.get(source_name, {})
        rendered = preprocess(text, variables)
        unresolved = sorted(set(re.findall(r"@[A-Za-z0-9_:+.-]+@", rendered)))
        if unresolved:
            raise ValueError(f"{source_name}: unresolved Workbench tokens: {unresolved}")
        output = bundle / bundle_name
        output.write_text(rendered, encoding="utf-8")
        materialization.append({
            "source": source_name,
            "bundle": bundle_name,
            "replacements": variables,
            "source_sha256": sha256_file(source_path),
            "bundle_sha256": sha256_file(output),
        })
    report = {
        "schema": "vela.templates_ldmos.materialization.v1",
        "policy": "Only Workbench path/node/output tokens are replaced.",
        "files": materialization,
    }
    write_json(bundle / "materialization.json", report)
    return report


def run(argv: Sequence[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv), check=check, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(destination, filter="data")


def fetch_applications_source(args: argparse.Namespace,
                              remote_dir: str,
                              run_dir: Path) -> tuple[Path, str]:
    original = run_dir / "sentaurus_original"
    archive = run_dir / "manifest" / "templates_ldmos_source.tgz"
    quoted_remote = shlex.quote(remote_dir)
    quoted_source = shlex.quote(args.applications_path)
    create = (
        f"umask 077; if test -e {quoted_remote}; then echo 'remote run exists' >&2; exit 17; fi; "
        f"mkdir -p {quoted_remote}/work; "
        f"tar -C {quoted_source} -czf {quoted_remote}/templates_ldmos_source.tgz ."
    )
    run([args.ssh_bin, args.ssh_target, create])
    archive.parent.mkdir(parents=True, exist_ok=True)
    run([
        args.scp_bin,
        f"{args.ssh_target}:{remote_dir}/templates_ldmos_source.tgz",
        str(archive),
    ])
    safe_extract(archive, original)
    return original, sha256_file(archive)


def copy_local_source(source_dir: Path, run_dir: Path) -> Path:
    original = run_dir / "sentaurus_original"
    shutil.copytree(source_dir, original)
    return original


def capture_host_metadata(args: argparse.Namespace) -> dict[str, Any]:
    probe = (
        "printf 'hostname='; hostname; printf 'user='; whoami; printf 'pwd='; pwd; "
        "printf 'uname='; uname -a; printf 'sprocess='; command -v sprocess; "
        "printf 'sdevice='; command -v sdevice; printf 'svisual='; command -v svisual; "
        "printf 'tdx='; command -v tdx; "
        "printf 'license='; printf '%s\\n' \"$SNPSLMD_LICENSE_FILE\"; "
        "printf 'stdb='; printf '%s\\n' \"$STDB\""
    )
    environment = run([args.ssh_bin, args.ssh_target, probe]).stdout
    device_banner = run([
        args.ssh_bin, args.ssh_target,
        "timeout 15s sdevice -h 2>&1 | sed -n '1,6p'",
    ]).stdout.strip()
    process_banner = run([
        args.ssh_bin, args.ssh_target,
        "timeout 15s sprocess -v 2>&1 | sed -n '1,6p'",
    ]).stdout.strip()
    visual_banner = run([
        args.ssh_bin, args.ssh_target,
        "timeout 15s svisual -h 2>&1 | sed -n '1,6p'",
    ]).stdout.strip()
    banners = {
        "sprocess": process_banner,
        "sdevice": device_banner,
        "svisual": visual_banner,
    }
    missing = [name for name, banner in banners.items() if args.sentaurus_version not in banner]
    if missing:
        raise RuntimeError(
            f"Sentaurus banner(s) do not contain {args.sentaurus_version!r}: {missing}"
        )
    values: dict[str, str] = {}
    for line in environment.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return {"environment": values, "banners": banners}


def upload_bundle(args: argparse.Namespace, bundle: Path, remote_dir: str) -> None:
    files = [str(path) for path in sorted(bundle.iterdir()) if path.is_file()]
    run([args.scp_bin, *files, f"{args.ssh_target}:{remote_dir}/work/"])


def execute_stages(args: argparse.Namespace,
                   remote_dir: str,
                   stages: Sequence[str]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    process_ok = True
    for stage in stages:
        if stage != "sprocess" and not process_ok:
            results.append({"stage": stage, "status": "skipped", "reason": "sprocess_failed"})
            continue
        executable_name, deck = STAGE_COMMANDS[stage]
        work = shlex.quote(f"{remote_dir}/work")
        command = (
            f"cd {work}; /usr/bin/time -v -o timing_{stage}.txt "
            f"{executable_name} {shlex.quote(deck)} > run_{stage}.out 2>&1; "
            f"rc=$?; printf '%s\\n' \"$rc\" > {stage}.exitcode; exit $rc"
        )
        started = time.perf_counter()
        completed = run([args.ssh_bin, args.ssh_target, command], check=False)
        elapsed = time.perf_counter() - started
        status = "pass" if completed.returncode == 0 else "fail"
        results.append({
            "stage": stage,
            "status": status,
            "exit_code": completed.returncode,
            "wall_clock_seconds_host": elapsed,
            "command": command,
            "ssh_output": completed.stdout.strip(),
        })
        if stage == "sprocess":
            process_ok = completed.returncode == 0
    archive_command = (
        f"cd {shlex.quote(f'{remote_dir}/work')}; "
        "find . -maxdepth 1 -type f ! -name templates_ldmos_results.tgz -print0 | "
        "tar --null -czf templates_ldmos_results.tgz --files-from=-"
    )
    run([args.ssh_bin, args.ssh_target, archive_command])
    return results


def fetch_results(args: argparse.Namespace, remote_dir: str, run_dir: Path) -> Path:
    archive = run_dir / "raw" / "templates_ldmos_results.tgz"
    archive.parent.mkdir(parents=True, exist_ok=True)
    run([
        args.scp_bin,
        f"{args.ssh_target}:{remote_dir}/work/templates_ldmos_results.tgz",
        str(archive),
    ])
    safe_extract(archive, archive.parent)
    return archive


def normalize_plt_files(raw: Path, normalized: Path) -> list[dict[str, Any]]:
    from sentaurus_import import parse_quoted_list, parse_values_block

    normalized.mkdir(parents=True, exist_ok=True)
    outputs: list[dict[str, Any]] = []
    for path in sorted(raw.glob("*.plt")):
        text = path.read_text(encoding="utf-8", errors="ignore")
        datasets = parse_quoted_list(text, "datasets")
        if not datasets:
            continue
        rows = parse_values_block(text, len(datasets))
        output = normalized / f"{path.stem}.csv"
        with output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(datasets)
            writer.writerows(rows)
        lower = path.name.lower()
        if "idvg" in lower or lower.startswith("n2"):
            stage, bias_candidates = "idvg", ("gate OuterVoltage", "gate Voltage")
        elif "idvd" in lower or lower.startswith("n4"):
            stage, bias_candidates = "idvd", ("drain InnerVoltage", "drain OuterVoltage", "drain Voltage")
        elif "bv" in lower or lower.startswith("n6"):
            stage, bias_candidates = "bv", ("drain InnerVoltage", "drain OuterVoltage", "drain Voltage")
        else:
            stage, bias_candidates = "unknown", ()
        bias_name = next((name for name in bias_candidates if name in datasets), None)
        current_name = "drain TotalCurrent" if "drain TotalCurrent" in datasets else None
        curve_output = None
        if bias_name and current_name:
            curve_output = normalized / f"{path.stem}_drain_curve.csv"
            bias_index, current_index = datasets.index(bias_name), datasets.index(current_name)
            with curve_output.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["bias_V", "current_total_A_per_um"])
                writer.writerows((row[bias_index], row[current_index]) for row in rows)
        outputs.append({
            "source": path.name,
            "stage": stage,
            "datasets": datasets,
            "row_count": len(rows),
            "normalized_csv": output.name,
            "drain_curve_csv": curve_output.name if curve_output else None,
            "bias_dataset": bias_name,
            "current_dataset": current_name,
        })
    write_json(normalized / "plt_normalization_manifest.json", {
        "schema": "vela.templates_ldmos.plt_normalization.v1",
        "current_unit": "A/um",
        "current_sign": "Sentaurus native drain TotalCurrent",
        "files": outputs,
    })
    return outputs


def git_commit() -> str:
    return run([
        r"D:\msys64\usr\bin\git.exe" if os.name == "nt" else "git",
        "-C", str(REPO), "rev-parse", "HEAD",
    ]).stdout.strip()


def parse_stages(raw: str) -> list[str]:
    requested = [item.strip().lower() for item in raw.split(",") if item.strip()]
    invalid = [item for item in requested if item not in STAGE_ORDER]
    if invalid or not requested:
        raise ValueError(f"invalid --stages: {invalid or requested}")
    return [stage for stage in STAGE_ORDER if stage in requested]


def build_summary(live: bool,
                  stage_results: Sequence[dict[str, Any]],
                  manifest_paths: Sequence[Path],
                  limitations: list[str]) -> dict[str, Any]:
    stage_gates = []
    for item in stage_results:
        stage_gates.append({
            "id": f"sentaurus_{item['stage']}",
            "status": item["status"] if item["status"] in {"pass", "fail"} else "not_run",
            "summary": (
                f"exit_code={item.get('exit_code')}"
                if "exit_code" in item else str(item.get("reason", "not run"))
            ),
            "metrics": {"wall_clock_seconds": item.get("wall_clock_seconds_host")},
        })
    all_pass = bool(stage_results) and all(item["status"] == "pass" for item in stage_results)
    status = "pass" if live and all_pass else ("fail" if live else "not_run")
    evidence = [
        {"path": str(path), "sha256": sha256_file(path), "description": path.stem}
        for path in manifest_paths if path.is_file()
    ]
    return {
        "schema": "vela.templates_ldmos.validation_summary.v1",
        "benchmark": BENCHMARK,
        "generated_at": utc_now(),
        "status": status,
        "highest_level": "L0" if status == "pass" else "none",
        "gates": [
            {"id": "source_integrity", "status": "pass", "summary": "Required source files and SHA-256 records are present."},
            {"id": "materialization", "status": "pass", "summary": "All Workbench tokens were resolved by declared replacements."},
            *stage_gates,
            {"id": "exact_mesh_import", "status": "not_run", "summary": "Stage 1 has not been evaluated by this runner."},
        ],
        "evidence": evidence,
        "limitations": limitations,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--staging-root", type=Path, default=DEFAULT_STAGING)
    parser.add_argument("--applications-path", default=DEFAULT_APPLICATIONS_PATH)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--run-id")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--sentaurus-version", default=EXPECTED_VERSION)
    parser.add_argument("--stages", default=",".join(STAGE_ORDER))
    parser.add_argument("--live", action="store_true")
    args = parser.parse_args(argv)

    stages = parse_stages(args.stages)
    if any(stage != "sprocess" for stage in stages) and "sprocess" not in stages:
        raise ValueError("device stages require sprocess in the same isolated run")
    run_id = args.run_id or datetime.now().strftime("templates_ldmos_%Y%m%d_%H%M%S")
    run_dir = args.staging_root.resolve() / run_id
    if run_dir.exists():
        raise FileExistsError(f"run directory already exists: {run_dir}")
    (run_dir / "manifest").mkdir(parents=True)
    remote_dir = f"{args.remote_root.rstrip('/')}/{run_id}"
    host_metadata: dict[str, Any] = {}
    source_archive_sha256 = None

    if args.live:
        host_metadata = capture_host_metadata(args)
        source, source_archive_sha256 = fetch_applications_source(args, remote_dir, run_dir)
    else:
        if args.source_dir is None:
            raise ValueError("--source-dir is required without --live")
        source = copy_local_source(args.source_dir.resolve(), run_dir)

    bundle = run_dir / "bundle"
    materialization = prepare_bundle(source, bundle)
    source_manifest_path = run_dir / "manifest" / "source_manifest.json"
    write_json(source_manifest_path, manifest("source", run_id, records(source), {
        "applications_path": args.applications_path,
        "source_archive_sha256": source_archive_sha256,
        "required_files": list(REQUIRED_SOURCE_FILES),
    }), validate=True)

    stage_results: list[dict[str, Any]] = []
    normalization: list[dict[str, Any]] = []
    if args.live:
        upload_bundle(args, bundle, remote_dir)
        stage_results = execute_stages(args, remote_dir, stages)
        fetch_results(args, remote_dir, run_dir)
        normalization = normalize_plt_files(run_dir / "raw", run_dir / "normalized")
    else:
        stage_results = [{"stage": stage, "status": "skipped", "reason": "dry_run"} for stage in stages]

    run_manifest_path = run_dir / "manifest" / "run_manifest.json"
    write_json(run_manifest_path, manifest("run", run_id, records(bundle), {
        "live": args.live,
        "vela_commit": git_commit(),
        "ssh_target": args.ssh_target,
        "remote_dir": remote_dir,
        "expected_sentaurus_version": args.sentaurus_version,
        "host": host_metadata,
        "stages": stage_results,
        "materialization": materialization,
    }), validate=True)
    artifact_manifest_path = run_dir / "manifest" / "artifact_manifest.json"
    artifact_records = records(run_dir / "raw") + [
        {**item, "path": f"normalized/{item['path']}"}
        for item in records(run_dir / "normalized")
    ]
    write_json(artifact_manifest_path, manifest("artifact", run_id, artifact_records, {
        "normalization": normalization,
        "raw_root": "raw",
        "normalized_root": "normalized",
    }), validate=True)

    source_manifest_hash = sha256_file(source_manifest_path)
    governance = draft_governance_contracts(git_commit(), source_manifest_hash)
    for name, payload in governance.items():
        write_json(run_dir / "manifest" / name, payload, validate=True)

    limitations = [] if args.live else ["Dry run only; no VM execution or oracle result exists."]
    summary = build_summary(
        args.live, stage_results,
        [source_manifest_path, run_manifest_path, artifact_manifest_path],
        limitations,
    )
    summary_path = run_dir / "reports" / "validation_summary.json"
    write_json(summary_path, summary, validate=True)
    report_path = run_dir / "reports" / "validation_summary.md"
    report_path.write_text(render_summary(summary), encoding="utf-8")
    print(json.dumps({
        "run_dir": str(run_dir),
        "remote_dir": remote_dir,
        "live": args.live,
        "status": summary["status"],
        "stage_results": stage_results,
    }, indent=2))
    return 0 if not args.live or summary["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
