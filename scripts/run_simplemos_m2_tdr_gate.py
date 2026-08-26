#!/usr/bin/env python3
"""Qualify the frozen SimpleMOS nominal TDR before any device workflow uses it."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = (
    REPO
    / "reference_tcad"
    / "simplemos_sentaurus2022"
    / "simplemos_sdevice_validation_contract_v1.json"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def gate_report(source: Path,
                contract_path: Path,
                expected_sha256: str,
                actual_sha256: str,
                importer_report: dict[str, Any] | None,
                importer_returncode: int | None) -> dict[str, Any]:
    hash_passed = actual_sha256.lower() == expected_sha256.lower()
    checks: list[dict[str, Any]] = [{
        "code": "source.sha256",
        "passed": hash_passed,
        "message": (
            "nominal TDR SHA-256 matches the frozen contract"
            if hash_passed
            else f"expected {expected_sha256}, observed {actual_sha256}"
        ),
    }]
    if importer_report is not None:
        checks.extend(importer_report.get("checks", []))
    elif hash_passed:
        checks.append({
            "code": "importer.execution",
            "passed": False,
            "message": f"sentaurus_import exited with code {importer_returncode}",
        })
    return {
        "schema": "vela.simplemos.sdevice.m2_tdr_gate_report.v1",
        "source": str(source),
        "contract": str(contract_path),
        "expected_sha256": expected_sha256,
        "actual_sha256": actual_sha256,
        "passed": bool(checks) and all(bool(check.get("passed")) for check in checks),
        "checks": checks,
        "importer_report": importer_report,
    }


def default_importer() -> Path:
    executable = "sentaurus_import.exe" if sys.platform.startswith("win") else "sentaurus_import"
    return REPO / "build" / executable


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdr", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--tdr-importer", type=Path, default=default_importer())
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    tdr = args.tdr.resolve()
    contract_path = args.contract.resolve()
    output_dir = args.output_dir.resolve()
    report_path = output_dir / "simplemos_m2_tdr_gate_report.json"
    core_report_path = output_dir / "tdr_qualification_report.json"

    if not tdr.is_file():
        raise FileNotFoundError(f"missing nominal TDR: {tdr}")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    expected_sha256 = str(contract["nominal_tdr"]["tdr_sha256"])
    actual_sha256 = sha256_file(tdr)
    if actual_sha256.lower() != expected_sha256.lower():
        report = gate_report(
            tdr, contract_path, expected_sha256, actual_sha256, None, None)
        write_json(report_path, report)
        print(json.dumps(report, indent=2))
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(args.tdr_importer.resolve()),
        "--tdr", str(tdr),
        "--inventory-json", str(output_dir / "tdr_inventory.json"),
        "--qualification-contract", str(contract_path),
        "--qualification-report", str(core_report_path),
        "--export-dir", str(output_dir),
    ]
    completed = subprocess.run(command, cwd=REPO, check=False)
    importer_report = (
        json.loads(core_report_path.read_text(encoding="utf-8"))
        if core_report_path.is_file()
        else None
    )
    report = gate_report(
        tdr,
        contract_path,
        expected_sha256,
        actual_sha256,
        importer_report,
        completed.returncode,
    )
    write_json(report_path, report)
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] and completed.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
