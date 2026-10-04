#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M12 terminal-sensitivity evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022/terminal_sensitivity"
          / "m12_terminal_sensitivity_report.json")
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m12_terminal_sensitivity_evidence.json")
FIGURE_DIR = REPO / "docs/validation/figures/simplemos_m12"
FIGURES = [
    "simplemos_m12_residual_spectrum.png",
    "simplemos_m12_h1_contact_ratio.png",
    "simplemos_m12_terminal_shift.png",
    "simplemos_m12_contact_conditioning.png",
]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m12_terminal_sensitivity_contract_v1.json",
    "scripts/run_simplemos_m12_terminal_sensitivity.py",
    "scripts/plot_simplemos_m12_terminal_sensitivity.py",
    "scripts/freeze_simplemos_m12_terminal_sensitivity_evidence.py",
    "tests/regression/test_simplemos_m12_terminal_sensitivity.py",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--output", type=Path, default=EVIDENCE)
    args = parser.parse_args()
    report_path = args.report.resolve()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report["status"] != "complete":
        raise ValueError("M12 report is not complete")
    if report["execution"]["state_count"] != 16:
        raise ValueError("M12 report does not contain 16 states")
    if report["findings"]["drain_replay_identity_pass_count"] != 16:
        raise ValueError("M12 drain-cut identity did not pass in all states")
    artifacts = dict(report["artifacts"])
    artifacts["report"] = {"path": portable(report_path),
                           "sha256": sha256(report_path)}
    figures = [{"path": portable(FIGURE_DIR / name),
                "sha256": sha256(FIGURE_DIR / name)} for name in FIGURES]
    implementation = {relative: sha256(REPO / relative)
                      for relative in IMPLEMENTATION}
    evidence = {
        "schema": "vela.simplemos.sdevice.m12_terminal_sensitivity_evidence.v1",
        "status": "complete",
        "scope": report["scope"],
        "sentaurus_release": "T-2022.03-SP2",
        "execution": report["execution"],
        "findings": report["findings"],
        "conclusions": {
            "h1": (
                "accepted" if report["findings"]["h1_all_states_pass"]
                else "rejected"),
            "supported": [
                "The drain-cut delta-current reconstruction exactly matches the M10 frozen replay in all 16 states.",
                "The source/drain contact-cut mobility response is far smaller than the all-active-edge P95 improvement.",
                "The M10 mobility improvement and frozen terminal-current response are spatially disconnected.",
                "The high-NWell low-drain residual remains a smooth, same-sign weak-current bias."
            ],
            "not_supported": [
                "Attributing the weak-current Id mismatch directly to the improved internal-edge HFS mobility residual.",
                "Attributing the deep-off terminal mismatch to large cancellation within the drain contact cut.",
                "Treating the frozen mapped-state contact sum as a self-consistent SRH continuity closure."
            ]
        },
        "continuity_scope": {
            "frozen_state": "contact-flux residual only",
            "self_consistent_context": "M8 SRH/KCL rows kept separate"
        },
        "source_integrity": report["inputs"],
        "artifacts": artifacts,
        "figures": figures,
        "implementation_sha256": implementation,
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"status": "frozen", "output": portable(output),
                      "h1": evidence["conclusions"]["h1"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
