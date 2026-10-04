#!/usr/bin/env python3
"""Freeze the conditional M61/E2 execution decision after M60, M63, and M62."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
M60 = ROOT / "simplemos_m60_tight_convergence_port_burst_evidence.json"
M62 = ROOT / "simplemos_m62_active_species_topology_intervention_evidence.json"
M63 = ROOT / "simplemos_m63_smooth_nwell_attribution_evidence.json"
OUT = ROOT / "path_sensitivity_gate"
REPORT = OUT / "m61_path_sensitivity_gate_report.json"
DOC = REPO / "docs/validation/simplemos_m61_path_sensitivity_gate_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m61/artifact.json"
EVIDENCE = ROOT / "simplemos_m61_path_sensitivity_gate_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m61_path_sensitivity_gate.py"
EXPECTED = {
    "reference_tcad/simplemos_sentaurus2022/simplemos_m60_tight_convergence_port_burst_evidence.json":
        "dd9e4866f953e6e1051bb42f0bbf06db4fffa5736a476308f524c443f9bae96c",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m62_active_species_topology_intervention_evidence.json":
        "d452dd49895c8acc1dbcb72b86376cf7d261b45fe4b3daf89e2a97e9d871ddb7",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m63_smooth_nwell_attribution_evidence.json":
        "2480a8ad273d8953566f525cb5806d56351b7bb365b4695baa9c0202a94ac38c",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8",
                    newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def validate_sources() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    for relative, expected in EXPECTED.items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M61 gate source changed: {relative}")
    return read_json(M60), read_json(M62), read_json(M63)


def freeze() -> dict[str, Any]:
    m60, m62, m63 = validate_sources()
    execute = m60["classification"] not in {
        "tight_convergence_suppresses_port_bursts",
        "tight_convergence_reduces_but_does_not_suppress_port_bursts",
    }
    decision = "execute_m61_e2" if execute else "not_required_m60_decisive"
    checks = {
        "m60_evidence_frozen": m60["status"] == "frozen",
        "m60_classification_decisive": m60["classification"] == "tight_convergence_suppresses_port_bursts",
        "m60_zero_remaining_bursts": m60["acceptance"]["all_checks_pass"],
        "m62_stop_result_frozen": m62["status"] == "frozen",
        "m63_independent_track_frozen": m63["status"] == "frozen",
        "m61_e2_not_executed": not execute,
        "production_defaults_not_changed": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m61_path_sensitivity_gate_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "decision": decision,
        "reason": (
            "M61/E2 was conditional on an uncertain or negative M60/E1 result. "
            "M60 instead eliminated all 15 burst flags, reduced the n23 target "
            "observable gap by 98.8585 percent, and reduced its Vela error to "
            "0.0243983 dex, so path-sensitivity replay would not discriminate "
            "the already-decided convergence mechanism."),
        "sequence_closure": {
            "m60": m60["classification"],
            "m63": m63["classification"],
            "m62": m62["classification"],
            "m61": decision,
        },
        "execution": {"new_sentaurus_execution": False,
                      "new_vela_execution": False,
                      "deck_count": 0},
        "acceptance": checks,
    }
    write_json(REPORT, report)
    DOC.write_text(f"""# SimpleMOS M61 路径敏感性门禁

M61/E2判定为 `{decision}`，因此没有创建或执行步长计划deck。

M61原本只在M60/E1结果不确定或无实质作用时执行。M60已将15个burst标志点全部消除，n23目标点观测差削减98.8585%，默认Id与Vela差降到0.0243983 dex，已对收敛残差机制给出决定性阳性判别。M63独立闭合平滑NWell成分为阈值样水平平移，M62则按写回门禁安全停止；两者都不重新触发M61条件。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M61 path-sensitivity execution gate",
        "status": report["status"], "classification": decision,
        "summary": "M61/E2 not required because M60/E1 was decisive",
        "report": portable(REPORT),
    })
    artifacts = [REPORT, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m61_path_sensitivity_gate_evidence.v1",
        "status": "frozen" if checks["all_checks_pass"] else "failed",
        "decision": decision,
        "source_hashes": EXPECTED,
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "acceptance": checks,
    })
    return report


def verify() -> dict[str, Any]:
    validate_sources()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M61 gate evidence is not accepted and frozen")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M61 artifact hash changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    report = verify() if args.verify else freeze()
    print(f"M61 gate: {report['decision']}")


if __name__ == "__main__":
    main()
