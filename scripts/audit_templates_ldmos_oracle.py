#!/usr/bin/env python3
"""Audit a sealed Templates/LDMOS live run without changing its raw evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

from run_templates_ldmos_sentaurus_vm import (
    EXPECTED_VERSION,
    capture_host_metadata,
    executable,
    records,
    sha256_file,
    utc_now,
    write_json,
)


# Sentaurus convergence logs routinely emit the success phrase
# ``Error smaller than 1 (...)``.  Only diagnostic headers with an explicit
# colon are treated as strong errors here; process exit codes remain a
# separate fail-closed gate.
STRONG_ERROR = re.compile(r"^\s*(?:error|fatal)\s*:", re.IGNORECASE)


def curve_metrics(path: Path) -> dict[str, Any]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    points = [
        (float(row["bias_V"]), float(row["current_total_A_per_um"]))
        for row in rows
    ]
    finite = all(math.isfinite(value) for point in points for value in point)
    monotonic = all(b >= a for (a, _), (b, _) in zip(points, points[1:]))
    duplicate_biases = sum(b == a for (a, _), (b, _) in zip(points, points[1:]))
    backtracks = sum(b < a for (a, _), (b, _) in zip(points, points[1:]))
    return {
        "path": path.name,
        "sha256": sha256_file(path),
        "point_count": len(points),
        "finite": finite,
        "bias_nondecreasing": monotonic,
        "duplicate_adjacent_bias_count": duplicate_biases,
        "bias_backtrack_count": backtracks,
        "bias_min_V": min((point[0] for point in points), default=None),
        "bias_max_V": max((point[0] for point in points), default=None),
        "max_abs_current_A_per_um": max((abs(point[1]) for point in points), default=None),
    }


def strong_log_errors(path: Path) -> list[dict[str, Any]]:
    result = []
    for number, line in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        if STRONG_ERROR.search(line):
            result.append({"line": number, "text": line.strip()[:500]})
    return result


def classify_tdr(name: str) -> str:
    lower = name.lower()
    if lower == "n1_fps.tdr":
        return "final_sprocess_structure"
    if "state_" in lower:
        return "derived_representative_state"
    if lower.startswith("n2"):
        return "idvg_final"
    if lower.startswith("n4_vg1"):
        return "idvd_vg4_saved_or_final"
    if lower.startswith("n4_vg2") or lower.startswith("n4"):
        return "idvd_saved_or_final"
    if lower.startswith("n6"):
        return "bv_final_or_maxcurrent"
    return "intermediate_process_or_unknown"


def parse_deck_contract(bundle: Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for name in ("IdVg.cmd", "IdVd.cmd", "BVdss.cmd"):
        path = bundle / name
        text = path.read_text(encoding="utf-8")
        entry: dict[str, Any] = {
            "sha256": sha256_file(path),
            "currentplot_interval_counts": [
                int(value) for value in re.findall(r"CurrentPlot\s*\([^\n]*Intervals\s*=\s*([0-9]+)", text)
            ],
            "couples_temperature": bool(re.search(r"Coupled\s*\{[^}]*\bTemperature\b", text)),
            "uses_okuto": "Avalanche(Okuto)" in text,
            "uses_thermode": "Thermode" in text,
        }
        if name == "BVdss.cmd":
            continuation = {}
            for key in (
                "Increment", "Decrement", "InitialVstep", "MinStep", "MaxVoltage",
                "MinVoltage", "MaxCurrent", "MinCurrent", "Iadapt",
            ):
                match = re.search(rf"\b{key}\s*=\s*([-+0-9.eE]+)", text, re.IGNORECASE)
                continuation[key] = float(match.group(1)) if match else None
            thermode_match = re.search(r"Thermode\s*\{(.*?)\}\s*File\s*\{", text, re.DOTALL | re.IGNORECASE)
            entry["continuation"] = continuation
            entry["thermal_contacts"] = (
                re.findall(r"\bName\s*=\s*\"([^\"]+)\"", thermode_match.group(1), re.IGNORECASE)
                if thermode_match else []
            )
        result[name] = entry
    return result


def audit(run_dir: Path) -> dict[str, Any]:
    raw = run_dir / "raw"
    normalized = run_dir / "normalized"
    curves = [curve_metrics(path) for path in sorted(normalized.glob("*_drain_curve.csv"))]
    logs = []
    log_paths = set(raw.glob("run_*.out")) | set(raw.glob("*.log"))
    for path in sorted(log_paths):
        errors = strong_log_errors(path)
        logs.append({
            "path": path.name,
            "sha256": sha256_file(path),
            "strong_error_count": len(errors),
            "strong_errors": errors,
        })
    exit_codes = {}
    for path in sorted(raw.glob("*.exitcode")):
        exit_codes[path.stem] = int(path.read_text(encoding="utf-8").strip())
    state_paths = list(raw.glob("*.tdr")) + list(
        (run_dir / "representative_states" / "raw").glob("*.tdr")
    )
    state_capture_manifest_path = run_dir / "representative_states" / "state_capture_manifest.json"
    state_capture_manifest = (
        json.loads(state_capture_manifest_path.read_text(encoding="utf-8"))
        if state_capture_manifest_path.is_file() else None
    )
    states = [
        {
            "path": path.relative_to(run_dir).as_posix(),
            "name": path.name,
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
            "classification": classify_tdr(path.name),
        }
        for path in sorted(state_paths)
    ]
    state_name_counts = {
        prefix: sum(path.name.lower().startswith(prefix) for path in state_paths)
        for prefix in (
            "state_idvg_eq_0v", "state_idvg_", "state_idvd_vg4_",
            "state_idvd_vg8_", "state_bv_path_",
        )
    }
    representative_counts_complete = (
        state_name_counts["state_idvg_eq_0v"] >= 1
        and state_name_counts["state_idvg_"] >= 4
        and state_name_counts["state_idvd_vg4_"] >= 6
        and state_name_counts["state_idvd_vg8_"] >= 6
        and state_name_counts["state_bv_path_"] >= 4
    )
    stage_curve_counts = {
        stage: sum(stage in curve["path"].lower() or prefix in curve["path"].lower() for curve in curves)
        for stage, prefix in (("idvg", "n2"), ("idvd", "n4"), ("bv", "n6"))
    }
    iv_curves = [
        item for item in curves
        if item["path"].lower().startswith(("idvg_", "idvd_"))
    ]
    bv_curves = [item for item in curves if item["path"].lower().startswith("n6_")]
    gates = {
        "all_exit_codes_zero": bool(exit_codes) and all(value == 0 for value in exit_codes.values()),
        "no_strong_log_errors": all(item["strong_error_count"] == 0 for item in logs),
        "three_curve_families_present": all(value > 0 for value in stage_curve_counts.values()),
        "curves_finite": bool(curves) and all(item["finite"] for item in curves),
        "iv_curve_bias_nondecreasing": bool(iv_curves) and all(
            item["bias_nondecreasing"] for item in iv_curves
        ),
        # The BV deck switches from voltage to current continuation at Iadapt;
        # voltage backtracking is a valid part of that ordered path.  Never
        # sort it before comparison.
        "bv_continuation_path_preserved": bool(bv_curves) and all(
            item["point_count"] > 1 and item["finite"] for item in bv_curves
        ),
        "final_sprocess_tdr_present": any(item["name"] == "n1_fps.tdr" for item in states),
        "representative_state_set_complete": any(
            item["classification"] == "derived_representative_state" for item in states
        ) and representative_counts_complete and state_capture_manifest is not None and all(
            item["status"] == "pass" for item in state_capture_manifest["stage_results"]
        ),
    }
    return {
        "schema": "vela.templates_ldmos.oracle_audit.v1",
        "generated_at": utc_now(),
        "run_id": run_dir.name,
        "exit_codes": exit_codes,
        "curve_family_counts": stage_curve_counts,
        "deck_contract": parse_deck_contract(run_dir / "bundle"),
        "curves": curves,
        "logs": logs,
        "tdr_states": states,
        "representative_state_prefix_counts": state_name_counts,
        "gates": gates,
        "limitations": [
            "Original and output-only derivative state-capture decks are classified separately.",
            "BV voltage backtracking after the Iadapt transition is retained in original row order and must never be silently sorted.",
        ],
    }


def render(audit_result: dict[str, Any]) -> str:
    lines = [
        "# Templates/LDMOS stage 0 oracle audit",
        "",
        f"- Run: `{audit_result['run_id']}`",
        "",
        "## Gates",
        "",
        "| Gate | Result |",
        "| --- | --- |",
    ]
    lines.extend(f"| `{name}` | `{value}` |" for name, value in audit_result["gates"].items())
    lines.extend([
        "",
        "## Curves",
        "",
        "| File | Points | Bias range (V) | Max current (A/um) | Finite | Monotonic |",
        "| --- | ---: | --- | ---: | --- | --- |",
    ])
    for curve in audit_result["curves"]:
        lines.append(
            f"| {curve['path']} | {curve['point_count']} | "
            f"{curve['bias_min_V']} to {curve['bias_max_V']} | "
            f"{curve['max_abs_current_A_per_um']} | {curve['finite']} | "
            f"{curve['bias_nondecreasing']} |"
        )
    lines.extend(["", "## State inventory", ""])
    lines.extend(f"- `{item['path']}`: {item['classification']}" for item in audit_result["tdr_states"])
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--ssh-target", help="Refresh and seal sprocess/sdevice/svisual banners.")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--sentaurus-version", default=EXPECTED_VERSION)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = audit(args.run_dir.resolve())
    if args.ssh_target:
        result["vm_metadata_refresh"] = capture_host_metadata(args)
    reports = args.run_dir.resolve() / "reports"
    write_json(reports / "oracle_audit.json", result)
    (reports / "oracle_audit.md").write_text(render(result), encoding="utf-8")
    write_json(args.run_dir.resolve() / "manifest" / "oracle_audit_artifacts.json", {
        "schema": "vela.templates_ldmos.oracle_audit_artifacts.v1",
        "records": records(reports),
    })
    print(json.dumps({"gates": result["gates"]}, indent=2))
    return 0 if all(result["gates"].values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
