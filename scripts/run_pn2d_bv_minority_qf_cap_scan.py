#!/usr/bin/env python
"""Task 3 validation: minority quasi-Fermi update-cap BV scan.

Clones the PN2D Sentaurus2018 BV deck (default half-box avalanche source
support; `source_volume_factor` removed) and varies only the Newton
quasi-Fermi update caps:

  * solver.quasi_fermi_update_limit_V            (global / majority cap)
  * solver.quasi_fermi_update_limit_minority_V   (tighter minority-only cap)

For each candidate it runs the release runner and records the deepest reverse
bias reached, the failure reason (if any), the max accepted p95 density jump,
and the Sentaurus knee markers / max log10 current error over the -20..-10 V
window. This isolates whether a per-node minority cap unblocks the low-bias
gate and/or preserves the knee window without regression.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_knee_shape as knee


REPO = Path(__file__).resolve().parents[1]
DEFAULT_BASE = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reports"
    / "qflim0p05_source_factor_scan"
    / "factor_0p875.json"
)
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
DEFAULT_REFERENCE = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reference_curves"
    / "pn2d_sentaurus2018_bv_reference.csv"
)


def build_deck(base: dict[str, Any], label: str, global_v: float,
               minority_v: float, stop_v: float, out_dir: Path) -> Path:
    deck = json.loads(json.dumps(base))  # deep copy
    solver = deck.setdefault("solver", {})
    solver["quasi_fermi_update_limit_V"] = global_v
    solver["quasi_fermi_update_limit_minority_V"] = minority_v
    impact = solver.get("impact_ionization")
    if isinstance(impact, dict):
        impact.pop("source_volume_factor", None)  # default half-box support
    sweep = deck.setdefault("sweep", {})
    sweep["stop"] = stop_v
    sweep["write_vtk"] = False
    sweep["write_state_file"] = f"{label}_last_state.csv"
    deck["output_csv"] = f"{label}.csv"
    config_path = out_dir / f"{label}.json"
    config_path.write_text(json.dumps(deck, indent=2) + "\n", encoding="utf-8")
    return config_path


def run_candidate(runner: Path, config_path: Path, out_dir: Path) -> dict[str, Any]:
    proc = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=str(out_dir),
        capture_output=True,
        text=True,
    )
    return {
        "exit_code": proc.returncode,
        "stdout_tail": proc.stdout.strip().splitlines()[-1:] or [""],
        "stderr_tail": proc.stderr.strip().splitlines()[-3:],
    }


def _float(value: str | None) -> float | None:
    return knee._float_or_none(value)


def parse_curve(csv_path: Path, reference: list[tuple[float, float]],
                bias_min: float, bias_max: float) -> dict[str, Any]:
    rows: list[dict[str, str]] = []
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))

    converged = [r for r in rows
                 if str(r.get("converged", "")).strip() in {"1", "true", "True"}]
    biases = [_float(r.get("bias_V")) for r in converged]
    biases = [b for b in biases if b is not None]
    deepest = min(biases) if biases else None
    reached_minus20 = deepest is not None and deepest <= -20.0 + 1.0e-2

    failure_reason = ""
    failed_bias = ""
    max_p95 = 0.0
    for r in rows:
        fr = (r.get("failure_reason") or "").strip()
        br = (r.get("branch_acceptance_reason") or "").strip()
        if fr and not failure_reason:
            failure_reason = fr
            failed_bias = (r.get("bias_V") or "").strip()
        if br and not failure_reason:
            failure_reason = br
            failed_bias = (r.get("bias_V") or "").strip()
        p95 = _float(r.get("electron_density_jump_p95_abs_dex"))
        if p95 is not None and str(r.get("converged", "")).strip() in {"1", "true", "True"}:
            max_p95 = max(max_p95, abs(p95))

    points = knee.load_curve(csv_path) if csv_path.exists() else []
    summary = None
    max_error = None
    if points:
        try:
            summary = knee.summarize_curve(points, bias_min, bias_max)
            max_error = knee.max_abs_log10_error(points, reference, bias_min, bias_max)
        except ValueError:
            # No overlap with the knee window (e.g. a low-bias-only sweep).
            summary = None
            max_error = None

    return {
        "rows": len(rows),
        "converged_rows": len(converged),
        "deepest_reverse_bias_V": -deepest if deepest is not None else None,
        "reached_minus20V": reached_minus20,
        "failure_reason": failure_reason,
        "failure_bias_V": failed_bias,
        "max_accepted_p95_jump": max_p95,
        "first_growth_over_1p5_V": summary.first_growth_over_1p5 if summary else None,
        "first_growth_over_2p0_V": summary.first_growth_over_2p0 if summary else None,
        "max_abs_log10_current_error_decades": max_error,
    }


def main() -> int:
    args = parse_args()
    args.runner = args.runner.resolve()
    args.base = args.base.resolve()
    args.reference = args.reference.resolve()
    base = json.loads(args.base.read_text(encoding="utf-8"))
    reference = knee.load_curve(args.reference)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir = args.out_dir.resolve()

    candidates = parse_candidates(args.candidate)
    results: list[dict[str, Any]] = []
    for label, global_v, minority_v, stop_v in candidates:
        config_path = build_deck(base, label, global_v, minority_v, stop_v, args.out_dir)
        run_info = run_candidate(args.runner, config_path, args.out_dir)
        curve = parse_curve(args.out_dir / f"{label}.csv", reference,
                            args.bias_min, args.bias_max)
        row = {
            "label": label,
            "global_qf_limit_V": global_v,
            "minority_qf_limit_V": minority_v,
            "sweep_stop_V": stop_v,
            **curve,
            "exit_code": run_info["exit_code"],
        }
        results.append(row)
        print(json.dumps(row, sort_keys=True))

    columns = [
        "label", "global_qf_limit_V", "minority_qf_limit_V", "sweep_stop_V",
        "converged_rows", "deepest_reverse_bias_V", "reached_minus20V",
        "failure_reason", "failure_bias_V", "max_accepted_p95_jump",
        "first_growth_over_1p5_V", "first_growth_over_2p0_V",
        "max_abs_log10_current_error_decades", "exit_code",
    ]
    summary_csv = args.out_dir / "minority_qf_cap_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in results:
            writer.writerow({c: row.get(c, "") for c in columns})
    (args.out_dir / "minority_qf_cap_summary.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {summary_csv}")
    return 0


def parse_candidates(entries: list[str]) -> list[tuple[str, float, float, float]]:
    parsed: list[tuple[str, float, float, float]] = []
    for entry in entries:
        parts = entry.split(":")
        if len(parts) != 4:
            raise SystemExit(
                f"invalid --candidate entry (label:global:minority:stop): {entry}")
        parsed.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-10.0)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018"
        / "reports" / "minority_qf_cap_scan",
    )
    parser.add_argument(
        "--candidate",
        action="append",
        required=True,
        metavar="LABEL:GLOBAL:MINORITY:STOP",
        help="Repeatable 'label:global_V:minority_V:stop_V' candidate entry.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main())
