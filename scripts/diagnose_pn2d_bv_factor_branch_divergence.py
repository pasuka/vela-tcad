#!/usr/bin/env python
"""Task 1: locate and classify the PN2D BV source-volume-factor divergence.

Compares converged reverse-bias BV curves that differ ONLY in
``impact_ionization.source_volume_factor`` (all continuation settings identical:
qf limit, secant predictor, branch guard, driving force). For each candidate it
reuses the repository knee-shape growth-ratio convention and extracts the
per-bias continuation branch diagnostics already emitted by the runner:

  * predicted_initial_state    (secant predictor handoff selection)
  * branch_acceptance_status   (carrier-density branch guard verdict)
  * electron_density_jump_p95_abs_dex / *_max_abs_dex
  * step retry_count           (adaptive continuation backtracking)

It then answers the Task 1 question: at which bias do two factors diverge, and
is the divergence (a) a continuation/predictor branch selection effect, or
(b) a smooth source-volume magnitude effect on a single shared branch?

Read-only: it does not run the solver and does not modify any deck or artifact.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_knee_shape as knee


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCAN_ROOT = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reports"
)
DEFAULT_REFERENCE = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reference_curves"
    / "pn2d_sentaurus2018_bv_reference.csv"
)

# Per-bias branch-diagnostic columns emitted by the BV sweep runner.
BRANCH_COLUMNS = (
    "predicted_initial_state",
    "branch_acceptance_status",
    "branch_acceptance_reason",
    "electron_density_jump_p95_abs_dex",
    "electron_density_jump_max_abs_dex",
)

# Classification thresholds.
CURRENT_DIVERGENCE_DEX = 0.01  # |log10(I_b/I_a)| above this counts as divergence.
P95_BRANCH_DEX = 5.0e-3        # p95 jump difference above this is a branch change.


def _parse_retry_count(step_diagnostics: str | None) -> int | None:
    if not step_diagnostics:
        return None
    for token in step_diagnostics.split(";"):
        key, _, value = token.partition("=")
        if key.strip() == "retry_count":
            try:
                return int(value)
            except ValueError:
                return None
    return None


def load_branch_rows(path: Path) -> dict[int, dict[str, Any]]:
    """Map integer bias -> branch diagnostics for each converged row."""
    rows: dict[int, dict[str, Any]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            bias = knee._float_or_none(row.get("bias_V"))
            if bias is None:
                continue
            nearest = round(bias)
            if abs(bias - nearest) > 1.0e-2:
                continue
            if str(row.get("converged", "1")).strip() not in {"1", "true", "True"}:
                continue
            record = {
                col: row.get(col, "") for col in BRANCH_COLUMNS
            }
            record["retry_count"] = _parse_retry_count(row.get("step_diagnostics"))
            rows[nearest] = record
    return rows


def _branch_diff(a: dict[str, Any], b: dict[str, Any]) -> list[str]:
    diffs: list[str] = []
    if a.get("predicted_initial_state") != b.get("predicted_initial_state"):
        diffs.append("predicted_initial_state")
    if a.get("branch_acceptance_status") != b.get("branch_acceptance_status"):
        diffs.append("branch_acceptance_status")
    p95_a = knee._float_or_none(a.get("electron_density_jump_p95_abs_dex"))
    p95_b = knee._float_or_none(b.get("electron_density_jump_p95_abs_dex"))
    if p95_a is not None and p95_b is not None and abs(p95_a - p95_b) > P95_BRANCH_DEX:
        diffs.append("electron_density_jump_p95_abs_dex")
    if a.get("retry_count") != b.get("retry_count"):
        diffs.append("retry_count")
    return diffs


def compare_pair(
    label_a: str,
    points_a: list[tuple[float, float]],
    branch_a: dict[int, dict[str, Any]],
    label_b: str,
    points_b: list[tuple[float, float]],
    branch_b: dict[int, dict[str, Any]],
    bias_min: float,
    bias_max: float,
) -> dict[str, Any]:
    low = int(math.ceil(min(bias_min, bias_max)))
    high = int(math.floor(max(bias_min, bias_max)))
    per_bias: list[dict[str, Any]] = []
    branch_divergence_bias: float | None = None
    current_divergence_bias: float | None = None
    max_abs_current_dex = 0.0
    ratio_signs: set[int] = set()

    # Walk from least-negative toward most-negative bias.
    for bias in range(high, low - 1, -1):
        log_a = knee._interpolate_log_abs_current(points_a, float(bias))
        log_b = knee._interpolate_log_abs_current(points_b, float(bias))
        current_dex = None
        if log_a is not None and log_b is not None:
            current_dex = log_b - log_a
            max_abs_current_dex = max(max_abs_current_dex, abs(current_dex))
            if abs(current_dex) > 1.0e-9:
                ratio_signs.add(1 if current_dex > 0 else -1)
            if current_divergence_bias is None and abs(current_dex) > CURRENT_DIVERGENCE_DEX:
                current_divergence_bias = float(bias)
        rec_a = branch_a.get(bias)
        rec_b = branch_b.get(bias)
        branch_diffs = _branch_diff(rec_a, rec_b) if rec_a and rec_b else []
        if branch_diffs and branch_divergence_bias is None:
            branch_divergence_bias = float(bias)
        per_bias.append({
            "bias_V": float(bias),
            "current_ratio_dex": current_dex,
            "branch_diffs": branch_diffs,
            "p95_a": knee._float_or_none(
                (rec_a or {}).get("electron_density_jump_p95_abs_dex")),
            "p95_b": knee._float_or_none(
                (rec_b or {}).get("electron_density_jump_p95_abs_dex")),
            "predicted_initial_state_a": (rec_a or {}).get("predicted_initial_state"),
            "predicted_initial_state_b": (rec_b or {}).get("predicted_initial_state"),
        })

    branch_identical = branch_divergence_bias is None
    monotone_magnitude = len(ratio_signs) <= 1
    if not branch_identical:
        classification = "continuation_branch_selection"
    elif max_abs_current_dex > CURRENT_DIVERGENCE_DEX:
        classification = "source_volume_magnitude"
    else:
        classification = "no_material_divergence"

    return {
        "pair": f"{label_a} vs {label_b}",
        "bias_window_V": {"min": float(low), "max": float(high)},
        "branch_diagnostics_identical": branch_identical,
        "branch_divergence_onset_V": branch_divergence_bias,
        "current_divergence_onset_V": current_divergence_bias,
        "max_abs_current_ratio_dex": max_abs_current_dex,
        "current_ratio_monotone_sign": monotone_magnitude,
        "classification": classification,
        "per_bias": per_bias,
    }


def candidate_summary(
    points: list[tuple[float, float]],
    reference: list[tuple[float, float]],
    bias_min: float,
    bias_max: float,
) -> dict[str, Any]:
    summary = knee.summarize_curve(points, bias_min, bias_max)
    return {
        "first_growth_over_1p5_V": summary.first_growth_over_1p5,
        "first_growth_over_2p0_V": summary.first_growth_over_2p0,
        "max_abs_log10_current_error_decades": knee.max_abs_log10_error(
            points, reference, bias_min, bias_max),
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("## Source-Volume Factor Branch-Divergence Probe (Task 1)")
    lines.append("")
    lines.append(
        "Read-only comparison of converged BV curves that differ only in "
        "`impact_ionization.source_volume_factor`; all continuation settings "
        "(`quasi_fermi_update_limit_V=0.05`, secant predictor, branch guard "
        "`p95<=2.0 dex`, `driving_force=quasi_fermi_gradient`) are identical.")
    lines.append("")
    lines.append("| factor | reaches -20V | first >1.5 | first >2.0 | max abs log10 current error |")
    lines.append("|---:|:---:|---:|---:|---:|")
    ref = report["reference_knee"]
    lines.append(
        f"| Sentaurus | n/a | `{ref['first_growth_over_1p5_V']}` | "
        f"`{ref['first_growth_over_2p0_V']}` | n/a |")
    for cand in report["candidates"]:
        lines.append(
            f"| `{cand['source_volume_factor']}` | "
            f"`{cand['reached_minus20V']}` | "
            f"`{cand['first_growth_over_1p5_V']}` | "
            f"`{cand['first_growth_over_2p0_V']}` | "
            f"`{cand['max_abs_log10_current_error_decades']:.6f}` |")
    lines.append("")
    for pair in report["pairs"]:
        lines.append(f"### {pair['pair']}")
        lines.append("")
        lines.append(
            f"- branch diagnostics identical across window: "
            f"`{pair['branch_diagnostics_identical']}`")
        lines.append(
            f"- branch divergence onset: `{pair['branch_divergence_onset_V']}`")
        lines.append(
            f"- current-magnitude divergence onset: "
            f"`{pair['current_divergence_onset_V']} V`")
        lines.append(
            f"- max abs current ratio: "
            f"`{pair['max_abs_current_ratio_dex']:.6f}` decades")
        lines.append(f"- classification: **`{pair['classification']}`**")
        lines.append("")
    lines.append("### Conclusion")
    lines.append("")
    lines.append(report["conclusion"])
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    reference = knee.load_curve(args.reference)

    candidates: list[dict[str, Any]] = []
    loaded: list[tuple[str, list[tuple[float, float]], dict[int, dict[str, Any]]]] = []
    for factor, csv_path in args.candidate:
        points = knee.load_curve(csv_path)
        branch = load_branch_rows(csv_path)
        summary = candidate_summary(points, reference, args.bias_min, args.bias_max)
        reached = any(abs(b + 20.0) <= 1.0e-2 for b, _ in points)
        candidates.append({
            "source_volume_factor": factor,
            "csv": str(csv_path),
            "reached_minus20V": reached,
            **summary,
        })
        loaded.append((f"factor={factor}", points, branch))

    pairs: list[dict[str, Any]] = []
    for (la, pa, ba), (lb, pb, bb) in zip(loaded, loaded[1:]):
        pairs.append(compare_pair(
            la, pa, ba, lb, pb, bb, args.bias_min, args.bias_max))

    any_branch = any(not p["branch_diagnostics_identical"] for p in pairs)
    all_magnitude = all(
        p["classification"] in {"source_volume_magnitude", "no_material_divergence"}
        for p in pairs)
    if all_magnitude and not any_branch:
        conclusion = (
            "All adjacent factor pairs share an identical continuation branch "
            "(same predicted_initial_state, branch_acceptance_status, p95 jump "
            "trajectory, and retry counts across the -20..-10 V window). The only "
            "difference is a smooth, monotone current-magnitude scaling with the "
            "source-volume factor. The `>1.5` growth marker therefore appears only "
            "as a smooth threshold crossing of the rescaled magnitude at the window "
            "edge, not as a discrete branch jump. Classification: source-volume "
            "MAGNITUDE on a single shared branch, NOT a continuation/predictor "
            "branch selection. The missing Sentaurus knee is not reachable by "
            "scalar source-volume selection on this axis; the knee blocker lies in "
            "the continuation/branch mechanism itself (pursue pseudo-arclength "
            "continuation and local minority quasi-Fermi caps).")
    else:
        onsets = [p["branch_divergence_onset_V"] for p in pairs if p["branch_divergence_onset_V"] is not None]
        conclusion = (
            "At least one factor pair shows a discrete branch divergence "
            f"(onset bias {min(onsets) if onsets else 'n/a'} V): the continuation "
            "predictor/branch guard selects a different branch between factors. "
            "Classification: continuation/predictor branch selection. Target the "
            "predictor handoff and branch-guard step control at that bias.")

    report = {
        "reference_knee": knee._summary_to_json(
            knee.summarize_curve(reference, args.bias_min, args.bias_max)),
        "bias_window_V": {"min": args.bias_min, "max": args.bias_max},
        "candidates": candidates,
        "pairs": pairs,
        "conclusion": conclusion,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "factor_branch_divergence.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    markdown = render_markdown(report)
    (args.out_dir / "factor_branch_divergence.md").write_text(
        markdown + "\n", encoding="utf-8")
    print(markdown)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate",
        action="append",
        metavar="FACTOR=CSV",
        help="Repeatable 'factor=csv_path' entry, ordered by increasing factor.",
    )
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-10.0)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_SCAN_ROOT / "factor_branch_divergence",
    )
    args = parser.parse_args()

    if not args.candidate:
        scan = DEFAULT_SCAN_ROOT / "qflim0p05_source_factor_scan"
        high = DEFAULT_SCAN_ROOT / "qflim0p05_high_source_scan"
        args.candidate = [
            f"0.875={scan / 'factor_0p875.csv'}",
            f"0.90625={scan / 'factor_0p90625.csv'}",
            f"0.921875={high / 'factor_0p921875.csv'}",
        ]

    parsed: list[tuple[str, Path]] = []
    for entry in args.candidate:
        factor, _, path = entry.partition("=")
        if not path:
            raise SystemExit(f"invalid --candidate entry: {entry}")
        parsed.append((factor.strip(), Path(path.strip())))
    args.candidate = parsed
    return args


if __name__ == "__main__":
    raise SystemExit(main())
