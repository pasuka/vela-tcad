#!/usr/bin/env python3
"""Replay three archived G3 Sentaurus states with the contact HFS fallback."""

from __future__ import annotations

import argparse
import json
import math
import statistics
from copy import deepcopy
from pathlib import Path
from typing import Any

from audit_templates_ldmos_g3_idvg_shift_kcl import run_sg_probe


BIAS_LABELS = (
    (1.0 / 6.0, "vg_0p166667"),
    (0.5, "vg_0p500000"),
    (5.0 / 6.0, "vg_0p833333"),
)


def replay(
    runner: Path,
    baseline_config: Path,
    archived_replay: Path,
    archived_summary: Path,
    output_dir: Path,
) -> dict[str, Any]:
    baseline = json.loads(baseline_config.read_text(encoding="utf-8"))
    mobility = baseline["solver"]["mobility"]
    if mobility.get("high_field_driving_force") != "quasi_fermi_gradient":
        raise ValueError("G3 contact-HFS replay requires a GradQF baseline")
    if "ialmob" in json.dumps(mobility).lower():
        raise ValueError("G3 contact-HFS replay requires IALMob disabled")
    mobility.update({
        "contact_electric_field_fallback": True,
        "contact_electric_field_fallback_scope": "contact_node_cell",
        "contact_electric_field_fallback_mode": "cell_gradient_magnitude",
    })

    archived = json.loads(archived_summary.read_text(encoding="utf-8"))[
        "sentaurus_state_sg_replay"
    ]["points"]
    archived_by_bias = {
        round(float(point["bias_V"]), 12): point for point in archived
    }

    points: list[dict[str, Any]] = []
    for bias, label in BIAS_LABELS:
        old = archived_by_bias[round(bias, 12)]
        state = archived_replay / label / "sentaurus_state.csv"
        probe = run_sg_probe(
            runner, deepcopy(baseline), state, bias, output_dir / label)
        current = float(probe["drain_cut"]["total_A_per_um"])
        terminal = float(old["sentaurus_terminal_A_per_um"])
        ratio = current / terminal
        error = abs(
            math.log10(max(abs(current), 1.0e-300))
            - math.log10(max(abs(terminal), 1.0e-300))
        )
        old_error = float(old["magnitude_error_dex"])
        points.append({
            "bias_V": bias,
            "sentaurus_terminal_A_per_um": terminal,
            "baseline_vela_sg_A_per_um": float(
                old["vela_sg_on_sentaurus_state_A_per_um"]),
            "baseline_signed_ratio": float(old["signed_ratio"]),
            "baseline_magnitude_error_dex": old_error,
            "contact_fallback_vela_sg_A_per_um": current,
            "contact_fallback_signed_ratio": ratio,
            "contact_fallback_magnitude_error_dex": error,
            "improvement_dex": old_error - error,
            "probe": probe,
        })

    old_median = statistics.median(
        point["baseline_magnitude_error_dex"] for point in points)
    new_median = statistics.median(
        point["contact_fallback_magnitude_error_dex"] for point in points)
    high_points = points[1:]
    clear_improvement = all(
        point["improvement_dex"] >= 0.5 for point in high_points
    ) and all(
        0.5 <= abs(point["contact_fallback_signed_ratio"]) <= 2.0
        for point in high_points
    )
    return {
        "schema": "vela.templates_ldmos.g3_contact_hfs_replay.v1",
        "status": "pass" if clear_improvement else "fail",
        "decision": {
            "clear_improvement": clear_improvement,
            "run_31_point_curve": clear_improvement,
            "rule": (
                "both Vg=0.5 and 0.833333 V improve by at least 0.5 dex "
                "and finish within a factor of two of Sentaurus"
            ),
        },
        "contracts": {
            "interior_hfs_drive": "quasi_fermi_gradient",
            "contact_scope": "contact_node_cell",
            "contact_hfs_drive": "cell_gradient_magnitude",
            "ialmob": "disabled",
            "predictor": "disabled",
        },
        "baseline_median_magnitude_error_dex": old_median,
        "contact_fallback_median_magnitude_error_dex": new_median,
        "median_improvement_dex": old_median - new_median,
        "points": points,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--archived-replay", type=Path, required=True)
    parser.add_argument("--archived-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = replay(
        args.runner, args.baseline_config, args.archived_replay,
        args.archived_summary, args.output_dir)
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"],
        "summary": str(summary_path.resolve()),
        "run_31_point_curve": summary["decision"]["run_31_point_curve"],
    }))
    return 0 if summary["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
