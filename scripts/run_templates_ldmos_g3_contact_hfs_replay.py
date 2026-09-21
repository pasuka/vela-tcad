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

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import run_sg_probe
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_idvg_shift_kcl import run_sg_probe


BIAS_LABELS = (
    (1.0 / 6.0, "vg_0p166667"),
    (0.5, "vg_0p500000"),
    (5.0 / 6.0, "vg_0p833333"),
)


def contact_hfs_variants(baseline: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build a same-contract frozen-state pair differing only at contact HFS."""
    hfs_off = deepcopy(baseline)
    hfs_on = deepcopy(baseline)
    off_mobility = hfs_off["solver"]["mobility"]
    on_mobility = hfs_on["solver"]["mobility"]
    off_mobility["contact_electric_field_fallback"] = False
    off_mobility.pop("contact_electric_field_fallback_scope", None)
    off_mobility.pop("contact_electric_field_fallback_mode", None)
    on_mobility.update({
        "contact_electric_field_fallback": True,
        "contact_electric_field_fallback_scope": "contact_node_cell",
        "contact_electric_field_fallback_mode": "cell_gradient_magnitude",
    })
    return hfs_off, hfs_on


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
    hfs_off, hfs_on = contact_hfs_variants(baseline)

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
        baseline_probe = run_sg_probe(
            runner, deepcopy(hfs_off), state, bias,
            output_dir / label / "contact_hfs_off")
        probe = run_sg_probe(
            runner, deepcopy(hfs_on), state, bias,
            output_dir / label / "contact_hfs_on")
        baseline_current = float(
            baseline_probe["drain_cut"]["total_A_per_um"])
        current = float(probe["drain_cut"]["total_A_per_um"])
        terminal = float(old["sentaurus_terminal_A_per_um"])
        ratio = current / terminal
        error = abs(
            math.log10(max(abs(current), 1.0e-300))
            - math.log10(max(abs(terminal), 1.0e-300))
        )
        baseline_ratio = baseline_current / terminal
        old_error = abs(math.log10(max(abs(baseline_ratio), 1.0e-300)))
        points.append({
            "bias_V": bias,
            "sentaurus_terminal_A_per_um": terminal,
            "archived_mesh_default_baseline_vela_sg_A_per_um": float(
                old["vela_sg_on_sentaurus_state_A_per_um"]),
            "archived_mesh_default_baseline_signed_ratio": float(
                old["signed_ratio"]),
            "archived_mesh_default_baseline_magnitude_error_dex": float(
                old["magnitude_error_dex"]),
            "same_contract_hfs_off_vela_sg_A_per_um": baseline_current,
            "same_contract_hfs_off_signed_ratio": baseline_ratio,
            "same_contract_hfs_off_magnitude_error_dex": old_error,
            "contact_fallback_vela_sg_A_per_um": current,
            "contact_fallback_signed_ratio": ratio,
            "contact_fallback_magnitude_error_dex": error,
            "improvement_dex": old_error - error,
            "hfs_off_probe": baseline_probe,
            "probe": probe,
        })

    old_median = statistics.median(
        point["same_contract_hfs_off_magnitude_error_dex"] for point in points)
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
        "schema": "vela.templates_ldmos.g3_contact_hfs_replay.v2",
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
            "comparison": "same_frozen_state_and_same_non_hfs_operator",
            "carrier_transport_couple_profile": baseline.get(
                "mesh_geometry", {}
            ).get("carrier_transport_couple_profile", "mesh_default"),
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
