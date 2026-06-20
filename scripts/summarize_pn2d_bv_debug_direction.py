#!/usr/bin/env python3
"""Summarize PN2D BV debug evidence into the next investigation direction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--debug-ranking", type=Path, required=True)
    parser.add_argument("--terminal-crosscheck", type=Path, required=True)
    parser.add_argument("--qf-floor-stability", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def by_quantity(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(row.get("quantity")): row for row in rows}


def assess_terminal(
    terminal: dict[str, Any],
    qf_stability: dict[str, Any],
    ranking: dict[str, Any],
) -> str:
    status = str(terminal.get("status", ""))
    terminal_dex = abs(finite_float(terminal.get("log10_contact_flux_over_plt")) or 0.0)
    qf_error = finite_float(qf_stability.get("qf_floor_abs_log10_error"))
    curve_error = finite_float(ranking.get("curve", {}).get("max_abs_log10_error"))
    current_error = qf_error if qf_error is not None else curve_error
    if status.endswith("_mismatch") and current_error is not None:
        if current_error <= terminal_dex * 1.1:
            return "plt_gap_explained_by_sentaurus_terminal_convention"
        return "plt_gap_exceeds_sentaurus_terminal_convention"
    if status.endswith("_consistent"):
        return "sentaurus_terminal_definitions_consistent"
    return "sentaurus_terminal_crosscheck_incomplete"


def make_next_actions(fields: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    actions: list[dict[str, Any]] = []
    avalanche = fields.get("avalanche_generation", {})
    if avalanche:
        status = str(avalanche.get("avalanche_status", ""))
        field_error = finite_float(avalanche.get("field_error"))
        junction_error = finite_float(avalanche.get("junction_error"))
        if status == "thresholded_peak" or (field_error is not None and field_error > 1.0):
            actions.append({
                "id": "thresholded_avalanche_support",
                "priority": 1,
                "reason": (
                    "Avalanche mismatch is dominated by thresholded/near-floor support; "
                    "compare high-generation support near the junction before changing alpha formulas."
                ),
                "field_error": field_error,
                "junction_error": junction_error,
            })
    electric = fields.get("electric_field", {})
    if electric:
        field_error = finite_float(electric.get("field_error"))
        junction_error = finite_float(electric.get("junction_error"))
        actions.append({
            "id": "junction_electric_field_reconstruction",
            "priority": 2,
            "reason": (
                "Global electric-field relative error is denominator-sensitive; "
                "reconstruct and compare the junction/high-field stencil directly."
            ),
            "field_error": field_error,
            "junction_error": junction_error,
        })
    if fields.get("electron_density") or fields.get("electron_mobility"):
        actions.append({
            "id": "carrier_density_mobility_followup",
            "priority": 3,
            "reason": (
                "Carrier density and mobility remain secondary contributors after "
                "terminal convention and high-field avalanche support are separated."
            ),
        })
    return sorted(actions, key=lambda item: int(item["priority"]))


def secondary_quantities(fields: dict[str, dict[str, Any]]) -> list[str]:
    ordered = ["electron_density", "electron_mobility", "hole_density", "hole_mobility"]
    return [name for name in ordered if name in fields]


def build_summary(
    ranking: dict[str, Any],
    terminal: dict[str, Any],
    qf_stability: dict[str, Any],
) -> dict[str, Any]:
    fields = by_quantity(list(ranking.get("field_rankings", [])))
    return {
        "terminal_assessment": assess_terminal(terminal, qf_stability, ranking),
        "current_error": {
            "qf_floor_abs_log10_error": qf_stability.get("qf_floor_abs_log10_error"),
            "debug_ranking_max_abs_log10_error": ranking.get("curve", {}).get("max_abs_log10_error"),
            "sentaurus_contact_flux_vs_plt_log10": terminal.get("log10_contact_flux_over_plt"),
            "sentaurus_contact_flux_vs_plt_relative": terminal.get("relative_difference"),
        },
        "field_evidence": {
            "electric_field": fields.get("electric_field"),
            "avalanche_generation": fields.get("avalanche_generation"),
        },
        "secondary_quantities": secondary_quantities(fields),
        "next_actions": make_next_actions(fields),
    }


def main() -> int:
    args = parse_args()
    summary = build_summary(
        load_json(args.debug_ranking),
        load_json(args.terminal_crosscheck),
        load_json(args.qf_floor_stability),
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
