#!/usr/bin/env python3
"""Read-only Anode hole-QF-drop sensitivity for PN2D BV contact current."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Q = 1.602176634e-19
VT_300 = 0.025851999786435535


EDGE_FIELDS = [
    "bias_V",
    "current_contact",
    "edge_id",
    "baseline_hole_qf_drop_V",
    "sentaurus_hole_qf_drop_V",
    "baseline_hole_current_A_per_um",
    "sentaurus_drop_hole_current_A_per_um",
    "hole_current_delta_A_per_um",
    "baseline_total_current_A_per_um",
    "sentaurus_drop_total_current_A_per_um",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-edges", type=Path, required=True)
    parser.add_argument("--edge-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--thermal-voltage", type=float, default=VT_300)
    return parser.parse_args()


def optional_float(raw: object) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def bernoulli(x: float) -> float:
    if abs(x) < 1.0e-8:
        return 1.0 - x / 2.0 + x * x / 12.0
    if x > 100.0:
        return x * math.exp(-x)
    if x < -100.0:
        return -x
    return x / math.expm1(x)


def limited_exp(value: float) -> float:
    return math.exp(max(-500.0, min(500.0, value)))


def sg_hole_current_from_qf_drop_a_per_m(row: dict[str, str], qf_drop: float, vt: float) -> float | None:
    psi0 = optional_float(row.get("psi0"))
    psi1 = optional_float(row.get("psi1"))
    phip0 = optional_float(row.get("phip0"))
    ni0 = optional_float(row.get("ni0"))
    ni1 = optional_float(row.get("ni1"))
    mup = optional_float(row.get("mup"))
    edge_length = optional_float(row.get("edge_length_m"))
    edge_couple = optional_float(row.get("edge_couple_m"))
    outward_sign = optional_float(row.get("outward_sign"))
    if None in {psi0, psi1, phip0, ni0, ni1, mup, edge_length, edge_couple, outward_sign}:
        return None
    if edge_length <= 0.0 or ni0 <= 0.0 or ni1 <= 0.0 or mup <= 0.0:
        return None

    phip1 = phip0 + qf_drop
    if phip0 == phip1:
        continuity_flux = 0.0
    else:
        eta = (psi1 - psi0) / vt + math.log(ni0 / ni1)
        p0 = ni0 * limited_exp((phip0 - psi0) / vt)
        p1 = ni1 * limited_exp((phip1 - psi1) / vt)
        coef = mup * vt / edge_length
        continuity_flux = coef * (bernoulli(eta) * p0 - bernoulli(-eta) * p1)
    hole_flux = -continuity_flux
    return Q * outward_sign * hole_flux * edge_couple


def current_a_per_um(row: dict[str, str], key: str) -> float | None:
    direct = optional_float(row.get(f"{key}_A_per_um"))
    if direct is not None:
        return direct
    value = optional_float(row.get(key))
    if value is not None:
        return value / 1.0e6
    return None


def log10_error(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate == 0.0 or reference == 0.0:
        return None
    return abs(math.log10(abs(candidate) / abs(reference)))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    edge_rows = [
        row for row in load_csv(args.contact_edges)
        if row.get("current_contact") == args.contact
        and abs((optional_float(row.get("bias_V")) or math.inf) - args.bias) <= 1.0e-9
    ]
    summary_rows = [
        row for row in load_csv(args.edge_summary)
        if abs((optional_float(row.get("bias_V")) or math.inf) - args.bias) <= 1.0e-9
    ]
    reference_current = None
    if summary_rows:
        reference_current = (
            optional_float(summary_rows[0].get("sentaurus_plt_current_A"))
            or optional_float(summary_rows[0].get("sentaurus_current_A"))
        )

    out_rows: list[dict[str, Any]] = []
    baseline_total = 0.0
    sensitivity_total = 0.0
    baseline_hole_total = 0.0
    sensitivity_hole_total = 0.0
    for row in edge_rows:
        baseline_hole = current_a_per_um(row, "current_hole") or 0.0
        baseline_total_edge = current_a_per_um(row, "current_total") or 0.0
        electron = current_a_per_um(row, "current_electron") or 0.0
        baseline_drop = optional_float(row.get("vela_hole_qf_drop_V"))
        if baseline_drop is None:
            phip0 = optional_float(row.get("phip0"))
            phip1 = optional_float(row.get("phip1"))
            baseline_drop = (phip1 - phip0) if phip0 is not None and phip1 is not None else None
        sentaurus_drop = optional_float(row.get("sentaurus_hole_qf_drop_V"))
        sensitivity_hole_a_per_m = (
            sg_hole_current_from_qf_drop_a_per_m(row, sentaurus_drop, args.thermal_voltage)
            if sentaurus_drop is not None
            else None
        )
        sensitivity_hole = (
            baseline_hole
            if sensitivity_hole_a_per_m is None
            else sensitivity_hole_a_per_m / 1.0e6
        )
        sensitivity_total_edge = electron - sensitivity_hole

        baseline_total += baseline_total_edge
        sensitivity_total += sensitivity_total_edge
        baseline_hole_total += baseline_hole
        sensitivity_hole_total += sensitivity_hole
        out_rows.append({
            "bias_V": args.bias,
            "current_contact": args.contact,
            "edge_id": row.get("edge_id"),
            "baseline_hole_qf_drop_V": baseline_drop,
            "sentaurus_hole_qf_drop_V": sentaurus_drop,
            "baseline_hole_current_A_per_um": baseline_hole,
            "sentaurus_drop_hole_current_A_per_um": sensitivity_hole,
            "hole_current_delta_A_per_um": sensitivity_hole - baseline_hole,
            "baseline_total_current_A_per_um": baseline_total_edge,
            "sentaurus_drop_total_current_A_per_um": sensitivity_total_edge,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "hole_qf_sensitivity_edges.csv", EDGE_FIELDS, out_rows)
    baseline_error = log10_error(baseline_total, reference_current)
    sensitivity_error = log10_error(sensitivity_total, reference_current)
    summary = {
        "bias_V": args.bias,
        "contact": args.contact,
        "edge_count": len(out_rows),
        "reference_current_A": reference_current,
        "baseline_total_current_A_per_um": baseline_total,
        "sentaurus_drop_total_current_A_per_um": sensitivity_total,
        "baseline_hole_current_A_per_um": baseline_hole_total,
        "sentaurus_drop_hole_current_A_per_um": sensitivity_hole_total,
        "abs_log10_error_baseline": baseline_error,
        "abs_log10_error_sentaurus_drop": sensitivity_error,
    }
    (args.out_dir / "hole_qf_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
