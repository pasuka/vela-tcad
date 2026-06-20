#!/usr/bin/env python3
"""Compare contact-current reporting policies for tiny Anode hole-QF floors."""

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
    "node0",
    "node1",
    "baseline_hole_qf_drop_V",
    "sentaurus_hole_qf_drop_V",
    "restart_hole_qf_drop_V",
    "ulp_floor_hole_qf_drop_V",
    "baseline_hole_current_A_per_um",
    "sentaurus_drop_hole_current_A_per_um",
    "restart_drop_hole_current_A_per_um",
    "ulp_floor_hole_current_A_per_um",
    "baseline_total_current_A_per_um",
    "sentaurus_drop_total_current_A_per_um",
    "restart_drop_total_current_A_per_um",
    "ulp_floor_total_current_A_per_um",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-edges", type=Path, required=True)
    parser.add_argument("--edge-summary", type=Path, required=True)
    parser.add_argument("--restart-state", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--thermal-voltage", type=float, default=VT_300)
    parser.add_argument("--ulp-floor-multiplier", type=float, default=5.0)
    parser.add_argument(
        "--ulp-sign-source",
        choices=("restart", "sentaurus", "baseline", "terminal"),
        default="restart",
    )
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


def optional_int(raw: object) -> int | None:
    value = optional_float(raw)
    if value is None:
        return None
    return int(value)


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_restart_phip(path: Path | None) -> dict[int, float]:
    if path is None or not path.exists():
        return {}
    values: dict[int, float] = {}
    for row in load_csv(path):
        node_id = optional_int(row.get("node_id"))
        value = optional_float(row.get("phip"))
        if node_id is not None and value is not None:
            values[node_id] = value
    return values


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


def sg_hole_current_from_qf_drop_a_per_m(
    row: dict[str, str],
    qf_drop: float,
    vt: float,
) -> float | None:
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


def sign_of(value: float | None) -> float | None:
    if value is None or value == 0.0:
        return None
    return 1.0 if value > 0.0 else -1.0


def choose_ulp_sign(
    row: dict[str, str],
    source: str,
    restart_drop: float | None,
    sentaurus_drop: float | None,
    baseline_drop: float | None,
) -> float:
    candidates: list[float | None]
    if source == "restart":
        candidates = [restart_drop, sentaurus_drop, baseline_drop]
    elif source == "sentaurus":
        candidates = [sentaurus_drop, restart_drop, baseline_drop]
    elif source == "baseline":
        candidates = [baseline_drop, restart_drop, sentaurus_drop]
    else:
        outward = optional_float(row.get("outward_sign"))
        candidates = [outward, restart_drop, sentaurus_drop, baseline_drop]
    for candidate in candidates:
        sign = sign_of(candidate)
        if sign is not None:
            return sign
    return -1.0


def restart_edge_drop(row: dict[str, str], phip_by_node: dict[int, float]) -> float | None:
    node0 = optional_int(row.get("node0"))
    node1 = optional_int(row.get("node1"))
    if node0 is None or node1 is None:
        return None
    phip0 = phip_by_node.get(node0)
    phip1 = phip_by_node.get(node1)
    if phip0 is None or phip1 is None:
        return None
    return phip1 - phip0


def recompute_hole_current_a_per_um(
    row: dict[str, str],
    qf_drop: float | None,
    vt: float,
    fallback: float,
) -> float:
    if qf_drop is None:
        return fallback
    value = sg_hole_current_from_qf_drop_a_per_m(row, qf_drop, vt)
    if value is None:
        return fallback
    return value / 1.0e6


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    restart_phip = load_restart_phip(args.restart_state)
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
    totals = {
        "baseline": 0.0,
        "sentaurus": 0.0,
        "restart": 0.0,
        "ulp": 0.0,
    }
    hole_totals = {
        "baseline": 0.0,
        "sentaurus": 0.0,
        "restart": 0.0,
        "ulp": 0.0,
    }
    restart_drop_edges = 0

    for row in edge_rows:
        baseline_hole = current_a_per_um(row, "current_hole") or 0.0
        baseline_total = current_a_per_um(row, "current_total") or 0.0
        electron = current_a_per_um(row, "current_electron") or 0.0
        phip0 = optional_float(row.get("phip0"))
        phip1 = optional_float(row.get("phip1"))
        baseline_drop = optional_float(row.get("vela_hole_qf_drop_V"))
        if baseline_drop is None and phip0 is not None and phip1 is not None:
            baseline_drop = phip1 - phip0
        sentaurus_drop = optional_float(row.get("sentaurus_hole_qf_drop_V"))
        restart_drop = restart_edge_drop(row, restart_phip)
        if restart_drop is not None:
            restart_drop_edges += 1

        ulp_sign = choose_ulp_sign(
            row,
            args.ulp_sign_source,
            restart_drop,
            sentaurus_drop,
            baseline_drop,
        )
        ulp_drop = ulp_sign * abs(math.ulp(args.bias)) * args.ulp_floor_multiplier

        sentaurus_hole = recompute_hole_current_a_per_um(
            row, sentaurus_drop, args.thermal_voltage, baseline_hole)
        restart_hole = recompute_hole_current_a_per_um(
            row, restart_drop, args.thermal_voltage, baseline_hole)
        ulp_hole = recompute_hole_current_a_per_um(
            row, ulp_drop, args.thermal_voltage, baseline_hole)

        sentaurus_total = electron - sentaurus_hole
        restart_total = electron - restart_hole
        ulp_total = electron - ulp_hole

        totals["baseline"] += baseline_total
        totals["sentaurus"] += sentaurus_total
        totals["restart"] += restart_total
        totals["ulp"] += ulp_total
        hole_totals["baseline"] += baseline_hole
        hole_totals["sentaurus"] += sentaurus_hole
        hole_totals["restart"] += restart_hole
        hole_totals["ulp"] += ulp_hole

        out_rows.append({
            "bias_V": args.bias,
            "current_contact": args.contact,
            "edge_id": row.get("edge_id"),
            "node0": row.get("node0"),
            "node1": row.get("node1"),
            "baseline_hole_qf_drop_V": baseline_drop,
            "sentaurus_hole_qf_drop_V": sentaurus_drop,
            "restart_hole_qf_drop_V": restart_drop,
            "ulp_floor_hole_qf_drop_V": ulp_drop,
            "baseline_hole_current_A_per_um": baseline_hole,
            "sentaurus_drop_hole_current_A_per_um": sentaurus_hole,
            "restart_drop_hole_current_A_per_um": restart_hole,
            "ulp_floor_hole_current_A_per_um": ulp_hole,
            "baseline_total_current_A_per_um": baseline_total,
            "sentaurus_drop_total_current_A_per_um": sentaurus_total,
            "restart_drop_total_current_A_per_um": restart_total,
            "ulp_floor_total_current_A_per_um": ulp_total,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "contact_qf_floor_reporting_edges.csv", EDGE_FIELDS, out_rows)
    summary = {
        "bias_V": args.bias,
        "contact": args.contact,
        "edge_count": len(out_rows),
        "restart_drop_edges": restart_drop_edges,
        "reference_current_A": reference_current,
        "baseline_total_current_A_per_um": totals["baseline"],
        "sentaurus_drop_total_current_A_per_um": totals["sentaurus"],
        "restart_drop_total_current_A_per_um": totals["restart"],
        "ulp_floor_total_current_A_per_um": totals["ulp"],
        "baseline_hole_current_A_per_um": hole_totals["baseline"],
        "sentaurus_drop_hole_current_A_per_um": hole_totals["sentaurus"],
        "restart_drop_hole_current_A_per_um": hole_totals["restart"],
        "ulp_floor_hole_current_A_per_um": hole_totals["ulp"],
        "abs_log10_error_baseline": log10_error(totals["baseline"], reference_current),
        "abs_log10_error_sentaurus_drop": log10_error(totals["sentaurus"], reference_current),
        "abs_log10_error_restart_drop": log10_error(totals["restart"], reference_current),
        "abs_log10_error_ulp_floor": log10_error(totals["ulp"], reference_current),
        "ulp_floor_multiplier": args.ulp_floor_multiplier,
        "ulp_sign_source": args.ulp_sign_source,
        "ulp_at_bias_V": math.ulp(args.bias),
    }
    (args.out_dir / "contact_qf_floor_reporting_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
