#!/usr/bin/env python3
"""Step 3: partition the PN2D BV avalanche-source integral by electric field.

Steps 1-2 established that the ~0.37x Vela/Sentaurus avalanche-source integral
ratio is geometry-independent and splits into two regimes:

* genuine high-field active edges at a clean ~0.75x, and
* a broad low-field tail where Sentaurus generates but Vela's quasi-Fermi-gradient
  path gives ~0, which drags the *global* integral ratio down to ~0.37x.

The terminal-current ratio at -13.2 V is ~0.698x, i.e. much closer to the active
0.75x than to the global 0.37x. This script quantifies that by partitioning the
edge-source integral at a sweep of electric-field thresholds and reporting, for
each partition, the Vela and Sentaurus partial integrals, their ratio, and the
generation-rate current proxy ``q * integral(G dV)``. The goal is to decide
whether closing the active-edge ~0.75x deficit alone would move the terminal
current to parity (lever = SG flux magnitude on active edges) or whether the
low-field tail must also be filled (lever = low-density driving-force
interpolation to ElectricField, Sentaurus manual p.439).

It consumes the Step 2 ``edge_geometry_audit_edges.csv`` so the geometry is the
already-validated box area and no mesh re-derivation is needed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

Q = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges-csv", type=Path, required=True,
                        help="Step 2 edge_geometry_audit_edges.csv")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=-13.2)
    parser.add_argument(
        "--terminal-current-ratio", type=float, default=0.698,
        help="Known Vela/Sentaurus terminal current ratio at this bias.")
    parser.add_argument(
        "--field-thresholds", type=str,
        default="0,1e6,5e6,1e7,1.5e7,2e7,2.4e7",
        help="Comma-separated E thresholds (V/m) for the cumulative partition.")
    return parser.parse_args()


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def read_edges(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            field = optional_float(row.get("electric_field_V_per_m"))
            vela = optional_float(row.get("edge_source_integral_s_inv"))
            sent = optional_float(row.get("sentaurus_edge_source_s_inv"))
            if field is None or vela is None or sent is None:
                continue
            rows.append({"E": field, "vela": vela, "sent": sent})
    return rows


def safe_ratio(num: float, den: float) -> float | None:
    if den == 0.0:
        return None
    return num / den


def partition_stats(rows: list[dict[str, float]],
                    lo: float, hi: float | None) -> dict[str, Any]:
    vela = 0.0
    sent = 0.0
    count = 0
    for r in rows:
        e = r["E"]
        if e < lo:
            continue
        if hi is not None and e >= hi:
            continue
        vela += r["vela"]
        sent += r["sent"]
        count += 1
    return {
        "edge_count": count,
        "vela_source_s_inv": vela,
        "sentaurus_source_s_inv": sent,
        "vela_over_sentaurus": safe_ratio(vela, sent),
        "vela_gen_current_A_per_um": Q * vela * 1.0e-6,
        "sentaurus_gen_current_A_per_um": Q * sent * 1.0e-6,
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_edges(args.edges_csv)

    total_vela = sum(r["vela"] for r in rows)
    total_sent = sum(r["sent"] for r in rows)

    thresholds = [float(x) for x in args.field_thresholds.split(",") if x != ""]
    thresholds = sorted(set(thresholds))

    # Cumulative "active" partitions: E >= threshold.
    cumulative: list[dict[str, Any]] = []
    for thr in thresholds:
        block = partition_stats(rows, thr, None)
        block["field_threshold_V_per_m"] = thr
        block["vela_fraction_of_total"] = safe_ratio(
            block["vela_source_s_inv"], total_vela)
        block["sentaurus_fraction_of_total"] = safe_ratio(
            block["sentaurus_source_s_inv"], total_sent)
        cumulative.append(block)

    # Disjoint bands between consecutive thresholds (plus the top open band).
    bands: list[dict[str, Any]] = []
    edges_seq = thresholds + [None]
    for lo, hi in zip(edges_seq[:-1], edges_seq[1:]):
        block = partition_stats(rows, lo, hi)
        block["field_lo_V_per_m"] = lo
        block["field_hi_V_per_m"] = hi
        bands.append(block)

    # Active vs tail split at the physically meaningful avalanche onset.
    # Avalanche alpha is negligible below ~1e7 V/m for silicon Van Overstraeten,
    # so treat E >= 1e7 as the active region and E < 1e7 as the low-field tail.
    active = partition_stats(rows, 1.0e7, None)
    tail = partition_stats(rows, 0.0, 1.0e7)

    report: dict[str, Any] = {
        "bias_V": args.bias,
        "terminal_current_ratio_vela_over_sentaurus": args.terminal_current_ratio,
        "edge_count": len(rows),
        "total_vela_source_s_inv": total_vela,
        "total_sentaurus_source_s_inv": total_sent,
        "global_source_ratio": safe_ratio(total_vela, total_sent),
        "active_region_E_ge_1e7": active,
        "low_field_tail_E_lt_1e7": tail,
        "active_vs_terminal": {
            "active_source_ratio": active["vela_over_sentaurus"],
            "terminal_current_ratio": args.terminal_current_ratio,
            "global_source_ratio": safe_ratio(total_vela, total_sent),
            "active_minus_terminal": (
                None if active["vela_over_sentaurus"] is None
                else active["vela_over_sentaurus"] - args.terminal_current_ratio),
            "global_minus_terminal": (
                safe_ratio(total_vela, total_sent) - args.terminal_current_ratio),
        },
        "tail_contribution": {
            "sentaurus_tail_fraction_of_total_gen": safe_ratio(
                tail["sentaurus_source_s_inv"], total_sent),
            "vela_tail_fraction_of_total_gen": safe_ratio(
                tail["vela_source_s_inv"], total_vela),
            "sentaurus_tail_minus_vela_tail_s_inv": (
                tail["sentaurus_source_s_inv"] - tail["vela_source_s_inv"]),
        },
        "cumulative_active_partitions": cumulative,
        "disjoint_field_bands": bands,
    }

    summary_path = args.out_dir / "source_field_partition_summary.json"
    summary_path.write_text(json.dumps(report, indent=2))

    bands_path = args.out_dir / "source_field_bands.csv"
    fields = [
        "field_lo_V_per_m", "field_hi_V_per_m", "edge_count",
        "vela_source_s_inv", "sentaurus_source_s_inv", "vela_over_sentaurus",
        "vela_gen_current_A_per_um", "sentaurus_gen_current_A_per_um",
    ]
    with bands_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(bands)

    print(json.dumps(report, indent=2))
    print(f"\nWrote {summary_path}")
    print(f"Wrote {bands_path}")


if __name__ == "__main__":
    main()
