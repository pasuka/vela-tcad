#!/usr/bin/env python3
"""Diagnose PN2D BV contact density differences via inferred BGN ni_eff."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


EDGE_FIELDS = [
    "bias_V",
    "contact",
    "edge_id",
    "carrier",
    "vela_density_cm3",
    "sentaurus_density_cm3",
    "log10_density_vela_over_sentaurus",
    "vela_drive_V",
    "sentaurus_drive_V",
    "drive_delta_V",
    "vela_inferred_ni_eff_cm3",
    "sentaurus_inferred_ni_eff_cm3",
    "log10_inferred_ni_eff_vela_over_sentaurus",
    "vela_electron_inferred_ni_eff_cm3",
    "sentaurus_electron_inferred_ni_eff_cm3",
    "vela_hole_inferred_ni_eff_cm3",
    "sentaurus_hole_inferred_ni_eff_cm3",
]

SUMMARY_FIELDS = [
    "bias_V",
    "contact",
    "carrier",
    "points",
    "log10_density_p50_vela_over_sentaurus",
    "log10_inferred_ni_eff_p50_vela_over_sentaurus",
    "drive_delta_V_p50",
    "drive_contribution_decade_p50",
    "classification",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-edge-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", required=True)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--thermal-voltage", type=float, default=0.025852)
    parser.add_argument("--mismatch-threshold-decade", type=float, default=0.1)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def value(row: dict[str, str], key: str) -> float | None:
    raw = row.get(key)
    if raw in (None, ""):
        return None
    try:
        result = float(raw)
    except ValueError:
        return None
    if not math.isfinite(result):
        return None
    return result


def log10_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate == 0.0 or reference == 0.0:
        return None
    return math.log10(abs(candidate) / abs(reference))


def limited_exp(raw: float) -> float:
    return math.exp(max(-500.0, min(500.0, raw)))


def carrier_drive(psi: float, qf: float, carrier: str) -> float:
    if carrier == "electron":
        return psi - qf
    if carrier == "hole":
        return qf - psi
    raise ValueError(carrier)


def inferred_ni_eff_cm3(density_cm3: float | None, drive_v: float | None, vt: float) -> float | None:
    if density_cm3 is None or drive_v is None or density_cm3 <= 0.0:
        return None
    return density_cm3 / limited_exp(drive_v / vt)


def percentile(values: list[float], p: float) -> float | None:
    clean = sorted(values)
    if not clean:
        return None
    idx = (len(clean) - 1) * p / 100.0
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return clean[lo]
    return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)


def classify(
    density_log: float | None,
    ni_log: float | None,
    drive_contribution_decade: float | None,
    threshold: float,
) -> str:
    if density_log is None:
        return "missing_density"
    if abs(density_log) <= threshold:
        return "density_close"
    if ni_log is not None and abs(ni_log) > threshold:
        return "ni_eff_mismatch"
    if drive_contribution_decade is not None and abs(drive_contribution_decade) > threshold:
        return "qf_potential_drive_mismatch"
    return "density_mismatch"


def edge_rows(
    rows: list[dict[str, str]],
    target_biases: set[float],
    contact: str,
    vt: float,
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for row in rows:
        bias = value(row, "bias_V")
        if bias is None or bias_key(bias) not in target_biases or row.get("current_contact") != contact:
            continue
        vela_psi = value(row, "vela_potential_avg_V")
        sent_psi = value(row, "sentaurus_potential_V")
        for carrier in ["electron", "hole"]:
            vela_density = value(row, f"vela_{carrier}_density_avg_cm3")
            sent_density = value(row, f"sentaurus_{carrier}_density_cm3")
            vela_qf = value(row, f"vela_{carrier}_qf_avg_V")
            sent_qf = value(row, f"sentaurus_{carrier}_qf_V")
            vela_drive = carrier_drive(vela_psi, vela_qf, carrier) if vela_psi is not None and vela_qf is not None else None
            sent_drive = carrier_drive(sent_psi, sent_qf, carrier) if sent_psi is not None and sent_qf is not None else None
            vela_ni = inferred_ni_eff_cm3(vela_density, vela_drive, vt)
            sent_ni = inferred_ni_eff_cm3(sent_density, sent_drive, vt)
            out = {
                "bias_V": bias_key(bias),
                "contact": contact,
                "edge_id": row.get("edge_id", ""),
                "carrier": carrier,
                "vela_density_cm3": vela_density,
                "sentaurus_density_cm3": sent_density,
                "log10_density_vela_over_sentaurus": log10_ratio(vela_density, sent_density),
                "vela_drive_V": vela_drive,
                "sentaurus_drive_V": sent_drive,
                "drive_delta_V": None if vela_drive is None or sent_drive is None else vela_drive - sent_drive,
                "vela_inferred_ni_eff_cm3": vela_ni,
                "sentaurus_inferred_ni_eff_cm3": sent_ni,
                "log10_inferred_ni_eff_vela_over_sentaurus": log10_ratio(vela_ni, sent_ni),
            }
            out[f"vela_{carrier}_inferred_ni_eff_cm3"] = vela_ni
            out[f"sentaurus_{carrier}_inferred_ni_eff_cm3"] = sent_ni
            result.append(out)
    return result


def summary_rows(edges: list[dict[str, Any]], biases: list[float], contact: str, vt: float, threshold: float) -> list[dict[str, Any]]:
    groups: dict[tuple[float, str], list[dict[str, Any]]] = defaultdict(list)
    for row in edges:
        groups[(bias_key(float(row["bias_V"])), str(row["carrier"]))].append(row)

    rows: list[dict[str, Any]] = []
    for bias in biases:
        for carrier in ["electron", "hole"]:
            items = groups.get((bias_key(bias), carrier), [])
            density_logs = [float(item["log10_density_vela_over_sentaurus"]) for item in items if item.get("log10_density_vela_over_sentaurus") is not None]
            ni_logs = [float(item["log10_inferred_ni_eff_vela_over_sentaurus"]) for item in items if item.get("log10_inferred_ni_eff_vela_over_sentaurus") is not None]
            drive_deltas = [float(item["drive_delta_V"]) for item in items if item.get("drive_delta_V") is not None]
            density_p50 = percentile(density_logs, 50)
            ni_p50 = percentile(ni_logs, 50)
            drive_p50 = percentile(drive_deltas, 50)
            drive_decade = None if drive_p50 is None else drive_p50 / (vt * math.log(10.0))
            rows.append({
                "bias_V": bias_key(bias),
                "contact": contact,
                "carrier": carrier,
                "points": len(items),
                "log10_density_p50_vela_over_sentaurus": density_p50,
                "log10_inferred_ni_eff_p50_vela_over_sentaurus": ni_p50,
                "drive_delta_V_p50": drive_p50,
                "drive_contribution_decade_p50": drive_decade,
                "classification": classify(density_p50, ni_p50, drive_decade, threshold),
            })
    return rows


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    biases = parse_biases(args.biases)
    selected_edges = edge_rows(
        read_csv_rows(args.contact_edge_csv),
        {bias_key(bias) for bias in biases},
        args.contact,
        args.thermal_voltage,
    )
    summaries = summary_rows(
        selected_edges,
        biases,
        args.contact,
        args.thermal_voltage,
        args.mismatch_threshold_decade,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "contact_bgn_density_edges.csv", EDGE_FIELDS, selected_edges)
    write_csv(args.out_dir / "contact_bgn_density_summary.csv", SUMMARY_FIELDS, summaries)
    (args.out_dir / "summary.json").write_text(
        json.dumps({
            "contact_edge_csv": str(args.contact_edge_csv),
            "contact": args.contact,
            "biases": biases,
            "thermal_voltage": args.thermal_voltage,
            "rows": summaries,
        }, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
