#!/usr/bin/env python3
"""Report signed PN2D field differences by spatial partition."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import compare_pn2d_bv_multibias_fields as field_compare


DENSITY_QUANTITIES = {"electron_density", "hole_density", "srh_recombination", "avalanche_generation"}


def main() -> int:
    args = parse_args()
    quantities = field_compare.parse_csv_list(args.quantities)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sx, sy, _ = field_compare.load_sentaurus_mesh(args.sentaurus_dir)
    node_ids = field_compare.load_sentaurus_node_ids(args.sentaurus_dir)
    masks = field_compare.load_region_masks(args.sentaurus_dir, sx, sy)
    masks = {
        "all": np.ones(len(node_ids), dtype=bool),
        "centerline": centerline_mask(sy),
        **masks,
    }

    vtk = field_compare.parse_vtk(args.vela_vtk)
    vx = np.asarray(vtk["x"], dtype=float)
    vy = np.asarray(vtk["y"], dtype=float)

    summary_rows: list[dict[str, Any]] = []
    top_rows: list[dict[str, Any]] = []

    for quantity in quantities:
        if quantity not in field_compare.FIELD_SPECS:
            raise SystemExit(f"unknown quantity: {quantity}")
        spec = field_compare.FIELD_SPECS[quantity]
        loaded = field_compare.load_sentaurus_scalar(args.sentaurus_dir, spec["sentaurus"])
        vraw = vtk["scalars"].get(spec["vela"])
        if loaded is None or vraw is None:
            summary_rows.append({
                "bias_V": args.bias,
                "quantity": quantity,
                "region": "all",
                "status": "missing_field",
            })
            continue
        sentaurus_field, svals = loaded
        vvals = field_compare.nearest_values(
            vx,
            vy,
            np.asarray(vraw, dtype=float) * float(spec["scale"]),
            sx,
            sy,
        )
        signed = signed_error(quantity, svals, vvals)
        abs_error = np.abs(signed)

        for region, mask in masks.items():
            row = summarize_region(
                args.bias,
                quantity,
                sentaurus_field,
                region,
                mask,
                svals,
                vvals,
                signed,
                abs_error,
            )
            summary_rows.append(row)

        top_rows.extend(top_error_rows(args.bias, quantity, node_ids, sx, sy, svals, vvals, signed, args.top_n))

    write_csv(args.out_dir / "signed_field_partition_summary.csv", summary_rows)
    write_csv(args.out_dir / "signed_field_partition_top_errors.csv", top_rows)
    summary = {
        "bias_V": args.bias,
        "sentaurus_dir": str(args.sentaurus_dir),
        "vela_vtk": str(args.vela_vtk),
        "summary_row_count": len(summary_rows),
        "top_error_row_count": len(top_rows),
        "quantities": quantities,
    }
    (args.out_dir / "signed_field_partition_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )

    print(json.dumps(summary, sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument(
        "--quantities",
        default="potential,electric_field,electron_density,hole_density",
    )
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def centerline_mask(y_um: np.ndarray) -> np.ndarray:
    indices = field_compare.centerline_indices(y_um)
    mask = np.zeros(len(y_um), dtype=bool)
    mask[indices] = True
    return mask


def signed_error(quantity: str, reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    signed = np.full_like(reference, math.nan, dtype=float)
    mask = np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(mask):
        return signed
    if quantity in DENSITY_QUANTITIES:
        eps = 1.0e-300
        signed[mask] = (
            np.log10(np.maximum(np.abs(candidate[mask]), eps))
            - np.log10(np.maximum(np.abs(reference[mask]), eps))
        )
    elif quantity == "potential":
        signed[mask] = candidate[mask] - reference[mask]
    else:
        denom = np.maximum(np.abs(reference[mask]), 1.0e-30)
        signed[mask] = (candidate[mask] - reference[mask]) / denom
    return signed


def summarize_region(
    bias: float,
    quantity: str,
    sentaurus_field: str,
    region: str,
    mask: np.ndarray,
    reference: np.ndarray,
    candidate: np.ndarray,
    signed: np.ndarray,
    abs_error: np.ndarray,
) -> dict[str, Any]:
    finite = mask & np.isfinite(signed)
    if not np.any(finite):
        return {
            "bias_V": bias,
            "quantity": quantity,
            "sentaurus_field": sentaurus_field,
            "region": region,
            "status": "empty",
            "count": 0,
        }
    signed_values = signed[finite]
    abs_values = abs_error[finite]
    ref_values = reference[finite]
    cand_values = candidate[finite]
    median_signed = float(np.nanmedian(signed_values))
    row = {
        "bias_V": bias,
        "quantity": quantity,
        "sentaurus_field": sentaurus_field,
        "region": region,
        "status": "ok",
        "count": int(np.sum(finite)),
        "signed_log10_median_vela_over_sentaurus": "",
        "signed_log10_p95_abs": "",
        "signed_delta_median": "",
        "signed_delta_p95_abs": "",
        "sentaurus_median": float(np.nanmedian(ref_values)),
        "vela_median": float(np.nanmedian(cand_values)),
        "direction": direction(median_signed),
    }
    if quantity in DENSITY_QUANTITIES:
        row["signed_log10_median_vela_over_sentaurus"] = median_signed
        row["signed_log10_p95_abs"] = float(np.nanpercentile(abs_values, 95))
        row["vela_over_sentaurus_median"] = 10.0 ** median_signed
    else:
        row["signed_delta_median"] = median_signed
        row["signed_delta_p95_abs"] = float(np.nanpercentile(abs_values, 95))
        row["vela_over_sentaurus_median"] = ""
    return row


def direction(value: float) -> str:
    if not math.isfinite(value):
        return "unknown"
    if value > 1.0e-12:
        return "vela_higher"
    if value < -1.0e-12:
        return "vela_lower"
    return "matched"


def top_error_rows(
    bias: float,
    quantity: str,
    node_ids: list[int],
    x_um: np.ndarray,
    y_um: np.ndarray,
    reference: np.ndarray,
    candidate: np.ndarray,
    signed: np.ndarray,
    top_n: int,
) -> list[dict[str, Any]]:
    finite_indices = [index for index, value in enumerate(signed) if math.isfinite(float(value))]
    finite_indices.sort(key=lambda index: abs(float(signed[index])), reverse=True)
    rows: list[dict[str, Any]] = []
    for rank, index in enumerate(finite_indices[:max(top_n, 0)], start=1):
        rows.append({
            "bias_V": bias,
            "quantity": quantity,
            "rank": rank,
            "node_id": node_ids[index],
            "x_um": float(x_um[index]),
            "y_um": float(y_um[index]),
            "sentaurus_value": float(reference[index]),
            "vela_value": float(candidate[index]),
            "signed_error": float(signed[index]),
            "abs_error": abs(float(signed[index])),
            "direction": direction(float(signed[index])),
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
