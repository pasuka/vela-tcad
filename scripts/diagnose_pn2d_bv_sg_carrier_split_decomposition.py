#!/usr/bin/env python3
"""Merge SG flux-form endpoints with edge current carrier-split rows."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


CSV_FIELDS = [
    "bias_V",
    "variant",
    "support_node_id",
    "support_class",
    "y_um",
    "edge_id",
    "edge_axis",
    "node0",
    "node1",
    "vela_psi_drop_V",
    "sentaurus_psi_drop_V",
    "mixed_psi_drop_V",
    "mixed_phin_drop_V",
    "mixed_phip_drop_V",
    "mixed_electron_eta",
    "mixed_hole_eta",
    "mixed_electron_bernoulli_plus",
    "mixed_electron_bernoulli_minus",
    "mixed_hole_bernoulli_plus",
    "mixed_hole_bernoulli_minus",
    "mixed_electron_density_factor_0",
    "mixed_electron_density_factor_1",
    "mixed_hole_density_factor_0",
    "mixed_hole_density_factor_1",
    "mixed_electron_density_factor_ratio_1_over_0",
    "mixed_hole_density_factor_ratio_1_over_0",
    "vela_electron_mobility_over_sentaurus",
    "vela_hole_mobility_over_sentaurus",
    "sentaurus_electron_fraction",
    "sentaurus_hole_fraction",
    "mixed_electron_fraction",
    "mixed_hole_fraction",
    "electron_fraction_delta",
    "hole_fraction_delta",
    "mixed_particle_over_sentaurus",
    "mixed_electron_over_sentaurus",
    "mixed_hole_over_sentaurus",
    "vela_qf_electron_over_sentaurus",
    "vela_qf_hole_over_sentaurus",
    "sentaurus_qf_electron_over_sentaurus",
    "sentaurus_qf_hole_over_sentaurus",
    "sentaurus_inferred_ni_electron_ratio_1_over_0",
    "sentaurus_inferred_ni_hole_ratio_1_over_0",
]

CORRELATION_FIELDS = [
    "electron_fraction_delta",
    "hole_fraction_delta",
    "mixed_particle_over_sentaurus",
    "mixed_electron_over_sentaurus",
    "mixed_hole_over_sentaurus",
    "mixed_electron_density_factor_ratio_1_over_0",
    "mixed_hole_density_factor_ratio_1_over_0",
    "mixed_electron_density_factor_0",
    "mixed_electron_density_factor_1",
    "mixed_hole_density_factor_0",
    "mixed_hole_density_factor_1",
    "vela_electron_mobility_over_sentaurus",
    "vela_hole_mobility_over_sentaurus",
    "sentaurus_inferred_ni_electron_ratio_1_over_0",
    "sentaurus_inferred_ni_hole_ratio_1_over_0",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-state-csv", type=Path, required=True)
    parser.add_argument("--sg-flux-form-csv", type=Path, required=True)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--support-class", default="overlap")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def read_csv(path: Path) -> list[dict[str, str]]:
    with open_path(path, "r", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(raw: Any, default: float = 0.0) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator in (None, 0.0):
        return None
    return numerator / denominator


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None or row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def flux_forms_by_edge(path: Path, bias: float | None) -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in read_csv(path):
        if not matches_bias(row, bias):
            continue
        result[(row.get("edge_id", ""), row.get("source", ""))] = row
    return result


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def pearson(xs: list[float], ys: list[float]) -> float | None:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return None
    mean_x = statistics.fmean(x for x, _ in pairs)
    mean_y = statistics.fmean(y for _, y in pairs)
    dx = [x - mean_x for x, _ in pairs]
    dy = [y - mean_y for _, y in pairs]
    denom = math.sqrt(sum(x * x for x in dx) * sum(y * y for y in dy))
    if denom == 0.0:
        return None
    return sum(x * y for x, y in zip(dx, dy)) / denom


def metric_summary(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"median": None, "min": None, "max": None}
    return {
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def density_factor(ni: float, exponent: float) -> float:
    return ni * sgdiag.limited_exp(exponent)


def split_fraction(component: float, particle: float) -> float | None:
    return safe_ratio(component, particle)


def build_row(
    edge_state: dict[str, str],
    vela: dict[str, str],
    sentaurus: dict[str, str],
    vt: float,
) -> dict[str, Any]:
    sent_e = finite_float(edge_state.get("sentaurus_edgeavg_electron_current_flux_m2_s"))
    sent_h = finite_float(edge_state.get("sentaurus_edgeavg_hole_current_flux_m2_s"))
    sent_particle = finite_float(edge_state.get("sentaurus_edgeavg_particle_flux_m2_s"))
    mixed_e = finite_float(edge_state.get("qf_old_slotboom_electron_flux_m2_s"))
    mixed_h = finite_float(edge_state.get("qf_old_slotboom_hole_flux_m2_s"))
    mixed_particle = finite_float(edge_state.get("qf_old_slotboom_particle_flux_m2_s"))

    ni0 = finite_float(vela.get("ni_model_0_m3"))
    ni1 = finite_float(vela.get("ni_model_1_m3"))
    psi0 = finite_float(vela.get("psi0_V"))
    psi1 = finite_float(vela.get("psi1_V"))
    phin0 = finite_float(sentaurus.get("phin0_V"))
    phin1 = finite_float(sentaurus.get("phin1_V"))
    phip0 = finite_float(sentaurus.get("phip0_V"))
    phip1 = finite_float(sentaurus.get("phip1_V"))
    electron_eta = (psi1 - psi0) / vt + math.log(ni1 / ni0) if ni0 > 0.0 and ni1 > 0.0 else None
    hole_eta = (psi1 - psi0) / vt + math.log(ni0 / ni1) if ni0 > 0.0 and ni1 > 0.0 else None
    electron_density_0 = density_factor(ni0, (psi0 - phin0) / vt)
    electron_density_1 = density_factor(ni1, (psi1 - phin1) / vt)
    hole_density_0 = density_factor(ni0, (phip0 - psi0) / vt)
    hole_density_1 = density_factor(ni1, (phip1 - psi1) / vt)

    sent_e_fraction = split_fraction(sent_e, sent_particle)
    sent_h_fraction = split_fraction(sent_h, sent_particle)
    mixed_e_fraction = split_fraction(mixed_e, mixed_particle)
    mixed_h_fraction = split_fraction(mixed_h, mixed_particle)
    return {
        "bias_V": edge_state.get("bias_V", ""),
        "variant": edge_state.get("variant", ""),
        "support_node_id": edge_state.get("support_node_id", ""),
        "support_class": edge_state.get("support_class", ""),
        "y_um": finite_float(edge_state.get("y_um")),
        "edge_id": edge_state.get("edge_id", ""),
        "edge_axis": edge_state.get("edge_axis", ""),
        "node0": vela.get("node0", ""),
        "node1": vela.get("node1", ""),
        "vela_psi_drop_V": finite_float(vela.get("psi1_V")) - finite_float(vela.get("psi0_V")),
        "sentaurus_psi_drop_V": finite_float(sentaurus.get("psi1_V")) - finite_float(sentaurus.get("psi0_V")),
        "mixed_psi_drop_V": psi1 - psi0,
        "mixed_phin_drop_V": phin1 - phin0,
        "mixed_phip_drop_V": phip1 - phip0,
        "mixed_electron_eta": electron_eta,
        "mixed_hole_eta": hole_eta,
        "mixed_electron_bernoulli_plus": sgdiag.bernoulli(electron_eta) if electron_eta is not None else None,
        "mixed_electron_bernoulli_minus": sgdiag.bernoulli(-electron_eta) if electron_eta is not None else None,
        "mixed_hole_bernoulli_plus": sgdiag.bernoulli(hole_eta) if hole_eta is not None else None,
        "mixed_hole_bernoulli_minus": sgdiag.bernoulli(-hole_eta) if hole_eta is not None else None,
        "mixed_electron_density_factor_0": electron_density_0,
        "mixed_electron_density_factor_1": electron_density_1,
        "mixed_hole_density_factor_0": hole_density_0,
        "mixed_hole_density_factor_1": hole_density_1,
        "mixed_electron_density_factor_ratio_1_over_0": safe_ratio(electron_density_1, electron_density_0),
        "mixed_hole_density_factor_ratio_1_over_0": safe_ratio(hole_density_1, hole_density_0),
        "vela_electron_mobility_over_sentaurus": safe_ratio(
            finite_float(vela.get("electron_mobility_m2_V_s")),
            finite_float(sentaurus.get("electron_mobility_m2_V_s")),
        ),
        "vela_hole_mobility_over_sentaurus": safe_ratio(
            finite_float(vela.get("hole_mobility_m2_V_s")),
            finite_float(sentaurus.get("hole_mobility_m2_V_s")),
        ),
        "sentaurus_electron_fraction": sent_e_fraction,
        "sentaurus_hole_fraction": sent_h_fraction,
        "mixed_electron_fraction": mixed_e_fraction,
        "mixed_hole_fraction": mixed_h_fraction,
        "electron_fraction_delta": (
            mixed_e_fraction - sent_e_fraction
            if mixed_e_fraction is not None and sent_e_fraction is not None else None
        ),
        "hole_fraction_delta": (
            mixed_h_fraction - sent_h_fraction
            if mixed_h_fraction is not None and sent_h_fraction is not None else None
        ),
        "mixed_particle_over_sentaurus": safe_ratio(mixed_particle, sent_particle),
        "mixed_electron_over_sentaurus": safe_ratio(mixed_e, sent_e),
        "mixed_hole_over_sentaurus": safe_ratio(mixed_h, sent_h),
        "vela_qf_electron_over_sentaurus": safe_ratio(
            finite_float(vela.get("electron_flux_qf_model_abs")), sent_e),
        "vela_qf_hole_over_sentaurus": safe_ratio(
            finite_float(vela.get("hole_flux_qf_model_abs")), sent_h),
        "sentaurus_qf_electron_over_sentaurus": safe_ratio(
            finite_float(sentaurus.get("electron_flux_qf_model_abs")), sent_e),
        "sentaurus_qf_hole_over_sentaurus": safe_ratio(
            finite_float(sentaurus.get("hole_flux_qf_model_abs")), sent_h),
        "sentaurus_inferred_ni_electron_ratio_1_over_0": safe_ratio(
            finite_float(sentaurus.get("ni_electron_inferred_1_m3")),
            finite_float(sentaurus.get("ni_electron_inferred_0_m3")),
        ),
        "sentaurus_inferred_ni_hole_ratio_1_over_0": safe_ratio(
            finite_float(sentaurus.get("ni_hole_inferred_1_m3")),
            finite_float(sentaurus.get("ni_hole_inferred_0_m3")),
        ),
    }


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    forms = flux_forms_by_edge(args.sg_flux_form_csv, args.bias)
    rows: list[dict[str, Any]] = []
    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    for edge_state in read_csv(args.edge_state_csv):
        if edge_state.get("variant") != args.variant:
            continue
        if edge_state.get("support_class") != args.support_class:
            continue
        if not matches_bias(edge_state, args.bias):
            continue
        edge_id = edge_state.get("edge_id", "")
        vela = forms.get((edge_id, "vela"))
        sentaurus = forms.get((edge_id, "sentaurus"))
        if vela is None or sentaurus is None:
            continue
        rows.append(build_row(edge_state, vela, sentaurus, vt))
    rows.sort(key=lambda row: (finite_float(row.get("y_um")), int(row["edge_id"])))
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    y = finite_values(rows, "y_um")
    correlations: dict[str, dict[str, Any]] = {}
    for key in CORRELATION_FIELDS:
        values = finite_values(rows, key)
        correlations[key] = {
            "corr_y": pearson(y, values),
            **metric_summary(values),
        }
    axes = sorted({str(row.get("edge_axis", "")) for row in rows if row.get("edge_axis", "")})
    return {
        "row_count": len(rows),
        "support_node_count": len({str(row.get("support_node_id", "")) for row in rows}),
        "edge_count": len({str(row.get("edge_id", "")) for row in rows}),
        "edge_axes": axes,
        "correlations": correlations,
    }


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    args = parse_args()
    rows = build_rows(args)
    if not rows:
        raise SystemExit("no merged edge rows matched the requested filters")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sg_carrier_split_decomposition_edges.csv", rows)
    with open_path(args.out_dir / "sg_carrier_split_decomposition_summary.json", "w") as handle:
        json.dump(clean_json(summarize(rows)), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
