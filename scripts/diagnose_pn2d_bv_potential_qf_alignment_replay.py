#!/usr/bin/env python3
"""Replay PN2D BV SG currents with swapped potential, qF, and mobility sources."""

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
    "input_variant",
    "replay_variant",
    "psi_shift_V",
    "support_node_id",
    "support_class",
    "y_um",
    "edge_id",
    "edge_axis",
    "node0",
    "node1",
    "psi_source",
    "qf_source",
    "mobility_source",
    "ni_source",
    "psi0_V",
    "psi1_V",
    "phin0_V",
    "phin1_V",
    "phip0_V",
    "phip1_V",
    "ni0_m3",
    "ni1_m3",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
    "electron_flux_m2_s",
    "hole_flux_m2_s",
    "particle_flux_m2_s",
    "sentaurus_electron_flux_m2_s",
    "sentaurus_hole_flux_m2_s",
    "sentaurus_particle_flux_m2_s",
    "electron_over_sentaurus",
    "hole_over_sentaurus",
    "particle_over_sentaurus",
    "sentaurus_electron_fraction",
    "sentaurus_hole_fraction",
    "electron_fraction",
    "hole_fraction",
    "electron_fraction_delta",
    "hole_fraction_delta",
]

NODE_SHIFT_FIELDS = [
    "support_node_id",
    "support_class",
    "y_um",
    "edge_count",
    "best_psi_shift_V",
    "combined_fraction_l2",
    "direct_sentaurus_minus_vela_psi_shift_median_V",
    "direct_sentaurus_minus_vela_psi_shift_min_V",
    "direct_sentaurus_minus_vela_psi_shift_max_V",
]

CORRELATION_FIELDS = [
    "electron_over_sentaurus",
    "hole_over_sentaurus",
    "particle_over_sentaurus",
    "electron_fraction_delta",
    "hole_fraction_delta",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-state-csv", type=Path, required=True)
    parser.add_argument("--sg-flux-form-csv", type=Path, required=True)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--support-class", default="overlap")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--shift-min", type=float, default=-0.2)
    parser.add_argument("--shift-max", type=float, default=0.2)
    parser.add_argument("--shift-steps", type=int, default=801)
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


def write_node_shift_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_SHIFT_FIELDS, extrasaction="ignore")
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
        if matches_bias(row, bias):
            result[(row.get("edge_id", ""), row.get("source", ""))] = row
    return result


def filtered_edge_states(args: argparse.Namespace) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_csv(args.edge_state_csv):
        if row.get("variant") != args.variant:
            continue
        if row.get("support_class") != args.support_class:
            continue
        if not matches_bias(row, args.bias):
            continue
        rows.append(row)
    return rows


def source_rows(
    edge_state: dict[str, str],
    forms: dict[tuple[str, str], dict[str, str]],
) -> tuple[dict[str, str], dict[str, str]] | None:
    edge_id = edge_state.get("edge_id", "")
    vela = forms.get((edge_id, "vela"))
    sentaurus = forms.get((edge_id, "sentaurus"))
    if vela is None or sentaurus is None:
        return None
    return vela, sentaurus


def pick(source: str, vela: dict[str, str], sentaurus: dict[str, str]) -> dict[str, str]:
    return sentaurus if source == "sentaurus" else vela


def component_fluxes(
    psi: dict[str, str],
    qf: dict[str, str],
    mobility: dict[str, str],
    ni: dict[str, str],
    vt: float,
    psi_shift: float = 0.0,
) -> tuple[float, float, float]:
    length = finite_float(mobility.get("length_m"))
    if length <= 0.0:
        raise ValueError("edge length_m must be positive")
    psi0 = finite_float(psi.get("psi0_V")) + psi_shift
    psi1 = finite_float(psi.get("psi1_V")) + psi_shift
    phin0 = finite_float(qf.get("phin0_V"))
    phin1 = finite_float(qf.get("phin1_V"))
    phip0 = finite_float(qf.get("phip0_V"))
    phip1 = finite_float(qf.get("phip1_V"))
    ni0 = finite_float(ni.get("ni_model_0_m3"))
    ni1 = finite_float(ni.get("ni_model_1_m3"))
    mun = finite_float(mobility.get("electron_mobility_m2_V_s"))
    mup = finite_float(mobility.get("hole_mobility_m2_V_s"))
    electron = abs(sgdiag.sg_electron_flux_qf_variable_ni(
        ni0, ni1, psi0, psi1, phin0, phin1, vt, mun * vt / length))
    hole = abs(sgdiag.sg_hole_flux_qf_variable_ni(
        ni0, ni1, psi0, psi1, phip0, phip1, vt, mup * vt / length))
    return electron, hole, electron + hole


def split_fraction(component: float, particle: float) -> float | None:
    return safe_ratio(component, particle)


def build_replay_row(
    edge_state: dict[str, str],
    vela: dict[str, str],
    sentaurus: dict[str, str],
    vt: float,
    replay_variant: str,
    psi_source: str,
    qf_source: str,
    mobility_source: str,
    ni_source: str,
    psi_shift: float = 0.0,
) -> dict[str, Any]:
    psi = pick(psi_source, vela, sentaurus)
    qf = pick(qf_source, vela, sentaurus)
    mobility = pick(mobility_source, vela, sentaurus)
    ni = pick(ni_source, vela, sentaurus)
    electron, hole, particle = component_fluxes(psi, qf, mobility, ni, vt, psi_shift)
    sent_e = finite_float(edge_state.get("sentaurus_edgeavg_electron_current_flux_m2_s"))
    sent_h = finite_float(edge_state.get("sentaurus_edgeavg_hole_current_flux_m2_s"))
    sent_particle = finite_float(edge_state.get("sentaurus_edgeavg_particle_flux_m2_s"))
    sent_e_fraction = split_fraction(sent_e, sent_particle)
    sent_h_fraction = split_fraction(sent_h, sent_particle)
    electron_fraction = split_fraction(electron, particle)
    hole_fraction = split_fraction(hole, particle)
    psi0 = finite_float(psi.get("psi0_V")) + psi_shift
    psi1 = finite_float(psi.get("psi1_V")) + psi_shift
    return {
        "bias_V": edge_state.get("bias_V", ""),
        "input_variant": edge_state.get("variant", ""),
        "replay_variant": replay_variant,
        "psi_shift_V": psi_shift,
        "support_node_id": edge_state.get("support_node_id", ""),
        "support_class": edge_state.get("support_class", ""),
        "y_um": finite_float(edge_state.get("y_um")),
        "edge_id": edge_state.get("edge_id", ""),
        "edge_axis": edge_state.get("edge_axis", ""),
        "node0": psi.get("node0", ""),
        "node1": psi.get("node1", ""),
        "psi_source": psi_source,
        "qf_source": qf_source,
        "mobility_source": mobility_source,
        "ni_source": ni_source,
        "psi0_V": psi0,
        "psi1_V": psi1,
        "phin0_V": finite_float(qf.get("phin0_V")),
        "phin1_V": finite_float(qf.get("phin1_V")),
        "phip0_V": finite_float(qf.get("phip0_V")),
        "phip1_V": finite_float(qf.get("phip1_V")),
        "ni0_m3": finite_float(ni.get("ni_model_0_m3")),
        "ni1_m3": finite_float(ni.get("ni_model_1_m3")),
        "electron_mobility_m2_V_s": finite_float(mobility.get("electron_mobility_m2_V_s")),
        "hole_mobility_m2_V_s": finite_float(mobility.get("hole_mobility_m2_V_s")),
        "electron_flux_m2_s": electron,
        "hole_flux_m2_s": hole,
        "particle_flux_m2_s": particle,
        "sentaurus_electron_flux_m2_s": sent_e,
        "sentaurus_hole_flux_m2_s": sent_h,
        "sentaurus_particle_flux_m2_s": sent_particle,
        "electron_over_sentaurus": safe_ratio(electron, sent_e),
        "hole_over_sentaurus": safe_ratio(hole, sent_h),
        "particle_over_sentaurus": safe_ratio(particle, sent_particle),
        "sentaurus_electron_fraction": sent_e_fraction,
        "sentaurus_hole_fraction": sent_h_fraction,
        "electron_fraction": electron_fraction,
        "hole_fraction": hole_fraction,
        "electron_fraction_delta": (
            electron_fraction - sent_e_fraction
            if electron_fraction is not None and sent_e_fraction is not None else None
        ),
        "hole_fraction_delta": (
            hole_fraction - sent_h_fraction
            if hole_fraction is not None and sent_h_fraction is not None else None
        ),
    }


REPLAY_VARIANTS = [
    (
        "vela_psi_sentaurus_qf_vela_mobility",
        "vela",
        "sentaurus",
        "vela",
        "vela",
        0.0,
    ),
    (
        "sentaurus_psi_sentaurus_qf_vela_mobility",
        "sentaurus",
        "sentaurus",
        "vela",
        "vela",
        0.0,
    ),
    (
        "vela_psi_sentaurus_qf_sentaurus_mobility",
        "vela",
        "sentaurus",
        "sentaurus",
        "vela",
        0.0,
    ),
    (
        "sentaurus_psi_sentaurus_qf_sentaurus_mobility",
        "sentaurus",
        "sentaurus",
        "sentaurus",
        "sentaurus",
        0.0,
    ),
]


def replay_rows_for_shift(
    matched: list[tuple[dict[str, str], dict[str, str], dict[str, str]]],
    vt: float,
    shift: float,
) -> list[dict[str, Any]]:
    return [
        build_replay_row(
            edge_state,
            vela,
            sentaurus,
            vt,
            "vela_psi_shifted_const_sentaurus_qf_vela_mobility",
            "vela",
            "sentaurus",
            "vela",
            "vela",
            shift,
        )
        for edge_state, vela, sentaurus in matched
    ]


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
    return {"median": statistics.median(values), "min": min(values), "max": max(values)}


def combined_fraction_l2(rows: list[dict[str, Any]]) -> float | None:
    terms: list[float] = []
    for row in rows:
        e_delta = optional_float(row.get("electron_fraction_delta"))
        h_delta = optional_float(row.get("hole_fraction_delta"))
        if e_delta is not None and h_delta is not None:
            terms.append(e_delta * e_delta + h_delta * h_delta)
    if not terms:
        return None
    return math.sqrt(statistics.fmean(terms))


def summarize_variant(rows: list[dict[str, Any]]) -> dict[str, Any]:
    y = finite_values(rows, "y_um")
    correlations: dict[str, dict[str, Any]] = {}
    for key in CORRELATION_FIELDS:
        values = finite_values(rows, key)
        correlations[key] = {"corr_y": pearson(y, values), **metric_summary(values)}
    shifts = sorted({finite_float(row.get("psi_shift_V")) for row in rows})
    return {
        "row_count": len(rows),
        "psi_shift_V": shifts[0] if len(shifts) == 1 else None,
        "combined_fraction_l2": combined_fraction_l2(rows),
        "correlations": correlations,
    }


def shift_candidates(minimum: float, maximum: float, steps: int) -> list[float]:
    if steps < 2:
        return [minimum]
    step = (maximum - minimum) / float(steps - 1)
    return [minimum + index * step for index in range(steps)]


def best_constant_shift(
    matched: list[tuple[dict[str, str], dict[str, str], dict[str, str]]],
    vt: float,
    minimum: float,
    maximum: float,
    steps: int,
) -> tuple[float, list[dict[str, Any]]]:
    best_shift = minimum
    best_rows = replay_rows_for_shift(matched, vt, best_shift)
    best_l2 = combined_fraction_l2(best_rows)
    for shift in shift_candidates(minimum, maximum, steps):
        rows = replay_rows_for_shift(matched, vt, shift)
        l2 = combined_fraction_l2(rows)
        if l2 is not None and (best_l2 is None or l2 < best_l2):
            best_shift = shift
            best_rows = rows
            best_l2 = l2
    return best_shift, best_rows


def matched_edge_inputs(args: argparse.Namespace) -> list[tuple[dict[str, str], dict[str, str], dict[str, str]]]:
    forms = flux_forms_by_edge(args.sg_flux_form_csv, args.bias)
    matched: list[tuple[dict[str, str], dict[str, str], dict[str, str]]] = []
    for edge_state in filtered_edge_states(args):
        pair = source_rows(edge_state, forms)
        if pair is not None:
            vela, sentaurus = pair
            matched.append((edge_state, vela, sentaurus))
    return matched


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    matched = matched_edge_inputs(args)
    if not matched:
        return []
    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    rows: list[dict[str, Any]] = []
    for edge_state, vela, sentaurus in matched:
        for replay in REPLAY_VARIANTS:
            rows.append(build_replay_row(edge_state, vela, sentaurus, vt, *replay))
    _best_shift, shifted_rows = best_constant_shift(
        matched, vt, args.shift_min, args.shift_max, args.shift_steps)
    rows.extend(shifted_rows)
    rows.sort(key=lambda row: (
        str(row.get("replay_variant", "")),
        finite_float(row.get("y_um")),
        int(str(row.get("edge_id", "0"))),
    ))
    return rows


def direct_psi_shifts(vela: dict[str, str], sentaurus: dict[str, str]) -> list[float]:
    return [
        finite_float(sentaurus.get("psi0_V")) - finite_float(vela.get("psi0_V")),
        finite_float(sentaurus.get("psi1_V")) - finite_float(vela.get("psi1_V")),
    ]


def build_node_shift_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    matched = matched_edge_inputs(args)
    if not matched:
        return []
    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    by_node: dict[str, list[tuple[dict[str, str], dict[str, str], dict[str, str]]]] = {}
    for item in matched:
        edge_state, _vela, _sentaurus = item
        by_node.setdefault(edge_state.get("support_node_id", ""), []).append(item)
    rows: list[dict[str, Any]] = []
    for node_id, node_items in by_node.items():
        best_shift, best_rows = best_constant_shift(
            node_items, vt, args.shift_min, args.shift_max, args.shift_steps)
        direct_shifts: list[float] = []
        for _edge_state, vela, sentaurus in node_items:
            direct_shifts.extend(direct_psi_shifts(vela, sentaurus))
        y_values = [
            finite_float(edge_state.get("y_um"))
            for edge_state, _vela, _sentaurus in node_items
        ]
        support_classes = {
            edge_state.get("support_class", "")
            for edge_state, _vela, _sentaurus in node_items
        }
        rows.append({
            "support_node_id": node_id,
            "support_class": ",".join(sorted(support_classes)),
            "y_um": statistics.median(y_values) if y_values else None,
            "edge_count": len(node_items),
            "best_psi_shift_V": best_shift,
            "combined_fraction_l2": combined_fraction_l2(best_rows),
            "direct_sentaurus_minus_vela_psi_shift_median_V": (
                statistics.median(direct_shifts) if direct_shifts else None
            ),
            "direct_sentaurus_minus_vela_psi_shift_min_V": min(direct_shifts) if direct_shifts else None,
            "direct_sentaurus_minus_vela_psi_shift_max_V": max(direct_shifts) if direct_shifts else None,
        })
    rows.sort(key=lambda row: (finite_float(row.get("y_um")), str(row.get("support_node_id", ""))))
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_variant: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_variant.setdefault(str(row.get("replay_variant", "")), []).append(row)
    return {
        "row_count": len(rows),
        "edge_count": len({str(row.get("edge_id", "")) for row in rows}),
        "support_node_count": len({str(row.get("support_node_id", "")) for row in rows}),
        "replay_variant_count": len(by_variant),
        "variants": {
            name: summarize_variant(variant_rows)
            for name, variant_rows in sorted(by_variant.items())
        },
    }


def summarize_node_shifts(rows: list[dict[str, Any]]) -> dict[str, Any]:
    y = finite_values(rows, "y_um")
    return {
        "row_count": len(rows),
        "best_psi_shift_V": {
            "corr_y": pearson(y, finite_values(rows, "best_psi_shift_V")),
            **metric_summary(finite_values(rows, "best_psi_shift_V")),
        },
        "direct_sentaurus_minus_vela_psi_shift_median_V": {
            "corr_y": pearson(y, finite_values(rows, "direct_sentaurus_minus_vela_psi_shift_median_V")),
            **metric_summary(finite_values(rows, "direct_sentaurus_minus_vela_psi_shift_median_V")),
        },
        "combined_fraction_l2": metric_summary(finite_values(rows, "combined_fraction_l2")),
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
        raise SystemExit("no replay edge rows matched the requested filters")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "potential_qf_alignment_replay_edges.csv", rows)
    node_shift_rows = build_node_shift_rows(args)
    write_node_shift_csv(args.out_dir / "potential_qf_alignment_node_shift_scan.csv", node_shift_rows)
    with open_path(args.out_dir / "potential_qf_alignment_replay_summary.json", "w") as handle:
        summary = summarize(rows)
        summary["node_shift_scan"] = summarize_node_shifts(node_shift_rows)
        json.dump(clean_json(summary), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
