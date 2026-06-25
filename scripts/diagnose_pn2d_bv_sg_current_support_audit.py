#!/usr/bin/env python3
"""Audit PN2D SG current support against Sentaurus edge current density.

This diagnostic keeps the Vela impact-ionization alpha and edge area fixed, then
compares the raw SG particle-flux support against Sentaurus |Jn|/q and |Jp|/q.
It also emits simple unit-factor checks based on mu * carrier * |grad(qF)|.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import pn2d_fixed_state_manifest as fixed_manifest


EDGE_FIELDS = [
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "sent_high_field_p95",
    "sent_high_h_current_p95",
    "edge_length_m",
    "edge_area_proxy_m2",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
    "electron_density_avg_m3",
    "hole_density_avg_m3",
    "psi_drop_V",
    "electron_qf_drop_V",
    "hole_qf_drop_V",
    "electron_qf_grad_abs_V_m",
    "hole_qf_grad_abs_V_m",
    "vela_sg_signed_e_flux_m2_s",
    "vela_sg_signed_h_flux_m2_s",
    "vela_sg_e_flux_m2_s",
    "vela_sg_h_flux_m2_s",
    "sent_edge_e_flux_m2_s",
    "sent_edge_h_flux_m2_s",
    "electron_flux_over_sentaurus",
    "hole_flux_over_sentaurus",
    "electron_linear_qf_flux_m2_s",
    "hole_linear_qf_flux_m2_s",
    "electron_sg_over_linear_qf",
    "hole_sg_over_linear_qf",
    "electron_sentaurus_over_linear_qf",
    "hole_sentaurus_over_linear_qf",
    "vela_e_source_integral",
    "vela_h_source_integral",
    "sent_e_source_integral",
    "sent_h_source_integral",
    "electron_source_over_sentaurus",
    "hole_source_over_sentaurus",
]

SUMMARY_FIELDS = [
    "bias_V",
    "edge_count",
    "high_field_edge_count",
    "vela_e_source_integral",
    "vela_h_source_integral",
    "vela_eh_source_integral",
    "sent_e_source_integral",
    "sent_h_source_integral",
    "sent_eh_source_integral",
    "vela_e_over_sent_e",
    "vela_h_over_sent_h",
    "vela_eh_over_sent_eh",
    "electron_missing_fraction",
    "hole_missing_fraction",
    "electron_flux_ratio_area_weighted",
    "hole_flux_ratio_area_weighted",
    "electron_sg_over_linear_qf_median",
    "hole_sg_over_linear_qf_median",
    "electron_sentaurus_over_linear_qf_median",
    "hole_sentaurus_over_linear_qf_median",
]


def parse_args() -> argparse.Namespace:
    root = Path("build-release/reference_tcad/pn2d_sentaurus2018")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overlap-csv",
        type=Path,
        default=root / "reports/high_field_overlap_compare/high_field_overlap_edges.csv",
    )
    parser.add_argument("--sg-edge-manifest", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, default=Path("build/pn2d_bv_gradqf_full20/mesh.json"))
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-16.0)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def finite_float(raw: Any, default: float = 0.0) -> float:
    value = optional_float(raw)
    return value if value is not None else default


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def bool_cell(raw: Any) -> bool:
    if isinstance(raw, bool):
        return raw
    return str(raw).strip().lower() in {"1", "true", "yes", "y"}


def in_bias_window(bias: float, lo: float, hi: float) -> bool:
    return min(lo, hi) - 1.0e-9 <= bias <= max(lo, hi) + 1.0e-9


def bias_key(bias: float) -> str:
    return f"{bias:.9f}"


def load_overlap_rows(path: Path, bias_lo: float, bias_hi: float) -> dict[tuple[str, int], dict[str, str]]:
    result: dict[tuple[str, int], dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            bias = float(row["bias_V"])
            if in_bias_window(bias, bias_lo, bias_hi):
                result[(bias_key(bias), int(row["edge_id"]))] = row
    return result


def load_nodes(path: Path) -> dict[int, tuple[float, float]]:
    data = read_json(path)
    return {
        int(node["id"]): (float(node["x"]) * 1.0e-6, float(node["y"]) * 1.0e-6)
        for node in data.get("nodes", [])
    }


def load_state(path: Path) -> dict[int, dict[str, float]]:
    rows = read_csv(path)
    state: dict[int, dict[str, float]] = {}
    for row in rows:
        state[int(row["node_id"])] = {
            "psi": finite_float(row.get("psi")),
            "phin": finite_float(row.get("phin")),
            "phip": finite_float(row.get("phip")),
            "electrons_m3": finite_float(row.get("electrons_m3")),
            "holes_m3": finite_float(row.get("holes_m3")),
        }
    return state


def state_path_from_manifest_item(item: dict[str, Any]) -> Path | None:
    if item.get("state_csv"):
        return Path(item["state_csv"])
    config_path = item.get("config")
    if not config_path:
        return None
    config = read_json(Path(config_path))
    sweep = config.get("sweep") or {}
    raw = sweep.get("write_state_file")
    return Path(raw) if raw else None


def load_sg_edges_and_states(path: Path, bias_lo: float, bias_hi: float) -> tuple[dict[tuple[str, int], dict[str, str]], dict[str, dict[int, dict[str, float]]]]:
    manifest = read_json(path)
    sg_edges: dict[tuple[str, int], dict[str, str]] = {}
    states: dict[str, dict[int, dict[str, float]]] = {}
    for item in manifest:
        bias = float(item["bias_V"])
        if not in_bias_window(bias, bias_lo, bias_hi):
            continue
        bkey = bias_key(bias)
        fixed_manifest.require_materials_file(item, path)
        for row in read_csv(Path(item["sg_edges_csv"])):
            sg_edges[(bias_key(float(row["bias_V"])), int(row["edge_id"]))] = row
        state_path = state_path_from_manifest_item(item)
        if state_path is not None and state_path.exists():
            states[bkey] = load_state(state_path)
    return sg_edges, states


def edge_length(row: dict[str, str], nodes: dict[int, tuple[float, float]]) -> float:
    length = optional_float(row.get("edge_length_m"))
    if length is not None and length > 0.0:
        return length
    n0 = int(row["node0"])
    n1 = int(row["node1"])
    if n0 in nodes and n1 in nodes:
        x0, y0 = nodes[n0]
        x1, y1 = nodes[n1]
        return math.hypot(x1 - x0, y1 - y0)
    return 0.0


def raw_signed_flux(row: dict[str, str], carrier: str) -> float:
    signed = optional_float(row.get(f"{carrier}_raw_signed_flux_proxy"))
    if signed is not None:
        return signed
    raw = optional_float(row.get(f"{carrier}_raw_flux_proxy"))
    if raw is not None:
        return raw
    return finite_float(row.get(f"{carrier}_flux_proxy"))


def endpoint_state(states: dict[int, dict[str, float]], node0: int, node1: int, key: str) -> tuple[float | None, float | None]:
    s0 = states.get(node0)
    s1 = states.get(node1)
    if s0 is None or s1 is None:
        return None, None
    return s0.get(key), s1.get(key)


def avg_pair(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return 0.5 * (a + b)


def abs_drop(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return b - a


def linear_qf_flux(mobility: float, density_avg: float | None, qf_drop: float | None, length: float) -> float | None:
    if density_avg is None or qf_drop is None or mobility <= 0.0 or length <= 0.0:
        return None
    return mobility * max(density_avg, 0.0) * abs(qf_drop) / length


def make_edge_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    overlap = load_overlap_rows(args.overlap_csv, args.bias_min, args.bias_max)
    sg_edges, states_by_bias = load_sg_edges_and_states(args.sg_edge_manifest, args.bias_min, args.bias_max)
    nodes = load_nodes(args.mesh)
    rows: list[dict[str, Any]] = []
    for key, sg in sorted(sg_edges.items(), key=lambda item: (float(item[0][0]), item[0][1])):
        bkey, edge_id = key
        ov = overlap.get(key)
        if ov is None:
            continue
        bias = float(sg.get("bias_V", ov["bias_V"]))
        node0 = int(sg.get("node0", ov["node0"]))
        node1 = int(sg.get("node1", ov["node1"]))
        state = states_by_bias.get(bkey, {})
        length = edge_length(sg, nodes)
        edge_area = finite_float(sg.get("edge_area_proxy_m2"), finite_float(ov.get("edge_area_proxy_m2")))
        alpha_e = finite_float(sg.get("electron_alpha_m_inv"))
        alpha_h = finite_float(sg.get("hole_alpha_m_inv"))
        mu_e = finite_float(sg.get("electron_mobility_m2_V_s"))
        mu_h = finite_float(sg.get("hole_mobility_m2_V_s"))
        e_signed = raw_signed_flux(sg, "electron")
        h_signed = raw_signed_flux(sg, "hole")
        e_flux = abs(e_signed)
        h_flux = abs(h_signed)
        sent_e_flux = abs(finite_float(ov.get("sent_e_particle_flux_m2_s")))
        sent_h_flux = abs(finite_float(ov.get("sent_h_particle_flux_m2_s")))

        psi0, psi1 = endpoint_state(state, node0, node1, "psi")
        phin0, phin1 = endpoint_state(state, node0, node1, "phin")
        phip0, phip1 = endpoint_state(state, node0, node1, "phip")
        n0, n1 = endpoint_state(state, node0, node1, "electrons_m3")
        p0, p1 = endpoint_state(state, node0, node1, "holes_m3")
        psi_drop = abs_drop(psi0, psi1)
        e_qf_drop = abs_drop(phin0, phin1)
        h_qf_drop = abs_drop(phip0, phip1)
        n_avg = avg_pair(n0, n1)
        p_avg = avg_pair(p0, p1)
        e_linear = linear_qf_flux(mu_e, n_avg, e_qf_drop, length)
        h_linear = linear_qf_flux(mu_h, p_avg, h_qf_drop, length)

        vela_e_source = alpha_e * e_flux * edge_area
        vela_h_source = alpha_h * h_flux * edge_area
        sent_e_source = alpha_e * sent_e_flux * edge_area
        sent_h_source = alpha_h * sent_h_flux * edge_area
        rows.append({
            "bias_V": bias,
            "edge_id": edge_id,
            "node0": node0,
            "node1": node1,
            "sent_high_field_p95": bool_cell(ov.get("sent_high_field_p95")),
            "sent_high_h_current_p95": bool_cell(ov.get("sent_high_h_current_p95")),
            "edge_length_m": length,
            "edge_area_proxy_m2": edge_area,
            "electron_alpha_m_inv": alpha_e,
            "hole_alpha_m_inv": alpha_h,
            "electron_mobility_m2_V_s": mu_e,
            "hole_mobility_m2_V_s": mu_h,
            "electron_density_avg_m3": n_avg,
            "hole_density_avg_m3": p_avg,
            "psi_drop_V": psi_drop,
            "electron_qf_drop_V": e_qf_drop,
            "hole_qf_drop_V": h_qf_drop,
            "electron_qf_grad_abs_V_m": abs(e_qf_drop) / length if e_qf_drop is not None and length > 0.0 else None,
            "hole_qf_grad_abs_V_m": abs(h_qf_drop) / length if h_qf_drop is not None and length > 0.0 else None,
            "vela_sg_signed_e_flux_m2_s": e_signed,
            "vela_sg_signed_h_flux_m2_s": h_signed,
            "vela_sg_e_flux_m2_s": e_flux,
            "vela_sg_h_flux_m2_s": h_flux,
            "sent_edge_e_flux_m2_s": sent_e_flux,
            "sent_edge_h_flux_m2_s": sent_h_flux,
            "electron_flux_over_sentaurus": ratio(e_flux, sent_e_flux),
            "hole_flux_over_sentaurus": ratio(h_flux, sent_h_flux),
            "electron_linear_qf_flux_m2_s": e_linear,
            "hole_linear_qf_flux_m2_s": h_linear,
            "electron_sg_over_linear_qf": ratio(e_flux, e_linear),
            "hole_sg_over_linear_qf": ratio(h_flux, h_linear),
            "electron_sentaurus_over_linear_qf": ratio(sent_e_flux, e_linear),
            "hole_sentaurus_over_linear_qf": ratio(sent_h_flux, h_linear),
            "vela_e_source_integral": vela_e_source,
            "vela_h_source_integral": vela_h_source,
            "sent_e_source_integral": sent_e_source,
            "sent_h_source_integral": sent_h_source,
            "electron_source_over_sentaurus": ratio(vela_e_source, sent_e_source),
            "hole_source_over_sentaurus": ratio(vela_h_source, sent_h_source),
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def weighted_flux_ratio(rows: list[dict[str, Any]], carrier: str) -> float | None:
    alpha_key = f"{carrier}_alpha_m_inv"
    area_key = "edge_area_proxy_m2"
    vela_key = f"vela_sg_{'e' if carrier == 'electron' else 'h'}_flux_m2_s"
    sent_key = f"sent_edge_{'e' if carrier == 'electron' else 'h'}_flux_m2_s"
    numerator = sum(float(row[alpha_key]) * float(row[area_key]) * float(row[vela_key]) for row in rows)
    denominator = sum(float(row[alpha_key]) * float(row[area_key]) * float(row[sent_key]) for row in rows)
    return ratio(numerator, denominator)


def summarize_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for bias in sorted({float(row["bias_V"]) for row in rows}):
        subset = [row for row in rows if abs(float(row["bias_V"]) - bias) <= 1.0e-9]
        vela_e = sum(float(row["vela_e_source_integral"]) for row in subset)
        vela_h = sum(float(row["vela_h_source_integral"]) for row in subset)
        sent_e = sum(float(row["sent_e_source_integral"]) for row in subset)
        sent_h = sum(float(row["sent_h_source_integral"]) for row in subset)
        missing_e = sent_e - vela_e
        missing_h = sent_h - vela_h
        missing_total = missing_e + missing_h
        result.append({
            "bias_V": bias,
            "edge_count": len(subset),
            "high_field_edge_count": sum(1 for row in subset if row["sent_high_field_p95"]),
            "vela_e_source_integral": vela_e,
            "vela_h_source_integral": vela_h,
            "vela_eh_source_integral": vela_e + vela_h,
            "sent_e_source_integral": sent_e,
            "sent_h_source_integral": sent_h,
            "sent_eh_source_integral": sent_e + sent_h,
            "vela_e_over_sent_e": ratio(vela_e, sent_e),
            "vela_h_over_sent_h": ratio(vela_h, sent_h),
            "vela_eh_over_sent_eh": ratio(vela_e + vela_h, sent_e + sent_h),
            "electron_missing_fraction": ratio(missing_e, missing_total),
            "hole_missing_fraction": ratio(missing_h, missing_total),
            "electron_flux_ratio_area_weighted": weighted_flux_ratio(subset, "electron"),
            "hole_flux_ratio_area_weighted": weighted_flux_ratio(subset, "hole"),
            "electron_sg_over_linear_qf_median": statistics.median(finite_values(subset, "electron_sg_over_linear_qf")) if finite_values(subset, "electron_sg_over_linear_qf") else None,
            "hole_sg_over_linear_qf_median": statistics.median(finite_values(subset, "hole_sg_over_linear_qf")) if finite_values(subset, "hole_sg_over_linear_qf") else None,
            "electron_sentaurus_over_linear_qf_median": statistics.median(finite_values(subset, "electron_sentaurus_over_linear_qf")) if finite_values(subset, "electron_sentaurus_over_linear_qf") else None,
            "hole_sentaurus_over_linear_qf_median": statistics.median(finite_values(subset, "hole_sentaurus_over_linear_qf")) if finite_values(subset, "hole_sentaurus_over_linear_qf") else None,
        })
    return result


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
    rows = make_edge_rows(args)
    summary = summarize_rows(rows)
    high_field_rows = [row for row in rows if row["sent_high_field_p95"]]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sg_current_support_edge_audit.csv", EDGE_FIELDS, rows)
    write_csv(args.out_dir / "sg_current_support_high_field_edges.csv", EDGE_FIELDS, high_field_rows)
    write_csv(args.out_dir / "sg_current_support_summary.csv", SUMMARY_FIELDS, summary)
    payload = {
        "schema": "pn2d.sg_current_support_audit.v1",
        "edge_row_count": len(rows),
        "high_field_edge_row_count": len(high_field_rows),
        "bias_count": len(summary),
        "summary": summary,
    }
    (args.out_dir / "sg_current_support_summary.json").write_text(
        json.dumps(clean_json(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "bias_count": len(summary),
        "edge_row_count": len(rows),
        "high_field_edge_row_count": len(high_field_rows),
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
