#!/usr/bin/env python3
"""Matched Vela/Sentaurus current-support factor audit for PN2D BV.

For each Vela fixed-state SG edge, match its endpoints to nearest Sentaurus
nodes at the same bias and decompose the local linear current-support ratio into
carrier density, mobility, and quasi-Fermi-gradient factors.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any, Iterable

import pn2d_fixed_state_manifest as fixed_manifest


EDGE_FIELDS = [
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "sentaurus_node0",
    "sentaurus_node1",
    "sentaurus_distance0_um",
    "sentaurus_distance1_um",
    "sent_high_field_p95",
    "edge_area_proxy_m2",
    "edge_length_m",
    "electron_density_ratio_vela_over_sentaurus",
    "electron_mobility_ratio_vela_over_sentaurus",
    "electron_qf_grad_ratio_vela_over_sentaurus",
    "electron_linear_support_ratio_vela_over_sentaurus",
    "hole_density_ratio_vela_over_sentaurus",
    "hole_mobility_ratio_vela_over_sentaurus",
    "hole_qf_grad_ratio_vela_over_sentaurus",
    "hole_linear_support_ratio_vela_over_sentaurus",
    "vela_e_density_avg_cm3",
    "sent_e_density_avg_cm3",
    "vela_h_density_avg_cm3",
    "sent_h_density_avg_cm3",
    "vela_e_mobility_cm2_V_s",
    "sent_e_mobility_avg_cm2_V_s",
    "vela_h_mobility_cm2_V_s",
    "sent_h_mobility_avg_cm2_V_s",
    "vela_e_qf_grad_abs_V_m",
    "sent_e_qf_grad_abs_V_m",
    "vela_h_qf_grad_abs_V_m",
    "sent_h_qf_grad_abs_V_m",
    "vela_sg_e_flux_m2_s",
    "vela_sg_h_flux_m2_s",
    "sent_e_flux_m2_s",
    "sent_h_flux_m2_s",
    "vela_e_flux_over_sentaurus_current",
    "vela_h_flux_over_sentaurus_current",
]

SUMMARY_FIELDS = [
    "bias_V",
    "edge_count",
    "high_field_edge_count",
    "electron_density_ratio_median",
    "electron_mobility_ratio_median",
    "electron_qf_grad_ratio_median",
    "electron_linear_support_ratio_median",
    "hole_density_ratio_median",
    "hole_mobility_ratio_median",
    "hole_qf_grad_ratio_median",
    "hole_linear_support_ratio_median",
    "electron_density_ratio_high_field_median",
    "electron_mobility_ratio_high_field_median",
    "electron_qf_grad_ratio_high_field_median",
    "electron_linear_support_ratio_high_field_median",
    "hole_density_ratio_high_field_median",
    "hole_mobility_ratio_high_field_median",
    "hole_qf_grad_ratio_high_field_median",
    "hole_linear_support_ratio_high_field_median",
    "electron_linear_support_ratio_source_weighted",
    "hole_linear_support_ratio_source_weighted",
    "electron_flux_over_sentaurus_current_source_weighted",
    "hole_flux_over_sentaurus_current_source_weighted",
]

Q_C = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    root = Path("build-release/reference_tcad/pn2d_sentaurus2018")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlap-csv", type=Path, default=root / "reports/high_field_overlap_compare/high_field_overlap_edges.csv")
    parser.add_argument("--sg-edge-manifest", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, default=root / "sentaurus_multibias_0p25v")
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


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
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


def in_bias_window(bias: float, lo: float, hi: float) -> bool:
    return min(lo, hi) - 1.0e-9 <= bias <= max(lo, hi) + 1.0e-9


def bias_key(bias: float) -> str:
    return f"{bias:.9f}"


def sentaurus_bias_dir(root: Path, bias: float) -> Path:
    text = f"{abs(bias):g}".replace(".", "p")
    return root / f"sentaurus_-{text}v"


def load_mesh_nodes(path: Path) -> dict[int, tuple[float, float]]:
    data = read_json(path)
    return {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in data.get("nodes", [])
    }


def load_sentaurus_nodes(case_dir: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(case_dir / "nodes.csv")
    }


def nearest_node(nodes: dict[int, tuple[float, float]], x_um: float, y_um: float) -> tuple[int, float]:
    best_id = -1
    best_d2 = math.inf
    for node_id, (x, y) in nodes.items():
        d2 = (x - x_um) * (x - x_um) + (y - y_um) * (y - y_um)
        if d2 < best_d2:
            best_id = node_id
            best_d2 = d2
    return best_id, math.sqrt(best_d2)


def numeric_components(row: dict[str, str]) -> list[float]:
    values: list[float] = []
    for key, raw in row.items():
        if key == "node_id":
            continue
        value = optional_float(raw)
        if value is not None:
            values.append(value)
    return values


def load_scalar_field(fields_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted(fields_dir.glob(f"{name}_region*.csv")):
        for row in read_csv(path):
            values = numeric_components(row)
            result[int(row["node_id"])] = values[0] if values else 0.0
    return result


def avg(field: dict[int, float], a: int, b: int) -> float | None:
    if a not in field or b not in field:
        return None
    return 0.5 * (field[a] + field[b])


def drop(field: dict[int, float], a: int, b: int) -> float | None:
    if a not in field or b not in field:
        return None
    return field[b] - field[a]


def state_path_from_manifest_item(item: dict[str, Any]) -> Path | None:
    if item.get("state_csv"):
        return Path(item["state_csv"])
    config_path = item.get("config")
    if not config_path:
        return None
    config = read_json(Path(config_path))
    raw = (config.get("sweep") or {}).get("write_state_file")
    return Path(raw) if raw else None


def load_state(path: Path) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    for row in read_csv(path):
        result[int(row["node_id"])] = {
            "phin": finite_float(row.get("phin")),
            "phip": finite_float(row.get("phip")),
            "electrons_m3": finite_float(row.get("electrons_m3")),
            "holes_m3": finite_float(row.get("holes_m3")),
        }
    return result


def load_manifest(path: Path, bias_lo: float, bias_hi: float) -> tuple[dict[tuple[str, int], dict[str, str]], dict[str, dict[int, dict[str, float]]]]:
    sg_edges: dict[tuple[str, int], dict[str, str]] = {}
    states: dict[str, dict[int, dict[str, float]]] = {}
    for item in read_json(path):
        bias = float(item["bias_V"])
        if not in_bias_window(bias, bias_lo, bias_hi):
            continue
        bkey = bias_key(bias)
        fixed_manifest.require_materials_file(item, path)
        for row in read_csv(Path(item["sg_edges_csv"])):
            sg_edges[(bias_key(float(row["bias_V"])), int(row["edge_id"]))] = row
        state_path = state_path_from_manifest_item(item)
        if state_path and state_path.exists():
            states[bkey] = load_state(state_path)
    return sg_edges, states


def load_overlap(path: Path, bias_lo: float, bias_hi: float) -> dict[tuple[str, int], dict[str, str]]:
    result: dict[tuple[str, int], dict[str, str]] = {}
    for row in read_csv(path):
        bias = float(row["bias_V"])
        if in_bias_window(bias, bias_lo, bias_hi):
            result[(bias_key(bias), int(row["edge_id"]))] = row
    return result


def load_sentaurus_case(case_dir: Path) -> tuple[dict[int, tuple[float, float]], dict[str, dict[int, float]]]:
    fields = case_dir / "fields"
    return load_sentaurus_nodes(case_dir), {
        "eDensity": load_scalar_field(fields, "eDensity"),
        "hDensity": load_scalar_field(fields, "hDensity"),
        "eMobility": load_scalar_field(fields, "eMobility"),
        "hMobility": load_scalar_field(fields, "hMobility"),
        "eQF": load_scalar_field(fields, "eQuasiFermiPotential"),
        "hQF": load_scalar_field(fields, "hQuasiFermiPotential"),
    }


def raw_flux(row: dict[str, str], carrier: str) -> float:
    signed = optional_float(row.get(f"{carrier}_raw_signed_flux_proxy"))
    if signed is not None:
        return abs(signed)
    raw = optional_float(row.get(f"{carrier}_raw_flux_proxy"))
    if raw is not None:
        return abs(raw)
    return abs(finite_float(row.get(f"{carrier}_flux_proxy")))


def endpoint_avg_state(state: dict[int, dict[str, float]], a: int, b: int, key: str) -> float | None:
    if a not in state or b not in state:
        return None
    return 0.5 * (state[a][key] + state[b][key])


def endpoint_drop_state(state: dict[int, dict[str, float]], a: int, b: int, key: str) -> float | None:
    if a not in state or b not in state:
        return None
    return state[b][key] - state[a][key]


def qf_grad(drop_v: float | None, length_m: float) -> float | None:
    if drop_v is None or length_m <= 0.0:
        return None
    return abs(drop_v) / length_m


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    mesh_nodes = load_mesh_nodes(args.mesh)
    overlap = load_overlap(args.overlap_csv, args.bias_min, args.bias_max)
    sg_edges, states = load_manifest(args.sg_edge_manifest, args.bias_min, args.bias_max)
    sent_cache: dict[str, tuple[dict[int, tuple[float, float]], dict[str, dict[int, float]]]] = {}
    nearest_cache: dict[tuple[str, int], tuple[int, float]] = {}
    rows: list[dict[str, Any]] = []

    for (bkey, edge_id), sg in sorted(sg_edges.items(), key=lambda item: (float(item[0][0]), item[0][1])):
        ov = overlap.get((bkey, edge_id))
        state = states.get(bkey)
        if ov is None or state is None:
            continue
        bias = float(sg["bias_V"])
        if bkey not in sent_cache:
            sent_cache[bkey] = load_sentaurus_case(sentaurus_bias_dir(args.sentaurus_root, bias))
        sent_nodes, sent_fields = sent_cache[bkey]
        node0 = int(sg["node0"])
        node1 = int(sg["node1"])
        if node0 not in mesh_nodes or node1 not in mesh_nodes:
            continue
        s0, d0 = nearest_cache.setdefault((bkey, node0), nearest_node(sent_nodes, *mesh_nodes[node0]))
        s1, d1 = nearest_cache.setdefault((bkey, node1), nearest_node(sent_nodes, *mesh_nodes[node1]))
        length = finite_float(sg.get("edge_length_m"))
        if length <= 0.0:
            x0, y0 = mesh_nodes[node0]
            x1, y1 = mesh_nodes[node1]
            length = math.hypot((x1 - x0) * 1.0e-6, (y1 - y0) * 1.0e-6)
        if length <= 0.0:
            continue

        vela_n_cm3 = (endpoint_avg_state(state, node0, node1, "electrons_m3") or 0.0) * 1.0e-6
        vela_p_cm3 = (endpoint_avg_state(state, node0, node1, "holes_m3") or 0.0) * 1.0e-6
        sent_n_cm3 = avg(sent_fields["eDensity"], s0, s1)
        sent_p_cm3 = avg(sent_fields["hDensity"], s0, s1)
        vela_mu_n_cm2 = finite_float(sg.get("electron_mobility_m2_V_s")) * 1.0e4
        vela_mu_p_cm2 = finite_float(sg.get("hole_mobility_m2_V_s")) * 1.0e4
        sent_mu_n_cm2 = avg(sent_fields["eMobility"], s0, s1)
        sent_mu_p_cm2 = avg(sent_fields["hMobility"], s0, s1)
        vela_e_grad = qf_grad(endpoint_drop_state(state, node0, node1, "phin"), length)
        vela_h_grad = qf_grad(endpoint_drop_state(state, node0, node1, "phip"), length)
        sent_e_grad = qf_grad(drop(sent_fields["eQF"], s0, s1), length)
        sent_h_grad = qf_grad(drop(sent_fields["hQF"], s0, s1), length)

        e_density_ratio = ratio(vela_n_cm3, sent_n_cm3)
        e_mu_ratio = ratio(vela_mu_n_cm2, sent_mu_n_cm2)
        e_grad_ratio = ratio(vela_e_grad, sent_e_grad)
        h_density_ratio = ratio(vela_p_cm3, sent_p_cm3)
        h_mu_ratio = ratio(vela_mu_p_cm2, sent_mu_p_cm2)
        h_grad_ratio = ratio(vela_h_grad, sent_h_grad)
        e_linear_ratio = None if None in (e_density_ratio, e_mu_ratio, e_grad_ratio) else e_density_ratio * e_mu_ratio * e_grad_ratio
        h_linear_ratio = None if None in (h_density_ratio, h_mu_ratio, h_grad_ratio) else h_density_ratio * h_mu_ratio * h_grad_ratio
        sent_e_flux = abs(finite_float(ov.get("sent_e_particle_flux_m2_s")))
        sent_h_flux = abs(finite_float(ov.get("sent_h_particle_flux_m2_s")))
        vela_e_flux = raw_flux(sg, "electron")
        vela_h_flux = raw_flux(sg, "hole")
        rows.append({
            "bias_V": bias,
            "edge_id": edge_id,
            "node0": node0,
            "node1": node1,
            "sentaurus_node0": s0,
            "sentaurus_node1": s1,
            "sentaurus_distance0_um": d0,
            "sentaurus_distance1_um": d1,
            "sent_high_field_p95": str(ov.get("sent_high_field_p95", "")).lower() in {"1", "true", "yes"},
            "edge_area_proxy_m2": finite_float(sg.get("edge_area_proxy_m2"), finite_float(ov.get("edge_area_proxy_m2"))),
            "edge_length_m": length,
            "electron_density_ratio_vela_over_sentaurus": e_density_ratio,
            "electron_mobility_ratio_vela_over_sentaurus": e_mu_ratio,
            "electron_qf_grad_ratio_vela_over_sentaurus": e_grad_ratio,
            "electron_linear_support_ratio_vela_over_sentaurus": e_linear_ratio,
            "hole_density_ratio_vela_over_sentaurus": h_density_ratio,
            "hole_mobility_ratio_vela_over_sentaurus": h_mu_ratio,
            "hole_qf_grad_ratio_vela_over_sentaurus": h_grad_ratio,
            "hole_linear_support_ratio_vela_over_sentaurus": h_linear_ratio,
            "vela_e_density_avg_cm3": vela_n_cm3,
            "sent_e_density_avg_cm3": sent_n_cm3,
            "vela_h_density_avg_cm3": vela_p_cm3,
            "sent_h_density_avg_cm3": sent_p_cm3,
            "vela_e_mobility_cm2_V_s": vela_mu_n_cm2,
            "sent_e_mobility_avg_cm2_V_s": sent_mu_n_cm2,
            "vela_h_mobility_cm2_V_s": vela_mu_p_cm2,
            "sent_h_mobility_avg_cm2_V_s": sent_mu_p_cm2,
            "vela_e_qf_grad_abs_V_m": vela_e_grad,
            "sent_e_qf_grad_abs_V_m": sent_e_grad,
            "vela_h_qf_grad_abs_V_m": vela_h_grad,
            "sent_h_qf_grad_abs_V_m": sent_h_grad,
            "vela_sg_e_flux_m2_s": vela_e_flux,
            "vela_sg_h_flux_m2_s": vela_h_flux,
            "sent_e_flux_m2_s": sent_e_flux,
            "sent_h_flux_m2_s": sent_h_flux,
            "vela_e_flux_over_sentaurus_current": ratio(vela_e_flux, sent_e_flux),
            "vela_h_flux_over_sentaurus_current": ratio(vela_h_flux, sent_h_flux),
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    out: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            out.append(value)
    return out


def median(rows: list[dict[str, Any]], key: str) -> float | None:
    values = finite_values(rows, key)
    return statistics.median(values) if values else None


def weighted_mean(rows: list[dict[str, Any]], key: str, weight_key: str) -> float | None:
    num = 0.0
    den = 0.0
    for row in rows:
        value = optional_float(row.get(key))
        weight = optional_float(row.get(weight_key))
        if value is None or weight is None or weight <= 0.0:
            continue
        num += value * weight
        den += weight
    return ratio(num, den)


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for bias in sorted({float(row["bias_V"]) for row in rows}):
        subset = [row for row in rows if abs(float(row["bias_V"]) - bias) <= 1.0e-9]
        high = [row for row in subset if row["sent_high_field_p95"]]
        result.append({
            "bias_V": bias,
            "edge_count": len(subset),
            "high_field_edge_count": len(high),
            "electron_density_ratio_median": median(subset, "electron_density_ratio_vela_over_sentaurus"),
            "electron_mobility_ratio_median": median(subset, "electron_mobility_ratio_vela_over_sentaurus"),
            "electron_qf_grad_ratio_median": median(subset, "electron_qf_grad_ratio_vela_over_sentaurus"),
            "electron_linear_support_ratio_median": median(subset, "electron_linear_support_ratio_vela_over_sentaurus"),
            "hole_density_ratio_median": median(subset, "hole_density_ratio_vela_over_sentaurus"),
            "hole_mobility_ratio_median": median(subset, "hole_mobility_ratio_vela_over_sentaurus"),
            "hole_qf_grad_ratio_median": median(subset, "hole_qf_grad_ratio_vela_over_sentaurus"),
            "hole_linear_support_ratio_median": median(subset, "hole_linear_support_ratio_vela_over_sentaurus"),
            "electron_density_ratio_high_field_median": median(high, "electron_density_ratio_vela_over_sentaurus"),
            "electron_mobility_ratio_high_field_median": median(high, "electron_mobility_ratio_vela_over_sentaurus"),
            "electron_qf_grad_ratio_high_field_median": median(high, "electron_qf_grad_ratio_vela_over_sentaurus"),
            "electron_linear_support_ratio_high_field_median": median(high, "electron_linear_support_ratio_vela_over_sentaurus"),
            "hole_density_ratio_high_field_median": median(high, "hole_density_ratio_vela_over_sentaurus"),
            "hole_mobility_ratio_high_field_median": median(high, "hole_mobility_ratio_vela_over_sentaurus"),
            "hole_qf_grad_ratio_high_field_median": median(high, "hole_qf_grad_ratio_vela_over_sentaurus"),
            "hole_linear_support_ratio_high_field_median": median(high, "hole_linear_support_ratio_vela_over_sentaurus"),
            "electron_linear_support_ratio_source_weighted": weighted_mean(subset, "electron_linear_support_ratio_vela_over_sentaurus", "sent_e_flux_m2_s"),
            "hole_linear_support_ratio_source_weighted": weighted_mean(subset, "hole_linear_support_ratio_vela_over_sentaurus", "sent_h_flux_m2_s"),
            "electron_flux_over_sentaurus_current_source_weighted": weighted_mean(subset, "vela_e_flux_over_sentaurus_current", "sent_e_flux_m2_s"),
            "hole_flux_over_sentaurus_current_source_weighted": weighted_mean(subset, "vela_h_flux_over_sentaurus_current", "sent_h_flux_m2_s"),
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
    rows = make_rows(args)
    summary = summarize(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "matched_current_factor_edges.csv", EDGE_FIELDS, rows)
    write_csv(args.out_dir / "matched_current_factor_high_field_edges.csv", EDGE_FIELDS, [r for r in rows if r["sent_high_field_p95"]])
    write_csv(args.out_dir / "matched_current_factor_summary.csv", SUMMARY_FIELDS, summary)
    payload = {
        "schema": "pn2d.matched_current_factor_audit.v1",
        "edge_row_count": len(rows),
        "bias_count": len(summary),
        "summary": summary,
    }
    (args.out_dir / "matched_current_factor_summary.json").write_text(
        json.dumps(clean_json(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"edge_row_count": len(rows), "bias_count": len(summary), "out_dir": str(args.out_dir)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
