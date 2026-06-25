#!/usr/bin/env python3
"""Decompose matched PN2D carrier-density deficit into level and ni-eff terms."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import pn2d_fixed_state_manifest as fixed_manifest


K_B_OVER_Q = 8.617333262145179e-5

EDGE_FIELDS = [
    "bias_V", "edge_id", "node0", "node1", "sentaurus_node0", "sentaurus_node1",
    "sentaurus_distance0_um", "sentaurus_distance1_um", "sent_high_field_p95",
    "electron_density_ratio_vela_over_sentaurus", "electron_level_delta_V",
    "electron_level_density_ratio", "electron_inferred_ni_eff_ratio",
    "hole_density_ratio_vela_over_sentaurus", "hole_level_delta_V",
    "hole_level_density_ratio", "hole_inferred_ni_eff_ratio",
    "vela_psi_minus_phin_avg_V", "sentaurus_psi_minus_phin_avg_V",
    "vela_phip_minus_psi_avg_V", "sentaurus_phip_minus_psi_avg_V",
    "vela_e_density_avg_cm3", "sent_e_density_avg_cm3",
    "vela_h_density_avg_cm3", "sent_h_density_avg_cm3",
]

SUMMARY_FIELDS = [
    "bias_V", "edge_count", "high_field_edge_count",
    "electron_density_ratio_median", "electron_level_density_ratio_median", "electron_inferred_ni_eff_ratio_median",
    "hole_density_ratio_median", "hole_level_density_ratio_median", "hole_inferred_ni_eff_ratio_median",
    "electron_density_ratio_high_field_median", "electron_level_density_ratio_high_field_median", "electron_inferred_ni_eff_ratio_high_field_median",
    "hole_density_ratio_high_field_median", "hole_level_density_ratio_high_field_median", "hole_inferred_ni_eff_ratio_high_field_median",
    "electron_level_delta_high_field_median_mV", "hole_level_delta_high_field_median_mV",
]


def parse_args() -> argparse.Namespace:
    root = Path("build-release/reference_tcad/pn2d_sentaurus2018")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlap-csv", type=Path, default=root / "reports/high_field_overlap_compare/high_field_overlap_edges.csv")
    parser.add_argument("--sg-edge-manifest", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, default=root / "sentaurus_multibias_0p25v")
    parser.add_argument("--mesh", type=Path, default=Path("build/pn2d_bv_gradqf_full20/mesh.json"))
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-16.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
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


def ratio(a: float | None, b: float | None) -> float | None:
    if a is None or b is None or b == 0.0:
        return None
    return a / b


def safe_exp(x: float | None) -> float | None:
    if x is None:
        return None
    if x > 700.0:
        return math.inf
    if x < -745.0:
        return 0.0
    return math.exp(x)


def in_bias_window(bias: float, lo: float, hi: float) -> bool:
    return min(lo, hi) - 1.0e-9 <= bias <= max(lo, hi) + 1.0e-9


def bias_key(bias: float) -> str:
    return f"{bias:.9f}"


def sentaurus_bias_dir(root: Path, bias: float) -> Path:
    return root / f"sentaurus_-{str(abs(bias)).rstrip('0').rstrip('.').replace('.', 'p')}v"


def load_mesh_nodes(path: Path) -> dict[int, tuple[float, float]]:
    data = read_json(path)
    return {int(node["id"]): (float(node["x"]), float(node["y"])) for node in data.get("nodes", [])}


def load_sentaurus_nodes(case_dir: Path) -> dict[int, tuple[float, float]]:
    return {int(row["id"]): (float(row["x_um"]), float(row["y_um"])) for row in read_csv(case_dir / "nodes.csv")}


def nearest_node(nodes: dict[int, tuple[float, float]], x_um: float, y_um: float) -> tuple[int, float]:
    best = (-1, math.inf)
    for node_id, (x, y) in nodes.items():
        d2 = (x - x_um) ** 2 + (y - y_um) ** 2
        if d2 < best[1]:
            best = (node_id, d2)
    return best[0], math.sqrt(best[1])


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


def load_sentaurus_case(case_dir: Path) -> tuple[dict[int, tuple[float, float]], dict[str, dict[int, float]]]:
    fields = case_dir / "fields"
    return load_sentaurus_nodes(case_dir), {
        "psi": load_scalar_field(fields, "ElectrostaticPotential"),
        "phin": load_scalar_field(fields, "eQuasiFermiPotential"),
        "phip": load_scalar_field(fields, "hQuasiFermiPotential"),
        "n": load_scalar_field(fields, "eDensity"),
        "p": load_scalar_field(fields, "hDensity"),
    }


def avg_field(field: dict[int, float], a: int, b: int) -> float | None:
    if a not in field or b not in field:
        return None
    return 0.5 * (field[a] + field[b])


def load_state(path: Path) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    for row in read_csv(path):
        result[int(row["node_id"])] = {
            "psi": finite_float(row.get("psi")),
            "phin": finite_float(row.get("phin")),
            "phip": finite_float(row.get("phip")),
            "n_cm3": finite_float(row.get("electrons_m3")) * 1.0e-6,
            "p_cm3": finite_float(row.get("holes_m3")) * 1.0e-6,
        }
    return result


def state_path_from_manifest_item(item: dict[str, Any]) -> Path | None:
    if item.get("state_csv"):
        return Path(item["state_csv"])
    if not item.get("config"):
        return None
    config = read_json(Path(item["config"]))
    raw = (config.get("sweep") or {}).get("write_state_file")
    return Path(raw) if raw else None


def load_manifest(path: Path, lo: float, hi: float) -> tuple[dict[tuple[str, int], dict[str, str]], dict[str, dict[int, dict[str, float]]]]:
    edges: dict[tuple[str, int], dict[str, str]] = {}
    states: dict[str, dict[int, dict[str, float]]] = {}
    for item in read_json(path):
        bias = float(item["bias_V"])
        if not in_bias_window(bias, lo, hi):
            continue
        bkey = bias_key(bias)
        fixed_manifest.require_materials_file(item, path)
        for row in read_csv(Path(item["sg_edges_csv"])):
            edges[(bias_key(float(row["bias_V"])), int(row["edge_id"]))] = row
        state_path = state_path_from_manifest_item(item)
        if state_path and state_path.exists():
            states[bkey] = load_state(state_path)
    return edges, states


def load_overlap(path: Path, lo: float, hi: float) -> dict[tuple[str, int], dict[str, str]]:
    result: dict[tuple[str, int], dict[str, str]] = {}
    for row in read_csv(path):
        bias = float(row["bias_V"])
        if in_bias_window(bias, lo, hi):
            result[(bias_key(bias), int(row["edge_id"]))] = row
    return result


def avg_state(state: dict[int, dict[str, float]], a: int, b: int, key: str) -> float | None:
    if a not in state or b not in state:
        return None
    return 0.5 * (state[a][key] + state[b][key])


def bool_cell(raw: Any) -> bool:
    return str(raw).strip().lower() in {"1", "true", "yes", "y"}


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    vt = K_B_OVER_Q * args.temperature_k
    mesh_nodes = load_mesh_nodes(args.mesh)
    overlap = load_overlap(args.overlap_csv, args.bias_min, args.bias_max)
    edges, states = load_manifest(args.sg_edge_manifest, args.bias_min, args.bias_max)
    sent_cache: dict[str, tuple[dict[int, tuple[float, float]], dict[str, dict[int, float]]]] = {}
    nearest: dict[tuple[str, int], tuple[int, float]] = {}
    rows: list[dict[str, Any]] = []

    for (bkey, edge_id), edge in sorted(edges.items(), key=lambda item: (float(item[0][0]), item[0][1])):
        ov = overlap.get((bkey, edge_id))
        state = states.get(bkey)
        if ov is None or state is None:
            continue
        bias = float(edge["bias_V"])
        if bkey not in sent_cache:
            sent_cache[bkey] = load_sentaurus_case(sentaurus_bias_dir(args.sentaurus_root, bias))
        sent_nodes, sent = sent_cache[bkey]
        node0 = int(edge["node0"])
        node1 = int(edge["node1"])
        if node0 not in mesh_nodes or node1 not in mesh_nodes:
            continue
        nkey0 = (bkey, node0)
        nkey1 = (bkey, node1)
        if nkey0 not in nearest:
            nearest[nkey0] = nearest_node(sent_nodes, *mesh_nodes[node0])
        if nkey1 not in nearest:
            nearest[nkey1] = nearest_node(sent_nodes, *mesh_nodes[node1])
        s0, d0 = nearest[nkey0]
        s1, d1 = nearest[nkey1]

        vela_psi = avg_state(state, node0, node1, "psi")
        vela_phin = avg_state(state, node0, node1, "phin")
        vela_phip = avg_state(state, node0, node1, "phip")
        vela_n = avg_state(state, node0, node1, "n_cm3")
        vela_p = avg_state(state, node0, node1, "p_cm3")
        sent_psi = avg_field(sent["psi"], s0, s1)
        sent_phin = avg_field(sent["phin"], s0, s1)
        sent_phip = avg_field(sent["phip"], s0, s1)
        sent_n = avg_field(sent["n"], s0, s1)
        sent_p = avg_field(sent["p"], s0, s1)
        if None in (vela_psi, vela_phin, vela_phip, vela_n, vela_p, sent_psi, sent_phin, sent_phip, sent_n, sent_p):
            continue

        vela_e_level = vela_psi - vela_phin
        sent_e_level = sent_psi - sent_phin
        vela_h_level = vela_phip - vela_psi
        sent_h_level = sent_phip - sent_psi
        e_delta = vela_e_level - sent_e_level
        h_delta = vela_h_level - sent_h_level
        e_density_ratio = ratio(vela_n, sent_n)
        h_density_ratio = ratio(vela_p, sent_p)
        e_level_ratio = safe_exp(e_delta / vt)
        h_level_ratio = safe_exp(h_delta / vt)
        rows.append({
            "bias_V": bias,
            "edge_id": edge_id,
            "node0": node0,
            "node1": node1,
            "sentaurus_node0": s0,
            "sentaurus_node1": s1,
            "sentaurus_distance0_um": d0,
            "sentaurus_distance1_um": d1,
            "sent_high_field_p95": bool_cell(ov.get("sent_high_field_p95")),
            "electron_density_ratio_vela_over_sentaurus": e_density_ratio,
            "electron_level_delta_V": e_delta,
            "electron_level_density_ratio": e_level_ratio,
            "electron_inferred_ni_eff_ratio": ratio(e_density_ratio, e_level_ratio),
            "hole_density_ratio_vela_over_sentaurus": h_density_ratio,
            "hole_level_delta_V": h_delta,
            "hole_level_density_ratio": h_level_ratio,
            "hole_inferred_ni_eff_ratio": ratio(h_density_ratio, h_level_ratio),
            "vela_psi_minus_phin_avg_V": vela_e_level,
            "sentaurus_psi_minus_phin_avg_V": sent_e_level,
            "vela_phip_minus_psi_avg_V": vela_h_level,
            "sentaurus_phip_minus_psi_avg_V": sent_h_level,
            "vela_e_density_avg_cm3": vela_n,
            "sent_e_density_avg_cm3": sent_n,
            "vela_h_density_avg_cm3": vela_p,
            "sent_h_density_avg_cm3": sent_p,
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def median(rows: list[dict[str, Any]], key: str) -> float | None:
    values = finite_values(rows, key)
    return statistics.median(values) if values else None


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for bias in sorted({float(row["bias_V"]) for row in rows}):
        subset = [row for row in rows if abs(float(row["bias_V"]) - bias) <= 1.0e-9]
        high = [row for row in subset if row["sent_high_field_p95"]]
        out.append({
            "bias_V": bias,
            "edge_count": len(subset),
            "high_field_edge_count": len(high),
            "electron_density_ratio_median": median(subset, "electron_density_ratio_vela_over_sentaurus"),
            "electron_level_density_ratio_median": median(subset, "electron_level_density_ratio"),
            "electron_inferred_ni_eff_ratio_median": median(subset, "electron_inferred_ni_eff_ratio"),
            "hole_density_ratio_median": median(subset, "hole_density_ratio_vela_over_sentaurus"),
            "hole_level_density_ratio_median": median(subset, "hole_level_density_ratio"),
            "hole_inferred_ni_eff_ratio_median": median(subset, "hole_inferred_ni_eff_ratio"),
            "electron_density_ratio_high_field_median": median(high, "electron_density_ratio_vela_over_sentaurus"),
            "electron_level_density_ratio_high_field_median": median(high, "electron_level_density_ratio"),
            "electron_inferred_ni_eff_ratio_high_field_median": median(high, "electron_inferred_ni_eff_ratio"),
            "hole_density_ratio_high_field_median": median(high, "hole_density_ratio_vela_over_sentaurus"),
            "hole_level_density_ratio_high_field_median": median(high, "hole_level_density_ratio"),
            "hole_inferred_ni_eff_ratio_high_field_median": median(high, "hole_inferred_ni_eff_ratio"),
            "electron_level_delta_high_field_median_mV": (median(high, "electron_level_delta_V") or 0.0) * 1.0e3,
            "hole_level_delta_high_field_median_mV": (median(high, "hole_level_delta_V") or 0.0) * 1.0e3,
        })
    return out


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
    write_csv(args.out_dir / "density_deficit_origin_edges.csv", EDGE_FIELDS, rows)
    write_csv(args.out_dir / "density_deficit_origin_high_field_edges.csv", EDGE_FIELDS, [row for row in rows if row["sent_high_field_p95"]])
    write_csv(args.out_dir / "density_deficit_origin_summary.csv", SUMMARY_FIELDS, summary)
    payload = {"schema": "pn2d.density_deficit_origin.v1", "edge_row_count": len(rows), "bias_count": len(summary), "summary": summary}
    (args.out_dir / "density_deficit_origin_summary.json").write_text(json.dumps(clean_json(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"bias_count": len(summary), "edge_row_count": len(rows), "out_dir": str(args.out_dir)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
