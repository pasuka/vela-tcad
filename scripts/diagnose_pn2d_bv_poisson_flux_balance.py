#!/usr/bin/env python3
"""Decompose PN2D BV Poisson edge-flux and charge balance for external states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


Q = 1.602176634e-19
EPS0 = 8.8541878128e-12

STATE_FIELDS = {
    "psi": ("ElectrostaticPotential", "Potential"),
    "n": ("eDensity", "Electrons"),
    "p": ("hDensity", "Holes"),
}

NODE_FIELDS = [
    "state",
    "bias_V",
    "rank",
    "node_id",
    "x_um",
    "y_um",
    "band",
    "contact_names",
    "volume_m2",
    "psi_V",
    "electron_density_cm3",
    "hole_density_cm3",
    "net_doping_cm3",
    "charge_density_cm3",
    "flux_term_C_per_m",
    "charge_term_C_per_m",
    "residual_C_per_m",
    "abs_residual_C_per_m",
    "residual_per_epsV",
]

BAND_FIELDS = [
    "state",
    "bias_V",
    "band",
    "x_min_um",
    "x_max_um",
    "node_count",
    "median_psi_V",
    "median_charge_density_cm3",
    "median_flux_term_C_per_m",
    "median_charge_term_C_per_m",
    "median_residual_C_per_m",
    "median_abs_residual_C_per_m",
    "max_abs_residual_C_per_m",
    "max_abs_residual_node",
    "median_residual_per_epsV",
    "max_abs_residual_per_epsV",
]

COMPARE_FIELDS = [
    "bias_V",
    "band",
    "x_min_um",
    "x_max_um",
    "node_count",
    "vela_median_residual_C_per_m",
    "sentaurus_median_residual_C_per_m",
    "delta_median_residual_C_per_m",
    "vela_median_abs_residual_C_per_m",
    "sentaurus_median_abs_residual_C_per_m",
    "delta_median_abs_residual_C_per_m",
    "vela_max_abs_residual_C_per_m",
    "sentaurus_max_abs_residual_C_per_m",
    "delta_max_abs_residual_C_per_m",
    "vela_median_residual_per_epsV",
    "sentaurus_median_residual_per_epsV",
    "delta_median_residual_per_epsV",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument("--relative-permittivity", type=float, default=11.7)
    parser.add_argument("--top", type=int, default=120)
    parser.add_argument(
        "--include-contact-nodes",
        action="store_true",
        help="Include Dirichlet contact nodes in rankings and band aggregates.",
    )
    parser.add_argument(
        "--bands",
        default="left_p:0.0:0.25,pre_junction_p:0.7:0.9,junction:0.95:1.05,post_junction_n:1.1:1.3,right_n:1.5:1.8",
        help="Comma-separated name:xmin:xmax bands in um.",
    )
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_bands(raw: str) -> list[tuple[str, float, float]]:
    result = []
    for item in raw.split(","):
        if not item.strip():
            continue
        name, x_min, x_max = item.split(":", 2)
        result.append((name.strip(), float(x_min), float(x_max)))
    return result


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if match:
            result[bias_key(float(match.group("bias")))] = path
    return result


def sentaurus_dir(root: Path, bias: float) -> Path:
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "fields").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def read_mesh(path: Path) -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]], dict[int, set[str]]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {
            "id": int(node["id"]),
            "x_um": float(node["x"]),
            "y_um": float(node["y"]),
        }
        for node in data["nodes"]
    }
    triangles = [
        {
            "id": int(triangle["id"]),
            "node_ids": [int(node_id) for node_id in triangle["node_ids"]],
        }
        for triangle in data["triangles"]
    ]
    contact_by_node: dict[int, set[str]] = defaultdict(set)
    for contact in data.get("contacts", []):
        name = str(contact.get("name", "contact"))
        for node_id in contact.get("node_ids", []):
            contact_by_node[int(node_id)].add(name)
    return nodes, triangles, contact_by_node


def triangle_area_m2(nodes: dict[int, dict[str, Any]], node_ids: list[int]) -> float:
    a, b, c = [nodes[node_id] for node_id in node_ids[:3]]
    ax, ay = a["x_um"] * 1.0e-6, a["y_um"] * 1.0e-6
    bx, by = b["x_um"] * 1.0e-6, b["y_um"] * 1.0e-6
    cx, cy = c["x_um"] * 1.0e-6, c["y_um"] * 1.0e-6
    return 0.5 * abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))


def node_volumes(nodes: dict[int, dict[str, Any]], triangles: list[dict[str, Any]]) -> list[float]:
    volumes = [0.0 for _ in range(max(nodes) + 1)]
    for triangle in triangles:
        area = triangle_area_m2(nodes, triangle["node_ids"])
        for node_id in triangle["node_ids"]:
            volumes[node_id] += area / 3.0
    return volumes


def load_csv_scalar(path: Path, node_count: int) -> list[float]:
    values = [math.nan for _ in range(node_count)]
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            values[int(row["node_id"])] = float(row["component0"])
    missing = [idx for idx, value in enumerate(values) if not math.isfinite(value)]
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return values


def load_doping_csv(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    seen = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            donors[node_id] = float(row["donors_cm3"]) * 1.0e6
            acceptors[node_id] = float(row["acceptors_cm3"]) * 1.0e6
            seen.add(node_id)
    missing = sorted(set(range(node_count)) - seen)
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return donors, acceptors


def load_state(sentaurus: Path, vtk: Path, node_count: int) -> dict[str, dict[str, list[float]]]:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    state: dict[str, dict[str, list[float]]] = {"vela": {}, "sentaurus": {}}
    for short_name, (sentaurus_name, vela_name) in STATE_FIELDS.items():
        if vela_name not in scalars:
            raise RuntimeError(f"{vtk} missing scalar {vela_name}")
        values = scalars[vela_name]
        if len(values) != node_count:
            raise RuntimeError(f"{vela_name} length mismatch in {vtk}")
        state["vela"][short_name] = values
        sentaurus_values = load_csv_scalar(sentaurus / "fields" / f"{sentaurus_name}_region0.csv", node_count)
        if short_name in {"n", "p"}:
            sentaurus_values = [value * 1.0e6 for value in sentaurus_values]
        state["sentaurus"][short_name] = sentaurus_values
    return state


def band_for_x(x_um: float, bands: list[tuple[str, float, float]]) -> str:
    for name, x_min, x_max in bands:
        if x_min <= x_um <= x_max:
            return name
    return "outside"


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def safe_max_abs(rows: list[dict[str, Any]], field: str) -> tuple[float | None, int | None]:
    if not rows:
        return None, None
    row = max(rows, key=lambda item: abs(float(item[field])))
    return abs(float(row[field])), int(row["node_id"])


def compute_flux_balance_rows(
    state_name: str,
    bias: float,
    nodes: dict[int, dict[str, Any]],
    triangles: list[dict[str, Any]],
    contact_by_node: dict[int, set[str]],
    volumes: list[float],
    edges: list[dict[str, Any]],
    psi: list[float],
    n: list[float],
    p: list[float],
    net_doping: list[float],
    eps: float,
    eps_v_scale: float,
    bands: list[tuple[str, float, float]],
) -> list[dict[str, Any]]:
    del triangles
    flux = [0.0 for _ in nodes]
    for edge in edges:
        i = int(edge["node0"])
        j = int(edge["node1"])
        length = float(edge["length_m"])
        if length <= 1.0e-30:
            continue
        conductance = eps * float(edge["couple_m"]) / length
        edge_flux = conductance * (psi[i] - psi[j])
        flux[i] += edge_flux
        flux[j] -= edge_flux

    rows = []
    for node_id in sorted(nodes):
        charge_density = p[node_id] - n[node_id] + net_doping[node_id]
        charge_term = Q * charge_density * volumes[node_id]
        residual = flux[node_id] - charge_term
        node = nodes[node_id]
        rows.append({
            "state": state_name,
            "bias_V": bias,
            "rank": 0,
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "band": band_for_x(float(node["x_um"]), bands),
            "contact_names": ";".join(sorted(contact_by_node.get(node_id, set()))),
            "volume_m2": volumes[node_id],
            "psi_V": psi[node_id],
            "electron_density_cm3": n[node_id] / 1.0e6,
            "hole_density_cm3": p[node_id] / 1.0e6,
            "net_doping_cm3": net_doping[node_id] / 1.0e6,
            "charge_density_cm3": charge_density / 1.0e6,
            "flux_term_C_per_m": flux[node_id],
            "charge_term_C_per_m": charge_term,
            "residual_C_per_m": residual,
            "abs_residual_C_per_m": abs(residual),
            "residual_per_epsV": residual / eps_v_scale if eps_v_scale > 0.0 else math.nan,
        })
    rows.sort(key=lambda row: float(row["abs_residual_C_per_m"]), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def summarize_bands(
    rows: list[dict[str, Any]],
    bands: list[tuple[str, float, float]],
) -> list[dict[str, Any]]:
    result = []
    for name, x_min, x_max in bands:
        selected = [row for row in rows if row["band"] == name]
        max_abs, max_node = safe_max_abs(selected, "residual_C_per_m")
        max_scaled, _ = safe_max_abs(selected, "residual_per_epsV")
        result.append({
            "state": selected[0]["state"] if selected else "",
            "bias_V": selected[0]["bias_V"] if selected else "",
            "band": name,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "node_count": len(selected),
            "median_psi_V": median([float(row["psi_V"]) for row in selected]),
            "median_charge_density_cm3": median([float(row["charge_density_cm3"]) for row in selected]),
            "median_flux_term_C_per_m": median([float(row["flux_term_C_per_m"]) for row in selected]),
            "median_charge_term_C_per_m": median([float(row["charge_term_C_per_m"]) for row in selected]),
            "median_residual_C_per_m": median([float(row["residual_C_per_m"]) for row in selected]),
            "median_abs_residual_C_per_m": median([float(row["abs_residual_C_per_m"]) for row in selected]),
            "max_abs_residual_C_per_m": max_abs,
            "max_abs_residual_node": max_node,
            "median_residual_per_epsV": median([float(row["residual_per_epsV"]) for row in selected]),
            "max_abs_residual_per_epsV": max_scaled,
        })
    return result


def compare_bands(vela_rows: list[dict[str, Any]],
                  sentaurus_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    sentaurus_by_key = {
        (float(row["bias_V"]), row["band"]): row
        for row in sentaurus_rows
    }
    result = []
    for vela in vela_rows:
        key = (float(vela["bias_V"]), vela["band"])
        sentaurus = sentaurus_by_key.get(key)
        if sentaurus is None:
            continue
        row = {
            "bias_V": vela["bias_V"],
            "band": vela["band"],
            "x_min_um": vela["x_min_um"],
            "x_max_um": vela["x_max_um"],
            "node_count": vela["node_count"],
        }
        pairs = [
            ("median_residual_C_per_m", "median_residual_C_per_m"),
            ("median_abs_residual_C_per_m", "median_abs_residual_C_per_m"),
            ("max_abs_residual_C_per_m", "max_abs_residual_C_per_m"),
            ("median_residual_per_epsV", "median_residual_per_epsV"),
        ]
        for source_field, output_field in pairs:
            v_value = vela[source_field]
            s_value = sentaurus[source_field]
            short = output_field
            row[f"vela_{short}"] = v_value
            row[f"sentaurus_{short}"] = s_value
            row[f"delta_{short}"] = (
                v_value - s_value
                if v_value is not None and s_value is not None
                else None
            )
        result.append(row)
    return result


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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
    nodes, triangles, contact_by_node = read_mesh(args.mesh)
    node_count = len(nodes)
    bands = parse_bands(args.bands)
    vtks = discover_vela_vtks(args.vela_vtk_root)
    volumes = node_volumes(nodes, triangles)
    edges = sgdiag.build_edges(nodes, triangles)
    vela_donors, vela_acceptors = load_doping_csv(args.doping_csv, node_count)
    vela_net = [donor - acceptor for donor, acceptor in zip(vela_donors, vela_acceptors)]
    eps = EPS0 * args.relative_permittivity
    eps_v_scale = eps

    all_top_rows: list[dict[str, Any]] = []
    vela_band_rows: list[dict[str, Any]] = []
    sentaurus_band_rows: list[dict[str, Any]] = []
    summaries = []
    for bias in parse_biases(args.biases):
        vtk = vtks.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        sentaurus = sentaurus_dir(args.sentaurus_root, bias)
        state = load_state(sentaurus, vtk, node_count)
        sentaurus_donors = load_csv_scalar(sentaurus / "fields" / "DonorConcentration_region0.csv", node_count)
        sentaurus_acceptors = load_csv_scalar(sentaurus / "fields" / "AcceptorConcentration_region0.csv", node_count)
        sentaurus_net = [
            (donor - acceptor) * 1.0e6
            for donor, acceptor in zip(sentaurus_donors, sentaurus_acceptors)
        ]
        bias_summary: dict[str, Any] = {
            "bias_V": bias,
            "vtk": str(vtk),
            "sentaurus_dir": str(sentaurus),
            "states": {},
        }
        for source, net in (("vela", vela_net), ("sentaurus", sentaurus_net)):
            rows = compute_flux_balance_rows(
                source,
                bias,
                nodes,
                triangles,
                contact_by_node,
                volumes,
                edges,
                state[source]["psi"],
                state[source]["n"],
                state[source]["p"],
                net,
                eps,
                eps_v_scale,
                bands,
            )
            diagnostic_rows = rows if args.include_contact_nodes else [
                row for row in rows if not row["contact_names"]
            ]
            band_rows = summarize_bands(diagnostic_rows, bands)
            all_top_rows.extend(diagnostic_rows[:args.top])
            if source == "vela":
                vela_band_rows.extend(band_rows)
            else:
                sentaurus_band_rows.extend(band_rows)
            bias_summary["states"][source] = {
                "top_abs_residual": diagnostic_rows[0] if diagnostic_rows else None,
                "contact_nodes_included": args.include_contact_nodes,
                "bands": band_rows,
            }
        summaries.append(bias_summary)

    comparison = compare_bands(vela_band_rows, sentaurus_band_rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "poisson_flux_balance_top_nodes.csv", all_top_rows, NODE_FIELDS)
    write_csv(args.out_dir / "poisson_flux_balance_bands.csv", vela_band_rows + sentaurus_band_rows, BAND_FIELDS)
    write_csv(args.out_dir / "poisson_flux_balance_band_compare.csv", comparison, COMPARE_FIELDS)
    (args.out_dir / "poisson_flux_balance_summary.json").write_text(
        json.dumps(clean_json({
            "biases": summaries,
            "top_rows": len(all_top_rows),
            "band_rows": len(vela_band_rows) + len(sentaurus_band_rows),
            "compare_rows": len(comparison),
            "relative_permittivity": args.relative_permittivity,
            "bands": [
                {"name": name, "x_min_um": x_min, "x_max_um": x_max}
                for name, x_min, x_max in bands
            ],
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
