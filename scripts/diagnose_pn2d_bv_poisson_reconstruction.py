#!/usr/bin/env python3
"""Solve Vela-style Poisson reconstructions with frozen PN2D BV charge states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any

import numpy as np

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


Q = 1.602176634e-19
EPS0 = 8.8541878128e-12
K_B_OVER_Q = 8.617333262145e-5

FIELD_MAP = {
    "psi": ("ElectrostaticPotential", "Potential"),
    "n": ("eDensity", "Electrons"),
    "p": ("hDensity", "Holes"),
}

BAND_FIELDS = [
    "bias_V",
    "bc_source",
    "charge_source",
    "band",
    "x_min_um",
    "x_max_um",
    "node_count",
    "median_charge_density_cm3",
    "median_reconstructed_psi_V",
    "median_vela_psi_V",
    "median_sentaurus_psi_V",
    "median_reconstructed_minus_vela_V",
    "median_reconstructed_minus_sentaurus_V",
    "median_vela_minus_sentaurus_V",
    "max_abs_reconstructed_minus_vela_V",
    "max_abs_reconstructed_minus_sentaurus_V",
]

NODE_FIELDS = [
    "bias_V",
    "bc_source",
    "charge_source",
    "node_id",
    "x_um",
    "y_um",
    "contact_names",
    "charge_density_cm3",
    "reconstructed_psi_V",
    "vela_psi_V",
    "sentaurus_psi_V",
    "reconstructed_minus_vela_V",
    "reconstructed_minus_sentaurus_V",
    "vela_minus_sentaurus_V",
]

CONTACT_FIELDS = [
    "bias_V",
    "bc_source",
    "contact",
    "node_count",
    "applied_bias_V",
    "net_doping_mean_cm3",
    "builtin_psi_V",
    "dirichlet_psi_V",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument("--focus-nodes", default="955,1089,351,986")
    parser.add_argument("--relative-permittivity", type=float, default=11.7)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.0e16)
    parser.add_argument("--bias-contact", default="Anode")
    parser.add_argument("--bc-sources", default="vela_expected,sentaurus_state")
    parser.add_argument("--charge-sources", default="depletion,vela_frozen,sentaurus_frozen")
    parser.add_argument(
        "--bands",
        default="left_p:0.0:0.25,pre_junction_p:0.7:0.9,junction:0.95:1.05,post_junction_n:1.1:1.3,right_n:1.5:1.8",
        help="Comma-separated name:xmin:xmax bands in um.",
    )
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def parse_token_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


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


def read_mesh(path: Path) -> tuple[dict[int, dict[str, float]], list[dict[str, Any]], list[dict[str, Any]]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {
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
    contacts = [
        {
            "name": str(contact.get("name", f"contact_{idx}")),
            "node_ids": [int(node_id) for node_id in contact.get("node_ids", [])],
        }
        for idx, contact in enumerate(data.get("contacts", []))
    ]
    return nodes, triangles, contacts


def triangle_area_m2(nodes: dict[int, dict[str, float]], node_ids: list[int]) -> float:
    a, b, c = [nodes[node_id] for node_id in node_ids[:3]]
    ax, ay = a["x_um"] * 1.0e-6, a["y_um"] * 1.0e-6
    bx, by = b["x_um"] * 1.0e-6, b["y_um"] * 1.0e-6
    cx, cy = c["x_um"] * 1.0e-6, c["y_um"] * 1.0e-6
    return 0.5 * abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))


def node_volumes(nodes: dict[int, dict[str, float]], triangles: list[dict[str, Any]]) -> list[float]:
    volumes = [0.0 for _ in range(max(nodes) + 1)]
    for triangle in triangles:
        area = triangle_area_m2(nodes, triangle["node_ids"])
        for node_id in triangle["node_ids"]:
            volumes[node_id] += area / 3.0
    return volumes


def contact_by_node(contacts: list[dict[str, Any]]) -> dict[int, set[str]]:
    result: dict[int, set[str]] = {}
    for contact in contacts:
        for node_id in contact["node_ids"]:
            result.setdefault(node_id, set()).add(contact["name"])
    return result


def load_doping_csv(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            donors[node_id] = float(row["donors_cm3"]) * 1.0e6
            acceptors[node_id] = float(row["acceptors_cm3"]) * 1.0e6
    return donors, acceptors


def load_csv_scalar(path: Path, node_count: int) -> list[float]:
    values = [math.nan for _ in range(node_count)]
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            values[int(row["node_id"])] = float(row["component0"])
    missing = [idx for idx, value in enumerate(values) if not math.isfinite(value)]
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return values


def load_state(sentaurus: Path, vtk: Path, node_count: int) -> dict[str, dict[str, list[float]]]:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    state: dict[str, dict[str, list[float]]] = {"vela": {}, "sentaurus": {}}
    for short_name, (sentaurus_name, vela_name) in FIELD_MAP.items():
        if vela_name not in scalars:
            raise RuntimeError(f"{vtk} missing scalar {vela_name}")
        values = list(scalars[vela_name])
        if len(values) != node_count:
            raise RuntimeError(f"{vela_name} has {len(values)} values, expected {node_count}")
        state["vela"][short_name] = values
        sentaurus_values = load_csv_scalar(sentaurus / "fields" / f"{sentaurus_name}_region0.csv", node_count)
        if short_name in {"n", "p"}:
            sentaurus_values = [value * 1.0e6 for value in sentaurus_values]
        state["sentaurus"][short_name] = sentaurus_values
    sentaurus_donors = load_csv_scalar(sentaurus / "fields" / "DonorConcentration_region0.csv", node_count)
    sentaurus_acceptors = load_csv_scalar(sentaurus / "fields" / "AcceptorConcentration_region0.csv", node_count)
    state["sentaurus"]["net_doping"] = [
        (donor - acceptor) * 1.0e6
        for donor, acceptor in zip(sentaurus_donors, sentaurus_acceptors)
    ]
    return state


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def max_abs(values: list[float]) -> float | None:
    return max((abs(value) for value in values), default=None)


def equilibrium_electron_density(net_m3: float, ni_m3: float) -> float:
    return 0.5 * (net_m3 + math.sqrt(net_m3 * net_m3 + 4.0 * ni_m3 * ni_m3))


def builtin_potential(net_cm3: float, temperature_k: float, ni_m3: float) -> float:
    net_m3 = net_cm3 * 1.0e6
    neq = equilibrium_electron_density(net_m3, ni_m3)
    vt = K_B_OVER_Q * temperature_k
    return vt * math.log(neq / ni_m3)


def median_at(values: list[float], node_ids: list[int]) -> float:
    return statistics.median(values[node_id] for node_id in node_ids)


def contact_dirichlet_values(
    bc_source: str,
    bias: float,
    contacts: list[dict[str, Any]],
    net_doping_m3: list[float],
    state: dict[str, dict[str, list[float]]],
    bias_contact: str,
    temperature_k: float,
    ni_m3: float,
) -> tuple[dict[int, float], list[dict[str, Any]]]:
    values: dict[int, float] = {}
    rows = []
    for contact in contacts:
        node_ids = contact["node_ids"]
        if not node_ids:
            continue
        name = contact["name"]
        applied_bias = bias if name == bias_contact else 0.0
        net_mean_cm3 = statistics.fmean(net_doping_m3[node_id] for node_id in node_ids) / 1.0e6
        builtin = builtin_potential(net_mean_cm3, temperature_k, ni_m3)
        if bc_source == "vela_expected":
            psi = applied_bias + builtin
        elif bc_source == "vela_state":
            psi = median_at(state["vela"]["psi"], node_ids)
        elif bc_source == "sentaurus_state":
            psi = median_at(state["sentaurus"]["psi"], node_ids)
        else:
            raise ValueError(f"unknown bc source: {bc_source}")
        for node_id in node_ids:
            values[node_id] = psi
        rows.append({
            "bias_V": bias,
            "bc_source": bc_source,
            "contact": name,
            "node_count": len(node_ids),
            "applied_bias_V": applied_bias,
            "net_doping_mean_cm3": net_mean_cm3,
            "builtin_psi_V": builtin,
            "dirichlet_psi_V": psi,
        })
    return values, rows


def charge_density_m3(charge_source: str,
                      net_doping_m3: list[float],
                      state: dict[str, dict[str, list[float]]]) -> list[float]:
    if charge_source == "depletion":
        return list(net_doping_m3)
    if charge_source == "vela_frozen":
        return [
            p - n + net
            for p, n, net in zip(state["vela"]["p"], state["vela"]["n"], net_doping_m3)
        ]
    if charge_source == "sentaurus_frozen":
        sentaurus_net = state["sentaurus"].get("net_doping", net_doping_m3)
        return [
            p - n + net
            for p, n, net in zip(state["sentaurus"]["p"], state["sentaurus"]["n"], sentaurus_net)
        ]
    raise ValueError(f"unknown charge source: {charge_source}")


def assemble_poisson_matrix(
    node_count: int,
    edges: list[dict[str, Any]],
    eps: float,
) -> np.ndarray:
    matrix = np.zeros((node_count, node_count), dtype=float)
    for edge in edges:
        i = int(edge["node0"])
        j = int(edge["node1"])
        length = float(edge["length_m"])
        if length <= 1.0e-30:
            continue
        conductance = eps * float(edge["couple_m"]) / length
        if conductance == 0.0:
            continue
        matrix[i, i] += conductance
        matrix[i, j] -= conductance
        matrix[j, j] += conductance
        matrix[j, i] -= conductance
    return matrix


def solve_poisson(
    base_matrix: np.ndarray,
    charge_m3: list[float],
    volumes_m2: list[float],
    dirichlet: dict[int, float],
) -> np.ndarray:
    matrix = np.array(base_matrix, copy=True)
    rhs = np.array([Q * rho * vol for rho, vol in zip(charge_m3, volumes_m2)], dtype=float)
    for node_id, value in dirichlet.items():
        matrix[node_id, :] = 0.0
        matrix[node_id, node_id] = 1.0
        rhs[node_id] = value
    return np.linalg.solve(matrix, rhs)


def band_rows(
    bias: float,
    bc_source: str,
    charge_source: str,
    nodes: dict[int, dict[str, float]],
    charge: list[float],
    reconstructed: np.ndarray,
    state: dict[str, dict[str, list[float]]],
    bands: list[tuple[str, float, float]],
) -> list[dict[str, Any]]:
    rows = []
    for name, x_min, x_max in bands:
        ids = [
            node_id for node_id, node in nodes.items()
            if x_min <= node["x_um"] <= x_max
        ]
        recon = [float(reconstructed[node_id]) for node_id in ids]
        vela = [state["vela"]["psi"][node_id] for node_id in ids]
        sentaurus = [state["sentaurus"]["psi"][node_id] for node_id in ids]
        diff_v = [r - v for r, v in zip(recon, vela)]
        diff_s = [r - s for r, s in zip(recon, sentaurus)]
        diff_vs = [v - s for v, s in zip(vela, sentaurus)]
        rows.append({
            "bias_V": bias,
            "bc_source": bc_source,
            "charge_source": charge_source,
            "band": name,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "node_count": len(ids),
            "median_charge_density_cm3": median([charge[node_id] / 1.0e6 for node_id in ids]),
            "median_reconstructed_psi_V": median(recon),
            "median_vela_psi_V": median(vela),
            "median_sentaurus_psi_V": median(sentaurus),
            "median_reconstructed_minus_vela_V": median(diff_v),
            "median_reconstructed_minus_sentaurus_V": median(diff_s),
            "median_vela_minus_sentaurus_V": median(diff_vs),
            "max_abs_reconstructed_minus_vela_V": max_abs(diff_v),
            "max_abs_reconstructed_minus_sentaurus_V": max_abs(diff_s),
        })
    return rows


def node_rows(
    bias: float,
    bc_source: str,
    charge_source: str,
    nodes: dict[int, dict[str, float]],
    focus_nodes: list[int],
    contact_names: dict[int, set[str]],
    charge: list[float],
    reconstructed: np.ndarray,
    state: dict[str, dict[str, list[float]]],
) -> list[dict[str, Any]]:
    rows = []
    for node_id in focus_nodes:
        if node_id not in nodes:
            continue
        recon = float(reconstructed[node_id])
        vela = state["vela"]["psi"][node_id]
        sentaurus = state["sentaurus"]["psi"][node_id]
        rows.append({
            "bias_V": bias,
            "bc_source": bc_source,
            "charge_source": charge_source,
            "node_id": node_id,
            "x_um": nodes[node_id]["x_um"],
            "y_um": nodes[node_id]["y_um"],
            "contact_names": ";".join(sorted(contact_names.get(node_id, set()))),
            "charge_density_cm3": charge[node_id] / 1.0e6,
            "reconstructed_psi_V": recon,
            "vela_psi_V": vela,
            "sentaurus_psi_V": sentaurus,
            "reconstructed_minus_vela_V": recon - vela,
            "reconstructed_minus_sentaurus_V": recon - sentaurus,
            "vela_minus_sentaurus_V": vela - sentaurus,
        })
    return rows


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
    if isinstance(value, np.floating):
        return float(value)
    return value


def main() -> int:
    args = parse_args()
    nodes, triangles, contacts = read_mesh(args.mesh)
    node_count = len(nodes)
    volumes = node_volumes(nodes, triangles)
    edges = sgdiag.build_edges(nodes, triangles)
    net_doping_m3 = [
        donor - acceptor
        for donor, acceptor in zip(*load_doping_csv(args.doping_csv, node_count))
    ]
    eps = EPS0 * args.relative_permittivity
    base_matrix = assemble_poisson_matrix(node_count, edges, eps)
    vtks = discover_vela_vtks(args.vela_vtk_root)
    bands = parse_bands(args.bands)
    focus_nodes = parse_int_list(args.focus_nodes)
    bc_sources = parse_token_list(args.bc_sources)
    charge_sources = parse_token_list(args.charge_sources)
    contacts_by_node = contact_by_node(contacts)

    all_band_rows: list[dict[str, Any]] = []
    all_node_rows: list[dict[str, Any]] = []
    all_contact_rows: list[dict[str, Any]] = []
    summaries = []
    for bias in parse_float_list(args.biases):
        vtk = vtks.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        sentaurus = sentaurus_dir(args.sentaurus_root, bias)
        state = load_state(sentaurus, vtk, node_count)
        bias_summary: dict[str, Any] = {
            "bias_V": bias,
            "vtk": str(vtk),
            "sentaurus_dir": str(sentaurus),
            "solutions": [],
        }
        for bc_source in bc_sources:
            dirichlet, contact_rows = contact_dirichlet_values(
                bc_source,
                bias,
                contacts,
                net_doping_m3,
                state,
                args.bias_contact,
                args.temperature_k,
                args.material_ni_m3,
            )
            all_contact_rows.extend(contact_rows)
            for charge_source in charge_sources:
                charge = charge_density_m3(charge_source, net_doping_m3, state)
                reconstructed = solve_poisson(base_matrix, charge, volumes, dirichlet)
                rows_b = band_rows(
                    bias, bc_source, charge_source, nodes, charge, reconstructed, state, bands)
                rows_n = node_rows(
                    bias, bc_source, charge_source, nodes, focus_nodes,
                    contacts_by_node, charge, reconstructed, state)
                all_band_rows.extend(rows_b)
                all_node_rows.extend(rows_n)
                bias_summary["solutions"].append({
                    "bc_source": bc_source,
                    "charge_source": charge_source,
                    "bands": rows_b,
                    "focus_nodes": rows_n,
                })
        summaries.append(bias_summary)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "poisson_reconstruction_bands.csv", all_band_rows, BAND_FIELDS)
    write_csv(args.out_dir / "poisson_reconstruction_focus_nodes.csv", all_node_rows, NODE_FIELDS)
    write_csv(args.out_dir / "poisson_reconstruction_contacts.csv", all_contact_rows, CONTACT_FIELDS)
    (args.out_dir / "poisson_reconstruction_summary.json").write_text(
        json.dumps(clean_json({
            "mesh": str(args.mesh),
            "biases": summaries,
            "band_rows": len(all_band_rows),
            "focus_node_rows": len(all_node_rows),
            "contact_rows": len(all_contact_rows),
            "bc_sources": bc_sources,
            "charge_sources": charge_sources,
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
