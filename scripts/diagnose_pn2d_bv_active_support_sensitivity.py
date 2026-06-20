#!/usr/bin/env python3
"""Finite-difference active-support continuity terms with QF perturbations."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from copy import deepcopy
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "carrier",
    "qf_variable",
    "delta_V",
    "baseline_flux_s_inv",
    "baseline_srh_s_inv",
    "baseline_impact_s_inv",
    "baseline_residual_s_inv",
    "forward_flux_s_inv",
    "backward_flux_s_inv",
    "d_flux_dqf_s_inv_per_V",
    "forward_srh_s_inv",
    "backward_srh_s_inv",
    "d_srh_dqf_s_inv_per_V",
    "forward_impact_s_inv",
    "backward_impact_s_inv",
    "d_impact_dqf_s_inv_per_V",
    "forward_residual_s_inv",
    "backward_residual_s_inv",
    "d_residual_dqf_s_inv_per_V",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--nodes", default="", help="Comma-separated node ids; default uses active support nodes.")
    parser.add_argument("--delta-v", type=float, default=1.0e-3)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    parser.add_argument("--taun-s", type=float, default=1.0e-5)
    parser.add_argument("--taup-s", type=float, default=3.0e-6)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def selected_support_nodes(path: Path, node_filter: set[int] | None) -> list[dict[str, str]]:
    rows = []
    for row in read_csv_rows(path):
        if row.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(row["node_id"])
        if node_filter is not None and node_id not in node_filter:
            continue
        rows.append(row)
    return rows


def parse_node_filter(raw: str) -> set[int] | None:
    if not raw.strip():
        return None
    return {int(item.strip()) for item in raw.split(",") if item.strip()}


def limited_exp(value: float) -> float:
    return math.exp(max(-100.0, min(100.0, value)))


def with_reconstructed_carriers(
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
) -> dict[str, list[float]]:
    result = deepcopy(scalars)
    psi = result["Potential"]
    phin = result["ElectronQuasiFermi"]
    phip = result["HoleQuasiFermi"]
    result["Electrons"] = [
        ni[i] * limited_exp((psi[i] - phin[i]) / vt)
        for i in range(len(psi))
    ]
    result["Holes"] = [
        ni[i] * limited_exp((phip[i] - psi[i]) / vt)
        for i in range(len(psi))
    ]
    return result


def srh_rate(n: float, p: float, ni: float, taun: float, taup: float) -> float:
    denominator = taup * (n + ni) + taun * (p + ni)
    if denominator == 0.0:
        return 0.0
    return (n * p - ni * ni) / denominator


def term_state(
    nodes: dict[int, dict[str, float]],
    triangles: list[dict[str, Any]],
    contacts: dict[int, set[str]],
    edges: list[dict[str, Any]],
    volumes: list[float],
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    taun: float,
    taup: float,
) -> dict[str, dict[int, float]]:
    state = with_reconstructed_carriers(scalars, ni, vt)
    source_rows, node_avalanche = cont.source_rows_with_signed_flux(
        nodes,
        edges,
        contacts,
        state,
        ni,
        vt,
    )
    electron_transport, hole_transport = cont.compute_signed_terms(source_rows)
    electron_srh: dict[int, float] = {}
    hole_srh: dict[int, float] = {}
    for node_id in range(len(volumes)):
        rate = srh_rate(
            state["Electrons"][node_id],
            state["Holes"][node_id],
            ni[node_id],
            taun,
            taup,
        )
        electron_srh[node_id] = rate * volumes[node_id]
        hole_srh[node_id] = rate * volumes[node_id]
    return {
        "electron_flux": dict(electron_transport),
        "hole_flux": dict(hole_transport),
        "electron_srh": electron_srh,
        "hole_srh": hole_srh,
        "electron_impact": {idx: node_avalanche[idx] for idx in range(len(node_avalanche))},
        "hole_impact": {idx: node_avalanche[idx] for idx in range(len(node_avalanche))},
    }


def value(terms: dict[str, dict[int, float]], key: str, node_id: int) -> float:
    return float(terms[key].get(node_id, 0.0))


def residual(terms: dict[str, dict[int, float]], carrier: str, node_id: int) -> float:
    return (
        value(terms, f"{carrier}_flux", node_id)
        + value(terms, f"{carrier}_srh", node_id)
        - value(terms, f"{carrier}_impact", node_id)
    )


def perturb_scalars(
    scalars: dict[str, list[float]],
    carrier: str,
    node_id: int,
    delta: float,
) -> dict[str, list[float]]:
    result = deepcopy(scalars)
    key = "ElectronQuasiFermi" if carrier == "electron" else "HoleQuasiFermi"
    result[key][node_id] += delta
    return result


def derivative(forward: float, backward: float, delta: float) -> float:
    return (forward - backward) / (2.0 * delta)


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.delta_v <= 0.0:
        raise SystemExit("--delta-v must be positive")
    nodes, triangles, contacts = cont.sgdiag.read_mesh(args.mesh)
    edges = cont.sgdiag.build_edges(nodes, triangles)
    volumes = cont.sgdiag.node_volumes(nodes, triangles)
    vt = cont.sgdiag.K_B_OVER_Q * args.temperature_k
    vtk = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    scalars = cont.sgdiag.parse_vtk_scalars(vtk)
    ni = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    support = selected_support_nodes(args.support_csv, parse_node_filter(args.nodes))
    baseline = term_state(
        nodes,
        triangles,
        contacts,
        edges,
        volumes,
        scalars,
        ni,
        vt,
        args.taun_s,
        args.taup_s,
    )
    rows: list[dict[str, Any]] = []
    for support_row in support:
        node_id = int(support_row["node_id"])
        for carrier in ["electron", "hole"]:
            forward = term_state(
                nodes,
                triangles,
                contacts,
                edges,
                volumes,
                perturb_scalars(scalars, carrier, node_id, args.delta_v),
                ni,
                vt,
                args.taun_s,
                args.taup_s,
            )
            backward = term_state(
                nodes,
                triangles,
                contacts,
                edges,
                volumes,
                perturb_scalars(scalars, carrier, node_id, -args.delta_v),
                ni,
                vt,
                args.taun_s,
                args.taup_s,
            )
            flux_key = f"{carrier}_flux"
            srh_key = f"{carrier}_srh"
            impact_key = f"{carrier}_impact"
            f_flux = value(forward, flux_key, node_id)
            b_flux = value(backward, flux_key, node_id)
            f_srh = value(forward, srh_key, node_id)
            b_srh = value(backward, srh_key, node_id)
            f_impact = value(forward, impact_key, node_id)
            b_impact = value(backward, impact_key, node_id)
            f_residual = residual(forward, carrier, node_id)
            b_residual = residual(backward, carrier, node_id)
            rows.append({
                "node_id": node_id,
                "x_um": support_row.get("x_um", nodes[node_id]["x_um"]),
                "y_um": support_row.get("y_um", nodes[node_id]["y_um"]),
                "support_class": support_row.get("support_class", ""),
                "carrier": carrier,
                "qf_variable": "phin" if carrier == "electron" else "phip",
                "delta_V": args.delta_v,
                "baseline_flux_s_inv": value(baseline, flux_key, node_id),
                "baseline_srh_s_inv": value(baseline, srh_key, node_id),
                "baseline_impact_s_inv": value(baseline, impact_key, node_id),
                "baseline_residual_s_inv": residual(baseline, carrier, node_id),
                "forward_flux_s_inv": f_flux,
                "backward_flux_s_inv": b_flux,
                "d_flux_dqf_s_inv_per_V": derivative(f_flux, b_flux, args.delta_v),
                "forward_srh_s_inv": f_srh,
                "backward_srh_s_inv": b_srh,
                "d_srh_dqf_s_inv_per_V": derivative(f_srh, b_srh, args.delta_v),
                "forward_impact_s_inv": f_impact,
                "backward_impact_s_inv": b_impact,
                "d_impact_dqf_s_inv_per_V": derivative(f_impact, b_impact, args.delta_v),
                "forward_residual_s_inv": f_residual,
                "backward_residual_s_inv": b_residual,
                "d_residual_dqf_s_inv_per_V": derivative(f_residual, b_residual, args.delta_v),
            })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize(rows: list[dict[str, Any]], args: argparse.Namespace) -> dict[str, Any]:
    by_carrier: dict[str, Any] = {}
    for carrier in sorted({str(row["carrier"]) for row in rows}):
        subset = [row for row in rows if row["carrier"] == carrier]
        item: dict[str, Any] = {"count": len(subset)}
        for key in [
            "d_flux_dqf_s_inv_per_V",
            "d_srh_dqf_s_inv_per_V",
            "d_impact_dqf_s_inv_per_V",
            "d_residual_dqf_s_inv_per_V",
        ]:
            values = finite_values(subset, key)
            abs_values = [abs(value) for value in values]
            if values:
                item[f"{key}_median"] = statistics.median(values)
                item[f"max_abs_{key}"] = max(abs_values)
        by_carrier[carrier] = item
    return {
        "bias": args.bias,
        "delta_V": args.delta_v,
        "row_count": len(rows),
        "by_carrier": by_carrier,
    }


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_support_sensitivity_rows.csv", rows)
    summary = summarize(rows, args)
    (args.out_dir / "active_support_sensitivity_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
