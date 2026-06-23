#!/usr/bin/env python3
"""Compare active-support continuity residual proxies under state/source substitutions."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "variant",
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "active_edge_count",
    "active_edge_ids",
    "node_volume_m2",
    "transport_model",
    "source_policy",
    "srh_policy",
    "electron_transport_s_inv",
    "hole_transport_s_inv",
    "impact_source_s_inv",
    "srh_source_s_inv",
    "electron_residual_proxy_s_inv",
    "hole_residual_proxy_s_inv",
    "electron_residual_over_impact",
    "hole_residual_over_impact",
    "impact_over_sentaurus_node_source",
    "electron_transport_over_sentaurus_node_source",
    "hole_transport_over_sentaurus_node_source",
]

SUMMARY_KEYS = [
    "electron_residual_over_impact",
    "hole_residual_over_impact",
    "impact_over_sentaurus_node_source",
    "electron_transport_over_sentaurus_node_source",
    "hole_transport_over_sentaurus_node_source",
]

VALID_TRANSPORT_MODES = {
    "density_sg",
    "qf_inferred_ni",
    "qf_old_slotboom_ni",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--edge-local-source-csv", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--active-axis", choices=["x", "y"], default="x")
    parser.add_argument("--active-source-relative-threshold", type=float, default=1.0e-12)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    parser.add_argument("--transport-modes", default="density_sg")
    parser.add_argument("--electron-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--hole-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--qf-shift-scope", choices=["all", "support_nodes"], default="all")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with open_path(path, "r", newline="") as handle:
        return list(csv.DictReader(handle))


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def support_nodes(path: Path) -> list[dict[str, str]]:
    return [
        row for row in read_csv_rows(path)
        if row.get("support_class") not in (None, "", "inactive")
    ]


def parse_transport_modes(raw: str) -> list[str]:
    modes = [item.strip() for item in raw.split(",") if item.strip()]
    if not modes:
        raise SystemExit("--transport-modes must name at least one mode")
    invalid = [mode for mode in modes if mode not in VALID_TRANSPORT_MODES]
    if invalid:
        raise SystemExit(f"unknown transport mode(s): {', '.join(invalid)}")
    return list(dict.fromkeys(modes))


def shifted(values: list[float], delta: float, node_ids: set[int] | None = None) -> list[float]:
    if node_ids is None:
        return [value + delta for value in values]
    return [value + delta if index in node_ids else value for index, value in enumerate(values)]


def inferred_ni_state(state: dict[str, list[float]], vt: float) -> dict[str, list[float]]:
    return {
        "nie": [
            fluxforms.inferred_ni(state["n"][index], state["psi"][index], state["phin"][index], vt, "electron")
            for index in range(len(state["psi"]))
        ],
        "nih": [
            fluxforms.inferred_ni(state["p"][index], state["psi"][index], state["phip"][index], vt, "hole")
            for index in range(len(state["psi"]))
        ],
    }


def reconstruct_density_state(
    ni: dict[str, list[float]],
    psi: list[float],
    phin: list[float],
    phip: list[float],
    mobility_source: dict[str, list[float]],
    vt: float,
) -> dict[str, list[float]]:
    return {
        "psi": list(psi),
        "phin": list(phin),
        "phip": list(phip),
        "n": [
            ni["nie"][index] * fluxforms.sgdiag.limited_exp((psi[index] - phin[index]) / vt)
            for index in range(len(psi))
        ],
        "p": [
            ni["nih"][index] * fluxforms.sgdiag.limited_exp((phip[index] - psi[index]) / vt)
            for index in range(len(psi))
        ],
        "mun": list(mobility_source["mun"]),
        "mup": list(mobility_source["mup"]),
    }


def make_states(
    vela: dict[str, list[float]],
    sentaurus: dict[str, list[float]],
    vt: float,
    electron_shift: float,
    hole_shift: float,
    shift_node_ids: set[int] | None = None,
) -> dict[str, dict[str, list[float]]]:
    vela_ni = inferred_ni_state(vela, vt)
    sentaurus_with_vela_mobility = deepcopy(sentaurus)
    sentaurus_with_vela_mobility["mun"] = list(vela["mun"])
    sentaurus_with_vela_mobility["mup"] = list(vela["mup"])
    return {
        "vela_state": deepcopy(vela),
        "sentaurus_state": deepcopy(sentaurus),
        "sentaurus_state_vela_mobility": sentaurus_with_vela_mobility,
        "shifted_vela_qf": reconstruct_density_state(
            vela_ni,
            vela["psi"],
            shifted(vela["phin"], electron_shift, shift_node_ids),
            shifted(vela["phip"], hole_shift, shift_node_ids),
            vela,
            vt,
        ),
    }


def active_edge_ids(active_by_node: dict[int, list[dict[str, Any]]], node_id: int) -> str:
    return ";".join(str(int(item["edge_id"])) for item in active_by_node.get(node_id, []))


def edge_source_integral(row: dict[str, str]) -> float:
    explicit = optional_float(row.get("edge_source_integral"))
    if explicit is not None:
        return explicit
    electron = optional_float(row.get("electron_source_integral")) or 0.0
    hole = optional_float(row.get("hole_source_integral")) or 0.0
    if electron != 0.0 or hole != 0.0:
        return electron + hole
    area = optional_float(row.get("edge_area_proxy_m2")) or 0.0
    electron_flux = abs(optional_float(row.get("electron_flux_proxy")) or 0.0)
    hole_flux = abs(optional_float(row.get("hole_flux_proxy")) or 0.0)
    electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
    hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
    return area * (electron_alpha * electron_flux + hole_alpha * hole_flux)


def load_active_edges_with_source(
    edge_csv: Path,
    bias: float,
    active_axis: str,
    threshold_scale: float,
) -> dict[int, list[dict[str, Any]]]:
    by_node: dict[int, list[dict[str, Any]]] = {}
    for row in read_csv_rows(edge_csv):
        if not fluxfac.matches_bias(row, bias):
            continue
        electron_flux = abs(optional_float(row.get("electron_flux_proxy")) or 0.0)
        hole_flux = abs(optional_float(row.get("hole_flux_proxy")) or 0.0)
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        alpha_flux = electron_alpha * electron_flux + hole_alpha * hole_flux
        edge = {
            "edge_id": int(row["edge_id"]),
            "axis": fluxfac.edge_axis(row),
            "endpoint_area_m2": 0.5 * (optional_float(row.get("edge_area_proxy_m2")) or 0.0),
            "endpoint_source_integral": 0.5 * edge_source_integral(row),
            "alpha_flux": alpha_flux,
        }
        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) not in (None, ""):
                by_node.setdefault(int(row[endpoint]), []).append(edge)
    active_by_node: dict[int, list[dict[str, Any]]] = {}
    for node_id, items in by_node.items():
        max_alpha_flux = max((abs(float(item["alpha_flux"])) for item in items), default=0.0)
        threshold = max_alpha_flux * threshold_scale
        active_by_node[node_id] = [
            item for item in items
            if item["axis"] == active_axis and abs(float(item["alpha_flux"])) > threshold
        ]
    return active_by_node


def transport_terms(
    state: dict[str, list[float]],
    edges: list[dict[str, Any]],
    vt: float,
    transport_model: str,
    ni_model: list[float] | None,
) -> tuple[dict[int, float], dict[int, float]]:
    electron: dict[int, float] = defaultdict(float)
    hole: dict[int, float] = defaultdict(float)
    for edge in edges:
        i = int(edge["node0"])
        j = int(edge["node1"])
        h = float(edge["length_m"])
        if h <= 0.0:
            continue
        couple = float(edge["couple_m"])
        electron_flux, hole_flux = edge_fluxes(state, edge, vt, transport_model, ni_model)
        electron[i] += electron_flux * couple
        electron[j] -= electron_flux * couple
        hole[i] += hole_flux * couple
        hole[j] -= hole_flux * couple
    return electron, hole


def edge_fluxes(
    state: dict[str, list[float]],
    edge: dict[str, Any],
    vt: float,
    transport_model: str,
    ni_model: list[float] | None,
) -> tuple[float, float]:
    i = int(edge["node0"])
    j = int(edge["node1"])
    h = float(edge["length_m"])
    dpsi = state["psi"][j] - state["psi"][i]
    coef_n = fluxforms.avg(state["mun"], i, j) * vt / h if h > 0.0 else 0.0
    coef_p = fluxforms.avg(state["mup"], i, j) * vt / h if h > 0.0 else 0.0
    if transport_model == "density_sg":
        return (
            fluxforms.sgdiag.sg_electron_flux(state["n"][i], state["n"][j], dpsi, vt, coef_n),
            fluxforms.sgdiag.sg_hole_flux(state["p"][i], state["p"][j], dpsi, vt, coef_p),
        )
    if transport_model == "qf_inferred_ni":
        ni_e_i = fluxforms.inferred_ni(state["n"][i], state["psi"][i], state["phin"][i], vt, "electron")
        ni_e_j = fluxforms.inferred_ni(state["n"][j], state["psi"][j], state["phin"][j], vt, "electron")
        ni_h_i = fluxforms.inferred_ni(state["p"][i], state["psi"][i], state["phip"][i], vt, "hole")
        ni_h_j = fluxforms.inferred_ni(state["p"][j], state["psi"][j], state["phip"][j], vt, "hole")
    elif transport_model == "qf_old_slotboom_ni":
        if ni_model is None:
            raise RuntimeError("qf_old_slotboom_ni transport requires --doping-csv")
        ni_e_i = ni_h_i = ni_model[i]
        ni_e_j = ni_h_j = ni_model[j]
    else:
        raise RuntimeError(f"unsupported transport mode: {transport_model}")
    return (
        fluxforms.sgdiag.sg_electron_flux_qf_variable_ni(
            ni_e_i, ni_e_j, state["psi"][i], state["psi"][j], state["phin"][i], state["phin"][j], vt, coef_n),
        fluxforms.sgdiag.sg_hole_flux_qf_variable_ni(
            ni_h_i, ni_h_j, state["psi"][i], state["psi"][j], state["phip"][i], state["phip"][j], vt, coef_p),
    )


def replayed_edge_source(
    state: dict[str, list[float]],
    active_items: list[dict[str, Any]],
    mesh_edges: dict[int, dict[str, Any]],
    vt: float,
    transport_model: str,
    ni_model: list[float] | None,
) -> float:
    total = 0.0
    for item in active_items:
        edge = mesh_edges.get(int(item["edge_id"]))
        if edge is None:
            continue
        i = int(edge["node0"])
        j = int(edge["node1"])
        h = float(edge["length_m"])
        if h <= 0.0:
            continue
        electron_flux, hole_flux = edge_fluxes(state, edge, vt, transport_model, ni_model)
        electron_flux = abs(electron_flux)
        hole_flux = abs(hole_flux)
        electron_field = fluxforms.sgdiag.field_between(state["phin"], i, j, h)
        hole_field = fluxforms.sgdiag.field_between(state["phip"], i, j, h)
        alpha_e = fluxforms.sgdiag.van_overstraeten_alpha(
            electron_field, 7.03e7, 7.03e7, 1.231e8, 1.231e8)
        alpha_h = fluxforms.sgdiag.van_overstraeten_alpha(
            hole_field, 1.582e8, 6.71e7, 2.036e8, 1.693e8)
        total += (alpha_e * electron_flux + alpha_h * hole_flux) * float(item["endpoint_area_m2"])
    return total


def cxx_edge_source(active_items: list[dict[str, Any]]) -> float:
    return sum(float(item["endpoint_source_integral"]) for item in active_items)


def load_edge_local_source_ratios(path: Path | None) -> dict[int, float]:
    if path is None:
        return {}
    values: dict[int, float] = {}
    for row in read_csv_rows(path):
        if row.get("node_id") in (None, ""):
            continue
        source = optional_float(row.get("sentaurus_edgeavg_current_scaled_source"))
        cxx_source = optional_float(row.get("current_cxx_endpoint_source"))
        if source is not None and cxx_source not in (None, 0.0):
            values[int(row["node_id"])] = source / cxx_source
    return values


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize_group(rows: list[dict[str, Any]], include_transport_models: bool = False) -> dict[str, Any]:
    summary: dict[str, Any] = {"count": len(rows)}
    for key in SUMMARY_KEYS:
        values = finite_values(rows, key)
        summary[f"{key}_median"] = statistics.median(values) if values else None
        summary[f"{key}_mean"] = statistics.fmean(values) if values else None
    if include_transport_models:
        by_transport: dict[str, list[dict[str, Any]]] = {}
        for row in rows:
            by_transport.setdefault(str(row["transport_model"]), []).append(row)
        summary["transport_models"] = {
            mode: summarize_group(items)
            for mode, items in sorted(by_transport.items())
        }
    return summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
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


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, _ = fluxforms.sgdiag.read_mesh(args.mesh)
    edges = fluxforms.sgdiag.build_edges(nodes, triangles)
    mesh_edges = {int(edge["edge_id"]): edge for edge in edges}
    volumes = fluxforms.sgdiag.node_volumes(nodes, triangles)
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    transport_modes = parse_transport_modes(args.transport_modes)
    if "qf_old_slotboom_ni" in transport_modes and args.doping_csv is None:
        raise SystemExit("qf_old_slotboom_ni transport requires --doping-csv")
    ni_model = (
        fluxforms.effective_ni_from_doping(
            args.doping_csv,
            len(nodes),
            args.material_ni_m3,
            vt,
            args.bandgap_narrowing,
        )
        if args.doping_csv is not None
        else None
    )
    vela_vtk = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    vela = fluxforms.load_vela_state(vela_vtk, len(nodes))
    sentaurus = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
    support = support_nodes(args.support_csv)
    shift_node_ids = {
        int(row["node_id"]) for row in support
    } if args.qf_shift_scope == "support_nodes" else None
    states = make_states(
        vela, sentaurus, vt, args.electron_qf_shift_v, args.hole_qf_shift_v, shift_node_ids)
    transports = {
        (state_name, transport_model): transport_terms(state, edges, vt, transport_model, ni_model)
        for state_name, state in states.items()
        for transport_model in transport_modes
    }
    active_by_node = load_active_edges_with_source(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    sent_nodes = cont.load_sentaurus_nodes(args.sentaurus_dir)
    sentaurus_impact = fluxforms.read_scalar_csv(
        args.sentaurus_dir / "fields" / "ImpactIonization_region0.csv")
    sentaurus_srh = fluxforms.read_scalar_csv(
        args.sentaurus_dir / "fields" / "srhRecombination_region0.csv")
    vela_scalars = fluxforms.sgdiag.parse_vtk_scalars(vela_vtk)
    vela_srh = vela_scalars.get("SRHRecombination", [0.0 for _ in nodes])
    edge_local_source_ratios = load_edge_local_source_ratios(args.edge_local_source_csv)

    variant_specs = [
        ("vela_state_vela_edge_source", "vela_state", "vela_edge_source", "vela_srh"),
        ("vela_state_sentaurus_node_source", "vela_state", "sentaurus_node_source", "vela_srh"),
        ("sentaurus_state_sentaurus_node_source", "sentaurus_state", "sentaurus_node_source", "sentaurus_srh"),
        ("sentaurus_state_replayed_edge_source", "sentaurus_state", "replayed_edge_source", "sentaurus_srh"),
        ("shifted_vela_qf_sentaurus_node_source", "shifted_vela_qf", "sentaurus_node_source", "vela_srh"),
        ("shifted_vela_qf_replayed_edge_source", "shifted_vela_qf", "replayed_edge_source", "vela_srh"),
        (
            "probe_vela_density_vela_mobility_vela_source",
            "vela_state",
            "vela_edge_source",
            "sentaurus_srh",
        ),
        (
            "probe_sentaurus_density_vela_mobility_vela_source",
            "sentaurus_state_vela_mobility",
            "vela_edge_source",
            "sentaurus_srh",
        ),
        (
            "probe_sentaurus_density_sentaurus_mobility_vela_source",
            "sentaurus_state",
            "vela_edge_source",
            "sentaurus_srh",
        ),
    ]
    if edge_local_source_ratios:
        variant_specs.append((
            "probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source",
            "sentaurus_state",
            "edgeavg_current_source",
            "sentaurus_srh",
        ))

    rows: list[dict[str, Any]] = []
    for support_row in support:
        node_id = int(support_row["node_id"])
        node = nodes[node_id]
        sent_id, sent_distance = cont.nearest_node(sent_nodes, node["x_um"], node["y_um"])
        volume = volumes[node_id]
        active_items = active_by_node.get(node_id, [])
        vela_edge_source = cxx_edge_source(active_items)
        sentaurus_node_source = sentaurus_impact[sent_id] * 1.0e6 * volume
        edgeavg_source_ratio = edge_local_source_ratios.get(node_id)
        source_values = {
            "vela_edge_source": vela_edge_source,
            "sentaurus_node_source": sentaurus_node_source,
            "edgeavg_current_source": (
                vela_edge_source * edgeavg_source_ratio
                if edgeavg_source_ratio is not None
                else None
            ),
        }
        srh_values = {
            "vela_srh": vela_srh[node_id] * volume,
            "sentaurus_srh": sentaurus_srh[sent_id] * 1.0e6 * volume,
        }
        for transport_model in transport_modes:
            for variant, state_name, source_policy, srh_policy in variant_specs:
                state = states[state_name]
                electron_transport, hole_transport = transports[(state_name, transport_model)]
                impact_source = (
                    replayed_edge_source(state, active_items, mesh_edges, vt, transport_model, ni_model)
                    if source_policy == "replayed_edge_source"
                    else source_values[source_policy]
                )
                if impact_source is None:
                    continue
                srh_source = srh_values[srh_policy]
                electron_residual = electron_transport[node_id] + srh_source - impact_source
                hole_residual = hole_transport[node_id] + srh_source - impact_source
                rows.append({
                    "variant": variant,
                    "node_id": node_id,
                    "x_um": node["x_um"],
                    "y_um": node["y_um"],
                    "support_class": support_row.get("support_class", ""),
                    "sentaurus_node_id": sent_id,
                    "sentaurus_distance_um": sent_distance,
                    "active_edge_count": len(active_items),
                    "active_edge_ids": active_edge_ids(active_by_node, node_id),
                    "node_volume_m2": volume,
                    "transport_model": transport_model,
                    "source_policy": source_policy,
                    "srh_policy": srh_policy,
                    "electron_transport_s_inv": electron_transport[node_id],
                    "hole_transport_s_inv": hole_transport[node_id],
                    "impact_source_s_inv": impact_source,
                    "srh_source_s_inv": srh_source,
                    "electron_residual_proxy_s_inv": electron_residual,
                    "hole_residual_proxy_s_inv": hole_residual,
                    "electron_residual_over_impact": ratio(electron_residual, impact_source),
                    "hole_residual_over_impact": ratio(hole_residual, impact_source),
                    "impact_over_sentaurus_node_source": ratio(impact_source, sentaurus_node_source),
                    "electron_transport_over_sentaurus_node_source": ratio(
                        electron_transport[node_id], sentaurus_node_source),
                    "hole_transport_over_sentaurus_node_source": ratio(
                        hole_transport[node_id], sentaurus_node_source),
                })
    return rows


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_support_residual_proxy_nodes.csv", rows)
    by_class: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_class.setdefault(str(row["support_class"]), []).append(row)
    summary = {
        "bias": args.bias,
        "row_count": len(rows),
        "variants": sorted({str(row["variant"]) for row in rows}),
        "transport_modes": sorted({str(row["transport_model"]) for row in rows}),
        "support_classes": {},
    }
    for support_class, class_rows in by_class.items():
        by_variant: dict[str, list[dict[str, Any]]] = {}
        for row in class_rows:
            by_variant.setdefault(str(row["variant"]), []).append(row)
        summary["support_classes"][support_class] = {
            "count": len(class_rows),
            "variants": {
                variant: summarize_group(items, include_transport_models=True)
                for variant, items in sorted(by_variant.items())
            },
        }
    payload = clean_json(summary)
    with open_path(args.out_dir / "active_support_residual_proxy_summary.json", "w", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
