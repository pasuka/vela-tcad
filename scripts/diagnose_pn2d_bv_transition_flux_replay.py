#!/usr/bin/env python3
"""Replay electron SG flux on transition edges with Vela/Sentaurus hybrid states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


FIELDS = [
    "bias_V",
    "state",
    "start_node",
    "target_contact",
    "step",
    "edge_id",
    "node_from",
    "node_to",
    "edge_mid_x_um",
    "length_m",
    "couple_m",
    "psi_source",
    "phin_source",
    "psi_from_V",
    "psi_to_V",
    "phin_from_V",
    "phin_to_V",
    "psi_drop_V",
    "phin_drop_V",
    "electron_density_exp_from",
    "electron_density_exp_to",
    "electron_qf_field_V_m",
    "electron_flux_signed_m2_s",
    "electron_flux_integral_s_inv",
    "electron_flux_integral_delta_vs_sentaurus_s_inv",
    "electron_flux_integral_ratio_vs_sentaurus",
    "d_flux_from_s_inv_per_V",
    "d_flux_to_s_inv_per_V",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--start-node", type=int)
    parser.add_argument("--target-contact", default="")
    parser.add_argument("--edge-ids", default="")
    parser.add_argument("--x-min-um", type=float)
    parser.add_argument("--x-max-um", type=float)
    parser.add_argument("--delta-v", type=float, default=1.0e-4)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_edge_ids(raw: str) -> set[int] | None:
    if not raw.strip():
        return None
    return {int(item.strip()) for item in raw.split(",") if item.strip()}


def selected_path_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    edge_ids = parse_edge_ids(args.edge_ids)
    rows: list[dict[str, str]] = []
    for row in read_csv_rows(args.path_edge_csv):
        if abs(float(row["bias_V"]) - args.bias) > 1.0e-9:
            continue
        if row.get("carrier") != "electron":
            continue
        if args.start_node is not None and int(row["start_node"]) != args.start_node:
            continue
        if args.target_contact and row.get("target_contact") != args.target_contact:
            continue
        if edge_ids is not None and int(row["edge_id"]) not in edge_ids:
            continue
        xmin = min(float(row["x_from_um"]), float(row["x_to_um"]))
        xmax = max(float(row["x_from_um"]), float(row["x_to_um"]))
        if args.x_min_um is not None and xmax < args.x_min_um:
            continue
        if args.x_max_um is not None and xmin > args.x_max_um:
            continue
        rows.append(row)
    rows.sort(key=lambda item: (int(item.get("start_node", 0)), int(item.get("step", 0))))
    return rows


def require_selected_path_rows(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise SystemExit(
            "no transition path edges selected; check bias/contact/start-node/edge-id/x-window filters"
        )


def require_vtk_scalars(scalars: dict[str, list[float]], required_names: tuple[str, ...], min_len: int, source: Path) -> None:
    missing = [name for name in required_names if name not in scalars]
    if missing:
        raise RuntimeError(f"VTK is missing required scalars {', '.join(missing)}: {source}")
    too_short = [name for name in required_names if len(scalars[name]) < min_len]
    if too_short:
        raise RuntimeError(
            f"VTK scalar length is too short for mesh node ids ({min_len} required): "
            f"{', '.join(too_short)} in {source}"
        )


def nearest_cache(sentaurus_nodes: dict[int, dict[str, float]], nodes: dict[int, dict[str, float]]) -> dict[int, tuple[int, float]]:
    return {
        node_id: cont.nearest_node(sentaurus_nodes, node["x_um"], node["y_um"])
        for node_id, node in nodes.items()
    }


def load_sentaurus_scalar(root: Path, name: str) -> dict[int, float]:
    raw = cont.load_components(root, name)
    return {node_id: values[0] for node_id, values in raw.items()}


def avg2(a: float, b: float) -> float:
    return 0.5 * (a + b)


def field_abs(drop: float, length_m: float) -> float:
    return abs(drop) / length_m if length_m > 1.0e-30 else 0.0


def safe_ratio(candidate: float, reference: float) -> float | None:
    if reference == 0.0:
        return None
    return candidate / reference


def flux_value(
    i: int,
    j: int,
    length_m: float,
    mobility: float,
    psi: list[float],
    phin: list[float],
    ni: list[float],
    vt: float,
) -> float:
    coef = mobility * vt / length_m if length_m > 1.0e-30 else 0.0
    return cont.sgdiag.sg_electron_flux_qf_variable_ni(
        ni[i], ni[j], psi[i], psi[j], phin[i], phin[j], vt, coef)


def derivative_at(
    node_id: int,
    i: int,
    j: int,
    length_m: float,
    couple_m: float,
    mobility: float,
    psi: list[float],
    phin: list[float],
    ni: list[float],
    vt: float,
    delta_v: float,
) -> float:
    plus = list(phin)
    minus = list(phin)
    plus[node_id] += delta_v
    minus[node_id] -= delta_v
    f_plus = flux_value(i, j, length_m, mobility, psi, plus, ni, vt)
    f_minus = flux_value(i, j, length_m, mobility, psi, minus, ni, vt)
    return (f_plus - f_minus) * couple_m / (2.0 * delta_v)


def make_state_vectors(
    nodes: dict[int, dict[str, float]],
    scalars: dict[str, list[float]],
    sent_nearest: dict[int, tuple[int, float]],
    sent_psi: dict[int, float],
    sent_phin: dict[int, float],
) -> dict[str, dict[str, list[float] | str]]:
    count = max(nodes) + 1
    psi_sent = [0.0 for _ in range(count)]
    phin_sent = [0.0 for _ in range(count)]
    for node_id in nodes:
        sent_id, _ = sent_nearest[node_id]
        psi_sent[node_id] = sent_psi[sent_id]
        phin_sent[node_id] = sent_phin[sent_id]
    psi_vela = scalars["Potential"]
    phin_vela = scalars["ElectronQuasiFermi"]
    return {
        "vela": {"psi": psi_vela, "phin": phin_vela, "psi_source": "vela", "phin_source": "vela"},
        "vela_psi_sentaurus_phin": {
            "psi": psi_vela,
            "phin": phin_sent,
            "psi_source": "vela",
            "phin_source": "sentaurus",
        },
        "sentaurus_psi_vela_phin": {
            "psi": psi_sent,
            "phin": phin_vela,
            "psi_source": "sentaurus",
            "phin_source": "vela",
        },
        "sentaurus": {
            "psi": psi_sent,
            "phin": phin_sent,
            "psi_source": "sentaurus",
            "phin_source": "sentaurus",
        },
    }


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, _, _ = cont.sgdiag.read_mesh(args.mesh)
    vt = cont.sgdiag.K_B_OVER_Q * args.temperature_k
    vtk_path = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    scalars = cont.sgdiag.parse_vtk_scalars(vtk_path)
    require_vtk_scalars(
        scalars,
        ("Potential", "ElectronQuasiFermi", "ElectronMobility"),
        max(nodes) + 1 if nodes else 0,
        vtk_path,
    )
    sentaurus_nodes = cont.load_sentaurus_nodes(args.sentaurus_dir)
    sent_nearest = nearest_cache(sentaurus_nodes, nodes)
    sent_psi = load_sentaurus_scalar(args.sentaurus_dir, "ElectrostaticPotential")
    sent_phin = load_sentaurus_scalar(args.sentaurus_dir, "eQuasiFermiPotential")
    ni = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    states = make_state_vectors(nodes, scalars, sent_nearest, sent_psi, sent_phin)
    selected = selected_path_rows(args)
    require_selected_path_rows(selected)
    rows: list[dict[str, Any]] = []
    for path_row in selected:
        i = int(path_row["node_from"])
        j = int(path_row["node_to"])
        length_m = float(path_row.get("length_m") or 0.0)
        if length_m <= 0.0:
            length_m = math.hypot(
                (nodes[j]["x_um"] - nodes[i]["x_um"]) * 1.0e-6,
                (nodes[j]["y_um"] - nodes[i]["y_um"]) * 1.0e-6,
            )
        couple_m = float(path_row.get("couple_m") or 1.0)
        mobility = avg2(scalars["ElectronMobility"][i], scalars["ElectronMobility"][j])
        flux_by_state: dict[str, float] = {}
        for name, state in states.items():
            psi = state["psi"]
            phin = state["phin"]
            assert isinstance(psi, list)
            assert isinstance(phin, list)
            flux_by_state[name] = flux_value(i, j, length_m, mobility, psi, phin, ni, vt) * couple_m
        reference = flux_by_state["sentaurus"]
        for name, state in states.items():
            psi = state["psi"]
            phin = state["phin"]
            assert isinstance(psi, list)
            assert isinstance(phin, list)
            flux = flux_value(i, j, length_m, mobility, psi, phin, ni, vt)
            integral = flux * couple_m
            psi_drop = psi[j] - psi[i]
            phin_drop = phin[j] - phin[i]
            rows.append({
                "bias_V": args.bias,
                "state": name,
                "start_node": path_row.get("start_node", ""),
                "target_contact": path_row.get("target_contact", ""),
                "step": path_row.get("step", ""),
                "edge_id": path_row["edge_id"],
                "node_from": i,
                "node_to": j,
                "edge_mid_x_um": avg2(nodes[i]["x_um"], nodes[j]["x_um"]),
                "length_m": length_m,
                "couple_m": couple_m,
                "psi_source": state["psi_source"],
                "phin_source": state["phin_source"],
                "psi_from_V": psi[i],
                "psi_to_V": psi[j],
                "phin_from_V": phin[i],
                "phin_to_V": phin[j],
                "psi_drop_V": psi_drop,
                "phin_drop_V": phin_drop,
                "electron_density_exp_from": (psi[i] - phin[i]) / vt,
                "electron_density_exp_to": (psi[j] - phin[j]) / vt,
                "electron_qf_field_V_m": field_abs(phin_drop, length_m),
                "electron_flux_signed_m2_s": flux,
                "electron_flux_integral_s_inv": integral,
                "electron_flux_integral_delta_vs_sentaurus_s_inv": integral - reference,
                "electron_flux_integral_ratio_vs_sentaurus": safe_ratio(integral, reference),
                "d_flux_from_s_inv_per_V": derivative_at(
                    i, i, j, length_m, couple_m, mobility, psi, phin, ni, vt, args.delta_v),
                "d_flux_to_s_inv_per_V": derivative_at(
                    j, i, j, length_m, couple_m, mobility, psi, phin, ni, vt, args.delta_v),
            })
    return rows


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_state: dict[str, list[float]] = {}
    for row in rows:
        if row["state"] == "sentaurus":
            continue
        value = optional_float(row.get("electron_flux_integral_delta_vs_sentaurus_s_inv"))
        if value is not None:
            by_state.setdefault(str(row["state"]), []).append(abs(value))
    median_by_state = {
        state: statistics.median(values)
        for state, values in by_state.items()
        if values
    }
    best = min(median_by_state, key=median_by_state.get) if median_by_state else ""
    return {
        "edge_count": len({(row["edge_id"], row["start_node"]) for row in rows}),
        "row_count": len(rows),
        "median_abs_flux_delta_vs_sentaurus_s_inv": median_by_state,
        "best_flux_match_state": best,
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
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "transition_flux_replay_edges.csv", rows)
    summary = summarize(rows)
    summary.update({
        "bias_V": args.bias,
        "path_edge_csv": str(args.path_edge_csv),
        "x_min_um": args.x_min_um,
        "x_max_um": args.x_max_um,
        "target_contact": args.target_contact,
        "start_node": args.start_node,
    })
    (args.out_dir / "transition_flux_replay_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
