#!/usr/bin/env python3
"""Compare Vela and Sentaurus states on selected BV transition path edges."""

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
    "carrier",
    "start_node",
    "target_contact",
    "step",
    "edge_id",
    "node_from",
    "node_to",
    "x_from_um",
    "x_to_um",
    "edge_mid_x_um",
    "length_m",
    "sentaurus_node_from",
    "sentaurus_node_to",
    "sentaurus_distance_from_um",
    "sentaurus_distance_to_um",
    "vela_psi_from_V",
    "vela_psi_to_V",
    "sentaurus_psi_from_V",
    "sentaurus_psi_to_V",
    "delta_psi_from_V",
    "delta_psi_to_V",
    "vela_phin_from_V",
    "vela_phin_to_V",
    "sentaurus_phin_from_V",
    "sentaurus_phin_to_V",
    "delta_phin_from_V",
    "delta_phin_to_V",
    "vela_phip_from_V",
    "vela_phip_to_V",
    "sentaurus_phip_from_V",
    "sentaurus_phip_to_V",
    "delta_phip_from_V",
    "delta_phip_to_V",
    "vela_psi_minus_phin_from_V",
    "vela_psi_minus_phin_to_V",
    "sentaurus_psi_minus_phin_from_V",
    "sentaurus_psi_minus_phin_to_V",
    "delta_psi_minus_phin_avg_V",
    "vela_phip_minus_psi_from_V",
    "vela_phip_minus_psi_to_V",
    "sentaurus_phip_minus_psi_from_V",
    "sentaurus_phip_minus_psi_to_V",
    "delta_phip_minus_psi_avg_V",
    "hybrid_vela_psi_sentaurus_phin_delta_exp_avg_V",
    "hybrid_sentaurus_psi_vela_phin_delta_exp_avg_V",
    "hybrid_vela_phip_sentaurus_psi_delta_exp_avg_V",
    "hybrid_sentaurus_phip_vela_psi_delta_exp_avg_V",
    "vela_psi_drop_V",
    "sentaurus_psi_drop_V",
    "delta_psi_drop_V",
    "vela_phin_drop_V",
    "sentaurus_phin_drop_V",
    "delta_phin_drop_V",
    "vela_phip_drop_V",
    "sentaurus_phip_drop_V",
    "delta_phip_drop_V",
    "vela_electric_field_V_m",
    "sentaurus_electric_field_V_m",
    "delta_electric_field_V_m",
    "vela_electron_qf_field_V_m",
    "sentaurus_electron_qf_field_V_m",
    "delta_electron_qf_field_V_m",
    "vela_hole_qf_field_V_m",
    "sentaurus_hole_qf_field_V_m",
    "delta_hole_qf_field_V_m",
    "model_ni_eff_avg_cm3",
    "vela_electron_density_avg_cm3",
    "sentaurus_electron_density_avg_cm3",
    "log10_vela_over_sentaurus_electron_density",
    "vela_hole_density_avg_cm3",
    "sentaurus_hole_density_avg_cm3",
    "log10_vela_over_sentaurus_hole_density",
    "dominant_edge_delta_metric",
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
    parser.add_argument("--carrier", choices=["all", "electron", "hole"], default="all")
    parser.add_argument("--start-node", type=int)
    parser.add_argument("--target-contact", default="")
    parser.add_argument("--edge-ids", default="")
    parser.add_argument("--x-min-um", type=float)
    parser.add_argument("--x-max-um", type=float)
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


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def log10_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate <= 0.0 or reference <= 0.0:
        return None
    return math.log10(candidate / reference)


def avg2(a: float, b: float) -> float:
    return 0.5 * (a + b)


def field(drop: float, length_m: float) -> float:
    if length_m <= 1.0e-30:
        return 0.0
    return abs(drop) / length_m


def nearest_cache(sentaurus_nodes: dict[int, dict[str, float]], nodes: dict[int, dict[str, float]]) -> dict[int, tuple[int, float]]:
    return {
        node_id: cont.nearest_node(sentaurus_nodes, node["x_um"], node["y_um"])
        for node_id, node in nodes.items()
    }


def load_sentaurus_fields(root: Path) -> dict[str, dict[int, list[float]]]:
    fields = {
        name: cont.load_components(root, name)
        for name in ("ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential")
    }
    for name in ("eDensity", "hDensity"):
        try:
            fields[name] = cont.load_components(root, name)
        except FileNotFoundError:
            fields[name] = {}
    return fields


def selected_path_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    edge_ids = parse_edge_ids(args.edge_ids)
    rows: list[dict[str, str]] = []
    for row in read_csv_rows(args.path_edge_csv):
        if abs(float(row["bias_V"]) - args.bias) > 1.0e-9:
            continue
        if args.carrier != "all" and row.get("carrier") != args.carrier:
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
    rows.sort(key=lambda item: (item.get("carrier", ""), int(item.get("start_node", 0)), int(item.get("step", 0))))
    return rows


def scalar_pair(values: list[float], i: int, j: int) -> tuple[float, float]:
    return values[i], values[j]


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


def sent_scalar(fields: dict[str, dict[int, list[float]]], name: str, node_id: int) -> float | None:
    if node_id not in fields.get(name, {}):
        return None
    return cont.scalar_at(fields, name, node_id)


def dominant_metric(row: dict[str, Any]) -> str:
    metrics = {
        "psi_drop": abs(float(row["delta_psi_drop_V"])),
        "electron_qf_drop": abs(float(row["delta_phin_drop_V"])),
        "hole_qf_drop": abs(float(row["delta_phip_drop_V"])),
        "electron_exponent": abs(float(row["delta_psi_minus_phin_avg_V"])),
        "hole_exponent": abs(float(row["delta_phip_minus_psi_avg_V"])),
    }
    return max(metrics, key=metrics.get)


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, _, _ = cont.sgdiag.read_mesh(args.mesh)
    vt = cont.sgdiag.K_B_OVER_Q * args.temperature_k
    vtk_path = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    scalars = cont.sgdiag.parse_vtk_scalars(vtk_path)
    require_vtk_scalars(
        scalars,
        ("Potential", "ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"),
        max(nodes) + 1 if nodes else 0,
        vtk_path,
    )
    sentaurus_nodes = cont.load_sentaurus_nodes(args.sentaurus_dir)
    sent_nearest = nearest_cache(sentaurus_nodes, nodes)
    sent_fields = load_sentaurus_fields(args.sentaurus_dir)
    ni = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )

    result: list[dict[str, Any]] = []
    selected = selected_path_rows(args)
    if not selected:
        raise SystemExit(
            "no transition path edges selected; check bias/contact/start-node/edge-id/x-window filters"
        )
    for path_row in selected:
        i = int(path_row["node_from"])
        j = int(path_row["node_to"])
        si, si_dist = sent_nearest[i]
        sj, sj_dist = sent_nearest[j]
        length_m = optional_float(path_row.get("length_m"))
        if length_m is None:
            length_m = math.hypot(
                (nodes[j]["x_um"] - nodes[i]["x_um"]) * 1.0e-6,
                (nodes[j]["y_um"] - nodes[i]["y_um"]) * 1.0e-6,
            )

        vela_psi_i, vela_psi_j = scalar_pair(scalars["Potential"], i, j)
        vela_phin_i, vela_phin_j = scalar_pair(scalars["ElectronQuasiFermi"], i, j)
        vela_phip_i, vela_phip_j = scalar_pair(scalars["HoleQuasiFermi"], i, j)
        sent_psi_i = cont.scalar_at(sent_fields, "ElectrostaticPotential", si)
        sent_psi_j = cont.scalar_at(sent_fields, "ElectrostaticPotential", sj)
        sent_phin_i = cont.scalar_at(sent_fields, "eQuasiFermiPotential", si)
        sent_phin_j = cont.scalar_at(sent_fields, "eQuasiFermiPotential", sj)
        sent_phip_i = cont.scalar_at(sent_fields, "hQuasiFermiPotential", si)
        sent_phip_j = cont.scalar_at(sent_fields, "hQuasiFermiPotential", sj)

        vela_psi_drop = vela_psi_j - vela_psi_i
        sent_psi_drop = sent_psi_j - sent_psi_i
        vela_phin_drop = vela_phin_j - vela_phin_i
        sent_phin_drop = sent_phin_j - sent_phin_i
        vela_phip_drop = vela_phip_j - vela_phip_i
        sent_phip_drop = sent_phip_j - sent_phip_i
        avg_delta_psi = avg2(vela_psi_i - sent_psi_i, vela_psi_j - sent_psi_j)
        avg_delta_phin = avg2(vela_phin_i - sent_phin_i, vela_phin_j - sent_phin_j)
        avg_delta_phip = avg2(vela_phip_i - sent_phip_i, vela_phip_j - sent_phip_j)

        vela_n_avg = avg2(scalars["Electrons"][i], scalars["Electrons"][j]) * 1.0e-6
        vela_p_avg = avg2(scalars["Holes"][i], scalars["Holes"][j]) * 1.0e-6
        sent_n_i = sent_scalar(sent_fields, "eDensity", si)
        sent_n_j = sent_scalar(sent_fields, "eDensity", sj)
        sent_p_i = sent_scalar(sent_fields, "hDensity", si)
        sent_p_j = sent_scalar(sent_fields, "hDensity", sj)
        sent_n_avg = avg2(sent_n_i, sent_n_j) if sent_n_i is not None and sent_n_j is not None else None
        sent_p_avg = avg2(sent_p_i, sent_p_j) if sent_p_i is not None and sent_p_j is not None else None

        row: dict[str, Any] = {
            "bias_V": args.bias,
            "carrier": path_row.get("carrier", ""),
            "start_node": path_row.get("start_node", ""),
            "target_contact": path_row.get("target_contact", ""),
            "step": path_row.get("step", ""),
            "edge_id": path_row["edge_id"],
            "node_from": i,
            "node_to": j,
            "x_from_um": nodes[i]["x_um"],
            "x_to_um": nodes[j]["x_um"],
            "edge_mid_x_um": avg2(nodes[i]["x_um"], nodes[j]["x_um"]),
            "length_m": length_m,
            "sentaurus_node_from": si,
            "sentaurus_node_to": sj,
            "sentaurus_distance_from_um": si_dist,
            "sentaurus_distance_to_um": sj_dist,
            "vela_psi_from_V": vela_psi_i,
            "vela_psi_to_V": vela_psi_j,
            "sentaurus_psi_from_V": sent_psi_i,
            "sentaurus_psi_to_V": sent_psi_j,
            "delta_psi_from_V": vela_psi_i - sent_psi_i,
            "delta_psi_to_V": vela_psi_j - sent_psi_j,
            "vela_phin_from_V": vela_phin_i,
            "vela_phin_to_V": vela_phin_j,
            "sentaurus_phin_from_V": sent_phin_i,
            "sentaurus_phin_to_V": sent_phin_j,
            "delta_phin_from_V": vela_phin_i - sent_phin_i,
            "delta_phin_to_V": vela_phin_j - sent_phin_j,
            "vela_phip_from_V": vela_phip_i,
            "vela_phip_to_V": vela_phip_j,
            "sentaurus_phip_from_V": sent_phip_i,
            "sentaurus_phip_to_V": sent_phip_j,
            "delta_phip_from_V": vela_phip_i - sent_phip_i,
            "delta_phip_to_V": vela_phip_j - sent_phip_j,
            "vela_psi_minus_phin_from_V": vela_psi_i - vela_phin_i,
            "vela_psi_minus_phin_to_V": vela_psi_j - vela_phin_j,
            "sentaurus_psi_minus_phin_from_V": sent_psi_i - sent_phin_i,
            "sentaurus_psi_minus_phin_to_V": sent_psi_j - sent_phin_j,
            "delta_psi_minus_phin_avg_V": (
                avg2(vela_psi_i - vela_phin_i, vela_psi_j - vela_phin_j)
                - avg2(sent_psi_i - sent_phin_i, sent_psi_j - sent_phin_j)
            ),
            "vela_phip_minus_psi_from_V": vela_phip_i - vela_psi_i,
            "vela_phip_minus_psi_to_V": vela_phip_j - vela_psi_j,
            "sentaurus_phip_minus_psi_from_V": sent_phip_i - sent_psi_i,
            "sentaurus_phip_minus_psi_to_V": sent_phip_j - sent_psi_j,
            "delta_phip_minus_psi_avg_V": (
                avg2(vela_phip_i - vela_psi_i, vela_phip_j - vela_psi_j)
                - avg2(sent_phip_i - sent_psi_i, sent_phip_j - sent_psi_j)
            ),
            "hybrid_vela_psi_sentaurus_phin_delta_exp_avg_V": avg_delta_psi,
            "hybrid_sentaurus_psi_vela_phin_delta_exp_avg_V": -avg_delta_phin,
            "hybrid_vela_phip_sentaurus_psi_delta_exp_avg_V": avg_delta_phip,
            "hybrid_sentaurus_phip_vela_psi_delta_exp_avg_V": -avg_delta_psi,
            "vela_psi_drop_V": vela_psi_drop,
            "sentaurus_psi_drop_V": sent_psi_drop,
            "delta_psi_drop_V": vela_psi_drop - sent_psi_drop,
            "vela_phin_drop_V": vela_phin_drop,
            "sentaurus_phin_drop_V": sent_phin_drop,
            "delta_phin_drop_V": vela_phin_drop - sent_phin_drop,
            "vela_phip_drop_V": vela_phip_drop,
            "sentaurus_phip_drop_V": sent_phip_drop,
            "delta_phip_drop_V": vela_phip_drop - sent_phip_drop,
            "vela_electric_field_V_m": field(vela_psi_drop, length_m),
            "sentaurus_electric_field_V_m": field(sent_psi_drop, length_m),
            "delta_electric_field_V_m": field(vela_psi_drop, length_m) - field(sent_psi_drop, length_m),
            "vela_electron_qf_field_V_m": field(vela_phin_drop, length_m),
            "sentaurus_electron_qf_field_V_m": field(sent_phin_drop, length_m),
            "delta_electron_qf_field_V_m": field(vela_phin_drop, length_m) - field(sent_phin_drop, length_m),
            "vela_hole_qf_field_V_m": field(vela_phip_drop, length_m),
            "sentaurus_hole_qf_field_V_m": field(sent_phip_drop, length_m),
            "delta_hole_qf_field_V_m": field(vela_phip_drop, length_m) - field(sent_phip_drop, length_m),
            "model_ni_eff_avg_cm3": avg2(ni[i], ni[j]) * 1.0e-6,
            "vela_electron_density_avg_cm3": vela_n_avg,
            "sentaurus_electron_density_avg_cm3": sent_n_avg,
            "log10_vela_over_sentaurus_electron_density": log10_ratio(vela_n_avg, sent_n_avg),
            "vela_hole_density_avg_cm3": vela_p_avg,
            "sentaurus_hole_density_avg_cm3": sent_p_avg,
            "log10_vela_over_sentaurus_hole_density": log10_ratio(vela_p_avg, sent_p_avg),
        }
        row["dominant_edge_delta_metric"] = dominant_metric(row)
        result.append(row)
    return result


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def median_abs(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [abs(value) for value in finite_values(rows, key)]
    return statistics.median(values) if values else None


def make_summary(args: argparse.Namespace, rows: list[dict[str, Any]]) -> dict[str, Any]:
    metrics = {
        "psi_drop": median_abs(rows, "delta_psi_drop_V"),
        "electron_qf_drop": median_abs(rows, "delta_phin_drop_V"),
        "hole_qf_drop": median_abs(rows, "delta_phip_drop_V"),
        "electron_exponent_avg": median_abs(rows, "delta_psi_minus_phin_avg_V"),
        "hole_exponent_avg": median_abs(rows, "delta_phip_minus_psi_avg_V"),
        "hybrid_vela_psi_sentaurus_phin": median_abs(rows, "hybrid_vela_psi_sentaurus_phin_delta_exp_avg_V"),
        "hybrid_sentaurus_psi_vela_phin": median_abs(rows, "hybrid_sentaurus_psi_vela_phin_delta_exp_avg_V"),
        "hybrid_vela_phip_sentaurus_psi": median_abs(rows, "hybrid_vela_phip_sentaurus_psi_delta_exp_avg_V"),
        "hybrid_sentaurus_phip_vela_psi": median_abs(rows, "hybrid_sentaurus_phip_vela_psi_delta_exp_avg_V"),
        "electric_field": median_abs(rows, "delta_electric_field_V_m"),
        "electron_qf_field": median_abs(rows, "delta_electron_qf_field_V_m"),
        "hole_qf_field": median_abs(rows, "delta_hole_qf_field_V_m"),
    }
    finite_metrics = {key: value for key, value in metrics.items() if value is not None}
    dominant = max(finite_metrics, key=finite_metrics.get) if finite_metrics else ""
    return {
        "bias_V": args.bias,
        "path_edge_csv": str(args.path_edge_csv),
        "edge_count": len(rows),
        "carrier_filter": args.carrier,
        "start_node_filter": args.start_node,
        "target_contact_filter": args.target_contact,
        "x_min_um": args.x_min_um,
        "x_max_um": args.x_max_um,
        "median_abs_metrics": finite_metrics,
        "dominant_drop_mismatch_metric": dominant,
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
    write_csv(args.out_dir / "transition_edge_state_compare.csv", rows)
    summary = make_summary(args, rows)
    (args.out_dir / "transition_edge_state_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
