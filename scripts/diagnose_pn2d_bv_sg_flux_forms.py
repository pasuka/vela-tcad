#!/usr/bin/env python3
"""Compare PN2D BV SG flux forms on selected edges for Vela and Sentaurus states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


FIELDS = [
    "bias_V",
    "edge_id",
    "source",
    "node0",
    "node1",
    "x0_um",
    "y0_um",
    "x1_um",
    "y1_um",
    "length_m",
    "couple_m",
    "edge_area_m2",
    "edge_class",
    "electric_field_V_m",
    "electron_qf_field_V_m",
    "hole_qf_field_V_m",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
    "electron_density_0_m3",
    "electron_density_1_m3",
    "hole_density_0_m3",
    "hole_density_1_m3",
    "psi0_V",
    "psi1_V",
    "phin0_V",
    "phin1_V",
    "phip0_V",
    "phip1_V",
    "ni_model_0_m3",
    "ni_model_1_m3",
    "ni_electron_inferred_0_m3",
    "ni_electron_inferred_1_m3",
    "ni_hole_inferred_0_m3",
    "ni_hole_inferred_1_m3",
    "electron_flux_density",
    "electron_flux_density_abs",
    "electron_flux_qf_model",
    "electron_flux_qf_model_abs",
    "electron_flux_qf_inferred",
    "electron_flux_qf_inferred_abs",
    "electron_qf_model_over_density_abs",
    "electron_qf_inferred_over_density_abs",
    "hole_flux_density",
    "hole_flux_density_abs",
    "hole_flux_qf_model",
    "hole_flux_qf_model_abs",
    "hole_flux_qf_inferred",
    "hole_flux_qf_inferred_abs",
    "electron_alpha_qf_m_inv",
    "electron_alpha_efield_m_inv",
    "hole_alpha_qf_m_inv",
    "hole_alpha_efield_m_inv",
    "electron_source_qf_model_qfdrive",
    "electron_source_qf_model_efielddrive",
    "electron_source_density_qfdrive",
    "electron_source_density_efielddrive",
]


def main() -> int:
    args = parse_args()
    biases = parse_float_list(args.biases)
    edge_ids = parse_int_list(args.edge_ids)
    nodes, triangles, contact_by_node = sgdiag.read_mesh(args.mesh)
    edges = {int(edge["edge_id"]): edge for edge in sgdiag.build_edges(nodes, triangles)}
    missing_edges = [edge_id for edge_id in edge_ids if edge_id not in edges]
    if missing_edges:
        raise SystemExit(f"missing edge ids: {missing_edges}")
    ni_model = effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        sgdiag.K_B_OVER_Q * args.temperature_k,
        args.bandgap_narrowing,
    )
    vela_vtks = discover_vela_vtks(args.vela_vtk_root)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        vtk = vela_vtks.get(bias_key(bias))
        if vtk is None:
            raise SystemExit(f"missing Vela VTK for {bias:g} V")
        sentaurus = sentaurus_dir(args.sentaurus_root, bias)
        states = [
            ("vela", load_vela_state(vtk, len(nodes))),
            ("sentaurus", load_sentaurus_state(sentaurus, len(nodes))),
        ]
        for source, state in states:
            for edge_id in edge_ids:
                rows.append(edge_row(
                    bias,
                    source,
                    edges[edge_id],
                    nodes,
                    contact_by_node,
                    state,
                    ni_model,
                    sgdiag.K_B_OVER_Q * args.temperature_k,
                ))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sg_flux_form_edges.csv", rows, FIELDS)
    summary = {
        "biases": biases,
        "edge_ids": edge_ids,
        "row_count": len(rows),
        "mesh": str(args.mesh),
        "sentaurus_root": str(args.sentaurus_root),
        "vela_vtk_root": str(args.vela_vtk_root),
        "ni_model": {
            "source": str(args.doping_csv),
            "material_ni_m3": args.material_ni_m3,
            "bandgap_narrowing": args.bandgap_narrowing,
        },
    }
    (args.out_dir / "sg_flux_form_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", required=True)
    parser.add_argument("--edge-ids", required=True)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def bias_token(bias: float) -> str:
    return f"{bias:g}".replace("-", "m").replace("+", "").replace(".", "p")


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    mtimes: dict[float, float] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if not match:
            continue
        key = bias_key(float(match.group("bias")))
        mtime = path.stat().st_mtime
        if key not in result or mtime >= mtimes[key]:
            result[key] = path
            mtimes[key] = mtime
    return result


def sentaurus_dir(root: Path, bias: float) -> Path:
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias_token(bias)}v",
        root / f"sentaurus_{bias_token(abs(bias))}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "fields").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def read_total_impurity_m3(path: Path, node_count: int) -> list[float]:
    values = [0.0 for _ in range(node_count)]
    seen = [False for _ in range(node_count)]
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            node_id = int(row["node_id"])
            if "donors_m3" in row or "acceptors_m3" in row:
                donors = float(row.get("donors_m3", 0.0) or 0.0)
                acceptors = float(row.get("acceptors_m3", 0.0) or 0.0)
            else:
                donors = float(row.get("donors_cm3", 0.0) or 0.0) * 1.0e6
                acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0) * 1.0e6
            values[node_id] = abs(donors) + abs(acceptors)
            seen[node_id] = True
    missing = [index for index, present in enumerate(seen) if not present]
    if missing:
        raise RuntimeError(f"doping CSV is missing node ids: {missing[:8]}")
    return values


def effective_ni_from_doping(
    path: Path,
    node_count: int,
    material_ni_m3: float,
    vt: float,
    bandgap_narrowing: str,
) -> list[float]:
    impurity = read_total_impurity_m3(path, node_count)
    values: list[float] = []
    for total in impurity:
        delta = sgdiag.old_slotboom_delta_eg_eV(total) if bandgap_narrowing == "old_slotboom" else 0.0
        values.append(material_ni_m3 * math.exp(delta / (2.0 * vt)) if delta > 0.0 else material_ni_m3)
    return values


def load_vela_state(path: Path, node_count: int) -> dict[str, list[float]]:
    scalars = sgdiag.parse_vtk_scalars(path)
    mapping = {
        "psi": ("Potential", 1.0),
        "phin": ("ElectronQuasiFermi", 1.0),
        "phip": ("HoleQuasiFermi", 1.0),
        "n": ("Electrons", 1.0),
        "p": ("Holes", 1.0),
        "mun": ("ElectronMobility", 1.0),
        "mup": ("HoleMobility", 1.0),
    }
    return mapped_state(scalars, mapping, node_count)


def load_sentaurus_state(root: Path, node_count: int) -> dict[str, list[float]]:
    fields = {
        "ElectrostaticPotential": read_scalar_csv(root / "fields" / "ElectrostaticPotential_region0.csv"),
        "eQuasiFermiPotential": read_scalar_csv(root / "fields" / "eQuasiFermiPotential_region0.csv"),
        "hQuasiFermiPotential": read_scalar_csv(root / "fields" / "hQuasiFermiPotential_region0.csv"),
        "eDensity": read_scalar_csv(root / "fields" / "eDensity_region0.csv"),
        "hDensity": read_scalar_csv(root / "fields" / "hDensity_region0.csv"),
        "eMobility": read_scalar_csv(root / "fields" / "eMobility_region0.csv"),
        "hMobility": read_scalar_csv(root / "fields" / "hMobility_region0.csv"),
    }
    mapping = {
        "psi": ("ElectrostaticPotential", 1.0),
        "phin": ("eQuasiFermiPotential", 1.0),
        "phip": ("hQuasiFermiPotential", 1.0),
        "n": ("eDensity", 1.0e6),
        "p": ("hDensity", 1.0e6),
        "mun": ("eMobility", 1.0e-4),
        "mup": ("hMobility", 1.0e-4),
    }
    return mapped_state(fields, mapping, node_count)


def read_scalar_csv(path: Path) -> list[float]:
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        value_col = next((name for name in reader.fieldnames or [] if name != "node_id"), None)
        if value_col is None:
            raise RuntimeError(f"{path} has no value column")
        for row in reader:
            values[int(row["node_id"])] = float(row[value_col])
    return [values[index] for index in range(max(values) + 1)]


def mapped_state(raw: dict[str, list[float]], mapping: dict[str, tuple[str, float]], node_count: int) -> dict[str, list[float]]:
    state: dict[str, list[float]] = {}
    for output, (source, scale) in mapping.items():
        if source not in raw:
            raise RuntimeError(f"missing state field {source}")
        values = [value * scale for value in raw[source]]
        if len(values) != node_count:
            raise RuntimeError(f"{source} has {len(values)} values, expected {node_count}")
        state[output] = values
    return state


def edge_row(
    bias: float,
    source: str,
    edge: dict[str, Any],
    nodes: dict[int, dict[str, float]],
    contact_by_node: dict[int, set[str]],
    state: dict[str, list[float]],
    ni_model: list[float],
    vt: float,
) -> dict[str, Any]:
    i = int(edge["node0"])
    j = int(edge["node1"])
    h = float(edge["length_m"])
    couple = float(edge["couple_m"])
    coef_n = avg(state["mun"], i, j) * vt / h if h > 0.0 else 0.0
    coef_p = avg(state["mup"], i, j) * vt / h if h > 0.0 else 0.0
    dpsi = state["psi"][j] - state["psi"][i]
    n_flux_density = sgdiag.sg_electron_flux(state["n"][i], state["n"][j], dpsi, vt, coef_n)
    p_flux_density = sgdiag.sg_hole_flux(state["p"][i], state["p"][j], dpsi, vt, coef_p)
    ni_e = [inferred_ni(state["n"][node], state["psi"][node], state["phin"][node], vt, "electron") for node in (i, j)]
    ni_h = [inferred_ni(state["p"][node], state["psi"][node], state["phip"][node], vt, "hole") for node in (i, j)]
    n_flux_qf_model = sgdiag.sg_electron_flux_qf_variable_ni(
        ni_model[i], ni_model[j], state["psi"][i], state["psi"][j], state["phin"][i], state["phin"][j], vt, coef_n)
    n_flux_qf_inferred = sgdiag.sg_electron_flux_qf_variable_ni(
        ni_e[0], ni_e[1], state["psi"][i], state["psi"][j], state["phin"][i], state["phin"][j], vt, coef_n)
    p_flux_qf_model = sgdiag.sg_hole_flux_qf_variable_ni(
        ni_model[i], ni_model[j], state["psi"][i], state["psi"][j], state["phip"][i], state["phip"][j], vt, coef_p)
    p_flux_qf_inferred = sgdiag.sg_hole_flux_qf_variable_ni(
        ni_h[0], ni_h[1], state["psi"][i], state["psi"][j], state["phip"][i], state["phip"][j], vt, coef_p)
    efield = sgdiag.field_between(state["psi"], i, j, h)
    nqf_field = sgdiag.field_between(state["phin"], i, j, h)
    pqf_field = sgdiag.field_between(state["phip"], i, j, h)
    alpha_n_qf = sgdiag.van_overstraeten_alpha(nqf_field, 7.03e7, 7.03e7, 1.231e8, 1.231e8)
    alpha_n_ef = sgdiag.van_overstraeten_alpha(efield, 7.03e7, 7.03e7, 1.231e8, 1.231e8)
    alpha_p_qf = sgdiag.van_overstraeten_alpha(pqf_field, 1.582e8, 6.71e7, 2.036e8, 1.693e8)
    alpha_p_ef = sgdiag.van_overstraeten_alpha(efield, 1.582e8, 6.71e7, 2.036e8, 1.693e8)
    edge_area = 0.5 * h * couple
    edge_class, _ = sgdiag.classify_edge(edge, sgdiag.contact_names_for_edge(edge, contact_by_node), None)
    n0 = nodes[i]
    n1 = nodes[j]
    return {
        "bias_V": bias,
        "edge_id": int(edge["edge_id"]),
        "source": source,
        "node0": i,
        "node1": j,
        "x0_um": n0["x_um"],
        "y0_um": n0["y_um"],
        "x1_um": n1["x_um"],
        "y1_um": n1["y_um"],
        "length_m": h,
        "couple_m": couple,
        "edge_area_m2": edge_area,
        "edge_class": edge_class,
        "electric_field_V_m": efield,
        "electron_qf_field_V_m": nqf_field,
        "hole_qf_field_V_m": pqf_field,
        "electron_mobility_m2_V_s": avg(state["mun"], i, j),
        "hole_mobility_m2_V_s": avg(state["mup"], i, j),
        "electron_density_0_m3": state["n"][i],
        "electron_density_1_m3": state["n"][j],
        "hole_density_0_m3": state["p"][i],
        "hole_density_1_m3": state["p"][j],
        "psi0_V": state["psi"][i],
        "psi1_V": state["psi"][j],
        "phin0_V": state["phin"][i],
        "phin1_V": state["phin"][j],
        "phip0_V": state["phip"][i],
        "phip1_V": state["phip"][j],
        "ni_model_0_m3": ni_model[i],
        "ni_model_1_m3": ni_model[j],
        "ni_electron_inferred_0_m3": ni_e[0],
        "ni_electron_inferred_1_m3": ni_e[1],
        "ni_hole_inferred_0_m3": ni_h[0],
        "ni_hole_inferred_1_m3": ni_h[1],
        "electron_flux_density": n_flux_density,
        "electron_flux_density_abs": abs(n_flux_density),
        "electron_flux_qf_model": n_flux_qf_model,
        "electron_flux_qf_model_abs": abs(n_flux_qf_model),
        "electron_flux_qf_inferred": n_flux_qf_inferred,
        "electron_flux_qf_inferred_abs": abs(n_flux_qf_inferred),
        "electron_qf_model_over_density_abs": ratio_abs(n_flux_qf_model, n_flux_density),
        "electron_qf_inferred_over_density_abs": ratio_abs(n_flux_qf_inferred, n_flux_density),
        "hole_flux_density": p_flux_density,
        "hole_flux_density_abs": abs(p_flux_density),
        "hole_flux_qf_model": p_flux_qf_model,
        "hole_flux_qf_model_abs": abs(p_flux_qf_model),
        "hole_flux_qf_inferred": p_flux_qf_inferred,
        "hole_flux_qf_inferred_abs": abs(p_flux_qf_inferred),
        "electron_alpha_qf_m_inv": alpha_n_qf,
        "electron_alpha_efield_m_inv": alpha_n_ef,
        "hole_alpha_qf_m_inv": alpha_p_qf,
        "hole_alpha_efield_m_inv": alpha_p_ef,
        "electron_source_qf_model_qfdrive": alpha_n_qf * abs(n_flux_qf_model) * edge_area,
        "electron_source_qf_model_efielddrive": alpha_n_ef * abs(n_flux_qf_model) * edge_area,
        "electron_source_density_qfdrive": alpha_n_qf * abs(n_flux_density) * edge_area,
        "electron_source_density_efielddrive": alpha_n_ef * abs(n_flux_density) * edge_area,
    }


def avg(values: list[float], i: int, j: int) -> float:
    return 0.5 * (values[i] + values[j])


def inferred_ni(density: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if density <= 0.0:
        return 0.0
    if carrier == "electron":
        return density / sgdiag.limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return density / sgdiag.limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def ratio_abs(candidate: float, reference: float) -> float | None:
    if reference == 0.0:
        return None
    return abs(candidate) / abs(reference)


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


if __name__ == "__main__":
    raise SystemExit(main())
