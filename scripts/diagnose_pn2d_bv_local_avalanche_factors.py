#!/usr/bin/env python3
"""Compare local PN2D BV avalanche factors around a selected Vela SG edge."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


Q = 1.602176634e-19


CSV_FIELDS = [
    "bias_V",
    "edge_id",
    "vela_node0",
    "vela_node1",
    "edge_mid_x_um",
    "edge_mid_y_um",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "vela_source_density_m3_s",
    "sentaurus_generation_m3_s",
    "log10_vela_over_sentaurus_generation",
    "vela_electric_field_V_m",
    "sentaurus_electric_field_V_m",
    "vela_electron_qf_field_V_m",
    "vela_hole_qf_field_V_m",
    "vela_electron_alpha_m_inv",
    "vela_hole_alpha_m_inv",
    "sentaurus_alpha_from_electric_field_electron_m_inv",
    "sentaurus_alpha_from_electric_field_hole_m_inv",
    "sentaurus_weighted_alpha_m_inv",
    "vela_electron_flux_abs_m2_s",
    "vela_hole_flux_abs_m2_s",
    "sentaurus_electron_flux_abs_m2_s",
    "sentaurus_hole_flux_abs_m2_s",
    "vela_electron_density_cm3",
    "sentaurus_electron_density_cm3",
    "vela_hole_density_cm3",
    "sentaurus_hole_density_cm3",
    "vela_electron_mobility_cm2_V_s",
    "sentaurus_electron_mobility_cm2_V_s",
    "vela_hole_mobility_cm2_V_s",
    "sentaurus_hole_mobility_cm2_V_s",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2,-20")
    parser.add_argument("--edge-id", type=int, default=2886)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


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
        if (candidate / "nodes.csv").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def load_sentaurus_nodes(root: Path) -> list[dict[str, float]]:
    rows = []
    with (root / "nodes.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append({
                "id": int(row["id"]),
                "x_um": float(row["x_um"]),
                "y_um": float(row["y_um"]),
            })
    return sorted(rows, key=lambda item: int(item["id"]))


def load_field(root: Path, name: str) -> dict[int, float]:
    path = root / "fields" / f"{name}_region0.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = [column for column in reader.fieldnames or [] if column != "node_id"]
        for row in reader:
            comps = [float(row[column]) for column in columns if row.get(column, "") != ""]
            values[int(row["node_id"])] = comps[0] if len(comps) == 1 else math.sqrt(
                sum(value * value for value in comps))
    return values


def nearest_node(nodes: list[dict[str, float]], x_um: float, y_um: float) -> tuple[dict[str, float], float]:
    best = min(
        nodes,
        key=lambda item: (item["x_um"] - x_um) ** 2 + (item["y_um"] - y_um) ** 2,
    )
    distance = math.hypot(best["x_um"] - x_um, best["y_um"] - y_um)
    return best, distance


def log10_ratio(candidate: float, reference: float) -> float | None:
    if candidate == 0.0 or reference == 0.0:
        return None
    return math.log10(abs(candidate) / abs(reference))


def sentaurus_local_row(root: Path, x_um: float, y_um: float) -> dict[str, float]:
    nodes = load_sentaurus_nodes(root)
    node, distance_um = nearest_node(nodes, x_um, y_um)
    node_id = int(node["id"])
    fields = {
        name: load_field(root, name)
        for name in [
            "ElectricField",
            "ImpactIonization",
            "eDensity",
            "hDensity",
            "eMobility",
            "hMobility",
            "eCurrentDensity",
            "hCurrentDensity",
        ]
    }
    electric_field_m = fields["ElectricField"][node_id] * 100.0
    generation_m3_s = fields["ImpactIonization"][node_id] * 1.0e6
    electron_flux = abs(fields["eCurrentDensity"][node_id]) * 1.0e4 / Q
    hole_flux = abs(fields["hCurrentDensity"][node_id]) * 1.0e4 / Q
    weighted_alpha = (
        generation_m3_s / (electron_flux + hole_flux)
        if electron_flux + hole_flux > 0.0 else 0.0
    )
    return {
        "sentaurus_node_id": node_id,
        "sentaurus_distance_um": distance_um,
        "sentaurus_generation_m3_s": generation_m3_s,
        "sentaurus_electric_field_V_m": electric_field_m,
        "sentaurus_alpha_from_electric_field_electron_m_inv": sgdiag.van_overstraeten_alpha(
            electric_field_m, 7.03e7, 7.03e7, 1.231e8, 1.231e8),
        "sentaurus_alpha_from_electric_field_hole_m_inv": sgdiag.van_overstraeten_alpha(
            electric_field_m, 1.582e8, 6.71e7, 2.036e8, 1.693e8),
        "sentaurus_weighted_alpha_m_inv": weighted_alpha,
        "sentaurus_electron_flux_abs_m2_s": electron_flux,
        "sentaurus_hole_flux_abs_m2_s": hole_flux,
        "sentaurus_electron_density_cm3": fields["eDensity"][node_id],
        "sentaurus_hole_density_cm3": fields["hDensity"][node_id],
        "sentaurus_electron_mobility_cm2_V_s": fields["eMobility"][node_id],
        "sentaurus_hole_mobility_cm2_V_s": fields["hMobility"][node_id],
    }


def vela_edge_row(mesh: Path, doping_csv: Path, vtk: Path, edge_id: int, temperature_k: float) -> dict[str, Any]:
    nodes, triangles, contacts = sgdiag.read_mesh(mesh)
    edges = sgdiag.build_edges(nodes, triangles)
    scalars = sgdiag.parse_vtk_scalars(vtk)
    vt = sgdiag.K_B_OVER_Q * temperature_k
    ni = sgdiag.effective_ni_from_doping(
        doping_csv, len(nodes), 1.0e16, vt, "old_slotboom")
    rows, _ = sgdiag.source_rows(
        nodes, edges, contacts, scalars, vt, "quasi_fermi_gradient",
        "quasi_fermi_variable_ni", ni)
    row = next(item for item in rows if int(item["edge_id"]) == edge_id)
    n0 = int(row["node0"])
    n1 = int(row["node1"])
    edge_area = float(row["edge_area_m2"])
    x_mid = 0.5 * (float(row["x0_um"]) + float(row["x1_um"]))
    y_mid = 0.5 * (float(row["y0_um"]) + float(row["y1_um"]))
    return {
        "edge_id": edge_id,
        "vela_node0": n0,
        "vela_node1": n1,
        "edge_mid_x_um": x_mid,
        "edge_mid_y_um": y_mid,
        "vela_source_density_m3_s": (
            float(row["source_integral"]) / edge_area if edge_area > 0.0 else 0.0),
        "vela_electric_field_V_m": row["electric_field_V_m"],
        "vela_electron_qf_field_V_m": row["electron_field_V_m"],
        "vela_hole_qf_field_V_m": row["hole_field_V_m"],
        "vela_electron_alpha_m_inv": row["electron_alpha_m_inv"],
        "vela_hole_alpha_m_inv": row["hole_alpha_m_inv"],
        "vela_electron_flux_abs_m2_s": row["electron_flux_abs"],
        "vela_hole_flux_abs_m2_s": row["hole_flux_abs"],
        "vela_electron_density_cm3": 0.5 * (
            scalars["Electrons"][n0] + scalars["Electrons"][n1]) * 1.0e-6,
        "vela_hole_density_cm3": 0.5 * (
            scalars["Holes"][n0] + scalars["Holes"][n1]) * 1.0e-6,
        "vela_electron_mobility_cm2_V_s": row["electron_mobility_m2_V_s"] * 1.0e4,
        "vela_hole_mobility_cm2_V_s": row["hole_mobility_m2_V_s"] * 1.0e4,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
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


def main() -> int:
    args = parse_args()
    biases = parse_biases(args.biases)
    vtks = discover_vela_vtks(args.vela_vtk_root)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        vtk = vtks.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        vela = vela_edge_row(args.mesh, args.doping_csv, vtk, args.edge_id, args.temperature_k)
        sentaurus = sentaurus_local_row(
            sentaurus_dir(args.sentaurus_root, bias),
            float(vela["edge_mid_x_um"]),
            float(vela["edge_mid_y_um"]),
        )
        row = {"bias_V": bias, **vela, **sentaurus}
        row["log10_vela_over_sentaurus_generation"] = log10_ratio(
            float(row["vela_source_density_m3_s"]),
            float(row["sentaurus_generation_m3_s"]),
        )
        rows.append(row)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "local_avalanche_factors.csv", rows)
    (args.out_dir / "local_avalanche_factors.json").write_text(
        json.dumps(clean_json({"edge_id": args.edge_id, "rows": rows}), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
