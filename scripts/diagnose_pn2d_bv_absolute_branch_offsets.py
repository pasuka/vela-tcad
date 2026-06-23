#!/usr/bin/env python3
"""Summarize absolute quasi-Fermi branch offsets against Sentaurus exports."""

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


NODE_FIELDS = [
    "bias_V",
    "node_id",
    "node_class",
    "contact_names",
    "doping_class",
    "x_um",
    "y_um",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "impact_active",
    "vela_psi_V",
    "sentaurus_psi_V",
    "delta_psi_V",
    "vela_phin_V",
    "sentaurus_phin_V",
    "delta_phin_V",
    "vela_phip_V",
    "sentaurus_phip_V",
    "delta_phip_V",
    "vela_psi_minus_phin_V",
    "sentaurus_psi_minus_phin_V",
    "delta_psi_minus_phin_V",
    "vela_phip_minus_psi_V",
    "sentaurus_phip_minus_psi_V",
    "delta_phip_minus_psi_V",
    "vela_electron_density_cm3",
    "sentaurus_electron_density_cm3",
    "log10_vela_over_sentaurus_electron_density",
    "vela_hole_density_cm3",
    "sentaurus_hole_density_cm3",
    "log10_vela_over_sentaurus_hole_density",
    "sentaurus_impact_m3_s",
    "vela_avalanche_m3_s",
]

SUMMARY_FLOATS = [
    "delta_psi_V",
    "delta_phin_V",
    "delta_phip_V",
    "delta_psi_minus_phin_V",
    "delta_phip_minus_psi_V",
    "log10_vela_over_sentaurus_electron_density",
    "log10_vela_over_sentaurus_hole_density",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", required=True)
    parser.add_argument("--doping-csv", type=Path)
    parser.add_argument("--active-impact-relative-threshold", type=float, default=1.0e-6)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def bias_token(bias: float) -> str:
    return f"{bias:g}".replace("-", "m").replace("+", "").replace(".", "p")


def sentaurus_dir(root: Path, bias: float) -> Path:
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias_token(bias)}v",
        root / f"sentaurus_{bias_token(abs(bias))}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "nodes.csv").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    result: dict[float, Path] = {}
    mtimes: dict[float, float] = {}
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


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def load_sentaurus_nodes(root: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["id"]): {
            "id": int(row["id"]),
            "x_um": float(row["x_um"]),
            "y_um": float(row["y_um"]),
        }
        for row in read_csv_rows(root / "nodes.csv")
    }


def load_scalar(root: Path, name: str) -> dict[int, float]:
    path = root / "fields" / f"{name}_region0.csv"
    if not path.exists():
        return {}
    result: dict[int, float] = {}
    for row in read_csv_rows(path):
        if row.get("component0") not in (None, ""):
            result[int(row["node_id"])] = float(row["component0"])
    return result


def nearest_node(nodes: dict[int, dict[str, float]], x_um: float, y_um: float) -> tuple[int, float]:
    best = min(
        nodes.values(),
        key=lambda item: (item["x_um"] - x_um) ** 2 + (item["y_um"] - y_um) ** 2,
    )
    distance_um = math.hypot(best["x_um"] - x_um, best["y_um"] - y_um)
    return int(best["id"]), distance_um


def log10_ratio(candidate: float, reference: float) -> float | None:
    if candidate <= 0.0 or reference <= 0.0:
        return None
    return math.log10(candidate / reference)


def boundary_nodes(nodes: dict[int, dict[str, float]], triangles: list[dict[str, Any]]) -> set[int]:
    edges = sgdiag.build_edges(nodes, triangles)
    result: set[int] = set()
    for edge in edges:
        if len(edge["cell_ids"]) == 1:
            result.add(int(edge["node0"]))
            result.add(int(edge["node1"]))
    return result


def load_doping_classes(path: Path | None) -> dict[int, str]:
    if path is None or not path.exists():
        return {}
    result: dict[int, str] = {}
    for row in read_csv_rows(path):
        node_id = int(row["node_id"])
        donors = float(row.get("donors_cm3") or 0.0)
        acceptors = float(row.get("acceptors_cm3") or 0.0)
        if donors == 0.0 and acceptors == 0.0:
            donors = float(row.get("donors_m3") or 0.0) * 1.0e-6
            acceptors = float(row.get("acceptors_m3") or 0.0) * 1.0e-6
        net = donors - acceptors
        if abs(net) <= max(abs(donors), abs(acceptors), 1.0) * 1.0e-9:
            result[node_id] = "compensated"
        elif net > 0.0:
            result[node_id] = "n_type"
        else:
            result[node_id] = "p_type"
    return result


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def median(values: list[float]) -> float | None:
    clean = [value for value in values if math.isfinite(value)]
    return statistics.median(clean) if clean else None


def summarize_group(bias: float, group: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    item: dict[str, Any] = {
        "bias_V": bias,
        "group": group,
        "node_count": len(rows),
    }
    for key in SUMMARY_FLOATS:
        values = [float(row[key]) for row in rows if row.get(key) is not None]
        item[f"median_{key}"] = median(values)
    return item


def grouped_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_bias: dict[float, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_bias[bias_key(float(row["bias_V"]))].append(row)
    result: list[dict[str, Any]] = []
    for bias, items in sorted(by_bias.items()):
        groups: dict[str, list[dict[str, Any]]] = {
            "all": items,
            "contact": [row for row in items if row["node_class"] == "contact"],
            "noncontact": [row for row in items if row["node_class"] != "contact"],
            "impact_active": [row for row in items if row["impact_active"] is True],
            "impact_inactive": [row for row in items if row["impact_active"] is False],
        }
        for doping_class in sorted({str(row.get("doping_class", "unknown")) for row in items}):
            groups[f"doping_{doping_class}"] = [
                row for row in items if row.get("doping_class") == doping_class
            ]
        for group, group_rows in groups.items():
            if group_rows:
                result.append(summarize_group(bias, group, group_rows))
    return result


def append_bias_rows(
    bias: float,
    mesh: Path,
    sentaurus: Path,
    vtk: Path,
    doping_classes: dict[int, str],
    active_relative_threshold: float,
    rows: list[dict[str, Any]],
) -> None:
    nodes, triangles, contacts = sgdiag.read_mesh(mesh)
    boundary = boundary_nodes(nodes, triangles)
    scalars = sgdiag.parse_vtk_scalars(vtk)
    sent_nodes = load_sentaurus_nodes(sentaurus)
    sent_psi = load_scalar(sentaurus, "ElectrostaticPotential")
    sent_phin = load_scalar(sentaurus, "eQuasiFermiPotential")
    sent_phip = load_scalar(sentaurus, "hQuasiFermiPotential")
    sent_n = load_scalar(sentaurus, "eDensity")
    sent_p = load_scalar(sentaurus, "hDensity")
    sent_impact = load_scalar(sentaurus, "ImpactIonization")
    max_impact = max((abs(value) for value in sent_impact.values()), default=0.0)
    active_threshold = max_impact * active_relative_threshold

    for node_id in sorted(nodes):
        node = nodes[node_id]
        sent_id, distance_um = nearest_node(sent_nodes, node["x_um"], node["y_um"])
        contact_names = sorted(contacts.get(node_id, set()))
        if contact_names:
            node_class = "contact"
        elif node_id in boundary:
            node_class = "boundary"
        else:
            node_class = "interior"

        vela_psi = scalars["Potential"][node_id]
        vela_phin = scalars["ElectronQuasiFermi"][node_id]
        vela_phip = scalars["HoleQuasiFermi"][node_id]
        s_psi = sent_psi[sent_id]
        s_phin = sent_phin[sent_id]
        s_phip = sent_phip[sent_id]
        vela_psi_minus_phin = vela_psi - vela_phin
        sent_psi_minus_phin = s_psi - s_phin
        vela_phip_minus_psi = vela_phip - vela_psi
        sent_phip_minus_psi = s_phip - s_psi
        vela_n = scalars.get("Electrons", [0.0 for _ in nodes])[node_id] * 1.0e-6
        vela_p = scalars.get("Holes", [0.0 for _ in nodes])[node_id] * 1.0e-6
        s_n = sent_n.get(sent_id, 0.0)
        s_p = sent_p.get(sent_id, 0.0)
        s_impact = sent_impact.get(sent_id, 0.0) * 1.0e6
        v_impact = scalars.get("AvalancheGeneration", [0.0 for _ in nodes])[node_id]
        rows.append({
            "bias_V": bias,
            "node_id": node_id,
            "node_class": node_class,
            "contact_names": ";".join(contact_names),
            "doping_class": doping_classes.get(node_id, "unknown"),
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "sentaurus_node_id": sent_id,
            "sentaurus_distance_um": distance_um,
            "impact_active": abs(sent_impact.get(sent_id, 0.0)) >= active_threshold if max_impact > 0.0 else False,
            "vela_psi_V": vela_psi,
            "sentaurus_psi_V": s_psi,
            "delta_psi_V": vela_psi - s_psi,
            "vela_phin_V": vela_phin,
            "sentaurus_phin_V": s_phin,
            "delta_phin_V": vela_phin - s_phin,
            "vela_phip_V": vela_phip,
            "sentaurus_phip_V": s_phip,
            "delta_phip_V": vela_phip - s_phip,
            "vela_psi_minus_phin_V": vela_psi_minus_phin,
            "sentaurus_psi_minus_phin_V": sent_psi_minus_phin,
            "delta_psi_minus_phin_V": vela_psi_minus_phin - sent_psi_minus_phin,
            "vela_phip_minus_psi_V": vela_phip_minus_psi,
            "sentaurus_phip_minus_psi_V": sent_phip_minus_psi,
            "delta_phip_minus_psi_V": vela_phip_minus_psi - sent_phip_minus_psi,
            "vela_electron_density_cm3": vela_n,
            "sentaurus_electron_density_cm3": s_n,
            "log10_vela_over_sentaurus_electron_density": log10_ratio(vela_n, s_n),
            "vela_hole_density_cm3": vela_p,
            "sentaurus_hole_density_cm3": s_p,
            "log10_vela_over_sentaurus_hole_density": log10_ratio(vela_p, s_p),
            "sentaurus_impact_m3_s": s_impact,
            "vela_avalanche_m3_s": v_impact,
        })


def main() -> int:
    args = parse_args()
    vtk_by_bias = discover_vela_vtks(args.vela_vtk_root)
    doping_classes = load_doping_classes(args.doping_csv)
    rows: list[dict[str, Any]] = []
    for bias in parse_biases(args.biases):
        vtk = vtk_by_bias.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        append_bias_rows(
            bias,
            args.mesh,
            sentaurus_dir(args.sentaurus_root, bias),
            vtk,
            doping_classes,
            args.active_impact_relative_threshold,
            rows,
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "absolute_branch_offsets_nodes.csv", rows, NODE_FIELDS)
    summary = {
        "biases": parse_biases(args.biases),
        "node_rows": len(rows),
        "groups": grouped_summary(rows),
    }
    (args.out_dir / "absolute_branch_offsets_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())