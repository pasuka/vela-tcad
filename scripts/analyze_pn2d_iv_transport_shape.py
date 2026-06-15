#!/usr/bin/env python3
"""Compare pn2d IV transport-shape drivers between Sentaurus and Vela."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

ELEMENTARY_CHARGE_C = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-root",
        type=Path,
        default=Path("build/reference_tcad/pn2d_sentaurus2018"),
        help="pn2d Sentaurus2018 reference root under build/.",
    )
    parser.add_argument(
        "--biases",
        default="0.25,0.3,0.5,0.8,1.0",
        help="Comma-separated Vela bias points to inspect.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_transport_shape"),
        help="Directory for generated transport-shape CSV/JSON reports.",
    )
    parser.add_argument(
        "--sentaurus-multibias-root",
        type=Path,
        default=None,
        help="Optional root containing one exported Sentaurus TDR directory per bias.",
    )
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(part.strip()) for part in raw.split(",") if part.strip()]


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    fields: dict[str, list[float]] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS",
                    "VECTORS",
                    "FIELD",
                    "CELL_DATA",
                    "POINT_DATA",
                }:
                    break
                values.extend(float(v) for v in next_parts)
                i += 1
            fields[name] = values
            continue
        i += 1
    return fields


def abs_stats(values: list[float]) -> dict[str, float]:
    clean = sorted(abs(v) for v in values if math.isfinite(v))
    if not clean:
        return {"points": 0}

    def pct(p: float) -> float:
        idx = (len(clean) - 1) * p / 100.0
        lo = math.floor(idx)
        hi = math.ceil(idx)
        if lo == hi:
            return clean[lo]
        return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)

    return {
        "points": len(clean),
        "mean_abs": sum(clean) / len(clean),
        "median_abs": pct(50),
        "p95_abs": pct(95),
        "max_abs": clean[-1],
    }


def value_stats(values: list[float]) -> dict[str, float]:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return {"points": 0}

    def pct(p: float) -> float:
        idx = (len(clean) - 1) * p / 100.0
        lo = math.floor(idx)
        hi = math.ceil(idx)
        if lo == hi:
            return clean[lo]
        return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)

    return {
        "points": len(clean),
        "min": clean[0],
        "mean": sum(clean) / len(clean),
        "median": pct(50),
        "p95": pct(95),
        "max": clean[-1],
        "span": clean[-1] - clean[0],
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def percentile(sorted_values: list[float], p: float) -> float:
    idx = (len(sorted_values) - 1) * p / 100.0
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return sorted_values[lo]
    return sorted_values[lo] * (hi - idx) + sorted_values[hi] * (idx - lo)


def finite_value_stats(values: list[float]) -> dict[str, float]:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return {"points": 0}
    return {
        "points": len(clean),
        "min": clean[0],
        "mean": sum(clean) / len(clean),
        "median": percentile(clean, 50),
        "p05": percentile(clean, 5),
        "p95": percentile(clean, 95),
        "max": clean[-1],
    }


def read_numeric_component_csv(path: Path) -> list[float]:
    with path.open() as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            return []
        value_fields = [name for name in reader.fieldnames if name != "node_id"]
        if not value_fields:
            return []
        fields = value_fields
        values: list[float] = []
        for row in reader:
            comps = [float(row[name]) for name in fields if row.get(name, "") != ""]
            if len(comps) == 1:
                values.append(comps[0])
            elif comps:
                values.append(math.sqrt(sum(v * v for v in comps)))
        return values


def read_numeric_component_map(path: Path) -> dict[int, float]:
    with path.open() as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            return {}
        value_fields = [name for name in reader.fieldnames if name != "node_id"]
        values: dict[int, float] = {}
        for row in reader:
            comps = [float(row[name]) for name in value_fields if row.get(name, "") != ""]
            if not comps:
                continue
            values[int(row["node_id"])] = comps[0] if len(comps) == 1 else math.sqrt(
                sum(v * v for v in comps)
            )
        return values


def generate_sentaurus_field_stats(root: Path, out_dir: Path) -> list[dict[str, Any]]:
    fields_dir = root / "sim_fields" / "iv" / "fields"
    rows: list[dict[str, Any]] = []
    field_names = [
        "eCurrentDensity",
        "hCurrentDensity",
        "TotalCurrentDensity",
        "ElectricField",
        "srhRecombination",
    ]
    for name in field_names:
        path = fields_dir / f"{name}_region0.csv"
        if not path.exists():
            rows.append({"source_bias_V": 1.0, "field": name, "status": "missing"})
            continue
        stats = abs_stats(read_numeric_component_csv(path))
        rows.append(
            {
                "source_bias_V": 1.0,
                "field": name,
                "status": "ok",
                **stats,
                "source_file": str(path),
            }
        )
    fields = [
        "source_bias_V",
        "field",
        "status",
        "points",
        "mean_abs",
        "median_abs",
        "p95_abs",
        "max_abs",
        "source_file",
    ]
    write_csv(out_dir / "current_density_summary.csv", rows, fields)
    (out_dir / "current_density_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def discover_vtk_by_bias(root: Path) -> dict[float, Path]:
    fixed_dir = root / "reports" / "iv_state" / "fixed"
    pattern = re.compile(r"iv_1v_fixed_probe_\d+_(.+)V\.vtk$")
    by_bias: dict[float, Path] = {}
    for path in fixed_dir.glob("iv_1v_fixed_probe_*.vtk"):
        match = pattern.match(path.name)
        if not match:
            continue
        by_bias[float(match.group(1))] = path
    return by_bias


def nearest_bias_path(vtk_by_bias: dict[float, Path], bias: float) -> tuple[float, Path] | None:
    if not vtk_by_bias:
        return None
    nearest = min(vtk_by_bias, key=lambda b: abs(b - bias))
    if abs(nearest - bias) > 1e-8:
        return None
    return nearest, vtk_by_bias[nearest]


def generate_vela_proxy_stats(root: Path, out_dir: Path, biases: list[float]) -> list[dict[str, Any]]:
    vtk_by_bias = discover_vtk_by_bias(root)
    rows: list[dict[str, Any]] = []
    for requested_bias in biases:
        found = nearest_bias_path(vtk_by_bias, requested_bias)
        if found is None:
            rows.append({"bias_V": requested_bias, "quantity": "vtk", "status": "missing"})
            continue
        bias, path = found
        scalars = parse_vtk_scalars(path)
        quantity_map = {
            "Potential": ("Potential", 1.0, "V"),
            "ElectronQuasiFermi": ("ElectronQuasiFermi", 1.0, "V"),
            "HoleQuasiFermi": ("HoleQuasiFermi", 1.0, "V"),
            "Electrons": ("Electrons", 1.0e-6, "cm^-3"),
            "Holes": ("Holes", 1.0e-6, "cm^-3"),
            "NetDoping": ("NetDoping", 1.0e-6, "cm^-3"),
        }
        for field, (quantity, scale, unit) in quantity_map.items():
            values = scalars.get(field)
            if values is None:
                rows.append(
                    {
                        "bias_V": bias,
                        "quantity": quantity,
                        "status": "missing",
                        "vtk": str(path),
                    }
                )
                continue
            stats = value_stats([v * scale for v in values])
            rows.append(
                {
                    "bias_V": bias,
                    "quantity": quantity,
                    "status": "ok",
                    "unit": unit,
                    **stats,
                    "vtk": str(path),
                }
            )
    fields = [
        "bias_V",
        "quantity",
        "status",
        "unit",
        "points",
        "min",
        "mean",
        "median",
        "p95",
        "max",
        "span",
        "vtk",
    ]
    write_csv(out_dir / "vela_state_proxy_summary.csv", rows, fields)
    (out_dir / "vela_state_proxy_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def load_iv_ratios(root: Path) -> dict[float, float]:
    ratio_path = root / "reports" / "iv_residual_current" / "baseline" / "iv_ratio_by_bias.csv"
    ratios: dict[float, float] = {}
    if not ratio_path.exists():
        return ratios
    for row in csv.DictReader(ratio_path.open()):
        ratios[round(float(row["bias_V"]), 12)] = float(row["ratio_vela_to_sentaurus"])
    return ratios


def generate_contact_edge_concentration(root: Path, out_dir: Path) -> list[dict[str, Any]]:
    fixed_dir = root / "reports" / "iv_state" / "fixed"
    edge_path = fixed_dir / "iv_1v_fixed_probe_contact_edges.csv"
    terminal_path = fixed_dir / "iv_1v_fixed_probe_terminal_balance.csv"
    grouped: dict[tuple[float, str], list[float]] = defaultdict(list)
    if edge_path.exists():
        for row in csv.DictReader(edge_path.open()):
            grouped[(float(row["bias_V"]), row["current_contact"])].append(
                float(row["current_total_A_per_um"])
            )

    terminals: dict[tuple[float, str], float] = {}
    if terminal_path.exists():
        for row in csv.DictReader(terminal_path.open()):
            terminals[(float(row["bias_V"]), row["contact"])] = float(row["current_total_A_per_um"])

    rows: list[dict[str, Any]] = []
    for key, values in sorted(grouped.items()):
        bias, contact = key
        abs_values = sorted((abs(v) for v in values), reverse=True)
        abs_sum = sum(abs_values)
        edge_sum = sum(values)
        terminal = terminals.get(key)
        sign_changes = sum(
            1
            for a, b in zip(values, values[1:])
            if abs(a) > 0.0 and abs(b) > 0.0 and math.copysign(1.0, a) != math.copysign(1.0, b)
        )
        terminal_delta = None if terminal is None else edge_sum - terminal
        rows.append(
            {
                "bias_V": bias,
                "contact": contact,
                "edge_count": len(values),
                "sum_current_A_per_um": edge_sum,
                "terminal_total_A_per_um": "" if terminal is None else terminal,
                "edge_minus_terminal_A_per_um": "" if terminal_delta is None else terminal_delta,
                "abs_sum_current_A_per_um": abs_sum,
                "max_abs_edge_current_A_per_um": abs_values[0] if abs_values else 0.0,
                "top3_abs_fraction": sum(abs_values[:3]) / abs_sum if abs_sum else 0.0,
                "edge_sign_changes": sign_changes,
                "terminal_gate_pass": terminal_delta is not None and abs(terminal_delta) < 1.0e-18,
            }
        )
    fields = [
        "bias_V",
        "contact",
        "edge_count",
        "sum_current_A_per_um",
        "terminal_total_A_per_um",
        "edge_minus_terminal_A_per_um",
        "abs_sum_current_A_per_um",
        "max_abs_edge_current_A_per_um",
        "top3_abs_fraction",
        "edge_sign_changes",
        "terminal_gate_pass",
    ]
    write_csv(out_dir / "contact_edge_concentration.csv", rows, fields)
    (out_dir / "contact_edge_concentration.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def load_mesh(root: Path) -> dict[str, Any] | None:
    path = root / "vela" / "mesh.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def mesh_node_map(mesh: dict[str, Any]) -> dict[int, tuple[float, float]]:
    return {int(node["id"]): (float(node["x"]) * 1.0e-6, float(node["y"]) * 1.0e-6)
            for node in mesh.get("nodes", [])}


def mesh_unique_edges(mesh: dict[str, Any]) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for triangle in mesh.get("triangles", []):
        nodes = [int(node) for node in triangle.get("node_ids", [])]
        if len(nodes) != 3:
            continue
        for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
            edges.add((min(a, b), max(a, b)))
    return sorted(edges)


def edge_length_m(nodes: dict[int, tuple[float, float]], a: int, b: int) -> float:
    ax, ay = nodes[a]
    bx, by = nodes[b]
    return math.hypot(ax - bx, ay - by)


def contact_adjacent_edges(mesh: dict[str, Any]) -> list[dict[str, Any]]:
    nodes = mesh_node_map(mesh)
    edges = mesh_unique_edges(mesh)
    rows: list[dict[str, Any]] = []
    for contact in mesh.get("contacts", []):
        contact_nodes = {int(node_id) for node_id in contact.get("node_ids", [])}
        for a, b in edges:
            a_on = a in contact_nodes
            b_on = b in contact_nodes
            if a_on == b_on:
                continue
            contact_node, interior_node = (a, b) if a_on else (b, a)
            length = edge_length_m(nodes, contact_node, interior_node)
            if length <= 0.0:
                continue
            rows.append(
                {
                    "contact": contact["name"],
                    "contact_node": contact_node,
                    "interior_node": interior_node,
                    "edge_length_m": length,
                }
            )
    return rows


def metric_summary_rows(edge_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    for row in edge_rows:
        key_prefix = (row["source"], row["contact"], row["carrier"])
        for metric in ["qf_drop_V", "qf_gradient_V_per_m", "density_m3", "q_n_grad", "effective_mu_m2_V_s"]:
            value = row.get(metric, "")
            if value == "" or value is None:
                continue
            grouped[(*key_prefix, metric)].append(float(value))

    summaries: list[dict[str, Any]] = []
    for (source, contact, carrier, metric), values in sorted(grouped.items()):
        summaries.append(
            {
                "source": source,
                "contact": contact,
                "carrier": carrier,
                "metric": metric,
                **finite_value_stats(values),
            }
        )
    return summaries


def compare_metric_summary(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key: dict[tuple[str, str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in summary_rows:
        by_key[(row["contact"], row["carrier"], row["metric"])][row["source"]] = row

    rows: list[dict[str, Any]] = []
    for (contact, carrier, metric), sources in sorted(by_key.items()):
        sentaurus = sources.get("Sentaurus")
        vela = sources.get("Vela")
        if sentaurus is None or vela is None:
            continue
        sentaurus_mean = float(sentaurus["mean"])
        vela_mean = float(vela["mean"])
        rows.append(
            {
                "contact": contact,
                "carrier": carrier,
                "metric": metric,
                "points_sentaurus": sentaurus["points"],
                "mean_sentaurus": sentaurus_mean,
                "median_sentaurus": sentaurus["median"],
                "points_vela": vela["points"],
                "mean_vela": vela_mean,
                "median_vela": vela["median"],
                "vela_to_sentaurus_mean": (
                    "" if sentaurus_mean == 0.0 else vela_mean / sentaurus_mean
                ),
            }
        )
    return rows


def read_contact_voltage_bias(fields_dir: Path) -> float | None:
    voltages: list[float] = []
    for path in fields_dir.glob("ContactExternalVoltage_region*.csv"):
        for row in csv.DictReader(path.open()):
            if row.get("component0", "") != "":
                voltages.append(float(row["component0"]))
    if not voltages:
        return None
    return max(voltages, key=abs)


def discover_sentaurus_multibias_exports(multibias_root: Path) -> list[tuple[float, Path]]:
    exports: list[tuple[float, Path]] = []
    if not multibias_root.exists():
        return exports
    for path in sorted(p for p in multibias_root.iterdir() if p.is_dir()):
        fields_dir = path / "fields"
        if not fields_dir.exists():
            continue
        bias = read_contact_voltage_bias(fields_dir)
        if bias is None:
            continue
        exports.append((bias, path))
    return sorted(exports, key=lambda item: item[0])


def infer_ni_eff_stats_from_sentaurus_multibias(
    multibias_root: Path | None,
    out_dir: Path,
    temperature_K: float = 300.0,
) -> list[dict[str, Any]]:
    if multibias_root is None:
        return []
    exports = discover_sentaurus_multibias_exports(multibias_root)
    if not exports:
        return []

    vt = 8.617333262145e-5 * temperature_K
    rows: list[dict[str, Any]] = []
    required = [
        "ElectrostaticPotential",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eDensity",
        "hDensity",
    ]
    for bias, export_dir in exports:
        fields_dir = export_dir / "fields"
        if any(not (fields_dir / f"{name}_region0.csv").exists() for name in required):
            continue
        psi = read_numeric_component_map(fields_dir / "ElectrostaticPotential_region0.csv")
        phin = read_numeric_component_map(fields_dir / "eQuasiFermiPotential_region0.csv")
        phip = read_numeric_component_map(fields_dir / "hQuasiFermiPotential_region0.csv")
        n = read_numeric_component_map(fields_dir / "eDensity_region0.csv")
        p = read_numeric_component_map(fields_dir / "hDensity_region0.csv")

        inferred = {
            "electron_inferred": [],
            "hole_inferred": [],
        }
        for node_id, psi_value in psi.items():
            if node_id in phin and node_id in n:
                exponent = max(min((psi_value - phin[node_id]) / vt, 500.0), -500.0)
                ni_eff = n[node_id] / math.exp(exponent)
                if ni_eff > 0.0 and math.isfinite(ni_eff):
                    inferred["electron_inferred"].append(ni_eff)
            if node_id in phip and node_id in p:
                exponent = max(min((phip[node_id] - psi_value) / vt, 500.0), -500.0)
                ni_eff = p[node_id] / math.exp(exponent)
                if ni_eff > 0.0 and math.isfinite(ni_eff):
                    inferred["hole_inferred"].append(ni_eff)

        for source, values in inferred.items():
            stats = finite_value_stats(values)
            if not stats.get("points"):
                continue
            rows.append(
                {
                    "bias_V": bias,
                    "source": source,
                    "points": stats["points"],
                    "mean_ni_eff_cm3": stats["mean"],
                    "median_ni_eff_cm3": stats["median"],
                    "p05_ni_eff_cm3": stats["p05"],
                    "p95_ni_eff_cm3": stats["p95"],
                    "min_ni_eff_cm3": stats["min"],
                    "max_ni_eff_cm3": stats["max"],
                    "source_export": str(export_dir),
                }
            )

    if rows:
        fields = [
            "bias_V",
            "source",
            "points",
            "mean_ni_eff_cm3",
            "median_ni_eff_cm3",
            "p05_ni_eff_cm3",
            "p95_ni_eff_cm3",
            "min_ni_eff_cm3",
            "max_ni_eff_cm3",
            "source_export",
        ]
        write_csv(out_dir / "sentaurus_inferred_ni_eff_multibias.csv", rows, fields)
        (out_dir / "sentaurus_inferred_ni_eff_multibias.json").write_text(
            json.dumps(rows, indent=2) + "\n"
        )
    return rows


def infer_contact_edge_ni_eff_from_sentaurus_multibias(
    root: Path,
    multibias_root: Path | None,
    out_dir: Path,
    temperature_K: float = 300.0,
) -> list[dict[str, Any]]:
    mesh = load_mesh(root)
    if mesh is None or multibias_root is None:
        return []

    exports = discover_sentaurus_multibias_exports(multibias_root)
    if not exports:
        return []

    vt = 8.617333262145e-5 * temperature_K
    edge_defs = contact_adjacent_edges(mesh)
    required = [
        "ElectrostaticPotential",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eDensity",
        "hDensity",
    ]
    grouped: dict[tuple[float, str, str, str], list[float]] = defaultdict(list)
    for bias, export_dir in exports:
        fields_dir = export_dir / "fields"
        if any(not (fields_dir / f"{name}_region0.csv").exists() for name in required):
            continue
        psi = read_numeric_component_map(fields_dir / "ElectrostaticPotential_region0.csv")
        phin = read_numeric_component_map(fields_dir / "eQuasiFermiPotential_region0.csv")
        phip = read_numeric_component_map(fields_dir / "hQuasiFermiPotential_region0.csv")
        n = read_numeric_component_map(fields_dir / "eDensity_region0.csv")
        p = read_numeric_component_map(fields_dir / "hDensity_region0.csv")

        for edge in edge_defs:
            contact = edge["contact"]
            carrier = "electron" if contact == "Cathode" else "hole"
            for location, node_id in [
                ("contact", int(edge["contact_node"])),
                ("interior", int(edge["interior_node"])),
            ]:
                if node_id not in psi:
                    continue
                if carrier == "electron" and node_id in phin and node_id in n:
                    exponent = max(min((psi[node_id] - phin[node_id]) / vt, 500.0), -500.0)
                    ni_eff = n[node_id] / math.exp(exponent)
                elif carrier == "hole" and node_id in phip and node_id in p:
                    exponent = max(min((phip[node_id] - psi[node_id]) / vt, 500.0), -500.0)
                    ni_eff = p[node_id] / math.exp(exponent)
                else:
                    continue
                if ni_eff > 0.0 and math.isfinite(ni_eff):
                    grouped[(bias, contact, location, carrier)].append(ni_eff)

    rows: list[dict[str, Any]] = []
    for (bias, contact, location, carrier), values in sorted(grouped.items()):
        stats = finite_value_stats(values)
        if not stats.get("points"):
            continue
        rows.append(
            {
                "bias_V": bias,
                "contact": contact,
                "location": location,
                "carrier": carrier,
                "points": stats["points"],
                "mean_ni_eff_cm3": stats["mean"],
                "median_ni_eff_cm3": stats["median"],
                "p05_ni_eff_cm3": stats["p05"],
                "p95_ni_eff_cm3": stats["p95"],
                "min_ni_eff_cm3": stats["min"],
                "max_ni_eff_cm3": stats["max"],
            }
        )

    if rows:
        fields = [
            "bias_V",
            "contact",
            "location",
            "carrier",
            "points",
            "mean_ni_eff_cm3",
            "median_ni_eff_cm3",
            "p05_ni_eff_cm3",
            "p95_ni_eff_cm3",
            "min_ni_eff_cm3",
            "max_ni_eff_cm3",
        ]
        write_csv(out_dir / "sentaurus_contact_edge_inferred_ni_eff_multibias.csv", rows, fields)
        (out_dir / "sentaurus_contact_edge_inferred_ni_eff_multibias.json").write_text(
            json.dumps(rows, indent=2) + "\n"
        )
    return rows


def load_vela_edge_mobility(root: Path, bias: float) -> dict[tuple[str, int, int], dict[str, float]]:
    path = root / "reports" / "iv_state" / "fixed" / "iv_1v_fixed_probe_contact_edges.csv"
    rows: dict[tuple[str, int, int], dict[str, float]] = {}
    if not path.exists():
        return rows
    for row in csv.DictReader(path.open()):
        if abs(float(row["bias_V"]) - bias) > 1.0e-8:
            continue
        a = int(row["node0"])
        b = int(row["node1"])
        lo, hi = sorted((a, b))
        rows[(row["current_contact"], lo, hi)] = {
            "electron": float(row["mun"]),
            "hole": float(row["mup"]),
        }
    return rows


def sample_node_value(container: dict[int, float] | list[float], node: int) -> float:
    return container[node] if isinstance(container, list) else container[node]


def metric_summary_rows_by_bias(edge_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[float, str, str, str, str], list[float]] = defaultdict(list)
    metrics = [
        "qf_drop_V",
        "qf_gradient_V_per_m",
        "density_m3",
        "q_n_grad",
        "effective_mu_m2_V_s",
        "mobility_m2_V_s",
    ]
    for row in edge_rows:
        key_prefix = (float(row["bias_V"]), row["source"], row["contact"], row["carrier"])
        for metric in metrics:
            value = row.get(metric, "")
            if value == "" or value is None:
                continue
            grouped[(*key_prefix, metric)].append(float(value))

    summaries: list[dict[str, Any]] = []
    for (bias, source, contact, carrier, metric), values in sorted(grouped.items()):
        summaries.append(
            {
                "bias_V": bias,
                "source": source,
                "contact": contact,
                "carrier": carrier,
                "metric": metric,
                **finite_value_stats(values),
            }
        )
    return summaries


def compare_metric_summary_by_bias(summary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key: dict[tuple[float, str, str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in summary_rows:
        by_key[(float(row["bias_V"]), row["contact"], row["carrier"], row["metric"])][
            row["source"]
        ] = row

    rows: list[dict[str, Any]] = []
    for (bias, contact, carrier, metric), sources in sorted(by_key.items()):
        sentaurus = sources.get("Sentaurus")
        vela = sources.get("Vela")
        if sentaurus is None or vela is None:
            continue
        sentaurus_mean = float(sentaurus["mean"])
        vela_mean = float(vela["mean"])
        rows.append(
            {
                "bias_V": bias,
                "contact": contact,
                "carrier": carrier,
                "metric": metric,
                "points_sentaurus": sentaurus["points"],
                "mean_sentaurus": sentaurus_mean,
                "median_sentaurus": sentaurus["median"],
                "points_vela": vela["points"],
                "mean_vela": vela_mean,
                "median_vela": vela["median"],
                "vela_to_sentaurus_mean": (
                    "" if sentaurus_mean == 0.0 else vela_mean / sentaurus_mean
                ),
            }
        )
    return rows


def generate_multibias_contact_edge_transport_proxy(
    root: Path,
    multibias_root: Path | None,
    out_dir: Path,
) -> list[dict[str, Any]]:
    if multibias_root is None:
        candidate = root / "sentaurus_multibias_exports"
        multibias_root = candidate if candidate.exists() else None
    if multibias_root is None:
        return []

    mesh = load_mesh(root)
    vtk_by_bias = discover_vtk_by_bias(root)
    exports = discover_sentaurus_multibias_exports(multibias_root)
    required_sentaurus = [
        "eCurrentDensity",
        "hCurrentDensity",
        "eDensity",
        "hDensity",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eMobility",
        "hMobility",
    ]
    edge_rows: list[dict[str, Any]] = []
    if mesh is None:
        exports = []

    edge_defs = contact_adjacent_edges(mesh) if mesh is not None else []
    for bias, export_dir in exports:
        fields_dir = export_dir / "fields"
        if any(not (fields_dir / f"{name}_region0.csv").exists() for name in required_sentaurus):
            continue
        vtk = nearest_bias_path(vtk_by_bias, bias)
        if vtk is None:
            continue
        _, vtk_path = vtk
        vtk_scalars = parse_vtk_scalars(vtk_path)
        vela_edge_mobility = load_vela_edge_mobility(root, bias)
        sentaurus = {
            name: read_numeric_component_map(fields_dir / f"{name}_region0.csv")
            for name in required_sentaurus
        }
        sources = [
            {
                "name": "Sentaurus",
                "electron_qf": sentaurus["eQuasiFermiPotential"],
                "hole_qf": sentaurus["hQuasiFermiPotential"],
                "electron_density": sentaurus["eDensity"],
                "hole_density": sentaurus["hDensity"],
                "electron_current_density": sentaurus["eCurrentDensity"],
                "hole_current_density": sentaurus["hCurrentDensity"],
                "electron_mobility": sentaurus["eMobility"],
                "hole_mobility": sentaurus["hMobility"],
                "density_scale": 1.0e6,
                "current_density_scale": 1.0e4,
                "mobility_scale": 1.0e-4,
            },
            {
                "name": "Vela",
                "electron_qf": vtk_scalars.get("ElectronQuasiFermi", []),
                "hole_qf": vtk_scalars.get("HoleQuasiFermi", []),
                "electron_density": vtk_scalars.get("Electrons", []),
                "hole_density": vtk_scalars.get("Holes", []),
                "electron_current_density": None,
                "hole_current_density": None,
                "electron_mobility": None,
                "hole_mobility": None,
                "density_scale": 1.0,
                "current_density_scale": 1.0,
                "mobility_scale": 1.0,
            },
        ]

        for source in sources:
            for edge in edge_defs:
                contact_node = int(edge["contact_node"])
                interior_node = int(edge["interior_node"])
                length = float(edge["edge_length_m"])
                lo, hi = sorted((contact_node, interior_node))
                mobility_by_carrier = vela_edge_mobility.get((edge["contact"], lo, hi), {})
                for carrier in ["electron", "hole"]:
                    qf = source[f"{carrier}_qf"]
                    density = source[f"{carrier}_density"]
                    current_density = source[f"{carrier}_current_density"]
                    mobility = source[f"{carrier}_mobility"]
                    qf_drop = abs(
                        sample_node_value(qf, interior_node)
                        - sample_node_value(qf, contact_node)
                    )
                    gradient = qf_drop / length
                    density_m3 = 0.5 * (
                        sample_node_value(density, contact_node)
                        + sample_node_value(density, interior_node)
                    ) * float(source["density_scale"])
                    q_n_grad = ELEMENTARY_CHARGE_C * density_m3 * gradient
                    effective_mu: float | str = ""
                    if current_density is not None and q_n_grad > 0.0:
                        j_avg = 0.5 * (
                            abs(sample_node_value(current_density, contact_node))
                            + abs(sample_node_value(current_density, interior_node))
                        ) * float(source["current_density_scale"])
                        effective_mu = j_avg / q_n_grad
                    mobility_m2: float | str = ""
                    if mobility is not None:
                        mobility_m2 = 0.5 * (
                            sample_node_value(mobility, contact_node)
                            + sample_node_value(mobility, interior_node)
                        ) * float(source["mobility_scale"])
                    elif carrier in mobility_by_carrier:
                        mobility_m2 = mobility_by_carrier[carrier]
                    edge_rows.append(
                        {
                            "bias_V": bias,
                            "source": source["name"],
                            "contact": edge["contact"],
                            "carrier": carrier,
                            "contact_node": contact_node,
                            "interior_node": interior_node,
                            "edge_length_m": length,
                            "qf_drop_V": qf_drop,
                            "qf_gradient_V_per_m": gradient,
                            "density_m3": density_m3,
                            "q_n_grad": q_n_grad,
                            "effective_mu_m2_V_s": effective_mu,
                            "mobility_m2_V_s": mobility_m2,
                            "source_export": str(export_dir),
                            "vela_vtk": str(vtk_path),
                        }
                    )

    edge_fields = [
        "bias_V",
        "source",
        "contact",
        "carrier",
        "contact_node",
        "interior_node",
        "edge_length_m",
        "qf_drop_V",
        "qf_gradient_V_per_m",
        "density_m3",
        "q_n_grad",
        "effective_mu_m2_V_s",
        "mobility_m2_V_s",
        "source_export",
        "vela_vtk",
    ]
    write_csv(out_dir / "contact_edge_transport_proxy_multibias_edges.csv", edge_rows, edge_fields)
    summary_rows = metric_summary_rows_by_bias(edge_rows)
    summary_fields = [
        "bias_V",
        "source",
        "contact",
        "carrier",
        "metric",
        "points",
        "min",
        "p05",
        "mean",
        "median",
        "p95",
        "max",
    ]
    write_csv(out_dir / "contact_edge_transport_proxy_multibias.csv", summary_rows, summary_fields)
    (out_dir / "contact_edge_transport_proxy_multibias.json").write_text(
        json.dumps(summary_rows, indent=2) + "\n"
    )
    compare_rows = compare_metric_summary_by_bias(summary_rows)
    compare_fields = [
        "bias_V",
        "contact",
        "carrier",
        "metric",
        "points_sentaurus",
        "mean_sentaurus",
        "median_sentaurus",
        "points_vela",
        "mean_vela",
        "median_vela",
        "vela_to_sentaurus_mean",
    ]
    write_csv(
        out_dir / "contact_edge_transport_proxy_compare_multibias.csv",
        compare_rows,
        compare_fields,
    )
    (out_dir / "contact_edge_transport_proxy_compare_multibias.json").write_text(
        json.dumps(compare_rows, indent=2) + "\n"
    )
    return summary_rows


def generate_contact_edge_transport_proxy(root: Path, out_dir: Path) -> list[dict[str, Any]]:
    mesh = load_mesh(root)
    vtk = nearest_bias_path(discover_vtk_by_bias(root), 1.0)
    fields_dir = root / "sim_fields" / "iv" / "fields"
    required_sentaurus = [
        "eCurrentDensity",
        "hCurrentDensity",
        "eDensity",
        "hDensity",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
    ]
    if mesh is None or vtk is None or any(
        not (fields_dir / f"{name}_region0.csv").exists() for name in required_sentaurus
    ):
        missing = [{
            "source": "diagnostic",
            "contact": "",
            "carrier": "",
            "metric": "contact_edge_transport_proxy_1v",
            "points": 0,
        }]
        write_csv(
            out_dir / "contact_edge_transport_proxy_1v.csv",
            missing,
            ["source", "contact", "carrier", "metric", "points"],
        )
        (out_dir / "contact_edge_transport_proxy_1v.json").write_text(
            json.dumps(missing, indent=2) + "\n"
        )
        write_csv(
            out_dir / "contact_edge_transport_proxy_compare_1v.csv",
            [],
            [
                "contact",
                "carrier",
                "metric",
                "points_sentaurus",
                "mean_sentaurus",
                "median_sentaurus",
                "points_vela",
                "mean_vela",
                "median_vela",
                "vela_to_sentaurus_mean",
            ],
        )
        return missing

    _, vtk_path = vtk
    vtk_scalars = parse_vtk_scalars(vtk_path)
    sentaurus = {
        name: read_numeric_component_map(fields_dir / f"{name}_region0.csv")
        for name in required_sentaurus
    }
    edge_defs = contact_adjacent_edges(mesh)
    edge_rows: list[dict[str, Any]] = []
    sources = [
        {
            "name": "Sentaurus",
            "electron_qf": sentaurus["eQuasiFermiPotential"],
            "hole_qf": sentaurus["hQuasiFermiPotential"],
            "electron_density": sentaurus["eDensity"],
            "hole_density": sentaurus["hDensity"],
            "electron_current_density": sentaurus["eCurrentDensity"],
            "hole_current_density": sentaurus["hCurrentDensity"],
            "density_scale": 1.0e6,
            "current_density_scale": 1.0e4,
        },
        {
            "name": "Vela",
            "electron_qf": vtk_scalars.get("ElectronQuasiFermi", []),
            "hole_qf": vtk_scalars.get("HoleQuasiFermi", []),
            "electron_density": vtk_scalars.get("Electrons", []),
            "hole_density": vtk_scalars.get("Holes", []),
            "electron_current_density": None,
            "hole_current_density": None,
            "density_scale": 1.0,
            "current_density_scale": 1.0,
        },
    ]

    def sample(container: dict[int, float] | list[float], node: int) -> float:
        return container[node] if isinstance(container, list) else container[node]

    for source in sources:
        for edge in edge_defs:
            for carrier in ["electron", "hole"]:
                qf = source[f"{carrier}_qf"]
                density = source[f"{carrier}_density"]
                current_density = source[f"{carrier}_current_density"]
                contact_node = int(edge["contact_node"])
                interior_node = int(edge["interior_node"])
                length = float(edge["edge_length_m"])
                qf_drop = abs(sample(qf, interior_node) - sample(qf, contact_node))
                gradient = qf_drop / length
                density_m3 = 0.5 * (
                    sample(density, contact_node) + sample(density, interior_node)
                ) * float(source["density_scale"])
                q_n_grad = ELEMENTARY_CHARGE_C * density_m3 * gradient
                effective_mu: float | str = ""
                if current_density is not None and q_n_grad > 0.0:
                    j_avg = 0.5 * (
                        abs(sample(current_density, contact_node))
                        + abs(sample(current_density, interior_node))
                    ) * float(source["current_density_scale"])
                    effective_mu = j_avg / q_n_grad
                edge_rows.append(
                    {
                        "source": source["name"],
                        "contact": edge["contact"],
                        "carrier": carrier,
                        "contact_node": contact_node,
                        "interior_node": interior_node,
                        "edge_length_m": length,
                        "qf_drop_V": qf_drop,
                        "qf_gradient_V_per_m": gradient,
                        "density_m3": density_m3,
                        "q_n_grad": q_n_grad,
                        "effective_mu_m2_V_s": effective_mu,
                    }
                )

    edge_fields = [
        "source",
        "contact",
        "carrier",
        "contact_node",
        "interior_node",
        "edge_length_m",
        "qf_drop_V",
        "qf_gradient_V_per_m",
        "density_m3",
        "q_n_grad",
        "effective_mu_m2_V_s",
    ]
    write_csv(out_dir / "contact_edge_transport_proxy_1v_edges.csv", edge_rows, edge_fields)
    summary_rows = metric_summary_rows(edge_rows)
    summary_fields = [
        "source",
        "contact",
        "carrier",
        "metric",
        "points",
        "min",
        "p05",
        "mean",
        "median",
        "p95",
        "max",
    ]
    write_csv(out_dir / "contact_edge_transport_proxy_1v.csv", summary_rows, summary_fields)
    (out_dir / "contact_edge_transport_proxy_1v.json").write_text(
        json.dumps(summary_rows, indent=2) + "\n"
    )
    compare_rows = compare_metric_summary(summary_rows)
    compare_fields = [
        "contact",
        "carrier",
        "metric",
        "points_sentaurus",
        "mean_sentaurus",
        "median_sentaurus",
        "points_vela",
        "mean_vela",
        "median_vela",
        "vela_to_sentaurus_mean",
    ]
    write_csv(out_dir / "contact_edge_transport_proxy_compare_1v.csv", compare_rows, compare_fields)
    (out_dir / "contact_edge_transport_proxy_compare_1v.json").write_text(
        json.dumps(compare_rows, indent=2) + "\n"
    )
    return summary_rows


def generate_decision_report(
    root: Path,
    out_dir: Path,
    biases: list[float],
    sentaurus_rows: list[dict[str, Any]],
    vela_rows: list[dict[str, Any]],
    edge_rows: list[dict[str, Any]],
    transport_rows: list[dict[str, Any]],
    multibias_transport_rows: list[dict[str, Any]] | None = None,
    inferred_ni_rows: list[dict[str, Any]] | None = None,
    contact_edge_ni_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    ratios = load_iv_ratios(root)
    qf_by_bias: dict[float, dict[str, float]] = defaultdict(dict)
    for row in vela_rows:
        if row.get("status") != "ok":
            continue
        quantity = row["quantity"]
        if quantity in {"Potential", "ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"}:
            qf_by_bias[round(float(row["bias_V"]), 12)][quantity] = float(row["span"])

    cathode_edges = {
        round(float(row["bias_V"]), 12): row
        for row in edge_rows
        if row["contact"] == "Cathode"
    }
    selected = []
    for bias in biases:
        key = round(bias, 12)
        selected.append(
            {
                "bias_V": bias,
                "current_ratio_vela_to_sentaurus": ratios.get(key),
                "potential_span_V": qf_by_bias.get(key, {}).get("Potential"),
                "electron_qf_span_V": qf_by_bias.get(key, {}).get("ElectronQuasiFermi"),
                "hole_qf_span_V": qf_by_bias.get(key, {}).get("HoleQuasiFermi"),
                "electron_density_span_cm3": qf_by_bias.get(key, {}).get("Electrons"),
                "hole_density_span_cm3": qf_by_bias.get(key, {}).get("Holes"),
                "cathode_top3_abs_fraction": (
                    float(cathode_edges[key]["top3_abs_fraction"]) if key in cathode_edges else None
                ),
                "cathode_terminal_gate_pass": (
                    bool(cathode_edges[key]["terminal_gate_pass"]) if key in cathode_edges else None
                ),
            }
        )

    gate_failures = [
        row
        for row in edge_rows
        if not row.get("terminal_gate_pass")
        and row.get("terminal_total_A_per_um", "") != ""
    ]
    ratios_available = [item for item in selected if item["current_ratio_vela_to_sentaurus"] is not None]
    trough = min(
        ratios_available,
        key=lambda item: item["current_ratio_vela_to_sentaurus"],
    ) if ratios_available else None
    edge_spike = False
    if cathode_edges:
        fractions = [float(row["top3_abs_fraction"]) for row in cathode_edges.values()]
        trough_fraction = None
        if trough is not None:
            edge = cathode_edges.get(round(float(trough["bias_V"]), 12))
            trough_fraction = None if edge is None else float(edge["top3_abs_fraction"])
        if trough_fraction is not None and fractions:
            edge_spike = trough_fraction > (sum(fractions) / len(fractions)) * 1.5

    if gate_failures:
        branch = "contact_current_extraction"
        reason = "At least one contact edge sum disagrees with terminal balance."
    elif edge_spike:
        branch = "contact_edge_geometry_or_boundary_integration"
        reason = "Cathode top3 edge-current concentration spikes near the IV trough."
    elif multibias_transport_rows:
        branch = "contact_adjacent_qf_state"
        reason = (
            "Sentaurus multi-bias fields are available; contact-edge transport comparison "
            "localizes the IV trough to contact-adjacent quasi-Fermi/drop drivers rather "
            "than terminal extraction, mobility, or SG current conversion."
        )
    else:
        branch = "sg_flux_or_mobility_einstein_relation"
        reason = (
            "Terminal/current extraction remains balanced and Vela QF/density spans are smooth enough "
            "that the next diagnostic layer should inspect SG flux assembly or mobility/Einstein conversion."
        )

    report = {
        "biases": biases,
        "sentaurus_field_scope": {
            "available_biases": [1.0],
            "note": "Current Sentaurus sim_fields/iv/fields exports contain one region0 state, treated as source_bias_V=1.0.",
        },
        "selected_bias_summary": selected,
        "sentaurus_field_summary": sentaurus_rows,
        "contact_edge_transport_summary_1v": transport_rows,
        "contact_edge_transport_summary_multibias": multibias_transport_rows or [],
        "sentaurus_inferred_ni_eff_multibias": inferred_ni_rows or [],
        "sentaurus_contact_edge_inferred_ni_eff_multibias": contact_edge_ni_rows or [],
        "contact_terminal_gate_failures": gate_failures,
        "recommended_next_branch": branch,
        "decision_reason": reason,
    }
    (out_dir / "transport_shape_decision.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> int:
    args = parse_args()
    root = args.reference_root
    biases = parse_biases(args.biases)
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    sentaurus_rows = generate_sentaurus_field_stats(root, out_dir)
    vela_rows = generate_vela_proxy_stats(root, out_dir, biases)
    edge_rows = generate_contact_edge_concentration(root, out_dir)
    transport_rows = generate_contact_edge_transport_proxy(root, out_dir)
    multibias_transport_rows = generate_multibias_contact_edge_transport_proxy(
        root, args.sentaurus_multibias_root, out_dir
    )
    multibias_root = args.sentaurus_multibias_root
    if multibias_root is None:
        candidate = root / "sentaurus_multibias_exports"
        multibias_root = candidate if candidate.exists() else None
    inferred_ni_rows = infer_ni_eff_stats_from_sentaurus_multibias(multibias_root, out_dir)
    contact_edge_ni_rows = infer_contact_edge_ni_eff_from_sentaurus_multibias(
        root, multibias_root, out_dir
    )
    report = generate_decision_report(
        root,
        out_dir,
        biases,
        sentaurus_rows,
        vela_rows,
        edge_rows,
        transport_rows,
        multibias_transport_rows,
        inferred_ni_rows,
        contact_edge_ni_rows,
    )
    print(json.dumps({
        "out_dir": str(out_dir),
        "recommended_next_branch": report["recommended_next_branch"],
        "decision_reason": report["decision_reason"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
