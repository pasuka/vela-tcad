#!/usr/bin/env python3
"""Diagnose PN2D BV contact-edge current extraction on existing outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


EDGE_SUMMARY_FIELDS = [
    "bias_V",
    "contact",
    "edge_rows",
    "skipped_incomplete_rows",
    "edge_current_A_per_um",
    "edge_electron_current_A_per_um",
    "edge_hole_current_A_per_um",
    "vela_iv_current_A_per_um",
    "edge_vs_iv_relative_error",
    "sentaurus_current_A",
    "log10_vela_over_sentaurus_current",
    "sentaurus_plt_current_A",
    "sentaurus_reference_vs_plt_relative_error",
    "sentaurus_contact_current_flux_A",
    "sentaurus_contact_current_flux_region",
    "sentaurus_contact_flux_vs_plt_relative_error",
    "sentaurus_terminal_crosscheck_status",
    "sentaurus_current_field_status",
    "classification",
    "current_gap_branch",
]


PASSTHROUGH_FIELDS = [
    "bias_V",
    "current_contact",
    "edge_id",
    "node0",
    "node1",
    "psi0",
    "psi1",
    "phin0",
    "phin1",
    "phip0",
    "phip1",
    "n0",
    "n1",
    "p0",
    "p1",
    "mun",
    "mup",
    "electron_continuity_flux",
    "hole_continuity_flux",
    "current_electron",
    "current_hole",
    "current_total_A_per_um",
    "edge_mid_x_um",
    "edge_mid_y_um",
    "sentaurus_nearest_node_id",
    "sentaurus_nearest_distance_um",
    "sentaurus_total_current_density",
    "sentaurus_electron_current_density",
    "sentaurus_hole_current_density",
    "sentaurus_electron_density_cm3",
    "sentaurus_hole_density_cm3",
    "sentaurus_electron_mobility_cm2_V_s",
    "sentaurus_hole_mobility_cm2_V_s",
    "sentaurus_potential_V",
    "sentaurus_electron_qf_V",
    "sentaurus_hole_qf_V",
    "sentaurus_potential_drop_V",
    "sentaurus_electron_qf_drop_V",
    "sentaurus_hole_qf_drop_V",
    "vela_potential_avg_V",
    "vela_electron_qf_avg_V",
    "vela_hole_qf_avg_V",
    "vela_potential_drop_V",
    "vela_electron_qf_drop_V",
    "vela_hole_qf_drop_V",
    "sentaurus_over_vela_hole_qf_drop",
    "vela_electron_density_avg_cm3",
    "vela_hole_density_avg_cm3",
    "vela_electron_mobility_cm2_V_s",
    "vela_hole_mobility_cm2_V_s",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-edge-csv", type=Path, required=True)
    parser.add_argument("--curve-candidate", type=Path, required=True)
    parser.add_argument("--curve-reference", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", required=True, help="Comma-separated biases in volts.")
    parser.add_argument("--contact", default="Anode")
    parser.add_argument(
        "--parity-rtol",
        type=float,
        default=0.01,
        help="Relative tolerance for contact-edge sum versus Vela IV current.",
    )
    parser.add_argument(
        "--sentaurus-root",
        type=Path,
        default=None,
        help="Root containing sentaurus_<bias>v exports for current-density field alignment.",
    )
    parser.add_argument("--sentaurus-plt", type=Path, default=None)
    parser.add_argument("--plt-bias-column", default=None)
    parser.add_argument("--plt-current-column", default=None)
    parser.add_argument("--mesh", type=Path, default=None, help="Vela mesh JSON for edge midpoint alignment.")
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def near_bias(value: float, targets: set[float]) -> bool:
    key = bias_key(value)
    return any(abs(key - target) <= 1.0e-9 for target in targets)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def optional_float(raw: str | None) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def load_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def parse_quoted_list(text: str, key: str) -> list[str]:
    match = re.search(rf"{re.escape(key)}\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def parse_plt_values_block(text: str, column_count: int) -> list[list[float]]:
    match = re.search(r"Values\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        match = re.search(r"Data\s*\{(.*?)\}", text, re.DOTALL | re.IGNORECASE)
    if not match:
        return []
    numbers = [float(value) for value in re.findall(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        match.group(1),
    )]
    if column_count <= 0 or len(numbers) % column_count != 0:
        return []
    return [
        numbers[index:index + column_count]
        for index in range(0, len(numbers), column_count)
    ]


def load_plt_points(path: Path | None, bias_column: str, current_column: str) -> dict[float, float]:
    if path is None or not path.exists():
        return {}
    text = path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if bias_column not in datasets or current_column not in datasets:
        return {}
    bias_index = datasets.index(bias_column)
    current_index = datasets.index(current_column)
    points: dict[float, float] = {}
    for row in parse_plt_values_block(text, len(datasets)):
        points[bias_key(row[bias_index])] = row[current_index]
    return points


def load_mesh_nodes(path: Path | None) -> dict[int, tuple[float, float]]:
    if path is None or not path.exists():
        return {}
    data = json.loads(path.read_text())
    nodes: dict[int, tuple[float, float]] = {}
    for node in data.get("nodes", []):
        nodes[int(node.get("id", len(nodes)))] = (float(node["x"]), float(node["y"]))
    return nodes


def discover_sentaurus_dir(root: Path | None, bias: float) -> Path | None:
    if root is None:
        return None
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "nodes.csv").exists():
            return candidate
    return None


def load_scalar_field(export_dir: Path, name: str) -> dict[int, float]:
    path = export_dir / "fields" / f"{name}_region0.csv"
    if not path.exists():
        return {}
    _, rows = load_csv(path)
    values: dict[int, float] = {}
    for row in rows:
        node_id = optional_float(row.get("node_id"))
        value = optional_float(row.get("component0"))
        if node_id is not None and value is not None:
            values[int(node_id)] = value
    return values


def load_contact_current_flux(export_dir: Path, contact: str) -> dict[str, Any]:
    manifest_path = export_dir / "field_manifest.json"
    if not manifest_path.exists():
        return {"status": "contact_current_flux_manifest_missing"}
    try:
        manifest = json.loads(manifest_path.read_text())
    except json.JSONDecodeError:
        return {"status": "contact_current_flux_manifest_invalid"}
    datasets = manifest.get("datasets") or manifest.get("fields") or (
        manifest if isinstance(manifest, list) else []
    )
    for item in datasets:
        if item.get("name") != "ContactCurrentFlux":
            continue
        if str(item.get("region_name", "")).lower() != contact.lower():
            continue
        region = item.get("region")
        path = export_dir / "fields" / f"ContactCurrentFlux_region{region}.csv"
        if not path.exists():
            return {
                "status": "contact_current_flux_csv_missing",
                "region": region,
            }
        _, rows = load_csv(path)
        for row in rows:
            value = optional_float(row.get("component0"))
            if value is not None:
                return {
                    "status": "contact_current_flux_available",
                    "region": region,
                    "value_A": value,
                }
        return {
            "status": "contact_current_flux_value_missing",
            "region": region,
        }
    return {"status": "contact_current_flux_contact_missing"}


def load_sentaurus_export(export_dir: Path | None, contact: str = "Anode") -> dict[str, Any] | None:
    if export_dir is None:
        return None
    _, node_rows = load_csv(export_dir / "nodes.csv")
    nodes: dict[int, tuple[float, float]] = {}
    for row in node_rows:
        node_id = optional_float(row.get("id"))
        x_um = optional_float(row.get("x_um"))
        y_um = optional_float(row.get("y_um"))
        if node_id is not None and x_um is not None and y_um is not None:
            nodes[int(node_id)] = (x_um, y_um)
    fields = {
        "total_current_density": load_scalar_field(export_dir, "TotalCurrentDensity"),
        "electron_current_density": load_scalar_field(export_dir, "eCurrentDensity"),
        "hole_current_density": load_scalar_field(export_dir, "hCurrentDensity"),
        "electron_density": load_scalar_field(export_dir, "eDensity"),
        "hole_density": load_scalar_field(export_dir, "hDensity"),
        "electron_mobility": load_scalar_field(export_dir, "eMobility"),
        "hole_mobility": load_scalar_field(export_dir, "hMobility"),
        "potential": load_scalar_field(export_dir, "ElectrostaticPotential"),
        "electron_qf": load_scalar_field(export_dir, "eQuasiFermiPotential"),
        "hole_qf": load_scalar_field(export_dir, "hQuasiFermiPotential"),
    }
    return {
        "nodes": nodes,
        "fields": fields,
        "contact_current_flux": load_contact_current_flux(export_dir, contact),
    }


def nearest_sentaurus_node(
    export: dict[str, Any],
    x_um: float,
    y_um: float,
) -> tuple[int | None, float | None]:
    nodes: dict[int, tuple[float, float]] = export["nodes"]
    if not nodes:
        return None, None
    best_id = min(
        nodes,
        key=lambda node_id: (nodes[node_id][0] - x_um) ** 2 + (nodes[node_id][1] - y_um) ** 2,
    )
    best = nodes[best_id]
    return best_id, math.hypot(best[0] - x_um, best[1] - y_um)


def load_curve_points(path: Path, prefer_per_um: bool) -> dict[float, float]:
    if not path.exists():
        return {}
    fields, rows = load_csv(path)
    if not rows:
        return {}
    bias_col = next((name for name in fields if name in {"bias_V", "voltage", "V"}), fields[0])
    if prefer_per_um:
        preferred = ["current_total_A_per_um", "current_total", "current"]
    else:
        preferred = ["current_total", "current", "current_total_A_per_um"]
    current_col = next((name for name in preferred if name in fields), fields[-1])
    points: dict[float, float] = {}
    for row in rows:
        converged = str(row.get("converged", "1")).strip()
        if converged not in {"1", "true", "True", ""}:
            continue
        bias = optional_float(row.get(bias_col))
        current = optional_float(row.get(current_col))
        if bias is None or current is None:
            continue
        points.setdefault(bias_key(bias), current)
    return points


def current_a_per_um(row: dict[str, str], total_key: str) -> float | None:
    value = optional_float(row.get(f"{total_key}_A_per_um"))
    if value is not None:
        return value
    value = optional_float(row.get(total_key))
    if value is None:
        return None
    return value / 1.0e6


def row_float(row: dict[str, str], key: str) -> float | None:
    return optional_float(row.get(key))


def annotate_edge_row(
    row: dict[str, str],
    mesh_nodes: dict[int, tuple[float, float]],
    sentaurus_exports: dict[float, dict[str, Any]],
) -> dict[str, Any]:
    annotated: dict[str, Any] = dict(row)
    bias = row_float(row, "bias_V")
    node0 = row_float(row, "node0")
    node1 = row_float(row, "node1")
    if bias is None or node0 is None or node1 is None:
        return annotated
    p0 = mesh_nodes.get(int(node0))
    p1 = mesh_nodes.get(int(node1))
    if p0 is None or p1 is None:
        return annotated
    x_um = 0.5 * (p0[0] + p1[0])
    y_um = 0.5 * (p0[1] + p1[1])
    annotated["edge_mid_x_um"] = x_um
    annotated["edge_mid_y_um"] = y_um

    n0 = row_float(row, "n0")
    n1 = row_float(row, "n1")
    p0_density = row_float(row, "p0")
    p1_density = row_float(row, "p1")
    psi0 = row_float(row, "psi0")
    psi1 = row_float(row, "psi1")
    phin0 = row_float(row, "phin0")
    phin1 = row_float(row, "phin1")
    phip0 = row_float(row, "phip0")
    phip1 = row_float(row, "phip1")
    mun = row_float(row, "mun")
    mup = row_float(row, "mup")
    if psi0 is not None and psi1 is not None:
        annotated["vela_potential_avg_V"] = 0.5 * (psi0 + psi1)
        annotated["vela_potential_drop_V"] = psi1 - psi0
    if phin0 is not None and phin1 is not None:
        annotated["vela_electron_qf_avg_V"] = 0.5 * (phin0 + phin1)
        annotated["vela_electron_qf_drop_V"] = phin1 - phin0
    if phip0 is not None and phip1 is not None:
        annotated["vela_hole_qf_avg_V"] = 0.5 * (phip0 + phip1)
        annotated["vela_hole_qf_drop_V"] = phip1 - phip0
    if n0 is not None and n1 is not None:
        annotated["vela_electron_density_avg_cm3"] = 0.5 * (n0 + n1) * 1.0e-6
    if p0_density is not None and p1_density is not None:
        annotated["vela_hole_density_avg_cm3"] = 0.5 * (p0_density + p1_density) * 1.0e-6
    if mun is not None:
        annotated["vela_electron_mobility_cm2_V_s"] = mun * 1.0e4
    if mup is not None:
        annotated["vela_hole_mobility_cm2_V_s"] = mup * 1.0e4

    export = sentaurus_exports.get(bias_key(bias))
    if export is None:
        return annotated
    nearest_id, distance_um = nearest_sentaurus_node(export, x_um, y_um)
    if nearest_id is None:
        return annotated
    annotated["sentaurus_nearest_node_id"] = nearest_id
    annotated["sentaurus_nearest_distance_um"] = distance_um
    fields = export["fields"]
    mapping = {
        "sentaurus_total_current_density": ("total_current_density", nearest_id),
        "sentaurus_electron_current_density": ("electron_current_density", nearest_id),
        "sentaurus_hole_current_density": ("hole_current_density", nearest_id),
        "sentaurus_electron_density_cm3": ("electron_density", nearest_id),
        "sentaurus_hole_density_cm3": ("hole_density", nearest_id),
        "sentaurus_electron_mobility_cm2_V_s": ("electron_mobility", nearest_id),
        "sentaurus_hole_mobility_cm2_V_s": ("hole_mobility", nearest_id),
        "sentaurus_potential_V": ("potential", nearest_id),
        "sentaurus_electron_qf_V": ("electron_qf", nearest_id),
        "sentaurus_hole_qf_V": ("hole_qf", nearest_id),
    }
    for out_key, (field_key, node_id) in mapping.items():
        value = fields[field_key].get(node_id)
        if value is not None:
            annotated[out_key] = value
    endpoint_drops = {
        "sentaurus_potential_drop_V": "potential",
        "sentaurus_electron_qf_drop_V": "electron_qf",
        "sentaurus_hole_qf_drop_V": "hole_qf",
    }
    for out_key, field_key in endpoint_drops.items():
        v0 = fields[field_key].get(int(node0))
        v1 = fields[field_key].get(int(node1))
        if v0 is not None and v1 is not None:
            annotated[out_key] = v1 - v0
    sent_hole_drop = optional_float(str(annotated.get("sentaurus_hole_qf_drop_V", "")))
    vela_hole_drop = optional_float(str(annotated.get("vela_hole_qf_drop_V", "")))
    if sent_hole_drop is not None and vela_hole_drop not in (None, 0.0):
        annotated["sentaurus_over_vela_hole_qf_drop"] = sent_hole_drop / vela_hole_drop
    return annotated


def complete_edge_row(row: dict[str, str]) -> bool:
    if optional_float(row.get("bias_V")) is None:
        return False
    if not row.get("current_contact"):
        return False
    return current_a_per_um(row, "current_total") is not None


def relative_error(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None:
        return None
    scale = max(abs(candidate), abs(reference), 1.0e-300)
    return abs(candidate - reference) / scale


def log10_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None:
        return None
    if candidate == 0.0 or reference == 0.0:
        return None
    return math.log10(abs(candidate) / abs(reference))


def classify(edge_current: float | None, iv_current: float | None, rel: float | None, rtol: float) -> str:
    if edge_current is None:
        return "missing_contact_edges"
    if iv_current is None:
        return "vela_iv_current_missing"
    if rel is not None and rel <= rtol:
        return "terminal_extraction_consistent"
    return "terminal_extraction_mismatch"


def terminal_crosscheck_status(
    plt_current: float | None,
    contact_flux: float | None,
    contact_flux_vs_plt: float | None,
) -> str:
    if plt_current is None and contact_flux is None:
        return "sentaurus_terminal_crosscheck_unavailable"
    if plt_current is None:
        return "sentaurus_plt_unavailable"
    if contact_flux is None:
        return "sentaurus_contact_flux_unavailable"
    if contact_flux_vs_plt is not None and contact_flux_vs_plt <= 0.1:
        return "sentaurus_plt_contact_flux_consistent"
    return "sentaurus_plt_contact_flux_mismatch"


def classify_current_gap(
    extraction_classification: str,
    log_ratio: float | None,
    crosscheck_status: str,
) -> str:
    if extraction_classification != "terminal_extraction_consistent":
        return "vela_terminal_extraction_mismatch"
    if crosscheck_status == "sentaurus_plt_contact_flux_mismatch":
        return "sentaurus_terminal_definition_mismatch"
    if log_ratio is None:
        return "current_gap_unclassified"
    if abs(log_ratio) <= 0.1:
        return "current_parity_close"
    if crosscheck_status == "sentaurus_plt_contact_flux_consistent":
        return "transport_or_field_mismatch"
    return "needs_sentaurus_terminal_crosscheck"


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    biases = parse_biases(args.biases)
    target_keys = {bias_key(bias) for bias in biases}
    fields, edge_rows = load_csv(args.contact_edge_csv)
    candidate_curve = load_curve_points(args.curve_candidate, prefer_per_um=True)
    reference_curve = load_curve_points(args.curve_reference, prefer_per_um=False)
    plt_bias_column = args.plt_bias_column or f"{args.contact} OuterVoltage"
    plt_current_column = args.plt_current_column or f"{args.contact} TotalCurrent"
    sentaurus_plt_curve = load_plt_points(args.sentaurus_plt, plt_bias_column, plt_current_column)
    mesh_nodes = load_mesh_nodes(args.mesh)
    sentaurus_exports = {
        bias_key(bias): export
        for bias in biases
        if (export := load_sentaurus_export(discover_sentaurus_dir(args.sentaurus_root, bias), args.contact)) is not None
    }

    skipped_incomplete = 0
    selected: list[dict[str, str]] = []
    grouped: dict[tuple[float, str], list[dict[str, str]]] = defaultdict(list)
    for row in edge_rows:
        bias = optional_float(row.get("bias_V"))
        if bias is None or not near_bias(bias, target_keys):
            continue
        if not complete_edge_row(row):
            skipped_incomplete += 1
            continue
        if row.get("current_contact") != args.contact:
            continue
        row = annotate_edge_row(row, mesh_nodes, sentaurus_exports)
        key = (bias_key(bias), row["current_contact"])
        grouped[key].append(row)
        selected.append(row)

    filtered_fields = list(dict.fromkeys([*fields, *PASSTHROUGH_FIELDS]))
    summary_rows: list[dict[str, Any]] = []
    json_rows: list[dict[str, Any]] = []
    for bias in biases:
        key = (bias_key(bias), args.contact)
        rows = grouped.get(key, [])
        edge_current = sum(current_a_per_um(row, "current_total") or 0.0 for row in rows) if rows else None
        electron_current = sum(current_a_per_um(row, "current_electron") or 0.0 for row in rows) if rows else None
        hole_current = sum(current_a_per_um(row, "current_hole") or 0.0 for row in rows) if rows else None
        iv_current = candidate_curve.get(bias_key(bias))
        sentaurus_current = reference_curve.get(bias_key(bias))
        sentaurus_plt_current = sentaurus_plt_curve.get(bias_key(bias))
        export = sentaurus_exports.get(bias_key(bias))
        contact_flux_info = export.get("contact_current_flux", {}) if export is not None else {}
        contact_flux = contact_flux_info.get("value_A")
        contact_flux_region = contact_flux_info.get("region")
        rel = relative_error(edge_current, iv_current)
        log_ratio = log10_ratio(iv_current, sentaurus_current)
        reference_vs_plt = relative_error(sentaurus_current, sentaurus_plt_current)
        flux_vs_plt = relative_error(contact_flux, sentaurus_plt_current)
        classification = classify(edge_current, iv_current, rel, args.parity_rtol)
        crosscheck_status = terminal_crosscheck_status(
            sentaurus_plt_current,
            contact_flux,
            flux_vs_plt,
        )
        sentaurus_status = (
            "sentaurus_current_field_available"
            if rows and any(row.get("sentaurus_total_current_density") not in (None, "") for row in rows)
            else "sentaurus_current_field_unavailable"
        )
        current_gap_branch = classify_current_gap(classification, log_ratio, crosscheck_status)
        summary = {
            "bias_V": bias,
            "contact": args.contact,
            "edge_rows": len(rows),
            "skipped_incomplete_rows": skipped_incomplete,
            "edge_current_A_per_um": edge_current,
            "edge_electron_current_A_per_um": electron_current,
            "edge_hole_current_A_per_um": hole_current,
            "vela_iv_current_A_per_um": iv_current,
            "edge_vs_iv_relative_error": rel,
            "sentaurus_current_A": sentaurus_current,
            "log10_vela_over_sentaurus_current": log_ratio,
            "sentaurus_plt_current_A": sentaurus_plt_current,
            "sentaurus_reference_vs_plt_relative_error": reference_vs_plt,
            "sentaurus_contact_current_flux_A": contact_flux,
            "sentaurus_contact_current_flux_region": contact_flux_region,
            "sentaurus_contact_flux_vs_plt_relative_error": flux_vs_plt,
            "sentaurus_terminal_crosscheck_status": crosscheck_status,
            "sentaurus_current_field_status": sentaurus_status,
            "classification": classification,
            "current_gap_branch": current_gap_branch,
        }
        summary_rows.append(summary)
        json_rows.append(summary)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "contact_edges_filtered.csv", filtered_fields, selected)
    write_csv(args.out_dir / "edge_summary.csv", EDGE_SUMMARY_FIELDS, summary_rows)
    (args.out_dir / "summary.json").write_text(
        json.dumps({
            "contact_edge_csv": str(args.contact_edge_csv),
            "curve_candidate": str(args.curve_candidate),
            "curve_reference": str(args.curve_reference),
            "sentaurus_plt": str(args.sentaurus_plt) if args.sentaurus_plt else None,
            "plt_bias_column": plt_bias_column,
            "plt_current_column": plt_current_column,
            "contact": args.contact,
            "biases": biases,
            "skipped_incomplete_rows": skipped_incomplete,
            "rows": json_rows,
        }, indent=2)
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
