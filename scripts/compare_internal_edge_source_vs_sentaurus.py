#!/usr/bin/env python3
"""Compare Vela internal avalanche edge source rows with Sentaurus edge-equivalent sources."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median


Q = 1.602176634e-19
RATIO_FLOOR = 1.0e-300


WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}


@dataclass
class EdgeTopology:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    edge_length_um: float
    source_area_cm2: float

    @property
    def x_mid_um(self) -> float:
        return 0.5 * (self.x0_um + self.x1_um)

    @property
    def y_mid_um(self) -> float:
        return 0.5 * (self.y0_um + self.y1_um)


@dataclass
class SentaurusBiasFields:
    bias: float
    electron_qf_V: dict[int, float]
    hole_qf_V: dict[int, float]
    electron_alpha_cm_inv: dict[int, float]
    hole_alpha_cm_inv: dict[int, float]
    electron_current_A_cm2: dict[int, tuple[float, float]]
    hole_current_A_cm2: dict[int, tuple[float, float]]
    gava_cm3_s: dict[int, float]


@dataclass
class WindowStats:
    rows: int = 0
    qg_self: float = 0.0
    qg_sentaurus_a: float = 0.0
    qg_sentaurus_b: float = 0.0
    g_self_max: float = 0.0
    g_sentaurus_b_max: float = 0.0
    ratio_logs_f: list[float] = None
    ratio_logs_alpha: list[float] = None
    ratio_logs_j: list[float] = None
    ratio_logs_g: list[float] = None

    def __post_init__(self) -> None:
        self.ratio_logs_f = []
        self.ratio_logs_alpha = []
        self.ratio_logs_j = []
        self.ratio_logs_g = []


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    default_internal = repo / "build/diagnostics/avalanche_internal_source_current_audit.csv"
    default_summary = repo / "build/diagnostics/avalanche_internal_source_current_audit_summary.md"
    default_sentaurus = repo / "build/diagnostics/pn2d_bv_codex_compare/reports/sentaurus_multibias"
    default_elements = repo / "build/diagnostics/pn2d_bv_codex_compare/imported_reference/elements.csv"
    default_nodes = repo / "build/diagnostics/pn2d_bv_codex_compare/imported_reference/nodes.csv"
    default_mesh_json = repo / "build/diagnostics/pn2d_bv_codex_compare/imported_reference/vela/mesh.json"
    default_node_compare = (
        repo / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
        / "coarse_vm_vector_compare/coarse_node_field_compare_aligned.csv"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-audit", type=Path, default=default_internal)
    parser.add_argument("--internal-summary", type=Path, default=default_summary)
    parser.add_argument("--self-edge-topology", type=Path)
    parser.add_argument("--sentaurus-multibias-dir", type=Path, default=default_sentaurus)
    parser.add_argument("--node-compare", type=Path, default=default_node_compare)
    parser.add_argument("--elements", type=Path, default=default_elements)
    parser.add_argument("--nodes", type=Path, default=default_nodes)
    parser.add_argument("--mesh-json", type=Path, default=default_mesh_json)
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=repo / "build/diagnostics/internal_edge_source_vs_sentaurus.csv",
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=repo / "build/diagnostics/internal_edge_source_vs_sentaurus_summary.md",
    )
    return parser.parse_args()


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {key.strip(): value.strip() for key, value in row.items() if key is not None}


def parse_float(value: str, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def parse_biases(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def nearest_bias(value: float, targets: list[float], tol: float) -> float | None:
    if not targets:
        return value
    best = min(targets, key=lambda bias: abs(value - bias))
    return best if abs(value - best) <= tol else None


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if abs(denominator) <= RATIO_FLOOR or not math.isfinite(denominator):
        return None
    ratio = numerator / denominator
    return ratio if math.isfinite(ratio) else None


def ratio_text(value: float | None) -> str:
    return "" if value is None else f"{value:.17g}"


def add_ratio_log(bucket: list[float], ratio: float | None) -> None:
    if ratio is not None and ratio > 0.0 and math.isfinite(ratio):
        bucket.append(abs(math.log10(ratio)))


def vector_magnitude(vector: tuple[float, float]) -> float:
    return math.hypot(vector[0], vector[1])


def load_vector_field(path: Path) -> dict[int, tuple[float, float]]:
    values: dict[int, tuple[float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = clean_row(raw)
            values[int(row["node_id"])] = (
                parse_float(row.get("component0", "0")),
                parse_float(row.get("component1", "0")),
            )
    return values


def load_scalar_field(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = clean_row(raw)
            values[int(row["node_id"])] = parse_float(row.get("component0", "0"))
    return values


def sentaurus_bias_from_dir(path: Path) -> float | None:
    match = re.match(r"sentaurus_([+-]?\d+(?:\.\d+)?)v$", path.name, re.IGNORECASE)
    if match:
        return float(match.group(1))
    return None


def find_sentaurus_dirs(root: Path) -> dict[float, Path]:
    dirs: dict[float, Path] = {}
    if not root.exists():
        return dirs
    for child in root.iterdir():
        if not child.is_dir():
            continue
        bias = sentaurus_bias_from_dir(child)
        if bias is not None:
            dirs[bias] = child
    return dirs


def load_sentaurus_bias_fields(root: Path, requested_biases: list[float]) -> dict[float, SentaurusBiasFields]:
    dirs = find_sentaurus_dirs(root)
    fields: dict[float, SentaurusBiasFields] = {}
    for target in requested_biases:
        if not dirs:
            continue
        actual = min(dirs, key=lambda bias: abs(bias - target))
        field_dir = dirs[actual] / "fields"
        fields[target] = SentaurusBiasFields(
            bias=actual,
            electron_qf_V=load_scalar_field(field_dir / "eQuasiFermiPotential_region0.csv"),
            hole_qf_V=load_scalar_field(field_dir / "hQuasiFermiPotential_region0.csv"),
            electron_alpha_cm_inv=load_scalar_field(field_dir / "eAlphaAvalanche_region0.csv"),
            hole_alpha_cm_inv=load_scalar_field(field_dir / "hAlphaAvalanche_region0.csv"),
            electron_current_A_cm2=load_vector_field(field_dir / "eCurrentDensity_region0.csv"),
            hole_current_A_cm2=load_vector_field(field_dir / "hCurrentDensity_region0.csv"),
            gava_cm3_s=load_scalar_field(field_dir / "ImpactIonization_region0.csv"),
        )
    return fields


def load_topology_from_sg_edges(path: Path) -> dict[int, EdgeTopology]:
    topology: dict[int, EdgeTopology] = {}
    if not path.exists():
        return topology
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = clean_row(raw)
            edge_id = int(row["edge_id"])
            if edge_id in topology:
                continue
            length_um = parse_float(row["edge_length_m"]) * 1.0e6
            area_cm2 = parse_float(row.get("edge_area_proxy_m2", "0")) * 1.0e4
            topology[edge_id] = EdgeTopology(
                edge_id=edge_id,
                node0=int(row["node0"]),
                node1=int(row["node1"]),
                x0_um=parse_float(row["x0_um"]),
                y0_um=parse_float(row["y0_um"]),
                x1_um=parse_float(row["x1_um"]),
                y1_um=parse_float(row["y1_um"]),
                edge_length_um=length_um,
                source_area_cm2=area_cm2,
            )
    return topology


def load_nodes_csv(path: Path) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = clean_row(raw)
            nodes[int(row["id"])] = (parse_float(row["x_um"]), parse_float(row["y_um"]))
    return nodes


def load_nodes_mesh_json(path: Path) -> dict[int, tuple[float, float]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {int(node["id"]): (float(node["x"]), float(node["y"])) for node in data["nodes"]}


def add_edge_from_nodes(
    topology: dict[int, EdgeTopology],
    seen: set[tuple[int, int]],
    nodes: dict[int, tuple[float, float]],
    a: int,
    b: int,
) -> None:
    key = tuple(sorted((a, b)))
    if key in seen:
        return
    seen.add(key)
    edge_id = len(topology)
    node0, node1 = key
    x0, y0 = nodes[node0]
    x1, y1 = nodes[node1]
    length_um = math.hypot(x1 - x0, y1 - y0)
    topology[edge_id] = EdgeTopology(
        edge_id=edge_id,
        node0=node0,
        node1=node1,
        x0_um=x0,
        y0_um=y0,
        x1_um=x1,
        y1_um=y1,
        edge_length_um=length_um,
        source_area_cm2=0.0,
    )


def load_topology_from_mesh_json(path: Path) -> dict[int, EdgeTopology]:
    data = json.loads(path.read_text(encoding="utf-8"))
    nodes = {int(node["id"]): (float(node["x"]), float(node["y"])) for node in data["nodes"]}
    topology: dict[int, EdgeTopology] = {}
    seen: set[tuple[int, int]] = set()
    for tri in data.get("triangles", []):
        tri_nodes = [int(node_id) for node_id in tri["node_ids"]]
        for a, b in ((tri_nodes[0], tri_nodes[1]), (tri_nodes[1], tri_nodes[2]), (tri_nodes[2], tri_nodes[0])):
            add_edge_from_nodes(topology, seen, nodes, a, b)
    return topology


def load_topology_from_elements(elements: Path, nodes_path: Path | None, mesh_json: Path | None) -> dict[int, EdgeTopology]:
    if nodes_path is not None and nodes_path.exists():
        nodes = load_nodes_csv(nodes_path)
    elif mesh_json is not None and mesh_json.exists():
        nodes = load_nodes_mesh_json(mesh_json)
    else:
        return {}

    topology: dict[int, EdgeTopology] = {}
    seen: set[tuple[int, int]] = set()
    with elements.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = clean_row(raw)
            tri_nodes = [int(row["node0"]), int(row["node1"]), int(row["node2"])]
            for a, b in ((tri_nodes[0], tri_nodes[1]), (tri_nodes[1], tri_nodes[2]), (tri_nodes[2], tri_nodes[0])):
                add_edge_from_nodes(topology, seen, nodes, a, b)
    return topology


def load_topology(args: argparse.Namespace) -> dict[int, EdgeTopology]:
    if args.self_edge_topology is not None and args.self_edge_topology.exists():
        topology = load_topology_from_sg_edges(args.self_edge_topology)
        if topology:
            return topology
    if args.mesh_json is not None and args.mesh_json.exists():
        topology = load_topology_from_mesh_json(args.mesh_json)
        if topology:
            return topology
    if args.elements is not None and args.elements.exists():
        topology = load_topology_from_elements(args.elements, args.nodes, args.mesh_json)
        return topology
    return {}


def sentaurus_edge_terms(edge: EdgeTopology, fields: SentaurusBiasFields) -> dict[str, float]:
    node0 = edge.node0
    node1 = edge.node1
    length_cm = edge.edge_length_um * 1.0e-4
    fn = abs(fields.electron_qf_V[node1] - fields.electron_qf_V[node0]) / length_cm if length_cm > 0.0 else 0.0
    fp = abs(fields.hole_qf_V[node1] - fields.hole_qf_V[node0]) / length_cm if length_cm > 0.0 else 0.0
    alpha_n = 0.5 * (fields.electron_alpha_cm_inv[node0] + fields.electron_alpha_cm_inv[node1])
    alpha_p = 0.5 * (fields.hole_alpha_cm_inv[node0] + fields.hole_alpha_cm_inv[node1])
    jn_vec = (
        0.5 * (fields.electron_current_A_cm2[node0][0] + fields.electron_current_A_cm2[node1][0]),
        0.5 * (fields.electron_current_A_cm2[node0][1] + fields.electron_current_A_cm2[node1][1]),
    )
    jp_vec = (
        0.5 * (fields.hole_current_A_cm2[node0][0] + fields.hole_current_A_cm2[node1][0]),
        0.5 * (fields.hole_current_A_cm2[node0][1] + fields.hole_current_A_cm2[node1][1]),
    )
    jn = vector_magnitude(jn_vec)
    jp = vector_magnitude(jp_vec)
    g_a = (alpha_n * jn + alpha_p * jp) / Q
    g_b = 0.5 * (fields.gava_cm3_s[node0] + fields.gava_cm3_s[node1])
    qg_a = Q * g_a * edge.source_area_cm2 / 1.0e4
    qg_b = Q * g_b * edge.source_area_cm2 / 1.0e4
    return {
        "Fn_sentaurus_edge_V_per_cm": fn,
        "Fp_sentaurus_edge_V_per_cm": fp,
        "alpha_n_sentaurus_edge_cm_inv": alpha_n,
        "alpha_p_sentaurus_edge_cm_inv": alpha_p,
        "Jn_sentaurus_edge_A_per_cm2": jn,
        "Jp_sentaurus_edge_A_per_cm2": jp,
        "Gava_sentaurus_edge_method_A_cm_minus3_s_minus1": g_a,
        "Gava_sentaurus_edge_method_B_cm_minus3_s_minus1": g_b,
        "qG_sentaurus_edge_method_A_A_per_um": qg_a,
        "qG_sentaurus_edge_method_B_A_per_um": qg_b,
    }


def window_for_x(x_um: float) -> list[str]:
    return [name for name, (lo, hi) in WINDOWS.items() if lo <= x_um <= hi]


def dominant_factor(stats: WindowStats) -> str:
    errors = {
        "F": median(stats.ratio_logs_f) if stats.ratio_logs_f else 0.0,
        "alpha": median(stats.ratio_logs_alpha) if stats.ratio_logs_alpha else 0.0,
        "J": median(stats.ratio_logs_j) if stats.ratio_logs_j else 0.0,
    }
    qg_ratio = safe_ratio(stats.qg_self, stats.qg_sentaurus_a)
    g_error = median(stats.ratio_logs_g) if stats.ratio_logs_g else 0.0
    qg_error = abs(math.log10(qg_ratio)) if qg_ratio is not None and qg_ratio > 0.0 else 0.0
    component_name, component_error = max(errors.items(), key=lambda item: item[1])
    if qg_error > max(component_error, g_error) + 0.25 and g_error <= component_error + 0.25:
        return "source_area_or_edge_mapping"
    return component_name


def parse_internal_qg_relative_error(path: Path) -> float | None:
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("- qG_internal_vs_solver_relative_error:"):
            return float(line.split(":", 1)[1].strip())
    return None


def write_summary(
    path: Path,
    *,
    internal_summary_error: float | None,
    stats_by_bias_window: dict[tuple[float, str], WindowStats],
    max_distance_um: float | None,
    sentaurus_method_a_over_b: float | None,
    used_node_compare: Path | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    qg_ok = internal_summary_error is None or internal_summary_error <= 1.0e-10
    minus20 = -20.0
    full = stats_by_bias_window.get((minus20, "full"))
    left = stats_by_bias_window.get((minus20, "left_shoulder"))
    right = stats_by_bias_window.get((minus20, "right_shoulder"))

    def consistency_text(values: list[float]) -> str:
        if not values:
            return "insufficient_data"
        med = median(values)
        return "yes" if med <= math.log10(2.0) else "no"

    f_consistency = consistency_text(full.ratio_logs_f if full else [])
    alpha_consistency = consistency_text(full.ratio_logs_alpha if full else [])
    j_consistency = consistency_text(full.ratio_logs_j if full else [])
    main_factor = dominant_factor(full) if full and full.rows else "insufficient_data"
    left_factor = dominant_factor(left) if left and left.rows else "insufficient_data"
    right_factor = dominant_factor(right) if right and right.rows else "insufficient_data"
    def qg_direction(stats: WindowStats | None) -> tuple[str, str]:
        if stats is None or not stats.rows:
            return "insufficient_data", ""
        ratio = safe_ratio(stats.qg_self, stats.qg_sentaurus_a)
        if ratio is None:
            return "insufficient_data", ""
        if ratio > 1.0:
            return "overstrong", f"self/Sentaurus={ratio:.17g}"
        return "low", f"self/Sentaurus={ratio:.17g}"

    left_direction, left_ratio_text = qg_direction(left)
    right_direction, right_ratio_text = qg_direction(right)
    nodal_ok = (
        sentaurus_method_a_over_b is not None and
        0.5 <= sentaurus_method_a_over_b <= 2.0
    )

    with path.open("w", encoding="utf-8", newline="\n") as output:
        output.write("# Internal Edge Source vs Sentaurus Comparison\n\n")
        output.write("- model: VanOverstraeten/de Man\n")
        output.write("- driving_force: GradQuasiFermi\n")
        output.write("- parameter_set: default\n")
        output.write("- A_scale: 1\n")
        output.write("- B_scale: 1\n")
        output.write("- self current source: internal edge source audit rows\n")
        output.write("- self exported node current density used: no\n")
        output.write("- current unit: A/cm^2\n")
        output.write("- field unit: V/cm\n")
        output.write("- alpha unit: cm^-1\n")
        output.write("- Gava unit: cm^-3 s^-1\n")
        output.write("- qG contribution unit: A/um\n")
        if used_node_compare is not None:
            output.write(f"- node_compare_reference: {used_node_compare}\n")
        if internal_summary_error is not None:
            output.write(f"- internal_qG_vs_solver_relative_error: {internal_summary_error:.17g}\n")
        if sentaurus_method_a_over_b is not None:
            output.write(
                "- sentaurus_method_A_qG_over_method_B_qG_full_minus20V: "
                f"{sentaurus_method_a_over_b:.17g}\n"
            )
        if max_distance_um is not None:
            output.write(f"- max_Gava_edge_distance_minus20V_um: {max_distance_um:.17g}\n")

        output.write("\n## Window Integrals\n\n")
        output.write("| bias_V | window | rows | qG_self_A_per_um | qG_sentaurus_method_A_A_per_um | qG_sentaurus_method_B_A_per_um | self_over_sentaurus_A | self_over_sentaurus_B | dominant_factor |\n")
        output.write("| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for (bias, window), stats in sorted(stats_by_bias_window.items()):
            ratio_a = safe_ratio(stats.qg_self, stats.qg_sentaurus_a)
            ratio_b = safe_ratio(stats.qg_self, stats.qg_sentaurus_b)
            output.write(
                f"| {bias:.17g} | {window} | {stats.rows} | {stats.qg_self:.17g} | "
                f"{stats.qg_sentaurus_a:.17g} | {stats.qg_sentaurus_b:.17g} | "
                f"{ratio_text(ratio_a)} | {ratio_text(ratio_b)} | {dominant_factor(stats)} |\n"
            )

        output.write("\n## Required Answers\n\n")
        output.write(f"1. self internal edge qG explains solver qG: {'yes' if qg_ok else 'no'}.\n")
        output.write(f"2. self edge Fn/Fp matches Sentaurus edge-equivalent Fn/Fp: {f_consistency}.\n")
        output.write(f"3. self edge alpha matches Sentaurus edge-equivalent alpha: {alpha_consistency}.\n")
        output.write(f"4. self edge Jn/Jp matches Sentaurus edge-equivalent Jn/Jp: {j_consistency}.\n")
        output.write(f"5. Main Gava difference factor at -20V full window: {main_factor}.\n")
        output.write(
            "6. -20V left_shoulder is "
            f"{left_direction} ({left_ratio_text}), dominant factor: {left_factor}; "
            "right_shoulder is "
            f"{right_direction} ({right_ratio_text}), dominant factor: {right_factor}.\n"
        )
        output.write(
            "7. max Gava edge distance at -20V: "
            f"{'insufficient_data' if max_distance_um is None else format(max_distance_um, '.17g')} um.\n"
        )
        if nodal_ok:
            output.write("8. Sentaurus nodal averaging can explain Sentaurus nodal Gava at edge level: yes.\n")
        else:
            output.write(
                "8. Sentaurus nodal averaging cannot explain its Gava; use cell/box-face reconstruction rather than edge nodal averaging.\n"
            )
        if j_consistency == "no":
            output.write("\nNext step: self edge J differs most; investigate edge current formula, SG current, and quasi-Fermi sign convention.\n")
        elif alpha_consistency == "no":
            output.write("\nNext step: self edge alpha differs most; investigate alpha location and VanOverstraeten branch.\n")
        elif main_factor == "source_area_or_edge_mapping":
            output.write("\nNext step: F/alpha/J are close but qG window is not; investigate source area and edge-to-control-volume weighting.\n")


def compare(args: argparse.Namespace) -> None:
    requested_biases = parse_biases(args.biases)
    topology = load_topology(args)
    if not topology:
        raise SystemExit("No edge topology available; pass --self-edge-topology or --elements with --nodes/--mesh-json.")
    sentaurus_by_bias = load_sentaurus_bias_fields(args.sentaurus_multibias_dir, requested_biases)
    if not sentaurus_by_bias:
        raise SystemExit("No Sentaurus multibias field directories were found.")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    stats_by_bias_window: dict[tuple[float, str], WindowStats] = defaultdict(WindowStats)
    self_max = (-math.inf, None)
    sentaurus_max = (-math.inf, None)
    sentaurus_a_qg_full_minus20 = 0.0
    sentaurus_b_qg_full_minus20 = 0.0

    header = [
        "bias_V",
        "edge_id",
        "node0",
        "node1",
        "x_mid_um",
        "y_mid_um",
        "edge_length_um",
        "source_area_cm2",
        "Fn_self_edge_V_per_cm",
        "Fp_self_edge_V_per_cm",
        "alpha_n_self_edge_cm_inv",
        "alpha_p_self_edge_cm_inv",
        "Jn_self_edge_A_per_cm2",
        "Jp_self_edge_A_per_cm2",
        "Gava_n_self_edge_cm_minus3_s_minus1",
        "Gava_p_self_edge_cm_minus3_s_minus1",
        "Gava_self_edge_cm_minus3_s_minus1",
        "qG_self_edge_A_per_um",
        "Fn_sentaurus_edge_V_per_cm",
        "Fp_sentaurus_edge_V_per_cm",
        "alpha_n_sentaurus_edge_cm_inv",
        "alpha_p_sentaurus_edge_cm_inv",
        "Jn_sentaurus_edge_A_per_cm2",
        "Jp_sentaurus_edge_A_per_cm2",
        "Gava_sentaurus_edge_method_A_cm_minus3_s_minus1",
        "Gava_sentaurus_edge_method_B_cm_minus3_s_minus1",
        "qG_sentaurus_edge_method_A_A_per_um",
        "qG_sentaurus_edge_method_B_A_per_um",
        "Fn_ratio_self_over_sentaurus",
        "Fp_ratio_self_over_sentaurus",
        "alpha_n_ratio_self_over_sentaurus",
        "alpha_p_ratio_self_over_sentaurus",
        "Jn_ratio_self_over_sentaurus",
        "Jp_ratio_self_over_sentaurus",
        "Gava_ratio_self_over_sentaurus_method_A",
        "Gava_ratio_self_over_sentaurus_method_B",
        "qG_ratio_self_over_sentaurus_method_A",
        "qG_ratio_self_over_sentaurus_method_B",
    ]

    with args.out_csv.open("w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=header)
        writer.writeheader()
        with args.internal_audit.open(newline="", encoding="utf-8") as in_handle:
            for raw in csv.DictReader(in_handle):
                row = clean_row(raw)
                if row.get("source_location_type", "edge") != "edge":
                    continue
                raw_bias = parse_float(row["bias_V"])
                bias = nearest_bias(raw_bias, requested_biases, args.bias_tol)
                if bias is None or bias not in sentaurus_by_bias:
                    continue
                edge_id = int(row["source_entity_id"])
                edge = topology.get(edge_id)
                if edge is None:
                    continue
                fields = sentaurus_by_bias[bias]
                try:
                    sent = sentaurus_edge_terms(edge, fields)
                except KeyError:
                    continue
                source_area_cm2 = parse_float(
                    row.get("source_weight_or_volume_cm2_for_2D")
                    or row.get("contribution_volume_cm3_or_area_cm2_for_2D")
                    or "0"
                )
                if source_area_cm2 > 0.0 and edge.source_area_cm2 <= 0.0:
                    edge.source_area_cm2 = source_area_cm2
                    sent = sentaurus_edge_terms(edge, fields)

                self_values = {
                    "Fn_self_edge_V_per_cm": parse_float(row["Fn_used_V_per_cm"]),
                    "Fp_self_edge_V_per_cm": parse_float(row["Fp_used_V_per_cm"]),
                    "alpha_n_self_edge_cm_inv": parse_float(row["alpha_n_used_cm_inv"]),
                    "alpha_p_self_edge_cm_inv": parse_float(row["alpha_p_used_cm_inv"]),
                    "Jn_self_edge_A_per_cm2": parse_float(row["Jn_mag_used_A_per_cm2"]),
                    "Jp_self_edge_A_per_cm2": parse_float(row["Jp_mag_used_A_per_cm2"]),
                    "Gava_n_self_edge_cm_minus3_s_minus1": parse_float(row["Gava_n_used_cm_minus3_s_minus1"]),
                    "Gava_p_self_edge_cm_minus3_s_minus1": parse_float(row["Gava_p_used_cm_minus3_s_minus1"]),
                    "Gava_self_edge_cm_minus3_s_minus1": parse_float(row["Gava_total_used_cm_minus3_s_minus1"]),
                    "qG_self_edge_A_per_um": parse_float(row["qG_contribution_A_per_um"]),
                }
                ratios = {
                    "Fn_ratio_self_over_sentaurus": safe_ratio(self_values["Fn_self_edge_V_per_cm"], sent["Fn_sentaurus_edge_V_per_cm"]),
                    "Fp_ratio_self_over_sentaurus": safe_ratio(self_values["Fp_self_edge_V_per_cm"], sent["Fp_sentaurus_edge_V_per_cm"]),
                    "alpha_n_ratio_self_over_sentaurus": safe_ratio(self_values["alpha_n_self_edge_cm_inv"], sent["alpha_n_sentaurus_edge_cm_inv"]),
                    "alpha_p_ratio_self_over_sentaurus": safe_ratio(self_values["alpha_p_self_edge_cm_inv"], sent["alpha_p_sentaurus_edge_cm_inv"]),
                    "Jn_ratio_self_over_sentaurus": safe_ratio(self_values["Jn_self_edge_A_per_cm2"], sent["Jn_sentaurus_edge_A_per_cm2"]),
                    "Jp_ratio_self_over_sentaurus": safe_ratio(self_values["Jp_self_edge_A_per_cm2"], sent["Jp_sentaurus_edge_A_per_cm2"]),
                    "Gava_ratio_self_over_sentaurus_method_A": safe_ratio(self_values["Gava_self_edge_cm_minus3_s_minus1"], sent["Gava_sentaurus_edge_method_A_cm_minus3_s_minus1"]),
                    "Gava_ratio_self_over_sentaurus_method_B": safe_ratio(self_values["Gava_self_edge_cm_minus3_s_minus1"], sent["Gava_sentaurus_edge_method_B_cm_minus3_s_minus1"]),
                    "qG_ratio_self_over_sentaurus_method_A": safe_ratio(self_values["qG_self_edge_A_per_um"], sent["qG_sentaurus_edge_method_A_A_per_um"]),
                    "qG_ratio_self_over_sentaurus_method_B": safe_ratio(self_values["qG_self_edge_A_per_um"], sent["qG_sentaurus_edge_method_B_A_per_um"]),
                }

                output_row = {
                    "bias_V": f"{bias:.17g}",
                    "edge_id": str(edge.edge_id),
                    "node0": str(edge.node0),
                    "node1": str(edge.node1),
                    "x_mid_um": f"{edge.x_mid_um:.17g}",
                    "y_mid_um": f"{edge.y_mid_um:.17g}",
                    "edge_length_um": f"{edge.edge_length_um:.17g}",
                    "source_area_cm2": f"{source_area_cm2:.17g}",
                }
                output_row.update({key: f"{value:.17g}" for key, value in self_values.items()})
                output_row.update({key: f"{value:.17g}" for key, value in sent.items()})
                output_row.update({key: ratio_text(value) for key, value in ratios.items()})
                writer.writerow(output_row)

                for window in window_for_x(edge.x_mid_um):
                    stats = stats_by_bias_window[(bias, window)]
                    stats.rows += 1
                    stats.qg_self += self_values["qG_self_edge_A_per_um"]
                    stats.qg_sentaurus_a += sent["qG_sentaurus_edge_method_A_A_per_um"]
                    stats.qg_sentaurus_b += sent["qG_sentaurus_edge_method_B_A_per_um"]
                    stats.g_self_max = max(stats.g_self_max, self_values["Gava_self_edge_cm_minus3_s_minus1"])
                    stats.g_sentaurus_b_max = max(stats.g_sentaurus_b_max, sent["Gava_sentaurus_edge_method_B_cm_minus3_s_minus1"])
                    add_ratio_log(stats.ratio_logs_f, ratios["Fn_ratio_self_over_sentaurus"])
                    add_ratio_log(stats.ratio_logs_f, ratios["Fp_ratio_self_over_sentaurus"])
                    add_ratio_log(stats.ratio_logs_alpha, ratios["alpha_n_ratio_self_over_sentaurus"])
                    add_ratio_log(stats.ratio_logs_alpha, ratios["alpha_p_ratio_self_over_sentaurus"])
                    add_ratio_log(stats.ratio_logs_j, ratios["Jn_ratio_self_over_sentaurus"])
                    add_ratio_log(stats.ratio_logs_j, ratios["Jp_ratio_self_over_sentaurus"])
                    add_ratio_log(stats.ratio_logs_g, ratios["Gava_ratio_self_over_sentaurus_method_A"])
                if abs(bias + 20.0) <= args.bias_tol:
                    if self_values["Gava_self_edge_cm_minus3_s_minus1"] > self_max[0]:
                        self_max = (self_values["Gava_self_edge_cm_minus3_s_minus1"], (edge.x_mid_um, edge.y_mid_um))
                    if sent["Gava_sentaurus_edge_method_B_cm_minus3_s_minus1"] > sentaurus_max[0]:
                        sentaurus_max = (sent["Gava_sentaurus_edge_method_B_cm_minus3_s_minus1"], (edge.x_mid_um, edge.y_mid_um))
                    sentaurus_a_qg_full_minus20 += sent["qG_sentaurus_edge_method_A_A_per_um"]
                    sentaurus_b_qg_full_minus20 += sent["qG_sentaurus_edge_method_B_A_per_um"]

    max_distance_um = None
    if self_max[1] is not None and sentaurus_max[1] is not None:
        max_distance_um = math.hypot(self_max[1][0] - sentaurus_max[1][0], self_max[1][1] - sentaurus_max[1][1])
    method_a_over_b = safe_ratio(sentaurus_a_qg_full_minus20, sentaurus_b_qg_full_minus20)
    write_summary(
        args.out_summary,
        internal_summary_error=parse_internal_qg_relative_error(args.internal_summary),
        stats_by_bias_window=stats_by_bias_window,
        max_distance_um=max_distance_um,
        sentaurus_method_a_over_b=method_a_over_b,
        used_node_compare=args.node_compare if args.node_compare.exists() else None,
    )


def main() -> int:
    args = parse_args()
    compare(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
