#!/usr/bin/env python3
"""Decompose internal edge current support against Sentaurus edge-equivalent current.

This diagnostic is intentionally post-processing only.  The self current used
here comes from the internal edge source audit rows, not from exported self node
current-density fields.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19
RATIO_FLOOR = 1.0e-300

DEFAULT_NODE_COMPARE = (
    REPO
    / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
    / "coarse_previous_full20_fields_20260629_155228/coarse_node_field_compare_aligned.csv"
)
DEFAULT_INTERNAL_AUDIT = REPO / "build/diagnostics/avalanche_internal_source_current_audit.csv"
DEFAULT_SG_EDGE_TOPOLOGY = (
    REPO
    / "build/diagnostics/vanoverstraeten_A2_B105_spatial_validation_case"
    / "BV-A2-B1p05-spatial/BV-A2-B1p05-spatial_sg_avalanche_edges.csv"
)
DEFAULT_OUT_CSV = REPO / "build/diagnostics/edge_current_decomposition_audit.csv"
DEFAULT_SUMMARY = REPO / "build/diagnostics/edge_current_decomposition_audit_summary.md"

WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}


@dataclass(frozen=True)
class EdgeTopology:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    edge_length_um: float

    @property
    def x_mid_um(self) -> float:
        return 0.5 * (self.x0_um + self.x1_um)

    @property
    def y_mid_um(self) -> float:
        return 0.5 * (self.y0_um + self.y1_um)


@dataclass
class DecompStats:
    rows: int = 0
    jn_ratios: list[float] = None
    jp_ratios: list[float] = None
    jn_recon_ratios: list[float] = None
    jp_recon_ratios: list[float] = None
    jn_self_closure: list[float] = None
    jp_self_closure: list[float] = None
    jn_sent_closure: list[float] = None
    jp_sent_closure: list[float] = None
    fn_ratios: list[float] = None
    fp_ratios: list[float] = None
    n_ratios: list[float] = None
    p_ratios: list[float] = None
    mu_n_ratios: list[float] = None
    mu_p_ratios: list[float] = None
    jn_factors: Counter[str] = None
    jp_factors: Counter[str] = None

    def __post_init__(self) -> None:
        self.jn_ratios = []
        self.jp_ratios = []
        self.jn_recon_ratios = []
        self.jp_recon_ratios = []
        self.jn_self_closure = []
        self.jp_self_closure = []
        self.jn_sent_closure = []
        self.jp_sent_closure = []
        self.fn_ratios = []
        self.fp_ratios = []
        self.n_ratios = []
        self.p_ratios = []
        self.mu_n_ratios = []
        self.mu_p_ratios = []
        self.jn_factors = Counter()
        self.jp_factors = Counter()


HEADER = [
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "x_mid_um",
    "y_mid_um",
    "edge_length_um",
    "Fn_self_edge_V_per_cm",
    "Fp_self_edge_V_per_cm",
    "n_self_edge_cm_minus3",
    "p_self_edge_cm_minus3",
    "mu_n_self_edge_cm2_per_Vs",
    "mu_p_self_edge_cm2_per_Vs",
    "Jn_self_used_A_per_cm2",
    "Jp_self_used_A_per_cm2",
    "Jn_self_reconstructed_q_mu_n_F_A_per_cm2",
    "Jp_self_reconstructed_q_mu_p_F_A_per_cm2",
    "Jn_self_closure_ratio",
    "Jp_self_closure_ratio",
    "Fn_sent_edge_V_per_cm",
    "Fp_sent_edge_V_per_cm",
    "n_sent_edge_cm_minus3",
    "p_sent_edge_cm_minus3",
    "mu_n_sent_edge_cm2_per_Vs",
    "mu_p_sent_edge_cm2_per_Vs",
    "Jn_sent_edge_A_per_cm2",
    "Jp_sent_edge_A_per_cm2",
    "Jn_sent_reconstructed_q_mu_n_F_A_per_cm2",
    "Jp_sent_reconstructed_q_mu_p_F_A_per_cm2",
    "Fn_ratio_self_over_sent",
    "Fp_ratio_self_over_sent",
    "n_ratio_self_over_sent",
    "p_ratio_self_over_sent",
    "mu_n_ratio_self_over_sent",
    "mu_p_ratio_self_over_sent",
    "Jn_ratio_self_over_sent",
    "Jp_ratio_self_over_sent",
    "Jn_reconstructed_ratio_self_over_sent",
    "Jp_reconstructed_ratio_self_over_sent",
    "dominant_Jn_difference_factor",
    "dominant_Jp_difference_factor",
    "window_label",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-audit", type=Path, default=DEFAULT_INTERNAL_AUDIT)
    parser.add_argument("--self-edge-topology", type=Path, default=DEFAULT_SG_EDGE_TOPOLOGY)
    parser.add_argument("--node-compare", type=Path, default=DEFAULT_NODE_COMPARE)
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--out-csv", type=Path, default=DEFAULT_OUT_CSV)
    parser.add_argument("--out-summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items() if key is not None}


def as_float(raw: str | None, default: float | None = None) -> float | None:
    if raw is None:
        return default
    text = str(raw).strip()
    if text == "":
        return default
    try:
        value = float(text)
    except ValueError:
        return default
    return value if math.isfinite(value) else default


def parse_float(raw: str | None, default: float = 0.0) -> float:
    value = as_float(raw, default)
    return float(default if value is None else value)


def parse_biases(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def nearest_bias(value: float, targets: list[float], tol: float) -> float | None:
    best = min(targets, key=lambda bias: abs(value - bias))
    return best if abs(value - best) <= tol else None


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None:
        return None
    if abs(denominator) <= RATIO_FLOOR or not math.isfinite(denominator):
        return None
    ratio = numerator / denominator
    return ratio if math.isfinite(ratio) else None


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.17g}"


def median_or_none(values: list[float]) -> float | None:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    return median(finite) if finite else None


def add_ratio(bucket: list[float], ratio: float | None) -> None:
    if ratio is not None and ratio > 0.0 and math.isfinite(ratio):
        bucket.append(ratio)


def load_topology_from_sg_edges(path: Path) -> dict[int, EdgeTopology]:
    topology: dict[int, EdgeTopology] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            edge_id = int(float(row["edge_id"]))
            if edge_id in topology:
                continue
            topology[edge_id] = EdgeTopology(
                edge_id=edge_id,
                node0=int(float(row["node0"])),
                node1=int(float(row["node1"])),
                x0_um=parse_float(row["x0_um"]),
                y0_um=parse_float(row["y0_um"]),
                x1_um=parse_float(row["x1_um"]),
                y1_um=parse_float(row["y1_um"]),
                edge_length_um=parse_float(row["edge_length_m"]) * 1.0e6,
            )
    return topology


def read_node_compare(path: Path) -> dict[tuple[float, int], dict[str, Any]]:
    nodes: dict[tuple[float, int], dict[str, Any]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            bias = as_float(row.get("bias_V"))
            node_id_raw = row.get("node_id")
            quantity = row.get("quantity", "")
            if bias is None or node_id_raw is None or quantity == "":
                continue
            node_id = int(float(node_id_raw))
            key = (bias, node_id)
            entry = nodes.setdefault(key, {
                "bias_V": bias,
                "node_id": node_id,
                "x_um": as_float(row.get("x_um")),
                "y_um": as_float(row.get("y_um")),
                "values": {},
            })
            entry["values"][quantity] = {
                "sent": as_float(row.get("sentaurus_value")),
                "self": as_float(row.get("vela_value_scaled_to_sentaurus_units")),
            }
    return nodes


def node_value(nodes: dict[tuple[float, int], dict[str, Any]], bias: float, node_id: int, quantity: str, side: str) -> float | None:
    node = nodes.get((bias, node_id))
    if node is None:
        return None
    item = node["values"].get(quantity)
    if item is None:
        return None
    return item.get(side)


def edge_average(
    nodes: dict[tuple[float, int], dict[str, Any]],
    bias: float,
    edge: EdgeTopology,
    quantity: str,
    side: str,
) -> float | None:
    a = node_value(nodes, bias, edge.node0, quantity, side)
    b = node_value(nodes, bias, edge.node1, quantity, side)
    if a is None or b is None:
        return None
    return 0.5 * (a + b)


def edge_qf_gradient(
    nodes: dict[tuple[float, int], dict[str, Any]],
    bias: float,
    edge: EdgeTopology,
    quantity: str,
    side: str,
) -> float | None:
    a = node_value(nodes, bias, edge.node0, quantity, side)
    b = node_value(nodes, bias, edge.node1, quantity, side)
    length_cm = edge.edge_length_um * 1.0e-4
    if a is None or b is None or length_cm <= 0.0:
        return None
    return abs(b - a) / length_cm


def sentaurus_edge_current(
    nodes: dict[tuple[float, int], dict[str, Any]],
    bias: float,
    edge: EdgeTopology,
    prefix: str,
) -> float | None:
    jx0 = node_value(nodes, bias, edge.node0, f"{prefix}_x", "sent")
    jy0 = node_value(nodes, bias, edge.node0, f"{prefix}_y", "sent")
    jx1 = node_value(nodes, bias, edge.node1, f"{prefix}_x", "sent")
    jy1 = node_value(nodes, bias, edge.node1, f"{prefix}_y", "sent")
    if None in (jx0, jy0, jx1, jy1):
        mag0 = node_value(nodes, bias, edge.node0, f"{prefix}_mag", "sent")
        mag1 = node_value(nodes, bias, edge.node1, f"{prefix}_mag", "sent")
        if mag0 is None or mag1 is None:
            return None
        return 0.5 * (mag0 + mag1)
    return math.hypot(0.5 * (jx0 + jx1), 0.5 * (jy0 + jy1))


def reconstructed_current(mu_cm2_per_vs: float | None, density_cm3: float | None, field_v_cm: float | None) -> float | None:
    if mu_cm2_per_vs is None or density_cm3 is None or field_v_cm is None:
        return None
    return Q_C * abs(mu_cm2_per_vs) * abs(density_cm3) * abs(field_v_cm)


def window_label(x_um: float) -> str:
    for name in ("left_shoulder", "center", "right_shoulder"):
        lo, hi = WINDOWS[name]
        if lo <= x_um <= hi:
            return name
    lo, hi = WINDOWS["junction"]
    if lo <= x_um <= hi:
        return "junction"
    return "outside_junction"


def window_members(x_um: float) -> list[str]:
    return [name for name, (lo, hi) in WINDOWS.items() if lo <= x_um <= hi]


def log_error(ratio: float | None) -> float:
    if ratio is None or ratio <= 0.0 or not math.isfinite(ratio):
        return 0.0
    return abs(math.log10(ratio))


def dominant_factor(
    *,
    actual_ratio: float | None,
    reconstructed_ratio: float | None,
    self_closure: float | None,
    sent_closure: float | None,
    field_ratio: float | None,
    density_ratio: float | None,
    mobility_ratio: float | None,
) -> str:
    if self_closure is not None and log_error(self_closure) > math.log10(2.0):
        return "SG_formula_or_density_averaging"
    if sent_closure is not None and log_error(sent_closure) > math.log10(2.0):
        return "sentaurus_q_mu_n_F_not_closing_node_current"
    if actual_ratio is not None and reconstructed_ratio is not None:
        residual = safe_ratio(actual_ratio, reconstructed_ratio)
        if residual is not None and log_error(residual) > math.log10(2.0):
            return "unit_conversion_or_sign_convention"
    factors = {
        "quasi_fermi_gradient": field_ratio,
        "carrier_density": density_ratio,
        "mobility": mobility_ratio,
    }
    return max(factors.items(), key=lambda item: log_error(item[1]))[0]


def write_summary(path: Path, rows: list[dict[str, Any]], stats: dict[tuple[float, str], DecompStats]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    full_stats = stats.get((-20.0, "full"), DecompStats())

    def closure_answer(values: list[float]) -> str:
        med = median_or_none([log_error(value) for value in values])
        if med is None:
            return "insufficient_data"
        return "yes" if med <= math.log10(2.0) else "no"

    self_closes = closure_answer(full_stats.jn_self_closure + full_stats.jp_self_closure)
    sent_closes = closure_answer(full_stats.jn_sent_closure + full_stats.jp_sent_closure)
    jn_error = median_or_none([log_error(value) for value in full_stats.jn_ratios]) or 0.0
    jp_error = median_or_none([log_error(value) for value in full_stats.jp_ratios]) or 0.0
    worse_carrier = "electron" if jn_error >= jp_error else "hole"

    factor_counter = full_stats.jn_factors + full_stats.jp_factors
    main_factor = factor_counter.most_common(1)[0][0] if factor_counter else "insufficient_data"

    def med_ratio(values: list[float]) -> str:
        value = median_or_none(values)
        return "insufficient_data" if value is None else f"{value:.6g}"

    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Edge Current Decomposition Audit\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- self current source: internal edge source rows\n")
        out.write("- self exported node current density used: no\n")
        out.write("- Sentaurus edge current: nodal vector average magnitude\n")
        out.write("- F unit: V/cm\n")
        out.write("- n,p unit: cm^-3\n")
        out.write("- mobility unit: cm^2/V/s\n")
        out.write("- J unit: A/cm^2\n")
        out.write(f"- rows: {len(rows)}\n\n")

        out.write("## Window Median Ratios\n\n")
        out.write("| bias_V | window | rows | Jn self/Sent | Jp self/Sent | Jn qmuF self/Sent | Jp qmuF self/Sent | dominant Jn | dominant Jp |\n")
        out.write("| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for (bias, window), item in sorted(stats.items()):
            jn_factor = item.jn_factors.most_common(1)[0][0] if item.jn_factors else "insufficient_data"
            jp_factor = item.jp_factors.most_common(1)[0][0] if item.jp_factors else "insufficient_data"
            out.write(
                f"| {bias:.17g} | {window} | {item.rows} | {med_ratio(item.jn_ratios)} | "
                f"{med_ratio(item.jp_ratios)} | {med_ratio(item.jn_recon_ratios)} | "
                f"{med_ratio(item.jp_recon_ratios)} | {jn_factor} | {jp_factor} |\n"
            )

        out.write("\n## Required Answers\n\n")
        out.write(f"1. self Jn/Jp can be reconstructed from self q*mu*n*F: {self_closes}.\n")
        out.write(f"2. Sentaurus edge Jn/Jp can be reconstructed from Sentaurus q*mu*n*F: {sent_closes}.\n")
        out.write(
            "3. Main self/Sentaurus J difference source at -20V full window: "
            f"{main_factor}; median Fn ratio={med_ratio(full_stats.fn_ratios)}, "
            f"Fp ratio={med_ratio(full_stats.fp_ratios)}, n ratio={med_ratio(full_stats.n_ratios)}, "
            f"p ratio={med_ratio(full_stats.p_ratios)}, mu_n ratio={med_ratio(full_stats.mu_n_ratios)}, "
            f"mu_p ratio={med_ratio(full_stats.mu_p_ratios)}.\n"
        )
        out.write("4. -20V full/junction/center/shoulder low-J main cause: see Window Median Ratios dominant columns; ")
        out.write(f"full-window dominant factor is {main_factor}.\n")
        out.write(f"5. Carrier with larger J discrepancy at -20V full window: {worse_carrier}.\n")
        if sent_closes == "yes" and self_closes == "no":
            out.write("6. Sentaurus q*mu*n*F closes while self does not; investigate self SG current formula first.\n")
        else:
            out.write("6. Sentaurus/self closure does not isolate a self-only SG formula failure by itself.\n")
        if self_closes == "yes" and main_factor in {"carrier_density", "mobility", "quasi_fermi_gradient"}:
            out.write(
                "7. self q*mu*n*F and self J agree; current is low because the listed n/p, mobility, or qF-gradient factor is low.\n"
            )
        else:
            out.write("7. If self q*mu*n*F later closes but J remains low, inspect which n/p/mu/F term is low.\n")
        if main_factor in {"unit_conversion_or_sign_convention", "SG_formula_or_density_averaging"}:
            out.write(
                "8. n/p/mu/F do not fully explain J; inspect SG Bernoulli/current formula, edge orientation, and sign convention.\n"
            )
        else:
            out.write(
                "8. If future runs show n/p/mu/F close but J low, inspect SG Bernoulli/current formula, edge orientation, and sign convention.\n"
            )


def run(args: argparse.Namespace) -> None:
    requested_biases = parse_biases(args.biases)
    if not args.internal_audit.exists():
        raise SystemExit(f"Internal audit CSV not found: {args.internal_audit}")
    if not args.self_edge_topology.exists():
        raise SystemExit(f"Edge topology CSV not found: {args.self_edge_topology}")
    if not args.node_compare.exists():
        raise SystemExit(f"Node compare CSV not found: {args.node_compare}")

    topology = load_topology_from_sg_edges(args.self_edge_topology)
    nodes = read_node_compare(args.node_compare)
    rows_out: list[dict[str, Any]] = []
    stats: dict[tuple[float, str], DecompStats] = defaultdict(DecompStats)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.internal_audit.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            if row.get("source_location_type", "edge") != "edge":
                continue
            raw_bias = parse_float(row.get("bias_V"))
            bias = nearest_bias(raw_bias, requested_biases, args.bias_tol)
            if bias is None:
                continue
            edge_id = int(float(row["source_entity_id"]))
            edge = topology.get(edge_id)
            if edge is None:
                continue

            fn_self = parse_float(row.get("Fn_used_V_per_cm"))
            fp_self = parse_float(row.get("Fp_used_V_per_cm"))
            jn_self = parse_float(row.get("Jn_mag_used_A_per_cm2"))
            jp_self = parse_float(row.get("Jp_mag_used_A_per_cm2"))

            n_self = edge_average(nodes, bias, edge, "electron_density", "self")
            p_self = edge_average(nodes, bias, edge, "hole_density", "self")
            mu_n_self = edge_average(nodes, bias, edge, "electron_mobility", "self")
            mu_p_self = edge_average(nodes, bias, edge, "hole_mobility", "self")
            jn_self_recon = reconstructed_current(mu_n_self, n_self, fn_self)
            jp_self_recon = reconstructed_current(mu_p_self, p_self, fp_self)
            jn_self_closure = safe_ratio(jn_self, jn_self_recon)
            jp_self_closure = safe_ratio(jp_self, jp_self_recon)

            fn_sent = edge_qf_gradient(nodes, bias, edge, "electron_qf", "sent")
            fp_sent = edge_qf_gradient(nodes, bias, edge, "hole_qf", "sent")
            n_sent = edge_average(nodes, bias, edge, "electron_density", "sent")
            p_sent = edge_average(nodes, bias, edge, "hole_density", "sent")
            mu_n_sent = edge_average(nodes, bias, edge, "electron_mobility", "sent")
            mu_p_sent = edge_average(nodes, bias, edge, "hole_mobility", "sent")
            jn_sent = sentaurus_edge_current(nodes, bias, edge, "electron_current_density")
            jp_sent = sentaurus_edge_current(nodes, bias, edge, "hole_current_density")
            jn_sent_recon = reconstructed_current(mu_n_sent, n_sent, fn_sent)
            jp_sent_recon = reconstructed_current(mu_p_sent, p_sent, fp_sent)
            jn_sent_closure = safe_ratio(jn_sent, jn_sent_recon)
            jp_sent_closure = safe_ratio(jp_sent, jp_sent_recon)

            ratios = {
                "Fn_ratio_self_over_sent": safe_ratio(fn_self, fn_sent),
                "Fp_ratio_self_over_sent": safe_ratio(fp_self, fp_sent),
                "n_ratio_self_over_sent": safe_ratio(n_self, n_sent),
                "p_ratio_self_over_sent": safe_ratio(p_self, p_sent),
                "mu_n_ratio_self_over_sent": safe_ratio(mu_n_self, mu_n_sent),
                "mu_p_ratio_self_over_sent": safe_ratio(mu_p_self, mu_p_sent),
                "Jn_ratio_self_over_sent": safe_ratio(jn_self, jn_sent),
                "Jp_ratio_self_over_sent": safe_ratio(jp_self, jp_sent),
                "Jn_reconstructed_ratio_self_over_sent": safe_ratio(jn_self_recon, jn_sent_recon),
                "Jp_reconstructed_ratio_self_over_sent": safe_ratio(jp_self_recon, jp_sent_recon),
            }
            dominant_jn = dominant_factor(
                actual_ratio=ratios["Jn_ratio_self_over_sent"],
                reconstructed_ratio=ratios["Jn_reconstructed_ratio_self_over_sent"],
                self_closure=jn_self_closure,
                sent_closure=jn_sent_closure,
                field_ratio=ratios["Fn_ratio_self_over_sent"],
                density_ratio=ratios["n_ratio_self_over_sent"],
                mobility_ratio=ratios["mu_n_ratio_self_over_sent"],
            )
            dominant_jp = dominant_factor(
                actual_ratio=ratios["Jp_ratio_self_over_sent"],
                reconstructed_ratio=ratios["Jp_reconstructed_ratio_self_over_sent"],
                self_closure=jp_self_closure,
                sent_closure=jp_sent_closure,
                field_ratio=ratios["Fp_ratio_self_over_sent"],
                density_ratio=ratios["p_ratio_self_over_sent"],
                mobility_ratio=ratios["mu_p_ratio_self_over_sent"],
            )

            out_row: dict[str, Any] = {
                "bias_V": bias,
                "edge_id": edge.edge_id,
                "node0": edge.node0,
                "node1": edge.node1,
                "x_mid_um": edge.x_mid_um,
                "y_mid_um": edge.y_mid_um,
                "edge_length_um": edge.edge_length_um,
                "Fn_self_edge_V_per_cm": fn_self,
                "Fp_self_edge_V_per_cm": fp_self,
                "n_self_edge_cm_minus3": n_self,
                "p_self_edge_cm_minus3": p_self,
                "mu_n_self_edge_cm2_per_Vs": mu_n_self,
                "mu_p_self_edge_cm2_per_Vs": mu_p_self,
                "Jn_self_used_A_per_cm2": jn_self,
                "Jp_self_used_A_per_cm2": jp_self,
                "Jn_self_reconstructed_q_mu_n_F_A_per_cm2": jn_self_recon,
                "Jp_self_reconstructed_q_mu_p_F_A_per_cm2": jp_self_recon,
                "Jn_self_closure_ratio": jn_self_closure,
                "Jp_self_closure_ratio": jp_self_closure,
                "Fn_sent_edge_V_per_cm": fn_sent,
                "Fp_sent_edge_V_per_cm": fp_sent,
                "n_sent_edge_cm_minus3": n_sent,
                "p_sent_edge_cm_minus3": p_sent,
                "mu_n_sent_edge_cm2_per_Vs": mu_n_sent,
                "mu_p_sent_edge_cm2_per_Vs": mu_p_sent,
                "Jn_sent_edge_A_per_cm2": jn_sent,
                "Jp_sent_edge_A_per_cm2": jp_sent,
                "Jn_sent_reconstructed_q_mu_n_F_A_per_cm2": jn_sent_recon,
                "Jp_sent_reconstructed_q_mu_p_F_A_per_cm2": jp_sent_recon,
                **ratios,
                "dominant_Jn_difference_factor": dominant_jn,
                "dominant_Jp_difference_factor": dominant_jp,
                "window_label": window_label(edge.x_mid_um),
            }
            rows_out.append(out_row)

            for window in window_members(edge.x_mid_um):
                item = stats[(bias, window)]
                item.rows += 1
                add_ratio(item.jn_ratios, ratios["Jn_ratio_self_over_sent"])
                add_ratio(item.jp_ratios, ratios["Jp_ratio_self_over_sent"])
                add_ratio(item.jn_recon_ratios, ratios["Jn_reconstructed_ratio_self_over_sent"])
                add_ratio(item.jp_recon_ratios, ratios["Jp_reconstructed_ratio_self_over_sent"])
                add_ratio(item.jn_self_closure, jn_self_closure)
                add_ratio(item.jp_self_closure, jp_self_closure)
                add_ratio(item.jn_sent_closure, jn_sent_closure)
                add_ratio(item.jp_sent_closure, jp_sent_closure)
                add_ratio(item.fn_ratios, ratios["Fn_ratio_self_over_sent"])
                add_ratio(item.fp_ratios, ratios["Fp_ratio_self_over_sent"])
                add_ratio(item.n_ratios, ratios["n_ratio_self_over_sent"])
                add_ratio(item.p_ratios, ratios["p_ratio_self_over_sent"])
                add_ratio(item.mu_n_ratios, ratios["mu_n_ratio_self_over_sent"])
                add_ratio(item.mu_p_ratios, ratios["mu_p_ratio_self_over_sent"])
                item.jn_factors[dominant_jn] += 1
                item.jp_factors[dominant_jp] += 1

    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for row in rows_out:
            writer.writerow({key: fmt(row.get(key)) if isinstance(row.get(key), float) or row.get(key) is None else row.get(key) for key in HEADER})

    write_summary(args.out_summary, rows_out, stats)


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
