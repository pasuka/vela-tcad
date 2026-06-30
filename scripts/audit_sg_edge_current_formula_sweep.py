#!/usr/bin/env python3
"""Sweep SG edge-current formula candidates against Sentaurus edge currents."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19
VT_300K = 0.025851999786435535
RATIO_FLOOR = 1.0e-300

DEFAULT_INTERNAL_AUDIT = REPO / "build/diagnostics/avalanche_internal_source_current_audit.csv"
DEFAULT_SG_EDGE_TOPOLOGY = (
    REPO
    / "build/diagnostics/avalanche_internal_source_current_audit_case"
    / "BV-A2-B1p05-internal-source-current-audit"
    / "BV-A2-B1p05-internal-source-current-audit_sg_avalanche_edges.csv"
)
DEFAULT_NODE_COMPARE = (
    REPO
    / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
    / "coarse_previous_full20_fields_20260629_155228/coarse_node_field_compare_aligned.csv"
)
DEFAULT_OUT_CSV = REPO / "build/diagnostics/sg_edge_current_formula_sweep.csv"
DEFAULT_SUMMARY = REPO / "build/diagnostics/sg_edge_current_formula_sweep_summary.md"

WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}

FORMULAS = [
    "self_current_impl",
    "sg_standard",
    "qf_arithmetic_density",
    "qf_geometric_density",
    "qf_logarithmic_density",
    "qf_upwind_density",
    "qf_sg_effective_density",
]
RANKED_FORMULAS = [name for name in FORMULAS if name != "self_current_impl"]


@dataclass(frozen=True)
class EdgeTopology:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    edge_length_cm: float

    @property
    def x_mid_um(self) -> float:
        return 0.5 * (self.x0_um + self.x1_um)

    @property
    def y_mid_um(self) -> float:
        return 0.5 * (self.y0_um + self.y1_um)


@dataclass
class FormulaStats:
    vector_errors: list[float]
    projected_errors: list[float]
    self_impl_errors: list[float]
    signed_projection_rel_errors: list[float]
    flipped_signed_projection_rel_errors: list[float]

    def __init__(self) -> None:
        self.vector_errors = []
        self.projected_errors = []
        self.self_impl_errors = []
        self.signed_projection_rel_errors = []
        self.flipped_signed_projection_rel_errors = []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-audit", type=Path, default=DEFAULT_INTERNAL_AUDIT)
    parser.add_argument("--self-edge-topology", type=Path, default=DEFAULT_SG_EDGE_TOPOLOGY)
    parser.add_argument("--node-compare", type=Path, default=DEFAULT_NODE_COMPARE)
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--thermal-voltage", type=float, default=VT_300K)
    parser.add_argument("--min-sent-current", type=float, default=1.0e-30)
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
    return default if value is None else value


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


def log_error(candidate: float | None, reference: float | None) -> float | None:
    ratio = safe_ratio(abs(candidate) if candidate is not None else None, abs(reference) if reference is not None else None)
    if ratio is None or ratio <= 0.0:
        return None
    return abs(math.log10(ratio))


def signed_relative_error(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None:
        return None
    scale = max(abs(reference), RATIO_FLOOR)
    value = abs(candidate - reference) / scale
    return value if math.isfinite(value) else None


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.17g}"
    return str(value)


def median_or_none(values: list[float]) -> float | None:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    return median(finite) if finite else None


def bernoulli(x: float) -> float:
    if abs(x) < 1.0e-8:
        return 1.0 - 0.5 * x + x * x / 12.0
    if x > 700.0:
        return 0.0
    if x < -700.0:
        return -x
    return x / math.expm1(x)


def logarithmic_mean(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.5 * (a + b)
    if abs(a - b) <= 1.0e-14 * max(abs(a), abs(b), 1.0):
        return 0.5 * (a + b)
    return (b - a) / (math.log(b) - math.log(a))


def geometric_mean(a: float, b: float) -> float:
    if a < 0.0 or b < 0.0:
        return 0.5 * (a + b)
    return math.sqrt(a * b)


def upwind_density(dphi: float, d0: float, d1: float) -> float:
    if dphi > 0.0:
        return d1
    return d0


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
                edge_length_cm=parse_float(row["edge_length_m"]) * 100.0,
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
            entry = nodes.setdefault((bias, node_id), {"values": {}})
            entry["values"][quantity] = {
                "sent": as_float(row.get("sentaurus_value")),
                "self": as_float(row.get("vela_value_scaled_to_sentaurus_units")),
            }
    return nodes


def node_value(nodes: dict[tuple[float, int], dict[str, Any]], bias: float, node_id: int, quantity: str) -> float | None:
    node = nodes.get((bias, node_id))
    if node is None:
        return None
    item = node["values"].get(quantity)
    if item is None:
        return None
    return item.get("sent")


def endpoint_pair(nodes: dict[tuple[float, int], dict[str, Any]], bias: float, edge: EdgeTopology, quantity: str) -> tuple[float | None, float | None]:
    return (
        node_value(nodes, bias, edge.node0, quantity),
        node_value(nodes, bias, edge.node1, quantity),
    )


def edge_vector_reference(
    nodes: dict[tuple[float, int], dict[str, Any]],
    bias: float,
    edge: EdgeTopology,
    prefix: str,
) -> tuple[float | None, float | None, float | None]:
    jx0, jx1 = endpoint_pair(nodes, bias, edge, f"{prefix}_x")
    jy0, jy1 = endpoint_pair(nodes, bias, edge, f"{prefix}_y")
    if None in (jx0, jx1, jy0, jy1):
        return None, None, None
    jx = 0.5 * (float(jx0) + float(jx1))
    jy = 0.5 * (float(jy0) + float(jy1))
    dx = edge.x1_um - edge.x0_um
    dy = edge.y1_um - edge.y0_um
    length_um = math.hypot(dx, dy)
    if length_um <= 0.0:
        return math.hypot(jx, jy), None, None
    projected_signed = (jx * dx + jy * dy) / length_um
    return math.hypot(jx, jy), abs(projected_signed), projected_signed


def qf_current(mu: float, density: float, dphi: float, length_cm: float) -> float:
    if length_cm <= 0.0:
        return 0.0
    return -Q_C * mu * density * dphi / length_cm


def standard_sg_electron(mu: float, n0: float, n1: float, psi0: float, psi1: float, vt: float, length_cm: float) -> float:
    if length_cm <= 0.0:
        return 0.0
    u = (psi1 - psi0) / vt
    return Q_C * mu * vt / length_cm * (bernoulli(u) * n1 - bernoulli(-u) * n0)


def standard_sg_hole(mu: float, p0: float, p1: float, psi0: float, psi1: float, vt: float, length_cm: float) -> float:
    if length_cm <= 0.0:
        return 0.0
    u = (psi1 - psi0) / vt
    return Q_C * mu * vt / length_cm * (bernoulli(-u) * p1 - bernoulli(u) * p0)


def sg_effective_density_current(standard_current: float, mu: float, dphi: float, length_cm: float) -> float:
    denom = Q_C * abs(mu) * abs(dphi) / length_cm if length_cm > 0.0 else 0.0
    if denom <= RATIO_FLOOR:
        return 0.0
    density_eff = abs(standard_current) / denom
    return qf_current(mu, density_eff, dphi, length_cm)


def formula_values(
    *,
    carrier: str,
    self_impl: float,
    psi0: float,
    psi1: float,
    qf0: float,
    qf1: float,
    density0: float,
    density1: float,
    mu0: float,
    mu1: float,
    length_cm: float,
    vt: float,
) -> dict[str, float]:
    mu = 0.5 * (mu0 + mu1)
    dphi = qf1 - qf0
    standard = (
        standard_sg_electron(mu, density0, density1, psi0, psi1, vt, length_cm)
        if carrier == "Jn"
        else standard_sg_hole(mu, density0, density1, psi0, psi1, vt, length_cm)
    )
    return {
        "self_current_impl": self_impl,
        "sg_standard": standard,
        "qf_arithmetic_density": qf_current(mu, 0.5 * (density0 + density1), dphi, length_cm),
        "qf_geometric_density": qf_current(mu, geometric_mean(density0, density1), dphi, length_cm),
        "qf_logarithmic_density": qf_current(mu, logarithmic_mean(density0, density1), dphi, length_cm),
        "qf_upwind_density": qf_current(mu, upwind_density(dphi, density0, density1), dphi, length_cm),
        "qf_sg_effective_density": sg_effective_density_current(standard, mu, dphi, length_cm),
    }


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


def candidate_headers() -> list[str]:
    headers: list[str] = []
    for carrier in ("Jn", "Jp"):
        for formula in FORMULAS:
            prefix = f"{carrier}_{formula}"
            headers.extend([
                f"{prefix}_A_per_cm2",
                f"{prefix}_ratio_to_sent_vector_mag",
                f"{prefix}_ratio_to_sent_projected",
                f"{prefix}_log_error_to_sent_vector_mag",
                f"{prefix}_log_error_to_sent_projected",
                f"{prefix}_signed_rel_error_to_sent_projected_signed",
                f"{prefix}_sign_flipped_rel_error_to_sent_projected_signed",
                f"{prefix}_sign_flipped_log_error_to_sent_projected_signed",
            ])
    return headers


BASE_HEADERS = [
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "x0_um",
    "y0_um",
    "x1_um",
    "y1_um",
    "x_mid_um",
    "y_mid_um",
    "edge_length_cm",
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
    "mu_n0",
    "mu_n1",
    "mu_p0",
    "mu_p1",
    "Jn_sent_edge_vector_mag_A_per_cm2",
    "Jn_sent_edge_projected_A_per_cm2",
    "Jn_sent_edge_projected_signed_A_per_cm2",
    "Jp_sent_edge_vector_mag_A_per_cm2",
    "Jp_sent_edge_projected_A_per_cm2",
    "Jp_sent_edge_projected_signed_A_per_cm2",
    "Jn_signal_label",
    "Jp_signal_label",
    "window_label",
]


def add_formula_outputs(
    out_row: dict[str, Any],
    stats: dict[tuple[float, str, str, str], FormulaStats],
    *,
    bias: float,
    windows: list[str],
    carrier: str,
    formulas: dict[str, float],
    sent_vector: float | None,
    sent_projected: float | None,
    sent_projected_signed: float | None,
    min_current: float,
) -> None:
    high_signal = sent_vector is not None and abs(sent_vector) >= min_current
    for formula, value in formulas.items():
        prefix = f"{carrier}_{formula}"
        ratio_vector = safe_ratio(abs(value), sent_vector)
        ratio_projected = safe_ratio(abs(value), sent_projected)
        err_vector = log_error(value, sent_vector)
        err_projected = log_error(value, sent_projected)
        signed_rel_error = signed_relative_error(value, sent_projected_signed)
        flipped_signed_rel_error = signed_relative_error(-value, sent_projected_signed)
        flipped_signed_err = log_error(-value, sent_projected_signed)
        out_row[f"{prefix}_A_per_cm2"] = value
        out_row[f"{prefix}_ratio_to_sent_vector_mag"] = ratio_vector
        out_row[f"{prefix}_ratio_to_sent_projected"] = ratio_projected
        out_row[f"{prefix}_log_error_to_sent_vector_mag"] = err_vector
        out_row[f"{prefix}_log_error_to_sent_projected"] = err_projected
        out_row[f"{prefix}_signed_rel_error_to_sent_projected_signed"] = signed_rel_error
        out_row[f"{prefix}_sign_flipped_rel_error_to_sent_projected_signed"] = flipped_signed_rel_error
        out_row[f"{prefix}_sign_flipped_log_error_to_sent_projected_signed"] = flipped_signed_err
        if high_signal and (formula in RANKED_FORMULAS or formula == "self_current_impl"):
            for window in windows:
                item = stats[(bias, window, carrier, formula)]
                if formula in RANKED_FORMULAS and err_vector is not None:
                    item.vector_errors.append(err_vector)
                if formula in RANKED_FORMULAS and err_projected is not None:
                    item.projected_errors.append(err_projected)
                if formula in RANKED_FORMULAS:
                    self_value = formulas["self_current_impl"]
                    self_err = log_error(value, self_value)
                    if self_err is not None:
                        item.self_impl_errors.append(self_err)
                if signed_rel_error is not None:
                    item.signed_projection_rel_errors.append(signed_rel_error)
                if flipped_signed_rel_error is not None:
                    item.flipped_signed_projection_rel_errors.append(flipped_signed_rel_error)


def best_formula(stats: dict[tuple[float, str, str, str], FormulaStats], bias: float, window: str, carrier: str, metric: str) -> tuple[str, float | None]:
    best_name = "insufficient_data"
    best_value: float | None = None
    for formula in RANKED_FORMULAS:
        item = stats.get((bias, window, carrier, formula))
        if item is None:
            continue
        values = getattr(item, metric)
        value = median_or_none(values)
        if value is None:
            continue
        if best_value is None or value < best_value:
            best_name = formula
            best_value = value
    return best_name, best_value


def sign_status(stats: dict[tuple[float, str, str, str], FormulaStats], carrier: str) -> str:
    item = stats.get((-20.0, "full", carrier, "self_current_impl"))
    if item is None:
        return "insufficient_data"
    normal = median_or_none(item.signed_projection_rel_errors)
    flipped = median_or_none(item.flipped_signed_projection_rel_errors)
    if normal is None or flipped is None:
        return "insufficient_data"
    if flipped < 0.5 * normal:
        return f"sign_flip_likely(normal_rel={normal:.6g}, flipped_rel={flipped:.6g})"
    if normal < 0.5 * flipped:
        return f"orientation_consistent(normal_rel={normal:.6g}, flipped_rel={flipped:.6g})"
    return f"inconclusive(normal_rel={normal:.6g}, flipped_rel={flipped:.6g})"
def write_summary(path: Path, stats: dict[tuple[float, str, str, str], FormulaStats], row_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    jn_vec, jn_vec_err = best_formula(stats, -20.0, "full", "Jn", "vector_errors")
    jp_vec, jp_vec_err = best_formula(stats, -20.0, "full", "Jp", "vector_errors")
    jn_proj, _ = best_formula(stats, -20.0, "full", "Jn", "projected_errors")
    jp_proj, _ = best_formula(stats, -20.0, "full", "Jp", "projected_errors")
    jn_self, _ = best_formula(stats, -20.0, "full", "Jn", "self_impl_errors")
    jp_self, _ = best_formula(stats, -20.0, "full", "Jp", "self_impl_errors")
    jn_sign_status = sign_status(stats, "Jn")
    jp_sign_status = sign_status(stats, "Jp")

    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# SG Edge Current Formula Sweep\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver default path changed: no\n")
        out.write("- self current source: internal edge source rows\n")
        out.write("- small-signal edges are marked and excluded from formula ranking\n")
        out.write("- F unit: V/cm\n")
        out.write("- Potential/qF unit: V\n")
        out.write("- length unit: cm\n")
        out.write("- density unit: cm^-3\n")
        out.write("- mobility unit: cm^2/V/s\n")
        out.write("- current-density unit: A/cm^2\n")
        out.write(f"- rows: {row_count}\n")
        out.write(f"- thermal_voltage_V: {VT_300K:.17g}\n\n")

        out.write("## -20V Window Winners\n\n")
        out.write("| window | carrier | best vector-mag formula | median log10 error | best projected formula | self-impl closest candidate |\n")
        out.write("| --- | --- | --- | --- | --- | --- |\n")
        for window in ("full", "junction", "center", "left_shoulder", "right_shoulder"):
            for carrier in ("Jn", "Jp"):
                best_vec, best_vec_error = best_formula(stats, -20.0, window, carrier, "vector_errors")
                best_projected, _ = best_formula(stats, -20.0, window, carrier, "projected_errors")
                best_self, _ = best_formula(stats, -20.0, window, carrier, "self_impl_errors")
                err_text = "" if best_vec_error is None else f"{best_vec_error:.6g}"
                out.write(f"| {window} | {carrier} | {best_vec} | {err_text} | {best_projected} | {best_self} |\n")

        out.write("\n## Required Answers\n\n")
        out.write(f"1. Best electron formula vs Sentaurus edge vector magnitude: {jn_vec}.\n")
        out.write(f"2. Best hole formula vs Sentaurus edge vector magnitude: {jp_vec}.\n")
        out.write(f"3. Best formulas vs Sentaurus edge projection: electron={jn_proj}, hole={jp_proj}.\n")
        out.write(f"4. Current self implementation closest candidate: electron={jn_self}, hole={jp_self}.\n")
        if "qf_sg_effective_density" in {jn_vec, jp_vec, jn_proj, jp_proj} or "sg_standard" in {jn_vec, jp_vec, jn_proj, jp_proj}:
            out.write("5. Bernoulli/SG exponential density averaging may matter in at least one main ranking; inspect sg_standard and qf_sg_effective_density columns by window.\n")
        else:
            out.write("5. No -20V full-window main ranking says the self implementation is simply missing Bernoulli/SG exponential density averaging.\n")
        density_mean_winners = {jn_vec, jp_vec, jn_proj, jp_proj, jn_self, jp_self} & {"qf_arithmetic_density", "qf_geometric_density", "qf_logarithmic_density", "qf_upwind_density", "qf_sg_effective_density"}
        if density_mean_winners:
            out.write(f"6. Density mean choice is active; winning/closest density means include: {', '.join(sorted(density_mean_winners))}.\n")
        else:
            out.write("6. Density mean choice is not isolated by the full-window winners.\n")
        out.write(f"7. Edge orientation/sign convention check from self implementation vs signed projection: electron={jn_sign_status}, hole={jp_sign_status}.\n")
        out.write("8. See -20V Window Winners for full/junction/center/left/right formula winners.\n")
        if jn_vec == "sg_standard" or jp_vec == "sg_standard":
            out.write("9. If standard SG reproduces Sentaurus while self does not, self SG formula is suspect.\n")
        else:
            out.write("9. Standard SG does not alone prove the self SG formula wrong in this ranking.\n")
        if "qf_logarithmic_density" in {jn_vec, jp_vec} or "qf_sg_effective_density" in {jn_vec, jp_vec}:
            out.write("10. A qF logarithmic/SG-effective winner would point to edge density averaging; inspect those columns per window.\n")
        else:
            out.write("10. qF logarithmic/SG-effective density is not the -20V full-window vector-mag winner.\n")
        proj_note = ""
        if jn_vec != jn_proj or jp_vec != jp_proj:
            proj_note = " Ranking differs between projected current and vector magnitude."
        out.write(
            "11. If Sentaurus projected current and vector magnitude prefer different formulas, distinguish edge scalar current from vector-current magnitude."
            f"{proj_note}\n"
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
    stats: dict[tuple[float, str, str, str], FormulaStats] = defaultdict(FormulaStats)
    output_rows: list[dict[str, Any]] = []

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
            if edge is None or edge.edge_length_cm <= 0.0:
                continue

            required_pairs = {
                "psi": endpoint_pair(nodes, bias, edge, "potential"),
                "phin": endpoint_pair(nodes, bias, edge, "electron_qf"),
                "phip": endpoint_pair(nodes, bias, edge, "hole_qf"),
                "n": endpoint_pair(nodes, bias, edge, "electron_density"),
                "p": endpoint_pair(nodes, bias, edge, "hole_density"),
                "mu_n": endpoint_pair(nodes, bias, edge, "electron_mobility"),
                "mu_p": endpoint_pair(nodes, bias, edge, "hole_mobility"),
            }
            if any(value is None for pair in required_pairs.values() for value in pair):
                continue
            psi0, psi1 = required_pairs["psi"]
            phin0, phin1 = required_pairs["phin"]
            phip0, phip1 = required_pairs["phip"]
            n0, n1 = required_pairs["n"]
            p0, p1 = required_pairs["p"]
            mu_n0, mu_n1 = required_pairs["mu_n"]
            mu_p0, mu_p1 = required_pairs["mu_p"]
            assert None not in (psi0, psi1, phin0, phin1, phip0, phip1, n0, n1, p0, p1, mu_n0, mu_n1, mu_p0, mu_p1)

            jn_sent_mag, jn_sent_proj, jn_sent_proj_signed = edge_vector_reference(
                nodes, bias, edge, "electron_current_density")
            jp_sent_mag, jp_sent_proj, jp_sent_proj_signed = edge_vector_reference(
                nodes, bias, edge, "hole_current_density")
            jn_self = parse_float(row.get("Jn_mag_used_A_per_cm2"))
            jp_self = parse_float(row.get("Jp_mag_used_A_per_cm2"))
            jn_formulas = formula_values(
                carrier="Jn",
                self_impl=jn_self,
                psi0=float(psi0),
                psi1=float(psi1),
                qf0=float(phin0),
                qf1=float(phin1),
                density0=float(n0),
                density1=float(n1),
                mu0=float(mu_n0),
                mu1=float(mu_n1),
                length_cm=edge.edge_length_cm,
                vt=args.thermal_voltage,
            )
            jp_formulas = formula_values(
                carrier="Jp",
                self_impl=jp_self,
                psi0=float(psi0),
                psi1=float(psi1),
                qf0=float(phip0),
                qf1=float(phip1),
                density0=float(p0),
                density1=float(p1),
                mu0=float(mu_p0),
                mu1=float(mu_p1),
                length_cm=edge.edge_length_cm,
                vt=args.thermal_voltage,
            )

            out_row: dict[str, Any] = {
                "bias_V": bias,
                "edge_id": edge.edge_id,
                "node0": edge.node0,
                "node1": edge.node1,
                "x0_um": edge.x0_um,
                "y0_um": edge.y0_um,
                "x1_um": edge.x1_um,
                "y1_um": edge.y1_um,
                "x_mid_um": edge.x_mid_um,
                "y_mid_um": edge.y_mid_um,
                "edge_length_cm": edge.edge_length_cm,
                "psi0": psi0,
                "psi1": psi1,
                "phin0": phin0,
                "phin1": phin1,
                "phip0": phip0,
                "phip1": phip1,
                "n0": n0,
                "n1": n1,
                "p0": p0,
                "p1": p1,
                "mu_n0": mu_n0,
                "mu_n1": mu_n1,
                "mu_p0": mu_p0,
                "mu_p1": mu_p1,
                "Jn_sent_edge_vector_mag_A_per_cm2": jn_sent_mag,
                "Jn_sent_edge_projected_A_per_cm2": jn_sent_proj,
                "Jn_sent_edge_projected_signed_A_per_cm2": jn_sent_proj_signed,
                "Jp_sent_edge_vector_mag_A_per_cm2": jp_sent_mag,
                "Jp_sent_edge_projected_A_per_cm2": jp_sent_proj,
                "Jp_sent_edge_projected_signed_A_per_cm2": jp_sent_proj_signed,
                "Jn_signal_label": "high_signal" if jn_sent_mag is not None and abs(jn_sent_mag) >= args.min_sent_current else "small_signal",
                "Jp_signal_label": "high_signal" if jp_sent_mag is not None and abs(jp_sent_mag) >= args.min_sent_current else "small_signal",
                "window_label": window_label(edge.x_mid_um),
            }
            windows = window_members(edge.x_mid_um)
            add_formula_outputs(
                out_row,
                stats,
                bias=bias,
                windows=windows,
                carrier="Jn",
                formulas=jn_formulas,
                sent_vector=jn_sent_mag,
                sent_projected=jn_sent_proj,
                sent_projected_signed=jn_sent_proj_signed,
                min_current=args.min_sent_current,
            )
            add_formula_outputs(
                out_row,
                stats,
                bias=bias,
                windows=windows,
                carrier="Jp",
                formulas=jp_formulas,
                sent_vector=jp_sent_mag,
                sent_projected=jp_sent_proj,
                sent_projected_signed=jp_sent_proj_signed,
                min_current=args.min_sent_current,
            )
            output_rows.append(out_row)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    headers = BASE_HEADERS + candidate_headers()
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in output_rows:
            writer.writerow({key: fmt(row.get(key)) for key in headers})
    write_summary(args.out_summary, stats, len(output_rows))


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
