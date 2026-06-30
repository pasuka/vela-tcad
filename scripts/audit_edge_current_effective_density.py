#!/usr/bin/env python3
"""Audit edge-current effective carrier density inferred from J = q mu n_eff grad(qF)."""

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
DELTA_PHI_FLOOR_V = 1.0e-12
RATIO_FLOOR = 1.0e-300
WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}


@dataclass(frozen=True)
class Edge:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    length_cm: float

    @property
    def x_mid_um(self) -> float:
        return 0.5 * (self.x0_um + self.x1_um)

    @property
    def y_mid_um(self) -> float:
        return 0.5 * (self.y0_um + self.y1_um)

    @property
    def tangent(self) -> tuple[float, float]:
        dx = self.x1_um - self.x0_um
        dy = self.y1_um - self.y0_um
        norm = math.hypot(dx, dy)
        if norm <= 0.0:
            return 0.0, 0.0
        return dx / norm, dy / norm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/a1_b1_internal_source_current_audit/avalanche_internal_source_current_audit.csv")
    parser.add_argument(
        "--self-edge-topology",
        type=Path,
        default=(
            REPO / "build/diagnostics/a1_b1_internal_source_current_audit"
            / "avalanche_internal_source_current_audit_case/BV-A1-B1p00-internal-source-current-audit"
            / "BV-A1-B1p00-internal-source-current-audit_sg_avalanche_edges.csv"
        ),
    )
    parser.add_argument(
        "--node-compare",
        type=Path,
        default=(
            REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
            / "coarse_previous_full20_fields_20260629_155228/coarse_node_field_compare_aligned.csv"
        ),
    )
    parser.add_argument("--sg-formula-sweep", type=Path, default=REPO / "build/diagnostics/sg_edge_current_formula_sweep.csv")
    parser.add_argument("--cell-vector-reconstruction", type=Path, default=REPO / "build/diagnostics/edge_to_cell_vector_current_reconstruction.csv")
    parser.add_argument("--elements-csv", type=Path, default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_fields_20260629_155228/elements.csv")
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--delta-phi-floor", type=float, default=DELTA_PHI_FLOOR_V)
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/edge_current_effective_density_audit.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/edge_current_effective_density_audit_summary.md")
    return parser.parse_args()


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items() if key is not None}


def f(value: Any, default: float = math.nan) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def parse_biases(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def nearest_bias(value: float, targets: list[float], tol: float) -> float | None:
    best = min(targets, key=lambda bias: abs(bias - value))
    return best if abs(best - value) <= tol else None


def safe_ratio(value: float | None, reference: float | None) -> float | None:
    if value is None or reference is None or not math.isfinite(value) or not math.isfinite(reference):
        return None
    if abs(reference) <= RATIO_FLOOR:
        return None
    ratio = value / reference
    return ratio if math.isfinite(ratio) else None


def log_error(value: float | None, reference: float | None) -> float:
    ratio = safe_ratio(value, reference)
    if ratio is None or ratio <= 0.0:
        return math.inf
    return abs(math.log10(abs(ratio)))


def load_edges(path: Path) -> dict[int, Edge]:
    edges: dict[int, Edge] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            edge_id = int(f(row["edge_id"]))
            if edge_id in edges:
                continue
            edges[edge_id] = Edge(
                edge_id=edge_id,
                node0=int(f(row["node0"])),
                node1=int(f(row["node1"])),
                x0_um=f(row["x0_um"]),
                y0_um=f(row["y0_um"]),
                x1_um=f(row["x1_um"]),
                y1_um=f(row["y1_um"]),
                length_cm=f(row["edge_length_m"]) * 100.0,
            )
    return edges


def load_internal(path: Path, biases: list[float], tol: float) -> dict[tuple[float, int], dict[str, float]]:
    rows: dict[tuple[float, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            if row.get("source_location_type", "edge") != "edge":
                continue
            bias = nearest_bias(f(row.get("bias_V")), biases, tol)
            if bias is None:
                continue
            edge_id = int(f(row.get("source_entity_id")))
            rows[(bias, edge_id)] = {
                "Jn": f(row.get("Jn_mag_used_A_per_cm2"), 0.0),
                "Jp": f(row.get("Jp_mag_used_A_per_cm2"), 0.0),
            }
    return rows


def load_node_compare(path: Path) -> dict[tuple[float, int], dict[str, dict[str, float]]]:
    nodes: dict[tuple[float, int], dict[str, dict[str, float]]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            bias = f(row.get("bias_V"))
            node_id = int(f(row.get("node_id")))
            quantity = row.get("quantity", "")
            nodes.setdefault((bias, node_id), {})[quantity] = {
                "sent": f(row.get("sentaurus_value")),
                "self": f(row.get("vela_value_scaled_to_sentaurus_units")),
            }
    return nodes


def node_value(nodes: dict[tuple[float, int], dict[str, dict[str, float]]], bias: float, node_id: int, quantity: str, side: str) -> float:
    return nodes.get((bias, node_id), {}).get(quantity, {}).get(side, math.nan)


def vector_ref(nodes: dict[tuple[float, int], dict[str, dict[str, float]]], bias: float, edge: Edge, prefix: str) -> tuple[float, float]:
    jx0 = node_value(nodes, bias, edge.node0, f"{prefix}_x", "sent")
    jy0 = node_value(nodes, bias, edge.node0, f"{prefix}_y", "sent")
    jx1 = node_value(nodes, bias, edge.node1, f"{prefix}_x", "sent")
    jy1 = node_value(nodes, bias, edge.node1, f"{prefix}_y", "sent")
    jx = 0.5 * (jx0 + jx1)
    jy = 0.5 * (jy0 + jy1)
    tx, ty = edge.tangent
    return math.hypot(jx, jy), abs(jx * tx + jy * ty)


def harmonic(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0 or a + b == 0.0:
        return math.nan
    return 2.0 * a * b / (a + b)


def geometric(a: float, b: float) -> float:
    if a < 0.0 or b < 0.0:
        return math.nan
    return math.sqrt(a * b)


def logarithmic(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return math.nan
    if abs(a - b) <= 1.0e-14 * max(abs(a), abs(b), 1.0):
        return 0.5 * (a + b)
    return (b - a) / (math.log(b) - math.log(a))


def infer_density(current: float, mu: float, delta_phi: float, length_cm: float) -> float | None:
    denom = Q_C * abs(mu) * abs(delta_phi) / length_cm if length_cm > 0.0 else 0.0
    if denom <= RATIO_FLOOR:
        return None
    value = abs(current) / denom
    return value if math.isfinite(value) else None


def best_mean(target: float | None, means: dict[str, float]) -> str:
    if target is None or not math.isfinite(target) or target <= 0.0:
        return "insufficient_data"
    candidates = [(name, log_error(target, value)) for name, value in means.items() if value > 0.0 and math.isfinite(value)]
    if not candidates:
        return "insufficient_data"
    return min(candidates, key=lambda item: item[1])[0]


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


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    biases = parse_biases(args.biases)
    edges = load_edges(args.self_edge_topology)
    internal = load_internal(args.internal_audit, biases, args.bias_tol)
    nodes = load_node_compare(args.node_compare)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        for edge in edges.values():
            int_row = internal.get((bias, edge.edge_id))
            if int_row is None:
                continue
            for carrier, prefix, density_q, mobility_q, phi_q, self_j_key in (
                ("electron", "electron_current_density", "electron_density", "electron_mobility", "electron_qf", "Jn"),
                ("hole", "hole_current_density", "hole_density", "hole_mobility", "hole_qf", "Jp"),
            ):
                d0 = node_value(nodes, bias, edge.node0, density_q, "self")
                d1 = node_value(nodes, bias, edge.node1, density_q, "self")
                mu0 = node_value(nodes, bias, edge.node0, mobility_q, "self")
                mu1 = node_value(nodes, bias, edge.node1, mobility_q, "self")
                p0 = node_value(nodes, bias, edge.node0, phi_q, "self")
                p1 = node_value(nodes, bias, edge.node1, phi_q, "self")
                if any(not math.isfinite(value) for value in (d0, d1, mu0, mu1, p0, p1)):
                    continue
                delta = p1 - p0
                if abs(delta) < args.delta_phi_floor:
                    continue
                mu_edge = 0.5 * (mu0 + mu1)
                sent_mag, sent_proj = vector_ref(nodes, bias, edge, prefix)
                self_j = int_row[self_j_key]
                eff_self = infer_density(self_j, mu_edge, delta, edge.length_cm)
                eff_sent_mag = infer_density(sent_mag, mu_edge, delta, edge.length_cm)
                eff_sent_proj = infer_density(sent_proj, mu_edge, delta, edge.length_cm)
                means = {
                    "arithmetic": 0.5 * (d0 + d1),
                    "geometric": geometric(d0, d1),
                    "harmonic": harmonic(d0, d1),
                    "logarithmic": logarithmic(d0, d1),
                    "upwind_from_phi_direction": d0 if delta >= 0.0 else d1,
                    "sg_effective_if_available": math.nan,
                }
                rows.append({
                    "bias_V": bias,
                    "edge_id": edge.edge_id,
                    "node0": edge.node0,
                    "node1": edge.node1,
                    "x_mid_um": edge.x_mid_um,
                    "y_mid_um": edge.y_mid_um,
                    "edge_length_cm": edge.length_cm,
                    "carrier": carrier,
                    "density0_cm_minus3": d0,
                    "density1_cm_minus3": d1,
                    "mu0_cm2_per_Vs": mu0,
                    "mu1_cm2_per_Vs": mu1,
                    "phi0_V": p0,
                    "phi1_V": p1,
                    "delta_phi_V": delta,
                    "grad_phi_abs_V_per_cm": abs(delta) / edge.length_cm,
                    "density_arithmetic": means["arithmetic"],
                    "density_geometric": means["geometric"],
                    "density_harmonic": means["harmonic"],
                    "density_logarithmic": means["logarithmic"],
                    "density_upwind_from_phi_direction": means["upwind_from_phi_direction"],
                    "density_sg_effective_if_available": means["sg_effective_if_available"],
                    "J_self_internal_A_per_cm2": self_j,
                    "J_sent_vector_mag_A_per_cm2": sent_mag,
                    "J_sent_edge_projection_A_per_cm2": sent_proj,
                    "density_eff_self_from_J": eff_self,
                    "density_eff_sent_vector_from_J": eff_sent_mag,
                    "density_eff_sent_projection_from_J": eff_sent_proj,
                    "self_eff_over_arithmetic": safe_ratio(eff_self, means["arithmetic"]),
                    "self_eff_over_geometric": safe_ratio(eff_self, means["geometric"]),
                    "self_eff_over_logarithmic": safe_ratio(eff_self, means["logarithmic"]),
                    "sent_vector_eff_over_arithmetic": safe_ratio(eff_sent_mag, means["arithmetic"]),
                    "sent_vector_eff_over_geometric": safe_ratio(eff_sent_mag, means["geometric"]),
                    "sent_vector_eff_over_logarithmic": safe_ratio(eff_sent_mag, means["logarithmic"]),
                    "self_eff_over_sent_vector_eff": safe_ratio(eff_self, eff_sent_mag),
                    "self_J_over_sent_vector_J": safe_ratio(self_j, sent_mag),
                    "self_J_over_sent_projection_J": safe_ratio(self_j, sent_proj),
                    "best_density_mean_for_self": best_mean(eff_self, means),
                    "best_density_mean_for_sent_vector": best_mean(eff_sent_mag, means),
                    "best_density_mean_for_sent_projection": best_mean(eff_sent_proj, means),
                    "window_label": window_label(edge.x_mid_um),
                })
    return rows


HEADER = [
    "bias_V", "edge_id", "node0", "node1", "x_mid_um", "y_mid_um", "edge_length_cm", "carrier",
    "density0_cm_minus3", "density1_cm_minus3", "mu0_cm2_per_Vs", "mu1_cm2_per_Vs",
    "phi0_V", "phi1_V", "delta_phi_V", "grad_phi_abs_V_per_cm",
    "density_arithmetic", "density_geometric", "density_harmonic", "density_logarithmic",
    "density_upwind_from_phi_direction", "density_sg_effective_if_available",
    "J_self_internal_A_per_cm2", "J_sent_vector_mag_A_per_cm2", "J_sent_edge_projection_A_per_cm2",
    "density_eff_self_from_J", "density_eff_sent_vector_from_J", "density_eff_sent_projection_from_J",
    "self_eff_over_arithmetic", "self_eff_over_geometric", "self_eff_over_logarithmic",
    "sent_vector_eff_over_arithmetic", "sent_vector_eff_over_geometric", "sent_vector_eff_over_logarithmic",
    "self_eff_over_sent_vector_eff", "self_J_over_sent_vector_J", "self_J_over_sent_projection_J",
    "best_density_mean_for_self", "best_density_mean_for_sent_vector", "best_density_mean_for_sent_projection",
    "window_label",
]


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.17g}"
    return str(value)


def median_or_none(values: list[float]) -> float | None:
    clean = [value for value in values if math.isfinite(value)]
    return median(clean) if clean else None


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    by_carrier = defaultdict(list)
    by_window_carrier = defaultdict(list)
    mean_counts_self = Counter()
    mean_counts_sent_vec = Counter()
    mean_counts_sent_proj = Counter()
    for row in rows:
        if float(row["bias_V"]) != -20.0:
            continue
        carrier = str(row["carrier"])
        by_carrier[carrier].append(row)
        for window in window_members(float(row["x_mid_um"])):
            by_window_carrier[(window, carrier)].append(row)
        mean_counts_self[str(row["best_density_mean_for_self"])] += 1
        mean_counts_sent_vec[str(row["best_density_mean_for_sent_vector"])] += 1
        mean_counts_sent_proj[str(row["best_density_mean_for_sent_projection"])] += 1

    def med(rows_in: list[dict[str, Any]], key: str) -> float | None:
        return median_or_none([float(row[key]) for row in rows_in if row.get(key) not in (None, "")])

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Edge Current Effective Density Audit\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver default path changed: no\n")
        out.write("- current-density unit: A/cm^2\n")
        out.write("- density unit: cm^-3\n")
        out.write("- mobility unit: cm^2/V/s\n")
        out.write("- potential/qF unit: V\n")
        out.write("- length unit: cm\n")
        out.write("- filtered edges: |delta_phi| below threshold are skipped\n")
        out.write(f"- rows: {len(rows)}\n\n")
        out.write("## -20V Median Ratios\n\n")
        out.write("| carrier | self_eff/sent_vector_eff | self_J/sent_vector_J | self_J/sent_projection_J |\n")
        out.write("| --- | --- | --- | --- |\n")
        for carrier in ("electron", "hole"):
            subset = by_carrier.get(carrier, [])
            out.write(
                f"| {carrier} | {format_ratio(med(subset, 'self_eff_over_sent_vector_eff'))} | "
                f"{format_ratio(med(subset, 'self_J_over_sent_vector_J'))} | "
                f"{format_ratio(med(subset, 'self_J_over_sent_projection_J'))} |\n"
            )
        out.write("\n## Required Answers\n\n")
        electron_gap = med(by_carrier.get("electron", []), "self_eff_over_sent_vector_eff")
        hole_gap = med(by_carrier.get("hole", []), "self_eff_over_sent_vector_eff")
        out.write(f"1. At -20V, self inferred effective density vs Sentaurus: electron={format_ratio(electron_gap)}, hole={format_ratio(hole_gap)}.\n")
        out.write(f"2. Dominant best density mean for self current: {mean_counts_self.most_common(1)[0][0] if mean_counts_self else 'insufficient_data'}.\n")
        out.write(f"3. Dominant best density mean for Sentaurus vector current: {mean_counts_sent_vec.most_common(1)[0][0] if mean_counts_sent_vec else 'insufficient_data'}.\n")
        out.write(f"4. Dominant best density mean for Sentaurus projected current: {mean_counts_sent_proj.most_common(1)[0][0] if mean_counts_sent_proj else 'insufficient_data'}.\n")
        out.write("5. Compare best-mean counters above to see whether self clusters near geometric/upwind while Sentaurus clusters near arithmetic/log/SG.\n")
        j_gap_e = med(by_carrier.get("electron", []), "self_J_over_sent_vector_J")
        j_gap_h = med(by_carrier.get("hole", []), "self_J_over_sent_vector_J")
        out.write(f"6. Inferred density gap tracks J ratio when qF/mu/length are shared; median J ratios are electron={format_ratio(j_gap_e)}, hole={format_ratio(j_gap_h)}.\n")
        worse = "electron" if (j_gap_e or 1.0) < (j_gap_h or 1.0) else "hole"
        out.write(f"7. Stronger carrier current problem by median self/Sentaurus J ratio: {worse}.\n")
        window_lines = []
        for window in ("center", "left_shoulder", "right_shoulder"):
            vals = []
            for carrier in ("electron", "hole"):
                vals.append(med(by_window_carrier.get((window, carrier), []), "self_J_over_sent_vector_J"))
            vals_clean = [value for value in vals if value is not None]
            if vals_clean:
                window_lines.append((window, min(vals_clean)))
        strongest_window = min(window_lines, key=lambda item: item[1])[0] if window_lines else "insufficient_data"
        out.write(f"8. Strongest spatial problem among center/shoulders: {strongest_window}.\n")
        if mean_counts_self and mean_counts_sent_vec and mean_counts_self.most_common(1)[0][0] != mean_counts_sent_vec.most_common(1)[0][0]:
            out.write("9. Density averaging plausibly explains a meaningful part of the gap; recommend a guarded avalanche-source-only edge density-average option for testing.\n")
        else:
            out.write("9. Density averaging alone is not isolated; do not change source density averaging from this evidence alone.\n")
        out.write("10. If density averaging is insufficient, continue edge orientation / box-face normal weighting audit.\n")


def format_ratio(value: float | None) -> str:
    return "insufficient_data" if value is None else f"{value:.6g}"


def run(args: argparse.Namespace) -> None:
    rows = make_rows(args)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in HEADER})
    write_summary(args.out_summary, rows)


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
