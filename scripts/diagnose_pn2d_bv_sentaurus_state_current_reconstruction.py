#!/usr/bin/env python3
"""Reconstruct Sentaurus-state current support from exported qF, density, and mobility.

For each Sentaurus multibias export, this script builds edges on the Sentaurus
mesh and compares exported |Jn|/q, |Jp|/q against local edge reconstructions
mu * carrier * |grad quasi-Fermi|.  It reports the equivalent mobility needed to
explain Sentaurus current on each edge.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any, Iterable


Q_C = 1.602176634e-19

EDGE_FIELDS = [
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "edge_length_m",
    "sent_high_field_p95",
    "sent_high_e_current_p95",
    "sent_high_h_current_p95",
    "electron_density_avg_cm3",
    "hole_density_avg_cm3",
    "electron_mobility_avg_cm2_V_s",
    "hole_mobility_avg_cm2_V_s",
    "electron_qf_drop_V",
    "hole_qf_drop_V",
    "electron_qf_grad_abs_V_m",
    "hole_qf_grad_abs_V_m",
    "sent_e_current_A_cm2",
    "sent_h_current_A_cm2",
    "sent_e_flux_m2_s",
    "sent_h_flux_m2_s",
    "electron_linear_qf_flux_m2_s",
    "hole_linear_qf_flux_m2_s",
    "electron_current_over_linear_qf",
    "hole_current_over_linear_qf",
    "electron_mu_eff_cm2_V_s",
    "hole_mu_eff_cm2_V_s",
    "electron_mu_eff_over_mobility",
    "hole_mu_eff_over_mobility",
]

SUMMARY_FIELDS = [
    "bias_V",
    "edge_count",
    "high_field_edge_count",
    "electron_current_over_linear_qf_median",
    "hole_current_over_linear_qf_median",
    "electron_current_over_linear_qf_high_field_median",
    "hole_current_over_linear_qf_high_field_median",
    "electron_mu_eff_over_mobility_median",
    "hole_mu_eff_over_mobility_median",
    "electron_mu_eff_over_mobility_high_field_median",
    "hole_mu_eff_over_mobility_high_field_median",
    "electron_flux_sum",
    "hole_flux_sum",
    "electron_linear_qf_flux_sum",
    "hole_linear_qf_flux_sum",
    "electron_flux_over_linear_qf_sum",
    "hole_flux_over_linear_qf_sum",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sentaurus-root",
        type=Path,
        default=Path("build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias_0p25v"),
    )
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-16.0)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def finite_float(raw: Any, default: float = 0.0) -> float:
    value = optional_float(raw)
    return value if value is not None else default


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def in_bias_window(bias: float, lo: float, hi: float) -> bool:
    return min(lo, hi) - 1.0e-9 <= bias <= max(lo, hi) + 1.0e-9


def bias_from_dir(path: Path) -> float | None:
    name = path.name
    if not name.startswith("sentaurus_-") or not name.endswith("v"):
        return None
    text = name[len("sentaurus_-"):-1].replace("p", ".")
    try:
        return -float(text)
    except ValueError:
        return None


def load_nodes(case_dir: Path) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    for row in read_csv(case_dir / "nodes.csv"):
        nodes[int(row["id"])] = (float(row["x_um"]) * 1.0e-6, float(row["y_um"]) * 1.0e-6)
    return nodes


def load_edges(case_dir: Path) -> list[tuple[int, int]]:
    seen: set[tuple[int, int]] = set()
    edges: list[tuple[int, int]] = []
    for row in read_csv(case_dir / "elements.csv"):
        ids = [int(row["node0"]), int(row["node1"]), int(row["node2"])]
        for i in range(3):
            a = ids[i]
            b = ids[(i + 1) % 3]
            key = (a, b) if a <= b else (b, a)
            if key not in seen:
                seen.add(key)
                edges.append(key)
    return edges


def numeric_components(row: dict[str, str]) -> list[float]:
    values: list[float] = []
    for key, raw in row.items():
        if key == "node_id":
            continue
        value = optional_float(raw)
        if value is not None:
            values.append(value)
    return values


def load_scalar_field(fields_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted(fields_dir.glob(f"{name}_region*.csv")):
        for row in read_csv(path):
            values = numeric_components(row)
            result[int(row["node_id"])] = values[0] if values else 0.0
    return result


def load_magnitude_field(fields_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted(fields_dir.glob(f"{name}_region*.csv")):
        for row in read_csv(path):
            values = numeric_components(row)
            result[int(row["node_id"])] = math.sqrt(sum(v * v for v in values)) if values else 0.0
    return result


def avg(field: dict[int, float], a: int, b: int) -> float | None:
    if a not in field or b not in field:
        return None
    return 0.5 * (field[a] + field[b])


def drop(field: dict[int, float], a: int, b: int) -> float | None:
    if a not in field or b not in field:
        return None
    return field[b] - field[a]


def particle_flux_from_current_a_cm2(current_a_cm2: float | None) -> float | None:
    if current_a_cm2 is None:
        return None
    return abs(current_a_cm2) * 1.0e4 / Q_C


def linear_qf_flux_m2_s(mobility_cm2: float | None,
                         density_cm3: float | None,
                         qf_drop_v: float | None,
                         length_m: float) -> float | None:
    if mobility_cm2 is None or density_cm3 is None or qf_drop_v is None or length_m <= 0.0:
        return None
    mobility_m2 = max(mobility_cm2, 0.0) * 1.0e-4
    density_m3 = max(density_cm3, 0.0) * 1.0e6
    return mobility_m2 * density_m3 * abs(qf_drop_v) / length_m


def mu_eff_cm2_v_s(flux_m2_s: float | None,
                    density_cm3: float | None,
                    qf_drop_v: float | None,
                    length_m: float) -> float | None:
    if flux_m2_s is None or density_cm3 is None or qf_drop_v is None:
        return None
    if density_cm3 <= 0.0 or abs(qf_drop_v) <= 0.0 or length_m <= 0.0:
        return None
    density_m3 = density_cm3 * 1.0e6
    mu_m2 = flux_m2_s * length_m / (density_m3 * abs(qf_drop_v))
    return mu_m2 * 1.0e4


def percentile_cut(values: list[float], q: float) -> float | None:
    finite = sorted(v for v in values if math.isfinite(v))
    if not finite:
        return None
    index = int(round((len(finite) - 1) * q))
    return finite[index]


def load_case_fields(case_dir: Path) -> dict[str, dict[int, float]]:
    fields = case_dir / "fields"
    return {
        "e_density": load_scalar_field(fields, "eDensity"),
        "h_density": load_scalar_field(fields, "hDensity"),
        "e_mobility": load_scalar_field(fields, "eMobility"),
        "h_mobility": load_scalar_field(fields, "hMobility"),
        "e_qf": load_scalar_field(fields, "eQuasiFermiPotential"),
        "h_qf": load_scalar_field(fields, "hQuasiFermiPotential"),
        "e_current": load_magnitude_field(fields, "eCurrentDensity"),
        "h_current": load_magnitude_field(fields, "hCurrentDensity"),
        "electric_field": load_magnitude_field(fields, "ElectricField"),
    }


def make_case_rows(case_dir: Path, bias: float) -> list[dict[str, Any]]:
    nodes = load_nodes(case_dir)
    edges = load_edges(case_dir)
    fields = load_case_fields(case_dir)
    edge_e_fields: list[float] = []
    edge_e_currents: list[float] = []
    edge_h_currents: list[float] = []
    for a, b in edges:
        efield = avg(fields["electric_field"], a, b)
        ecurr = avg(fields["e_current"], a, b)
        hcurr = avg(fields["h_current"], a, b)
        if efield is not None:
            edge_e_fields.append(efield)
        if ecurr is not None:
            edge_e_currents.append(ecurr)
        if hcurr is not None:
            edge_h_currents.append(hcurr)
    efield_p95 = percentile_cut(edge_e_fields, 0.95)
    ecurr_p95 = percentile_cut(edge_e_currents, 0.95)
    hcurr_p95 = percentile_cut(edge_h_currents, 0.95)

    rows: list[dict[str, Any]] = []
    for edge_id, (a, b) in enumerate(edges):
        if a not in nodes or b not in nodes:
            continue
        x0, y0 = nodes[a]
        x1, y1 = nodes[b]
        length = math.hypot(x1 - x0, y1 - y0)
        if length <= 1.0e-30:
            continue
        n_avg = avg(fields["e_density"], a, b)
        p_avg = avg(fields["h_density"], a, b)
        mun = avg(fields["e_mobility"], a, b)
        mup = avg(fields["h_mobility"], a, b)
        e_qf_drop = drop(fields["e_qf"], a, b)
        h_qf_drop = drop(fields["h_qf"], a, b)
        e_current = avg(fields["e_current"], a, b)
        h_current = avg(fields["h_current"], a, b)
        e_flux = particle_flux_from_current_a_cm2(e_current)
        h_flux = particle_flux_from_current_a_cm2(h_current)
        e_linear = linear_qf_flux_m2_s(mun, n_avg, e_qf_drop, length)
        h_linear = linear_qf_flux_m2_s(mup, p_avg, h_qf_drop, length)
        e_mu_eff = mu_eff_cm2_v_s(e_flux, n_avg, e_qf_drop, length)
        h_mu_eff = mu_eff_cm2_v_s(h_flux, p_avg, h_qf_drop, length)
        efield = avg(fields["electric_field"], a, b)
        rows.append({
            "bias_V": bias,
            "edge_id": edge_id,
            "node0": a,
            "node1": b,
            "edge_length_m": length,
            "sent_high_field_p95": bool(efield_p95 is not None and efield is not None and efield >= efield_p95),
            "sent_high_e_current_p95": bool(ecurr_p95 is not None and e_current is not None and e_current >= ecurr_p95),
            "sent_high_h_current_p95": bool(hcurr_p95 is not None and h_current is not None and h_current >= hcurr_p95),
            "electron_density_avg_cm3": n_avg,
            "hole_density_avg_cm3": p_avg,
            "electron_mobility_avg_cm2_V_s": mun,
            "hole_mobility_avg_cm2_V_s": mup,
            "electron_qf_drop_V": e_qf_drop,
            "hole_qf_drop_V": h_qf_drop,
            "electron_qf_grad_abs_V_m": abs(e_qf_drop) / length if e_qf_drop is not None else None,
            "hole_qf_grad_abs_V_m": abs(h_qf_drop) / length if h_qf_drop is not None else None,
            "sent_e_current_A_cm2": e_current,
            "sent_h_current_A_cm2": h_current,
            "sent_e_flux_m2_s": e_flux,
            "sent_h_flux_m2_s": h_flux,
            "electron_linear_qf_flux_m2_s": e_linear,
            "hole_linear_qf_flux_m2_s": h_linear,
            "electron_current_over_linear_qf": ratio(e_flux, e_linear),
            "hole_current_over_linear_qf": ratio(h_flux, h_linear),
            "electron_mu_eff_cm2_V_s": e_mu_eff,
            "hole_mu_eff_cm2_V_s": h_mu_eff,
            "electron_mu_eff_over_mobility": ratio(e_mu_eff, mun),
            "hole_mu_eff_over_mobility": ratio(h_mu_eff, mup),
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def median(rows: list[dict[str, Any]], key: str) -> float | None:
    values = finite_values(rows, key)
    return statistics.median(values) if values else None


def sum_values(rows: list[dict[str, Any]], key: str) -> float:
    return sum(finite_values(rows, key))


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for bias in sorted({float(row["bias_V"]) for row in rows}):
        subset = [row for row in rows if abs(float(row["bias_V"]) - bias) <= 1.0e-9]
        high_field = [row for row in subset if row["sent_high_field_p95"]]
        e_flux_sum = sum_values(subset, "sent_e_flux_m2_s")
        h_flux_sum = sum_values(subset, "sent_h_flux_m2_s")
        e_linear_sum = sum_values(subset, "electron_linear_qf_flux_m2_s")
        h_linear_sum = sum_values(subset, "hole_linear_qf_flux_m2_s")
        result.append({
            "bias_V": bias,
            "edge_count": len(subset),
            "high_field_edge_count": len(high_field),
            "electron_current_over_linear_qf_median": median(subset, "electron_current_over_linear_qf"),
            "hole_current_over_linear_qf_median": median(subset, "hole_current_over_linear_qf"),
            "electron_current_over_linear_qf_high_field_median": median(high_field, "electron_current_over_linear_qf"),
            "hole_current_over_linear_qf_high_field_median": median(high_field, "hole_current_over_linear_qf"),
            "electron_mu_eff_over_mobility_median": median(subset, "electron_mu_eff_over_mobility"),
            "hole_mu_eff_over_mobility_median": median(subset, "hole_mu_eff_over_mobility"),
            "electron_mu_eff_over_mobility_high_field_median": median(high_field, "electron_mu_eff_over_mobility"),
            "hole_mu_eff_over_mobility_high_field_median": median(high_field, "hole_mu_eff_over_mobility"),
            "electron_flux_sum": e_flux_sum,
            "hole_flux_sum": h_flux_sum,
            "electron_linear_qf_flux_sum": e_linear_sum,
            "hole_linear_qf_flux_sum": h_linear_sum,
            "electron_flux_over_linear_qf_sum": ratio(e_flux_sum, e_linear_sum),
            "hole_flux_over_linear_qf_sum": ratio(h_flux_sum, h_linear_sum),
        })
    return result


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def selected_case_dirs(root: Path, bias_min: float, bias_max: float) -> list[tuple[float, Path]]:
    cases: list[tuple[float, Path]] = []
    for path in sorted(root.glob("sentaurus_-*v")):
        if not path.is_dir():
            continue
        bias = bias_from_dir(path)
        if bias is None or not in_bias_window(bias, bias_min, bias_max):
            continue
        cases.append((bias, path))
    return sorted(cases)


def main() -> int:
    args = parse_args()
    rows: list[dict[str, Any]] = []
    for bias, case_dir in selected_case_dirs(args.sentaurus_root, args.bias_min, args.bias_max):
        rows.extend(make_case_rows(case_dir, bias))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = summarize(rows)
    high_field_rows = [row for row in rows if row["sent_high_field_p95"]]
    write_csv(args.out_dir / "sentaurus_state_current_edges.csv", EDGE_FIELDS, rows)
    write_csv(args.out_dir / "sentaurus_state_current_high_field_edges.csv", EDGE_FIELDS, high_field_rows)
    write_csv(args.out_dir / "sentaurus_state_current_summary.csv", SUMMARY_FIELDS, summary_rows)
    payload = {
        "schema": "pn2d.sentaurus_state_current_reconstruction.v1",
        "bias_count": len(summary_rows),
        "edge_row_count": len(rows),
        "high_field_edge_row_count": len(high_field_rows),
        "summary": summary_rows,
    }
    (args.out_dir / "sentaurus_state_current_summary.json").write_text(
        json.dumps(clean_json(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "bias_count": len(summary_rows),
        "edge_row_count": len(rows),
        "high_field_edge_row_count": len(high_field_rows),
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
