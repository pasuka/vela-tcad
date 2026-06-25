#!/usr/bin/env python3
"""Diagnose PN2D BV ionization magnitude from field, alpha, seed, and source evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
DEFAULT_DECK = REPO / "build" / "pn2d_bv_gradqf_full20" / "simulation_bv_genius_cell_gradient_full20.json"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
DEFAULT_REFERENCE_CURVE = (
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reference_curves" /
    "pn2d_sentaurus2018_bv_reference.csv"
)
DEFAULT_SENT_FIELD_ROOTS = [
    REPO / "reference_tcad" / "pn2d_sentaurus2018" / "sim_fields",
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "sim_fields",
    REPO / "build" / "reference_tcad" / "pn2d_sentaurus2018" / "sim_fields",
    REPO / "build" / "reference_tcad" / "pn2d_sentaurus2018_full" / "sim_fields",
]
DEFAULT_OUT_DIR = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reports" / "ionization_magnitude"
TARGET_BIASES = [round(-16.0 - 0.25 * i, 12) for i in range(17)]
FIELD_RATIO_THRESHOLD = 0.9
SEED_RATIO_THRESHOLD = 0.75
Q = 1.602176634e-19
KB = 1.380649e-23


DEFAULT_IMPACT = {
    "model": "none",
    "electron_A_m_inv": 7.03e7,
    "electron_B_V_m": 1.231e8,
    "hole_A_m_inv": 1.582e8,
    "hole_B_V_m": 2.036e8,
    "electron_a_low_m_inv": 7.03e7,
    "electron_a_high_m_inv": 7.03e7,
    "electron_b_low_V_m": 1.231e8,
    "electron_b_high_V_m": 1.231e8,
    "hole_a_low_m_inv": 1.582e8,
    "hole_a_high_m_inv": 6.71e7,
    "hole_b_low_V_m": 2.036e8,
    "hole_b_high_V_m": 1.693e8,
    "switch_field_V_m": 4.0e7,
    "phonon_energy_eV": 0.063,
    "reference_temperature_K": 300.0,
    "temperature_K": 300.0,
    "minimum_field_V_m": 0.0,
}


def float_or_none(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def absolutize(value: str, base: Path) -> str:
    p = Path(value)
    if p.is_absolute():
        return str(p)
    return str((base / p).resolve())


def bias_token(value: float) -> str:
    token = f"{abs(value):.6f}".replace(".", "p")
    return ("m" if value < 0.0 else "") + token


def sorted_unique(values: Iterable[float]) -> list[float]:
    seen: dict[float, float] = {}
    for value in values:
        seen[round(value, 12)] = value
    return sorted(seen.values(), reverse=True)


def prepare_deck(base_deck: Path, out_dir: Path, contacts: list[str]) -> Path:
    base_deck = base_deck.resolve()
    deck_dir = base_deck.parent
    data = json.loads(base_deck.read_text(encoding="utf-8"))
    for key in ["mesh_file", "node_doping_file", "materials_file"]:
        if key in data and isinstance(data[key], str):
            data[key] = absolutize(data[key], deck_dir)

    sweep = dict(data.get("sweep", {}))
    base_points = [float(v) for v in sweep.get("bias_points", [])]
    path_points = [v for v in base_points if v >= min(TARGET_BIASES) or abs(v - min(TARGET_BIASES)) < 1.0e-9]
    sweep["bias_points"] = sorted_unique(path_points + TARGET_BIASES)
    sweep["start"] = 0.0
    sweep["stop"] = min(TARGET_BIASES)
    sweep["step"] = -0.25
    sweep["write_vtk"] = False
    sweep["write_state_every_point_prefix"] = str((out_dir / "states" / "accepted_state").resolve())
    sweep["output_csv"] = str((out_dir / "ionization_magnitude_curve.csv").resolve())
    data["output_csv"] = str((out_dir / "ionization_magnitude_curve.csv").resolve())

    diagnostics = dict(sweep.get("diagnostics", {}))
    diagnostics["sg_avalanche_edges"] = {
        "enabled": True,
        "csv_file": str((out_dir / "sg_avalanche_edges.csv").resolve()),
    }
    diagnostics["terminal_current_method_compare"] = {
        "enabled": True,
        "contacts": contacts,
        "csv_file": str((out_dir / "terminal_current_method_compare.csv").resolve()),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": contacts,
        "csv_file": str((out_dir / "continuity_balance.csv").resolve()),
    }
    sweep["diagnostics"] = diagnostics
    data["sweep"] = sweep

    out_dir.mkdir(parents=True, exist_ok=True)
    deck_path = out_dir / "simulation_ionization_magnitude.json"
    deck_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return deck_path


def load_nodes(path: Path) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    for row in load_csv(path):
        node_id = int(row["id"])
        nodes[node_id] = (float(row["x_um"]) * 1.0e-6, float(row["y_um"]) * 1.0e-6)
    return nodes


def load_elements(path: Path) -> list[tuple[int, int, int]]:
    elements: list[tuple[int, int, int]] = []
    for row in load_csv(path):
        elements.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
    return elements


def load_doping(path: Path) -> dict[int, float]:
    doping: dict[int, float] = {}
    for row in load_csv(path):
        donors = float(row.get("donors_cm3", "0") or 0.0)
        acceptors = float(row.get("acceptors_cm3", "0") or 0.0)
        doping[int(row["node_id"])] = donors - acceptors
    return doping


def load_scalar_field(path: Path, column: str = "component0") -> dict[int, float]:
    values: dict[int, float] = {}
    for row in load_csv(path):
        node_key = "node_id" if "node_id" in row else "id"
        values[int(row[node_key])] = float(row[column])
    return values


def load_state_potential(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["psi"]) for row in load_csv(path)}


def triangle_gradient(coords: list[tuple[float, float]], values: list[float]) -> tuple[float, float] | None:
    (x0, y0), (x1, y1), (x2, y2) = coords
    v0, v1, v2 = values
    det = x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1)
    if abs(det) <= 1.0e-40:
        return None
    gx = (v0 * (y1 - y2) + v1 * (y2 - y0) + v2 * (y0 - y1)) / det
    gy = (v0 * (x2 - x1) + v1 * (x0 - x2) + v2 * (x1 - x0)) / det
    return gx, gy


def compute_junction_geometry(nodes: dict[int, tuple[float, float]], elements: list[tuple[int, int, int]], doping: dict[int, float]) -> dict[str, Any]:
    edges: set[tuple[int, int]] = set()
    for a, b, c in elements:
        for i, j in [(a, b), (b, c), (c, a)]:
            edges.add((min(i, j), max(i, j)))
    points: list[tuple[float, float]] = []
    normal_lengths: list[float] = []
    for a, b in edges:
        da = doping.get(a, 0.0)
        db = doping.get(b, 0.0)
        if da == 0.0 or db == 0.0 or da * db > 0.0:
            continue
        xa, ya = nodes[a]
        xb, yb = nodes[b]
        t = abs(da) / (abs(da) + abs(db))
        points.append((xa + t * (xb - xa), ya + t * (yb - ya)))
        normal_lengths.append(math.hypot(xb - xa, yb - ya))
    if not points:
        zero_nodes = [node for node, value in doping.items() if value == 0.0 and node in nodes]
        points = [nodes[node] for node in zero_nodes]
    if not points:
        xs = [xy[0] for xy in nodes.values()]
        ys = [xy[1] for xy in nodes.values()]
        center = ((min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0)
        return {"center": center, "tangent": (0.0, 1.0), "normal": (1.0, 0.0), "dx_junction_m": None, "junction_points": 0}
    cx = sum(p[0] for p in points) / len(points)
    cy = sum(p[1] for p in points) / len(points)
    sxx = sum((p[0] - cx) ** 2 for p in points)
    syy = sum((p[1] - cy) ** 2 for p in points)
    sxy = sum((p[0] - cx) * (p[1] - cy) for p in points)
    if sxx == 0.0 and syy == 0.0 and sxy == 0.0:
        tangent = (0.0, 1.0)
    else:
        delta = math.sqrt((sxx - syy) * (sxx - syy) + 4.0 * sxy * sxy)
        lambda_max = 0.5 * (sxx + syy + delta)
        if abs(sxy) > 1.0e-40:
            tangent = (sxy, lambda_max - sxx)
        elif sxx >= syy:
            tangent = (1.0, 0.0)
        else:
            tangent = (0.0, 1.0)
        norm = math.hypot(tangent[0], tangent[1])
        tangent = (tangent[0] / norm, tangent[1] / norm)
    normal = (-tangent[1], tangent[0])
    normal_lengths.sort()
    dx = normal_lengths[len(normal_lengths) // 2] if normal_lengths else None
    if dx is None:
        normal_positions = sorted(set(round(project(point, (cx, cy), normal), 18) for point in nodes.values()))
        spacings = [abs(b - a) for a, b in zip(normal_positions, normal_positions[1:]) if abs(b - a) > 1.0e-15]
        dx = min(spacings) if spacings else None
    return {"center": (cx, cy), "tangent": tangent, "normal": normal, "dx_junction_m": dx, "junction_points": len(points)}


def project(point: tuple[float, float], origin: tuple[float, float], axis: tuple[float, float]) -> float:
    return (point[0] - origin[0]) * axis[0] + (point[1] - origin[1]) * axis[1]


def field_profile(nodes: dict[int, tuple[float, float]], elements: list[tuple[int, int, int]], potentials: dict[int, float], junction: dict[str, Any]) -> list[dict[str, float]]:
    center = junction["center"]
    tangent = junction["tangent"]
    normal = junction["normal"]
    band = max((junction.get("dx_junction_m") or 3.125e-8) * 2.5, 6.25e-8)
    rows: list[dict[str, float]] = []
    for tri in elements:
        if any(node not in potentials for node in tri):
            continue
        coords = [nodes[node] for node in tri]
        grad = triangle_gradient(coords, [potentials[node] for node in tri])
        if grad is None:
            continue
        centroid = (sum(x for x, _ in coords) / 3.0, sum(y for _, y in coords) / 3.0)
        tangent_pos = project(centroid, center, tangent)
        if abs(tangent_pos) > band:
            continue
        normal_pos = project(centroid, center, normal)
        field = math.hypot(grad[0], grad[1])
        rows.append({
            "normal_um": normal_pos * 1.0e6,
            "tangent_um": tangent_pos * 1.0e6,
            "x_um": centroid[0] * 1.0e6,
            "y_um": centroid[1] * 1.0e6,
            "electric_field_V_m": field,
        })
    rows.sort(key=lambda row: row["normal_um"])
    return rows


def gamma_factor(cfg: dict[str, Any]) -> float:
    phonon = float(cfg.get("phonon_energy_eV", DEFAULT_IMPACT["phonon_energy_eV"]))
    tref = float(cfg.get("reference_temperature_K", DEFAULT_IMPACT["reference_temperature_K"]))
    temp = float(cfg.get("temperature_K", DEFAULT_IMPACT["temperature_K"]))
    kb_ev = KB / Q
    denom = math.tanh(phonon / (2.0 * kb_ev * temp))
    if abs(denom) <= 0.0:
        return 1.0
    return math.tanh(phonon / (2.0 * kb_ev * tref)) / denom


def alpha_coefficients(field: float, cfg: dict[str, Any]) -> tuple[float, float]:
    model = str(cfg.get("model", "none"))
    minimum = float(cfg.get("minimum_field_V_m", 0.0))
    f = abs(field)
    if f <= 0.0 or f < minimum or model == "none":
        return 0.0, 0.0
    if model == "selberherr":
        an = float(cfg.get("electron_A_m_inv", DEFAULT_IMPACT["electron_A_m_inv"]))
        bn = float(cfg.get("electron_B_V_m", DEFAULT_IMPACT["electron_B_V_m"]))
        ap = float(cfg.get("hole_A_m_inv", DEFAULT_IMPACT["hole_A_m_inv"]))
        bp = float(cfg.get("hole_B_V_m", DEFAULT_IMPACT["hole_B_V_m"]))
        return an * math.exp(max(-700.0, -bn / f)), ap * math.exp(max(-700.0, -bp / f))
    if model == "van_overstraeten":
        g = gamma_factor(cfg)
        switch = float(cfg.get("switch_field_V_m", DEFAULT_IMPACT["switch_field_V_m"]))
        suffix = "low" if f < switch else "high"
        an = float(cfg.get(f"electron_a_{suffix}_m_inv", DEFAULT_IMPACT[f"electron_a_{suffix}_m_inv"]))
        bn = float(cfg.get(f"electron_b_{suffix}_V_m", DEFAULT_IMPACT[f"electron_b_{suffix}_V_m"]))
        ap = float(cfg.get(f"hole_a_{suffix}_m_inv", DEFAULT_IMPACT[f"hole_a_{suffix}_m_inv"]))
        bp = float(cfg.get(f"hole_b_{suffix}_V_m", DEFAULT_IMPACT[f"hole_b_{suffix}_V_m"]))
        return g * an * math.exp(max(-700.0, -bn * g / f)), g * ap * math.exp(max(-700.0, -bp * g / f))
    return 0.0, 0.0


def alpha_integral(profile: list[dict[str, float]], cfg: dict[str, Any]) -> float | None:
    if len(profile) < 2:
        return None
    enriched = []
    for row in profile:
        an, ap = alpha_coefficients(row["electric_field_V_m"], cfg)
        enriched.append((row["normal_um"] * 1.0e-6, an + ap))
    total = 0.0
    for (x0, a0), (x1, a1) in zip(enriched, enriched[1:]):
        dx = abs(x1 - x0)
        total += 0.5 * (a0 + a1) * dx
    return total


def point_at(points: list[tuple[float, float]], bias: float, log_abs: bool = True) -> float | None:
    ordered = sorted(points)
    for b, v in ordered:
        if abs(b - bias) <= 1.0e-8:
            return v
    for (b0, v0), (b1, v1) in zip(ordered, ordered[1:]):
        if min(b0, b1) <= bias <= max(b0, b1) and b0 != b1:
            t = (bias - b0) / (b1 - b0)
            if log_abs and v0 != 0.0 and v1 != 0.0:
                logv = (1.0 - t) * math.log10(abs(v0)) + t * math.log10(abs(v1))
                return math.copysign(10.0 ** logv, v1)
            return (1.0 - t) * v0 + t * v1
    return None


def load_curve(path: Path) -> list[tuple[float, float]]:
    points: dict[float, float] = {}
    for row in load_csv(path):
        bias = float_or_none(row.get("bias_V"))
        current = float_or_none(row.get("current_total_A_per_um"))
        if current is None:
            current = float_or_none(row.get("current_total"))
        if bias is None or current is None or current == 0.0:
            continue
        previous = points.get(bias)
        if previous is None or abs(current) > abs(previous):
            points[bias] = current
    return sorted(points.items())


def extract_impact_config(deck: dict[str, Any]) -> dict[str, Any]:
    solver = deck.get("solver", {})
    raw = solver.get("impact_ionization", {})
    cfg = dict(DEFAULT_IMPACT)
    if isinstance(raw, str):
        cfg["model"] = raw
    elif isinstance(raw, dict):
        cfg.update(raw)
    if "temperature_K" in solver and "temperature_K" not in raw:
        cfg["temperature_K"] = solver["temperature_K"]
    return cfg


def state_path_for_bias(out_dir: Path, bias: float) -> Path:
    return out_dir / "states" / f"accepted_state_bias_{bias_token(bias)}.csv"


def read_metadata_bias(path: Path) -> float | None:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    values: list[float] = []
    for field in metadata.get("fields", []):
        if field.get("name") == "ContactExternalVoltage":
            raw = field.get("raw_values", [])
            if raw:
                values.append(float(raw[0]))
    negative = [v for v in values if v < 0.0]
    if negative:
        return min(negative)
    return min(values) if values else None


def find_sentaurus_fields(roots: list[Path], target_biases: list[float]) -> tuple[dict[float, Path], list[str]]:
    found: dict[float, Path] = {}
    notes: list[str] = []
    candidates: list[tuple[float, Path]] = []
    for root in roots:
        if not root.exists():
            notes.append(f"missing sentaurus field root: {root}")
            continue
        for metadata in root.rglob("metadata.json"):
            fields_dir = metadata.parent / "fields"
            potential = fields_dir / "ElectrostaticPotential_region0.csv"
            if not potential.exists():
                continue
            bias = read_metadata_bias(metadata)
            if bias is not None:
                candidates.append((bias, fields_dir))
    for target in target_biases:
        exact = [(bias, fields) for bias, fields in candidates if abs(bias - target) <= 1.0e-8]
        if exact:
            found[target] = exact[0][1]
    return found, notes


def summarize_profile(profile: list[dict[str, float]], cfg: dict[str, Any]) -> dict[str, Any]:
    if not profile:
        return {
            "E_peak_V_m": None,
            "alpha_integral": None,
            "M_minus_1_estimate": None,
            "M_minus_1_note": "no profile",
        }
    peak = max(row["electric_field_V_m"] for row in profile)
    integral = alpha_integral(profile, cfg)
    driving_force = str(cfg.get("driving_force", "electric_field"))
    generation = str(cfg.get("generation", "local"))
    terminal_compatible = driving_force == "electric_field" and generation != "current_density"
    return {
        "E_peak_V_m": peak,
        "alpha_integral": integral,
        "M_minus_1_estimate": None if integral is None or not terminal_compatible else math.expm1(min(integral, 700.0)),
        "M_minus_1_note": (
            "junction-normal electrostatic alpha integral is an audit scalar, not a terminal M-1 estimate"
            if not terminal_compatible else "terminal-compatible local electric-field estimate"
        ),
    }


def source_distribution(rows: list[dict[str, str]], bias: float, junction: dict[str, Any]) -> dict[str, Any]:
    center = junction["center"]
    normal = junction["normal"]
    tangent = junction["tangent"]
    samples: list[tuple[float, float, float, float, float]] = []
    for row in rows:
        b = float_or_none(row.get("bias_V"))
        source = float_or_none(row.get("edge_source_integral"))
        if b is None or abs(b - bias) > 1.0e-8 or source is None:
            continue
        x = 0.5 * (float(row["x0_um"]) + float(row["x1_um"])) * 1.0e-6
        y = 0.5 * (float(row["y0_um"]) + float(row["y1_um"])) * 1.0e-6
        samples.append((abs(source), x, y, project((x, y), center, normal), project((x, y), center, tangent)))
    if not samples:
        return {"source_peak_integral": None, "source_peak_x_um": None, "source_peak_y_um": None, "source_hwhm_normal_um": None, "source_hwhm_tangent_um": None}
    peak = max(samples, key=lambda item: item[0])
    half = 0.5 * peak[0]
    above = [item for item in samples if item[0] >= half]
    normal_width = max(item[3] for item in above) - min(item[3] for item in above) if above else 0.0
    tangent_width = max(item[4] for item in above) - min(item[4] for item in above) if above else 0.0
    return {
        "source_peak_integral": peak[0],
        "source_peak_x_um": peak[1] * 1.0e6,
        "source_peak_y_um": peak[2] * 1.0e6,
        "source_hwhm_normal_um": 0.5 * normal_width * 1.0e6,
        "source_hwhm_tangent_um": 0.5 * tangent_width * 1.0e6,
    }


def safe_ratio(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or den == 0.0:
        return None
    return abs(num) / abs(den)


def classify_from_metrics(metrics: dict[str, Any]) -> str:
    e_ratio = metrics.get("e_peak_ratio")
    field_gain = metrics.get("sentaurus_field_alpha_integral_gain")
    seed_ratio = metrics.get("seed_current_ratio")
    m_ratio = metrics.get("vela_m_minus_1_over_reference")
    if e_ratio is not None and e_ratio < FIELD_RATIO_THRESHOLD:
        return "field_underestimate"
    if field_gain is not None and field_gain >= 1.5 and (e_ratio is None or e_ratio < 0.98):
        return "field_underestimate"
    if e_ratio is not None and e_ratio >= FIELD_RATIO_THRESHOLD and m_ratio is not None and m_ratio < 0.5:
        return "alpha_coefficient_gap"
    if seed_ratio is not None and seed_ratio < SEED_RATIO_THRESHOLD and m_ratio is not None and m_ratio >= 0.5:
        return "seed_leakage_gap"
    return "inconclusive"


def build_rows(out_dir: Path, deck: dict[str, Any], sentaurus_roots: list[Path], reference_curve_path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    mesh_path = Path(deck["mesh_file"])
    doping_path = Path(deck["node_doping_file"])
    nodes = load_nodes(mesh_path.with_name("nodes.csv") if mesh_path.name == "mesh.json" and (mesh_path.with_name("nodes.csv")).exists() else REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "nodes.csv")
    elements = load_elements(mesh_path.with_name("elements.csv") if mesh_path.name == "mesh.json" and (mesh_path.with_name("elements.csv")).exists() else REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "elements.csv")
    doping = load_doping(doping_path if doping_path.exists() else REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "doping.csv")
    junction = compute_junction_geometry(nodes, elements, doping)
    impact_cfg = extract_impact_config(deck)
    sentaurus_fields, notes = find_sentaurus_fields(sentaurus_roots, TARGET_BIASES)
    reference_curve = load_curve(reference_curve_path)
    vela_curve = load_curve(out_dir / "ionization_magnitude_curve.csv")
    source_rows = load_csv(out_dir / "sg_avalanche_edges.csv") if (out_dir / "sg_avalanche_edges.csv").exists() else []

    rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []
    for bias in TARGET_BIASES:
        state_path = state_path_for_bias(out_dir, bias)
        vela_profile: list[dict[str, float]] = []
        if state_path.exists():
            vela_profile = field_profile(nodes, elements, load_state_potential(state_path), junction)
        sent_profile: list[dict[str, float]] = []
        sentaurus_field_dir = sentaurus_fields.get(bias)
        if sentaurus_field_dir is not None:
            sent_profile = field_profile(nodes, elements, load_scalar_field(sentaurus_field_dir / "ElectrostaticPotential_region0.csv"), junction)
        vela_summary = summarize_profile(vela_profile, impact_cfg)
        sent_summary = summarize_profile(sent_profile, impact_cfg)
        for source, profile in [("vela", vela_profile), ("sentaurus", sent_profile)]:
            for row in profile:
                profile_rows.append({"bias_V": bias, "source": source, **row})
        i_vela = point_at(vela_curve, bias)
        i_ref = point_at(reference_curve, bias)
        i_ref_prev = point_at(reference_curve, bias + 0.25)
        ref_m = None if i_ref is None or i_ref_prev is None or i_ref_prev == 0.0 else max(abs(i_ref) / abs(i_ref_prev) - 1.0, 0.0)
        source_summary = source_distribution(source_rows, bias, junction)
        rows.append({
            "bias_V": bias,
            "vela_state_csv": str(state_path) if state_path.exists() else None,
            "sentaurus_field_dir": str(sentaurus_field_dir) if sentaurus_field_dir is not None else None,
            "E_peak_vela_V_m": vela_summary["E_peak_V_m"],
            "E_peak_sentaurus_V_m": sent_summary["E_peak_V_m"],
            "E_peak_vela_over_sentaurus": safe_ratio(vela_summary["E_peak_V_m"], sent_summary["E_peak_V_m"]),
            "alpha_integral_vela_field": vela_summary["alpha_integral"],
            "alpha_integral_sentaurus_field": sent_summary["alpha_integral"],
            "sentaurus_field_alpha_integral_gain": safe_ratio(sent_summary["alpha_integral"], vela_summary["alpha_integral"]),
            "M_minus_1_vela_field_estimate": vela_summary["M_minus_1_estimate"],
            "M_minus_1_sentaurus_field_estimate": sent_summary["M_minus_1_estimate"],
            "M_minus_1_vela_field_note": vela_summary["M_minus_1_note"],
            "M_minus_1_sentaurus_field_note": sent_summary["M_minus_1_note"],
            "M_minus_1_sentaurus_current_adjacent": ref_m,
            "vela_abs_current_A_per_um": None if i_vela is None else abs(i_vela),
            "sentaurus_abs_current_A_per_um": None if i_ref is None else abs(i_ref),
            "vela_over_sentaurus_current": safe_ratio(i_vela, i_ref),
            "dx_junction_um": None if junction["dx_junction_m"] is None else junction["dx_junction_m"] * 1.0e6,
            **source_summary,
        })
    write_csv(out_dir / "ionization_magnitude_points.csv", rows)
    write_csv(out_dir / "ionization_magnitude_profiles.csv", profile_rows)
    summary = {
        "junction": {
            "center_x_um": junction["center"][0] * 1.0e6,
            "center_y_um": junction["center"][1] * 1.0e6,
            "normal": junction["normal"],
            "tangent": junction["tangent"],
            "dx_junction_um": None if junction["dx_junction_m"] is None else junction["dx_junction_m"] * 1.0e6,
            "junction_crossing_edges": junction["junction_points"],
        },
        "impact_config": impact_cfg,
        "sentaurus_field_notes": notes,
        "missing_sentaurus_bias_fields_V": [bias for bias in TARGET_BIASES if bias not in sentaurus_fields],
    }
    return rows, summary


def choose_classification(rows: list[dict[str, Any]], seed_bias: float) -> tuple[str, dict[str, Any]]:
    comparable = [row for row in rows if row.get("E_peak_vela_over_sentaurus") is not None]
    focus = comparable[-1] if comparable else (rows[-1] if rows else {})
    seed_row = min(rows, key=lambda row: abs(row["bias_V"] - seed_bias)) if rows else {}
    ref_m = focus.get("M_minus_1_sentaurus_current_adjacent")
    vela_m = focus.get("M_minus_1_vela_field_estimate")
    metrics = {
        "bias_V": focus.get("bias_V"),
        "e_peak_ratio": focus.get("E_peak_vela_over_sentaurus"),
        "sentaurus_field_alpha_integral_gain": focus.get("sentaurus_field_alpha_integral_gain"),
        "seed_current_ratio": seed_row.get("vela_over_sentaurus_current"),
        "vela_m_minus_1_over_reference": safe_ratio(vela_m, ref_m),
    }
    return classify_from_metrics(metrics), metrics


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# PN2D BV ionization magnitude diagnostic",
        "",
        f"classification: `{payload['classification']}`",
        "",
        "## Key Metrics",
        "",
        "| metric | value |",
        "| --- | ---: |",
    ]
    for key, value in payload["classification_metrics"].items():
        lines.append(f"| {key} | {value} |")
    lines += [
        "",
        "## Bias Rows",
        "",
        "| bias V | E vela/sentaurus | alpha integral gain | Vela M-1 est | Sentaurus current M-1 | current ratio | dx junction um | source HWHM normal um |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["rows"]:
        lines.append(
            f"| {row['bias_V']} | {row.get('E_peak_vela_over_sentaurus')} | "
            f"{row.get('sentaurus_field_alpha_integral_gain')} | "
            f"{row.get('M_minus_1_vela_field_estimate')} | "
            f"{row.get('M_minus_1_sentaurus_current_adjacent')} | "
            f"{row.get('vela_over_sentaurus_current')} | "
            f"{row.get('dx_junction_um')} | {row.get('source_hwhm_normal_um')} |"
        )
    missing = payload["summary"].get("missing_sentaurus_bias_fields_V", [])
    lines += [
        "",
        "## Data Coverage",
        "",
        f"Missing Sentaurus bias fields: `{missing}`",
        "",
        "## Artifacts",
        "",
        f"- points_csv: `{payload['points_csv']}`",
        f"- profiles_csv: `{payload['profiles_csv']}`",
        f"- deck: `{payload['deck']}`",
        f"- curve_csv: `{payload['curve_csv']}`",
        f"- sg_avalanche_edges_csv: `{payload['sg_avalanche_edges_csv']}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def self_test() -> int:
    assert classify_from_metrics({
        "e_peak_ratio": 0.82,
        "sentaurus_field_alpha_integral_gain": 2.5,
        "seed_current_ratio": 0.9,
    }) == "field_underestimate"
    assert classify_from_metrics({
        "e_peak_ratio": 0.97,
        "vela_m_minus_1_over_reference": 0.2,
        "seed_current_ratio": 1.0,
    }) == "alpha_coefficient_gap"
    assert classify_from_metrics({
        "e_peak_ratio": 1.0,
        "vela_m_minus_1_over_reference": 0.8,
        "seed_current_ratio": 0.1,
    }) == "seed_leakage_gap"
    assert classify_from_metrics({
        "e_peak_ratio": 1.0,
        "seed_current_ratio": 0.1,
    }) == "inconclusive"
    grad = triangle_gradient([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], [1.0, 3.0, 4.0])
    assert grad is not None
    assert abs(grad[0] - 2.0) < 1.0e-12
    assert abs(grad[1] - 3.0) < 1.0e-12
    cfg = dict(DEFAULT_IMPACT)
    cfg["model"] = "van_overstraeten"
    an, ap = alpha_coefficients(5.0e7, cfg)
    assert an > 0.0 and ap > 0.0
    print("self-test passed")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--base-deck", type=Path, default=DEFAULT_DECK)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--reference-curve", type=Path, default=DEFAULT_REFERENCE_CURVE)
    parser.add_argument("--sentaurus-field-root", action="append", type=Path, default=[])
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--contacts", default="Anode,Cathode")
    parser.add_argument("--seed-bias", type=float, default=-16.0)
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        return self_test()
    out_dir = args.out_dir.resolve()
    if out_dir.exists() and not args.skip_run:
        shutil.rmtree(out_dir)
    contacts = [item.strip() for item in args.contacts.split(",") if item.strip()]
    deck_path = prepare_deck(args.base_deck, out_dir, contacts)
    if not args.skip_run:
        subprocess.run([str(args.runner.resolve()), "--config", str(deck_path)], check=True)
    deck = json.loads(deck_path.read_text(encoding="utf-8"))
    sentaurus_roots = [p.resolve() for p in args.sentaurus_field_root] or [p.resolve() for p in DEFAULT_SENT_FIELD_ROOTS]
    rows, summary = build_rows(out_dir, deck, sentaurus_roots, args.reference_curve.resolve())
    classification, metrics = choose_classification(rows, args.seed_bias)
    payload = {
        "classification": classification,
        "classification_metrics": metrics,
        "target_biases_V": TARGET_BIASES,
        "deck": str(deck_path),
        "curve_csv": str((out_dir / "ionization_magnitude_curve.csv").resolve()),
        "sg_avalanche_edges_csv": str((out_dir / "sg_avalanche_edges.csv").resolve()),
        "points_csv": str((out_dir / "ionization_magnitude_points.csv").resolve()),
        "profiles_csv": str((out_dir / "ionization_magnitude_profiles.csv").resolve()),
        "reference_curve": str(args.reference_curve.resolve()),
        "sentaurus_field_roots": [str(p) for p in sentaurus_roots],
        "summary": summary,
        "rows": rows,
    }
    summary_path = out_dir / "ionization_magnitude_report.json"
    md_path = out_dir / "ionization_magnitude_report.md"
    summary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(md_path, payload)
    print(json.dumps({
        "classification": classification,
        "classification_metrics": metrics,
        "report_json": str(summary_path),
        "report_md": str(md_path),
        "missing_sentaurus_bias_fields_V": summary.get("missing_sentaurus_bias_fields_V", []),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())