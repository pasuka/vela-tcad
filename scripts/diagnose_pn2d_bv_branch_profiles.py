#!/usr/bin/env python3
"""Compare PN2D BV carrier/potential profiles across avalanche/no-impact branches."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


FIELD_MAP = {
    "psi": ("ElectrostaticPotential", "Potential", 1.0),
    "phin": ("eQuasiFermiPotential", "ElectronQuasiFermi", 1.0),
    "phip": ("hQuasiFermiPotential", "HoleQuasiFermi", 1.0),
    "electron": ("eDensity", "Electrons", 1.0e-6),
    "hole": ("hDensity", "Holes", 1.0e-6),
}

BRANCH_FIELDS = [
    "bias_V",
    "band",
    "branch",
    "node_count",
    "psi_median_V",
    "phin_median_V",
    "phip_median_V",
    "psi_minus_phin_median_V",
    "phip_minus_psi_median_V",
    "electron_median_cm3",
    "hole_median_cm3",
    "electron_log10_median_cm3",
    "hole_log10_median_cm3",
]

COMPARISON_FIELDS = [
    "bias_V",
    "band",
    "node_count",
    "avalanche_minus_sentaurus_psi_median_V",
    "noimpact_minus_sentaurus_psi_median_V",
    "avalanche_minus_noimpact_psi_median_V",
    "avalanche_minus_sentaurus_psi_minus_phin_median_V",
    "noimpact_minus_sentaurus_psi_minus_phin_median_V",
    "avalanche_minus_noimpact_psi_minus_phin_median_V",
    "avalanche_minus_sentaurus_phip_minus_psi_median_V",
    "noimpact_minus_sentaurus_phip_minus_psi_median_V",
    "avalanche_minus_noimpact_phip_minus_psi_median_V",
    "avalanche_log10_electron_over_sentaurus_median",
    "noimpact_log10_electron_over_sentaurus_median",
    "avalanche_log10_electron_over_noimpact_median",
    "avalanche_log10_hole_over_sentaurus_median",
    "noimpact_log10_hole_over_sentaurus_median",
    "avalanche_log10_hole_over_noimpact_median",
]

FOCUS_FIELDS = [
    "bias_V",
    "node_id",
    "x_um",
    "y_um",
    "branch",
    "psi_V",
    "phin_V",
    "phip_V",
    "psi_minus_phin_V",
    "phip_minus_psi_V",
    "electron_cm3",
    "hole_cm3",
    "electron_log10_cm3",
    "hole_log10_cm3",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-avalanche-vtk-root", type=Path, required=True)
    parser.add_argument("--vela-no-impact-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument(
        "--bands",
        default=(
            "left_p:0.0:0.25,pre_junction_p:0.7:0.9,junction:0.95:1.05,"
            "post_junction_n:1.1:1.3,right_n:1.5:1.8"
        ),
        help="Comma-separated name:xmin:xmax plateau bands in um.",
    )
    parser.add_argument("--focus-nodes", default="955,1089,351,986")
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def parse_bands(raw: str) -> list[tuple[str, float, float]]:
    bands = []
    for item in raw.split(","):
        if not item.strip():
            continue
        name, x_min, x_max = item.split(":", 2)
        bands.append((name.strip(), float(x_min), float(x_max)))
    return bands


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
        if (candidate / "fields").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def load_sentaurus_scalar(root: Path, name: str, node_count: int) -> list[float]:
    values = [math.nan for _ in range(node_count)]
    path = root / "fields" / f"{name}_region0.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            values[int(row["node_id"])] = float(row["component0"])
    missing = [idx for idx, value in enumerate(values) if not math.isfinite(value)]
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return values


def read_nodes(mesh: Path) -> list[dict[str, float]]:
    data = json.loads(mesh.read_text())
    return [
        {"id": int(node["id"]), "x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data["nodes"]
    ]


def load_vela_state(vtk: Path, node_count: int) -> dict[str, list[float]]:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    state: dict[str, list[float]] = {}
    for short_name, (_sentaurus_name, vela_name, scale) in FIELD_MAP.items():
        if vela_name not in scalars:
            raise RuntimeError(f"{vtk} missing scalar {vela_name}")
        values = scalars[vela_name]
        if len(values) != node_count:
            raise RuntimeError(f"{vtk}:{vela_name} length mismatch")
        state[short_name] = [value * scale for value in values]
    return state


def load_sentaurus_state(root: Path, node_count: int) -> dict[str, list[float]]:
    state: dict[str, list[float]] = {}
    for short_name, (sentaurus_name, _vela_name, _scale) in FIELD_MAP.items():
        state[short_name] = load_sentaurus_scalar(root, sentaurus_name, node_count)
    return state


def finite_log10(value: float) -> float | None:
    if value > 0.0 and math.isfinite(value):
        return math.log10(value)
    return None


def median(values: list[float]) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.median(finite) if finite else None


def median_log_ratio(candidate: list[float], reference: list[float]) -> float | None:
    values = []
    for cand, ref in zip(candidate, reference):
        if cand > 0.0 and ref > 0.0 and math.isfinite(cand) and math.isfinite(ref):
            values.append(math.log10(cand / ref))
    return statistics.median(values) if values else None


def band_node_ids(nodes: list[dict[str, float]], x_min: float, x_max: float) -> list[int]:
    return [int(node["id"]) for node in nodes if x_min <= node["x_um"] <= x_max]


def take(values: list[float], ids: list[int]) -> list[float]:
    return [values[idx] for idx in ids]


def branch_row(bias: float,
               band: str,
               branch: str,
               ids: list[int],
               state: dict[str, list[float]]) -> dict[str, Any]:
    psi = take(state["psi"], ids)
    phin = take(state["phin"], ids)
    phip = take(state["phip"], ids)
    electrons = take(state["electron"], ids)
    holes = take(state["hole"], ids)
    psi_minus_phin = [p - q for p, q in zip(psi, phin)]
    phip_minus_psi = [q - p for q, p in zip(phip, psi)]
    electron_median = median(electrons)
    hole_median = median(holes)
    return {
        "bias_V": bias,
        "band": band,
        "branch": branch,
        "node_count": len(ids),
        "psi_median_V": median(psi),
        "phin_median_V": median(phin),
        "phip_median_V": median(phip),
        "psi_minus_phin_median_V": median(psi_minus_phin),
        "phip_minus_psi_median_V": median(phip_minus_psi),
        "electron_median_cm3": electron_median,
        "hole_median_cm3": hole_median,
        "electron_log10_median_cm3": finite_log10(electron_median or 0.0),
        "hole_log10_median_cm3": finite_log10(hole_median or 0.0),
    }


def diff_median(field: str,
                ids: list[int],
                candidate: dict[str, list[float]],
                reference: dict[str, list[float]]) -> float | None:
    if field == "psi_minus_phin":
        values = [
            candidate["psi"][idx] - candidate["phin"][idx] -
            (reference["psi"][idx] - reference["phin"][idx])
            for idx in ids
        ]
    elif field == "phip_minus_psi":
        values = [
            candidate["phip"][idx] - candidate["psi"][idx] -
            (reference["phip"][idx] - reference["psi"][idx])
            for idx in ids
        ]
    else:
        values = [candidate[field][idx] - reference[field][idx] for idx in ids]
    return median(values)


def comparison_row(bias: float,
                   band: str,
                   ids: list[int],
                   states: dict[str, dict[str, list[float]]]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "bias_V": bias,
        "band": band,
        "node_count": len(ids),
    }
    pairs = [
        ("avalanche", "sentaurus"),
        ("noimpact", "sentaurus"),
        ("avalanche", "noimpact"),
    ]
    for candidate, reference in pairs:
        prefix = f"{candidate}_minus_{reference}"
        row[f"{prefix}_psi_median_V"] = diff_median(
            "psi", ids, states[candidate], states[reference])
        row[f"{prefix}_psi_minus_phin_median_V"] = diff_median(
            "psi_minus_phin", ids, states[candidate], states[reference])
        row[f"{prefix}_phip_minus_psi_median_V"] = diff_median(
            "phip_minus_psi", ids, states[candidate], states[reference])
    for carrier in ["electron", "hole"]:
        row[f"avalanche_log10_{carrier}_over_sentaurus_median"] = median_log_ratio(
            take(states["avalanche"][carrier], ids), take(states["sentaurus"][carrier], ids))
        row[f"noimpact_log10_{carrier}_over_sentaurus_median"] = median_log_ratio(
            take(states["noimpact"][carrier], ids), take(states["sentaurus"][carrier], ids))
        row[f"avalanche_log10_{carrier}_over_noimpact_median"] = median_log_ratio(
            take(states["avalanche"][carrier], ids), take(states["noimpact"][carrier], ids))
    return row


def focus_rows(bias: float,
               nodes: list[dict[str, float]],
               focus_node_ids: list[int],
               states: dict[str, dict[str, list[float]]]) -> list[dict[str, Any]]:
    node_by_id = {int(node["id"]): node for node in nodes}
    rows = []
    for node_id in focus_node_ids:
        if node_id not in node_by_id:
            continue
        node = node_by_id[node_id]
        for branch in ["avalanche", "noimpact", "sentaurus"]:
            state = states[branch]
            electron = state["electron"][node_id]
            hole = state["hole"][node_id]
            rows.append({
                "bias_V": bias,
                "node_id": node_id,
                "x_um": node["x_um"],
                "y_um": node["y_um"],
                "branch": branch,
                "psi_V": state["psi"][node_id],
                "phin_V": state["phin"][node_id],
                "phip_V": state["phip"][node_id],
                "psi_minus_phin_V": state["psi"][node_id] - state["phin"][node_id],
                "phip_minus_psi_V": state["phip"][node_id] - state["psi"][node_id],
                "electron_cm3": electron,
                "hole_cm3": hole,
                "electron_log10_cm3": finite_log10(electron),
                "hole_log10_cm3": finite_log10(hole),
            })
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
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
    nodes = read_nodes(args.mesh)
    node_count = len(nodes)
    biases = parse_float_list(args.biases)
    bands = parse_bands(args.bands)
    focus_node_ids = parse_int_list(args.focus_nodes)
    avalanche_vtks = discover_vela_vtks(args.vela_avalanche_vtk_root)
    noimpact_vtks = discover_vela_vtks(args.vela_no_impact_vtk_root)
    branch_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    focus: list[dict[str, Any]] = []
    summaries = []
    for bias in biases:
        key = bias_key(bias)
        if key not in avalanche_vtks:
            raise FileNotFoundError(f"missing avalanche Vela VTK for {bias:g} V")
        if key not in noimpact_vtks:
            raise FileNotFoundError(f"missing no-impact Vela VTK for {bias:g} V")
        states = {
            "avalanche": load_vela_state(avalanche_vtks[key], node_count),
            "noimpact": load_vela_state(noimpact_vtks[key], node_count),
            "sentaurus": load_sentaurus_state(sentaurus_dir(args.sentaurus_root, bias), node_count),
        }
        bias_comparison = []
        for band, x_min, x_max in bands:
            ids = band_node_ids(nodes, x_min, x_max)
            for branch in ["avalanche", "noimpact", "sentaurus"]:
                branch_rows.append(branch_row(bias, band, branch, ids, states[branch]))
            row = comparison_row(bias, band, ids, states)
            comparison_rows.append(row)
            bias_comparison.append(row)
        focus.extend(focus_rows(bias, nodes, focus_node_ids, states))
        summaries.append({
            "bias_V": bias,
            "avalanche_vtk": str(avalanche_vtks[key]),
            "noimpact_vtk": str(noimpact_vtks[key]),
            "comparisons": bias_comparison,
        })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "branch_profile_medians.csv", branch_rows, BRANCH_FIELDS)
    write_csv(args.out_dir / "branch_profile_comparison.csv", comparison_rows, COMPARISON_FIELDS)
    write_csv(args.out_dir / "branch_profile_focus_nodes.csv", focus, FOCUS_FIELDS)
    (args.out_dir / "branch_profile_summary.json").write_text(
        json.dumps(clean_json({
            "biases": summaries,
            "branch_rows": len(branch_rows),
            "comparison_rows": len(comparison_rows),
            "focus_rows": len(focus),
            "bands": [
                {"name": name, "x_min_um": x_min, "x_max_um": x_max}
                for name, x_min, x_max in bands
            ],
            "focus_nodes": focus_node_ids,
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
