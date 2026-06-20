#!/usr/bin/env python3
"""Diagnose PN2D quasi-Fermi absolute anchoring against Sentaurus fields."""

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


K_B_OVER_Q = 8.617333262145e-5

FIELD_MAP = {
    "psi": ("ElectrostaticPotential", "Potential", 1.0),
    "phin": ("eQuasiFermiPotential", "ElectronQuasiFermi", 1.0),
    "phip": ("hQuasiFermiPotential", "HoleQuasiFermi", 1.0),
    "electron": ("eDensity", "Electrons", 1.0e-6),
    "hole": ("hDensity", "Holes", 1.0e-6),
}

CONTACT_FIELDS = [
    "bias_V",
    "contact",
    "node_count",
    "vela_psi_median_V",
    "sentaurus_psi_median_V",
    "delta_psi_median_V",
    "vela_phin_median_V",
    "sentaurus_phin_median_V",
    "delta_phin_median_V",
    "vela_phip_median_V",
    "sentaurus_phip_median_V",
    "delta_phip_median_V",
    "vela_psi_minus_phin_median_V",
    "sentaurus_psi_minus_phin_median_V",
    "delta_psi_minus_phin_median_V",
    "electron_log10_ratio_median",
    "hole_log10_ratio_median",
]

BAND_FIELDS = [
    "bias_V",
    "band",
    "node_count",
    "delta_psi_median_V",
    "delta_phin_median_V",
    "delta_phip_median_V",
    "delta_psi_minus_phin_median_V",
    "delta_phip_minus_psi_median_V",
    "electron_log10_ratio_median",
    "electron_log10_ratio_from_exponent_median",
    "hole_log10_ratio_median",
    "hole_log10_ratio_from_exponent_median",
    "uniform_phin_shift_to_match_electron_median_V",
    "contact_phin_violation_if_uniform_shift_V",
]

FOCUS_FIELDS = [
    "bias_V",
    "node_id",
    "x_um",
    "y_um",
    "vela_psi_V",
    "sentaurus_psi_V",
    "vela_phin_V",
    "sentaurus_phin_V",
    "vela_phip_V",
    "sentaurus_phip_V",
    "delta_psi_minus_phin_V",
    "delta_phip_minus_psi_V",
    "electron_log10_ratio",
    "electron_log10_ratio_from_exponent",
    "hole_log10_ratio",
    "hole_log10_ratio_from_exponent",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-13.2")
    parser.add_argument(
        "--bands",
        default=(
            "left_p:0.0:0.25,pre_junction_p:0.7:0.9,junction:0.95:1.05,"
            "post_junction_n:1.1:1.3,right_n:1.5:1.8"
        ),
    )
    parser.add_argument("--focus-nodes", default="955,1089,351,986")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def parse_bands(raw: str) -> list[tuple[str, float, float]]:
    result = []
    for item in raw.split(","):
        if not item.strip():
            continue
        name, x_min, x_max = item.split(":", 2)
        result.append((name.strip(), float(x_min), float(x_max)))
    return result


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    result: dict[float, Path] = {}
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


def read_mesh(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    nodes = [
        {"id": int(node["id"]), "x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data["nodes"]
    ]
    contacts = [
        {
            "name": str(contact.get("name", "contact")),
            "node_ids": [int(node_id) for node_id in contact.get("node_ids", [])],
        }
        for contact in data.get("contacts", [])
    ]
    return nodes, contacts


def load_sentaurus_scalar(root: Path, name: str, node_count: int) -> list[float]:
    path = root / "fields" / f"{name}_region0.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    values = [math.nan for _ in range(node_count)]
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            values[int(row["node_id"])] = float(row["component0"])
    missing = [idx for idx, value in enumerate(values) if not math.isfinite(value)]
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return values


def load_sentaurus_state(root: Path, node_count: int) -> dict[str, list[float]]:
    return {
        short_name: load_sentaurus_scalar(root, sentaurus_name, node_count)
        for short_name, (sentaurus_name, _vela_name, _scale) in FIELD_MAP.items()
    }


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


def median(values: list[float]) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.median(finite) if finite else None


def take(values: list[float], ids: list[int]) -> list[float]:
    return [values[idx] for idx in ids]


def log10_ratio(candidate: float, reference: float) -> float | None:
    if candidate > 0.0 and reference > 0.0 and math.isfinite(candidate) and math.isfinite(reference):
        return math.log10(candidate / reference)
    return None


def median_log10_ratio(candidate: list[float], reference: list[float]) -> float | None:
    return median([
        value for cand, ref in zip(candidate, reference)
        for value in [log10_ratio(cand, ref)]
        if value is not None
    ])


def exponent_log10_ratio(delta_exponent_v: list[float], vt: float) -> float | None:
    return median([value / (vt * math.log(10.0)) for value in delta_exponent_v])


def delta_values(field: str,
                 ids: list[int],
                 vela: dict[str, list[float]],
                 sentaurus: dict[str, list[float]]) -> list[float]:
    return [vela[field][idx] - sentaurus[field][idx] for idx in ids]


def delta_psi_minus_phin(ids: list[int],
                         vela: dict[str, list[float]],
                         sentaurus: dict[str, list[float]]) -> list[float]:
    return [
        (vela["psi"][idx] - vela["phin"][idx]) -
        (sentaurus["psi"][idx] - sentaurus["phin"][idx])
        for idx in ids
    ]


def delta_phip_minus_psi(ids: list[int],
                         vela: dict[str, list[float]],
                         sentaurus: dict[str, list[float]]) -> list[float]:
    return [
        (vela["phip"][idx] - vela["psi"][idx]) -
        (sentaurus["phip"][idx] - sentaurus["psi"][idx])
        for idx in ids
    ]


def contact_shift_violation(contacts: list[dict[str, Any]],
                            vela: dict[str, list[float]],
                            sentaurus: dict[str, list[float]],
                            phin_shift: float | None) -> float | None:
    if phin_shift is None:
        return None
    values = []
    for contact in contacts:
        for node_id in contact["node_ids"]:
            values.append(abs((vela["phin"][node_id] + phin_shift) - sentaurus["phin"][node_id]))
    return max(values) if values else None


def contact_rows(bias: float,
                 contacts: list[dict[str, Any]],
                 vela: dict[str, list[float]],
                 sentaurus: dict[str, list[float]]) -> list[dict[str, Any]]:
    rows = []
    for contact in contacts:
        ids = contact["node_ids"]
        if not ids:
            continue
        d_psi_minus_phin = delta_psi_minus_phin(ids, vela, sentaurus)
        rows.append({
            "bias_V": bias,
            "contact": contact["name"],
            "node_count": len(ids),
            "vela_psi_median_V": median(take(vela["psi"], ids)),
            "sentaurus_psi_median_V": median(take(sentaurus["psi"], ids)),
            "delta_psi_median_V": median(delta_values("psi", ids, vela, sentaurus)),
            "vela_phin_median_V": median(take(vela["phin"], ids)),
            "sentaurus_phin_median_V": median(take(sentaurus["phin"], ids)),
            "delta_phin_median_V": median(delta_values("phin", ids, vela, sentaurus)),
            "vela_phip_median_V": median(take(vela["phip"], ids)),
            "sentaurus_phip_median_V": median(take(sentaurus["phip"], ids)),
            "delta_phip_median_V": median(delta_values("phip", ids, vela, sentaurus)),
            "vela_psi_minus_phin_median_V": median([
                vela["psi"][idx] - vela["phin"][idx] for idx in ids
            ]),
            "sentaurus_psi_minus_phin_median_V": median([
                sentaurus["psi"][idx] - sentaurus["phin"][idx] for idx in ids
            ]),
            "delta_psi_minus_phin_median_V": median(d_psi_minus_phin),
            "electron_log10_ratio_median": median_log10_ratio(
                take(vela["electron"], ids), take(sentaurus["electron"], ids)),
            "hole_log10_ratio_median": median_log10_ratio(
                take(vela["hole"], ids), take(sentaurus["hole"], ids)),
        })
    return rows


def band_node_ids(nodes: list[dict[str, Any]], x_min: float, x_max: float) -> list[int]:
    return [int(node["id"]) for node in nodes if x_min <= float(node["x_um"]) <= x_max]


def band_row(bias: float,
             band: str,
             ids: list[int],
             contacts: list[dict[str, Any]],
             vela: dict[str, list[float]],
             sentaurus: dict[str, list[float]],
             vt: float) -> dict[str, Any]:
    d_psi_minus_phin = delta_psi_minus_phin(ids, vela, sentaurus)
    d_phip_minus_psi = delta_phip_minus_psi(ids, vela, sentaurus)
    phin_shift = median(d_psi_minus_phin)
    return {
        "bias_V": bias,
        "band": band,
        "node_count": len(ids),
        "delta_psi_median_V": median(delta_values("psi", ids, vela, sentaurus)),
        "delta_phin_median_V": median(delta_values("phin", ids, vela, sentaurus)),
        "delta_phip_median_V": median(delta_values("phip", ids, vela, sentaurus)),
        "delta_psi_minus_phin_median_V": phin_shift,
        "delta_phip_minus_psi_median_V": median(d_phip_minus_psi),
        "electron_log10_ratio_median": median_log10_ratio(
            take(vela["electron"], ids), take(sentaurus["electron"], ids)),
        "electron_log10_ratio_from_exponent_median": exponent_log10_ratio(d_psi_minus_phin, vt),
        "hole_log10_ratio_median": median_log10_ratio(
            take(vela["hole"], ids), take(sentaurus["hole"], ids)),
        "hole_log10_ratio_from_exponent_median": exponent_log10_ratio(d_phip_minus_psi, vt),
        "uniform_phin_shift_to_match_electron_median_V": phin_shift,
        "contact_phin_violation_if_uniform_shift_V": contact_shift_violation(
            contacts, vela, sentaurus, phin_shift),
    }


def focus_rows(bias: float,
               nodes: list[dict[str, Any]],
               focus_node_ids: list[int],
               vela: dict[str, list[float]],
               sentaurus: dict[str, list[float]],
               vt: float) -> list[dict[str, Any]]:
    node_by_id = {int(node["id"]): node for node in nodes}
    rows = []
    for node_id in focus_node_ids:
        if node_id not in node_by_id:
            continue
        d_e = (vela["psi"][node_id] - vela["phin"][node_id]) - (
            sentaurus["psi"][node_id] - sentaurus["phin"][node_id])
        d_h = (vela["phip"][node_id] - vela["psi"][node_id]) - (
            sentaurus["phip"][node_id] - sentaurus["psi"][node_id])
        node = node_by_id[node_id]
        rows.append({
            "bias_V": bias,
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "vela_psi_V": vela["psi"][node_id],
            "sentaurus_psi_V": sentaurus["psi"][node_id],
            "vela_phin_V": vela["phin"][node_id],
            "sentaurus_phin_V": sentaurus["phin"][node_id],
            "vela_phip_V": vela["phip"][node_id],
            "sentaurus_phip_V": sentaurus["phip"][node_id],
            "delta_psi_minus_phin_V": d_e,
            "delta_phip_minus_psi_V": d_h,
            "electron_log10_ratio": log10_ratio(vela["electron"][node_id], sentaurus["electron"][node_id]),
            "electron_log10_ratio_from_exponent": d_e / (vt * math.log(10.0)),
            "hole_log10_ratio": log10_ratio(vela["hole"][node_id], sentaurus["hole"][node_id]),
            "hole_log10_ratio_from_exponent": d_h / (vt * math.log(10.0)),
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
    vt = K_B_OVER_Q * args.temperature_k
    nodes, contacts = read_mesh(args.mesh)
    node_count = len(nodes)
    vtks = discover_vela_vtks(args.vela_vtk_root)
    biases = parse_float_list(args.biases)
    bands = parse_bands(args.bands)
    focus_node_ids = parse_int_list(args.focus_nodes)

    all_contact_rows: list[dict[str, Any]] = []
    all_band_rows: list[dict[str, Any]] = []
    all_focus_rows: list[dict[str, Any]] = []
    summaries = []

    for bias in biases:
        key = bias_key(bias)
        if key not in vtks:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        vela = load_vela_state(vtks[key], node_count)
        sentaurus = load_sentaurus_state(sentaurus_dir(args.sentaurus_root, bias), node_count)

        bias_contact_rows = contact_rows(bias, contacts, vela, sentaurus)
        bias_band_rows = []
        for band, x_min, x_max in bands:
            ids = band_node_ids(nodes, x_min, x_max)
            if not ids:
                continue
            row = band_row(bias, band, ids, contacts, vela, sentaurus, vt)
            all_band_rows.append(row)
            bias_band_rows.append(row)
        all_contact_rows.extend(bias_contact_rows)
        all_focus_rows.extend(focus_rows(bias, nodes, focus_node_ids, vela, sentaurus, vt))
        summaries.append({
            "bias_V": bias,
            "vela_vtk": str(vtks[key]),
            "contact_rows": len(bias_contact_rows),
            "band_rows": len(bias_band_rows),
            "max_contact_phin_violation_if_uniform_shift_V": max(
                [
                    abs(float(row["contact_phin_violation_if_uniform_shift_V"]))
                    for row in bias_band_rows
                    if row["contact_phin_violation_if_uniform_shift_V"] is not None
                ],
                default=None,
            ),
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "qf_anchor_contact_summary.csv", all_contact_rows, CONTACT_FIELDS)
    write_csv(args.out_dir / "qf_anchor_band_summary.csv", all_band_rows, BAND_FIELDS)
    write_csv(args.out_dir / "qf_anchor_focus_nodes.csv", all_focus_rows, FOCUS_FIELDS)
    (args.out_dir / "qf_anchor_summary.json").write_text(
        json.dumps(clean_json({
            "biases": summaries,
            "bands": [
                {"name": name, "x_min_um": x_min, "x_max_um": x_max}
                for name, x_min, x_max in bands
            ],
            "focus_nodes": focus_node_ids,
            "contact_rows": len(all_contact_rows),
            "band_rows": len(all_band_rows),
            "focus_rows": len(all_focus_rows),
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
