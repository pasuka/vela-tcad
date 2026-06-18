#!/usr/bin/env python3
"""Compare PN2D BV potential/quasi-Fermi spatial profiles and plateau offsets."""

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
    "psi": ("ElectrostaticPotential", "Potential"),
    "phin": ("eQuasiFermiPotential", "ElectronQuasiFermi"),
    "phip": ("hQuasiFermiPotential", "HoleQuasiFermi"),
}

PROFILE_FIELDS = [
    "bias_V",
    "axis",
    "axis_value_um",
    "coordinate_um",
    "node_id",
    "x_um",
    "y_um",
    "vela_psi_V",
    "sentaurus_psi_V",
    "delta_psi_V",
    "vela_phin_V",
    "sentaurus_phin_V",
    "delta_phin_V",
    "vela_phip_V",
    "sentaurus_phip_V",
    "delta_phip_V",
    "vela_psi_minus_phin_V",
    "sentaurus_psi_minus_phin_V",
    "delta_psi_minus_phin_V",
    "vela_phip_minus_psi_V",
    "sentaurus_phip_minus_psi_V",
    "delta_phip_minus_psi_V",
]

PLATEAU_FIELDS = [
    "bias_V",
    "band",
    "x_min_um",
    "x_max_um",
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
    "vela_phip_minus_psi_median_V",
    "sentaurus_phip_minus_psi_median_V",
    "delta_phip_minus_psi_median_V",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument("--horizontal-y-um", default="0.015625,0.25")
    parser.add_argument("--vertical-x-um", default="0.1875,1.0,1.53125")
    parser.add_argument("--profile-tolerance-um", type=float, default=1.0e-6)
    parser.add_argument(
        "--bands",
        default="left_p:0.0:0.25,junction:0.9:1.1,right_n:1.5:1.8",
        help="Comma-separated name:xmin:xmax plateau bands in um.",
    )
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_biases(raw: str) -> list[float]:
    return parse_float_list(raw)


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


def load_state(sentaurus: Path, vtk: Path, node_count: int) -> dict[str, dict[str, list[float]]]:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    state: dict[str, dict[str, list[float]]] = {"vela": {}, "sentaurus": {}}
    for short_name, (sentaurus_name, vela_name) in FIELD_MAP.items():
        if vela_name not in scalars:
            raise RuntimeError(f"{vtk} missing scalar {vela_name}")
        if len(scalars[vela_name]) != node_count:
            raise RuntimeError(f"{vela_name} length mismatch")
        state["vela"][short_name] = scalars[vela_name]
        state["sentaurus"][short_name] = load_sentaurus_scalar(sentaurus, sentaurus_name, node_count)
    return state


def read_nodes(mesh: Path) -> list[dict[str, float]]:
    data = json.loads(mesh.read_text())
    return [
        {"id": int(node["id"]), "x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data["nodes"]
    ]


def value_row(bias: float,
              node: dict[str, float],
              state: dict[str, dict[str, list[float]]],
              axis: str,
              axis_value_um: float,
              coordinate_um: float) -> dict[str, Any]:
    node_id = int(node["id"])
    vpsi = state["vela"]["psi"][node_id]
    spsi = state["sentaurus"]["psi"][node_id]
    vphin = state["vela"]["phin"][node_id]
    sphin = state["sentaurus"]["phin"][node_id]
    vphip = state["vela"]["phip"][node_id]
    sphip = state["sentaurus"]["phip"][node_id]
    v_exp_n = vpsi - vphin
    s_exp_n = spsi - sphin
    v_exp_p = vphip - vpsi
    s_exp_p = sphip - spsi
    return {
        "bias_V": bias,
        "axis": axis,
        "axis_value_um": axis_value_um,
        "coordinate_um": coordinate_um,
        "node_id": node_id,
        "x_um": node["x_um"],
        "y_um": node["y_um"],
        "vela_psi_V": vpsi,
        "sentaurus_psi_V": spsi,
        "delta_psi_V": vpsi - spsi,
        "vela_phin_V": vphin,
        "sentaurus_phin_V": sphin,
        "delta_phin_V": vphin - sphin,
        "vela_phip_V": vphip,
        "sentaurus_phip_V": sphip,
        "delta_phip_V": vphip - sphip,
        "vela_psi_minus_phin_V": v_exp_n,
        "sentaurus_psi_minus_phin_V": s_exp_n,
        "delta_psi_minus_phin_V": v_exp_n - s_exp_n,
        "vela_phip_minus_psi_V": v_exp_p,
        "sentaurus_phip_minus_psi_V": s_exp_p,
        "delta_phip_minus_psi_V": v_exp_p - s_exp_p,
    }


def profile_rows(bias: float,
                 nodes: list[dict[str, float]],
                 state: dict[str, dict[str, list[float]]],
                 horizontal_y: list[float],
                 vertical_x: list[float],
                 tolerance_um: float) -> list[dict[str, Any]]:
    rows = []
    for y_value in horizontal_y:
        selected = [
            node for node in nodes
            if abs(node["y_um"] - y_value) <= tolerance_um
        ]
        selected.sort(key=lambda node: (node["x_um"], node["y_um"]))
        for node in selected:
            rows.append(value_row(
                bias, node, state, "horizontal_y", y_value, node["x_um"]))
    for x_value in vertical_x:
        selected = [
            node for node in nodes
            if abs(node["x_um"] - x_value) <= tolerance_um
        ]
        selected.sort(key=lambda node: (node["y_um"], node["x_um"]))
        for node in selected:
            rows.append(value_row(
                bias, node, state, "vertical_x", x_value, node["y_um"]))
    return rows


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def plateau_rows(bias: float,
                 nodes: list[dict[str, float]],
                 state: dict[str, dict[str, list[float]]],
                 bands: list[tuple[str, float, float]]) -> list[dict[str, Any]]:
    rows = []
    for name, x_min, x_max in bands:
        ids = [
            int(node["id"]) for node in nodes
            if x_min <= node["x_um"] <= x_max
        ]
        vpsi = [state["vela"]["psi"][idx] for idx in ids]
        spsi = [state["sentaurus"]["psi"][idx] for idx in ids]
        vphin = [state["vela"]["phin"][idx] for idx in ids]
        sphin = [state["sentaurus"]["phin"][idx] for idx in ids]
        vphip = [state["vela"]["phip"][idx] for idx in ids]
        sphip = [state["sentaurus"]["phip"][idx] for idx in ids]
        v_exp_n = [psi - phin for psi, phin in zip(vpsi, vphin)]
        s_exp_n = [psi - phin for psi, phin in zip(spsi, sphin)]
        v_exp_p = [phip - psi for phip, psi in zip(vphip, vpsi)]
        s_exp_p = [phip - psi for phip, psi in zip(sphip, spsi)]
        row = {
            "bias_V": bias,
            "band": name,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "node_count": len(ids),
            "vela_psi_median_V": median(vpsi),
            "sentaurus_psi_median_V": median(spsi),
            "vela_phin_median_V": median(vphin),
            "sentaurus_phin_median_V": median(sphin),
            "vela_phip_median_V": median(vphip),
            "sentaurus_phip_median_V": median(sphip),
            "vela_psi_minus_phin_median_V": median(v_exp_n),
            "sentaurus_psi_minus_phin_median_V": median(s_exp_n),
            "vela_phip_minus_psi_median_V": median(v_exp_p),
            "sentaurus_phip_minus_psi_median_V": median(s_exp_p),
        }
        for field in ["psi", "phin", "phip"]:
            row[f"delta_{field}_median_V"] = (
                row[f"vela_{field}_median_V"] - row[f"sentaurus_{field}_median_V"]
                if row[f"vela_{field}_median_V"] is not None
                and row[f"sentaurus_{field}_median_V"] is not None
                else None
            )
        row["delta_psi_minus_phin_median_V"] = (
            row["vela_psi_minus_phin_median_V"] -
            row["sentaurus_psi_minus_phin_median_V"]
        )
        row["delta_phip_minus_psi_median_V"] = (
            row["vela_phip_minus_psi_median_V"] -
            row["sentaurus_phip_minus_psi_median_V"]
        )
        rows.append(row)
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
    vtks = discover_vela_vtks(args.vela_vtk_root)
    biases = parse_biases(args.biases)
    horizontal_y = parse_float_list(args.horizontal_y_um)
    vertical_x = parse_float_list(args.vertical_x_um)
    bands = parse_bands(args.bands)
    profile: list[dict[str, Any]] = []
    plateaus: list[dict[str, Any]] = []
    summaries = []
    for bias in biases:
        vtk = vtks.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        state = load_state(sentaurus_dir(args.sentaurus_root, bias), vtk, node_count)
        profile.extend(profile_rows(
            bias, nodes, state, horizontal_y, vertical_x, args.profile_tolerance_um))
        bias_plateaus = plateau_rows(bias, nodes, state, bands)
        plateaus.extend(bias_plateaus)
        summaries.append({
            "bias_V": bias,
            "vtk": str(vtk),
            "plateaus": bias_plateaus,
        })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "potential_profile_samples.csv", profile, PROFILE_FIELDS)
    write_csv(args.out_dir / "potential_plateau_offsets.csv", plateaus, PLATEAU_FIELDS)
    (args.out_dir / "potential_profile_summary.json").write_text(
        json.dumps(clean_json({
            "biases": summaries,
            "profile_rows": len(profile),
            "plateau_rows": len(plateaus),
            "horizontal_y_um": horizontal_y,
            "vertical_x_um": vertical_x,
            "bands": [
                {"name": name, "x_min_um": x_min, "x_max_um": x_max}
                for name, x_min, x_max in bands
            ],
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
