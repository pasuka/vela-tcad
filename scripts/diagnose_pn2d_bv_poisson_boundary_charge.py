#!/usr/bin/env python3
"""Compare PN2D BV doping, contact boundary values, and Poisson charge terms."""

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

STATE_FIELDS = {
    "psi": ("ElectrostaticPotential", "Potential"),
    "phin": ("eQuasiFermiPotential", "ElectronQuasiFermi"),
    "phip": ("hQuasiFermiPotential", "HoleQuasiFermi"),
    "n": ("eDensity", "Electrons"),
    "p": ("hDensity", "Holes"),
}

DOPING_FIELDS = [
    "source",
    "donors_cm3",
    "acceptors_cm3",
    "net_cm3",
    "node_count",
    "x_min_um",
    "x_max_um",
]

CONTACT_FIELDS = [
    "bias_V",
    "contact",
    "node_count",
    "applied_bias_V",
    "vela_net_doping_median_cm3",
    "vela_net_doping_mean_cm3",
    "expected_psi_builtin_V",
    "expected_psi_V",
    "expected_phin_V",
    "expected_phip_V",
    "vela_psi_median_V",
    "sentaurus_psi_median_V",
    "delta_psi_median_V",
    "vela_phin_median_V",
    "sentaurus_phin_median_V",
    "delta_phin_median_V",
    "vela_phip_median_V",
    "sentaurus_phip_median_V",
    "delta_phip_median_V",
    "vela_expected_psi_delta_median_V",
    "sentaurus_expected_psi_delta_median_V",
]

CHARGE_FIELDS = [
    "bias_V",
    "band",
    "x_min_um",
    "x_max_um",
    "node_count",
    "vela_donor_median_cm3",
    "sentaurus_donor_median_cm3",
    "delta_donor_median_cm3",
    "vela_acceptor_median_cm3",
    "sentaurus_acceptor_median_cm3",
    "delta_acceptor_median_cm3",
    "vela_net_doping_median_cm3",
    "sentaurus_net_doping_median_cm3",
    "delta_net_doping_median_cm3",
    "vela_electron_median_cm3",
    "sentaurus_electron_median_cm3",
    "log10_electron_ratio",
    "vela_hole_median_cm3",
    "sentaurus_hole_median_cm3",
    "log10_hole_ratio",
    "vela_charge_density_median_cm3",
    "sentaurus_charge_density_median_cm3",
    "delta_charge_density_median_cm3",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument(
        "--bands",
        default="left_p:0.0:0.25,pre_junction_p:0.7:0.9,junction:0.95:1.05,post_junction_n:1.1:1.3,right_n:1.5:1.8",
        help="Comma-separated name:xmin:xmax bands in um.",
    )
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.0e16)
    parser.add_argument("--bias-contact", default="Anode")
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


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


def read_mesh(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    data = json.loads(path.read_text())
    nodes = [
        {
            "id": int(node["id"]),
            "x_um": float(node["x"]),
            "y_um": float(node["y"]),
        }
        for node in data["nodes"]
    ]
    contacts = [
        {
            "name": str(contact.get("name", f"contact_{idx}")),
            "node_ids": [int(node_id) for node_id in contact.get("node_ids", [])],
        }
        for idx, contact in enumerate(data.get("contacts", []))
    ]
    return nodes, contacts


def load_doping_csv(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    seen = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            donors[node_id] = float(row["donors_cm3"])
            acceptors[node_id] = float(row["acceptors_cm3"])
            seen.add(node_id)
    missing = sorted(set(range(node_count)) - seen)
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return donors, acceptors


def load_state(sentaurus: Path, vtk: Path, node_count: int) -> dict[str, dict[str, list[float]]]:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    state: dict[str, dict[str, list[float]]] = {"vela": {}, "sentaurus": {}}
    for short_name, (sentaurus_name, vela_name) in STATE_FIELDS.items():
        if vela_name not in scalars:
            raise RuntimeError(f"{vtk} missing scalar {vela_name}")
        if len(scalars[vela_name]) != node_count:
            raise RuntimeError(f"{vela_name} length mismatch in {vtk}")
        values = scalars[vela_name]
        if short_name in {"n", "p"}:
            values = [value / 1.0e6 for value in values]
        state["vela"][short_name] = values
        state["sentaurus"][short_name] = load_sentaurus_scalar(sentaurus, sentaurus_name, node_count)
    return state


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def mean(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def log10_ratio(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or num <= 0.0 or den <= 0.0:
        return None
    return math.log10(num / den)


def equilibrium_electron_density(net_m3: float, ni_m3: float) -> float:
    return 0.5 * (net_m3 + math.sqrt(net_m3 * net_m3 + 4.0 * ni_m3 * ni_m3))


def builtin_potential(net_cm3: float, temperature_k: float, ni_m3: float) -> float:
    net_m3 = net_cm3 * 1.0e6
    neq = equilibrium_electron_density(net_m3, ni_m3)
    vt = K_B_OVER_Q * temperature_k
    if neq <= 0.0 or ni_m3 <= 0.0:
        return math.nan
    return vt * math.log(neq / ni_m3)


def values_at(values: list[float], node_ids: list[int]) -> list[float]:
    return [values[node_id] for node_id in node_ids]


def summarize_doping(source: str,
                     nodes: list[dict[str, Any]],
                     donors: list[float],
                     acceptors: list[float]) -> list[dict[str, Any]]:
    groups: dict[tuple[float, float], list[int]] = {}
    node_by_id = {int(node["id"]): node for node in nodes}
    for node in nodes:
        node_id = int(node["id"])
        key = (donors[node_id], acceptors[node_id])
        groups.setdefault(key, []).append(node_id)
    rows = []
    for (donor, acceptor), node_ids in sorted(groups.items(), key=lambda item: (item[0][0], item[0][1])):
        xs = [node_by_id[node_id]["x_um"] for node_id in node_ids]
        rows.append({
            "source": source,
            "donors_cm3": donor,
            "acceptors_cm3": acceptor,
            "net_cm3": donor - acceptor,
            "node_count": len(node_ids),
            "x_min_um": min(xs),
            "x_max_um": max(xs),
        })
    return rows


def load_sentaurus_contact_voltages(sentaurus: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    for path in sorted((sentaurus / "fields").glob("ContactExternalVoltage_region*.csv")):
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                values[int(row["node_id"])] = float(row["component0"])
    return values


def contact_rows(bias: float,
                 sentaurus: Path,
                 contacts: list[dict[str, Any]],
                 state: dict[str, dict[str, list[float]]],
                 vela_donors: list[float],
                 vela_acceptors: list[float],
                 bias_contact: str,
                 temperature_k: float,
                 ni_m3: float) -> list[dict[str, Any]]:
    sentaurus_contact_voltages = load_sentaurus_contact_voltages(sentaurus)
    rows = []
    for contact in contacts:
        node_ids = contact["node_ids"]
        if not node_ids:
            continue
        name = contact["name"]
        contact_voltage_values = [
            sentaurus_contact_voltages[node_id]
            for node_id in node_ids
            if node_id in sentaurus_contact_voltages
        ]
        applied_bias = (
            median(contact_voltage_values)
            if contact_voltage_values
            else (bias if name == bias_contact else 0.0)
        )
        net = [vela_donors[node_id] - vela_acceptors[node_id] for node_id in node_ids]
        net_mean = mean(net)
        net_median = median(net)
        expected_builtin = builtin_potential(net_mean or 0.0, temperature_k, ni_m3)
        expected_psi = (applied_bias or 0.0) + expected_builtin
        row: dict[str, Any] = {
            "bias_V": bias,
            "contact": name,
            "node_count": len(node_ids),
            "applied_bias_V": applied_bias,
            "vela_net_doping_median_cm3": net_median,
            "vela_net_doping_mean_cm3": net_mean,
            "expected_psi_builtin_V": expected_builtin,
            "expected_psi_V": expected_psi,
            "expected_phin_V": applied_bias,
            "expected_phip_V": applied_bias,
        }
        for short_name in ("psi", "phin", "phip"):
            v_med = median(values_at(state["vela"][short_name], node_ids))
            s_med = median(values_at(state["sentaurus"][short_name], node_ids))
            row[f"vela_{short_name}_median_V"] = v_med
            row[f"sentaurus_{short_name}_median_V"] = s_med
            row[f"delta_{short_name}_median_V"] = (
                v_med - s_med if v_med is not None and s_med is not None else None
            )
        row["vela_expected_psi_delta_median_V"] = (
            row["vela_psi_median_V"] - expected_psi
            if row["vela_psi_median_V"] is not None else None
        )
        row["sentaurus_expected_psi_delta_median_V"] = (
            row["sentaurus_psi_median_V"] - expected_psi
            if row["sentaurus_psi_median_V"] is not None else None
        )
        rows.append(row)
    return rows


def charge_rows(bias: float,
                nodes: list[dict[str, Any]],
                state: dict[str, dict[str, list[float]]],
                vela_donors: list[float],
                vela_acceptors: list[float],
                sentaurus_donors: list[float],
                sentaurus_acceptors: list[float],
                bands: list[tuple[str, float, float]]) -> list[dict[str, Any]]:
    rows = []
    for name, x_min, x_max in bands:
        ids = [
            int(node["id"]) for node in nodes
            if x_min <= float(node["x_um"]) <= x_max
        ]
        vdonor = values_at(vela_donors, ids)
        sdonor = values_at(sentaurus_donors, ids)
        vacc = values_at(vela_acceptors, ids)
        sacc = values_at(sentaurus_acceptors, ids)
        vnet = [donor - acceptor for donor, acceptor in zip(vdonor, vacc)]
        snet = [donor - acceptor for donor, acceptor in zip(sdonor, sacc)]
        vn = values_at(state["vela"]["n"], ids)
        sn = values_at(state["sentaurus"]["n"], ids)
        vp = values_at(state["vela"]["p"], ids)
        sp = values_at(state["sentaurus"]["p"], ids)
        vcharge = [p - n + net for p, n, net in zip(vp, vn, vnet)]
        scharge = [p - n + net for p, n, net in zip(sp, sn, snet)]
        row = {
            "bias_V": bias,
            "band": name,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "node_count": len(ids),
            "vela_donor_median_cm3": median(vdonor),
            "sentaurus_donor_median_cm3": median(sdonor),
            "vela_acceptor_median_cm3": median(vacc),
            "sentaurus_acceptor_median_cm3": median(sacc),
            "vela_net_doping_median_cm3": median(vnet),
            "sentaurus_net_doping_median_cm3": median(snet),
            "vela_electron_median_cm3": median(vn),
            "sentaurus_electron_median_cm3": median(sn),
            "vela_hole_median_cm3": median(vp),
            "sentaurus_hole_median_cm3": median(sp),
            "vela_charge_density_median_cm3": median(vcharge),
            "sentaurus_charge_density_median_cm3": median(scharge),
        }
        for field in ("donor", "acceptor", "net_doping", "charge_density"):
            v_key = f"vela_{field}_median_cm3"
            s_key = f"sentaurus_{field}_median_cm3"
            row[f"delta_{field}_median_cm3"] = (
                row[v_key] - row[s_key]
                if row[v_key] is not None and row[s_key] is not None
                else None
            )
        row["log10_electron_ratio"] = log10_ratio(
            row["vela_electron_median_cm3"], row["sentaurus_electron_median_cm3"])
        row["log10_hole_ratio"] = log10_ratio(
            row["vela_hole_median_cm3"], row["sentaurus_hole_median_cm3"])
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
    nodes, contacts = read_mesh(args.mesh)
    node_count = len(nodes)
    biases = parse_biases(args.biases)
    bands = parse_bands(args.bands)
    vtks = discover_vela_vtks(args.vela_vtk_root)
    vela_donors, vela_acceptors = load_doping_csv(args.doping_csv, node_count)

    first_sentaurus = sentaurus_dir(args.sentaurus_root, biases[0])
    sentaurus_donors = load_sentaurus_scalar(first_sentaurus, "DonorConcentration", node_count)
    sentaurus_acceptors = load_sentaurus_scalar(first_sentaurus, "AcceptorConcentration", node_count)
    doping_summary = summarize_doping("vela", nodes, vela_donors, vela_acceptors)
    doping_summary.extend(summarize_doping(
        "sentaurus", nodes, sentaurus_donors, sentaurus_acceptors))

    all_contact_rows: list[dict[str, Any]] = []
    all_charge_rows: list[dict[str, Any]] = []
    summaries = []
    for bias in biases:
        vtk = vtks.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        sentaurus = sentaurus_dir(args.sentaurus_root, bias)
        state = load_state(sentaurus, vtk, node_count)
        bias_contact_rows = contact_rows(
            bias,
            sentaurus,
            contacts,
            state,
            vela_donors,
            vela_acceptors,
            args.bias_contact,
            args.temperature_k,
            args.material_ni_m3,
        )
        bias_charge_rows = charge_rows(
            bias,
            nodes,
            state,
            vela_donors,
            vela_acceptors,
            sentaurus_donors,
            sentaurus_acceptors,
            bands,
        )
        all_contact_rows.extend(bias_contact_rows)
        all_charge_rows.extend(bias_charge_rows)
        summaries.append({
            "bias_V": bias,
            "vtk": str(vtk),
            "sentaurus_dir": str(sentaurus),
            "contacts": bias_contact_rows,
            "charge_bands": bias_charge_rows,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "doping_distribution_summary.csv", doping_summary, DOPING_FIELDS)
    write_csv(args.out_dir / "contact_boundary_summary.csv", all_contact_rows, CONTACT_FIELDS)
    write_csv(args.out_dir / "charge_balance_bands.csv", all_charge_rows, CHARGE_FIELDS)
    (args.out_dir / "poisson_boundary_charge_summary.json").write_text(
        json.dumps(clean_json({
            "biases": summaries,
            "doping_groups": doping_summary,
            "doping_group_count": len(doping_summary),
            "contact_rows": len(all_contact_rows),
            "charge_rows": len(all_charge_rows),
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
