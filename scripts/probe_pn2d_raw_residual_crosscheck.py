from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

Q = 1.602176634e-19
KB = 1.380649e-23
T0 = 300.0
VT = KB * T0 / Q
NI = 1.0e16
SI_MUN = 0.135
SI_MUP = 0.048


def bernoulli(x: float) -> float:
    ax = abs(x)
    if ax < 1.0e-8:
        return 1.0 - x / 2.0 + x * x / 12.0
    if x > 500.0:
        return x * math.exp(-x)
    if x < -500.0:
        return -x
    return x / math.expm1(x)


def caughey(mu_max: float, net_doping: float, mu_min: float, n_ref: float, alpha: float) -> float:
    mu_min = min(max(mu_min, 0.0), mu_max)
    return mu_min + (mu_max - mu_min) / (1.0 + (abs(net_doping) / n_ref) ** alpha)


def load_doping(path: Path) -> List[Dict[str, float]]:
    out: List[Dict[str, float]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            donors = float(row["donors_cm3"]) * 1.0e6
            acceptors = float(row["acceptors_cm3"]) * 1.0e6
            out.append({"donors": donors, "acceptors": acceptors, "net": donors - acceptors})
    return out


def edge_key(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def triangle_area(nodes: List[Dict[str, float]], ids: List[int]) -> float:
    a, b, c = (nodes[i] for i in ids)
    return 0.5 * abs((b["x"] - a["x"]) * (c["y"] - a["y"]) - (c["x"] - a["x"]) * (b["y"] - a["y"]))


def cotangent(nodes: List[Dict[str, float]], a: int, b: int, opp: int) -> float:
    na, nb, no = nodes[a], nodes[b], nodes[opp]
    ux, uy = na["x"] - no["x"], na["y"] - no["y"]
    vx, vy = nb["x"] - no["x"], nb["y"] - no["y"]
    cross = ux * vy - uy * vx
    if abs(cross) < 1.0e-30:
        return 0.0
    return (ux * vx + uy * vy) / abs(cross)


def build_geometry(mesh: Dict) -> Tuple[List[Dict[str, float]], List[float]]:
    nodes = [{"x": n["x"] * 1.0e-6, "y": n["y"] * 1.0e-6} for n in mesh["nodes"]]
    edges: List[Dict[str, float]] = []
    index: Dict[Tuple[int, int], int] = {}
    volumes = [0.0 for _ in nodes]

    for cell in mesh["triangles"]:
        ids = cell["node_ids"]
        for k in range(3):
            key = edge_key(ids[k], ids[(k + 1) % 3])
            if key not in index:
                a, b = key
                dx = nodes[b]["x"] - nodes[a]["x"]
                dy = nodes[b]["y"] - nodes[a]["y"]
                index[key] = len(edges)
                edges.append({"n0": a, "n1": b, "length": math.hypot(dx, dy), "couple": 0.0})

    for cell in mesh["triangles"]:
        ids = cell["node_ids"]
        area = triangle_area(nodes, ids)
        if area <= 1.0e-30:
            continue
        for node_id in ids:
            volumes[node_id] += area / 3.0
        for k in range(3):
            a, b, opp = ids[k], ids[(k + 1) % 3], ids[(k + 2) % 3]
            e = edges[index[edge_key(a, b)]]
            cot = cotangent(nodes, a, b, opp)
            local = 0.5 * cot * e["length"]
            if cot < 0.0:
                local = area / (3.0 * e["length"])
            e["couple"] += max(local, 0.0)

    return edges, volumes


def slotboom_ni(doping: List[Dict[str, float]], bgn: str) -> List[float]:
    if bgn == "none":
        return [NI for _ in doping]
    out: List[float] = []
    for row in doping:
        effective = max(row["donors"] + row["acceptors"], 0.0)
        if effective <= 0.0:
            out.append(NI)
            continue
        x = math.log(effective / 1.0e23)
        delta = max(9.0e-3 * (x + math.sqrt(x * x + 0.5)), 0.0)
        out.append(NI * math.exp(delta / (2.0 * VT)))
    return out


def recombination_rate(n: float, p: float, ni: float, mechanisms: set[str], taun: float, taup: float) -> float:
    excess = max(min(n * p - ni * ni, 1.0e60), -1.0e60)
    rate = 0.0
    if "srh" in mechanisms:
        den = taup * (n + ni) + taun * (p + ni)
        if abs(den) >= 1.0e-100:
            rate += excess / den
    if "auger" in mechanisms:
        rate += (2.8e-43 * max(min(n, 1.0e30), -1.0e30) + 9.9e-44 * max(min(p, 1.0e30), -1.0e30)) * excess
    return rate


def read_vtk_scalars(path: Path) -> Dict[str, List[float]]:
    scalars: Dict[str, List[float]] = {}
    lines = path.read_text().splitlines()
    npoints = 0
    point_data_start = 0
    for idx, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "POINTS":
            npoints = int(parts[1])
        if len(parts) >= 2 and parts[0] == "POINT_DATA":
            point_data_start = idx + 1
            break

    i = point_data_start
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            scalars[name] = [float(lines[i + k].strip()) for k in range(npoints)]
            i += npoints
        else:
            i += 1
    return scalars


def load_contact_currents(vela_csv: Path, contact: str, bias: float) -> Dict[str, float]:
    best = None
    with vela_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("current_contact", "") != contact:
                continue
            b = float(row["bias_V"])
            if abs(b - bias) > 1.0e-9:
                continue
            best = row
            break

    if best is None:
        raise RuntimeError(f"No row found in {vela_csv} for contact={contact} bias={bias}")

    return {
        "total": float(best["current_total_A_per_um"]),
        "electron": float(best["current_electron_A_per_um"]),
        "hole": float(best["current_hole_A_per_um"]),
    }


def top10_text(contact_nodes: set[int], residual: List[float]) -> str:
    contributions = [(-Q * residual[idx] * 1.0e-6, idx) for idx in contact_nodes]
    contributions.sort(key=lambda x: abs(x[0]), reverse=True)
    top = contributions[:10]
    return ";".join(f"{idx}:{val:.6e}" for val, idx in top)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--doping", required=True, type=Path)
    parser.add_argument("--vtk", required=True, type=Path)
    parser.add_argument("--bias", required=True, type=float)
    parser.add_argument("--iv-cathode", required=True, type=Path)
    parser.add_argument("--iv-anode", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bgn", choices=["none", "slotboom"], default="slotboom")
    parser.add_argument("--recombination", default="srh,auger")
    parser.add_argument("--taun", type=float, default=1.0e-7)
    parser.add_argument("--taup", type=float, default=1.0e-7)
    args = parser.parse_args()

    mobility = {
        "electron_mu_min": 0.00522,
        "electron_nref": 9.68e22,
        "electron_alpha": 0.68,
        "hole_mu_min": 0.00449,
        "hole_nref": 2.23e23,
        "hole_alpha": 0.70,
        "field": True,
    }

    mesh = json.loads(args.mesh.read_text())
    doping = load_doping(args.doping)
    fields = read_vtk_scalars(args.vtk)
    contacts = {c["name"]: set(c["node_ids"]) for c in mesh["contacts"]}

    psi = fields["Potential"]
    phin = fields["ElectronQuasiFermi"]
    phip = fields["HoleQuasiFermi"]
    n = fields["Electrons"]
    p = fields["Holes"]

    ni = slotboom_ni(doping, args.bgn)
    edges, volumes = build_geometry(mesh)

    residual_density_n = [0.0 for _ in n]
    residual_density_p = [0.0 for _ in p]
    residual_qf_n = [0.0 for _ in n]
    residual_qf_p = [0.0 for _ in p]
    residual_asm_n = [0.0 for _ in n]
    residual_asm_p = [0.0 for _ in p]

    for edge in edges:
        i, j = int(edge["n0"]), int(edge["n1"])
        if edge["length"] < 1.0e-30 or edge["couple"] <= 0.0:
            continue

        dpsi = psi[j] - psi[i]
        field = abs(dpsi / edge["length"])
        net = 0.5 * (doping[i]["net"] + doping[j]["net"])

        mun = caughey(SI_MUN, net, mobility["electron_mu_min"], mobility["electron_nref"], mobility["electron_alpha"])
        mup = caughey(SI_MUP, net, mobility["hole_mu_min"], mobility["hole_nref"], mobility["hole_alpha"])
        if mobility["field"]:
            mun = mun / math.sqrt(1.0 + (mun * field / 1.0e5) ** 2)
            mup = mup / math.sqrt(1.0 + (mup * field / 1.0e5) ** 2)

        u = dpsi / VT
        bp = bernoulli(u)
        bm = bernoulli(-u)

        coef_n = mun * VT / edge["length"]
        coef_p = mup * VT / edge["length"]
        coef_n_c = coef_n * edge["couple"]
        coef_p_c = coef_p * edge["couple"]

        density_n_res = coef_n_c * (bm * n[i] - bp * n[j])
        density_p_res = coef_p_c * (bp * p[i] - bm * p[j])

        qf_n_res = coef_n_c * bp * ni[i] * math.exp(max(min(psi[j] / VT, 500.0), -500.0)) * (
            math.exp(max(min(-phin[i] / VT, 500.0), -500.0)) - math.exp(max(min(-phin[j] / VT, 500.0), -500.0))
        )
        qf_p_res = coef_p_c * bp * ni[i] * math.exp(max(min(-psi[i] / VT, 500.0), -500.0)) * (
            math.exp(max(min(phip[i] / VT, 500.0), -500.0)) - math.exp(max(min(phip[j] / VT, 500.0), -500.0))
        )

        asm_n_res = qf_n_res if ni[i] == ni[j] else density_n_res
        asm_p_res = qf_p_res if ni[i] == ni[j] else density_p_res

        residual_density_n[i] += density_n_res
        residual_density_n[j] -= density_n_res
        residual_density_p[i] += density_p_res
        residual_density_p[j] -= density_p_res

        residual_qf_n[i] += qf_n_res
        residual_qf_n[j] -= qf_n_res
        residual_qf_p[i] += qf_p_res
        residual_qf_p[j] -= qf_p_res

        residual_asm_n[i] += asm_n_res
        residual_asm_n[j] -= asm_n_res
        residual_asm_p[i] += asm_p_res
        residual_asm_p[j] -= asm_p_res

    mechanisms = {item.strip() for item in args.recombination.split(",") if item.strip() and item.strip() != "none"}
    for node_id in range(len(n)):
        rate = recombination_rate(n[node_id], p[node_id], ni[node_id], mechanisms, args.taun, args.taup)
        source = rate * volumes[node_id]
        residual_density_n[node_id] += source
        residual_density_p[node_id] += source
        residual_qf_n[node_id] += source
        residual_qf_p[node_id] += source
        residual_asm_n[node_id] += source
        residual_asm_p[node_id] += source

    contact_rows = {
        "Cathode": load_contact_currents(args.iv_cathode, "Cathode", args.bias),
        "Anode": load_contact_currents(args.iv_anode, "Anode", args.bias),
    }

    out_rows: List[Dict[str, str]] = []

    for contact_name, node_set in contacts.items():
        if contact_name not in contact_rows:
            continue

        raw_e = -Q * sum(residual_asm_n[i] for i in node_set) * 1.0e-6
        raw_h = -Q * sum(residual_asm_p[i] for i in node_set) * 1.0e-6
        raw_t = raw_e + raw_h

        # Dirichlet-replaced continuity residual at Dirichlet nodes is zero by construction.
        replaced_e = 0.0
        replaced_h = 0.0
        replaced_t = 0.0

        curr = contact_rows[contact_name]
        entries = [
            ("electron", curr["electron"], raw_e, replaced_e, top10_text(node_set, residual_asm_n)),
            ("hole", curr["hole"], raw_h, replaced_h, top10_text(node_set, residual_asm_p)),
            ("total", curr["total"], raw_t, replaced_t, ""),
        ]

        for carrier, from_contact_current, from_raw, from_replaced, top10 in entries:
            out_rows.append(
                {
                    "bias_V": f"{args.bias:.12g}",
                    "contact": contact_name,
                    "carrier": carrier,
                    "current_from_contact_current": f"{from_contact_current:.16e}",
                    "current_from_raw_residual": f"{from_raw:.16e}",
                    "current_from_dirichlet_replaced_residual": f"{from_replaced:.16e}",
                    "difference": f"{(from_contact_current - from_raw):.16e}",
                    "node_contributions_top10": top10,
                }
            )

    # Two-terminal sums for quick conservation check.
    cath_total = next(float(r["current_from_raw_residual"]) for r in out_rows if r["contact"] == "Cathode" and r["carrier"] == "total")
    ano_total = next(float(r["current_from_raw_residual"]) for r in out_rows if r["contact"] == "Anode" and r["carrier"] == "total")
    out_rows.append(
        {
            "bias_V": f"{args.bias:.12g}",
            "contact": "Both",
            "carrier": "total_sum",
            "current_from_contact_current": f"{(contact_rows['Cathode']['total'] + contact_rows['Anode']['total']):.16e}",
            "current_from_raw_residual": f"{(cath_total + ano_total):.16e}",
            "current_from_dirichlet_replaced_residual": f"{0.0:.16e}",
            "difference": f"{(contact_rows['Cathode']['total'] + contact_rows['Anode']['total'] - (cath_total + ano_total)):.16e}",
            "node_contributions_top10": "",
        }
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "bias_V",
                "contact",
                "carrier",
                "current_from_contact_current",
                "current_from_raw_residual",
                "current_from_dirichlet_replaced_residual",
                "difference",
                "node_contributions_top10",
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)


if __name__ == "__main__":
    main()
