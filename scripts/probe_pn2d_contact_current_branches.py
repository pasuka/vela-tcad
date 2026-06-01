"""Compare density-SG and quasi-Fermi-SG contact current extraction.

This is a narrow pn2d diagnostic. It reads one generated Vela mesh and one
solution VTK, then recomputes contact currents with both SG branches used by
ContactCurrent. It intentionally implements only the silicon pn2d case.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

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


def load_doping(path: Path) -> list[dict[str, float]]:
    out: list[dict[str, float]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            donors = float(row["donors_cm3"]) * 1.0e6
            acceptors = float(row["acceptors_cm3"]) * 1.0e6
            out.append({"donors": donors, "acceptors": acceptors, "net": donors - acceptors})
    return out


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def triangle_area(nodes: list[dict], ids: list[int]) -> float:
    a, b, c = (nodes[i] for i in ids)
    return 0.5 * abs((b["x"] - a["x"]) * (c["y"] - a["y"]) - (c["x"] - a["x"]) * (b["y"] - a["y"]))


def cotangent(nodes: list[dict], a: int, b: int, opp: int) -> float:
    na, nb, no = nodes[a], nodes[b], nodes[opp]
    ux, uy = na["x"] - no["x"], na["y"] - no["y"]
    vx, vy = nb["x"] - no["x"], nb["y"] - no["y"]
    cross = ux * vy - uy * vx
    if abs(cross) < 1.0e-30:
        return 0.0
    return (ux * vx + uy * vy) / abs(cross)


def build_geometry(mesh: dict) -> tuple[list[dict], list[float]]:
    nodes = [{"x": n["x"] * 1.0e-6, "y": n["y"] * 1.0e-6} for n in mesh["nodes"]]
    edges: list[dict] = []
    index: dict[tuple[int, int], int] = {}
    volumes = [0.0 for _ in nodes]
    for cell in mesh["triangles"]:
        ids = cell["node_ids"]
        for k in range(3):
            key = edge_key(ids[k], ids[(k + 1) % 3])
            if key not in index:
                a, b = key
                dx, dy = nodes[b]["x"] - nodes[a]["x"], nodes[b]["y"] - nodes[a]["y"]
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


def slotboom_ni(doping: list[dict[str, float]], bgn: str) -> list[float]:
    if bgn == "none":
        return [NI for _ in doping]
    out: list[float] = []
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
        rate += (2.8e-43 * max(min(n, 1.0e30), -1.0e30) +
                 9.9e-44 * max(min(p, 1.0e30), -1.0e30)) * excess
    return rate


def read_vtk_scalars(path: Path) -> dict[str, list[float]]:
    scalars: dict[str, list[float]] = {}
    lines = path.read_text().splitlines()
    npoints = 0
    point_data_start = 0
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "POINTS":
            npoints = int(parts[1])
        if len(parts) >= 2 and parts[0] == "POINT_DATA":
            point_data_start = lines.index(line) + 1
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


def branch_currents(mesh: dict,
                    doping: list[dict[str, float]],
                    fields: dict[str, list[float]],
                    contact: str,
                    mobility: dict,
                    bgn: str,
                    recombination: set[str],
                    taun: float,
                    taup: float) -> dict[str, float]:
    contacts = {c["name"]: set(c["node_ids"]) for c in mesh["contacts"]}
    contact_nodes = contacts[contact]
    psi = fields["Potential"]
    phin = fields["ElectronQuasiFermi"]
    phip = fields["HoleQuasiFermi"]
    n = fields["Electrons"]
    p = fields["Holes"]
    ni = slotboom_ni(doping, bgn)
    edges, volumes = build_geometry(mesh)
    residual_density_n = [0.0 for _ in n]
    residual_density_p = [0.0 for _ in p]
    residual_qf_n = [0.0 for _ in n]
    residual_qf_p = [0.0 for _ in p]
    residual_asm_n = [0.0 for _ in n]
    residual_asm_p = [0.0 for _ in p]

    result = {k: 0.0 for k in ["density_e", "density_h", "qf_e", "qf_h"]}
    for edge in edges:
        n0_on = edge["n0"] in contact_nodes
        n1_on = edge["n1"] in contact_nodes
        if n0_on == n1_on or edge["length"] < 1.0e-30 or edge["couple"] <= 0.0:
            continue
        i, j = edge["n0"], edge["n1"]
        dpsi = psi[j] - psi[i]
        field = abs(dpsi / edge["length"])
        net = 0.5 * (doping[i]["net"] + doping[j]["net"])
        mun = caughey(SI_MUN, net, mobility["electron_mu_min"], mobility["electron_nref"], mobility["electron_alpha"])
        mup = caughey(SI_MUP, net, mobility["hole_mu_min"], mobility["hole_nref"], mobility["hole_alpha"])
        if mobility.get("field", False):
            mun = mun / math.sqrt(1.0 + (mun * field / 1.0e5) ** 2)
            mup = mup / math.sqrt(1.0 + (mup * field / 1.0e5) ** 2)
        u = dpsi / VT
        bp = bernoulli(u)
        bm = bernoulli(-u)
        coef_n = mun * VT / edge["length"]
        coef_p = mup * VT / edge["length"]
        coef_n_c = coef_n * edge["couple"]
        coef_p_c = coef_p * edge["couple"]
        density_n_flux = -(coef_n * (bm * n[i] - bp * n[j]))
        density_p_flux = -(coef_p * (bp * p[i] - bm * p[j]))
        qf_n_flux = -(coef_n * bp * ni[i] * math.exp(max(min(psi[j] / VT, 500.0), -500.0)) *
                       (math.exp(max(min(-phin[i] / VT, 500.0), -500.0)) -
                        math.exp(max(min(-phin[j] / VT, 500.0), -500.0))))
        qf_p_flux = -(coef_p * bp * ni[i] * math.exp(max(min(-psi[i] / VT, 500.0), -500.0)) *
                       (math.exp(max(min(phip[i] / VT, 500.0), -500.0)) -
                        math.exp(max(min(phip[j] / VT, 500.0), -500.0))))
        density_n_res = coef_n_c * (bm * n[i] - bp * n[j])
        density_p_res = coef_p_c * (bp * p[i] - bm * p[j])
        qf_n_res = coef_n_c * bp * ni[i] * math.exp(max(min(psi[j] / VT, 500.0), -500.0)) * (
            math.exp(max(min(-phin[i] / VT, 500.0), -500.0)) -
            math.exp(max(min(-phin[j] / VT, 500.0), -500.0)))
        qf_p_res = coef_p_c * bp * ni[i] * math.exp(max(min(-psi[i] / VT, 500.0), -500.0)) * (
            math.exp(max(min(phip[i] / VT, 500.0), -500.0)) -
            math.exp(max(min(phip[j] / VT, 500.0), -500.0)))
        asm_n_res = qf_n_res if ni[i] == ni[j] else density_n_res
        asm_p_res = qf_p_res if ni[i] == ni[j] else density_p_res
        for residual, flux in ((residual_density_n, density_n_res),
                               (residual_qf_n, qf_n_res),
                               (residual_asm_n, asm_n_res)):
            residual[i] += flux
            residual[j] -= flux
        for residual, flux in ((residual_density_p, density_p_res),
                               (residual_qf_p, qf_p_res),
                               (residual_asm_p, asm_p_res)):
            residual[i] += flux
            residual[j] -= flux
        sign = 1.0 if n0_on else -1.0
        scale = Q * sign * edge["couple"] * 1.0e-6
        result["density_e"] += scale * density_n_flux
        result["density_h"] += scale * density_p_flux
        result["qf_e"] += scale * qf_n_flux
        result["qf_h"] += scale * qf_p_flux

    for node_id in range(len(n)):
        rate = recombination_rate(n[node_id], p[node_id], ni[node_id], recombination, taun, taup)
        source = rate * volumes[node_id]
        residual_density_n[node_id] += source
        residual_density_p[node_id] += source
        residual_qf_n[node_id] += source
        residual_qf_p[node_id] += source
        residual_asm_n[node_id] += source
        residual_asm_p[node_id] += source

    result["density_total"] = result["density_e"] + result["density_h"]
    result["qf_total"] = result["qf_e"] + result["qf_h"]
    for name, rn, rp in (("density_residual", residual_density_n, residual_density_p),
                         ("qf_residual", residual_qf_n, residual_qf_p),
                         ("assembler_residual", residual_asm_n, residual_asm_p)):
        result[name + "_e"] = -Q * sum(rn[i] for i in contact_nodes) * 1.0e-6
        result[name + "_h"] = -Q * sum(rp[i] for i in contact_nodes) * 1.0e-6
        result[name + "_total"] = result[name + "_e"] + result[name + "_h"]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--doping", required=True, type=Path)
    parser.add_argument("--vtk", required=True, type=Path)
    parser.add_argument("--contact", required=True)
    parser.add_argument("--mobility", choices=["default_field", "promoted_bv"], default="default_field")
    parser.add_argument("--bgn", choices=["none", "slotboom"], default="none")
    parser.add_argument("--recombination", default="none")
    parser.add_argument("--taun", type=float, default=1.0e-7)
    parser.add_argument("--taup", type=float, default=1.0e-7)
    args = parser.parse_args()

    mobility = {
        "default_field": {
            "electron_mu_min": 0.00522, "electron_nref": 9.68e22, "electron_alpha": 0.68,
            "hole_mu_min": 0.00449, "hole_nref": 2.23e23, "hole_alpha": 0.70,
            "field": True,
        },
        "promoted_bv": {
            "electron_mu_min": 46.458e-4, "electron_nref": 9.68e22, "electron_alpha": 0.6052,
            "hole_mu_min": 39.961e-4, "hole_nref": 2.23e23, "hole_alpha": 0.623,
            "field": False,
        },
    }[args.mobility]
    mesh = json.loads(args.mesh.read_text())
    mechanisms = {item.strip() for item in args.recombination.split(",") if item.strip() and item.strip() != "none"}
    result = branch_currents(
        mesh,
        load_doping(args.doping),
        read_vtk_scalars(args.vtk),
        args.contact,
        mobility,
        args.bgn,
        mechanisms,
        args.taun,
        args.taup)
    print("branch,electron_A_per_um,hole_A_per_um,total_A_per_um")
    print(f"density,{result['density_e']:.16e},{result['density_h']:.16e},{result['density_total']:.16e}")
    print(f"qf,{result['qf_e']:.16e},{result['qf_h']:.16e},{result['qf_total']:.16e}")
    print(f"density_residual,{result['density_residual_e']:.16e},{result['density_residual_h']:.16e},{result['density_residual_total']:.16e}")
    print(f"qf_residual,{result['qf_residual_e']:.16e},{result['qf_residual_h']:.16e},{result['qf_residual_total']:.16e}")
    print(f"assembler_residual,{result['assembler_residual_e']:.16e},{result['assembler_residual_h']:.16e},{result['assembler_residual_total']:.16e}")


if __name__ == "__main__":
    main()
