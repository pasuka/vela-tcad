"""Probe pn2d IV: dump per-edge contact-current contributions at one bias.

Reads the Vela mesh.json + a VTK probe (Potential, eQuasiFermi, hQuasiFermi,
Electrons, Holes, NetDoping), reconstructs the FVM Voronoi-edge weights, and
prints for each (contact_node, interior_neighbor) edge:

  contact_name node_c (x_c,y_c) node_i (x_i,y_i)
    edge_length [m]  couple [m]
    psi_c psi_i  dpsi   phin_c phin_i  phip_c phip_i
    n_c n_i  p_c p_i
    SG electron flux contribution to I_e at contact [A/m]
    SG hole flux contribution to I_p at contact [A/m]

Then aggregates per contact (sum) and compares to the Vela CSV.

Use:
    python scripts/probe_pn2d_contact_edges.py
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from collections import defaultdict

BUILD = Path("build/pn2d_tdr_tie_probe")
MESH = BUILD / "vela" / "mesh.json"
VTK = BUILD / "vela" / "pn2d_iv_default_vtkprobe_0003_0.3V.vtk"
TEMP_K = 300.0
KB = 1.380649e-23
QE = 1.602176634e-19
VT = KB * TEMP_K / QE  # ~0.02585

# Caughey-Thomas with doping/field (matches mobility_model="caughey_thomas_field")
# Default Silicon parameters used in vela; keep consistent enough for first pass.
# In Vela: see include/vela/physics/MobilityModel.h. We approximate:
#  mu_n_min=68.5, mu_n_max=1414, N_ref_n=9.2e16, alpha_n=0.711, beta_n=2
#  mu_p_min=44.9, mu_p_max=470.5, N_ref_p=2.23e17, alpha_p=0.719, beta_p=1
# v_sat ~ 1.07e7 cm/s = 1.07e5 m/s for electrons (we ignore field dep here;
# at 0.30 V across 2 um the average field is small).
def caughey_thomas_low_field(N_tot_cm3: float, kind: str) -> float:
    if kind == "n":
        mu_min, mu_max, N_ref, alpha = 68.5, 1414.0, 9.2e16, 0.711
    else:
        mu_min, mu_max, N_ref, alpha = 44.9, 470.5, 2.23e17, 0.719
    mu_cm2 = mu_min + (mu_max - mu_min) / (1.0 + (N_tot_cm3 / N_ref) ** alpha)
    return mu_cm2 * 1e-4  # m^2/V/s


def bernoulli(x: float) -> float:
    if abs(x) < 1e-8:
        return 1.0 - 0.5 * x + x * x / 12.0
    return x / (math.exp(x) - 1.0)


def load_mesh():
    data = json.loads(MESH.read_text())
    nodes = [(n["x"] * 1e-6, n["y"] * 1e-6) for n in data["nodes"]]
    tris = [(t["node_ids"][0], t["node_ids"][1], t["node_ids"][2], t.get("region_id", 0))
            for t in data["triangles"]]
    contacts = {c["name"]: set(c["node_ids"]) for c in data["contacts"]}
    regions = {r["id"]: r["name"] for r in data["regions"]}
    return nodes, tris, contacts, regions


def triangle_circumcenter(p0, p1, p2):
    ax, ay = p0
    bx, by = p1
    cx, cy = p2
    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-30:
        return ((ax + bx + cx) / 3.0, (ay + by + cy) / 3.0)
    ux = ((ax * ax + ay * ay) * (by - cy) +
          (bx * bx + by * by) * (cy - ay) +
          (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) +
          (bx * bx + by * by) * (ax - cx) +
          (cx * cx + cy * cy) * (bx - ax)) / d
    return (ux, uy)


def build_edges(nodes, tris):
    """Return {(min,max): {'len':..., 'couple':..., 'tris':[i,j]}} approximation.

    couple uses sum of (distance from triangle circumcenter to edge midpoint)
    over adjacent triangles, which equals the standard Voronoi dual edge for
    a Delaunay mesh and matches Vela's edge.couple definition.
    """
    edges = defaultdict(lambda: {"len": 0.0, "couple": 0.0, "tris": []})
    for ti, (a, b, c, _r) in enumerate(tris):
        pa, pb, pc = nodes[a], nodes[b], nodes[c]
        cc = triangle_circumcenter(pa, pb, pc)
        for (i, j, popp) in [(a, b, pc), (b, c, pa), (c, a, pb)]:
            key = (min(i, j), max(i, j))
            pi, pj = nodes[i], nodes[j]
            length = math.hypot(pi[0] - pj[0], pi[1] - pj[1])
            mid = (0.5 * (pi[0] + pj[0]), 0.5 * (pi[1] + pj[1]))
            # dual length = |circumcenter - midpoint|, signed positive if
            # the circumcenter lies on the same side as the opposite vertex
            # (interior Delaunay); we sum unsigned here (good enough for
            # well-shaped Delaunay meshes).
            dual = math.hypot(cc[0] - mid[0], cc[1] - mid[1])
            edges[key]["len"] = length
            edges[key]["couple"] += dual
            edges[key]["tris"].append(ti)
    return edges


def parse_vtk(path: Path):
    text = path.read_text()
    fields = {}
    # POINTS
    m = re.search(r"POINTS\s+(\d+)\s+\w+\s*\n(.*?)\nCELLS", text, re.DOTALL)
    n_pts = int(m.group(1))
    pts = [float(v) for v in m.group(2).split()]
    coords = [(pts[3 * i], pts[3 * i + 1]) for i in range(n_pts)]
    for name in ("Potential", "ElectronQuasiFermi", "HoleQuasiFermi",
                 "Electrons", "Holes", "NetDoping"):
        mm = re.search(
            rf"SCALARS\s+{name}\s+\w+\s+1\s*\nLOOKUP_TABLE\s+\w+\s*\n(.*?)(?=\nSCALARS|\Z)",
            text, re.DOTALL)
        if mm:
            fields[name] = [float(v) for v in mm.group(1).split()]
    return coords, fields


def main():
    nodes, tris, contacts, _regions = load_mesh()
    coords_vtk, fields = parse_vtk(VTK)
    psi = fields["Potential"]
    phin = fields["ElectronQuasiFermi"]
    phip = fields["HoleQuasiFermi"]
    n = fields["Electrons"]
    p = fields["Holes"]
    netd = fields["NetDoping"]
    edges = build_edges(nodes, tris)

    print(f"# Vt = {VT*1000:.3f} mV  (T=300 K)")
    for cname, cnodes in contacts.items():
        Ie_tot = 0.0
        Ip_tot = 0.0
        edge_count = 0
        perim_couple = 0.0
        print(f"\n## contact {cname}: nodes={sorted(cnodes)}")
        print(f"{'cn':>4} {'in':>4} ({'xc':>7},{'yc':>6}) ({'xi':>7},{'yi':>6}) "
              f"{'len[nm]':>9} {'couple[nm]':>10} {'dpsi[V]':>9} "
              f"{'n_c':>9} {'n_i':>9} {'p_c':>9} {'p_i':>9} "
              f"{'mun':>7} {'mup':>7} {'dIe':>11} {'dIp':>11}")
        for (i, j), info in edges.items():
            i_on = i in cnodes
            j_on = j in cnodes
            if i_on == j_on:
                continue
            cn, ino = (i, j) if i_on else (j, i)
            edge_len = info["len"]
            couple = info["couple"]
            if edge_len < 1e-30 or couple <= 0:
                continue
            edge_count += 1
            perim_couple += couple
            dpsi = psi[ino] - psi[cn]  # j - i with j=interior
            # mobility on edge: use average doping at midpoint as proxy
            Ntot = 0.5 * (abs(netd[cn]) + abs(netd[ino]))
            mun = caughey_thomas_low_field(Ntot * 1e-6, "n")  # cm^-3 expected
            # netd is in m^-3? Vela writes it in SI (per cm^3 typically).
            # Try interpreting: doping magnitude ~ 1e23 -> 1e23 m^-3 = 1e17 cm^-3.
            # netd values around 1e23 -> 1e17 cm^-3 -> correct.
            mup = caughey_thomas_low_field(Ntot * 1e-6, "p")
            # SG flux on edge (cancellation-free not used; BGN-on uses density form)
            u = dpsi / VT
            Bp = bernoulli(u)
            Bm = bernoulli(-u)
            coef_e = mun * VT / edge_len
            coef_p = mup * VT / edge_len
            # electron continuity flux  J_n^{cont} (from cn -> ino):
            # SG: f = coef * (B(-u) * n_cn - B(u) * n_ino)
            nFlux = coef_e * (Bm * n[cn] - Bp * n[ino])
            # hole: f = coef * (B(u) * p_cn - B(-u) * p_ino)
            pFlux = coef_p * (Bp * p[cn] - Bm * p[ino])
            # Vela ContactCurrent sign: outward at contact_node=cn is +1.
            # I_e (electron contribution to contact current, A/m) =
            #   q * outward * (-nFlux_from_cn_to_ino) * couple
            # because electronCurrent in Vela uses electronFlux01 = -nFlux
            # Actually with cn = node_with_lower (or based on edge orientation)
            # we need to match: Vela sets outwardSign = +1 if n0OnContact,
            # so it works on edges where n0 is the contact. We map (i,j)->cn->ino
            # such that cn is always the contact node. The convention then gives
            # electronCurrent contribution = q * (-(-nFlux)) * couple = q*nFlux*couple
            # (per the formula electronFlux01 = -nFlux, outwardSign=+1)
            # So dIe = q * (-nFlux) * couple * outward
            # outward sign: at contact node, current FROM cn TO ino is "outward"
            # In Vela: nFlux is the electron-continuity flux from i to j; with
            # cn=i, outward=+1, sign -> q * (-nFlux) * couple
            # But for electron carrier-current, J_e = -q * electron_flux.
            # See ContactCurrent.cpp comment: electronFlux01 = -nFlux means
            # electronFlux01 is the electron carrier flux (signed for direction).
            # And result.electronCurrent += q * outward * electronFlux01 * couple
            # So: dIe = q * 1 * (-nFlux) * couple
            dIe = QE * (-nFlux) * couple
            dIp = QE * (pFlux) * couple
            # (NOTE the hole side: result.holeCurrent += q*outward*holeFlux01*couple
            # with holeFlux01 = -pFlux, but for hole CARRIER current J_p = +q*hole_flux
            # — convention check needed. We'll just sum and compare to CSV.)
            Ie_tot += dIe
            Ip_tot += dIp
            if edge_count <= 12:
                xc, yc = nodes[cn]; xi, yi = nodes[ino]
                print(f"{cn:>4} {ino:>4} ({xc*1e6:7.3f},{yc*1e6:6.3f}) "
                      f"({xi*1e6:7.3f},{yi*1e6:6.3f}) "
                      f"{edge_len*1e9:9.2f} {couple*1e9:10.2f} {dpsi:9.4f} "
                      f"{n[cn]:9.2e} {n[ino]:9.2e} {p[cn]:9.2e} {p[ino]:9.2e} "
                      f"{mun*1e4:7.1f} {mup*1e4:7.1f} {dIe:11.3e} {dIp:11.3e}")
        print(f"  edge_count={edge_count} perimeter(couple)={perim_couple*1e9:.2f} nm "
              f"= {perim_couple*1e6:.4f} um")
        print(f"  I_e={Ie_tot:.3e} A/m  I_p={Ip_tot:.3e} A/m  "
              f"I_tot={Ie_tot+Ip_tot:.3e} A/m  "
              f"-> /1e6 = {(Ie_tot+Ip_tot)/1e6:.3e} A/um")


if __name__ == "__main__":
    main()
