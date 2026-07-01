#!/usr/bin/env python3
"""Refine a Vela mesh.json 1->4 (add edge midpoints) + interpolate node doping.

Used for the (ii) BV mesh-refinement experiment: subdivide every triangle into
four by inserting edge midpoints, linearly interpolate per-node doping onto the
new midpoint nodes, and propagate contacts (a midpoint of a boundary edge whose
both endpoints belong to a contact joins that contact). Regions are preserved.

Inputs : <in-dir>/mesh.json, <in-dir>/doping.csv (node_id,donors_cm3,acceptors_cm3)
Outputs: <out-dir>/mesh.json, <out-dir>/doping.csv
Std lib only.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict


def load_doping(path):
    d = {}
    with open(path, newline="") as h:
        r = csv.DictReader(h)
        for row in r:
            d[int(row["node_id"])] = (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with open(os.path.join(args.in_dir, "mesh.json")) as h:
        mesh = json.load(h)
    doping = load_doping(os.path.join(args.in_dir, "doping.csv"))

    coord = {n["id"]: (n["x"], n["y"]) for n in mesh["nodes"]}
    next_id = max(coord) + 1

    # boundary edges = edges belonging to exactly one triangle
    edge_tri_count = defaultdict(int)
    for t in mesh["triangles"]:
        a, b, c = t["node_ids"]
        for u, v in ((a, b), (b, c), (c, a)):
            edge_tri_count[(min(u, v), max(u, v))] += 1
    boundary_edges = {e for e, n in edge_tri_count.items() if n == 1}

    midpoint = {}

    def get_mid(u, v):
        nonlocal next_id
        key = (min(u, v), max(u, v))
        if key in midpoint:
            return midpoint[key]
        xu, yu = coord[u]
        xv, yv = coord[v]
        mid = next_id
        next_id += 1
        coord[mid] = (0.5 * (xu + xv), 0.5 * (yu + yv))
        du, au = doping.get(u, (0.0, 0.0))
        dv, av = doping.get(v, (0.0, 0.0))
        doping[mid] = (0.5 * (du + dv), 0.5 * (au + av))
        midpoint[key] = mid
        return mid

    new_tris = []
    tid = 0
    for t in mesh["triangles"]:
        a, b, c = t["node_ids"]
        rid = t["region_id"]
        mab = get_mid(a, b)
        mbc = get_mid(b, c)
        mca = get_mid(c, a)
        for tri in ((a, mab, mca), (mab, b, mbc), (mca, mbc, c), (mab, mbc, mca)):
            new_tris.append({"id": tid, "region_id": rid, "node_ids": list(tri)})
            tid += 1

    # rebuild region cell_ids
    region_cells = defaultdict(list)
    for t in new_tris:
        region_cells[t["region_id"]].append(t["id"])
    new_regions = []
    for reg in mesh["regions"]:
        new_regions.append({
            "id": reg["id"],
            "name": reg["name"],
            "material": reg["material"],
            "cell_ids": region_cells[reg["id"]],
        })

    # propagate contacts: add midpoint of a boundary edge if both endpoints are in the contact
    new_contacts = []
    for ct in mesh["contacts"]:
        nodes = set(ct["node_ids"])
        add = []
        for e in boundary_edges:
            u, v = e
            if u in nodes and v in nodes and e in midpoint:
                add.append(midpoint[e])
        new_contacts.append({
            "id": ct["id"],
            "name": ct["name"],
            "region_id": ct["region_id"],
            "node_ids": sorted(nodes | set(add)),
        })

    new_nodes = [{"id": i, "x": coord[i][0], "y": coord[i][1]} for i in sorted(coord)]
    out_mesh = {
        "_comment": "1->4 refined from " + os.path.join(args.in_dir, "mesh.json"),
        "nodes": new_nodes,
        "triangles": new_tris,
        "regions": new_regions,
        "contacts": new_contacts,
    }

    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "mesh.json"), "w") as h:
        json.dump(out_mesh, h, indent=2)
    with open(os.path.join(args.out_dir, "doping.csv"), "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["node_id", "donors_cm3", "acceptors_cm3"])
        for i in sorted(coord):
            dn, ac = doping[i]
            w.writerow([i, repr(dn), repr(ac)])

    print(f"nodes {len(mesh['nodes'])} -> {len(new_nodes)}, "
          f"tris {len(mesh['triangles'])} -> {len(new_tris)}")
    for ct in new_contacts:
        print(f"  contact {ct['name']}: {len(ct['node_ids'])} nodes")


if __name__ == "__main__":
    main()
