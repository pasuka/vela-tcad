#!/usr/bin/env python3
"""Uniform 1->4 triangular mesh refinement for a Vela DeviceMesh JSON, with
linear node-doping interpolation. Used to test whether the PN2D BV -13.2V
current deficit (0.73x) is a depletion-region Scharfetter-Gummel discretization
artifact: if halving every edge moves Vela's converged current toward Sentaurus,
the residual is discretization; if the current is mesh-converged, it is not.

Each triangle (a,b,c) is split into four by inserting the three edge midpoints
mab, mbc, mca:  (a,mab,mca), (mab,b,mbc), (mca,mbc,c), (mab,mbc,mca).
Edge midpoints are shared between neighbouring triangles (deduplicated), so the
result is conforming. Node doping (donors/acceptors) at a midpoint is the linear
average of its two edge endpoints. A midpoint of an edge whose two endpoints both
belong to the same contact joins that contact.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", type=Path, required=True)
    p.add_argument("--doping-csv", type=Path, required=True)
    p.add_argument("--out-mesh", type=Path, required=True)
    p.add_argument("--out-doping-csv", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    mesh = json.loads(args.mesh.read_text())
    nodes = {int(n["id"]): (float(n["x"]), float(n["y"])) for n in mesh["nodes"]}

    doping: dict[int, tuple[float, float]] = {}
    with args.doping_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            doping[int(row["node_id"])] = (
                float(row["donors_cm3"]), float(row["acceptors_cm3"]))

    next_id = max(nodes) + 1
    midpoint_cache: dict[tuple[int, int], int] = {}

    def midpoint(a: int, b: int) -> int:
        nonlocal next_id
        key = (a, b) if a <= b else (b, a)
        if key in midpoint_cache:
            return midpoint_cache[key]
        ax, ay = nodes[a]
        bx, by = nodes[b]
        mid = next_id
        next_id += 1
        nodes[mid] = (0.5 * (ax + bx), 0.5 * (ay + by))
        da, aa = doping[a]
        db, ab = doping[b]
        doping[mid] = (0.5 * (da + db), 0.5 * (aa + ab))
        midpoint_cache[key] = mid
        return mid

    new_triangles: list[dict[str, Any]] = []

    def add_tri(region_id: int, a: int, b: int, c: int) -> None:
        new_triangles.append({
            "id": len(new_triangles),
            "region_id": region_id,
            "node_ids": [a, b, c],
        })

    for tri in mesh["triangles"]:
        a, b, c = (int(x) for x in tri["node_ids"])
        region_id = int(tri["region_id"])
        mab = midpoint(a, b)
        mbc = midpoint(b, c)
        mca = midpoint(c, a)
        add_tri(region_id, a, mab, mca)
        add_tri(region_id, mab, b, mbc)
        add_tri(region_id, mca, mbc, c)
        add_tri(region_id, mab, mbc, mca)

    # Contacts: a midpoint of an edge whose endpoints are both in a contact
    # joins that contact.
    new_contacts: list[dict[str, Any]] = []
    for contact in mesh["contacts"]:
        members = set(int(n) for n in contact["node_ids"])
        extra: set[int] = set()
        for (a, b), mid in midpoint_cache.items():
            if a in members and b in members:
                extra.add(mid)
        new_contacts.append({
            **{k: v for k, v in contact.items() if k != "node_ids"},
            "node_ids": sorted(members | extra),
        })

    all_tri_ids = [t["id"] for t in new_triangles]
    new_regions = [
        {**{k: v for k, v in r.items() if k != "cell_ids"}, "cell_ids": all_tri_ids}
        for r in mesh["regions"]
    ]

    out_mesh = {
        "nodes": [{"id": nid, "x": xy[0], "y": xy[1]}
                  for nid, xy in sorted(nodes.items())],
        "triangles": new_triangles,
        "regions": new_regions,
        "contacts": new_contacts,
    }
    args.out_mesh.parent.mkdir(parents=True, exist_ok=True)
    args.out_mesh.write_text(json.dumps(out_mesh))

    with args.out_doping_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "donors_cm3", "acceptors_cm3"])
        for nid in sorted(doping):
            d, a = doping[nid]
            writer.writerow([nid, repr(d), repr(a)])

    print(f"nodes {len(mesh['nodes'])} -> {len(out_mesh['nodes'])}")
    print(f"triangles {len(mesh['triangles'])} -> {len(out_mesh['triangles'])}")
    for c in out_mesh["contacts"]:
        print(f"contact {c['name']}: {len(c['node_ids'])} nodes")
    print(f"Wrote {args.out_mesh}")
    print(f"Wrote {args.out_doping_csv}")


if __name__ == "__main__":
    main()
