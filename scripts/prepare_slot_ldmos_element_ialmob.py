"""Prepare isolated SLOT element-IALMob controls from its own geometry oracle.

The local-geometry control preserves the existing assembled edge couples and
mixed-Voronoi source volumes. Positive oracle cell coefficients only distribute
each existing couple among its adjacent silicon elements. This is an explicit
diagnostic decomposition, not a claim of Sentaurus geometry equivalence.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

try:
    from .prepare_slot_ldmos_ialmob_probe import extract_parameters
    from .audit_templates_ldmos_averagebox_node4492 import (
        coefficient_edge, measure_in_tdr_vertex_order, parse_debug_block,
    )
except ImportError:
    from prepare_slot_ldmos_ialmob_probe import extract_parameters
    from audit_templates_ldmos_averagebox_node4492 import (
        coefficient_edge, measure_in_tdr_vertex_order, parse_debug_block,
    )


def native_parameters(path: Path) -> dict:
    original = extract_parameters(path)
    result = {}
    for carrier, mass, other in (("electron", "me_over_m0", "mh_over_m0"),
                                 ("hole", "mh_over_m0", "me_over_m0")):
        p = dict(original[carrier])
        if p.pop("EnormMinimum", 0.) != 0.:
            raise ValueError("Nonzero EnormMinimum requires a separately audited branch")
        for name in ("nref_D", "nref_A", "cref_D", "cref_A"):
            if name in original:
                p[name] = original[name]
        if mass in original:
            p["mass"] = original[mass]
        if other in original:
            p["other_mass"] = original[other]
        # The unnamed isotropic block must remain identical for every family.
        result[carrier + "_parameters_cm"] = {f: dict(p) for f in ("100", "110", "111")}
    return result


def triangle_geometry(points: list) -> tuple[float, list, list]:
    cross = abs((points[1][0]-points[0][0])*(points[2][1]-points[0][1])
                -(points[1][1]-points[0][1])*(points[2][0]-points[0][0]))
    if not math.isfinite(cross) or cross <= 0.:
        raise ValueError("Nonfinite or degenerate triangle")
    cot = []
    for i in range(3):
        a, b, o = points[i], points[(i+1)%3], points[(i+2)%3]
        cot.append(((a[0]-o[0])*(b[0]-o[0])+(a[1]-o[1])*(b[1]-o[1]))/cross)
    area = cross/2
    if min(cot) < 0:
        shares = [area/4]*3
        shares[(cot.index(min(cot))+2)%3] = area/2
    else:
        shares = [(math.dist(points[i], points[(i+1)%3])**2*cot[i]
                   + math.dist(points[i], points[(i+2)%3])**2*cot[(i+2)%3])/8
                  for i in range(3)]
    return area, shares, [c/2 for c in cot]


def geometry_controls(mesh: dict, debug: Path, edge_probe: Path) -> tuple[dict, dict, dict]:
    measures = parse_debug_block(debug, "Measure")
    coefficients = parse_debug_block(debug, "Coefficients")
    materials = {r["id"]: r["material"] for r in mesh["regions"]}
    xy = {n["id"]: (n["x"], n["y"]) for n in mesh["nodes"]}
    raw, fallback, native_sum = defaultdict(float), defaultdict(float), defaultdict(float)
    local, native = [], []
    silicon_edge_count = defaultdict(int)
    area_sum = native_measure_sum = 0.
    mapping_checks = 0
    for c in mesh["triangles"]:
        ns = c["node_ids"]
        pts = [xy[n] for n in ns]
        area, shares, weights = triangle_geometry(pts)
        for k in range(3):
            edge = tuple(sorted((ns[k], ns[(k+1)%3])))
            raw[edge] += weights[k]
            fallback[edge] += area/(3*math.dist(pts[k], pts[(k+1)%3])**2)
        if materials[c["region_id"]] != "Si":
            continue
        cid = c["id"]
        by_edge = {tuple(sorted(coefficient_edge(ns, k)[:2])): v
                   for k, v in enumerate(coefficients[cid])}
        w = [by_edge[tuple(sorted((ns[k], ns[(k+1)%3])))] for k in range(3)]
        if len(w) != 3 or any(not math.isfinite(v) or v < 0 for v in w):
            raise ValueError("Invalid native edge coefficient")
        m = measure_in_tdr_vertex_order(measures[cid])
        # Independently check the debug-slot convention on asymmetric acute cells.
        if min(weights) > .01 and max(shares)-min(shares) > area*.05:
            if all(math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-20) for a,b in zip(m,shares)):
                mapping_checks += 1
        for k in range(3):
            edge = tuple(sorted((ns[k], ns[(k+1)%3])))
            native_sum[edge] += w[k]
            silicon_edge_count[edge] += 1
        row = dict(cell_id=cid, node_ids=ns, vertex_measure_m2=[v*1e-12 for v in m], edge_coefficients=w)
        native.append(row)
        local.append(dict(row, vertex_measure_m2=[v*1e-12 for v in shares]))
        area_sum += area
        native_measure_sum += sum(m)
    if mapping_checks < 10:
        raise ValueError("SLOT debug vertex mapping has insufficient independent acute-cell checks")
    # Use production floating-point couples, especially when nearly cancelling
    # cotangents make a Python reconstruction lose relative accuracy.
    with edge_probe.open(newline="",encoding="utf-8") as stream:
        legacy = {tuple(sorted((int(r["node0"]),int(r["node1"])))):
                  float(r["couple_m"])/float(r["length_m"]) for r in csv.DictReader(stream)}
    if legacy.keys() != raw.keys() or any(not math.isfinite(v) or v<0 for v in legacy.values()):
        raise ValueError("Production edge probe topology or finite positive geometry mismatch")
    for e, v in raw.items():
        expected = max(v,0.) if v >= -5e-13 else fallback[e]
        if not math.isclose(legacy[e],expected,rel_tol=1e-9,abs_tol=1e-12):
            raise ValueError(f"Production edge probe differs at {e}: {legacy[e]} vs {expected}")
    zero_mismatch = 0
    for c in local:
        ns = c["node_ids"]
        adjusted = []
        for k, w in enumerate(c["edge_coefficients"]):
            e = tuple(sorted((ns[k], ns[(k+1)%3])))
            if native_sum[e] == 0 and legacy[e] > 1e-12:
                zero_mismatch += 1
            # Preserve legacy edges absent from native transport support. The
            # diagnostic uses equal adjacent-cell weights and records them;
            # this fallback is never applied to the native geometry file.
            adjusted.append(w*legacy[e]/native_sum[e] if native_sum[e]
                            else legacy[e]/silicon_edge_count[e])
        c["edge_coefficients"] = adjusted
    envelope = dict(schema="vela.ialmob.transport_geometry.v1", node_count=len(mesh["nodes"]), cell_count=len(mesh["triangles"]))
    return (dict(envelope,cells=local), dict(envelope,cells=native),
            dict(silicon_cells=len(local), silicon_area_um2=area_sum,
                 native_measure_um2=native_measure_sum, mapping_acute_checks=mapping_checks,
                 zero_native_positive_legacy_edges=sum(v==0 and legacy[e]>0 for e,v in native_sum.items()),
                 zero_native_legacy_fallback_cell_entries=zero_mismatch,
                 transport_edges=len(native_sum),
                 native_vs_legacy_edge_coefficient_l1=sum(abs(v-legacy[e]) for e,v in native_sum.items())))


def prepare(base_path: Path, parameter_path: Path, debug: Path, edge_probe: Path, output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise ValueError("Refusing to overwrite an element IALMob run")
    config = json.loads(base_path.read_text(encoding="utf-8"))
    if config["mesh_geometry"] != dict(edge_coupling_policy="assembled_delaunay_fallback", node_volume_policy="mixed_voronoi", require_non_obtuse=False):
        raise ValueError("Expected frozen SLOT local geometry policy")
    mesh_path = Path(config["mesh_file"])
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    local, native, audit = geometry_controls(mesh, debug, edge_probe)
    params = native_parameters(parameter_path)
    output.mkdir(parents=True)
    for name, value in (("local_geometry.json",local),("native_geometry.json",native),("parameters.json",params)):
        (output/name).write_text(json.dumps(value, allow_nan=False),encoding="utf-8")
    old_prefix = str(base_path.resolve().parent)
    def relocate(value):
        if isinstance(value,str): return value.replace(old_prefix,str(output/"element_local"))
        if isinstance(value,list): return [relocate(v) for v in value]
        if isinstance(value,dict): return {k:relocate(v) for k,v in value.items()}
        return value
    cfg = relocate(copy.deepcopy(config))
    cfg["_comment"] = "SLOT element IALMob only: original unnamed isotropic parameters; legacy edge couples, mixed volumes, Boltzmann and SRH/Auger unchanged."
    mobility = cfg["solver"]["mobility"]
    mobility.pop("surface",None)
    mobility.update(model="ialmob",contact_electric_field_fallback=False,ialmob=dict(
        geometry_file=str(output/"local_geometry.json"),effective_electrodes=["source","drain","substrate"],
        crystal_x=[1,0,0],crystal_y=[0,1,0],high_field=True,reference_density_m3=0.,**params))
    (output/"element_local").mkdir()
    (output/"element_local/config.json").write_text(json.dumps(cfg,indent=2),encoding="utf-8")
    audit.update(scope="Diagnostic geometry decomposition, not native coupled qualification",
                 sources={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (base_path,parameter_path,debug,edge_probe,mesh_path)})
    (output/"manifest.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    return audit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config",type=Path,required=True)
    parser.add_argument("--parameters",type=Path,required=True)
    parser.add_argument("--oracle-debug",type=Path,required=True)
    parser.add_argument("--edge-probe",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.base_config,args.parameters,args.oracle_debug,args.edge_probe,args.output),indent=2))
