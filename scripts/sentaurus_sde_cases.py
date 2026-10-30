#!/usr/bin/env python3
"""Run the versioned offline SDE corpus; expected rejections are not support."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from jsonschema import Draft7Validator
from sentaurus_device_ir import MATERIAL_ALIASES
from sentaurus_mesh_exporters import read_doping_csv
from vela_sdevice.mesh_import import import_grid

REPO = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO / "configs/reference_sde_cases.json"
PRIMARY = ("device.grd", "device.dat", "device.msh", "mesh.json", "doping.csv")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def material_digest(materials: dict) -> str:
    return hashlib.sha256(json.dumps(materials, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_manifest(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    schema = json.loads((REPO / "schemas/vela.reference_sde_cases.v1.schema.json").read_text())
    Draft7Validator(schema).validate(data)
    ids = [case["id"] for case in data["cases"]]
    require(len(ids) == len(set(ids)), "duplicate case ids")
    for case in data["cases"]:
        source = (REPO / case["source"]).resolve()
        require(source.is_relative_to(REPO.resolve()), "source outside repository")
        require(source.is_file() and hashlib.sha256(source.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == case["source_sha256_lf"],
                f"{case['id']}: source hash mismatch")
        if case["expected"] == "generated":
            materials = case["contract"]["materials"]
            require(material_digest(materials) == case["contract"]["materials_sha256"],
                    f"{case['id']}: material contract hash mismatch")
            require(all(MATERIAL_ALIASES.get(k) == v for k, v in materials.items()),
                    f"{case['id']}: material mapping changed")
    return data


def close(a: float, b: float, label: str, *, relative: float = 1e-12) -> None:
    require(math.isfinite(a) and math.isfinite(b) and
            math.isclose(a, b, rel_tol=relative, abs_tol=1e-12), label)


def topology(mesh: dict) -> Counter:
    regions = {r["id"]: r["name"] for r in mesh["regions"]}
    return Counter((regions[t["region_id"]], tuple(sorted(t["node_ids"])))
                   for t in mesh["triangles"])


def mesh_with_edges(output: Path) -> dict:
    mesh = json.loads((output / "mesh.json").read_text())
    bundle = json.loads((output / "mesh_field_bundle.json").read_text())
    tags = {c["name"]: c["physical_tag"] for c in bundle["contacts"]}
    for contact in mesh["contacts"]:
        contact["edge_node_ids"] = [line["node_ids"] for line in bundle["lines"]
            if line["kind"] == "contact" and tags[contact["name"]] in line["physical_curve_tags"]]
    return mesh


def contacts(mesh: dict) -> dict:
    names = {r["id"]: r["name"] for r in mesh["regions"]}
    return {c["name"]: (names[c["region_id"]], set(c["node_ids"]),
            Counter(tuple(sorted(edge)) for edge in c["edge_node_ids"])) for c in mesh["contacts"]}


def compare_native(original: Path, imported: Path) -> None:
    a = mesh_with_edges(original)
    b = json.loads((imported / "mesh.json").read_text())
    require(len(a["nodes"]) == len(b["nodes"]), "native node count mismatch")
    for x, y in zip(a["nodes"], b["nodes"]):
        require(x["id"] == y["id"], "native node identity mismatch")
        for axis in ("x", "y"):
            close(x[axis], y[axis], "native coordinate mismatch")
    coords = {n["id"]: n for n in b["nodes"]}
    for triangle in b["triangles"]:
        x, y, z = [coords[i] for i in triangle["node_ids"]]
        require((y["x"] - x["x"]) * (z["y"] - x["y"]) -
                (y["y"] - x["y"]) * (z["x"] - x["x"]) > 0,
                "native triangle orientation mismatch")
    require(topology(a) == topology(b), "native region/topology mismatch")
    require({r["name"]: r["material"] for r in a["regions"]} ==
            {r["name"]: r["material"] for r in b["regions"] if r["cell_ids"]}, "native material mismatch")
    bundle = json.loads((original / "mesh_field_bundle.json").read_text())
    require({r["name"] for r in b.get("interfaces", [])} ==
            {r["name"] for r in bundle["interfaces"]}, "native interface identity mismatch")
    for interface in bundle["interfaces"]:
        actual = next(i for i in b["interfaces"] if i["name"] == interface["name"])
        expected = Counter(tuple(sorted(line["node_ids"])) for line in bundle["lines"]
            if line["kind"] == "interface" and interface["physical_tag"] in line["physical_curve_tags"])
        require(expected == Counter(tuple(sorted(pair)) for pair in actual["edge_node_ids"]),
                "native interface edges mismatch")
    require(contacts(a) == contacts(b), "native contact edges/nodes mismatch")
    fields_a = read_doping_csv(original / "doping.csv")
    with (imported / "doping.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    require([int(r["node_id"]) for r in rows] == [n["id"] for n in a["nodes"]],
            "native doping node identity mismatch")
    fields_b = {"DonorConcentration": [float(r["donors_cm3"]) for r in rows],
                "AcceptorConcentration": [float(r["acceptors_cm3"]) for r in rows],
                "DopingConcentration": [float(r["donors_cm3"]) - float(r["acceptors_cm3"]) for r in rows]}
    for name in ("DonorConcentration", "AcceptorConcentration", "DopingConcentration"):
        require(len(fields_a[name]) == len(fields_b[name]), "native field length mismatch")
        for x, y in zip(fields_a[name], fields_b[name]):
            close(x, y, f"native {name} mismatch")


def check_contract(output: Path, contract: dict) -> dict:
    mesh = mesh_with_edges(output)
    bundle = json.loads((output / "mesh_field_bundle.json").read_text())
    require(bundle["coordinate_unit"] == "um" and
            bundle["field_units"]["concentration"] == "cm^-3", "unit mismatch")
    require({m["name"] for m in bundle["materials"]} == set(contract["materials"]),
            "source material identity mismatch")
    require({r["name"]: r["material"] for r in mesh["regions"]} ==
            contract["regions"], "region/material contract mismatch")
    require(set(contacts(mesh)) == set(contract["contacts"]), "contact names mismatch")
    areas = Counter()
    nodes = {n["id"]: n for n in mesh["nodes"]}
    region_names = {r["id"]: r["name"] for r in mesh["regions"]}
    for tri in mesh["triangles"]:
        a, b, c = [nodes[i] for i in tri["node_ids"]]
        area = ((b["x"] - a["x"]) * (c["y"] - a["y"]) -
                (b["y"] - a["y"]) * (c["x"] - a["x"])) / 2
        require(area > 0, "nonpositive triangle area")
        areas[region_names[tri["region_id"]]] += area
    for name, area in contract["region_area_um2"].items():
        if name in contract.get('region_area_error_um2',{}):
            require(abs(areas[name]-area)<=contract['region_area_error_um2'][name],f'curved area mismatch: {name}')
        else:close(areas[name], area, f"region area mismatch: {name}")
    for i, value in enumerate((min(n["x"] for n in nodes.values()), min(n["y"] for n in nodes.values()),
                               max(n["x"] for n in nodes.values()), max(n["y"] for n in nodes.values()))):
        close(value, contract["bounds_um"][i], "device bounds mismatch")
    for name, edge_contract in contract["contacts"].items():
        contact = next(c for c in mesh["contacts"] if c["name"] == name)
        length = 0.0
        for edge in contact["edge_node_ids"]:
            a, b = [nodes[i] for i in edge]
            for node in (a, b):
                close(node[edge_contract["axis"]], edge_contract["coordinate_um"], "contact location mismatch")
            length += math.hypot(a["x"] - b["x"], a["y"] - b["y"])
        close(length, edge_contract["length_um"], "contact length mismatch")
        if "span_um" in edge_contract:
            axis = "x" if edge_contract["axis"] == "y" else "y"
            positions = [nodes[i][axis] for i in contact["node_ids"]]
            close(min(positions), edge_contract["span_um"][0], "contact span start mismatch")
            close(max(positions), edge_contract["span_um"][1], "contact span end mismatch")
            require(len(positions) >= 2, "contact has too few nodes")
    # Analytic pointwise doping, independent of the mesh backend or the exporter.
    bjt_report = None
    if contract["doping"]["kind"] == "bjt":
        from sentaurus_bjt_contract import check
        bjt_report = check(bundle, local=contract["doping"]["local_refined"])
        require(bjt_report["status"] == "passed", f"BJT contract failed: {bjt_report}")
        require(q_adaptive := bundle["provenance"].get("maxtransdiff"), "missing adaptive audit")
        require(q_adaptive["status"] == "passed" and q_adaptive["rounds"][-1]["refineable_nodes"] == 0,
                "adaptive termination contract failed")
        require(all(v["default_evaluations"] == 0 for v in bundle["provenance"]["general_evaluation"].values()),
                "BJT General used a default concentration")
    mos_report=None
    if contract['doping']['kind']=='mos':
        from sentaurus_mos_contract import check
        mos_report=check(bundle)
        require(mos_report['status']=='passed',f'MOS contract failed: {mos_report}')
        # Read the authored half through the native consumer on every repeat.
        import shutil
        from sentaurus_mesh_exporters import bundle_to_vela_mesh, write_doping_csv
        half=json.loads((output/'half_mesh_field_bundle.json').read_text())
        count=len(half['nodes'])
        require(all(a['x']==b['x'] and abs(a['y']-b['y'])<=1e-8 for a,b in zip(half['nodes'],bundle['nodes'][:count])),
                'full device does not retain authored half nodes')
        with tempfile.TemporaryDirectory(prefix='mos-half-native-') as temp:
            stage=Path(temp)/'half';stage.mkdir()
            (stage/'mesh_field_bundle.json').write_text(json.dumps(half))
            (stage/'mesh.json').write_text(json.dumps(bundle_to_vela_mesh(half)))
            write_doping_csv(stage/'doping.csv',half)
            for extension in ('grd','dat'):shutil.copy2(output/('half_device.'+extension),stage/('device.'+extension))
            import_grid(stage/'device.grd',Path(temp)/'native')
            compare_native(stage,Path(temp)/'native')
        mos_report['half_native_import']='passed'
    for node in mesh["nodes"] if bjt_report is None and mos_report is None else []:
        x, y = node["x"], node["y"]
        rule = contract["doping"]
        if rule["kind"] == "pn":
            donor = rule["concentration_cm3"] if x >= rule["junction_x_um"] - 1e-12 else 0.0
            acceptor = rule["concentration_cm3"] if x <= rule["junction_x_um"] + 1e-12 else 0.0
        else:
            donor = rule["concentration_cm3"] if y <= rule["silicon_top_um"] + 1e-12 else 0.0
            acceptor = 0.0
        for field, expected in (("DonorConcentration", donor), ("AcceptorConcentration", acceptor),
                                ("DopingConcentration", donor - acceptor)):
            close(bundle["fields"][field]["values"][node["id"]], expected, "analytic doping mismatch")


    # Independent P95 component sizes on unique edges of each volume region.
    quality = contract["mesh_quality"]
    q = json.loads((output / "mesh_qualification.json").read_text())
    require(q["all_triangles_ccw"] and q["connected_components"] == quality["connected_components"],
            "mesh orientation/connectivity mismatch")
    region_edges = {name: set() for name in contract["regions"]}
    for tri in mesh["triangles"]:
        ids = tri["node_ids"]
        for i in range(3):
            region_edges[region_names[tri["region_id"]]].add(tuple(sorted((ids[i], ids[(i+1)%3]))))
    measured = {}
    for name, edges in region_edges.items():
        measured[name] = []
        for axis_index, axis in enumerate(("x", "y")):
            lengths = sorted(abs(nodes[a][axis] - nodes[b][axis]) for a, b in edges)
            p95 = lengths[max(0, math.ceil(0.95 * len(lengths)) - 1)]
            require(p95 <= quality["spacing_p95_factor"] * quality["region_max_spacing_um"][name][axis_index] + 1e-12,
                    f"mesh P95 spacing mismatch: {name}/{axis}")
            measured[name].append(p95)
    return {"node_count": len(nodes), "cell_count": len(mesh["triangles"]),
            "max_angle_degrees": q["max_angle_degrees"], "region_p95_spacing_um": measured,
            **({'mos':mos_report} if mos_report else {}),
            **({"bjt": bjt_report, "maxtransdiff": bundle["provenance"]["maxtransdiff"],
                "refinement_placements": q["refinement"]["placements"]} if bjt_report else {})}


def run_case(case: dict, output: Path, runtime: Path, repeat: int) -> dict:
    source = (REPO / case["source"]).resolve()
    row = {"id": case["id"], "source_sha256": digest(source), "expected": case["expected"],
           "stages": {name: "not_run" for name in ("parse", "model", "export", "native_import", "reference")}}
    first = None
    runs = []
    try:
        for _ in range(repeat if case["expected"] == "generated" else 1):
            proc = subprocess.run([sys.executable, str(REPO / "scripts/sentaurus_gmsh_pipeline.py"),
                    "--sde", str(source), "--output-dir", str(output), "--runtime-dir", str(runtime)]+[arg for key,value in case.get('parameters',{}).items() for arg in ('--parameter',key+'='+value)],
                    capture_output=True, text=True)
            if proc.returncode:
                require("output_safety" not in proc.stderr, proc.stderr.strip())
                failure = json.loads((output / "fail_report.json").read_text())
                row["failure"] = {k: failure[k] for k in ("stage", "message")}
                require(case["expected"] == "rejected", failure["message"])
                expected = case["rejection"]
                require(proc.returncode == 2 and failure["stage"] == expected["stage"] and
                        f":{expected['line']}:{expected['column']}:" in failure["message"] and
                        source.name in failure["message"] and expected["message_contains"] in failure["message"],
                        "rejection boundary changed")
                require(not (output / "run_manifest.json").exists() and
                        not any((output / name).exists() for name in PRIMARY), "stale success outputs")
                row["stages"]["parse"] = "expected_rejection"
                row["status"] = "expected_rejection"
                return row
            require(case["expected"] == "generated", "unsupported case unexpectedly generated; review contract")
            row["stages"].update(parse="passed", model="passed")
            manifest = json.loads((output / "run_manifest.json").read_text())
            require(manifest["status"] == "generated" and manifest["cross_format"]["status"] == "passed",
                    "cross-format check failed")
            quality = check_contract(output, case["contract"])
            row["stages"]["export"] = "passed"
            with tempfile.TemporaryDirectory(prefix="native-import-", dir=output.parent) as tmp:
                import_grid(output / "device.grd", Path(tmp))
                compare_native(output, Path(tmp))
            row["stages"]["native_import"] = "passed"
            hashes = {name: digest(output / name) for name in PRIMARY}
            runs.append({"sha256": hashes, "mesh_quality": quality,
                         **({'half_sha256':{name:digest(output/name) for name in
                            ('half_device.grd','half_device.dat','half_device.msh','half_mesh_field_bundle.json')}}
                            if case['contract']['doping']['kind']=='mos' else {})})
            if first is not None and case["repeat_policy"] == "byte_identical":
                require(first == hashes, "repeat output digest mismatch")
            if first is None:
                first = hashes
        row.update(status="generated", runs=runs, repeats=repeat,
                   repeat_policy=case["repeat_policy"],
                   byte_identical_observed=all(run["sha256"] == first for run in runs))
        return row
    except Exception as error:
        row.update(status="failed", error=str(error))
        return row
    finally:
        require(digest(source) == row["source_sha256"], "source modified by run")


def run_suite(manifest: Path, output: Path, runtime: Path, repeat: int = 2) -> dict:
    summary_path = output / "case_summary.json"
    require(summary_path.resolve() != manifest.resolve(), "manifest collides with summary")
    require(not summary_path.is_symlink(), "refusing linked summary")
    if summary_path.exists():
        require(json.loads(summary_path.read_text()).get("schema") == "vela.reference_sde_run.v1",
                "refusing to overwrite unrelated summary")
        summary_path.unlink()
    require(repeat >= 1, "repeat must be positive")
    data = load_manifest(manifest)
    for case in data["cases"]:
        require(not (REPO / case["source"]).resolve().is_relative_to(output.resolve()),
                "output directory contains an input source")
    output.mkdir(parents=True, exist_ok=True)
    rows = [run_case(case, output / case["id"], runtime, repeat) for case in data["cases"]]
    counts = Counter(row["status"] for row in rows)
    summary = {"schema": "vela.reference_sde_run.v1", "manifest_sha256": digest(manifest),
               "status": "passed" if not counts["failed"] else "failed", "counts": dict(counts),
               "reference": "not_run", "cases": rows}
    from sentaurus_gmsh_pipeline import _strict_json
    _strict_json(summary_path, summary)
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--runtime-dir", type=Path, default=REPO / "build/gmsh-python-runtime")
    parser.add_argument("--repeat", type=int, default=2)
    args = parser.parse_args(argv)
    try:
        summary = run_suite(args.manifest, args.output_dir, args.runtime_dir, args.repeat)
    except Exception as error:
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=2))
    return 0 if summary["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
