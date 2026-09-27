"""Mesh import for ``vela device`` (WP2, contract C2).

Converts the ``File{Grid}`` input into the Vela mesh bundle written next to the
generated decks:

* ``mesh.json``   -- the format read by the C++ ``JsonMeshReader`` (same keys as
  ``scripts/convert_tcad_export.py``: ``nodes[{id,x,y}]``, ``triangles[{id,
  region_id,node_ids}]``, ``regions[{id,name,material,cell_ids}]``,
  ``contacts[{id,name,region_id,node_ids,edge_node_ids}]``; coordinates in um);
* ``doping.csv``  -- ``node_id,donors_cm3,acceptors_cm3``;
* ``mesh_qualification.json`` -- contract C2 report
  (``schemas/vela.mesh_qualification.v1.schema.json``).

Two input paths:

* DF-ISE ``.grd`` + sibling ``.dat`` (pure Python, reuses ``parse_grd`` /
  ``parse_dat`` from ``scripts/compare_sentaurus_tdr_tdx.py``).  Elements are
  edge-indexed; triangles are rebuilt as counter-clockwise node triples,
  ``material=Contact`` regions become contacts.  Quadrilaterals and 3-D
  elements are rejected (E402) -- no deterministic split is attempted in
  stage 1.  Vertex datasets must be global (one value per vertex); region-local
  dataset ordering is rejected with E402.
* Sentaurus ``.tdr`` through the external ``sentaurus_import`` executable
  (``--tdr <file> --export-dir <dir>``), whose CSV export is converted with the
  same rules as ``convert_tcad_export.load_mesh``.  Without an importer the
  import fails closed (E402).

Material names follow ``src/io/SentaurusTdrReader.cpp``: ``Silicon -> Si``,
``Oxide -> SiO2``; anything else is passed through verbatim and recorded in
``regions[].source_material``.
"""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_SCRIPTS = Path(__file__).resolve().parents[2]
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from compare_sentaurus_tdr_tdx import parse_dat, parse_grd  # noqa: E402

SCHEMA_ID = "vela.mesh_qualification.v1"
MATERIAL_MAP = {"silicon": "Si", "si": "Si", "oxide": "SiO2", "sio2": "SiO2"}
DOPING_COLUMNS = ["node_id", "donors_cm3", "acceptors_cm3"]
PLAUSIBLE_EXTENT_UM = (1e-3, 1e4)
DONOR_DATASETS = ("DonorConcentration", "PhosphorusActiveConcentration", "ArsenicActiveConcentration",
                  "AntimonyActiveConcentration")
ACCEPTOR_DATASETS = ("AcceptorConcentration", "BoronActiveConcentration", "IndiumActiveConcentration")

__all__ = ["MeshImportError", "import_grid", "check_electrodes", "find_importer", "SCHEMA_ID"]


class MeshImportError(Exception):
    """Fail-closed import problem carrying the catalogue message id (E114/E304/E402)."""

    def __init__(self, message_id: str, message: str) -> None:
        super().__init__(f"{message_id}: {message}")
        self.message_id = message_id
        self.message = message


def map_material(name: str) -> str:
    return MATERIAL_MAP.get(name.strip().lower(), name.strip())


def find_importer(importer: Optional[str] = None) -> Optional[List[str]]:
    """Resolve the sentaurus_import executable (argument, env, or build directories)."""
    candidates: List[Optional[str]] = [importer, os.environ.get("VELA_SENTAURUS_IMPORT")]
    exe = "sentaurus_import.exe" if sys.platform.startswith("win") else "sentaurus_import"
    worktree = Path(__file__).resolve().parents[3]
    for root in (worktree / "build", worktree / "build-release", Path(r"D:\code-repo\vela-tcad\build")):
        candidates.append(str(root / exe))
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return [candidate]
    return None


# --- geometry helpers ------------------------------------------------------------------

def _signed_area(p0, p1, p2) -> float:
    return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))


def _angles_deg(p0, p1, p2) -> List[float]:
    def length(a, b):
        return math.hypot(a[0] - b[0], a[1] - b[1])
    a, b, c = length(p1, p2), length(p0, p2), length(p0, p1)
    def angle(opposite, s1, s2):
        cos_value = (s1 * s1 + s2 * s2 - opposite * opposite) / (2.0 * s1 * s2)
        return math.degrees(math.acos(max(-1.0, min(1.0, cos_value))))
    return [angle(a, b, c), angle(b, a, c), angle(c, a, b)]


def _ccw(node_ids: Sequence[int], coords: Dict[int, Tuple[float, float]]) -> List[int]:
    ids = list(node_ids)
    if _signed_area(coords[ids[0]], coords[ids[1]], coords[ids[2]]) < 0.0:
        ids[1], ids[2] = ids[2], ids[1]
    return ids


# --- DF-ISE path -------------------------------------------------------------------------

def _edge_nodes(ref: int, edges: List[List[int]]) -> Tuple[int, int]:
    if ref >= 0:
        start, end = edges[ref]
    else:
        end, start = edges[-ref - 1]
    return int(start), int(end)


def _triangle_from_edges(record: List[int], edges: List[List[int]]) -> List[int]:
    chain = [_edge_nodes(ref, edges) for ref in record[1:4]]
    nodes = [chain[0][0], chain[0][1]]
    for start, end in chain[1:]:
        if start == nodes[-1] and end not in nodes:
            nodes.append(end)
        elif end == nodes[-1] and start not in nodes:
            nodes.append(start)
    if len(nodes) != 3:
        nodes = sorted({node for pair in chain for node in pair})
        if len(nodes) != 3:
            raise MeshImportError("E402", f"DF-ISE element {record} does not describe a triangle")
    return nodes


def _load_dfise(grd_path: Path) -> dict:
    dat_path = grd_path.with_suffix(".dat")
    if not dat_path.is_file():
        raise MeshImportError("E304", f"DF-ISE companion file '{dat_path.name}' not found")
    grid = parse_grd(grd_path)
    datasets = parse_dat(dat_path)
    vertices = grid["vertices"]
    if not vertices:
        raise MeshImportError("E402", f"no vertices found in '{grd_path.name}'")
    dimension = len(vertices[0])
    if dimension != 2:
        raise MeshImportError("E402", f"grid dimension {dimension} is not 2-D")
    coords = {index: (float(v[0]), float(v[1])) for index, v in enumerate(vertices)}
    edges = grid["edges"]
    elements = grid["elements"]
    by_kind = {"point": 0, "line": 0, "triangle": 0, "quadrilateral": 0, "tetrahedron": 0, "other": 0}
    for record in elements:
        kind = record[0] if record else -1
        key = {0: "point", 1: "line", 2: "triangle", 3: "quadrilateral", 5: "tetrahedron"}.get(kind, "other")
        by_kind[key] += 1
    if by_kind["quadrilateral"] or by_kind["tetrahedron"] or by_kind["other"]:
        raise MeshImportError("E402", "grid contains non-triangle 2-D elements "
                              f"(quadrilaterals={by_kind['quadrilateral']}, tetrahedra={by_kind['tetrahedron']}, "
                              f"other={by_kind['other']}); only triangles are supported")

    triangles: List[dict] = []
    regions: List[dict] = []
    contacts: List[dict] = []
    element_region: Dict[int, int] = {}
    region_names = [name for name in grid["region_blocks"]]
    for name in region_names:
        info = grid["region_blocks"][name]
        element_ids = grid["region_elements"].get(name, [])
        material = str(info.get("material") or "")
        if material.lower() == "contact":
            node_ids: List[int] = []
            edge_pairs: List[List[int]] = []
            for element_id in element_ids:
                record = elements[element_id]
                if record[0] != 1:
                    raise MeshImportError("E402", f"contact region '{name}' contains a non-line element")
                if len(record) == 3:
                    pair = [int(record[1]), int(record[2])]
                else:
                    pair = list(_edge_nodes(int(record[1]), edges))
                edge_pairs.append(pair)
                for node in pair:
                    if node not in node_ids:
                        node_ids.append(node)
            contacts.append({"name": name, "node_ids": sorted(node_ids), "edge_node_ids": edge_pairs,
                             "adjacent_regions": set()})
            continue
        region_id = len(regions)
        cell_ids: List[int] = []
        for element_id in element_ids:
            record = elements[element_id]
            if record[0] != 2:
                continue
            node_ids = _ccw(_triangle_from_edges(record, edges), coords)
            cell_id = len(triangles)
            triangles.append({"id": cell_id, "region_id": region_id, "node_ids": node_ids})
            cell_ids.append(cell_id)
            element_region[element_id] = region_id
        regions.append({"id": region_id, "name": name, "material": map_material(material),
                        "source_material": material, "cell_ids": cell_ids})

    # contact adjacency: the region owning a triangle that contains both edge nodes
    node_sets = [(set(tri["node_ids"]), tri["region_id"]) for tri in triangles]
    for contact in contacts:
        for pair in contact["edge_node_ids"]:
            pair_set = set(pair)
            for nodes, region_id in node_sets:
                if pair_set <= nodes:
                    contact["adjacent_regions"].add(region_id)
                    break

    vertex_count = len(coords)
    doping, source_fields, net_only = _doping_from_datasets(datasets, vertex_count)
    return {
        "format": "dfise", "coords": coords, "triangles": triangles, "regions": regions,
        "contacts": contacts, "by_kind": by_kind, "doping": doping, "source_fields": source_fields,
        "net_only": net_only, "unit_declared": "um", "unit_scale": 1.0,
        "companion": dat_path,
    }


def _doping_from_datasets(datasets: dict, vertex_count: int):
    def values_of(name: str) -> Optional[List[float]]:
        entry = datasets.get(name)
        if entry is None:
            return None
        values = entry["values"]
        if len(values) != vertex_count:
            raise MeshImportError("E402", f"dataset '{name}' has {len(values)} values for {vertex_count} vertices; "
                                  "region-local dataset ordering is not supported")
        return [float(v) for v in values]

    donors = [0.0] * vertex_count
    acceptors = [0.0] * vertex_count
    used: List[str] = []
    for name in DONOR_DATASETS:
        values = values_of(name)
        if values is not None and name != "DonorConcentration" and "DonorConcentration" in datasets:
            continue
        if values is not None:
            donors = [a + b for a, b in zip(donors, values)]
            used.append(name)
    for name in ACCEPTOR_DATASETS:
        values = values_of(name)
        if values is not None and name != "AcceptorConcentration" and "AcceptorConcentration" in datasets:
            continue
        if values is not None:
            acceptors = [a + b for a, b in zip(acceptors, values)]
            used.append(name)
    net_only = False
    if not used:
        net = values_of("DopingConcentration")
        if net is None:
            raise MeshImportError("E402", "no doping dataset (DonorConcentration/AcceptorConcentration/"
                                  "DopingConcentration) found in the DF-ISE data file")
        donors = [max(v, 0.0) for v in net]
        acceptors = [max(-v, 0.0) for v in net]
        used = ["DopingConcentration"]
        net_only = True
    rows = [{"node_id": index, "donors_cm3": donors[index], "acceptors_cm3": acceptors[index]}
            for index in range(vertex_count)]
    return rows, used, net_only


# --- TDR path ------------------------------------------------------------------------------

def _load_tdr(tdr_path: Path, importer: Optional[str], work_dir: Path) -> dict:
    command = find_importer(importer)
    if command is None:
        raise MeshImportError("E402", "no sentaurus_import executable available for .tdr grids "
                              "(set VELA_SENTAURUS_IMPORT or build sentaurus_import)")
    export_dir = work_dir / "tdr_export"
    export_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([*command, "--tdr", str(tdr_path), "--export-dir", str(export_dir),
                           "--compensated-doping-policy", "reported"], capture_output=True, text=True)
    if proc.returncode != 0:
        raise MeshImportError("E402", f"sentaurus_import failed ({proc.returncode}): {proc.stderr.strip()[:400]}")
    return _load_export_dir(export_dir)


def _load_export_dir(export_dir: Path) -> dict:
    """Convert the sentaurus_import CSV export (nodes/elements/contacts/doping.csv)."""
    def read(name: str) -> List[dict]:
        with (export_dir / name).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    coords = {int(r["id"]): (float(r["x_um"]), float(r["y_um"])) for r in read("nodes.csv")}
    triangles: List[dict] = []
    regions: List[dict] = []
    region_index: Dict[str, int] = {}
    for row in read("elements.csv"):
        name = row["region"]
        if name not in region_index:
            region_index[name] = len(regions)
            regions.append({"id": len(regions), "name": name, "material": map_material(row["material"]),
                            "source_material": row["material"], "cell_ids": []})
        node_ids = _ccw([int(row["node0"]), int(row["node1"]), int(row["node2"])], coords)
        cell_id = len(triangles)
        triangles.append({"id": cell_id, "region_id": region_index[name], "node_ids": node_ids})
        regions[region_index[name]]["cell_ids"].append(cell_id)
    contacts: List[dict] = []
    contacts_path = export_dir / "contacts.csv"
    if contacts_path.is_file():
        for row in read("contacts.csv"):
            raw = str(row.get("node_ids", "")).replace("|", ";")
            node_ids = sorted(int(v) for v in raw.split(";") if v.strip())
            pairs = []
            for item in str(row.get("edge_node_ids", "")).replace("|", ";").split(";"):
                if item.strip():
                    pairs.append([int(part) for part in item.split("-")])
            adjacent = {region_index[row["region"]]} if row.get("region") in region_index else set()
            contacts.append({"name": row["name"], "node_ids": node_ids, "edge_node_ids": pairs,
                             "adjacent_regions": adjacent})
    doping_rows = []
    doping_path = export_dir / "doping.csv"
    source_fields = []
    if doping_path.is_file():
        for row in read("doping.csv"):
            doping_rows.append({"node_id": int(row["node_id"]), "donors_cm3": float(row["donors_cm3"]),
                                "acceptors_cm3": float(row["acceptors_cm3"])})
        source_fields = ["DonorConcentration", "AcceptorConcentration"]
    else:
        doping_rows = [{"node_id": i, "donors_cm3": 0.0, "acceptors_cm3": 0.0} for i in sorted(coords)]
    return {
        "format": "tdr", "coords": coords, "triangles": triangles, "regions": regions, "contacts": contacts,
        "by_kind": {"triangle": len(triangles)}, "doping": sorted(doping_rows, key=lambda r: r["node_id"]),
        "source_fields": source_fields, "net_only": False, "unit_declared": "um", "unit_scale": 1.0,
        "companion": None,
    }


# --- outputs ---------------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_outputs(model: dict, grid_path: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    coords = model["coords"]
    triangles = model["triangles"]
    regions = model["regions"]
    contacts_out = []
    for index, contact in enumerate(model["contacts"]):
        adjacent = sorted(contact["adjacent_regions"])
        contacts_out.append({"id": index, "name": contact["name"], "region_id": adjacent[0] if adjacent else 0,
                             "node_ids": contact["node_ids"], "edge_node_ids": contact["edge_node_ids"]})
    mesh = {
        "nodes": [{"id": node_id, "x": xy[0], "y": xy[1]} for node_id, xy in sorted(coords.items())],
        "triangles": triangles,
        "regions": [{"id": r["id"], "name": r["name"], "material": r["material"], "cell_ids": r["cell_ids"]}
                    for r in regions],
        "contacts": contacts_out,
    }
    (out_dir / "mesh.json").write_text(json.dumps(mesh, indent=1) + "\n", encoding="utf-8")
    with (out_dir / "doping.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DOPING_COLUMNS)
        writer.writeheader()
        for row in model["doping"]:
            writer.writerow({"node_id": row["node_id"], "donors_cm3": repr(float(row["donors_cm3"])),
                             "acceptors_cm3": repr(float(row["acceptors_cm3"]))})

    xs = [xy[0] for xy in coords.values()]
    ys = [xy[1] for xy in coords.values()]
    max_deg, min_deg, obtuse = 0.0, 180.0, 0
    for tri in triangles:
        angles = _angles_deg(*(coords[n] for n in tri["node_ids"]))
        max_deg = max(max_deg, max(angles))
        min_deg = min(min_deg, min(angles))
        if max(angles) > 90.0 + 1e-9:
            obtuse += 1
    extent = (max(xs) - min(xs), max(ys) - min(ys))
    plausible = all(PLAUSIBLE_EXTENT_UM[0] <= e <= PLAUSIBLE_EXTENT_UM[1] for e in extent)
    materials: List[str] = []
    for region in regions:
        if region["material"] not in materials:
            materials.append(region["material"])
    by_kind = {"triangle": int(model["by_kind"].get("triangle", len(triangles)))}
    for key in ("point", "line", "quadrilateral", "tetrahedron", "other"):
        if model["by_kind"].get(key):
            by_kind[key] = int(model["by_kind"][key])
    report = {
        "schema": SCHEMA_ID,
        "source": {"grid_file": grid_path.name, "format": model["format"], "sha256": _sha256(grid_path),
                   "importer": "vela_sdevice.mesh_import dfise_reader v1" if model["format"] == "dfise"
                   else "sentaurus_import --export-dir"},
        "dimension": 2,
        "coordinate_unit": {"declared": model["unit_declared"], "applied_scale_to_um": model["unit_scale"]},
        "bounding_box_um": {"min": [min(xs), min(ys)], "max": [max(xs), max(ys)]},
        "nodes": {"count": len(coords)},
        "elements": {"count": len(triangles), "by_kind": by_kind, "converted_quadrilaterals": 0},
        "regions": [{"name": r["name"], "material": r["material"], "element_count": len(r["cell_ids"]),
                     **({"source_material": r["source_material"]} if r["source_material"] != r["material"] else {})}
                    for r in regions],
        "materials": materials,
        "contacts": [{"name": c["name"], "node_count": len(c["node_ids"]), "edge_count": len(c["edge_node_ids"]),
                      "adjacent_materials": sorted({regions[r]["material"] for r in c["adjacent_regions"]})}
                     for c in model["contacts"]],
        "angles": {"max_deg": max_deg, "obtuse_element_count": obtuse,
                   "obtuse_fraction": (obtuse / len(triangles)) if triangles else 0.0, "min_deg": min_deg},
        "doping": {"file": "doping.csv", "columns": DOPING_COLUMNS, "source_fields": model["source_fields"],
                   "net_only": bool(model["net_only"])},
        "outputs": {"mesh_json": "mesh.json", "doping_csv": "doping.csv"},
        "qualification": {
            "two_dimensional": True,
            "triangles_only": not any(model["by_kind"].get(k) for k in ("quadrilateral", "tetrahedron", "other")),
            "coordinates_plausible": plausible,
            "non_obtuse": obtuse == 0,
            "recommended_node_volume_policy": "barycentric",
            "allow_obtuse": True,
        },
    }
    if model.get("companion") is not None:
        report["source"]["companion_file"] = Path(model["companion"]).name
        report["source"]["companion_sha256"] = _sha256(Path(model["companion"]))
    (out_dir / "mesh_qualification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def import_grid(grid_path: os.PathLike, out_dir: os.PathLike, *, importer: Optional[str] = None) -> dict:
    """Convert a Sentaurus grid into mesh.json / doping.csv / mesh_qualification.json in ``out_dir``."""
    grid_path = Path(grid_path)
    out_dir = Path(out_dir)
    suffix = grid_path.suffix.lower()
    if not grid_path.is_file():
        raise MeshImportError("E302", f"Grid file '{grid_path.name}' not found")
    if suffix == ".grd":
        model = _load_dfise(grid_path)
    elif suffix == ".tdr":
        model = _load_tdr(grid_path, importer, out_dir)
    else:
        raise MeshImportError("E114", f"Grid file '{grid_path.name}' has an unsupported extension (expected .tdr or .grd)")
    if not model["triangles"]:
        raise MeshImportError("E402", "grid contains no triangles")
    return _write_outputs(model, grid_path, out_dir)


def check_electrodes(report: dict, electrode_names: Sequence[str]) -> List[dict]:
    """Gate-5 helper: match cmd electrode names against mesh contacts (case-sensitive, like sdevice)."""
    contacts = {c["name"] for c in report.get("contacts", [])}
    result = [{"electrode": name, "matched": name in contacts} for name in electrode_names]
    report.setdefault("qualification", {})["electrodes_checked"] = result
    return result


def _main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="vela_sdevice.mesh_import",
                                     description="Convert a Sentaurus .grd/.dat or .tdr grid into the Vela mesh bundle.")
    parser.add_argument("grid", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--importer", help="sentaurus_import executable for .tdr input")
    args = parser.parse_args(argv)
    try:
        report = import_grid(args.grid, args.out_dir, importer=args.importer)
    except MeshImportError as error:
        print(f"{error.message_id} {error.message}", file=sys.stderr)
        return {"E114": 1, "E302": 3, "E304": 3, "E402": 4}.get(error.message_id, 4)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
