"""Stage the exact archived TransportModels mesh for isothermal DD sweeps."""

from __future__ import annotations

from collections import Counter
import csv
import json
import math
import shutil

from sentaurus_mesh_exporters import write_dfise


def stage_transportmodels_inputs(fixture, work, stem):
    work.mkdir()
    for name in (f"{stem}.cmd", f"{stem}.vela.json"):
        shutil.copyfile(fixture / name, work / name)
    source_par = (fixture.parent / "source" / f"{stem}.par").read_text(encoding="utf-8")
    poly_par = (fixture / "polysilicon.par").read_text(encoding="utf-8")
    (work / f"{stem}.par").write_text(source_par + "\n" + poly_par, encoding="utf-8")

    archived = fixture.parent / "vela"
    mesh = json.loads((archived / "mesh.json").read_text(encoding="utf-8"))
    material_names = sorted({region["material"] for region in mesh["regions"]})
    materials = [
        {"id": index, "name": "Silicon" if name == "Si" else name}
        for index, name in enumerate(material_names)
    ]
    regions = {r["id"]: material_names.index(r["material"]) for r in mesh["regions"]}
    triangles = [
        {**cell, "dfise_id": index, "material_id": regions[cell["region_id"]]}
        for index, cell in enumerate(mesh["triangles"])
    ]
    incidence = Counter()
    for triangle in triangles:
        a, b, c = triangle["node_ids"]
        incidence.update(tuple(sorted(edge)) for edge in ((a, b), (b, c), (c, a)))
    contacts = []
    lines = []
    for contact in mesh["contacts"]:
        nodes = set(contact["node_ids"])
        edges = sorted(edge for edge, count in incidence.items()
                       if count == 1 and nodes.issuperset(edge))
        if {node for edge in edges for node in edge} != nodes:
            raise ValueError(f"Incomplete archived boundary for {contact['name']}")
        tag = contact["id"]
        contacts.append({**contact, "physical_tag": tag})
        lines.extend({"physical_curve_tags": [tag], "node_ids": edge} for edge in edges)

    with (archived / "doping.csv").open(encoding="utf-8", newline="") as stream:
        doping = sorted(csv.DictReader(stream), key=lambda row: int(row["node_id"]))
    vertex_ids = list(range(len(mesh["nodes"])))
    if ([node["id"] for node in mesh["nodes"]] != vertex_ids
            or [int(row["node_id"]) for row in doping] != vertex_ids):
        raise ValueError("Archived geometry and doping must use the same global vertex numbering")
    donors = [float(row["donors_cm3"]) for row in doping]
    acceptors = [float(row["acceptors_cm3"]) for row in doping]
    net = [donor - acceptor for donor, acceptor in zip(donors, acceptors)]
    magnitudes = [abs(value) for value in net]
    fields = {
        "DonorConcentration": {"values": donors},
        "AcceptorConcentration": {"values": acceptors},
        "DopingConcentration": {"values": net},
        "AbsDopingConcentration": {"values": magnitudes},
        "Log10AbsDopingConcentration": {
            "values": [math.log10(value) if value else 0.0 for value in magnitudes]
        },
        # The source mesh has one donor species (arsenic) and one acceptor
        # species (boron). Sentaurus reconstructs their sum from these fields.
        "species": {
            "ArsenicActiveConcentration": {"values": donors},
            "BoronActiveConcentration": {"values": acceptors},
        },
    }
    write_dfise(work / "n1_msh.grd", work / "n1_msh.dat", {
        **mesh, "materials": materials, "triangles": triangles, "contacts": contacts,
        "interfaces": [], "lines": lines, "fields": fields,
    })
