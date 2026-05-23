#!/usr/bin/env python3
"""Convert neutral reference TCAD CSV exports into Vela unit_scaling decks."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any


DEFAULT_SWEEPS: dict[str, dict[str, Any]] = {
    "iv": {
        "mode": "iv",
        "start": 0.0,
        "stop": 0.5,
        "step": 0.25,
        "write_vtk": False,
    },
    "cv": {
        "mode": "cv_quasistatic",
        "start": 0.0,
        "stop": 0.5,
        "step": 0.25,
        "write_vtk": False,
        "terminal_charge": {"per_meter": True},
    },
    "bv": {
        "mode": "bv_reverse",
        "start": 0.0,
        "stop": -1.0,
        "step": -0.5,
        "write_vtk": False,
        "breakdown": {
            "max_electric_field_V_per_m": 1.0e12,
            "current_jump_ratio": 1.0e12,
            "non_convergence": True,
        },
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_node_ids(value: str) -> list[int]:
    for sep in (";", "|"):
        value = value.replace(sep, ",")
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def load_mesh(input_dir: Path) -> tuple[dict[str, Any], dict[str, set[int]]]:
    nodes_rows = read_csv(input_dir / "nodes.csv")
    element_rows = read_csv(input_dir / "elements.csv")
    contact_rows = read_csv(input_dir / "contacts.csv")

    nodes = [
        {
            "id": int(row["id"]),
            "x": float(row["x_um"]),
            "y": float(row["y_um"]),
        }
        for row in nodes_rows
    ]

    region_cells: dict[str, list[int]] = defaultdict(list)
    region_materials: dict[str, str] = {}
    region_nodes: dict[str, set[int]] = defaultdict(set)
    triangles: list[dict[str, Any]] = []
    for row in element_rows:
        cell_id = int(row["id"])
        node_ids = [int(row["node0"]), int(row["node1"]), int(row["node2"])]
        region = row["region"]
        material = row["material"]
        triangles.append({
            "id": cell_id,
            "region_id": None,
            "node_ids": node_ids,
        })
        region_cells[region].append(cell_id)
        region_materials.setdefault(region, material)
        region_nodes[region].update(node_ids)

    region_ids = {name: idx for idx, name in enumerate(region_cells)}
    for tri, row in zip(triangles, element_rows, strict=True):
        tri["region_id"] = region_ids[row["region"]]

    regions = [
        {
            "id": region_ids[name],
            "name": name,
            "material": region_materials[name],
            "cell_ids": cells,
        }
        for name, cells in region_cells.items()
    ]

    contacts = [
        {
            "id": idx,
            "name": row["name"],
            "region_id": region_ids[row["region"]],
            "node_ids": parse_node_ids(row["node_ids"]),
        }
        for idx, row in enumerate(contact_rows)
    ]

    mesh = {
        "nodes": nodes,
        "triangles": triangles,
        "regions": regions,
        "contacts": contacts,
    }
    return mesh, region_nodes


def load_region_doping(input_dir: Path, region_nodes: dict[str, set[int]]) -> list[dict[str, Any]]:
    rows = read_csv(input_dir / "doping.csv")
    by_node = {
        int(row["node_id"]): (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
        for row in rows
    }

    doping: list[dict[str, Any]] = []
    for region, node_ids in region_nodes.items():
        samples: list[tuple[float, float]] = []
        for node_id in sorted(node_ids):
            donor, acceptor = by_node.get(node_id, (0.0, 0.0))
            samples.append((donor, acceptor))
        dominant = dominant_net_sign(samples)
        if dominant != 0:
            selected = [
                sample for sample in samples
                if sign(sample[0] - sample[1]) == dominant
            ]
        else:
            selected = samples
        donors = [sample[0] for sample in selected]
        acceptors = [sample[1] for sample in selected]
        count = max(len(donors), 1)
        doping.append({
            "region": region,
            "donors": sum(donors) / count,
            "acceptors": sum(acceptors) / count,
        })
    return doping


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def dominant_net_sign(samples: list[tuple[float, float]]) -> int:
    counts = {1: 0, -1: 0}
    for donor, acceptor in samples:
        net_sign = sign(donor - acceptor)
        if net_sign != 0:
            counts[net_sign] += 1
    if counts[1] > counts[-1]:
        return 1
    if counts[-1] > counts[1]:
        return -1
    return 0


def choose_sweep_contact(contacts: list[dict[str, Any]], device: str) -> str:
    preferred = {
        "pn_diode": ["anode", "cathode"],
        "nmos2d": ["drain", "gate"],
        "pmos2d": ["drain", "gate"],
        "ldmos2d": ["drain", "gate"],
        "igbt2d": ["collector", "anode", "drain"],
    }.get(device, [])
    names = [contact["name"] for contact in contacts]
    for name in preferred:
        if name in names:
            return name
    if not names:
        raise ValueError("contacts.csv must contain at least one contact")
    return names[0]


def base_deck(mesh_file: str,
              output_csv: str,
              mesh: dict[str, Any],
              doping: list[dict[str, Any]],
              device: str,
              simulation_type: str) -> dict[str, Any]:
    sweep_contact = choose_sweep_contact(mesh["contacts"], device)
    contacts = [{"name": contact["name"], "bias": 0.0} for contact in mesh["contacts"]]
    sweep = dict(DEFAULT_SWEEPS[simulation_type])
    sweep["contact"] = sweep_contact
    sweep["current_contact"] = sweep_contact
    if simulation_type == "cv":
        sweep["terminal_charge"] = dict(sweep["terminal_charge"])
        sweep["terminal_charge"]["contact"] = sweep_contact

    return {
        "_comment": (
            "Generated from explicit reference TCAD CSV export. "
            "Coordinates, doping, mobility-like inputs, fields, and sheet densities "
            "use scaling.mode = unit_scaling conventions."
        ),
        "simulation_type": "dc_sweep",
        "mesh_file": mesh_file,
        "output_csv": output_csv,
        "scaling": {"mode": "unit_scaling"},
        "node_doping_file": "doping.csv",
        "doping": doping,
        "contacts": contacts,
        "solver": {
            "method": "gummel",
            "max_iter": 100,
            "reltol": 1.0e-5,
            "damping_psi": 0.5,
        },
        "sweep": sweep,
    }


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--device",
        choices=["pn_diode", "nmos2d", "pmos2d", "ldmos2d", "igbt2d"],
        required=True,
    )
    parser.add_argument(
        "--simulation-types",
        default="iv",
        help="Comma-separated subset of iv,cv,bv.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    simulation_types = [
        value.strip() for value in args.simulation_types.split(",") if value.strip()
    ]
    unknown = sorted(set(simulation_types) - set(DEFAULT_SWEEPS))
    if unknown:
        raise ValueError(f"unsupported simulation type(s): {', '.join(unknown)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    mesh, region_nodes = load_mesh(args.input_dir)
    doping = load_region_doping(args.input_dir, region_nodes)
    write_json(args.output_dir / "mesh.json", mesh)
    shutil.copyfile(args.input_dir / "doping.csv", args.output_dir / "doping.csv")

    generated = ["mesh.json", "doping.csv"]
    for sim_type in simulation_types:
        suffix = {"iv": "iv", "cv": "cv", "bv": "bv"}[sim_type]
        deck_name = f"simulation_{suffix}.json"
        csv_name = f"{args.device}_{suffix}.csv"
        write_json(
            args.output_dir / deck_name,
            base_deck("mesh.json", csv_name, mesh, doping, args.device, sim_type),
        )
        generated.append(deck_name)

    print(json.dumps({"generated": generated}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
