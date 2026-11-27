#!/usr/bin/env python3
"""Prepare SI point inputs for the production IALMob formula from an exported TDR.

This separates formula/parameter replay from Vela's cell/edge field reconstruction.
The default uses native eEnormal/hEnormal exports (case-sensitive metadata).
An electric-field projection is available as a separate reconstruction control.
Generated evidence belongs in an ignored build dir.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path


def extract_parameters(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    blocks = re.findall(r"^IALMob\s*:\s*\{(.*?)^\s*\}", text, re.M | re.S)
    if len(blocks) != 1:
        raise ValueError("Expected exactly one unnamed IALMob parameter block")
    result = {"electron": {}, "hole": {}}
    shared = {"nref_D", "nref_A", "cref_D", "cref_A", "me_over_m0", "mh_over_m0"}
    for line in blocks[0].splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"(\w+)\s*=\s*(.+)", line)
        if not match:
            raise ValueError(f"Unrecognized IALMob parameter line: {line}")
        key, payload = match.groups()
        values = [float(x.strip()) for x in payload.split(",")]
        if not all(math.isfinite(x) for x in values) or len(values) not in (1, 2):
            raise ValueError(f"Invalid IALMob parameter: {key}")
        if key in shared:
            if len(values) != 1 or key in result:
                raise ValueError(f"Invalid shared IALMob parameter: {key}")
            result[key] = values[0]
        else:
            if key in result["electron"]:
                raise ValueError(f"Duplicate IALMob parameter: {key}")
            result["electron"][key] = values[0]
            result["hole"][key] = values[-1]
    return result


def read_field(directory: Path, name: str, region: int = 0) -> dict[int, tuple[float, ...]]:
    result = {}
    with (directory / "fields" / f"{name}_region{region}.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            node = int(row["node_id"])
            values = tuple(float(row[f"component{i}"]) for i in range(len(row) - 1))
            if node in result or not all(math.isfinite(x) for x in values):
                raise ValueError(f"Duplicate or non-finite {name} node: {node}")
            result[node] = values
    return result


def nearest_segment(point: tuple[float, float], segments: list[tuple]) -> tuple[float, float, float]:
    nearest = math.inf
    normal = (0.0, 0.0)
    for ax, ay, bx, by in segments:
        dx, dy = bx - ax, by - ay
        squared = dx * dx + dy * dy
        if squared <= 0:
            raise ValueError("Degenerate interface segment")
        fraction = min(1.0, max(0.0, ((point[0] - ax) * dx + (point[1] - ay) * dy) / squared))
        distance = math.hypot(point[0] - ax - fraction * dx, point[1] - ay - fraction * dy)
        if distance < nearest:
            nearest = distance
            length = math.sqrt(squared)
            # Gradient of true geometric distance, including segment ends.
            normal = ((point[0] - ax - fraction * dx) / distance,
                      (point[1] - ay - fraction * dy) / distance) if distance > 1e-14 else (-dy / length, dx / length)
    if not math.isfinite(nearest):
        raise ValueError("No qualifying semiconductor/insulator interface")
    return nearest, normal[0], normal[1]


def prepare(export_dir: Path, parameter_file: Path, output_dir: Path,
            normal_source: str = "native_enormal") -> dict:
    if normal_source not in {"native_enormal", "electric_projection"}:
        raise ValueError("Unsupported normal field source")
    export_dir, parameter_file, output_dir = export_dir.resolve(), parameter_file.resolve(), output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((export_dir / "metadata.json").read_text(encoding="utf-8"))
    silicon = [r for r in metadata.get("regions", []) if r.get("name") == "Silicon_1"]
    if len(silicon) != 1 or silicon[0].get("index") != 0 or silicon[0].get("material") not in {"Si", "Silicon"}:
        raise ValueError("This Slot LDMOS replay requires Silicon_1 at region index 0")
    units = {(f["name"], f["region"]): f["unit"] for f in metadata["fields"]}
    required_units = {"eDensity": "cm^-3", "hDensity": "cm^-3", "DonorConcentration": "cm^-3",
        "AcceptorConcentration": "cm^-3", "ElectricField": "V*cm^-1", "eGradQuasiFermi": "V*cm^-1",
        "hGradQuasiFermi": "V*cm^-1", "eMobility": "cm^2*V^-1*s^-1", "hMobility": "cm^2*V^-1*s^-1",
        "LatticeTemperature": "K"}
    if normal_source == "native_enormal":
        required_units.update({"eEnormal": "V*cm^-1", "hEnormal": "V*cm^-1"})
    for name, unit in required_units.items():
        if units.get((name, 0)) != unit:
            raise ValueError(f"Unexpected unit for {name}: {units.get((name, 0))}")
    fields = {name: read_field(export_dir, name) for name in required_units}
    nodes = set(fields["eDensity"])
    if not nodes or any(set(field) != nodes for field in fields.values()):
        raise ValueError("The silicon field supports must match exactly")
    if any(len(v) != 1 or abs(v[0] - 300) > 1e-8 for v in fields["LatticeTemperature"].values()):
        raise ValueError("The original Slot IALMob high-field profile is defined at uniform 300 K")
    with (export_dir / "nodes.csv").open(newline="", encoding="utf-8") as stream:
        coordinates = {int(row["id"]): (float(row["x_um"]), float(row["y_um"])) for row in csv.DictReader(stream)}
    edges = defaultdict(set)
    with (export_dir / "elements.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            triangle = [int(row[f"node{i}"]) for i in range(3)]
            for i in range(3):
                edge = tuple(sorted((triangle[i], triangle[(i + 1) % 3])))
                edges[edge].add((row["region"], row["material"]))
    segments = []
    interface_pairs = defaultdict(int)
    for (a, b), regions in edges.items():
        if not any(name == "Silicon_1" for name, _ in regions):
            continue
        insulators = {name for name, material in regions if material not in {"Si", "Silicon", "PolySi", "PolySilicon"}}
        if insulators:
            segments.append((*coordinates[a], *coordinates[b]))
            for name in insulators:
                interface_pairs[name] += 1
    points, reference = [], []
    for node in sorted(nodes):
        distance, nx, ny = nearest_segment(coordinates[node], segments)
        ex, ey = fields["ElectricField"][node]
        normal = abs(ex * nx + ey * ny) * 100.0
        points.append({"node_id": node, "donors_m3": fields["DonorConcentration"][node][0] * 1e6,
            "acceptors_m3": fields["AcceptorConcentration"][node][0] * 1e6,
            "electrons_m3": fields["eDensity"][node][0] * 1e6, "holes_m3": fields["hDensity"][node][0] * 1e6,
            "electron_driving_field_V_per_m": math.hypot(*fields["eGradQuasiFermi"][node]) * 100,
            "hole_driving_field_V_per_m": math.hypot(*fields["hGradQuasiFermi"][node]) * 100,
            "electron_normal_field_V_per_m": abs(fields["eEnormal"][node][0]) * 100 if normal_source == "native_enormal" else normal,
            "hole_normal_field_V_per_m": abs(fields["hEnormal"][node][0]) * 100 if normal_source == "native_enormal" else normal,
            "interface_distance_m": distance * 1e-6})
        reference.append({"node_id": node, "electron_mobility_cm2_V_s": fields["eMobility"][node][0],
            "hole_mobility_cm2_V_s": fields["hMobility"][node][0]})
    mobility = {"model": "ialmob_field", "ialmob": extract_parameters(parameter_file)}
    config = {"simulation_type": "ialmob_point_probe", "temperature_K": 300,
        "mobility": mobility, "points_file": "points.json", "output_csv": "mobility.csv"}
    manifest = {"schema": "slot_ldmos_ialmob_point_replay.v1", "nodes": len(points),
        "export_directory": str(export_dir), "temperature_K": 300, "parameter_file": str(parameter_file),
        "parameter_sha256": hashlib.sha256(parameter_file.read_bytes()).hexdigest(),
        "interface_segment_counts": dict(interface_pairs), "normal_field_source": normal_source,
        "drive_source": "exported GradQuasiFermi vector magnitude",
        "limitations": ["A nodal model replay; not a self-consistent I-V point or an edge-discretization parity test.",
            "Nonlinear model evaluation and element/vertex output averaging do not commute.",
            "Uses isotropic unnamed IALMob, full PhuMob, default PhononCombination=1, no stress or ThinLayer."]}
    for name, value in [("config.json", config), ("points.json", points), ("reference.json", reference),
                        ("mobility_config.json", mobility), ("manifest.json", manifest)]:
        (output_dir / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", required=True, type=Path)
    parser.add_argument("--parameter-file", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--normal-source", choices=["native_enormal", "electric_projection"], default="native_enormal")
    args = parser.parse_args()
    print(json.dumps(prepare(args.export_dir, args.parameter_file, args.output_dir, args.normal_source), indent=2))


if __name__ == "__main__":
    main()
