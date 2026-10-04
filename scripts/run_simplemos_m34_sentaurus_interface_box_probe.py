#!/usr/bin/env python3
"""Run and analyze the SimpleMOS M34 Sentaurus interface box probe."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m34_sentaurus_interface_box_probe_contract_v1.json"
)
TDR = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8_upstream_tdrs/n23_fps.tdr"
)
OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m34_interface_box_probe"
)
IMPORTER = REPO / "build-release/sentaurus_import.exe"
MATERIALS = (
    REPO / "reference_tcad/transportmodels_sentaurus2022/vela"
    / "materials_sentaurus2022.json"
)
REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/"
    "simplemos_m34_interface_box_probe_20260831_v1"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO).as_posix()
    except ValueError:
        return str(resolved)


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (
            Path(os.environ.get("SystemRoot", r"C:\Windows"))
            / "System32" / "OpenSSH" / f"{name}.exe"
        )
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv),
        cwd=REPO,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    return completed.stdout or ""


def deck() -> str:
    return '''File {
  Grid="input_fps.tdr"
  Plot="m34_box_probe.tdr"
  Current="m34_box_probe"
  Output="m34_box_probe.log"
}

Electrode {
  { Name="source" Voltage=0.0 }
  { Name="drain" Voltage=0.0 }
  { Name="gate" Voltage=0.05 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  EffectiveIntrinsicDensity(OldSlotboom)
}

Plot {
  BM_AngleElements
  BM_CoeffIntersectionNonDelaunayElements
  BM_ElementVolume
  BM_IntersectionNonDelaunayElements
  BM_VolumeIntersectionNonDelaunayElements
}

Math {
  AverageBoxMethod
  BoxMeasureFromFile(GrdNumbering)
}

Solve {
  Coupled(Iterations=1) { Poisson }
  Plot(FilePrefix="m34_box_probe")
}
'''


def prepare(contract: dict[str, Any], output: Path, importer: Path) -> dict[str, Any]:
    if sha256(TDR) != contract["input_tdr_sha256"]:
        raise ValueError("n23 input TDR hash does not match the M34 contract")
    bundle = output / "sentaurus_bundle/n23"
    bundle.mkdir(parents=True, exist_ok=True)
    shutil.copy2(TDR, bundle / "input_fps.tdr")
    deck_path = bundle / "m34_box_probe_des.cmd"
    deck_path.write_text(deck(), encoding="utf-8", newline="\n")

    import_dir = output / "n23_import"
    if import_dir.exists():
        shutil.rmtree(import_dir)
    run([
        str(importer), "--tdr", str(TDR),
        "--inventory-json", str(output / "n23_inventory.json"),
        "--export-dir", str(import_dir),
    ])
    manifest = {
        "schema": "vela.simplemos.sdevice.m34_probe_bundle.v1",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "input_tdr": portable(TDR),
        "input_tdr_sha256": sha256(TDR),
        "deck": portable(deck_path),
        "deck_sha256": sha256(deck_path),
        "expected_outputs": [
            "MeasureCoefficients.debug",
            "m34_box_probe.console.log",
            "m34_box_probe.log_des.log",
            "m34_box_probe_des.tdr",
        ],
    }
    write_json(output / "probe_bundle_manifest.json", manifest)
    return manifest


def run_sentaurus(
    output: Path,
    ssh_target: str,
    ssh_bin: str,
    scp_bin: str,
    remote_root: str,
    expected_release: str,
) -> str:
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True,
    ).strip()
    if expected_release not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    probe_exists = run(
        [ssh_bin, ssh_target, f"if test -e {remote_root}; then echo EXISTS; fi"],
        capture=True,
    ).strip()
    if probe_exists:
        raise FileExistsError(f"remote M34 root already exists: {remote_root}")
    run([ssh_bin, ssh_target, f"mkdir -p {remote_root}"])
    run([
        scp_bin, "-r", str(output / "sentaurus_bundle"),
        f"{ssh_target}:{remote_root}/",
    ])
    remote_case = f"{remote_root}/sentaurus_bundle/n23"
    command = (
        f"set -eu; cd {remote_case}; "
        "test ! -e MeasureCoefficients.debug; "
        "sdevice m34_box_probe_des.cmd > m34_box_probe.console.log 2>&1"
    )
    run([ssh_bin, ssh_target, command])
    archive_name = "simplemos_m34_interface_box_probe_results.tgz"
    run([
        ssh_bin, ssh_target,
        f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle",
    ])
    raw = output / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (output / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n"
    )
    return banner


RECORD = re.compile(
    r"^\s*(?P<grd>\d+)\s+(?P<des>-?\d+)\s+(?P<type>\d+)\s+"
    r"(?P<values>[^#]+?)\s*$"
)


def parse_debug_info(text: str) -> dict[str, int]:
    result = {}
    for key in ("dimension", "nb_vertices", "nb_grd_elements", "nb_des_elements"):
        match = re.search(rf"\b{key}\s*=\s*(\d+)", text)
        if match is None:
            raise ValueError(f"debug Info is missing {key}")
        result[key] = int(match.group(1))
    return result


def parse_debug_block(text: str, name: str) -> dict[int, dict[str, Any]]:
    match = re.search(
        rf"\n\s*{re.escape(name)}\s*\{{.*?\n(.*?)\n\s*\}}",
        text,
        re.DOTALL,
    )
    if match is None:
        raise ValueError(f"MeasureCoefficients.debug is missing {name}")
    result: dict[int, dict[str, Any]] = {}
    for line in match.group(1).splitlines():
        parsed = RECORD.match(line)
        if parsed is None:
            continue
        values = [float(value) for value in parsed.group("values").split()]
        result[int(parsed.group("grd"))] = {
            "des": int(parsed.group("des")),
            "type": int(parsed.group("type")),
            "values": values,
        }
    return result


def parse_log_stats(text: str) -> dict[str, Any]:
    patterns = {
        "input_grid_points": r"Number of grid points is\s+(\d+)",
        "internal_edges": r"NumberOfEdges\s*=\s*(\d+)",
        "geometrical_edges": r"NumberOfGeometricalEdges\s*=\s*(\d+)",
        "double_edges": r"NumberOfDoubleEdges\s*=\s*(\d+)",
        "internal_vertices": r"NumberOfVertices\s*=\s*(\d+)",
        "internal_elements": r"NumberOfElements\s*=\s*(\d+)",
        "obtuse_elements": r"NumberOfObtuseElements\s*=\s*(\d+)",
        "non_delaunay_elements": r"NumberOfNonDelaunayElements\s*=\s*(\d+)",
    }
    result: dict[str, Any] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match is None:
            raise ValueError(f"Sentaurus log is missing {key}")
        result[key] = int(match.group(1))
    result["average_box_method"] = "CVPL_AverageBoxMethod = TRUE" in text
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def area(points: list[tuple[float, float]]) -> float:
    return 0.5 * abs(
        (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
        - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1])
    )


def raw_coefficient(points: list[tuple[float, float]], local: int) -> float:
    a = points[local]
    b = points[(local + 1) % 3]
    c = points[(local + 2) % 3]
    ux, uy = b[0] - a[0], b[1] - a[1]
    vx, vy = c[0] - a[0], c[1] - a[1]
    return 0.5 * (ux * vx + uy * vy) / abs(ux * vy - uy * vx)


def vela_positive_coefficient(points: list[tuple[float, float]], local: int) -> float:
    raw = raw_coefficient(points, local)
    if raw >= 0.0:
        return raw
    b = points[(local + 1) % 3]
    c = points[(local + 2) % 3]
    length2 = (b[0] - c[0]) ** 2 + (b[1] - c[1]) ** 2
    return area(points) / (3.0 * length2)


def signed_average_box_measures(
    points: list[tuple[float, float]],
) -> list[float]:
    def distance_squared(a: int, b: int) -> float:
        return sum((points[b][axis] - points[a][axis]) ** 2 for axis in (0, 1))

    result = []
    for local in range(3):
        j = (local + 1) % 3
        k = (local + 2) % 3
        result.append(0.25 * (
            distance_squared(local, k) * raw_coefficient(points, j)
            + distance_squared(local, j) * raw_coefficient(points, k)
        ))
    return result


def infer_input_to_debug_local_permutation(
    element_ids: set[int],
    elements: dict[int, dict[str, Any]],
    nodes: dict[int, tuple[float, float]],
    coefficients: dict[int, dict[str, Any]],
) -> tuple[tuple[int, int, int], float]:
    """Map each input-TDR local vertex to its debug-file local slot."""
    candidates: list[tuple[float, tuple[int, int, int]]] = []
    for permutation in itertools.permutations(range(3)):
        maximum = 0.0
        for element in element_ids:
            points = [nodes[node] for node in elements[element]["nodes"]]
            for local in range(3):
                maximum = max(
                    maximum,
                    abs(
                        float(coefficients[element]["values"][permutation[local]])
                        - raw_coefficient(points, local)
                    ),
                )
        candidates.append((maximum, permutation))
    maximum, permutation = min(candidates)
    return permutation, maximum


def infer_input_to_debug_measure_permutation(
    element_ids: set[int],
    elements: dict[int, dict[str, Any]],
    nodes: dict[int, tuple[float, float]],
    measures: dict[int, dict[str, Any]],
) -> tuple[tuple[int, int, int], float]:
    """Map input local vertices to Measure slots using signed AverageBox shares."""
    candidates: list[tuple[float, tuple[int, int, int]]] = []
    for permutation in itertools.permutations(range(3)):
        maximum = 0.0
        for element in element_ids:
            points = [nodes[node] for node in elements[element]["nodes"]]
            expected = signed_average_box_measures(points)
            for local in range(3):
                maximum = max(
                    maximum,
                    abs(float(measures[element]["values"][permutation[local]])
                        - expected[local]),
                )
        candidates.append((maximum, permutation))
    maximum, permutation = min(candidates)
    return permutation, maximum


def material_eps() -> dict[str, float]:
    values = {item["name"]: float(item["eps_r"])
              for item in read_json(MATERIALS)["materials"]}
    return {"Si": values["Si"], "Oxide": values["SiO2"]}


def analyze(contract: dict[str, Any], output: Path, banner: str) -> dict[str, Any]:
    raw = output / "sentaurus_raw/sentaurus_bundle/n23"
    debug_path = raw / "MeasureCoefficients.debug"
    log_path = raw / "m34_box_probe.console.log"
    debug_text = debug_path.read_text(encoding="utf-8", errors="replace")
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    info = parse_debug_info(debug_text)
    measures = parse_debug_block(debug_text, "Measure")
    coefficients = parse_debug_block(debug_text, "Coefficients")
    log_stats = parse_log_stats(log_text)

    imported = output / "n23_import"
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_rows(imported / "nodes.csv")
    }
    elements: dict[int, dict[str, Any]] = {}
    for row in read_rows(imported / "elements.csv"):
        element = int(row["id"])
        elements[element] = {
            "nodes": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
            "region": row["region"],
            "material": row["material"],
        }
    contacts = read_rows(imported / "contacts.csv")
    contact_nodes = {
        int(node)
        for row in contacts
        for node in row["node_ids"].split(";")
        if node
    }

    triangle_debug_ids = {
        grd for grd, row in measures.items()
        if row["type"] == 2 and row["des"] >= 0
    }
    if triangle_debug_ids != set(elements):
        missing = sorted(set(elements) - triangle_debug_ids)[:10]
        extra = sorted(triangle_debug_ids - set(elements))[:10]
        raise ValueError(f"incomplete triangle mapping: missing={missing}, extra={extra}")
    if set(measures) != set(coefficients):
        raise ValueError("Measure and Coefficients grid-element supports differ")

    edge_cells: dict[tuple[int, int], list[tuple[int, int]]] = {}
    node_elements: dict[int, list[tuple[int, int]]] = {node: [] for node in nodes}
    for element, item in elements.items():
        tri = item["nodes"]
        for local, node in enumerate(tri):
            node_elements[node].append((element, local))
            opposite = tuple(sorted((tri[(local + 1) % 3], tri[(local + 2) % 3])))
            edge_cells.setdefault(opposite, []).append((element, local))

    interface_adjacency = {
        edge: adjacent
        for edge, adjacent in edge_cells.items()
        if {elements[element]["material"] for element, _ in adjacent}
        == {"Si", "SiO2"}
    }
    interface_element_ids = {
        element for adjacent in interface_adjacency.values()
        for element, _ in adjacent
    }
    coefficient_permutation, coefficient_permutation_max_error = (
        infer_input_to_debug_local_permutation(
            interface_element_ids, elements, nodes, coefficients
        )
    )
    measure_permutation, measure_permutation_max_error = (
        infer_input_to_debug_measure_permutation(
            interface_element_ids, elements, nodes, measures
        )
    )

    eps = material_eps()
    edge_rows: list[dict[str, Any]] = []
    interface_nodes: set[int] = set()
    for edge, adjacent in sorted(interface_adjacency.items()):
        by_material = {elements[element]["material"]: (element, local)
                       for element, local in adjacent}
        interface_nodes.update(edge)
        values: dict[str, dict[str, float | int]] = {}
        for material, key in (("Si", "Si"), ("SiO2", "Oxide")):
            element, local = by_material[material]
            points = [nodes[node] for node in elements[element]["nodes"]]
            values[key] = {
                "element": element,
                "local": local,
                "sent": float(
                    coefficients[element]["values"][coefficient_permutation[local]]
                ),
                "raw": raw_coefficient(points, local),
                "vela": vela_positive_coefficient(points, local),
            }
        si = values["Si"]
        oxide = values["Oxide"]
        sent_transport = float(si["sent"])
        sent_poisson = eps["Si"] * sent_transport + eps["Oxide"] * float(oxide["sent"])
        vela_transport = float(si["vela"])
        vela_poisson = eps["Si"] * vela_transport + eps["Oxide"] * float(oxide["vela"])
        legacy_poisson = 0.5 * (eps["Si"] + eps["Oxide"]) * (
            float(si["vela"]) + float(oxide["vela"])
        )
        edge_rows.append({
            "node0": edge[0],
            "node1": edge[1],
            "si_element": si["element"],
            "oxide_element": oxide["element"],
            "sentaurus_si_local_coefficient": si["sent"],
            "sentaurus_oxide_local_coefficient": oxide["sent"],
            "raw_si_local_coefficient": si["raw"],
            "raw_oxide_local_coefficient": oxide["raw"],
            "vela_si_local_coefficient": si["vela"],
            "vela_oxide_local_coefficient": oxide["vela"],
            "sentaurus_total_over_si_transport_coefficient": (
                (float(si["sent"]) + float(oxide["sent"])) / sent_transport
            ),
            "sentaurus_region_poisson_coefficient": sent_poisson,
            "vela_region_poisson_coefficient": vela_poisson,
            "vela_legacy_poisson_coefficient": legacy_poisson,
            "vela_legacy_over_sentaurus_region_poisson": legacy_poisson / sent_poisson,
            "vela_region_over_sentaurus_region_poisson": vela_poisson / sent_poisson,
        })

    node_rows: list[dict[str, Any]] = []
    for node in sorted(interface_nodes):
        sums = {"Si": 0.0, "SiO2": 0.0, "other": 0.0}
        bary = {"Si": 0.0, "SiO2": 0.0, "other": 0.0}
        for element, local in node_elements[node]:
            material = elements[element]["material"]
            bucket = material if material in ("Si", "SiO2") else "other"
            sums[bucket] += float(
                measures[element]["values"][measure_permutation[local]]
            )
            points = [nodes[item] for item in elements[element]["nodes"]]
            bary[bucket] += area(points) / 3.0
        node_rows.append({
            "node": node,
            "x_um": nodes[node][0],
            "y_um": nodes[node][1],
            "sentaurus_si_measure_um2": sums["Si"],
            "sentaurus_oxide_measure_um2": sums["SiO2"],
            "sentaurus_other_measure_um2": sums["other"],
            "sentaurus_total_over_si_measure": sum(sums.values()) / sums["Si"],
            "barycentric_si_measure_um2": bary["Si"],
            "barycentric_oxide_measure_um2": bary["SiO2"],
            "barycentric_total_over_si_measure": sum(bary.values()) / bary["Si"],
        })

    coordinate_groups: dict[tuple[float, float], list[int]] = {}
    for node, point in nodes.items():
        coordinate_groups.setdefault(point, []).append(node)
    duplicate_groups = [group for group in coordinate_groups.values() if len(group) > 1]
    internal_vertex_delta = log_stats["internal_vertices"] - len(nodes)
    internal_edge_delta = log_stats["internal_edges"] - log_stats["geometrical_edges"]
    contact_explains_vertex_delta = internal_vertex_delta == len(contact_nodes)
    no_material_double_node_evidence = (
        contact_explains_vertex_delta
        and log_stats["double_edges"] == 0
        and not duplicate_groups
    )

    des_ids = [row["des"] for row in measures.values()
               if row["type"] == 2 and row["des"] >= 0]
    identity_count = sum(
        grd == row["des"] for grd, row in measures.items()
        if row["type"] == 2 and row["des"] >= 0
    )
    si_raw_delta = [
        abs(float(row["sentaurus_si_local_coefficient"])
            - float(row["raw_si_local_coefficient"]))
        for row in edge_rows
    ]
    oxide_raw_delta = [
        abs(float(row["sentaurus_oxide_local_coefficient"])
            - float(row["raw_oxide_local_coefficient"]))
        for row in edge_rows
    ]
    region_poisson_ratio = sorted(
        float(row["vela_region_over_sentaurus_region_poisson"])
        for row in edge_rows
    )
    legacy_poisson_ratio = sorted(
        float(row["vela_legacy_over_sentaurus_region_poisson"])
        for row in edge_rows
    )
    node_measure_ratio = sorted(
        float(row["sentaurus_total_over_si_measure"]) for row in node_rows
    )
    bary_measure_ratio = sorted(
        float(row["barycentric_total_over_si_measure"]) for row in node_rows
    )
    median = lambda values: values[len(values) // 2]

    mapping_rows = [{
        "input_grid_vertices": len(nodes),
        "debug_grid_vertices": info["nb_vertices"],
        "sdevice_internal_vertices": log_stats["internal_vertices"],
        "internal_vertex_delta": internal_vertex_delta,
        "unique_contact_nodes": len(contact_nodes),
        "contact_nodes_explain_delta": contact_explains_vertex_delta,
        "material_interface_nodes": len(interface_nodes),
        "input_duplicate_coordinate_groups": len(duplicate_groups),
        "geometrical_edges": log_stats["geometrical_edges"],
        "internal_edges": log_stats["internal_edges"],
        "internal_edge_delta": internal_edge_delta,
        "double_edges": log_stats["double_edges"],
        "no_material_double_node_evidence": no_material_double_node_evidence,
    }]
    report = {
        "schema": "vela.simplemos.sdevice.m34_sentaurus_interface_box_probe_report.v1",
        "status": "qualified",
        "sentaurus_banner": banner,
        "input_tdr_sha256": sha256(TDR),
        "debug_sha256": sha256(debug_path),
        "log_sha256": sha256(log_path),
        "probe": {
            "debug_info": info,
            "log_stats": log_stats,
            "measure_record_count": len(measures),
            "coefficient_record_count": len(coefficients),
            "triangle_record_count": len(triangle_debug_ids),
            "triangle_grd_des_identity_count": identity_count,
            "triangle_des_ids_complete": sorted(des_ids) == list(range(len(elements))),
            "input_to_debug_local_permutation": list(coefficient_permutation),
            "input_to_debug_coefficient_local_permutation": list(
                coefficient_permutation
            ),
            "input_to_debug_measure_local_permutation": list(
                measure_permutation
            ),
            "interface_permutation_max_abs_coefficient_error": (
                coefficient_permutation_max_error
            ),
            "interface_permutation_max_abs_measure_error": (
                measure_permutation_max_error
            ),
        },
        "vertex_mapping": mapping_rows[0],
        "si_oxide_interface": {
            "edge_count": len(edge_rows),
            "node_count": len(node_rows),
            "sentaurus_vs_raw_max_abs_si_coefficient": max(si_raw_delta),
            "sentaurus_vs_raw_max_abs_oxide_coefficient": max(oxide_raw_delta),
            "vela_region_over_sentaurus_region_poisson_min": min(region_poisson_ratio),
            "vela_region_over_sentaurus_region_poisson_median": median(region_poisson_ratio),
            "vela_region_over_sentaurus_region_poisson_max": max(region_poisson_ratio),
            "vela_legacy_over_sentaurus_region_poisson_min": min(legacy_poisson_ratio),
            "vela_legacy_over_sentaurus_region_poisson_median": median(legacy_poisson_ratio),
            "vela_legacy_over_sentaurus_region_poisson_max": max(legacy_poisson_ratio),
            "sentaurus_total_over_si_measure_min": min(node_measure_ratio),
            "sentaurus_total_over_si_measure_median": median(node_measure_ratio),
            "sentaurus_total_over_si_measure_max": max(node_measure_ratio),
            "barycentric_total_over_si_measure_min": min(bary_measure_ratio),
            "barycentric_total_over_si_measure_median": median(bary_measure_ratio),
            "barycentric_total_over_si_measure_max": max(bary_measure_ratio),
        },
        "conclusions": {
            "input_si_oxide_topology": "shared_node_ids",
            "internal_vertex_expansion": (
                "fully_explained_by_unique_contact_nodes"
                if contact_explains_vertex_delta else "not_closed"
            ),
            "explicit_si_oxide_double_nodes_supported": False,
            "region_local_box_contributions_supported": True,
            "production_default_changed": False,
        },
        "claim_policy": contract["claim_policy"],
    }
    if len(edge_rows) < contract["acceptance"]["minimum_interface_edges"]:
        raise ValueError("M34 did not find the contracted Si/Oxide interface")
    if not log_stats["average_box_method"]:
        raise ValueError("M34 probe did not use the original AverageBoxMethod")
    for row in itertools.chain(edge_rows, node_rows):
        for value in row.values():
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError("M34 produced a non-finite diagnostic value")

    analysis = output / "analysis"
    write_csv(analysis / "m34_interface_edge_coefficients.csv", edge_rows)
    write_csv(analysis / "m34_interface_node_measures.csv", node_rows)
    write_csv(analysis / "m34_vertex_mapping_summary.csv", mapping_rows)
    write_json(analysis / "m34_sentaurus_interface_box_probe_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--importer", type=Path, default=IMPORTER)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    if not any((args.prepare, args.live_sentaurus, args.analyze)):
        args.analyze = True

    contract = read_json(args.contract)
    output = args.output.resolve()
    manifest_path = output / "probe_bundle_manifest.json"
    if args.prepare:
        manifest = prepare(contract, output, args.importer.resolve())
    else:
        manifest = read_json(manifest_path)

    banner_path = output / "sentaurus_banner.txt"
    if args.live_sentaurus:
        banner = run_sentaurus(
            output,
            args.ssh_target,
            args.ssh_bin,
            args.scp_bin,
            args.remote_root,
            contract["sentaurus_release"],
        )
    elif args.analyze:
        banner = banner_path.read_text(encoding="utf-8").strip()
    else:
        banner = ""
    report = analyze(contract, output, banner) if args.analyze else None
    print(json.dumps({
        "manifest": manifest,
        "report": report,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
