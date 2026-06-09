#!/usr/bin/env python3
"""Compare Sentaurus TDR neutral exports with tdx-generated DF-ISE text files."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def read_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text())


def parse_count(text: str, name: str) -> int | None:
    match = re.search(rf"\b{name}\s*=\s*(\d+)", text)
    return int(match.group(1)) if match else None


def parse_ise_list(text: str, name: str) -> list[str]:
    match = re.search(rf"\b{name}\s*=\s*\[(.*?)\]", text, re.S)
    if not match:
        return []
    raw = match.group(1)
    return [token[0] or token[1] for token in re.findall(r'"([^"]+)"|([A-Za-z0-9_.:+-]+)', raw)]


def parse_grd(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace")
    regions = parse_ise_list(text, "regions")
    materials = parse_ise_list(text, "materials")
    region_materials = {
        region: materials[index]
        for index, region in enumerate(regions)
        if index < len(materials)
    }
    region_blocks: dict[str, dict[str, Any]] = {}
    matches = list(re.finditer(r'Region\s*\("([^"]+)"\)\s*\{', text))
    for index, match in enumerate(matches):
        name = match.group(1)
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        body = text[match.end():end]
        material_match = re.search(r"\bmaterial\s*=\s*([A-Za-z0-9_.:+-]+)", body)
        elements_match = re.search(r"\bElements\s*\((\d+)\)", body)
        material = material_match.group(1) if material_match else region_materials.get(name)
        region_blocks[name] = {
            "material": material,
            "elements": int(elements_match.group(1)) if elements_match else None,
        }
    contact_counts = {
        name: int(info["elements"])
        for name, info in region_blocks.items()
        if str(info.get("material", "")).lower() == "contact" and info.get("elements") is not None
    }
    return {
        "path": str(path),
        "vertex_count": parse_count(text, "nb_vertices"),
        "element_count": parse_count(text, "nb_elements"),
        "region_count": parse_count(text, "nb_regions"),
        "regions": regions,
        "materials": materials,
        "region_blocks": region_blocks,
        "contact_element_counts": contact_counts,
    }


def parse_dat(path: Path) -> dict[str, dict[str, Any]]:
    text = path.read_text(errors="replace")
    datasets: dict[str, dict[str, Any]] = {}
    pattern = re.compile(
        r'Dataset\s*\("([^"]+)"\)\s*\{(?P<header>.*?)\bValues\s*\((\d+)\)\s*\{(?P<values>.*?)\n\s*\}',
        re.S,
    )
    for match in pattern.finditer(text):
        name = match.group(1)
        header = match.group("header")
        dimension_match = re.search(r"\bdimension\s*=\s*(\d+)", header)
        dimension = int(dimension_match.group(1)) if dimension_match else 1
        values = [float(item) for item in NUMBER_RE.findall(match.group("values"))]
        count = int(match.group(3))
        datasets[name] = {
            "count": count,
            "dimension": dimension,
            "values": values,
        }
    return datasets


def count_csv_rows(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open(newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def load_field_csv(path: Path) -> list[list[float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        component_columns = sorted(
            (name for name in reader.fieldnames if name.startswith("component")),
            key=lambda name: int(name.removeprefix("component") or "0"),
        )
        if not component_columns:
            component_columns = [name for name in reader.fieldnames if name != "node_id"]
        rows = []
        for row in reader:
            rows.append([float(row[name]) for name in component_columns])
        return rows


def find_tdr_field(root: Path, field: str) -> Path | None:
    fields_root = root / "fields"
    exact = fields_root / f"{field}_region0.csv"
    if exact.is_file():
        return exact
    matches = sorted(fields_root.glob(f"{field}_region*.csv"))
    return matches[0] if matches else None


def compare_flat_values(tdr_rows: list[list[float]], tdx_dataset: dict[str, Any]) -> dict[str, Any]:
    tdr_values = [component for row in tdr_rows for component in row]
    tdx_values = list(tdx_dataset.get("values", []))
    points = min(len(tdr_values), len(tdx_values))
    max_abs = 0.0
    max_rel = 0.0
    worst: dict[str, Any] | None = None
    for index in range(points):
        abs_diff = abs(tdr_values[index] - tdx_values[index])
        rel_diff = abs_diff / max(abs(tdx_values[index]), 1.0e-300)
        if abs_diff > max_abs:
            max_abs = abs_diff
            worst = {
                "index": index,
                "tdr": tdr_values[index],
                "tdx": tdx_values[index],
                "abs_diff": abs_diff,
                "rel_diff": rel_diff,
            }
        max_rel = max(max_rel, rel_diff)
    return {
        "points_compared": points,
        "tdr_values": len(tdr_values),
        "tdx_values": len(tdx_values),
        "max_abs_diff": max_abs,
        "max_rel_diff": max_rel,
        "worst": worst,
    }


def tdr_contact_counts(root: Path) -> dict[str, int]:
    metadata = read_json_if_exists(root / "metadata.json")
    counts: dict[str, int] = {}
    for region in metadata.get("regions", []):
        edges = region.get("edges")
        material = str(region.get("material", ""))
        if edges is not None and (material.lower() == "contact" or int(region.get("triangles", 0)) == 0):
            counts[str(region.get("name"))] = int(edges)
    return counts


def status_entry(tdr: int | None, tdx: int | None) -> dict[str, Any]:
    return {"tdr": tdr, "tdx": tdx, "match": tdr == tdx if tdr is not None and tdx is not None else None}


def non_contact_element_count(grd: dict[str, Any]) -> int | None:
    total = 0
    found = False
    for info in grd.get("region_blocks", {}).values():
        if str(info.get("material", "")).lower() == "contact":
            continue
        elements = info.get("elements")
        if elements is None:
            continue
        total += int(elements)
        found = True
    return total if found else grd.get("element_count")


def tdr_total_element_count(root: Path, contact_counts: dict[str, int]) -> int | None:
    bulk_count = count_csv_rows(root / "elements.csv")
    if bulk_count is None:
        return None
    return bulk_count + sum(contact_counts.values())


def choose_first(root: Path, pattern: str, preferred: list[str] | None = None) -> Path:
    if preferred:
        for name in preferred:
            candidate = root / name
            if candidate.is_file():
                return candidate
    matches = sorted(root.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"no {pattern} found in {root}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdr-export", required=True, type=Path)
    parser.add_argument("--tdx-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-abs-diff", type=float, default=0.0)
    parser.add_argument("--max-rel-diff", type=float, default=1.0e-12)
    args = parser.parse_args()

    grd = parse_grd(choose_first(args.tdx_dir, "*.grd", ["pn2d_msh.grd"]))
    dat = parse_dat(choose_first(args.tdx_dir, "*.dat", ["pn2d_msh.dat"]))
    metadata = read_json_if_exists(args.tdr_export / "metadata.json")
    tdr_regions = [
        str(region.get("name"))
        for region in metadata.get("regions", [])
        if region.get("name") is not None
    ]
    if not tdr_regions and (args.tdr_export / "elements.csv").is_file():
        with (args.tdr_export / "elements.csv").open(newline="") as handle:
            reader = csv.DictReader(handle)
            tdr_regions = sorted({row.get("region", "") for row in reader if row.get("region")})
    if (args.tdr_export / "contacts.csv").is_file():
        with (args.tdr_export / "contacts.csv").open(newline="") as handle:
            reader = csv.DictReader(handle)
            contact_names = {row.get("name", "") for row in reader if row.get("name")}
        tdr_regions = sorted(set(tdr_regions) | contact_names)
    contact_counts = tdr_contact_counts(args.tdr_export)

    mesh = {
        "vertex_count": status_entry(
            int(metadata.get("vertex_count")) if metadata.get("vertex_count") is not None else count_csv_rows(args.tdr_export / "nodes.csv"),
            grd.get("vertex_count"),
        ),
        "element_count": status_entry(tdr_total_element_count(args.tdr_export, contact_counts), grd.get("element_count")),
        "bulk_element_count": status_entry(
            count_csv_rows(args.tdr_export / "elements.csv"),
            non_contact_element_count(grd),
        ),
        "region_names": {
            "tdr": sorted(tdr_regions),
            "tdx": sorted(grd.get("regions", [])),
            "match": sorted(tdr_regions) == sorted(grd.get("regions", [])) if tdr_regions else None,
        },
        "contact_element_counts": {},
    }
    for name in sorted(set(contact_counts) | set(grd.get("contact_element_counts", {}))):
        mesh["contact_element_counts"][name] = {
            "tdr": contact_counts.get(name),
            "tdx": grd.get("contact_element_counts", {}).get(name),
            "match": (
                contact_counts.get(name) == grd.get("contact_element_counts", {}).get(name)
                if name in contact_counts
                else None
            ),
        }

    fields: dict[str, Any] = {}
    for field, dataset in sorted(dat.items()):
        field_path = find_tdr_field(args.tdr_export, field)
        if field_path is None:
            continue
        report = compare_flat_values(load_field_csv(field_path), dataset)
        report["status"] = "pass"
        if (
            report["tdr_values"] != report["tdx_values"]
            or (
                report["max_abs_diff"] > args.max_abs_diff
                and report["max_rel_diff"] > args.max_rel_diff
            )
        ):
            report["status"] = "fail"
        fields[field] = report

    failed = any(
        entry.get("match") is False
        for key, entry in mesh.items()
        if isinstance(entry, dict) and key != "contact_element_counts"
    )
    failed = failed or any(report.get("status") == "fail" for report in fields.values())
    report = {
        "status": "fail" if failed else "pass",
        "mesh": mesh,
        "fields": fields,
        "tdx": {
            "grd": grd["path"],
            "datasets": sorted(dat),
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "tdr_tdx_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
