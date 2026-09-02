#!/usr/bin/env python3
"""Compare Genius NPN BJT current-density and recombination fields."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


Q_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sentaurus_scalar(path: Path) -> list[float]:
    return [float(row["component0"]) for row in read_csv(path)]


def sentaurus_vector(path: Path) -> list[tuple[float, float]]:
    return [
        (float(row["component0"]), float(row["component1"]))
        for row in read_csv(path)
    ]


def read_vtk_point_data(
    path: Path,
) -> tuple[int, dict[str, list[float]], dict[str, list[tuple[float, float, float]]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    point_line = next(line for line in lines if line.startswith("POINT_DATA "))
    count = int(point_line.split()[1])
    start = lines.index(point_line) + 1
    scalars: dict[str, list[float]] = {}
    vectors: dict[str, list[tuple[float, float, float]]] = {}
    index = start
    while index < len(lines):
        line = lines[index]
        if line.startswith("SCALARS "):
            name = line.split()[1]
            if index + 1 >= len(lines) or not lines[index + 1].startswith("LOOKUP_TABLE"):
                raise ValueError(f"VTK scalar {name} is missing LOOKUP_TABLE")
            scalars[name] = [float(value) for value in lines[index + 2 : index + 2 + count]]
            index += count + 2
            continue
        if line.startswith("VECTORS "):
            name = line.split()[1]
            values = []
            for row in lines[index + 1 : index + 1 + count]:
                components = [float(value) for value in row.split()]
                if len(components) != 3:
                    raise ValueError(f"VTK vector {name} does not have three components")
                values.append((components[0], components[1], components[2]))
            vectors[name] = values
            index += count + 1
            continue
        index += 1
    return count, scalars, vectors


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot calculate percentile of an empty population")
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def scalar_error_statistics(values: list[float]) -> dict[str, float | int]:
    absolute = [abs(value) for value in values]
    return {
        "node_count": len(absolute),
        "rmse": math.sqrt(sum(value * value for value in absolute) / len(absolute)),
        "p95_absolute_error": percentile(absolute, 0.95),
        "maximum_absolute_error": max(absolute),
    }


def vector_metrics(
    reference: list[tuple[float, float]],
    actual: list[tuple[float, float, float]],
    reference_fraction: float,
) -> dict[str, object]:
    if len(reference) != len(actual):
        raise ValueError("current-density vector length mismatch")
    reference_magnitude = [math.hypot(x, y) for x, y in reference]
    actual_magnitude = [math.hypot(x, y) for x, y, _ in actual]
    reference_peak = max(reference_magnitude)
    floor = reference_peak * reference_fraction
    selected = [value >= floor for value in reference_magnitude]
    log_errors = [
        math.log10(actual_value / reference_value)
        for actual_value, reference_value, keep in zip(
            actual_magnitude, reference_magnitude, selected, strict=True
        )
        if keep and actual_value > 0.0 and reference_value > 0.0
    ]
    squared_error = 0.0
    squared_reference = 0.0
    dot = 0.0
    squared_actual = 0.0
    for (rx, ry), (ax, ay, _), keep in zip(reference, actual, selected, strict=True):
        if not keep:
            continue
        squared_error += (ax - rx) ** 2 + (ay - ry) ** 2
        squared_reference += rx * rx + ry * ry
        squared_actual += ax * ax + ay * ay
        dot += ax * rx + ay * ry
    return {
        "units": "A/cm^2",
        "mask": {
            "type": "sentaurus_reference_magnitude_fraction_of_peak",
            "fraction": reference_fraction,
            "minimum_A_per_cm2": floor,
        },
        "selected_node_count": sum(selected),
        "reference_peak_A_per_cm2": reference_peak,
        "vela_peak_A_per_cm2": max(actual_magnitude),
        "log10_magnitude_error": scalar_error_statistics(log_errors),
        "normalized_vector_rmse": math.sqrt(squared_error / squared_reference),
        "global_vector_cosine_similarity": (
            dot / math.sqrt(squared_reference * squared_actual)
            if squared_reference > 0.0 and squared_actual > 0.0
            else None
        ),
    }


def lumped_node_areas_um2(
    nodes: list[tuple[float, float]], triangles: list[tuple[int, int, int]]
) -> list[float]:
    areas = [0.0] * len(nodes)
    for triangle in triangles:
        a, b, c = (nodes[index] for index in triangle)
        area = 0.5 * abs(
            (b[0] - a[0]) * (c[1] - a[1])
            - (c[0] - a[0]) * (b[1] - a[1])
        )
        for index in triangle:
            areas[index] += area / 3.0
    return areas


def source_metrics(
    reference: list[float],
    actual: list[float],
    node_areas_um2: list[float],
    reference_fraction: float,
) -> dict[str, object]:
    if not (len(reference) == len(actual) == len(node_areas_um2)):
        raise ValueError("source-rate field length mismatch")
    reference_peak = max(abs(value) for value in reference)
    floor = reference_peak * reference_fraction
    selected = [abs(value) >= floor for value in reference]
    log_errors = [
        math.log10(abs(value) / abs(target))
        for value, target, keep in zip(actual, reference, selected, strict=True)
        if keep and value != 0.0 and target != 0.0
    ]
    weights = [area * 1.0e-12 for area in node_areas_um2]
    reference_signed = sum(value * weight for value, weight in zip(reference, weights, strict=True))
    actual_signed = sum(value * weight for value, weight in zip(actual, weights, strict=True))
    reference_absolute = sum(abs(value) * weight for value, weight in zip(reference, weights, strict=True))
    actual_absolute = sum(abs(value) * weight for value, weight in zip(actual, weights, strict=True))
    normalized_l1 = sum(
        abs(value - target) * weight
        for value, target, weight in zip(actual, reference, weights, strict=True)
    ) / max(reference_absolute, 1.0e-300)
    shape_tv = 0.5 * sum(
        abs(
            abs(target) * weight / max(reference_absolute, 1.0e-300)
            - abs(value) * weight / max(actual_absolute, 1.0e-300)
        )
        for value, target, weight in zip(actual, reference, weights, strict=True)
    )
    return {
        "units": "cm^-3 s^-1",
        "mask": {
            "type": "sentaurus_reference_absolute_fraction_of_peak",
            "fraction": reference_fraction,
            "minimum_cm3_per_s": floor,
        },
        "selected_node_count": sum(selected),
        "reference_peak_absolute_cm3_per_s": reference_peak,
        "vela_peak_absolute_cm3_per_s": max(abs(value) for value in actual),
        "log10_magnitude_error": scalar_error_statistics(log_errors),
        "signed_integral_A_per_um": {
            "sentaurus": Q_C * reference_signed,
            "vela": Q_C * actual_signed,
            "ratio_vela_over_sentaurus": (
                actual_signed / reference_signed if reference_signed else None
            ),
        },
        "absolute_integral_A_per_um": {
            "sentaurus": Q_C * reference_absolute,
            "vela": Q_C * actual_absolute,
            "ratio_vela_over_sentaurus": (
                actual_absolute / reference_absolute if reference_absolute else None
            ),
        },
        "normalized_l1_error": normalized_l1,
        "absolute_shape_total_variation": shape_tv,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_markdown(path: Path, report: dict[str, object]) -> None:
    electron = report["current_density"]["electron"]
    hole = report["current_density"]["hole"]
    srh = report["recombination"]["srh"]
    auger = report["recombination"]["auger"]
    lines = [
        "# Genius NPN BJT transport/source spatial characterization",
        "",
        "These quantities are diagnostic characterization, not asserted acceptance gates.",
        "Current-density comparison uses Vela's Sentaurus-style nodal quasi-Fermi-gradient reconstruction in A/cm^2.",
        "Recombination integrals use the common triangular mesh and a 1 um out-of-plane depth.",
        "",
        "| Quantity | Selected nodes | P95 log-magnitude error [decade] | Normalized error |",
        "|---|---:|---:|---:|",
        f"| Electron current density | {electron['selected_node_count']} | {electron['log10_magnitude_error']['p95_absolute_error']:.6g} | vector RMSE {electron['normalized_vector_rmse']:.6g} |",
        f"| Hole current density | {hole['selected_node_count']} | {hole['log10_magnitude_error']['p95_absolute_error']:.6g} | vector RMSE {hole['normalized_vector_rmse']:.6g} |",
        f"| SRH recombination | {srh['selected_node_count']} | {srh['log10_magnitude_error']['p95_absolute_error']:.6g} | L1 {srh['normalized_l1_error']:.6g} |",
        f"| Auger recombination | {auger['selected_node_count']} | {auger['log10_magnitude_error']['p95_absolute_error']:.6g} | L1 {auger['normalized_l1_error']:.6g} |",
        "",
        f"Comparison status: **{report['status']}**.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--sentaurus-fields-root", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    threshold_document = json.loads(
        (args.reference_root / "contracts" / "comparison_thresholds.json").read_text(
            encoding="utf-8"
        )
    )
    contract = threshold_document["wp3_wp5_vela_comparison"]["transport_source_comparison"]
    if contract["status"] != "characterization_only":
        raise ValueError("transport/source comparison unexpectedly asserts a gate")
    fields_root = args.sentaurus_fields_root / "fields"
    nodes_rows = read_csv(args.sentaurus_fields_root / "nodes.csv")
    nodes = [(float(row["x_um"]), float(row["y_um"])) for row in nodes_rows]
    triangles = [
        (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        for row in read_csv(args.sentaurus_fields_root / "elements.csv")
    ]
    node_areas = lumped_node_areas_um2(nodes, triangles)
    count, vtk_scalars, vtk_vectors = read_vtk_point_data(args.vela_vtk)
    if count != len(nodes):
        raise ValueError(f"VTK node count {count} does not match Sentaurus {len(nodes)}")

    electron_field = contract["vela_current_density_fields"]["electron"]
    hole_field = contract["vela_current_density_fields"]["hole"]
    report = {
        "schema_version": 1,
        "device": "Genius NPN BJT",
        "bias": contract["bias"],
        "status": contract["status"],
        "reason": contract["reason"],
        "common_node_count": count,
        "current_density": {
            "electron": vector_metrics(
                sentaurus_vector(fields_root / "eCurrentDensity_region0.csv"),
                vtk_vectors[electron_field],
                float(contract["current_density_reference_fraction"]),
            ),
            "hole": vector_metrics(
                sentaurus_vector(fields_root / "hCurrentDensity_region0.csv"),
                vtk_vectors[hole_field],
                float(contract["current_density_reference_fraction"]),
            ),
        },
        "recombination": {
            "srh": source_metrics(
                sentaurus_scalar(fields_root / "srhRecombination_region0.csv"),
                vtk_scalars["SRHRecombinationCm3PerS"],
                node_areas,
                float(contract["source_rate_reference_fraction"]),
            ),
            "auger": source_metrics(
                sentaurus_scalar(fields_root / "AugerRecombination_region0.csv"),
                vtk_scalars["AugerRecombinationCm3PerS"],
                node_areas,
                float(contract["source_rate_reference_fraction"]),
            ),
        },
        "source_sha256": {
            "sentaurus_field_manifest": sha256(args.sentaurus_fields_root / "field_manifest.json"),
            "vela_vtk": sha256(args.vela_vtk),
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "transport_source_comparison.json"
    output.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_markdown(args.output_dir / "transport_source_comparison.md", report)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
