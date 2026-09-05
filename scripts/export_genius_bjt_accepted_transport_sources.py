#!/usr/bin/env python3
"""Export split recombination and current-density fields from accepted BJT states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from compare_genius_bjt_transport_fields import (
    lumped_node_areas_um2,
    read_csv,
    read_vtk_point_data,
    sentaurus_scalar,
    sentaurus_vector,
    source_metrics,
    vector_metrics,
)


Q_C = 1.602176634e-19
REPRESENTATIVE_INDICES = (0, 10, 20, 30)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vector_closure_error(
    total: tuple[float, float, float],
    drift: tuple[float, float, float],
    diffusion: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum((total[i] - drift[i] - diffusion[i]) ** 2 for i in range(3))
    )


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_point(
    vtk_path: Path,
    output_root: Path,
    index: int,
    nodes: list[tuple[float, float]],
    areas_um2: list[float],
    sentaurus_root: Path | None,
) -> dict[str, object]:
    count, scalars, vectors = read_vtk_point_data(vtk_path)
    if count != len(nodes):
        raise ValueError(f"{vtk_path}: VTK node count {count} != mesh node count {len(nodes)}")

    required_scalars = ("SRHRecombinationCm3PerS", "AugerRecombinationCm3PerS")
    required_vectors = (
        "J_n_drift",
        "J_n_diffusion",
        "J_n_total",
        "J_p_drift",
        "J_p_diffusion",
        "J_p_total",
        "SentaurusElectronCurrentDensityVector",
        "SentaurusHoleCurrentDensityVector",
        "DualFaceSgElectronCurrentDensityVector",
        "DualFaceSgHoleCurrentDensityVector",
    )
    missing = [name for name in required_scalars if name not in scalars]
    missing.extend(name for name in required_vectors if name not in vectors)
    if missing:
        raise ValueError(f"{vtk_path}: missing diagnostic fields: {missing}")

    token = f"vce_{index:03d}"
    source_rows: list[dict[str, object]] = []
    transport_rows: list[dict[str, object]] = []
    integrals = {"srh": 0.0, "auger": 0.0, "total": 0.0}
    max_rate_closure = 0.0
    max_current_closure = {"electron": 0.0, "hole": 0.0}
    peak = {
        "srh": {"absolute_rate_cm3_s": -1.0, "node_id": -1},
        "auger": {"absolute_rate_cm3_s": -1.0, "node_id": -1},
        "total": {"absolute_rate_cm3_s": -1.0, "node_id": -1},
    }

    for node_id, ((x_um, y_um), area_um2) in enumerate(zip(nodes, areas_um2, strict=True)):
        srh = scalars["SRHRecombinationCm3PerS"][node_id]
        auger = scalars["AugerRecombinationCm3PerS"][node_id]
        total = srh + auger
        area_cm2 = area_um2 * 1.0e-8
        depth_cm = 1.0e-4
        factor_A_per_um = Q_C * area_cm2 * depth_cm
        signed_currents = {
            "srh": srh * factor_A_per_um,
            "auger": auger * factor_A_per_um,
            "total": total * factor_A_per_um,
        }
        for name, value in signed_currents.items():
            integrals[name] += value
        for name, value in (("srh", srh), ("auger", auger), ("total", total)):
            if abs(value) > peak[name]["absolute_rate_cm3_s"]:
                peak[name] = {"absolute_rate_cm3_s": abs(value), "node_id": node_id}

        source_rows.append(
            {
                "node_id": node_id,
                "x_um": x_um,
                "y_um": y_um,
                "lumped_area_um2": area_um2,
                "srh_rate_cm3_s": srh,
                "auger_rate_cm3_s": auger,
                "total_recombination_rate_cm3_s": total,
                "srh_integrated_A_per_um": signed_currents["srh"],
                "auger_integrated_A_per_um": signed_currents["auger"],
                "total_integrated_A_per_um": signed_currents["total"],
                "electron_continuity_source_A_per_um": -signed_currents["total"],
                "hole_continuity_source_A_per_um": -signed_currents["total"],
            }
        )

        row: dict[str, object] = {"node_id": node_id, "x_um": x_um, "y_um": y_um}
        for carrier, prefix in (("electron", "J_n"), ("hole", "J_p")):
            drift = vectors[f"{prefix}_drift"][node_id]
            diffusion = vectors[f"{prefix}_diffusion"][node_id]
            total_vector = vectors[f"{prefix}_total"][node_id]
            max_current_closure[carrier] = max(
                max_current_closure[carrier],
                vector_closure_error(total_vector, drift, diffusion),
            )
            for component, axis in enumerate(("x", "y")):
                row[f"{carrier}_drift_{axis}_legacy_diagnostic"] = drift[component]
                row[f"{carrier}_diffusion_{axis}_legacy_diagnostic"] = diffusion[component]
                row[f"{carrier}_total_{axis}_legacy_diagnostic"] = total_vector[component]
            physical = vectors[
                "DualFaceSgElectronCurrentDensityVector"
                if carrier == "electron"
                else "DualFaceSgHoleCurrentDensityVector"
            ][node_id]
            row[f"{carrier}_physical_total_x_A_per_cm2"] = physical[0]
            row[f"{carrier}_physical_total_y_A_per_cm2"] = physical[1]
        transport_rows.append(row)

    source_path = output_root / "sources" / f"{token}.csv"
    transport_path = output_root / "transport" / f"{token}.csv"
    write_csv(source_path, source_rows)
    write_csv(transport_path, transport_rows)
    result: dict[str, object] = {
        "index": index,
        "VBE_V": 0.7,
        "VCE_V": index / 10.0,
        "node_count": count,
        "source_file": source_path.relative_to(output_root).as_posix(),
        "transport_file": transport_path.relative_to(output_root).as_posix(),
        "source_sha256": sha256(source_path),
        "transport_sha256": sha256(transport_path),
        "vtk_sha256": sha256(vtk_path),
        "signed_recombination_integral_A_per_um": integrals,
        "continuity_source_integral_A_per_um": {
            "electron": -integrals["total"],
            "hole": -integrals["total"],
        },
        "peak_locations": peak,
        "closure": {
            "maximum_total_rate_minus_components_cm3_s": max_rate_closure,
            "maximum_electron_current_vector_closure_legacy_diagnostic": max_current_closure["electron"],
            "maximum_hole_current_vector_closure_legacy_diagnostic": max_current_closure["hole"],
        },
    }
    sentaurus_point = sentaurus_root / token if sentaurus_root is not None else None
    if sentaurus_point is None or not (sentaurus_point / "field_manifest.json").is_file():
        result["sentaurus_comparison"] = {
            "status": "missing_reference",
            "reason": "No exported SDevice state exists at this representative bias.",
        }
        return result

    fields = sentaurus_point / "fields"
    comparison = {
        "status": "characterization_only" if index != 30 else "asserted_elsewhere",
        "sentaurus_field_manifest_sha256": sha256(sentaurus_point / "field_manifest.json"),
        "current_density": {
            "electron": vector_metrics(
                sentaurus_vector(fields / "eCurrentDensity_region0.csv"),
                vectors["DualFaceSgElectronCurrentDensityVector"],
                1.0e-6,
            ),
            "hole": vector_metrics(
                sentaurus_vector(fields / "hCurrentDensity_region0.csv"),
                vectors["DualFaceSgHoleCurrentDensityVector"],
                1.0e-6,
            ),
        },
        "recombination": {
            "srh": source_metrics(
                sentaurus_scalar(fields / "srhRecombination_region0.csv"),
                scalars["SRHRecombinationCm3PerS"],
                areas_um2,
                1.0e-6,
                nodes,
            ),
            "auger": source_metrics(
                sentaurus_scalar(fields / "AugerRecombination_region0.csv"),
                scalars["AugerRecombinationCm3PerS"],
                areas_um2,
                1.0e-6,
                nodes,
            ),
        },
    }
    result["sentaurus_comparison"] = comparison
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accepted-root", type=Path, required=True)
    parser.add_argument("--mesh-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--sentaurus-root",
        type=Path,
        help="Optional root containing exported vce_NNN SDevice states",
    )
    args = parser.parse_args()

    node_rows = read_csv(args.mesh_root / "nodes.csv")
    nodes = [(float(row["x_um"]), float(row["y_um"])) for row in node_rows]
    triangles = [
        (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        for row in read_csv(args.mesh_root / "elements.csv")
    ]
    areas = lumped_node_areas_um2(nodes, triangles)
    points = [
        export_point(
            args.accepted_root / "fields" / f"vce_{index:03d}.vtk",
            args.output_root,
            index,
            nodes,
            areas,
            args.sentaurus_root,
        )
        for index in REPRESENTATIVE_INDICES
    ]
    manifest = {
        "schema_version": 1,
        "device": "Genius NPN BJT",
        "state_contract": "derived only from unique fixed-bias accepted Vela states",
        "unit_contract": {
            "rate": "cm^-3 s^-1",
            "physical_current_density": "A/cm^2, dual-face-length-weighted nodal representation of production SG line flux",
            "drift_diffusion_decomposition": "legacy Vela diagnostic scale; algebraic closure only",
            "integrated_source": "A/um for 1 um out-of-plane depth",
            "continuity_sign": "positive recombination is a negative source in both carrier continuity residuals",
        },
        "mesh_sha256": sha256(args.mesh_root / "field_manifest.json"),
        "points": points,
        "sentaurus_comparison_coverage": {
            "available_VCE_V": [
                point["VCE_V"]
                for point in points
                if point["sentaurus_comparison"]["status"] != "missing_reference"
            ],
            "missing_VCE_V": [
                point["VCE_V"]
                for point in points
                if point["sentaurus_comparison"]["status"] == "missing_reference"
            ],
        },
        "overall_pass": all(
            point["node_count"] == len(nodes)
            and point["closure"]["maximum_total_rate_minus_components_cm3_s"] == 0.0
            and point["closure"]["maximum_electron_current_vector_closure_legacy_diagnostic"] <= 1.0e-9
            and point["closure"]["maximum_hole_current_vector_closure_legacy_diagnostic"] <= 1.0e-9
            for point in points
        ),
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, allow_nan=False))
    return 0 if manifest["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
