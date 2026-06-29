#!/usr/bin/env python3
"""Audit Vela exported node current-density vectors and their code path."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ELECTRON_FIELD = "ElectronCurrentDensityVector"
HOLE_FIELD = "HoleCurrentDensityVector"
TOTAL_FIELD = "TotalCurrentDensityVector"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtk", type=Path, required=True, help="ASCII VTK file written by writeDDSolutionVTK")
    parser.add_argument("--out-md", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    return parser.parse_args()


def parse_vtk_points_and_vectors(path: Path) -> tuple[list[tuple[float, float, float]], dict[str, list[tuple[float, float, float]]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    points: list[tuple[float, float, float]] = []
    vectors: dict[str, list[tuple[float, float, float]]] = {}
    i = 0
    point_count = 0
    while i < len(lines):
        tokens = lines[i].strip().split()
        if len(tokens) >= 3 and tokens[0] == "POINTS":
            point_count = int(tokens[1])
            i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < 3 * point_count:
                values.extend(float(part) for part in lines[i].strip().split())
                i += 1
            points = [
                (values[j], values[j + 1], values[j + 2])
                for j in range(0, 3 * point_count, 3)
            ]
            continue
        if len(tokens) >= 3 and tokens[0] == "VECTORS":
            name = tokens[1]
            i += 1
            values = []
            while i < len(lines) and len(values) < 3 * point_count:
                stripped = lines[i].strip()
                if stripped:
                    values.extend(float(part) for part in stripped.split())
                i += 1
            vectors[name] = [
                (values[j], values[j + 1], values[j + 2])
                for j in range(0, 3 * point_count, 3)
            ]
            continue
        i += 1
    if not points:
        raise ValueError(f"{path} does not contain POINTS")
    missing = [name for name in (ELECTRON_FIELD, HOLE_FIELD, TOTAL_FIELD) if name not in vectors]
    if missing:
        raise ValueError(f"{path} is missing VECTORS fields: {', '.join(missing)}")
    for name in (ELECTRON_FIELD, HOLE_FIELD, TOTAL_FIELD):
        if len(vectors[name]) != len(points):
            raise ValueError(f"{name} vector count does not match POINTS count")
    return points, vectors


def magnitude(vector: tuple[float, float, float]) -> float:
    return math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def write_consistency_csv(
    path: Path,
    points: list[tuple[float, float, float]],
    vectors: dict[str, list[tuple[float, float, float]]],
) -> tuple[int, float, float, bool]:
    path.parent.mkdir(parents=True, exist_ok=True)
    max_abs_error = 0.0
    max_relative_error = 0.0
    all_ok = True
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "node_id",
            "x_um",
            "y_um",
            "electron_current_density_x_A_per_cm2",
            "electron_current_density_y_A_per_cm2",
            "electron_current_density_mag_A_per_cm2",
            "hole_current_density_x_A_per_cm2",
            "hole_current_density_y_A_per_cm2",
            "hole_current_density_mag_A_per_cm2",
            "total_current_density_x_A_per_cm2",
            "total_current_density_y_A_per_cm2",
            "total_current_density_mag_A_per_cm2",
            "reconstructed_total_current_density_x_A_per_cm2",
            "reconstructed_total_current_density_y_A_per_cm2",
            "reconstructed_total_current_density_mag_A_per_cm2",
            "total_minus_electron_plus_hole_x_A_per_cm2",
            "total_minus_electron_plus_hole_y_A_per_cm2",
            "absolute_error_mag_A_per_cm2",
            "relative_error_mag",
            "total_equals_electron_plus_hole",
        ])
        for node_id, point in enumerate(points):
            electron = vectors[ELECTRON_FIELD][node_id]
            hole = vectors[HOLE_FIELD][node_id]
            total = vectors[TOTAL_FIELD][node_id]
            reconstructed = (
                electron[0] + hole[0],
                electron[1] + hole[1],
                electron[2] + hole[2],
            )
            delta = (
                total[0] - reconstructed[0],
                total[1] - reconstructed[1],
                total[2] - reconstructed[2],
            )
            abs_error = magnitude(delta)
            rel_error = relative_error(magnitude(total), magnitude(reconstructed))
            ok = abs_error <= 1.0e-9 * max(magnitude(total), magnitude(reconstructed), 1.0)
            max_abs_error = max(max_abs_error, abs_error)
            max_relative_error = max(max_relative_error, rel_error)
            all_ok = all_ok and ok
            writer.writerow([
                node_id,
                f"{point[0] * 1.0e6:.17g}",
                f"{point[1] * 1.0e6:.17g}",
                f"{electron[0]:.17g}",
                f"{electron[1]:.17g}",
                f"{magnitude(electron):.17g}",
                f"{hole[0]:.17g}",
                f"{hole[1]:.17g}",
                f"{magnitude(hole):.17g}",
                f"{total[0]:.17g}",
                f"{total[1]:.17g}",
                f"{magnitude(total):.17g}",
                f"{reconstructed[0]:.17g}",
                f"{reconstructed[1]:.17g}",
                f"{magnitude(reconstructed):.17g}",
                f"{delta[0]:.17g}",
                f"{delta[1]:.17g}",
                f"{abs_error:.17g}",
                f"{rel_error:.17g}",
                "yes" if ok else "no",
            ])
    return len(points), max_abs_error, max_relative_error, all_ok


def write_summary(
    path: Path,
    vtk: Path,
    node_count: int,
    max_abs_error: float,
    max_relative_error: float,
    all_ok: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    total_answer = "yes" if all_ok else "no"
    field_rows = [
        ("electron_current_density_x/y/mag", "ElectronCurrentDensityVector", "node recovery / least-squares gradient post-processing"),
        ("hole_current_density_x/y/mag", "HoleCurrentDensityVector", "node recovery / least-squares gradient post-processing"),
        ("total_current_density_x/y/mag", "TotalCurrentDensityVector", "electron node vector plus hole node vector"),
    ]
    with path.open("w", encoding="utf-8", newline="\n") as output:
        output.write("# Current Density Export Path Audit\n\n")
        output.write(f"- vtk_source: {vtk}\n")
        output.write("- code_path: src/solver/GummelSolver.cpp::writeDDSolutionVTK\n")
        output.write("- data source: node recovery / least-squares gradient post-processing\n")
        output.write("- gradient source: computeNodeWeightedLeastSquaresGradients(psi, phin, phip)\n")
        output.write("- internal unit: A/cm^2\n")
        output.write("- output unit: A/cm^2\n")
        output.write("- conventional current density: yes\n")
        output.write("- electron particle flux: no\n")
        output.write("- multiplied by edge length/cell area/control volume/2D depth: no\n")
        output.write("- node_count: " + str(node_count) + "\n")
        output.write(f"- max_total_sum_absolute_error_A_per_cm2: {max_abs_error:.17g}\n")
        output.write(f"- max_total_sum_relative_error: {max_relative_error:.17g}\n")
        output.write(
            "- TotalCurrentDensity = ElectronCurrentDensity + HoleCurrentDensity: "
            + total_answer
            + "\n\n"
        )

        output.write("## Field Mapping\n\n")
        output.write("| requested field | exported VTK vector | data source | internal unit | output unit | conventional current | electron particle flux | geometry factor |\n")
        output.write("| --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for requested, exported, source in field_rows:
            output.write(
                f"| {requested} | {exported} | {source} | A/cm^2 | A/cm^2 | yes | no | none |\n"
            )

        output.write("\n## Unit Scale Checks\n\n")
        output.write("- A/m^2 -> A/cm^2: explicit /1e4 conversion in q * mobility * density * field scale.\n")
        output.write("- A/cm -> A/cm^2: not used by this export path.\n")
        output.write("- A/um -> A/cm^2: not used by this export path.\n")
        output.write("- edge length, cell area, control volume, and 2D depth are not applied in this node-vector export path.\n")
        output.write(
            "- exported node current density is a VTK post-processing recovery; it is not the SG edge source current used by avalanche assembly.\n"
        )


def main() -> int:
    args = parse_args()
    points, vectors = parse_vtk_points_and_vectors(args.vtk)
    node_count, max_abs_error, max_relative_error, all_ok = write_consistency_csv(
        args.out_csv, points, vectors
    )
    write_summary(args.out_md, args.vtk, node_count, max_abs_error, max_relative_error, all_ok)
    print(args.out_csv)
    print(args.out_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
