#!/usr/bin/env python3
"""Audit BJT hole SG edge flux against Vela and SDevice nodal current fields."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path

from compare_genius_bjt_transport_fields import read_vtk_point_data, sentaurus_vector


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"


def absolute(path: Path) -> str:
    return str(path.resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def projection_metrics(reference: list[float], actual: list[float]) -> dict[str, float | int]:
    if len(reference) != len(actual) or not reference:
        raise ValueError("projection metrics require equal non-empty arrays")
    peak = max(abs(value) for value in reference)
    selected = [abs(value) >= peak * 1.0e-3 for value in reference]
    pairs = [(r, a) for r, a, keep in zip(reference, actual, selected, strict=True) if keep]
    rr = sum(r * r for r, _ in pairs)
    aa = sum(a * a for _, a in pairs)
    dot = sum(r * a for r, a in pairs)
    ratios = [abs(a / r) for r, a in pairs if r != 0.0]
    return {
        "selected_edge_count": len(pairs),
        "reference_peak_A_per_cm2": peak,
        "normalized_rmse": math.sqrt(sum((a - r) ** 2 for r, a in pairs) / rr),
        "signed_cosine_similarity": dot / math.sqrt(rr * aa) if rr > 0.0 and aa > 0.0 else 0.0,
        "sign_flipped_normalized_rmse": math.sqrt(sum((-a - r) ** 2 for r, a in pairs) / rr),
        "least_squares_scale_actual_to_reference": dot / aa if aa > 0.0 else 0.0,
        "median_absolute_ratio": statistics.median(ratios),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=REPO / "build-release" / "vela_example_runner.exe")
    parser.add_argument("--output-root", type=Path, default=BUILD_ROOT / "m1_hole_current_edge_audit")
    parser.add_argument(
        "--sentaurus-root",
        type=Path,
        default=BUILD_ROOT / "m1_current_diagnosis" / "sentaurus_vce3",
    )
    args = parser.parse_args()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)

    config = json.loads(
        (FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json").read_text(encoding="utf-8")
    )
    config.pop("sweep", None)
    config.update(
        {
            "simulation_type": "sg_edge_flux_probe",
            "mesh_file": absolute(FIXTURE / "vela" / "input" / "mesh.json"),
            "node_doping_file": absolute(FIXTURE / "vela" / "input" / "doping.csv"),
            "materials_file": absolute(FIXTURE / "vela" / "materials_sentaurus2022.json"),
            "state_file": absolute(BUILD_ROOT / "m1_accepted_states" / "states" / "vce_030.csv"),
            "output_csv": absolute(output / "sg_edges.csv"),
        }
    )
    config_path = output / "sg_edge_probe.json"
    write_json(config_path, config)
    process = subprocess.run(
        [absolute(args.runner), "--config", absolute(config_path)],
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    (output / "runner.stdout.log").write_text(process.stdout, encoding="utf-8")
    (output / "runner.stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode:
        raise RuntimeError(process.stderr or process.stdout)

    sentaurus = sentaurus_vector(args.sentaurus_root / "fields" / "hCurrentDensity_region0.csv")
    _, _, vectors = read_vtk_point_data(BUILD_ROOT / "m1_accepted_states" / "fields" / "vce_030.vtk")
    vela = vectors["SentaurusHoleCurrentDensityVector"]
    doping = read_csv(FIXTURE / "vela" / "input" / "doping.csv")
    net_doping = [float(row["donors_cm3"]) - float(row["acceptors_cm3"]) for row in doping]

    rows: list[dict[str, object]] = []
    for edge in read_csv(output / "sg_edges.csv"):
        n0, n1 = int(edge["node0"]), int(edge["node1"])
        x0, y0 = float(edge["x0"]), float(edge["y0"])
        x1, y1 = float(edge["x1"]), float(edge["y1"])
        length = float(edge["length_m"])
        couple = float(edge["couple_m"])
        if length <= 0.0 or couple <= 0.0:
            continue
        ux, uy = (x1 - x0) / length, (y1 - y0) / length
        midpoint_x_um = 0.5 * (x0 + x1) * 1.0e6
        midpoint_y_um = 0.5 * (y0 + y1) * 1.0e6
        in_bc_window = (
            2.5 <= midpoint_x_um <= 4.5
            and 0.55 <= midpoint_y_um <= 0.85
            and abs(uy) >= 0.5
        )
        if not in_bc_window:
            continue
        sent_x = 0.5 * (sentaurus[n0][0] + sentaurus[n1][0])
        sent_y = 0.5 * (sentaurus[n0][1] + sentaurus[n1][1])
        vela_x = 0.5 * (vela[n0][0] + vela[n1][0])
        vela_y = 0.5 * (vela[n0][1] + vela[n1][1])
        # The 2-D SG diagnostic is particles/(m of out-of-plane depth*s).
        # Divide by the dual-face length and multiply by q to obtain A/m^2.
        sg_A_cm2 = (
            Q_C * float(edge["hole_particle_line_flux_per_m_s"]) / couple / 1.0e4
        )
        rows.append(
            {
                "edge_id": int(edge["edge_id"]),
                "node0": n0,
                "node1": n1,
                "midpoint_x_um": midpoint_x_um,
                "midpoint_y_um": midpoint_y_um,
                "unit_x": ux,
                "unit_y": uy,
                "crosses_net_doping_sign": net_doping[n0] * net_doping[n1] < 0.0,
                "net_doping0_cm3": net_doping[n0],
                "net_doping1_cm3": net_doping[n1],
                "sg_hole_projection_A_cm2": sg_A_cm2,
                "vela_node_projection_A_cm2": vela_x * ux + vela_y * uy,
                "sentaurus_node_projection_A_cm2": sent_x * ux + sent_y * uy,
            }
        )
    if not rows:
        raise RuntimeError("base-collector edge selection is empty")
    write_csv(output / "base_collector_hole_edges.csv", rows)

    sent = [float(row["sentaurus_node_projection_A_cm2"]) for row in rows]
    sg = [float(row["sg_hole_projection_A_cm2"]) for row in rows]
    vela_node = [float(row["vela_node_projection_A_cm2"]) for row in rows]
    crossing = [row for row in rows if bool(row["crosses_net_doping_sign"])]
    summary = {
        "schema_version": 1,
        "bias": {"VBE_V": 0.7, "VCE_V": 3.0},
        "selection": {
            "midpoint_x_um": [2.5, 4.5],
            "midpoint_y_um": [0.55, 0.85],
            "minimum_absolute_vertical_edge_direction": 0.5,
            "edge_count": len(rows),
            "net_doping_sign_crossing_edge_count": len(crossing),
        },
        "semantics": {
            "sentaurus": "endpoint-mean exported nodal hCurrentDensity projected on the Vela primal edge",
            "vela_node": "endpoint-mean SentaurusHoleCurrentDensityVector projected on the primal edge",
            "vela_sg": "production residual SG particle line flux divided by dual-face length and multiplied by q",
            "limitation": "SDevice internal directed-edge flux is unavailable; its nodal projection is not a native SG oracle.",
        },
        "comparisons": {
            "vela_sg_vs_sentaurus_node_projection": projection_metrics(sent, sg),
            "vela_node_reconstruction_vs_sentaurus_node_projection": projection_metrics(sent, vela_node),
            "vela_sg_vs_vela_node_reconstruction": projection_metrics(vela_node, sg),
        },
        "source_sha256": {
            "accepted_state": sha256(BUILD_ROOT / "m1_accepted_states" / "states" / "vce_030.csv"),
            "accepted_vtk": sha256(BUILD_ROOT / "m1_accepted_states" / "fields" / "vce_030.vtk"),
            "sentaurus_field_manifest": sha256(args.sentaurus_root / "field_manifest.json"),
            "sg_edges": sha256(output / "sg_edges.csv"),
        },
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
