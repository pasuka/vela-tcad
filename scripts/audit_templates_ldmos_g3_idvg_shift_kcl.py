#!/usr/bin/env python3
"""Audit Templates/LDMOS G3 subthreshold shift and low-current KCL.

The audit keeps IALMob and continuation prediction disabled.  It compares
exact-node Sentaurus/Vela states at the seven shared 0--1 V gate biases and
replays the first nonzero Vela point through the production carrier operator.
Large generated artifacts belong in ignored ``reference_staging``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
BIAS_POINTS_V = tuple(index / 6.0 for index in range(7))
Q_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(value for value in values if math.isfinite(value))
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def distribution(values: Iterable[float]) -> dict[str, float | int]:
    samples = [value for value in values if math.isfinite(value)]
    if not samples:
        return {"count": 0, "median": math.nan, "p95": math.nan,
                "maximum_abs": math.nan}
    return {
        "count": len(samples),
        "median": statistics.median(samples),
        "p95": percentile([abs(value) for value in samples], 0.95),
        "maximum_abs": max(abs(value) for value in samples),
    }


def scalar_field(export_dir: Path, name: str, region: int = 0) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    }


def vector_field(
    export_dir: Path, name: str, region: int = 0,
) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    }


def state_rows(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            key: float(row[key])
            for key in ("psi", "phin", "phip", "electrons_m3", "holes_m3")
        }
        for row in read_csv(path)
    }


def silicon_interface_nodes(mesh_path: Path) -> set[int]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    by_region: dict[int, set[int]] = {}
    for triangle in mesh["triangles"]:
        by_region.setdefault(int(triangle["region_id"]), set()).update(
            int(node) for node in triangle["node_ids"]
        )
    silicon = next(
        int(region["id"]) for region in mesh["regions"]
        if region["material"].lower() == "si"
    )
    oxide_nodes: set[int] = set()
    for region in mesh["regions"]:
        if region["material"].lower() == "sio2":
            oxide_nodes.update(by_region[int(region["id"])])
    return by_region[silicon] & oxide_nodes


def vela_state_name(bias: float) -> str:
    return f"g3_idvg_point_bias_{bias:.6f}".replace(".", "p") + ".csv"


def sentaurus_state(export_dir: Path) -> dict[str, dict[int, float]]:
    return {
        "psi": scalar_field(export_dir, "ElectrostaticPotential"),
        "phin": scalar_field(export_dir, "eQuasiFermiPotential"),
        "phip": scalar_field(export_dir, "hQuasiFermiPotential"),
        "electrons_m3": {
            node: value * 1.0e6
            for node, value in scalar_field(export_dir, "eDensity").items()
        },
        "holes_m3": {
            node: value * 1.0e6
            for node, value in scalar_field(export_dir, "hDensity").items()
        },
    }


def write_sentaurus_state_csv(
    export_dir: Path, mesh_path: Path, output: Path,
) -> Path:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    fields: dict[str, dict[int, float]] = {}
    for key, name in (
        ("psi", "ElectrostaticPotential"),
        ("phin", "eQuasiFermiPotential"),
        ("phip", "hQuasiFermiPotential"),
    ):
        merged: dict[int, float] = {}
        for region in range(3):
            path = export_dir / "fields" / f"{name}_region{region}.csv"
            if path.is_file():
                for node, value in scalar_field(export_dir, name, region).items():
                    merged.setdefault(node, value)
        fields[key] = merged
    fields["electrons_m3"] = {
        node: value * 1.0e6
        for node, value in scalar_field(export_dir, "eDensity").items()
    }
    fields["holes_m3"] = {
        node: value * 1.0e6
        for node, value in scalar_field(export_dir, "hDensity").items()
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
        )
        for node in sorted(int(item["id"]) for item in mesh["nodes"]):
            writer.writerow([
                node,
                format(fields["psi"][node], ".17g"),
                format(fields["phin"][node], ".17g"),
                format(fields["phip"][node], ".17g"),
                format(fields["electrons_m3"].get(node, 0.0), ".17g"),
                format(fields["holes_m3"].get(node, 0.0), ".17g"),
            ])
    return output


def compare_nodes(
    sentaurus: dict[str, dict[int, float]],
    vela: dict[int, dict[str, float]],
    nodes: set[int],
) -> dict[str, Any]:
    common = sorted(nodes & set(vela) & set(sentaurus["psi"]))
    result: dict[str, Any] = {"nodes": len(common)}
    for key in ("psi", "phin", "phip"):
        result[f"vela_minus_sentaurus_{key}_V"] = distribution(
            vela[node][key] - sentaurus[key][node] for node in common
        )
    for key in ("electrons_m3", "holes_m3"):
        result[f"vela_minus_sentaurus_{key}_dex"] = distribution(
            math.log10(max(vela[node][key], 1.0))
            - math.log10(max(sentaurus[key][node], 1.0))
            for node in common
        )
    return result


def terminal_precision(rows: list[dict[str, str]]) -> dict[str, Any]:
    methods = {
        "native": "current_total_A_per_um",
        "compensated": "current_total_compensated_A_per_um",
        "long_double": "current_total_long_double_reference_A_per_um",
    }
    result: dict[str, Any] = {}
    for method, column in methods.items():
        values = [float(row[column]) for row in rows]
        result[method] = {
            "contact_sum_A_per_um": sum(values),
            "sum_abs_contacts_A_per_um": sum(abs(value) for value in values),
        }
    result["native_vs_long_double_sum_difference_A_per_um"] = abs(
        result["native"]["contact_sum_A_per_um"]
        - result["long_double"]["contact_sum_A_per_um"]
    )
    return result


def runner_environment() -> dict[str, str]:
    environment = os.environ.copy()
    if os.name == "nt":
        environment["PATH"] = os.pathsep.join([
            r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin",
            environment.get("PATH", ""),
        ])
    environment["VELA_LINEAR_SOLVER"] = "sparselu"
    return environment


def run_probe(
    runner: Path, baseline: dict[str, Any], state: Path,
    bias: float, output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    config = deepcopy(baseline)
    config["simulation_type"] = "newton_carrier_term_probe"
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str((output_dir / "carrier_terms.csv").resolve())
    config.pop("sweep", None)
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    config_path = output_dir / "carrier_probe.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output_dir / "carrier_probe.log").resolve())],
        cwd=REPO, text=True, capture_output=True,
        env=runner_environment(), check=False,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)
    rows = read_csv(Path(config["output_csv"]))
    return {
        "return_code": completed.returncode,
        "electron_residual_signed_sum": sum(
            float(row["electron_residual"]) for row in rows
        ),
        "electron_residual_abs_sum": sum(
            abs(float(row["electron_residual"])) for row in rows
        ),
        "hole_residual_signed_sum": sum(
            float(row["hole_residual"]) for row in rows
        ),
        "hole_residual_abs_sum": sum(
            abs(float(row["hole_residual"])) for row in rows
        ),
        "artifacts": {
            "config": str(config_path.resolve()),
            "carrier_terms": config["output_csv"],
        },
    }


def drain_nodes(mesh_path: Path) -> set[int]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    return next(
        {int(node) for node in contact["node_ids"]}
        for contact in mesh["contacts"]
        if contact["name"].lower() == "drain"
    )


def drain_sg_cut(rows: list[dict[str, str]], nodes: set[int]) -> dict[str, float]:
    electron_particles = 0.0
    hole_particles = 0.0
    for row in rows:
        node0, node1 = int(row["node0"]), int(row["node1"])
        at0, at1 = node0 in nodes, node1 in nodes
        if at0 == at1:
            continue
        outward = 1.0 if at0 else -1.0
        electron_particles += outward * float(
            row["electron_particle_line_flux_per_m_s"]
        )
        hole_particles += outward * float(row["hole_particle_line_flux_per_m_s"])
    electron = -Q_C * electron_particles / 1.0e6
    hole = -Q_C * hole_particles / 1.0e6
    edge_currents: list[float] = []
    for row in rows:
        node0, node1 = int(row["node0"]), int(row["node1"])
        at0, at1 = node0 in nodes, node1 in nodes
        if at0 == at1:
            continue
        outward = 1.0 if at0 else -1.0
        edge_currents.append(
            -Q_C * outward * float(row["electron_particle_line_flux_per_m_s"])
            / 1.0e6
        )
    return {
        "electron_A_per_um": electron,
        "hole_signed_particle_convention_A_per_um": hole,
        "total_A_per_um": electron - hole,
        "electron_sum_abs_edges_A_per_um": sum(abs(value) for value in edge_currents),
        "electron_edge_cancellation_condition": (
            sum(abs(value) for value in edge_currents) / max(abs(electron), 1.0e-300)
        ),
    }


def sentaurus_drain_vector_current(
    export_dir: Path, mesh_path: Path,
) -> dict[str, float | int]:
    """Trapezoid-integrate native Sentaurus Jn along oriented drain edges."""
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    drain = next(
        contact for contact in mesh["contacts"]
        if contact["name"].lower() == "drain"
    )
    current_density = vector_field(export_dir, "eCurrentDensity")
    current = 0.0
    used = 0
    for node0, node1 in drain["edge_node_ids"]:
        if node0 not in current_density or node1 not in current_density:
            continue
        x0, y0 = coordinates[node0]
        x1, y1 = coordinates[node1]
        jx0, jy0 = current_density[node0]
        jx1, jy1 = current_density[node1]
        # Dot J [A/cm2] with the oriented line normal (dy,-dx) [um], then
        # convert A/cm of device depth to A/um of device depth (1e-8).
        current += (
            0.5 * (jx0 + jx1) * (y1 - y0)
            - 0.5 * (jy0 + jy1) * (x1 - x0)
        ) * 1.0e-8
        used += 1
    return {"oriented_electron_A_per_um": current, "edges": used}


def cell_qf_gradient_drain_current(
    export_dir: Path, mesh_path: Path,
) -> dict[str, float | int]:
    """Reconstruct Jn=q*mu*n*grad(phin) in drain-adjacent Tri3 cells."""
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    drain = next(
        contact for contact in mesh["contacts"]
        if contact["name"].lower() == "drain"
    )
    silicon = next(
        int(region["id"]) for region in mesh["regions"]
        if region["material"].lower() == "si"
    )
    edge_cells: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for triangle in mesh["triangles"]:
        ids = [int(node) for node in triangle["node_ids"]]
        for index in range(3):
            key = tuple(sorted((ids[index], ids[(index + 1) % 3])))
            edge_cells.setdefault(key, []).append(triangle)
    phin = scalar_field(export_dir, "eQuasiFermiPotential")
    density_cm3 = {
        node: value for node, value in scalar_field(export_dir, "eDensity").items()
    }
    mobility_cm2 = scalar_field(export_dir, "eMobility")
    current = 0.0
    used = 0
    for node0, node1 in drain["edge_node_ids"]:
        adjacent = [
            triangle
            for triangle in edge_cells.get(tuple(sorted((node0, node1))), [])
            if int(triangle["region_id"]) == silicon
        ]
        if len(adjacent) != 1:
            continue
        ids = [int(node) for node in adjacent[0]["node_ids"]]
        x0, y0 = coordinates[ids[0]]
        x1, y1 = coordinates[ids[1]]
        x2, y2 = coordinates[ids[2]]
        determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if abs(determinant) <= 1.0e-300:
            continue
        values = [phin[node] for node in ids]
        grad_x_per_um = (
            values[0] * (y1 - y2)
            + values[1] * (y2 - y0)
            + values[2] * (y0 - y1)
        ) / determinant
        grad_y_per_um = (
            values[0] * (x2 - x1)
            + values[1] * (x0 - x2)
            + values[2] * (x1 - x0)
        ) / determinant
        density = sum(density_cm3[node] for node in ids) / 3.0
        mobility = sum(mobility_cm2[node] for node in ids) / 3.0
        jx_A_cm2 = Q_C * mobility * density * grad_x_per_um * 1.0e4
        jy_A_cm2 = Q_C * mobility * density * grad_y_per_um * 1.0e4
        ex0, ey0 = coordinates[node0]
        ex1, ey1 = coordinates[node1]
        current += (
            jx_A_cm2 * (ey1 - ey0) - jy_A_cm2 * (ex1 - ex0)
        ) * 1.0e-8
        used += 1
    return {"oriented_electron_A_per_um": current, "edges": used}


def drain_density_reconstruction(
    rows: list[dict[str, str]], state: dict[int, dict[str, float]],
    nodes: set[int],
) -> dict[str, Any]:
    """Compare supplied and SG-reconstructed densities on drain-cut endpoints."""
    contact_errors: list[float] = []
    interior_errors: list[float] = []
    qf_drops: list[float] = []
    crossing_edges = 0
    for row in rows:
        node0, node1 = int(row["node0"]), int(row["node1"])
        at0, at1 = node0 in nodes, node1 in nodes
        if at0 == at1:
            continue
        crossing_edges += 1
        qf_drops.append(abs(float(row["phin1_V"]) - float(row["phin0_V"])))
        for endpoint, node, is_contact in ((0, node0, at0), (1, node1, at1)):
            supplied = max(state[node]["electrons_m3"], 1.0)
            reconstructed = max(
                float(row[f"electron_density{endpoint}_m3"]), 1.0
            )
            error = math.log10(reconstructed) - math.log10(supplied)
            (contact_errors if is_contact else interior_errors).append(error)
    return {
        "crossing_edges": crossing_edges,
        "contact_reconstructed_minus_supplied_electron_density_dex": distribution(
            contact_errors
        ),
        "interior_reconstructed_minus_supplied_electron_density_dex": distribution(
            interior_errors
        ),
        "electron_qf_drop_V": distribution(qf_drops),
    }


def dominant_drain_edges(
    rows: list[dict[str, str]], nodes: set[int], limit: int = 10,
) -> list[dict[str, float | int]]:
    records = []
    for row in rows:
        node0, node1 = int(row["node0"]), int(row["node1"])
        at0, at1 = node0 in nodes, node1 in nodes
        if at0 == at1:
            continue
        outward = 1.0 if at0 else -1.0
        current = (
            -Q_C * outward * float(row["electron_particle_line_flux_per_m_s"])
            / 1.0e6
        )
        records.append({
            "edge_id": int(row["edge_id"]),
            "electron_A_per_um": current,
            "electron_qf_drop_V": abs(
                float(row["phin1_V"]) - float(row["phin0_V"])
            ),
            "electron_mobility_m2_V_s": float(row["electron_mobility_m2_V_s"]),
            "couple_over_length": float(row["couple_m"]) / max(
                float(row["length_m"]), 1.0e-300
            ),
        })
    return sorted(records, key=lambda item: abs(item["electron_A_per_um"]), reverse=True)[:limit]


def run_sg_probe(
    runner: Path, baseline: dict[str, Any], state: Path,
    bias: float, output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    config = deepcopy(baseline)
    config["simulation_type"] = "sg_edge_flux_probe"
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str((output_dir / "sg_edges.csv").resolve())
    config.pop("sweep", None)
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    config_path = output_dir / "sg_probe.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output_dir / "sg_probe.log").resolve())],
        cwd=REPO, text=True, capture_output=True,
        env=runner_environment(), check=False,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)
    rows = read_csv(Path(config["output_csv"]))
    nodes = drain_nodes(Path(baseline["mesh_file"]))
    return {
        "return_code": completed.returncode,
        "drain_cut": drain_sg_cut(rows, nodes),
        "drain_density_reconstruction": drain_density_reconstruction(
            rows, state_rows(state), nodes
        ),
        "dominant_drain_edges": dominant_drain_edges(rows, nodes),
        "artifacts": {
            "config": str(config_path.resolve()),
            "sg_edges": config["output_csv"],
        },
    }


def vela_state_sg_replay(
    runner: Path, baseline: dict[str, Any], state_root: Path,
    candidate_curve: Path, output_dir: Path,
) -> dict[str, Any]:
    """Close the drain-cut integration against Vela's own terminal current."""
    candidate = {
        round(float(row["bias_V"]), 12): float(row["current_total_A_per_um"])
        for row in read_csv(candidate_curve)
    }
    results = []
    for index in (1, 3, 5):
        bias = BIAS_POINTS_V[index]
        state = state_root / vela_state_name(bias)
        case_dir = output_dir / f"vg_{bias:.6f}".replace(".", "p")
        probe = run_sg_probe(runner, baseline, state, bias, case_dir)
        cut_current = float(probe["drain_cut"]["total_A_per_um"])
        terminal_current = candidate[round(bias, 12)]
        results.append({
            "bias_V": bias,
            "vela_terminal_A_per_um": terminal_current,
            "vela_sg_drain_cut_A_per_um": cut_current,
            "relative_error": abs(cut_current - terminal_current) / max(
                abs(terminal_current), 1.0e-300
            ),
            "probe": probe,
        })
    return {
        "points": results,
        "maximum_relative_error": max(row["relative_error"] for row in results),
    }


def strict_kcl_reclose(
    runner: Path, baseline: dict[str, Any], state: Path,
    bias: float, output_dir: Path,
) -> dict[str, Any]:
    """Reclose one low-current point with a tighter electron block ceiling."""
    output_dir.mkdir(parents=True, exist_ok=True)
    config = deepcopy(baseline)
    config["output_csv"] = str((output_dir / "curve.csv").resolve())
    config["solver"]["block_absolute_convergence"] = {
        "mode": "enforce",
        "psi_residual_ceiling": 5.0e-8,
        "electron_residual_ceiling": 1.0e-11,
        "hole_residual_ceiling": 3.0e-10,
    }
    contacts = [contact["name"] for contact in config["contacts"]]
    config["sweep"] = {
        "mode": "iv",
        "contact": "gate",
        "current_contact": "drain",
        "start": bias,
        "stop": bias,
        "step": 1.0,
        "bias_points": [bias],
        "initial_state_file": str(state.resolve()),
        "write_state_file": str((output_dir / "state.csv").resolve()),
        "write_vtk": False,
        "diagnostics": {
            "terminal_balance": {
                "enabled": True,
                "contacts": contacts,
                "csv_file": str((output_dir / "terminal_balance.csv").resolve()),
            },
            "srh_balance": {
                "enabled": True,
                "material": "Si",
                "drain_contact": "drain",
                "substrate_contact": "substrate",
                "kcl_contacts": contacts,
                "resolution_margin_ratio": 10.0,
                "csv_file": str((output_dir / "srh_balance.csv").resolve()),
            },
            "newton_history": {
                "enabled": True,
                "csv_file": str((output_dir / "newton_history.csv").resolve()),
                "attempts_csv_file": str((output_dir / "newton_attempts.csv").resolve()),
                "iterations_csv_file": str((output_dir / "newton_iterations.csv").resolve()),
                "rejected_state_directory": str(
                    (output_dir / "rejected_states").resolve()
                ),
            },
        },
    }
    config["_comment"] = (
        "Single-factor low-current KCL reclose: G3 physics unchanged; "
        "electron block ceiling tightened from 2e-9 to 1e-11."
    )
    config_path = output_dir / "strict_kcl_reclose.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output_dir / "strict_kcl_reclose.log").resolve())],
        cwd=REPO, text=True, capture_output=True,
        env=runner_environment(), check=False,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    result: dict[str, Any] = {
        "return_code": completed.returncode,
        "electron_residual_ceiling": 1.0e-11,
        "artifacts": {"config": str(config_path.resolve())},
    }
    curve = output_dir / "curve.csv"
    if curve.is_file() and read_csv(curve):
        result["curve"] = read_csv(curve)[-1]
    srh = output_dir / "srh_balance.csv"
    if srh.is_file() and read_csv(srh):
        row = read_csv(srh)[-1]
        result["kcl"] = {
            "drain_current_A_per_um": float(row["drain_total_current_A_per_um"]),
            "four_terminal_kcl_residual_A_per_um": float(
                row["four_terminal_kcl_residual_A_per_um"]
            ),
            "id_to_kcl_ratio": float(row["id_to_kcl_residual_ratio"]),
        }
    return result


def sentaurus_state_sg_replay(
    runner: Path, baseline: dict[str, Any], export_root: Path,
    reference_curve: Path, output_dir: Path,
) -> dict[str, Any]:
    reference = {
        round(float(row["bias_V"]), 12): float(row["current_total_A_per_um"])
        for row in read_csv(reference_curve)
    }
    results = []
    for index in (1, 3, 5):
        bias = BIAS_POINTS_V[index]
        case_dir = output_dir / f"vg_{bias:.6f}".replace(".", "p")
        state = write_sentaurus_state_csv(
            export_root / f"export_{index:04d}",
            Path(baseline["mesh_file"]), case_dir / "sentaurus_state.csv",
        )
        probe = run_sg_probe(runner, baseline, state, bias, case_dir)
        current = float(probe["drain_cut"]["total_A_per_um"])
        ref = reference[round(bias, 12)]
        results.append({
            "bias_V": bias,
            "sentaurus_terminal_A_per_um": ref,
            "vela_sg_on_sentaurus_state_A_per_um": current,
            "signed_ratio": current / ref,
            "magnitude_error_dex": abs(
                math.log10(max(abs(current), 1.0e-300))
                - math.log10(max(abs(ref), 1.0e-300))
            ),
            "sentaurus_native_drain_vector_integration": sentaurus_drain_vector_current(
                export_root / f"export_{index:04d}", Path(baseline["mesh_file"])
            ),
            "cell_qf_gradient_drain_reconstruction": cell_qf_gradient_drain_current(
                export_root / f"export_{index:04d}", Path(baseline["mesh_file"])
            ),
            "probe": probe,
        })
    return {
        "points": results,
        "median_magnitude_error_dex": statistics.median(
            row["magnitude_error_dex"] for row in results
        ),
    }


def run_high_field_off_control(
    runner: Path, baseline: dict[str, Any], initial_state: Path,
    reference_curve: Path, output_dir: Path,
) -> dict[str, Any]:
    """Run the exact seven-point G4 control with only HFS removed."""
    output_dir.mkdir(parents=True, exist_ok=True)
    config = deepcopy(baseline)
    config["_comment"] = (
        "Single-factor G4 control: G3 exact mesh/physics with high-field "
        "saturation disabled; IALMob and predictor remain disabled."
    )
    config["solver"]["mobility"] = {"model": "constant"}
    config["output_csv"] = str((output_dir / "g4_idvg.csv").resolve())
    sweep = config["sweep"]
    sweep.pop("predictor", None)
    sweep.pop("continuation", None)
    sweep["bias_points"] = list(BIAS_POINTS_V)
    sweep["start"] = BIAS_POINTS_V[0]
    sweep["stop"] = BIAS_POINTS_V[-1]
    sweep["step"] = BIAS_POINTS_V[-1]
    sweep["initial_state_file"] = str(initial_state.resolve())
    sweep["write_state_file"] = str((output_dir / "g4_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str(
        (output_dir / "g4_point").resolve()
    )
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": [contact["name"] for contact in config["contacts"]],
        "csv_file": str((output_dir / "terminal_balance.csv").resolve()),
    }
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": [contact["name"] for contact in config["contacts"]],
        "resolution_margin_ratio": 10.0,
        "csv_file": str((output_dir / "srh_balance.csv").resolve()),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((output_dir / "newton_history.csv").resolve()),
        "attempts_csv_file": str((output_dir / "newton_attempts.csv").resolve()),
        "iterations_csv_file": str((output_dir / "newton_iterations.csv").resolve()),
        "rejected_state_directory": str((output_dir / "rejected_states").resolve()),
    }
    config_path = output_dir / "g4_high_field_off.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output_dir / "g4_high_field_off.log").resolve())],
        cwd=REPO, text=True, capture_output=True,
        env=runner_environment(), check=False,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    candidate_rows = read_csv(Path(config["output_csv"]))
    reference = {
        round(float(row["bias_V"]), 12): float(row["current_total_A_per_um"])
        for row in read_csv(reference_curve)
    }
    errors = []
    point_metrics = []
    for row in candidate_rows:
        bias = float(row["bias_V"])
        key = round(bias, 12)
        if key not in reference:
            continue
        candidate = float(row["current_total_A_per_um"])
        ref = reference[key]
        error = abs(math.log10(max(abs(candidate), 1.0e-300))
                    - math.log10(max(abs(ref), 1.0e-300)))
        errors.append(error)
        point_metrics.append({
            "bias_V": bias,
            "reference_A_per_um": ref,
            "candidate_A_per_um": candidate,
            "abs_log_error_dex": error,
        })
    return {
        "return_code": completed.returncode,
        "points": len(candidate_rows),
        "median_abs_log_error_dex": (
            statistics.median(errors) if errors else math.nan
        ),
        "p95_abs_log_error_dex": percentile(errors, 0.95),
        "maximum_abs_log_error_dex": max(errors, default=math.nan),
        "point_metrics": point_metrics,
        "artifacts": {
            "config": str(config_path.resolve()),
            "curve": config["output_csv"],
            "state": sweep["write_state_file"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--sentaurus-export-root", type=Path, required=True)
    parser.add_argument("--vela-state-root", type=Path, required=True)
    parser.add_argument("--terminal-balance", type=Path, required=True)
    parser.add_argument("--srh-balance", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-high-field-off-control", action="store_true")
    parser.add_argument("--g4-reference-curve", type=Path)
    parser.add_argument("--run-sentaurus-state-sg-replay", action="store_true")
    parser.add_argument("--g3-reference-curve", type=Path)
    parser.add_argument("--run-vela-state-sg-replay", action="store_true")
    parser.add_argument("--g3-candidate-curve", type=Path)
    parser.add_argument("--run-strict-kcl-reclose", action="store_true")
    parser.add_argument(
        "--run-no-negative-cotangent-fallback-control", action="store_true"
    )
    args = parser.parse_args()

    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    if "predictor" in baseline.get("sweep", {}):
        raise ValueError("G3 shift audit requires predictor disabled")
    mobility = baseline["solver"]["mobility"]
    if "ialmob" in json.dumps(mobility).lower():
        raise ValueError("G3 shift audit requires IALMob disabled")

    interface = silicon_interface_nodes(Path(baseline["mesh_file"]))
    sentaurus_states = [
        sentaurus_state(args.sentaurus_export_root / f"export_{index:04d}")
        for index in range(len(BIAS_POINTS_V))
    ]
    vela_states = [
        state_rows(args.vela_state_root / vela_state_name(bias))
        for bias in BIAS_POINTS_V
    ]
    endpoint_log_gain = {
        node: math.log10(max(sentaurus_states[-1]["electrons_m3"][node], 1.0))
        - math.log10(max(sentaurus_states[0]["electrons_m3"][node], 1.0))
        for node in interface & set(sentaurus_states[0]["electrons_m3"])
    }
    responsive_interface = {
        node for node, gain in endpoint_log_gain.items() if gain >= 2.0
    }

    state_comparison = []
    equivalent_gate_shifts: list[float] = []
    for index, bias in enumerate(BIAS_POINTS_V):
        sentaurus = sentaurus_states[index]
        vela = vela_states[index]
        state_comparison.append({
            "bias_V": bias,
            "all_silicon": compare_nodes(
                sentaurus, vela, set(sentaurus["psi"])
            ),
            "responsive_interface": compare_nodes(
                sentaurus, vela, responsive_interface
            ),
        })
        if 0 < index < len(BIAS_POINTS_V) - 1:
            for node in responsive_interface:
                lower = math.log10(max(
                    sentaurus_states[index - 1]["electrons_m3"][node], 1.0))
                upper = math.log10(max(
                    sentaurus_states[index + 1]["electrons_m3"][node], 1.0))
                slope = (upper - lower) / (
                    BIAS_POINTS_V[index + 1] - BIAS_POINTS_V[index - 1])
                if abs(slope) < 1.0e-12:
                    continue
                log_error = (
                    math.log10(max(vela[node]["electrons_m3"], 1.0))
                    - math.log10(max(sentaurus["electrons_m3"][node], 1.0))
                )
                equivalent_gate_shifts.append(log_error / slope)

    terminal_rows = read_csv(args.terminal_balance)
    low_bias = BIAS_POINTS_V[1]
    low_terminal = [
        row for row in terminal_rows
        if abs(float(row["bias_V"]) - low_bias) <= 1.0e-12
    ]
    srh_row = next(
        row for row in read_csv(args.srh_balance)
        if abs(float(row["bias_V"]) - low_bias) <= 1.0e-12
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    probe = run_probe(
        args.runner, baseline, args.vela_state_root / vela_state_name(low_bias),
        low_bias, args.output_dir / "vg_0p166667_carrier_probe",
    )
    kcl = float(srh_row["four_terminal_kcl_residual_A_per_um"])
    electron_sum = float(probe["electron_residual_signed_sum"])
    report = {
        "schema": "vela.templates_ldmos.g3_idvg_shift_kcl_audit.v1",
        "status": "complete",
        "contracts": {
            "ialmob": "disabled",
            "predictor": "disabled",
            "bias_alignment": "seven exact shared points from 0 to 1 V",
        },
        "interface_nodes": len(interface),
        "gate_responsive_interface_nodes": len(responsive_interface),
        "state_comparison": state_comparison,
        "equivalent_gate_shift_from_interface_density_V": distribution(
            equivalent_gate_shifts
        ),
        "low_current_kcl": {
            "bias_V": low_bias,
            "drain_current_A_per_um": float(srh_row["drain_total_current_A_per_um"]),
            "four_terminal_kcl_residual_A_per_um": kcl,
            "id_to_kcl_ratio": float(srh_row["id_to_kcl_residual_ratio"]),
            "terminal_precision": terminal_precision(low_terminal),
            "carrier_probe": probe,
            "electron_residual_sum_to_kcl_scale": (
                kcl / electron_sum if electron_sum != 0.0 else math.nan
            ),
        },
    }
    if args.run_high_field_off_control:
        if args.g4_reference_curve is None:
            raise ValueError(
                "--g4-reference-curve is required with "
                "--run-high-field-off-control"
            )
        initial_state = Path(baseline["sweep"]["initial_state_file"])
        report["high_field_off_control"] = run_high_field_off_control(
            args.runner, baseline, initial_state, args.g4_reference_curve,
            args.output_dir / "g4_high_field_off_control",
        )
    if args.run_sentaurus_state_sg_replay:
        if args.g3_reference_curve is None:
            raise ValueError(
                "--g3-reference-curve is required with "
                "--run-sentaurus-state-sg-replay"
            )
        report["sentaurus_state_sg_replay"] = sentaurus_state_sg_replay(
            args.runner, baseline, args.sentaurus_export_root,
            args.g3_reference_curve,
            args.output_dir / "sentaurus_state_sg_replay",
        )
    if args.run_vela_state_sg_replay:
        if args.g3_candidate_curve is None:
            raise ValueError(
                "--g3-candidate-curve is required with "
                "--run-vela-state-sg-replay"
            )
        report["vela_state_sg_replay"] = vela_state_sg_replay(
            args.runner, baseline, args.vela_state_root,
            args.g3_candidate_curve,
            args.output_dir / "vela_state_sg_replay",
        )
    if args.run_strict_kcl_reclose:
        report["strict_kcl_reclose"] = strict_kcl_reclose(
            args.runner, baseline,
            args.vela_state_root / vela_state_name(low_bias), low_bias,
            args.output_dir / "strict_kcl_reclose",
        )
    if args.run_no_negative_cotangent_fallback_control:
        if args.g3_reference_curve is None:
            raise ValueError(
                "--g3-reference-curve is required with "
                "--run-no-negative-cotangent-fallback-control"
            )
        control = deepcopy(baseline)
        control.setdefault("mesh_geometry", {})[
            "fallback_negative_cotangent"
        ] = False
        report["no_negative_cotangent_fallback_control"] = {
            "single_factor": (
                "negative local cotangent contributions use zero instead of "
                "positive barycentric fallback"
            ),
            "replay": sentaurus_state_sg_replay(
                args.runner, control, args.sentaurus_export_root,
                args.g3_reference_curve,
                args.output_dir / "no_negative_cotangent_fallback_control",
            ),
        }
    output = args.output_dir / "summary.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": str(output.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
