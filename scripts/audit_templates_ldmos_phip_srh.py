#!/usr/bin/env python3
"""Classify the Templates/LDMOS hole-QF offset and its SRH sensitivity.

The audit separates an absolute phip offset from centered spatial residuals
and P1 cell gradients.  It then shifts only that measured uniform component,
uses the production SG probe to reconstruct mapping-consistent hole density,
and replays both the original and shifted states without a nonlinear solve.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import (
        BIAS_POINTS_V,
        read_csv,
        run_sg_probe,
        runner_environment,
        scalar_field,
        state_rows,
    )
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_idvg_shift_kcl import (
        BIAS_POINTS_V,
        read_csv,
        run_sg_probe,
        runner_environment,
        scalar_field,
        state_rows,
    )


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(value for value in values if math.isfinite(value))
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def stats(values: Iterable[float]) -> dict[str, float | int]:
    samples = [value for value in values if math.isfinite(value)]
    return {
        "count": len(samples),
        "median": statistics.median(samples) if samples else math.nan,
        "p95_abs": percentile([abs(value) for value in samples], 0.95),
        "maximum_abs": max((abs(value) for value in samples), default=math.nan),
    }


def silicon_nodes_and_cells(mesh: dict) -> tuple[set[int], list[dict]]:
    regions = {
        int(region["id"]) for region in mesh["regions"]
        if region["material"].lower() == "si"
    }
    cells = [
        cell for cell in mesh["triangles"]
        if int(cell["region_id"]) in regions
    ]
    nodes = {
        int(node) for cell in cells for node in cell["node_ids"]
    }
    return nodes, cells


def p1_gradient(
    cell: dict, coordinates: dict[int, tuple[float, float]],
    field: dict[int, float],
) -> tuple[float, float]:
    n0, n1, n2 = (int(node) for node in cell["node_ids"])
    x0, y0 = coordinates[n0]
    x1, y1 = coordinates[n1]
    x2, y2 = coordinates[n2]
    determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if determinant == 0.0:
        return math.nan, math.nan
    grad_x = (
        field[n0] * (y1 - y2)
        + field[n1] * (y2 - y0)
        + field[n2] * (y0 - y1)
    ) / determinant
    grad_y = (
        field[n0] * (x2 - x1)
        + field[n1] * (x0 - x2)
        + field[n2] * (x1 - x0)
    ) / determinant
    return grad_x, grad_y


def sentaurus_state(export_dir: Path) -> dict[str, dict[int, float]]:
    return {
        "phip": scalar_field(export_dir, "hQuasiFermiPotential"),
        "holes_m3": {
            node: value * 1.0e6
            for node, value in scalar_field(export_dir, "hDensity").items()
        },
    }


def sentaurus_srh_integral(
    export_dir: Path, silicon_cells: list[dict],
    coordinates_um: dict[int, tuple[float, float]],
) -> dict[str, float]:
    """Integrate signed Sentaurus srhRecombination over the 2-D mesh.

    TDR rate is cm^-3 s^-1 and mesh coordinates are um.  Conversion to A/um
    contributes 1e6 (cm^-3 to m^-3), 1e-12 (um2 to m2), and 1e-6
    (A/m to A/um), for a net factor of 1e-12 after multiplying by q.
    """
    rate = scalar_field(export_dir, "srhRecombination")
    net = generation = recombination = 0.0
    for cell in silicon_cells:
        ids = [int(node) for node in cell["node_ids"]]
        if not all(node in rate for node in ids):
            continue
        a, b, c = (coordinates_um[node] for node in ids)
        area_um2 = 0.5 * abs(
            (b[0] - a[0]) * (c[1] - a[1])
            - (c[0] - a[0]) * (b[1] - a[1])
        )
        value = statistics.mean(rate[node] for node in ids)
        net += value * area_um2
        generation += max(-value, 0.0) * area_um2
        recombination += max(value, 0.0) * area_um2
    factor = Q_C * 1.0e-12
    return {
        "srh_net_current_A_per_um": factor * net,
        "srh_generation_current_A_per_um": factor * generation,
        "srh_recombination_current_A_per_um": factor * recombination,
    }


def write_required_state(
    source: Path, output: Path, silicon_nodes: set[int], phip_shift_V: float,
    hole_override_m3: dict[int, float] | None = None,
) -> Path:
    rows = state_rows(source)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
        )
        for node in sorted(rows):
            row = rows[node]
            phip = row["phip"] - phip_shift_V if node in silicon_nodes else row["phip"]
            holes = (
                hole_override_m3[node]
                if hole_override_m3 is not None and node in hole_override_m3
                else row["holes_m3"]
            )
            writer.writerow([
                node, format(row["psi"], ".17g"), format(row["phin"], ".17g"),
                format(phip, ".17g"), format(row["electrons_m3"], ".17g"),
                format(holes, ".17g"),
            ])
    return output


def probe_node_holes(path: Path) -> dict[int, float]:
    values: dict[int, list[float]] = {}
    for row in read_csv(path):
        for endpoint in (0, 1):
            node = int(row[f"node{endpoint}"])
            values.setdefault(node, []).append(float(row[f"hole_density{endpoint}_m3"]))
    return {
        node: statistics.median(samples) for node, samples in values.items()
    }


def run_frozen(
    runner: Path, baseline: dict[str, Any], state: Path, bias: float,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    config = deepcopy(baseline)
    config["solver"]["method"] = "frozen_state"
    config["output_csv"] = str((output_dir / "curve.csv").resolve())
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    config["sweep"] = {
        "mode": "iv",
        "contact": "gate",
        "current_contact": "drain",
        "start": bias,
        "stop": bias,
        "step": 1.0,
        "bias_points": [bias],
        "initial_state_file": str(state.resolve()),
        "write_state_file": str((output_dir / "replayed_state.csv").resolve()),
        "frozen_state_compute_current": True,
        "write_vtk": False,
        "diagnostics": {
            "srh_balance": {
                "enabled": True,
                "material": "Si",
                "drain_contact": "drain",
                "substrate_contact": "substrate",
                "kcl_contacts": [
                    contact["name"] for contact in config["contacts"]
                ],
                "resolution_margin_ratio": 10.0,
                "csv_file": str((output_dir / "srh_balance.csv").resolve()),
            }
        },
    }
    config["_comment"] = (
        "Templates/LDMOS phip absolute-offset sensitivity; frozen-state only."
    )
    config_path = output_dir / "frozen.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output_dir / "frozen.log").resolve())],
        cwd=REPO, text=True, capture_output=True,
        env=runner_environment(), check=False,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)
    curve = read_csv(output_dir / "curve.csv")[-1]
    srh = read_csv(output_dir / "srh_balance.csv")[-1]
    return {
        "drain_total_current_A_per_um": float(curve["current_total_A_per_um"]),
        "drain_electron_current_A_per_um": float(curve["current_electron_A_per_um"]),
        "drain_hole_current_A_per_um": float(curve["current_hole_A_per_um"]),
        "srh_net_current_A_per_um": float(srh["srh_net_current_A_per_um"]),
        "srh_generation_current_A_per_um": float(
            srh["srh_generation_current_A_per_um"]
        ),
        "srh_recombination_current_A_per_um": float(
            srh["srh_recombination_current_A_per_um"]
        ),
        "artifacts": {
            "config": str(config_path.resolve()),
            "curve": str((output_dir / "curve.csv").resolve()),
            "srh_balance": str((output_dir / "srh_balance.csv").resolve()),
        },
    }


def state_name(bias: float) -> str:
    return f"g3_idvg_point_bias_{bias:.6f}".replace(".", "p") + ".csv"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vela-state-root", type=Path, required=True)
    parser.add_argument("--sentaurus-export-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hole-density-floor-m3", type=float, default=1.0e16)
    args = parser.parse_args()

    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    coordinates_um = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    coordinates = {
        node: (position[0] * 1.0e-6, position[1] * 1.0e-6)
        for node, position in coordinates_um.items()
    }
    silicon_nodes, silicon_cells = silicon_nodes_and_cells(mesh)
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {
        "schema": "vela.templates_ldmos_phip_srh_audit.v1",
        "mask_contract": {
            "all_silicon_nodes": len(silicon_nodes),
            "hole_populated_definition": (
                "max(Vela supplied p, Sentaurus supplied p) >= hole_density_floor_m3"
            ),
            "hole_density_floor_m3": args.hole_density_floor_m3,
            "acceptance_use": False,
        },
        "points": [],
    }
    for index, bias in enumerate(BIAS_POINTS_V):
        point_dir = output / f"vg_{bias:.6f}".replace(".", "p")
        vela_path = args.vela_state_root / state_name(bias)
        vela = state_rows(vela_path)
        sentaurus = sentaurus_state(args.sentaurus_export_root / f"export_{index:04d}")
        common = sorted(silicon_nodes & set(vela) & set(sentaurus["phip"]))
        delta = {node: vela[node]["phip"] - sentaurus["phip"][node] for node in common}
        offset = statistics.median(delta.values())
        centered = {node: delta[node] - offset for node in common}
        populated = {
            node for node in common
            if max(vela[node]["holes_m3"], sentaurus["holes_m3"][node])
            >= args.hole_density_floor_m3
        }
        populated_centered = {
            node: centered[node] for node in populated
        }
        gradient_delta = []
        gradient_relative = []
        populated_gradient_delta = []
        populated_gradient_relative = []
        for cell in silicon_cells:
            ids = {int(node) for node in cell["node_ids"]}
            if not ids.issubset(common):
                continue
            grad_vela = p1_gradient(
                cell, coordinates, {node: vela[node]["phip"] for node in ids}
            )
            grad_sentaurus = p1_gradient(cell, coordinates, sentaurus["phip"])
            difference = math.hypot(
                grad_vela[0] - grad_sentaurus[0],
                grad_vela[1] - grad_sentaurus[1],
            )
            scale = max(math.hypot(*grad_vela), math.hypot(*grad_sentaurus), 1.0)
            gradient_delta.append(difference)
            gradient_relative.append(difference / scale)
            if ids.issubset(populated):
                populated_gradient_delta.append(difference)
                populated_gradient_relative.append(difference / scale)

        qf_only = write_required_state(
            vela_path, point_dir / "phip_shift_qf_only.csv", silicon_nodes, offset,
        )
        shifted_probe = run_sg_probe(
            args.runner, baseline, qf_only, bias, point_dir / "mapping_probe",
        )
        reconstructed_holes = probe_node_holes(
            Path(shifted_probe["artifacts"]["sg_edges"])
        )
        mapping_state = write_required_state(
            vela_path, point_dir / "phip_shift_mapping_consistent.csv",
            silicon_nodes, offset, reconstructed_holes,
        )
        original_required = write_required_state(
            vela_path, point_dir / "original_required_columns.csv",
            silicon_nodes, 0.0,
        )
        hole_error_original = [
            math.log10(max(vela[node]["holes_m3"], 1.0))
            - math.log10(max(sentaurus["holes_m3"][node], 1.0))
            for node in populated
        ]
        hole_error_original_all = [
            math.log10(max(vela[node]["holes_m3"], 1.0))
            - math.log10(max(sentaurus["holes_m3"][node], 1.0))
            for node in common
        ]
        hole_error_shifted = [
            math.log10(max(reconstructed_holes[node], 1.0))
            - math.log10(max(sentaurus["holes_m3"][node], 1.0))
            for node in populated if node in reconstructed_holes
        ]
        hole_error_shifted_all = [
            math.log10(max(reconstructed_holes[node], 1.0))
            - math.log10(max(sentaurus["holes_m3"][node], 1.0))
            for node in common if node in reconstructed_holes
        ]
        mapping_closure = [
            math.log10(max(reconstructed_holes[node], 1.0))
            - math.log10(max(vela[node]["holes_m3"], 1.0))
            for node in populated if node in reconstructed_holes
        ]
        frozen = {
            "original": run_frozen(
                args.runner, baseline, original_required, bias,
                point_dir / "frozen_original",
            ),
            "phip_shift_qf_only": run_frozen(
                args.runner, baseline, qf_only, bias,
                point_dir / "frozen_phip_shift_qf_only",
            ),
            "phip_shift_mapping_consistent": run_frozen(
                args.runner, baseline, mapping_state, bias,
                point_dir / "frozen_phip_shift_mapping_consistent",
            ),
        }
        summary["points"].append({
            "bias_V": bias,
            "absolute_phip_offset_V": offset,
            "phip_delta_V": stats(delta.values()),
            "centered_phip_residual_V": stats(centered.values()),
            "hole_populated_centered_phip_residual_V": stats(
                populated_centered.values()
            ),
            "p1_phip_gradient_difference_V_per_m": stats(gradient_delta),
            "p1_phip_gradient_relative_difference": stats(gradient_relative),
            "hole_populated_p1_phip_gradient_difference_V_per_m": stats(
                populated_gradient_delta
            ),
            "hole_populated_p1_phip_gradient_relative_difference": stats(
                populated_gradient_relative
            ),
            "hole_populated_nodes": len(populated),
            "all_silicon_hole_density_original_minus_sentaurus_dex": stats(
                hole_error_original_all
            ),
            "hole_density_original_minus_sentaurus_dex": stats(hole_error_original),
            "all_silicon_hole_density_shifted_mapping_minus_sentaurus_dex": stats(
                hole_error_shifted_all
            ),
            "hole_density_shifted_mapping_minus_sentaurus_dex": stats(
                hole_error_shifted
            ),
            "shifted_mapping_minus_original_vela_holes_dex": stats(mapping_closure),
            "frozen_sensitivity": frozen,
            "sentaurus_native_srh_integral": sentaurus_srh_integral(
                args.sentaurus_export_root / f"export_{index:04d}",
                silicon_cells, coordinates_um,
            ),
        })

    offsets = [point["absolute_phip_offset_V"] for point in summary["points"]]
    summary["decision"] = {
        "offset_across_bias_V": stats(offsets),
        "classification": (
            "absolute_hole_qf_reference_offset_in_hole_populated_region; "
            "near-zero-hole QF remains unanchored"
        ),
        "ledger_status": "draft",
        "ledger_approval_allowed": False,
        "caveat": (
            "The shifted states are fixed-state sensitivities, not self-consistent solutions."
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8",
    )
    print(json.dumps(summary["decision"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
