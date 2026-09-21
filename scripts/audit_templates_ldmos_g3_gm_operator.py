#!/usr/bin/env python3
"""Decompose the exact LDMOS G3 maximum-gm segment on frozen Vela states.

The two endpoint states remain fixed.  Each replay changes one documented
operator family, so the reported current and slope changes are direct frozen-
state effects and must not be interpreted as nonlinear state feedback.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import (
        drain_nodes,
        run_sg_probe,
    )
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_idvg_shift_kcl import drain_nodes, run_sg_probe


BIASES = (1.0, 7.0 / 6.0)
VARIANTS = (
    "qualified",
    "mesh_default",
    "contact_hfs_off",
    "global_hfs_off",
    "mesh_default_global_hfs_off",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def curve_map(path: Path) -> dict[float, float]:
    return {
        round(float(row["bias_V"]), 12): float(row["current_total_A_per_um"])
        for row in read_csv(path)
    }


def state_path(root: Path, prefix: str, bias: float) -> Path:
    label = f"{bias:.6f}".replace(".", "p")
    return root / f"{prefix}_bias_{label}.csv"


def without_averagebox(config: dict[str, Any]) -> None:
    geometry = config.setdefault("mesh_geometry", {})
    for key in (
        "carrier_transport_couple_profile",
        "external_averagebox_couples_file",
        "external_averagebox_expected_edges",
    ):
        geometry.pop(key, None)


def without_contact_hfs(config: dict[str, Any]) -> None:
    mobility = config["solver"]["mobility"]
    mobility["contact_electric_field_fallback"] = False
    mobility.pop("contact_electric_field_fallback_scope", None)
    mobility.pop("contact_electric_field_fallback_mode", None)


def without_global_hfs(config: dict[str, Any]) -> None:
    # The qualified G3 model is constant low-field mobility followed by HFS.
    # Selecting constant preserves that low-field model while removing HFS.
    config["solver"]["mobility"] = {"model": "constant"}


def variant_config(base: dict[str, Any], name: str) -> dict[str, Any]:
    if name not in VARIANTS:
        raise ValueError(f"unknown gm operator variant: {name}")
    config = deepcopy(base)
    if name in ("mesh_default", "mesh_default_global_hfs_off"):
        without_averagebox(config)
    if name == "contact_hfs_off":
        without_contact_hfs(config)
    if name in ("global_hfs_off", "mesh_default_global_hfs_off"):
        without_global_hfs(config)
    config["_comment"] = (
        "Templates/LDMOS G3 exact maximum-gm frozen-state operator replay; "
        f"variant={name}."
    )
    return config


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    ucrt = r"D:\msys64\ucrt64\bin"
    env["PATH"] = ucrt + os.pathsep + env.get("PATH", "")
    return env


def run_edge_mobility_probe(
    runner: Path,
    config: dict[str, Any],
    state: Path,
    bias: float,
    output: Path,
) -> Path:
    probe = deepcopy(config)
    probe.pop("sweep", None)
    probe["simulation_type"] = "edge_mobility_probe"
    probe["state_file"] = str(state.resolve())
    path = output / "edge_mobility.csv"
    probe["output_csv"] = str(path.resolve())
    for contact in probe["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    config_path = output / "edge_mobility.json"
    config_path.write_text(json.dumps(probe, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (output / "edge_mobility.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise RuntimeError(f"edge mobility probe failed: {config_path}")
    return path


def flux_weighted_mobility(
    sg_path: Path,
    mobility_path: Path,
    contact_nodes: set[int] | None = None,
) -> dict[str, float | int]:
    sg = {int(row["edge_id"]): row for row in read_csv(sg_path)}
    mobility = {int(row["edge_id"]): row for row in read_csv(mobility_path)}
    rows: list[tuple[float, dict[str, str]]] = []
    for edge in set(sg) & set(mobility):
        srow = sg[edge]
        if contact_nodes is not None:
            crosses_cut = (
                (int(srow["node0"]) in contact_nodes)
                != (int(srow["node1"]) in contact_nodes)
            )
            if not crosses_cut:
                continue
        weight = abs(float(srow["electron_particle_line_flux_per_m_s"]))
        if weight > 0.0:
            rows.append((weight, mobility[edge]))
    total = sum(weight for weight, _ in rows)

    def mean(field: str) -> float:
        return sum(weight * float(row[field]) for weight, row in rows) / max(
            total, 1.0e-300
        )

    return {
        "active_flux_edges": len(rows),
        "electron_mobility_m2_V_s": mean(
            "electron_final_mobility_m2_V_s"
        ),
        "electron_low_field_mobility_m2_V_s": mean(
            "electron_low_field_mobility_m2_V_s"
        ),
        "electron_mobility_limiter": mean("electron_mobility_limiter"),
        "electron_mobility_field_V_m": mean("electron_mobility_field_V_m"),
        "electron_qf_field_V_m": mean("electron_qf_field_V_m"),
        "electric_field_V_m": mean("electric_field_V_m"),
    }


def segment_metrics(
    currents: dict[float, float], reference_gm: float,
) -> dict[str, float]:
    growth = currents[BIASES[1]] / currents[BIASES[0]]
    gm = (currents[BIASES[1]] - currents[BIASES[0]]) / (
        BIASES[1] - BIASES[0]
    )
    ratio = gm / reference_gm
    return {
        "gm_A_per_um_V": gm,
        "ratio_to_sentaurus_gm": ratio,
        "relative_error_to_sentaurus_gm": abs(ratio - 1.0),
        "endpoint_current_growth_ratio": growth,
        "log_current_slope_dex_per_V": math.log10(abs(growth)) / (
            BIASES[1] - BIASES[0]
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--state-prefix", default="g3_averagebox_point")
    parser.add_argument("--reference-curve", type=Path, required=True)
    parser.add_argument("--vela-curve", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    base = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    if base["solver"].get("contact_boundary_reconstruction") != "legacy_node_local":
        raise ValueError("gm operator audit requires legacy_node_local")
    geometry = base.get("mesh_geometry", {})
    if geometry.get("carrier_transport_couple_profile") != (
        "templates_ldmos_external_averagebox"
    ):
        raise ValueError("gm operator audit requires the qualified AverageBox profile")
    mobility = base["solver"].get("mobility", {})
    if mobility.get("model") != "constant_field":
        raise ValueError("gm operator audit requires the qualified HFS model")
    if not mobility.get("contact_electric_field_fallback", False):
        raise ValueError("gm operator audit requires qualified contact HFS")
    drain_contact_nodes = drain_nodes(Path(base["mesh_file"]))

    reference = curve_map(args.reference_curve)
    vela = curve_map(args.vela_curve)
    reference_currents = {
        bias: reference[round(bias, 12)] for bias in BIASES
    }
    reference_gm = segment_metrics(reference_currents, 1.0)["gm_A_per_um_V"]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    currents: dict[str, dict[float, float]] = {
        name: {} for name in VARIANTS
    }
    points: list[dict[str, Any]] = []

    for bias in BIASES:
        state = state_path(args.state_root, args.state_prefix, bias)
        if not state.is_file():
            raise FileNotFoundError(state)
        label = f"vg_{bias:.6f}".replace(".", "p")
        point: dict[str, Any] = {
            "bias_V": bias,
            "sentaurus_terminal_A_per_um": reference[round(bias, 12)],
            "vela_curve_A_per_um": vela[round(bias, 12)],
            "variants": {},
        }
        for name in VARIANTS:
            case = output / label / name
            config = variant_config(base, name)
            sg = run_sg_probe(args.runner, config, state, bias, case)
            mobility_path = run_edge_mobility_probe(
                args.runner, config, state, bias, case
            )
            current = float(sg["drain_cut"]["total_A_per_um"])
            currents[name][bias] = current
            point["variants"][name] = {
                "drain_current_A_per_um": current,
                "ratio_to_sentaurus": current / reference[round(bias, 12)],
                "flux_weighted_mobility": {
                    "all_active_edges": flux_weighted_mobility(
                        Path(sg["artifacts"]["sg_edges"]), mobility_path
                    ),
                    "drain_cut": flux_weighted_mobility(
                        Path(sg["artifacts"]["sg_edges"]),
                        mobility_path,
                        drain_contact_nodes,
                    ),
                },
                "artifacts": {
                    **sg["artifacts"],
                    "edge_mobility": str(mobility_path.resolve()),
                },
            }
        qualified = point["variants"]["qualified"]["drain_current_A_per_um"]
        point["qualified_probe_to_curve_relative_difference"] = (
            qualified / point["vela_curve_A_per_um"] - 1.0
        )
        points.append(point)

    segments = {
        name: segment_metrics(values, reference_gm)
        for name, values in currents.items()
    }
    segments["sentaurus"] = segment_metrics(reference_currents, reference_gm)
    segments["vela_curve"] = segment_metrics(
        {bias: vela[round(bias, 12)] for bias in BIASES}, reference_gm
    )
    qualified_gm = segments["qualified"]["gm_A_per_um_V"]
    deficit = reference_gm - qualified_gm
    direct = {
        "averagebox_delta_gm_A_per_um_V": qualified_gm
        - segments["mesh_default"]["gm_A_per_um_V"],
        "contact_hfs_delta_gm_A_per_um_V": qualified_gm
        - segments["contact_hfs_off"]["gm_A_per_um_V"],
        "global_hfs_delta_gm_A_per_um_V": qualified_gm
        - segments["global_hfs_off"]["gm_A_per_um_V"],
        "remaining_qualified_gm_deficit_A_per_um_V": deficit,
    }
    direct["averagebox_delta_over_deficit"] = (
        direct["averagebox_delta_gm_A_per_um_V"] / deficit
    )
    direct["contact_hfs_delta_over_deficit"] = (
        direct["contact_hfs_delta_gm_A_per_um_V"] / deficit
    )
    direct["global_hfs_delta_over_deficit"] = (
        direct["global_hfs_delta_gm_A_per_um_V"] / deficit
    )
    operator_log_slopes = [
        segments[name]["log_current_slope_dex_per_V"] for name in VARIANTS
    ]
    direct["qualified_log_slope_gap_to_sentaurus_dex_per_V"] = (
        segments["sentaurus"]["log_current_slope_dex_per_V"]
        - segments["qualified"]["log_current_slope_dex_per_V"]
    )
    direct["operator_variant_log_slope_span_dex_per_V"] = (
        max(operator_log_slopes) - min(operator_log_slopes)
    )
    direct["operator_span_over_qualified_log_slope_gap"] = (
        direct["operator_variant_log_slope_span_dex_per_V"]
        / direct["qualified_log_slope_gap_to_sentaurus_dex_per_V"]
    )
    direct["classification"] = "frozen_operator_slope_insensitive"
    report = {
        "schema": "vela.templates_ldmos.g3_gm_operator_audit.v1",
        "contracts": {
            "states": "frozen converged Vela node-local plus AverageBox states",
            "biases_V": list(BIASES),
            "state_feedback_included": False,
            "sentaurus_field_state_required_for_feedback": True,
            "variants": list(VARIANTS),
        },
        "points": points,
        "segments": segments,
        "direct_operator_attribution": direct,
    }
    summary = output / "summary.json"
    summary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "summary": str(summary),
        "segments": segments,
        "direct_operator_attribution": direct,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
