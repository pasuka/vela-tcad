#!/usr/bin/env python3
"""Prepare direct-bordered SLOT-LDMOS BVDS probes for the Delaunay profile.

The two cases are deliberately named ``masetti_field`` and
``enhanced_lombardi_proxy``.  Vela does not yet implement the self-contained
Sentaurus IALMob model, so these decks must not be presented as a true IALMob
on/off comparison.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


OUTPUT_ROOT = "outputs/delaunay_direct_bordered_20260825"
HOMOTOPY_SCALES: tuple[float | None, ...] = (
    None, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0,
)
POISSON_CAP_PROBES: tuple[float, ...] = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0)
ATOMIC_MESH = {
    "node_volume_policy": "mixed_voronoi",
    "require_non_obtuse": False,
    "edge_coupling_policy": "assembled_delaunay_fallback",
}
ATOMIC_IMPACT = {
    "current_approximation": "element_edge_sg_gss_laux",
    "source_mapping_mode": "element_vertex_box_measure",
    "cell_reconstructed_midpoint_density": "bernoulli",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_atomic_profile(template: dict[str, Any]) -> None:
    if template.get("mesh_geometry") != ATOMIC_MESH:
        raise ValueError("template is not the atomic assembled-Delaunay mesh profile")
    impact = template.get("solver", {}).get("impact_ionization", {})
    for key, expected in ATOMIC_IMPACT.items():
        if impact.get(key) != expected:
            raise ValueError(f"template impact_ionization.{key} is not {expected!r}")


def configure_case(template: dict[str, Any], *, case: str,
                   seed_state: Path, inner_voltage_v: float,
                   outer_voltage_v: float, probe_only: bool,
                   allow_missing_seed: bool = False) -> dict[str, Any]:
    document = copy.deepcopy(template)
    validate_atomic_profile(document)
    if not allow_missing_seed and not seed_state.is_file():
        raise FileNotFoundError(seed_state)
    model = {
        "masetti_field": "masetti_field",
        "enhanced_lombardi_proxy": "masetti_field_lombardi",
    }.get(case)
    if model is None:
        raise ValueError(f"unknown case {case!r}")

    phase = "probe" if probe_only else "continuation"
    output = f"{OUTPUT_ROOT}/{case}/{phase}"
    # CSVWriter does not create parent directories at runtime.
    (Path(template["_bundle_root"]) / output).mkdir(parents=True, exist_ok=True)
    document["_comment"] = (
        "Direct-bordered assembled-Delaunay SLOT-LDMOS BVDS run. Enhanced "
        "Lombardi is only a surface-mobility proxy and is not Sentaurus IALMob."
    )
    document["_surface_mobility_ablation"] = {
        "case": case,
        "true_sentaurus_ialmob_implemented": False,
        "mesh_profile": "element_edge_sg_gss_laux_delaunay",
        "source_jacobian": "local_ad",
    }
    document["output_csv"] = f"{output}/iv.csv"
    solver = document["solver"]
    solver["abstol"] = 1.0e-9
    solver["stall_residual_floor"] = 1.0e-8
    solver["carrier_regularization_scale"] = 1.0e-8
    solver["carrier_row_qualified_stall_acceptance"] = True
    solver["carrier_row_convergence"] = {
        "mode": "enforce",
        "eps_row": 1.0e-3,
        "scale_floor": 1.0e-30,
        "min_source_scale_fraction": 0.0,
        "min_source_scale": 1.0e-13,
        "min_newton_max_iter": 80,
    }
    solver["line_search"] = True
    solver["line_search_mode"] = "merit"
    solver["global_continuity_closure"] = {
        "mode": "enforce",
        "tolerance": 1.0e-2,
        "source_floor": 1.0e-11,
        "absolute_tolerance": 1.0e-12,
        "linear_system_stabilization": "replace_free_row_sum",
    }
    solver["quasi_fermi_update_limit_V"] = 0.1
    solver["quasi_fermi_update_limit_minority_V"] = 0.05
    solver["quasi_fermi_update_limit_mode"] = "componentwise_poisson_recorrect"
    solver["impact_ionization"]["source_jacobian"] = "local_ad"
    mobility = solver["mobility"]
    mobility["model"] = model
    mobility["high_field_gradient_discretization"] = "transport_cell_vector"
    mobility["jacobian_field_derivatives"] = True
    mobility.pop("surface", None)
    if case == "enhanced_lombardi_proxy":
        mobility["surface"] = {
            "surface_interface": ["Silicon_1", "Oxide_1"]
        }

    sweep = document["sweep"]
    final_outer = outer_voltage_v if probe_only else 110000.0
    sweep["bias_points"] = [outer_voltage_v] if probe_only else [outer_voltage_v, final_outer]
    sweep["start"] = outer_voltage_v
    sweep["stop"] = final_outer
    sweep["initial_state_file"] = str(seed_state.resolve())
    sweep["write_state_file"] = f"{output}/final_state.csv"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["write_vtk"] = False
    sweep["continuation"] = {"arclength": {"enabled": False}}
    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": f"{output}/newton_attempts.csv",
            "iterations_csv_file": f"{output}/newton_iterations.csv",
            "rejected_state_directory": f"{output}/rejected_states",
        }
    }
    boundary = sweep["boundary_control"]
    boundary.update({
        "checkpoint_directory": f"{output}/checkpoints",
        "evaluation_csv": f"{output}/boundary_evaluations.csv",
        "resume": True,
        "adaptive_device_continuation": True,
    })
    circuit = sweep["external_circuit"]
    circuit.update({
        "solver": "coupled_newton",
        "coupled_linear_solver": "direct_bordered",
        "initial_inner_voltage_V": inner_voltage_v,
        "resistance_ohm_um": 1.0e12,
        "coupled_apply_device_update_limit": True,
        "coupled_damping_factor": 0.5,
        "coupled_equation_tolerance": 1.0e-6,
        "coupled_line_search_mode": "residual_filter",
        "coupled_filter_gamma": 1.0e-4,
        "coupled_filter_envelope_factor": 2.0,
        "coupled_max_line_search_steps": 12,
        "coupled_initial_outer_step_V": 12.5,
        "coupled_min_outer_step_V": 0.1,
        "coupled_max_outer_step_V": 30000.0,
        "coupled_outer_growth_factor": 4.0,
        "coupled_outer_shrink_factor": 0.5,
        "coupled_max_step_retries": 16,
        "max_inner_voltage_step_V": 1.0,
        "max_iterations": 80,
    })
    return document


def configure_device_reclosure(template: dict[str, Any], *, case: str,
                               seed_state: Path,
                               inner_voltage_v: float) -> dict[str, Any]:
    document = configure_case(
        template, case=case, seed_state=seed_state,
        inner_voltage_v=inner_voltage_v, outer_voltage_v=inner_voltage_v,
        probe_only=True)
    output = f"{OUTPUT_ROOT}/{case}/device_reclosure"
    (Path(document["_bundle_root"]) / output).mkdir(parents=True, exist_ok=True)
    document["_comment"] = (
        "Fixed-inner-voltage device reclosure after changing from the legacy "
        "cell reconstruction to the assembled-Delaunay atomic profile."
    )
    document["output_csv"] = f"{output}/iv.csv"
    document["solver"].setdefault("handoff", {})["gummel_max_iter"] = 100
    sweep = document["sweep"]
    sweep["bias_points"] = [inner_voltage_v]
    sweep["start"] = inner_voltage_v
    sweep["stop"] = inner_voltage_v
    sweep["external_circuit"]["enabled"] = False
    sweep["write_state_file"] = f"{output}/final_state.csv"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["boundary_control"]["checkpoint_directory"] = f"{output}/checkpoints"
    sweep["boundary_control"]["evaluation_csv"] = f"{output}/boundary_evaluations.csv"
    sweep["diagnostics"]["newton_history"].update({
        "attempts_csv_file": f"{output}/newton_attempts.csv",
        "iterations_csv_file": f"{output}/newton_iterations.csv",
        "rejected_state_directory": f"{output}/rejected_states",
    })
    return document


def configure_homotopy_stage(template: dict[str, Any], *, case: str,
                             seed_state: Path, inner_voltage_v: float,
                             stage_index: int,
                             avalanche_scale: float | None) -> dict[str, Any]:
    label = "avalanche_off" if avalanche_scale is None else (
        f"avalanche_{avalanche_scale:.0e}".replace("+", "")
    )
    output = f"{OUTPUT_ROOT}/{case}/homotopy/{stage_index:02d}_{label}"
    document = configure_case(
        template, case=case, seed_state=seed_state,
        inner_voltage_v=inner_voltage_v, outer_voltage_v=inner_voltage_v,
        probe_only=True, allow_missing_seed=stage_index > 0)
    (Path(document["_bundle_root"]) / output).mkdir(parents=True, exist_ok=True)
    document["_comment"] = (
        "Fixed-inner-voltage assembled-Delaunay homotopy stage. Avalanche is "
        "disabled first, then its conservative element-box source and local-AD "
        "Jacobian are restored through explicit source scales."
    )
    document["_delaunay_homotopy"] = {
        "stage_index": stage_index,
        "avalanche_enabled": avalanche_scale is not None,
        "avalanche_source_scale": avalanche_scale,
        "fixed_inner_voltage_V": inner_voltage_v,
        "seed_state": str(seed_state.resolve()),
    }
    document["output_csv"] = f"{output}/iv.csv"
    solver = document["solver"]
    if avalanche_scale is None:
        solver.pop("impact_ionization", None)
    else:
        solver["impact_ionization"]["source_geometry_scale"] = avalanche_scale
        solver["impact_ionization"]["source_jacobian"] = "local_ad"
    solver.setdefault("handoff", {})["gummel_max_iter"] = (
        100 if avalanche_scale is None else 20)

    sweep = document["sweep"]
    sweep["bias_points"] = [inner_voltage_v]
    sweep["start"] = inner_voltage_v
    sweep["stop"] = inner_voltage_v
    sweep["external_circuit"]["enabled"] = False
    sweep["write_state_file"] = f"{output}/final_state.csv"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["boundary_control"]["checkpoint_directory"] = f"{output}/checkpoints"
    sweep["boundary_control"]["evaluation_csv"] = f"{output}/boundary_evaluations.csv"
    sweep["diagnostics"]["newton_history"].update({
        "attempts_csv_file": f"{output}/newton_attempts.csv",
        "iterations_csv_file": f"{output}/newton_iterations.csv",
        "rejected_state_directory": f"{output}/rejected_states",
    })
    document.pop("_bundle_root", None)
    return document


def configure_poisson_reclosure(template: dict[str, Any], *, case: str,
                                seed_state: Path,
                                inner_voltage_v: float) -> dict[str, Any]:
    output = f"{OUTPUT_ROOT}/{case}/homotopy/00_poisson_only"
    document = configure_case(
        template, case=case, seed_state=seed_state,
        inner_voltage_v=inner_voltage_v, outer_voltage_v=inner_voltage_v,
        probe_only=True)
    (Path(document["_bundle_root"]) / output).mkdir(parents=True, exist_ok=True)
    document["_comment"] = (
        "Poisson-only assembled-Delaunay profile projection while preserving "
        "the legacy converged electron and hole quasi-Fermi potentials."
    )
    document["_delaunay_homotopy"] = {
        "stage_index": 0,
        "mode": "poisson_only_qf_preserved",
        "fixed_inner_voltage_V": inner_voltage_v,
        "seed_state": str(seed_state.resolve()),
    }
    document["output_csv"] = f"{output}/iv.csv"
    solver = document["solver"]
    solver.pop("impact_ionization", None)
    solver["method"] = "poisson_only"
    solver["max_iter"] = 200
    # The exact SLOT mesh contains a 0.044-degree triangle.  The raw mixed-
    # Voronoi Poisson direction can therefore be millions of volts even when
    # the profile-switch residual is only O(1e4).  Limit the complete Poisson
    # direction uniformly; do not split the atomic geometry bundle.
    solver["max_update"] = 0.1
    solver["line_search"] = True
    solver["line_search_mode"] = "merit"
    solver["poisson_only_quasi_fermi_initialization"] = "preserve_initial"
    solver.setdefault("handoff", {})["gummel_max_iter"] = 0
    sweep = document["sweep"]
    sweep["bias_points"] = [inner_voltage_v]
    sweep["start"] = inner_voltage_v
    sweep["stop"] = inner_voltage_v
    sweep["external_circuit"]["enabled"] = False
    sweep["write_state_file"] = f"{output}/final_state.csv"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["boundary_control"]["checkpoint_directory"] = f"{output}/checkpoints"
    sweep["boundary_control"]["evaluation_csv"] = f"{output}/boundary_evaluations.csv"
    sweep["diagnostics"]["newton_history"].update({
        "attempts_csv_file": f"{output}/newton_attempts.csv",
        "iterations_csv_file": f"{output}/newton_iterations.csv",
        "rejected_state_directory": f"{output}/rejected_states",
    })
    document.pop("_bundle_root", None)
    return document


def configure_poisson_block_probe(poisson_deck: dict[str, Any], *,
                                  seed_state: Path,
                                  inner_voltage_v: float,
                                  output_csv: str,
                                  max_update: float) -> dict[str, Any]:
    document = copy.deepcopy(poisson_deck)
    document["_comment"] = (
        "Read-only raw Poisson block direction audit before homotopy."
    )
    document["simulation_type"] = "newton_block_step_probe"
    document["state_file"] = str(seed_state.resolve())
    document["output_csv"] = output_csv
    document["block_modes"] = ["poisson_only"]
    document["solver"]["max_update"] = max_update
    document.pop("sweep", None)
    for contact in document["contacts"]:
        if contact["name"] == "drain":
            contact["bias"] = inner_voltage_v
    return document


def prepare(bundle: Path, *, off_seed: Path, off_inner_v: float,
            off_outer_v: float, proxy_seed: Path, proxy_inner_v: float,
            proxy_outer_v: float) -> dict[str, Any]:
    bundle = bundle.resolve()
    template = read_json(bundle / "simulation_05_avalanche_on_60v.json")
    template["_bundle_root"] = str(bundle)
    cases = {
        "masetti_field": (off_seed, off_inner_v, off_outer_v),
        "enhanced_lombardi_proxy": (proxy_seed, proxy_inner_v, proxy_outer_v),
    }
    configs: dict[str, dict[str, str]] = {}
    for case, (seed, inner, outer) in cases.items():
        poisson = configure_poisson_reclosure(
            template, case=case, seed_state=seed, inner_voltage_v=inner)
        poisson_name = f"simulation_delaunay_{case}_homotopy_00_poisson.json"
        write_json(bundle / poisson_name, poisson)
        probe_dir = bundle / OUTPUT_ROOT / case / "diagnostics"
        probe_dir.mkdir(parents=True, exist_ok=True)
        poisson_probe_names: list[str] = []
        for max_update in (0.0, *POISSON_CAP_PROBES):
            cap_label = "raw" if max_update == 0.0 else str(max_update).replace(".", "p")
            poisson_probe = configure_poisson_block_probe(
                poisson, seed_state=seed, inner_voltage_v=inner,
                output_csv=(
                    f"{OUTPUT_ROOT}/{case}/diagnostics/"
                    f"poisson_block_step_cap_{cap_label}.csv"),
                max_update=max_update)
            poisson_probe_name = (
                f"simulation_delaunay_{case}_poisson_block_probe_cap_{cap_label}.json"
            )
            write_json(bundle / poisson_probe_name, poisson_probe)
            poisson_probe_names.append(poisson_probe_name)
        homotopy_configs: list[str] = [poisson_name]
        stage_seed = (
            bundle / OUTPUT_ROOT / case / "homotopy" /
            "00_poisson_only" / "final_state.csv"
        )
        for stage_index, avalanche_scale in enumerate(HOMOTOPY_SCALES, start=1):
            stage = configure_homotopy_stage(
                template, case=case, seed_state=stage_seed,
                inner_voltage_v=inner, stage_index=stage_index,
                avalanche_scale=avalanche_scale)
            stage_name = (
                f"simulation_delaunay_{case}_homotopy_{stage_index:02d}.json"
            )
            write_json(bundle / stage_name, stage)
            homotopy_configs.append(stage_name)
            label = "avalanche_off" if avalanche_scale is None else (
                f"avalanche_{avalanche_scale:.0e}".replace("+", "")
            )
            stage_seed = (
                bundle / OUTPUT_ROOT / case / "homotopy" /
                f"{stage_index:02d}_{label}" / "final_state.csv"
            )
        reclosed_state = stage_seed

        probe = configure_case(
            template, case=case,
            seed_state=reclosed_state,
            inner_voltage_v=inner,
            outer_voltage_v=outer, probe_only=True,
            allow_missing_seed=True)
        if not reclosed_state.is_file():
            probe["_requires_completed_homotopy"] = str(reclosed_state)
        probe_name = f"simulation_delaunay_direct_bordered_{case}_probe.json"
        probe.pop("_bundle_root", None)
        write_json(bundle / probe_name, probe)

        probe_state = bundle / OUTPUT_ROOT / case / "probe" / "final_state.csv"
        continuation = configure_case(
            template, case=case,
            seed_state=probe_state if probe_state.is_file() else seed,
            inner_voltage_v=inner, outer_voltage_v=outer, probe_only=False)
        if not probe_state.is_file():
            continuation["_requires_completed_probe_state"] = str(probe_state)
        continuation.pop("_bundle_root", None)
        continuation_name = f"simulation_delaunay_direct_bordered_{case}_bvds.json"
        write_json(bundle / continuation_name, continuation)
        configs[case] = {
            "homotopy": homotopy_configs,
            "poisson_block_probes": poisson_probe_names,
            "probe": probe_name,
            "continuation": continuation_name,
        }

    manifest = {
        "schema": "vela.slot_ldmos.delaunay_direct_bordered_bvds.v1",
        "output_root": OUTPUT_ROOT,
        "true_sentaurus_ialmob_implemented": False,
        "breakdown_current_A_per_um": 1.0e-7,
        "global_continuity_contract": {
            "tolerance": 1.0e-2,
            "source_floor": 1.0e-11,
            "absolute_tolerance": 1.0e-12,
            "linear_system_stabilization": "replace_free_row_sum",
            "stall_residual_floor": 1.0e-8,
            "carrier_row_qualified_stall_acceptance": True,
            "carrier_row_tolerance": 1.0e-3,
            "carrier_row_min_source_scale": 1.0e-13,
        },
        "configs": configs,
    }
    write_json(bundle / "delaunay_direct_bordered_bvds_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--off-seed", type=Path, required=True)
    parser.add_argument("--off-inner-v", type=float, required=True)
    parser.add_argument("--off-outer-v", type=float, required=True)
    parser.add_argument("--proxy-seed", type=Path, required=True)
    parser.add_argument("--proxy-inner-v", type=float, required=True)
    parser.add_argument("--proxy-outer-v", type=float, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(
        args.bundle, off_seed=args.off_seed, off_inner_v=args.off_inner_v,
        off_outer_v=args.off_outer_v, proxy_seed=args.proxy_seed,
        proxy_inner_v=args.proxy_inner_v, proxy_outer_v=args.proxy_outer_v),
        sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
