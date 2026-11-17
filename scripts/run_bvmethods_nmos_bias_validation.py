#!/usr/bin/env python3
"""Run postprocessed and self-consistent BVmethods NMOS bias validation."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REQUESTED_BIASES = [0.0, 1.0, 2.0, 4.0, 5.0, 6.0]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def resolve_input(config_path: Path, value: str) -> str:
    path = Path(value)
    if not path.is_absolute():
        path = config_path.parent / path
    return str(path.resolve())


def build_config(
    base_path: Path,
    base: dict[str, Any],
    case_dir: Path,
    initial_state: Path,
    coupling_mode: str,
    solver_method: str,
    biases: list[float],
    reltol: float,
    abstol: float,
    stall_residual_floor: float,
    poisson_stall_relative_increase: float,
    poisson_stall_contact_qf_drop_limit_V: float,
    carrier_statistics: str,
    lightweight: bool,
    lag_high_field_mobility: bool,
    continuity_row_scaling: bool,
    qf_bounds_mode: str,
    qf_bounds_margin_V: float,
    qf_bounds_min_carrier_density_m3: float,
    carrier_row_convergence_mode: str,
    carrier_row_eps: float,
    carrier_row_min_source_scale: float,
    carrier_row_recovery_mode: str,
    carrier_row_recovery_max_attempts: int,
    carrier_row_recovery_max_cycles: int,
    newton_history: bool,
    qf_update_limit_V: float,
    minority_qf_update_limit_V: float,
    impact_driving_force: str,
    impact_current_approximation: str,
    impact_current_magnitude_mode: str,
    path_ionization_integrals: bool,
    path_ionization_max_paths: int,
    path_ionization_break_rank: int,
    path_ionization_break_value: float,
    path_ionization_driving_force: str,
    path_ionization_stop_field_V_per_m: float,
    path_ionization_tracing_mode: str,
    path_ionization_seed_field_V_per_m: float,
) -> Path:
    cfg = deepcopy(base)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in cfg:
            cfg[key] = resolve_input(base_path, str(cfg[key]))
    cfg["doping"] = []
    cfg["output_csv"] = str((case_dir / "sweep.csv").resolve())

    solver = cfg.setdefault("solver", {})
    solver.update(
        {
            "method": solver_method,
            "max_iter": 120,
            "reltol": reltol,
            "abstol": abstol,
            "stall_residual_floor": stall_residual_floor,
            "poisson_line_search_stall_relative_increase":
                poisson_stall_relative_increase,
            "poisson_line_search_stall_contact_majority_qf_drop_limit_V":
                poisson_stall_contact_qf_drop_limit_V,
            "carrier_row_qualified_stall_acceptance": True,
            "residual_scales": {
                "psi": 1138.7290351540657,
                "phin": 1.0,
                "phip": 1.0,
            },
            "line_search": True,
            "verbose": False,
            "warm_start": True,
            "bandgap_narrowing": "old_slotboom",
            "recombination": ["srh"],
            "mobility": {
                "model": "masetti_field",
                "high_field_driving_force": "quasi_fermi_gradient",
                "jacobian_field_derivatives": not lag_high_field_mobility,
            },
            "impact_ionization": {
                "model": "van_overstraeten",
                "coupling_mode": coupling_mode,
                "driving_force": impact_driving_force,
                "generation": "current_density",
                "current_approximation": impact_current_approximation,
                "current_magnitude_mode": impact_current_magnitude_mode,
            },
            "quasi_fermi_update_limit_V": qf_update_limit_V,
            "quasi_fermi_update_limit_minority_V": minority_qf_update_limit_V,
            "carrier_statistics": carrier_statistics,
        }
    )
    if coupling_mode == "self_consistent":
        # Do not let a loose sweep tolerance accept the no-impact restart
        # without taking a feedback-coupled Newton step.
        solver["abstol"] = 1.0e-12
    if solver_method == "gummel_newton":
        solver["handoff"] = {
            "fallback": "none",
            "require_gummel_convergence": False,
            "gummel_max_iter": 80,
            "newton_max_iter": 120,
        }
    if carrier_row_convergence_mode != "off":
        solver["carrier_row_convergence"] = {
            "mode": carrier_row_convergence_mode,
            "eps_row": carrier_row_eps,
            "min_source_scale": carrier_row_min_source_scale,
            "diagnostic_csv": str((case_dir / "carrier_row_convergence.csv").resolve()),
            "trace_csv": str((case_dir / "carrier_row_trace.csv").resolve()),
            "trace_first_iterations": 8,
            "trace_every_iterations": 1,
            "recovery": {
                "mode": carrier_row_recovery_mode,
                "max_attempts": carrier_row_recovery_max_attempts,
                "max_cycles": carrier_row_recovery_max_cycles,
                "density_change_reltol": 1.0e-8,
            },
        }
    if continuity_row_scaling:
        solver["continuity_row_scaling"] = {
            "enabled": True,
            "flux_fraction": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale": carrier_row_min_source_scale,
            "min_weight": 1.0e-12,
            "max_weight": 1.0e12,
        }

    sweep = cfg.setdefault("sweep", {})
    sweep.pop("initialization", None)
    sweep.update(
        {
            "mode": "bv_reverse",
            "contact": "drain",
            "current_contact": "drain",
            "bias_points": biases,
            "start": biases[0],
            "stop": biases[-1],
            "step": 0.05,
            "initial_step": 0.05,
            "min_step": 1.0e-8,
            "max_step": 0.1,
            "growth_factor": 1.35,
            "shrink_factor": 0.5,
            "max_retries": 29,
            "initial_state_file": str(initial_state.resolve()),
            "write_state_file": str((case_dir / "last_state.csv").resolve()),
            "write_state_every_point_prefix": str(
                (case_dir / "states" / "accepted_state").resolve()
            ),
            "write_vtk": not lightweight,
            "vtk_prefix": str((case_dir / "vtk" / "state").resolve()),
        }
    )
    diagnostics = {
        "qf_bounds": {
            "enabled": True,
            "mode": qf_bounds_mode,
            "margin_V": qf_bounds_margin_V,
            "min_carrier_density_m3": qf_bounds_min_carrier_density_m3,
            "csv_file": str((case_dir / "qf_bounds.csv").resolve()),
        },
    }
    if newton_history:
        diagnostics["newton_history"] = {
            "enabled": True,
            "csv_file": str((case_dir / "newton_history.csv").resolve()),
        }
    if not lightweight:
        diagnostics.update(
            {
                "sg_avalanche_edges": {
                    "enabled": True,
                    "csv_file": str((case_dir / "sg_avalanche_edges.csv").resolve()),
                },
                "terminal_current_method_compare": {
                    "enabled": True,
                    "contacts": ["drain", "source", "substrate"],
                    "csv_file": str(
                        (case_dir / "terminal_current_method_compare.csv").resolve()
                    ),
                },
                "continuity_balance": {
                    "enabled": True,
                    "contacts": ["drain", "source", "substrate"],
                    "csv_file": str((case_dir / "continuity_balance.csv").resolve()),
                },
            }
        )
    if path_ionization_integrals:
        diagnostics["path_ionization_integrals"] = {
            "enabled": True,
            "csv_file": str(
                (case_dir / "path_ionization_integrals.csv").resolve()
            ),
            "segments_csv_file": str(
                (case_dir / "path_ionization_integral_segments.csv").resolve()
            ),
            "max_paths": path_ionization_max_paths,
            "break_rank": path_ionization_break_rank,
            "break_value": path_ionization_break_value,
            "driving_force": path_ionization_driving_force,
            "stop_field_V_per_m": path_ionization_stop_field_V_per_m,
            "tracing_mode": path_ionization_tracing_mode,
            "seed_field_V_per_m": path_ionization_seed_field_V_per_m,
        }
    sweep["diagnostics"] = diagnostics
    cfg["_validation_case"] = {
        "coupling_mode": coupling_mode,
        "requested_biases_V": biases,
        "sentaurus_reference_BV_V": 6.377494277837012,
    }
    path = case_dir / "simulation.json"
    write_json(path, cfg)
    return path


def summarize_edges(case_dir: Path) -> None:
    edge_path = case_dir / "sg_avalanche_edges.csv"
    if not edge_path.exists():
        return
    aggregate: dict[str, dict[str, float]] = {}
    with edge_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            bias = row["bias_V"]
            item = aggregate.setdefault(
                bias,
                {
                    "max_electron_alpha_m_inv": 0.0,
                    "max_hole_alpha_m_inv": 0.0,
                    "sum_edge_source_integral": 0.0,
                },
            )
            item["max_electron_alpha_m_inv"] = max(
                item["max_electron_alpha_m_inv"], float(row["electron_alpha_m_inv"])
            )
            item["max_hole_alpha_m_inv"] = max(
                item["max_hole_alpha_m_inv"], float(row["hole_alpha_m_inv"])
            )
            item["sum_edge_source_integral"] += float(row["edge_source_integral"])
    with (case_dir / "avalanche_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        fields = [
            "bias_V",
            "max_electron_alpha_m_inv",
            "max_hole_alpha_m_inv",
            "sum_edge_source_integral",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for bias, item in sorted(aggregate.items(), key=lambda pair: float(pair[0])):
            writer.writerow({"bias_V": bias, **item})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--initial-state", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--node-doping-file",
        type=Path,
        help="Override the base deck node_doping_file (for corrected TDR reimports)",
    )
    parser.add_argument("--runner", type=Path)
    parser.add_argument(
        "--modes", default="postprocess_only,self_consistent",
        help="Comma-separated impact coupling modes",
    )
    parser.add_argument(
        "--impact-driving-force",
        choices=(
            "electric_field", "quasi_fermi_gradient",
            "grad_potential_parallel_j", "effective_field_parallel_j",
            "eparallel",
        ),
        help="Override solver.impact_ionization.driving_force",
    )
    parser.add_argument(
        "--impact-current-approximation",
        help="Override solver.impact_ionization.current_approximation",
    )
    parser.add_argument(
        "--impact-current-magnitude-mode",
        choices=("edge_scalar_abs", "dual_face_vector_mag"),
        help="Override solver.impact_ionization.current_magnitude_mode",
    )
    parser.add_argument(
        "--impact-eparallel-field-recovery",
        choices=("edge_adjacent_cells", "nodal_vertex_star"),
        help="Override solver.impact_ionization.eparallel_field_recovery",
    )
    parser.add_argument(
        "--impact-source-mapping-mode",
        choices=(
            "node_F_node_alpha_node_G",
            "edge_F_edge_alpha_edge_G_to_node",
            "cell_F_cell_alpha_cell_G_to_node",
            "nodal_eparallel_p1",
        ),
        help="Override solver.impact_ionization.source_mapping_mode",
    )
    parser.add_argument(
        "--solver-method",
        choices=("newton", "gummel_newton"),
        default="newton",
    )
    parser.add_argument(
        "--biases", default=",".join(str(value) for value in REQUESTED_BIASES)
    )
    parser.add_argument(
        "--reltol", type=float, default=0.0,
        help="Disable initial-residual-relative early acceptance for BV branch validation",
    )
    parser.add_argument("--abstol", type=float, default=1.0e-11)
    parser.add_argument(
        "--stall-residual-floor",
        type=float,
        default=2.0e-9,
        help=(
            "Benchmark-specific normalized residual noise floor; local carrier "
            "row convergence remains independently enforced"
        ),
    )
    parser.add_argument(
        "--poisson-stall-relative-increase",
        type=float,
        default=5.0e-5,
    )
    parser.add_argument(
        "--poisson-stall-contact-qf-drop-limit-v",
        type=float,
        default=1.0e-7,
    )
    parser.add_argument(
        "--carrier-statistics",
        choices=("boltzmann", "fermi_dirac"),
        default="fermi_dirac",
    )
    parser.add_argument(
        "--lightweight",
        action="store_true",
        help="Disable per-step VTK and heavy edge/current diagnostics for trunk sweeps",
    )
    parser.add_argument(
        "--lag-high-field-mobility",
        action="store_true",
        help=(
            "Evaluate the same high-field mobility model while omitting its "
            "quasi-Fermi-field derivatives from the Newton Jacobian"
        ),
    )
    parser.add_argument(
        "--mobility-model",
        choices=(
            "masetti_field", "masetti_surface", "masetti_field_surface",
            "masetti_lombardi", "masetti_field_lombardi",
        ),
        default="masetti_field",
    )
    parser.add_argument("--surface-theta-electron-m-per-v", type=float)
    parser.add_argument("--surface-theta-hole-m-per-v", type=float)
    parser.add_argument("--surface-beta", type=float, default=1.0)
    parser.add_argument("--surface-min-factor", type=float, default=0.05)
    parser.add_argument("--surface-region", default="R.Substrate")
    parser.add_argument(
        "--surface-interface",
        default="R.Substrate,R.Gateox",
        help="Comma-separated semiconductor and insulator region names",
    )
    parser.add_argument(
        "--srh-doping-dependence",
        action="store_true",
        help=(
            "Enable the Sentaurus Scharfetter SRH(DopingDep) lifetime law "
            "with the BVmethods pp6 defaults"
        ),
    )
    parser.add_argument(
        "--continuity-row-scaling",
        action="store_true",
        help="Enable source/flux-aware left scaling of Newton continuity rows",
    )
    parser.add_argument(
        "--qf-bounds-mode",
        choices=("warn", "reject_and_recover"),
        default="reject_and_recover",
    )
    parser.add_argument("--qf-bounds-margin-v", type=float, default=1.0)
    parser.add_argument(
        "--qf-bounds-min-carrier-density-m3",
        type=float,
        default=1.0e6,
        help=(
            "Ignore finite quasi-Fermi excursions when the corresponding "
            "carrier density is below this numerical observability floor"
        ),
    )
    parser.add_argument(
        "--carrier-row-convergence-mode",
        choices=("off", "report", "enforce"),
        default="enforce",
    )
    parser.add_argument("--carrier-row-eps", type=float, default=1.0e-3)
    parser.add_argument(
        "--carrier-row-min-source-scale",
        type=float,
        default=1.0e-14,
        help="Ignore local source/flux balances below this native numerical-noise floor",
    )
    parser.add_argument(
        "--carrier-row-recovery-mode",
        choices=("off", "gummel_density"),
        default="gummel_density",
    )
    parser.add_argument("--carrier-row-recovery-max-attempts", type=int, default=4)
    parser.add_argument("--carrier-row-recovery-max-cycles", type=int, default=3)
    parser.add_argument(
        "--newton-history", action="store_true",
        help="Write detailed nonlinear iteration and line-search history",
    )
    parser.add_argument("--qf-update-limit-v", type=float, default=0.01)
    parser.add_argument("--minority-qf-update-limit-v", type=float, default=0.003)
    parser.add_argument(
        "--predictor-mode",
        choices=("none", "constant", "linear", "secant"),
        default="none",
        help="Continuation predictor for psi/phin/phip on high-voltage branches",
    )
    parser.add_argument("--path-ionization-integrals", action="store_true")
    parser.add_argument("--path-ionization-max-paths", type=int, default=3)
    parser.add_argument("--path-ionization-break-rank", type=int, default=0)
    parser.add_argument("--path-ionization-break-value", type=float, default=1.0)
    parser.add_argument(
        "--path-ionization-stop-field-v-per-m", type=float, default=0.0,
        help="Stop a traced path when |E| drops below this SI threshold",
    )
    parser.add_argument(
        "--path-ionization-electron-stop-field-v-per-m", type=float, default=0.0,
        help="Terminate electron-injection support on its own Eparallel threshold",
    )
    parser.add_argument(
        "--path-ionization-hole-stop-field-v-per-m", type=float, default=0.0,
        help="Terminate hole-injection support on its own Eparallel threshold",
    )
    parser.add_argument(
        "--path-ionization-mean-definition",
        choices=(
            "carrier_integral_arithmetic",
            "carrier_alpha_length_arithmetic",
        ),
        default="carrier_integral_arithmetic",
    )
    parser.add_argument(
        "--path-ionization-break-ordering",
        choices=("path_mean", "carrier_integrals"),
        default="path_mean",
        help="Select validated path-mean or experimental flattened-carrier ordering",
    )
    parser.add_argument(
        "--path-ionization-tracing-mode",
        choices=("edge_graph", "continuous_cell"),
        default="edge_graph",
    )
    parser.add_argument(
        "--path-ionization-path-retention",
        choices=(
            "all_seed_trajectories",
            "numbered_peak_groups",
            "distinct_local_maxima",
            "corridor_deduplicated",
        ),
        default="numbered_peak_groups",
    )
    parser.add_argument(
        "--path-ionization-seed-mode",
        choices=("cell_local_maxima", "nodal_local_maxima"),
        default="nodal_local_maxima",
    )
    parser.add_argument(
        "--path-ionization-tracing-vector",
        choices=(
            "electric_field",
            "electron_current",
            "hole_current",
            "electron_qf_gradient",
            "hole_qf_gradient",
            "cell_electric_field",
            "electric_field_rk4",
            "sentaurus_eparallel_adaptive",
        ),
        default="electric_field",
    )
    parser.add_argument(
        "--path-ionization-tracing-direction",
        choices=("bidirectional", "along_vector", "opposite_vector"),
        default="bidirectional",
    )
    parser.add_argument(
        "--path-ionization-tracing-current-relative-floor",
        type=float,
        default=1.0e-8,
    )
    parser.add_argument(
        "--path-ionization-tracing-qf-relative-floor",
        type=float,
        default=5.1e-3,
    )
    parser.add_argument(
        "--path-ionization-seed-field-v-per-m", type=float, default=0.0,
    )
    parser.add_argument(
        "--path-ionization-driving-force",
        choices=(
            "solver", "electric_field", "quasi_fermi_gradient",
            "grad_potential_parallel_j", "effective_field_parallel_j",
            "eparallel",
        ),
        default="solver",
    )
    args = parser.parse_args()

    base_path = args.base_config.resolve()
    initial_state = args.initial_state.resolve()
    out_dir = args.out_dir.resolve()
    biases = [float(value) for value in args.biases.split(",") if value.strip()]
    modes = [value.strip() for value in args.modes.split(",") if value.strip()]
    runner = args.runner
    if runner is None:
        name = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
        runner = REPO / "build-release" / name
    runner = runner.resolve()
    base = read_json(base_path)
    impact = base.setdefault("solver", {}).setdefault("impact_ionization", {})
    if args.impact_driving_force is not None:
        impact["driving_force"] = args.impact_driving_force
    if args.impact_current_approximation is not None:
        impact["current_approximation"] = args.impact_current_approximation
    if args.impact_current_magnitude_mode is not None:
        impact["current_magnitude_mode"] = args.impact_current_magnitude_mode
    if args.node_doping_file is not None:
        base["node_doping_file"] = str(args.node_doping_file.resolve())
    impact_driving_force = str(
        impact.get("driving_force", "quasi_fermi_gradient")
    )
    impact_current_approximation = str(
        impact.get("current_approximation", "density_gradient")
    )
    impact_current_magnitude_mode = str(
        impact.get("current_magnitude_mode", "edge_scalar_abs")
    )

    overall = 0
    for mode in modes:
        case_dir = out_dir / mode
        config = build_config(
            base_path, base, case_dir, initial_state, mode,
            args.solver_method, biases,
            args.reltol, args.abstol, args.stall_residual_floor,
            args.poisson_stall_relative_increase,
            args.poisson_stall_contact_qf_drop_limit_v,
            args.carrier_statistics,
            args.lightweight, args.lag_high_field_mobility,
            args.continuity_row_scaling,
            args.qf_bounds_mode, args.qf_bounds_margin_v,
            args.qf_bounds_min_carrier_density_m3,
            args.carrier_row_convergence_mode, args.carrier_row_eps,
            args.carrier_row_min_source_scale,
            args.carrier_row_recovery_mode,
            args.carrier_row_recovery_max_attempts,
            args.carrier_row_recovery_max_cycles,
            args.newton_history,
            args.qf_update_limit_v,
            args.minority_qf_update_limit_v,
            impact_driving_force,
            impact_current_approximation,
            impact_current_magnitude_mode,
            args.path_ionization_integrals,
            args.path_ionization_max_paths,
            args.path_ionization_break_rank,
            args.path_ionization_break_value,
            args.path_ionization_driving_force,
            args.path_ionization_stop_field_v_per_m,
            args.path_ionization_tracing_mode,
            args.path_ionization_seed_field_v_per_m,
        )
        if (
            args.impact_eparallel_field_recovery is not None
            or args.impact_source_mapping_mode is not None
        ):
            runtime = read_json(config)
            runtime_impact = runtime["solver"]["impact_ionization"]
            if args.impact_eparallel_field_recovery is not None:
                runtime_impact["eparallel_field_recovery"] = (
                    args.impact_eparallel_field_recovery
                )
            if args.impact_source_mapping_mode is not None:
                runtime_impact["source_mapping_mode"] = (
                    args.impact_source_mapping_mode
                )
            write_json(config, runtime)
        runtime = read_json(config)
        runtime_sweep = runtime.setdefault("sweep", {})
        runtime_sweep.pop("predictor", None)
        runtime_sweep.setdefault("continuation", {})["predictor"] = {
            "mode": args.predictor_mode,
        }
        write_json(config, runtime)
        if args.mobility_model != "masetti_field":
            lombardi = args.mobility_model in {
                "masetti_lombardi", "masetti_field_lombardi"
            }
            if not lombardi and (
                args.surface_theta_electron_m_per_v is None
                or args.surface_theta_hole_m_per_v is None
            ):
                parser.error(
                    "surface mobility requires both electron and hole theta values"
                )
            runtime = read_json(config)
            mobility = runtime.setdefault("solver", {}).setdefault("mobility", {})
            mobility["model"] = args.mobility_model
            mobility["surface"] = {
                "surface_region": args.surface_region,
                "surface_interface": [
                    value.strip()
                    for value in args.surface_interface.split(",")
                    if value.strip()
                ],
            }
            if not lombardi:
                mobility["surface"].update({
                    "theta_electron_m_per_V": args.surface_theta_electron_m_per_v,
                    "theta_hole_m_per_V": args.surface_theta_hole_m_per_v,
                    "beta": args.surface_beta,
                    "min_factor": args.surface_min_factor,
                    "max_factor": 1.0,
                })
            write_json(config, runtime)
        if args.srh_doping_dependence:
            runtime = read_json(config)
            runtime["solver"]["srh_doping_dependence"] = {
                "enabled": True,
                "concentration_basis": "total_impurity",
                "electron": {
                    "tau_min_s": 0.0,
                    "tau_max_s": 1.0e-7,
                    "reference_doping_m3": 1.0e16,
                    "gamma": 1.0,
                },
                "hole": {
                    "tau_min_s": 0.0,
                    "tau_max_s": 1.0e-7,
                    "reference_doping_m3": 1.0e16,
                    "gamma": 1.0,
                },
            }
            write_json(config, runtime)
        if args.path_ionization_tracing_vector != "electric_field":
            runtime = read_json(config)
            runtime["sweep"]["diagnostics"]["path_ionization_integrals"][
                "tracing_vector"
            ] = args.path_ionization_tracing_vector
            write_json(config, runtime)
        if args.path_ionization_tracing_current_relative_floor != 1.0e-8:
            runtime = read_json(config)
            runtime["sweep"]["diagnostics"]["path_ionization_integrals"][
                "tracing_current_relative_floor"
            ] = args.path_ionization_tracing_current_relative_floor
            write_json(config, runtime)
        if (
            args.path_ionization_tracing_vector == "sentaurus_eparallel_adaptive"
            or args.path_ionization_tracing_qf_relative_floor != 5.1e-3
        ):
            runtime = read_json(config)
            runtime["sweep"]["diagnostics"]["path_ionization_integrals"][
                "tracing_qf_relative_floor"
            ] = args.path_ionization_tracing_qf_relative_floor
            write_json(config, runtime)
        if args.path_ionization_tracing_direction != "bidirectional":
            runtime = read_json(config)
            runtime["sweep"]["diagnostics"]["path_ionization_integrals"][
                "tracing_direction"
            ] = args.path_ionization_tracing_direction
            write_json(config, runtime)
        if (
            args.path_ionization_integrals
            or args.path_ionization_electron_stop_field_v_per_m != 0.0
            or args.path_ionization_hole_stop_field_v_per_m != 0.0
            or args.path_ionization_mean_definition
                != "carrier_integral_arithmetic"
            or args.path_ionization_break_ordering != "path_mean"
            or args.path_ionization_path_retention
                != "numbered_peak_groups"
            or args.path_ionization_seed_mode != "nodal_local_maxima"
        ):
            runtime = read_json(config)
            path_diag = runtime["sweep"]["diagnostics"][
                "path_ionization_integrals"
            ]
            path_diag["electron_stop_field_V_per_m"] = (
                args.path_ionization_electron_stop_field_v_per_m
            )
            path_diag["hole_stop_field_V_per_m"] = (
                args.path_ionization_hole_stop_field_v_per_m
            )
            path_diag["mean_definition"] = args.path_ionization_mean_definition
            path_diag["break_ordering"] = args.path_ionization_break_ordering
            path_diag["seed_mode"] = args.path_ionization_seed_mode
            path_diag["path_retention"] = (
                args.path_ionization_path_retention
            )
            write_json(config, runtime)
        print(f":: running {mode}", flush=True)
        completed = subprocess.run(
            [str(runner), "--config", str(config)], cwd=case_dir, check=False
        )
        summarize_edges(case_dir)
        if completed.returncode != 0:
            overall = completed.returncode
            break
    return overall


if __name__ == "__main__":
    raise SystemExit(main())
