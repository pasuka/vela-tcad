"""Explicit numerical profile for the qualified BVMethods Tri3 operator."""
from __future__ import annotations
import math

STATE_FIELDS = {"ElectrostaticPotential", "eDensity", "hDensity", "SpaceCharge", "DopingConcentration"}
GEOMETRY = {"node_volume_policy": "barycentric", "require_non_obtuse": False,
            "poisson_permittivity_policy": "cell_material", "cell_box_policy": "delaunay_transfer"}


# These numerical coordinates do not exist in a concentration-primary deck.
POPULATION_QF_DEFAULTS = {
    "quasi_fermi_reference": "none", "quasi_fermi_update_limit_V": 0.,
    "quasi_fermi_update_limit_minority_V": 0., "max_update": 0.,
    "quasi_fermi_update_limit_mode": "componentwise_poisson_recorrect",
    "quasi_fermi_trust_region_growth_factor": 1., "quasi_fermi_trust_region_max_multiplier": 1.,
    "quasi_fermi_trust_region_expansion_threshold": .75,
    "quasi_fermi_trust_region_shrink_factor": 1., "quasi_fermi_trust_region_min_multiplier": 1.,
    "quasi_fermi_recenter_on_initial_state": False, "quasi_fermi_recenter_on_stall": False,
    "contact_majority_qf_branch_drop_limit_V": 0.,
}


def apply_bvmethods_profile(deck, ir, stage, sidecar, *, source_stage=None):
    solver, sweep = deck["solver"], deck["sweep"]
    mobility, avalanche = solver["mobility"], solver["impact_ionization"]
    if (mobility["model"] != "masetti_field_lombardi" or
            mobility.get("high_field_driving_force") != "quasi_fermi_gradient" or
            solver.get("carrier_statistics", {}).get("model") != "fermi_dirac" or
            solver.get("bandgap_narrowing", {}).get("model") != "old_slotboom" or
            avalanche["model"] != "van_overstraeten" or
            avalanche.get("driving_force") != "quasi_fermi_gradient" or
            solver.get("band_to_band", {}).get("model", "none") != "none" or
            set(solver.get("recombination", [])) - {"srh"} or
            solver["temperature_K"] != 300 or not ir["math"]["aval_derivatives"]):
        raise ValueError("BVMethods profile requires classical 300 K Fermi/OldSlotboom, Masetti/Enormal/GradQF, SRH and GradQF VanOverstraeten with AvalDerivatives; BTBT is unsupported")
    requested = sidecar.get("solver", {})
    if requested.get("quasi_fermi_reference", "none") != "none" or requested.get("continuity_row_scaling", {}).get("enabled", False):
        raise ValueError("BVMethods profile requires zero-reference QF coordinates without legacy row scaling")
    if sidecar.get("avalanche_current_support_profile"):
        raise ValueError("BVMethods and PN2D avalanche discretization profiles cannot be combined")
    for key, value in GEOMETRY.items():
        if key in sidecar.get("mesh_geometry", {}) and sidecar["mesh_geometry"][key] != value:
            raise ValueError("partial BVMethods geometry override")
    if not sweep["vtk_fields"] or set(sweep["vtk_fields"]) - STATE_FIELDS:
        raise ValueError("BVMethods CMD profile currently requires explicit basic-state Plot fields; current/source Plot fields are not qualified")
    # Fermi remains enabled. Select the explicit numerical approximation already
    # qualified with this operator, retaining all cmd/par model coefficients.
    solver["carrier_statistics"]["model"] = "fermi_dirac_joyce_dixon"
    solver["coupled_operator"] = "bvmethods_tri3_qf"
    variables = requested.get("state_variables", "potential_density")
    if variables not in ("potential_quasi_fermi", "potential_density"):
        raise ValueError("unknown state_variables")
    sg = requested.get("population_sg_current", "quasi_fermi")
    if sg not in ("quasi_fermi", "band_population") or (sg != "quasi_fermi" and variables != "potential_density"):
        raise ValueError("invalid population_sg_current")
    statistical = requested.get("population_sg_bootstrap_linearization", "unchanged")
    if statistical not in ("unchanged", "floor_adjacent_edges"):
        raise ValueError("invalid population_sg_bootstrap_linearization")
    limits = requested.get("population_update", {})
    if not isinstance(limits, dict):
        raise ValueError("population_update must be an object")
    limits = dict(limits)
    if not {"minimum_density_cm3", "minimum_density_m3"} & limits.keys():
        limits["minimum_density_cm3"] = 1e-100
    if statistical == "floor_adjacent_edges":
        values = [limits[k] for k in ("minimum_density_cm3", "minimum_density_m3") if k in limits]
        if (variables != "potential_density" or sg != "quasi_fermi" or len(values) != 1 or
                isinstance(values[0], bool) or not isinstance(values[0], (int, float)) or
                not math.isfinite(values[0]) or values[0] <= 0):
            raise ValueError("statistical bootstrap requires population quasi_fermi SG and a positive density floor")
    if "population_sg_bootstrap_linearization" in requested:
        solver["population_sg_bootstrap_linearization"] = statistical
    recovery = requested.get("population_derivative_recovery", "after_protection" if variables == "potential_density" else "full")
    if recovery not in ("full", "automatic", "after_protection") or (recovery != "full" and variables != "potential_density"):
        raise ValueError("invalid population_derivative_recovery")
    if variables == "potential_density":
        if requested.get("method", "newton") not in ("newton", "poisson_only"):
            raise ValueError("potential_density does not support Gummel, hybrid or frozen-state solver paths")
        solver["state_variables"] = variables
        if "population_sg_current" in requested:
            solver["population_sg_current"] = sg
        solver["population_derivative_recovery"] = recovery
        solver["primary_variable_coordinates"] = requested.get("primary_variable_coordinates", "unit")
        solver["population_stopping_rule"] = requested.get("population_stopping_rule", "native_approx")
        if solver["primary_variable_coordinates"] not in ("unit", "physical"):
            raise ValueError("invalid primary_variable_coordinates")
        if solver["population_stopping_rule"] not in ("legacy", "native_approx"):
            raise ValueError("invalid population_stopping_rule")
        solver["population_update"] = {"protection": "local_nonpositive",
            "local_retained_fraction": .01, **limits}
        for key, default in POPULATION_QF_DEFAULTS.items():
            if requested.get(key, default) != default:
                raise ValueError(f"potential_density does not support {key}")
        for key in ("stable_merit_comparison", "carrier_row_qualified_stall_acceptance"):
            if requested.get(key, False):
                raise ValueError(f"potential_density does not support {key}")

    solver["quasi_fermi_reference"] = "none"
    solver["continuity_row_scaling"] = {"enabled": False}
    # Do not inherit the legacy warm-start emitter's implicit 0.1 V cap.
    # Positive values remain explicit controls for reproducible cap comparisons.
    for key in ("quasi_fermi_update_limit_V", "quasi_fermi_update_limit_minority_V"):
        solver[key] = requested.get(key, 0.0)
    if variables == "potential_density":
        # Remove inherited inactive QF controls instead of presenting them as
        # population update or stopping policies. Explicit active values fail above.
        for key in POPULATION_QF_DEFAULTS:
            solver.pop(key, None)
    for key in ("carrier_row_convergence", "global_continuity_closure"):
        if key in requested:
            # Preserve the explicitly selected gate contract, rather than
            # retaining unrelated legacy warm-start defaults underneath it.
            solver[key] = dict(requested[key])
    if solver["method"] not in ("poisson_only", "frozen_state"):
        solver.update(method="newton", warm_start=True)
    avalanche.update(coupling_mode="self_consistent", source_jacobian="local_ad")
    solver["band_to_band"] = {"model": "none"}
    for key, value in (("coefficient_eV", .009), ("reference_doping_cm3", 1e17), ("smoothing", .5)):
        solver["bandgap_narrowing"].setdefault(key, value)
    solver["region_resolved_interface_assembly"] = {"poisson_charge_node_volume": "signed_transport"}
    deck["mesh_geometry"].update(GEOMETRY)
    if solver["method"] == "newton":
        sweep.setdefault("step_growth_mode", "newton_iterations")
    sweep.setdefault("stop_on_failure", True)
    numerical_circuit = sidecar.get("sweep", {}).get("external_circuit", {})
    if "coupled_max_step_retries" in numerical_circuit:
        sweep["max_retries"] = numerical_circuit["coupled_max_step_retries"]
    if circuit := sweep.get("external_circuit"):
        # Sentaurus settings are resolved per Coupled stage, overriding Math.
        numeric = source_stage if source_stage is not None else stage
        damping = numeric.get("line_search_damping", ir["math"].get("line_search_damping", 1.0))
        not_damped = numeric.get("not_damped", ir["math"].get("not_damped", 1000))
        # Resolve the budget at the actual Coupled stage. An unrelated seed
        # solve may have a larger solver.max_iter and must not enlarge the
        # resistor sweep budget or its adaptive step-growth denominator.
        iterations = numeric.get("iterations")
        if iterations is None:
            iterations = ir["math"].get("iterations", 50)
        circuit.setdefault("max_iterations", iterations)
        iterations = circuit["max_iterations"]
        if damping != 1 or not_damped < iterations:
            raise ValueError("BVMethods CMD currently supports undamped Sentaurus settings only: LineSearchDamping=1 and NotDamped>=the effective iteration budget; vendor line-search/Bank-Rose algorithms are not implemented")
        circuit.setdefault("coupled_line_search_mode", "full_step")
        circuit.setdefault("solver", "coupled_newton")
        if variables == "potential_density" and solver["population_stopping_rule"] == "native_approx":
            if not ir["math"].get("rel_err_control", True):
                raise ValueError("native_approx supports RelErrControl only; choose legacy stopping for absolute update control")
            metric = circuit.setdefault("coupled_update_convergence", {})
            metric.setdefault("mode", "report")
            metric.setdefault("digits", ir["math"].get("digits", 5))
            for key, source in (("electron_reference_cm3", "err_ref_electron_cm3"),
                                ("hole_reference_cm3", "err_ref_hole_cm3")):
                if ir["math"].get(source) is not None:
                    metric.setdefault(key, ir["math"][source])
        if circuit["solver"] != "coupled_newton":
            raise ValueError("BVMethods CMD profile requires coupled_newton")
        for target, source in (("coupled_initial_outer_step_V", "initial_step"),
                               ("coupled_min_outer_step_V", "min_step"),
                               ("coupled_max_outer_step_V", "max_step"),
                               ("coupled_outer_growth_factor", "growth_factor"),
                               ("coupled_outer_shrink_factor", "shrink_factor")):
            circuit.setdefault(target, sweep[source])
        mode = circuit.setdefault("coupled_step_control", "quasistationary")
        if mode not in ("quasistationary", "legacy_inner_guard"):
            raise ValueError("unsupported external circuit coupled_step_control")
        if mode == "quasistationary":
            extra = circuit.get("coupled_continuation", {})
            if ("max_inner_voltage_step_V" in circuit or extra.get("enabled", False)
                    or "inner_step_control" in extra or extra.get("retry_from_accepted", False)):
                raise ValueError("quasistationary does not accept inner voltage guards or retry_from_accepted; select legacy_inner_guard explicitly")
        else:
            circuit.setdefault("coupled_continuation", {"enabled": True, "retry_from_accepted": True})
        sweep["step_growth_mode"] = "newton_iterations"
        boundary = sweep.get("boundary_control", {})
        boundary.pop("adaptive_outer_continuation", None)
        if not boundary:
            sweep.pop("boundary_control", None)
        # MaxStep is an adaptive cap. Explicit Plot times remain mandatory;
        # an absent output grid does not invent additional voltage targets.
        if stage.get("sweep", {}).get("mode") == "voltage" and not stage["required_t"]:
            delta = sweep["stop"] - sweep["start"]
            sweep["step"] = delta
            stage["sweep"]["step_V"] = delta
    return {"operator": "bvmethods_tri3_qf", "statistics_approximation": "fermi_dirac_joyce_dixon",
            "coordinates": "potential_density" if variables == "potential_density" else "zero_reference", "geometry": GEOMETRY,
            "scope": "explicit numerical profile; cmd/par retain physical laws and coefficients"}
