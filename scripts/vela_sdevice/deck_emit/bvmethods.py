"""Explicit numerical profile for the qualified BVMethods Tri3 operator."""
from __future__ import annotations

STATE_FIELDS = {"ElectrostaticPotential", "eDensity", "hDensity", "SpaceCharge", "DopingConcentration"}
GEOMETRY = {"node_volume_policy": "barycentric", "require_non_obtuse": False,
            "poisson_permittivity_policy": "cell_material", "cell_box_policy": "delaunay_transfer"}


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
    solver["quasi_fermi_reference"] = "none"
    solver["continuity_row_scaling"] = {"enabled": False}
    # Do not inherit the legacy warm-start emitter's implicit 0.1 V cap.
    # Positive values remain explicit controls for reproducible cap comparisons.
    for key in ("quasi_fermi_update_limit_V", "quasi_fermi_update_limit_minority_V"):
        solver[key] = requested.get(key, 0.0)
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
        iterations = numeric.get("iterations") or ir["math"].get("iterations", 50)
        iterations = max(iterations, solver.get("max_iter", iterations), circuit.get("max_iterations", iterations))
        if damping != 1 or not_damped < iterations:
            raise ValueError("BVMethods CMD currently supports undamped Sentaurus settings only: LineSearchDamping=1 and NotDamped>=the effective iteration budget; vendor line-search/Bank-Rose algorithms are not implemented")
        circuit.setdefault("coupled_line_search_mode", "full_step")
        circuit.setdefault("solver", "coupled_newton")
        if circuit["solver"] != "coupled_newton":
            raise ValueError("BVMethods CMD profile requires coupled_newton")
        for target, source in (("coupled_initial_outer_step_V", "initial_step"),
                               ("coupled_min_outer_step_V", "min_step"),
                               ("coupled_max_outer_step_V", "max_step"),
                               ("coupled_outer_growth_factor", "growth_factor"),
                               ("coupled_outer_shrink_factor", "shrink_factor")):
            circuit.setdefault(target, sweep[source])
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
            "coordinates": "zero_reference", "geometry": GEOMETRY,
            "scope": "explicit numerical profile; cmd/par retain physical laws and coefficients"}
