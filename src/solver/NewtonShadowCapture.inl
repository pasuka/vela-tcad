// Included inside solveClassicalWithFrozenElectronQuantumPotential.
// Experiment-only environment switch; no calibrated convergence policy.
const bool shadowCapture = [] {
    const char* value = std::getenv("VELA_P1_SHADOW_CAPTURE");
    return value && std::string(value) == "1";
}();
nlohmann::json shadowLinear = vela::shadow_convergence::unavailable("initial state has no solve");
nlohmann::json shadowSolvedLinear = vela::shadow_convergence::unavailable("initial state has no solve");
const auto appendShadow = [&](nlohmann::json& doc, int iteration,
    const VectorXd& raw, const VectorXd& processed, double alpha,
    const NewtonIterationInfo& info, const NewtonGlobalContinuityClosureEvaluation& closure) {
    if (!shadowCapture) return;
    using namespace vela::shadow_convergence;
    Json s = {{"schema", "vela.p1_shadow.capture.v1"}, {"affects_acceptance", false},
        {"iteration", iteration}, {"thermal_voltage_V", Vt},
        {"potential_scale_V", potentialScale}, {"statistics", cfg_.carrierStatistics.model},
        {"temperature_K", cfg_.temperature_K}, {"alpha", iteration > 0
            ? uncalibrated(alpha) : unavailable("initial state has no update")},
        {"s0_physical_direction", unavailable("independent unmodified physical-system direction not captured")},
        {"s1_identity", "configured solve after exact-Dirichlet update; before cap/restoration/damping"},
        {"s2_identity", "QF-coordinate finite map is identity for this frozen isothermal branch"},
        {"s3_identity", "alpha times processed direction; split writeback is in adjacent HDF5 states"},
        {"physical_residual_identity", "pre_closure_replacement; post_constraint_replacement"},
        {"residual_internal", vector(r)}, {"linear_quality", shadowLinear},
        {"solved_equation_quality", shadowSolvedLinear},
        {"fixed_reference_scales", {{"status", "not_calibrated"},
            {"basis", "unit-valued reporting references; not calibrated physical accuracy tolerances"},
            {"poisson_C_per_m", 1.0}, {"continuity_A_per_um", 1.0}, {"constraint_V", 1.0}}},
        {"shadow_numerically_converged", {{"status", "not_calibrated"}}},
        {"shadow_device_audit_passed", unavailable("device tolerances and independent reference not bound")},
        {"shadow_continuation_eligible", {{"status", "not_calibrated"}}}};
    const auto gate = [](bool passed) { return Json{{"status", passed ? "pass" : "fail"}}; };
    s["legacy_row_gate_passed"] = gate(carrierRowsAcceptConvergence(info.carrierRowConvergence));
    s["retained_guards"]["global_continuity_closure"] = gate(globalClosureAcceptsConvergence(closure));
    s["retained_guards"]["finite_packed_and_residual"] = gate(x.allFinite() && r.allFinite());
    s["retained_guards"]["contact_majority_qf_branch"] = gate(
        cfg_.contactMajorityQfBranchDropLimit_V <= 0.0 ||
        (std::isfinite(info.contactMajorityQfDrop) &&
         info.contactMajorityQfDrop <= cfg_.contactMajorityQfBranchDropLimit_V));
    s["retained_guards"]["contact_majority_qf_branch"]["drop_V"] = info.contactMajorityQfDrop;
    s["retained_guards"]["contact_majority_qf_branch"]["limit_V"] = cfg_.contactMajorityQfBranchDropLimit_V;
    if (iteration == 0)
        s["retained_guards"]["contact_majority_qf_branch"] = unavailable("initial progress record does not bind a branch verdict");
    for (const char* key : {"population_positivity", "inverse_mapping", "boundary_identity",
                           "archive_consistency", "contact_floor_guard"})
        s["retained_guards"][key] = unavailable("no per-iteration original verdict available at this hook");
    s["legacy_full_qualification"] = unavailable("full success-path verdict is not implied by partial guards");
    const double currentFactor = constants::q * assembler.continuityResidualScale() *
        (scaling.enabled ? scaling.currentDensityLineIntegralFactor : 1.0) * 1e-6;
    s["continuity_residual_to_A_per_um"] = currentFactor;
    s["poisson_residual_to_C_per_m"] = scaling.enabled
        ? scaling.permittivityReference_F_per_m * scaling.V0 : 1.0;
    s["constraint_residual_to_V"] = potentialScale;
    s["continuity_units_source"] = "production evaluateCarrierTermDiagnostics conversion";
    const auto closureLedger = [&](const NewtonGlobalContinuityCarrierClosure& c) {
        return Json{{"contact_flux_internal", c.contactFlux},
            {"integrated_source_internal", c.integratedSource}, {"mismatch_internal", c.mismatch},
            {"ratio", c.ratio}, {"qualified", c.qualified}, {"satisfied", c.satisfied},
            {"identity", "original global closure evaluation; not a regional conservation proof"}};
    };
    s["global_ledger"] = {{"enabled", closure.enabled}, {"enforced", closure.enforced},
        {"electron", closureLedger(closure.electron)}, {"hole", closureLedger(closure.hole)}};
    if (iteration > 0) {
        s["s1_direction_internal"] = vector(raw);
        s["s3_direction_before_alpha_internal"] = vector(processed);
    } else {
        s["E_adv"] = unavailable("initial state has no direction");
        s["E_act"] = unavailable("initial state has no previous update");
    }
    s["constrained_indices"] = Json::array();
    for (const auto& [node, value] : bcs.psi) s["constrained_indices"].push_back(node);
    for (const auto& [node, value] : bcs.phin) s["constrained_indices"].push_back(N + node);
    for (const auto& [node, value] : bcs.phip) s["constrained_indices"].push_back(2 * N + node);
    std::vector<std::string> states{"not_calibrated"};
    for (const auto& entry : s["retained_guards"].items())
        states.push_back(entry.value().at("status").get<std::string>());
    s["shadow_numerically_converged"]["status"] = combine(states);
    s["shadow_continuation_eligible"]["status"] =
        combine({s["shadow_numerically_converged"]["status"].get<std::string>(), "not_evaluated"});
    doc["p1_shadow"] = std::move(s);
};
