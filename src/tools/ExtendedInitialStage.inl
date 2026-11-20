// Included after the normal dc_sweep problem loaders.
bool isExtendedSlotStage(const nlohmann::json& deck)
{
    if (!deck.contains("solver") || !deck.contains("sweep"))
        return false;
    const auto& solver = deck.at("solver");
    const auto& sweep = deck.at("sweep");
    return deck.value("simulation_type", "") == "dc_sweep"
        && solver.value("arithmetic_precision", "") == "extended80"
        && solver.value("state_variables", "") == "potential_density"
        && solver.value("coupled_operator", "") == "legacy"
        && solver.value("primary_variable_coordinates", "") == "physical"
        && solver.value("population_stopping_rule", "") == "native_approx"
        && solver.value("population_derivative_recovery", "full") == "full"
        && solver.value("population_writeback", "independent") == "independent"
        && solver.value("carrier_statistics", nlohmann::json::object()).value("model", "") == "boltzmann"
        && solver.value("mobility", nlohmann::json::object()).value("model", "") == "ialmob"
        && deck.contains("decimal_stage_controls") && sweep.value("mode", "") == "iv"
        && !sweep.contains("voltage_to_current") && !sweep.contains("continuation")
        && !sweep.contains("boundary_control");
}

std::string extendedStageDecimal(long double value)
{
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::setprecision(std::numeric_limits<long double>::max_digits10) << value;
    return text.str();
}

std::string extendedStageFileHash(const std::filesystem::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
        throw std::runtime_error("Cannot hash extended stage file: " + path.string());
    return vela::stateSha256(std::string(std::istreambuf_iterator<char>(stream), {}));
}

nlohmann::json runExtendedSlotStage(const std::string& configFile, const nlohmann::json& deck, bool publish)
{
    using namespace vela;
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::invalid_argument("SLOT stage requires the 80-bit runner");
    using WideJson = nlohmann::basic_json<std::map, std::vector, std::string, bool,
                                         std::int64_t, std::uint64_t, long double>;
    std::ifstream raw(configFile);
    if (!raw)
        throw std::invalid_argument("Cannot open the SLOT stage configuration");
    const auto wideDeck = WideJson::parse(raw);
    const auto& rawInput = wideDeck.at("solver");
    const auto& input = deck.at("solver");
    const auto& sweep = deck.at("sweep");
    if (!isExtendedSlotStage(deck))
        throw std::invalid_argument("Unsupported extended80 SLOT stage profile");
    const auto mode = input.at("method").get<std::string>();
    const bool coldStage = mode == "newton" && sweep.contains("initialization");
    const bool loadOnly = mode == "frozen_state";
    const bool circuitStage = mode == "newton" && sweep.contains("external_circuit");
    if (publish) {
        const char* backend = std::getenv("VELA_LINEAR_SOLVER");
        if (backend && *backend && std::string_view(backend) != "sparse_lu")
            throw std::invalid_argument("Extended SLOT stage requires the long-double SparseLU backend");
        if (input.value("compensated_population_state", false) ||
            input.value("population_sg_current", "quasi_fermi") != "quasi_fermi" ||
            input.value("population_sg_bootstrap_linearization", "unchanged") != "unchanged")
            throw std::invalid_argument("Unsupported extended SLOT population settings");
    }
    const auto start = stageControl(deck, "sweep", "start", sweep.at("start").get<double>());
    const auto stop = stageControl(deck, "sweep", "stop", sweep.at("stop").get<double>());
    if (coldStage ? (start != 0 || stop != 0 || sweep.contains("initial_state_file") ||
                     sweep.contains("external_circuit") ||
                     sweep.at("initialization").at("mode") != "poisson_block")
                  : (!publish || !sweep.contains("initial_state_file") || !(loadOnly || circuitStage) ||
                     (loadOnly && (start != stop || sweep.contains("external_circuit")))))
        throw std::invalid_argument("Invalid extended80 SLOT stage ordering or starting state");

    const auto directory = configDirectory(configFile);
    const auto inputHashes = stateInputProvenance(deck, directory);
    const auto deckHash = extendedStageFileHash(configFile);
    std::vector<long double> externalCouples;
    auto problem = loadNewtonProblem(configFile, deck, &externalCouples);
    auto& config = problem.newton;
    const auto& scaling = config.inputScaling;
    const auto& refs = config.unitScalingRefs;
    if (!scaling.isUnitScaling() || !refs.decimalTemperature_K ||
        !deck.contains("decimal_stage_controls") || !problem.fixedCharges.empty() ||
        !problem.sheetCharges.empty() || !parseBoundarySegmentSpecs(deck).empty() ||
        config.mobility.model != "ialmob" ||
        config.impactIonization.model != "van_overstraeten" ||
        config.electronQuantumPotential.enabled || config.splitDDState ||
        config.recombinationVolumePolicy != "material_local_box" ||
        config.carrierRowConvergence.mode == "enforce" ||
        config.globalContinuityClosure.mode == "enforce" ||
        config.contactMajorityQfBranchDropLimit_V > 0)
        throw std::invalid_argument("Unsupported physics or acceptance options in extended80 initial stage");
    const auto temperature = refs.decimalTemperature_K->kelvin;
    const auto vt = arithmeticThermalVoltageForInputs(config.temperature_K, refs);
    const auto& units = scaling.unitSystem();
    constexpr long double length = 1e-6L;
    constexpr long double unitConcentration_m3 = 1e6L;
    auto contacts = parseContactBoundarySpecs(deck);
    const auto currentContact = sweep.at("contact").get<std::string>();
    if (sweep.value("current_contact", currentContact) != currentContact)
        throw std::invalid_argument("Extended SLOT requires current and driven contact identity");
    for (auto& contact : contacts) {
        if (contact.name != currentContact || coldStage)
            continue;
        contact.bias = 0;
        contact.decimalInitialBias_V = 0;
    }
    const auto& solverText = deck.at("decimal_stage_controls").at("solver");
    const auto& sourceText = solverText.at("slot_source_coefficients");
    const auto coefficient = [&](const std::string& key, Real view) {
        return decimalStageControl(sourceText.at(key), key, view);
    };
    const auto narrowingCoefficient = [&](const std::string& key, Real compatibility, long double fallback) {
        const auto path = "bandgap_narrowing." + key;
        if (sourceText.contains(path))
            return coefficient(path, compatibility);
        if (static_cast<Real>(fallback) != compatibility)
            throw std::invalid_argument("Missing decimal OldSlotboom coefficient: " + key);
        return fallback;
    };
    if (config.bandgapNarrowing.model != "old_slotboom")
        throw std::invalid_argument("Extended SLOT stage requires OldSlotboom");
    const auto& narrowing = config.bandgapNarrowing;
    const ExtendedSlotboomCoefficients narrowingCoefficients{
        narrowingCoefficient("reference_doping_m3", narrowing.referenceDoping, 1e17L),
        narrowingCoefficient("coefficient_eV", narrowing.coefficient, 9.0e-3L),
        narrowingCoefficient("smoothing", narrowing.smoothing, 0.5L),
        narrowingCoefficient("offset_eV", narrowing.offset, 0.L)};
    const auto material = buildExtendedBoltzmannMaterials(
        problem.mesh, problem.matdb, problem.doping, narrowing, temperature, vt, narrowingCoefficients);
    const auto targets = buildExtendedZeroBiasOhmicTargets(problem.mesh, problem.doping, material,
        contacts, vt, config.contactBoundaryReconstruction == "dominant_signed_contact_mean");
    const auto port = std::ranges::find(problem.mesh.contacts(), currentContact, &Contact::name);
    if (port == problem.mesh.contacts().end() || port->node_ids.empty())
        throw std::invalid_argument("Unknown extended80 initial stage contact");

    PoissonScalingSpec po;
    po.enabled = true;
    po.potentialScale_V = static_cast<Real>(vt);
    Real epsr = 1;
    for (const auto& region : problem.mesh.regions())
        epsr = std::max(epsr, problem.matdb.getMaterial(region.material).eps_r);
    po.permittivityReference_F_per_m = constants::eps0 * epsr;
    po.chargeAreaFactor = units.chargeAreaFactor();
    po.chargeLineFactor = units.chargeLineFactor();
    po.boundaryLengthFactor = units.lengthMPerInternal();
    po.poissonChargeVolumePolicy = config.poissonChargeVolumePolicy;
    const auto box = parseBoxGeometryOptions(deck);
    if (box.cellBoxPolicy != BoxGeometryBuilder::CellBoxPolicy::AverageIntersection ||
        config.poissonChargeVolumePolicy != PoissonChargeVolumePolicy::MaterialLocalBox)
        throw std::invalid_argument("Extended80 initial stage requires material-local AverageBox");
    const auto permittivities = readExtendedPermittivities(
        resolvePath(configDirectory(configFile), deck.at("materials_file").get<std::string>()), problem.matdb);
    const auto geometry = buildExtendedAveragePoissonGeometry(
        problem.mesh, problem.matdb, temperature, permittivities);
    const auto normalizationInput = scaling.deckFormatVersion == 2
        ? rawInput.value("normalization", WideJson::object()) : wideDeck.at("scaling");
    std::optional<ExtendedPopulationNormalization> normalization;
    long double concentration_cm3 = 1.L;
    if (normalizationInput.contains("reference_concentration_cm3") ||
        normalizationInput.contains("characteristic_length_um") ||
        normalizationInput.contains("reference_mobility_cm2_V_s")) {
        Real niFloor = 1.;
        for (const auto& region : problem.mesh.regions())
            niFloor = std::max(niFloor, problem.matdb.getMaterial(region.material, config.temperature_K).ni);
        const auto inputs = UnitScalingSystem::autoInputsFrom(problem.mesh, problem.doping, problem.matdb, niFloor);
        const auto scales = UnitScalingSystem::fromInputs(config.temperature_K,
            static_cast<Real>(geometry.permittivityReference_F_per_m), inputs, refs, units);
        const auto& calculated = scales.arithmeticValues();
        const auto exact = [&](const char* name, long double fallback) {
            const auto value = normalizationInput.find(name);
            return value != normalizationInput.end() && value->is_number() ? value->get<long double>() : fallback;
        };
        concentration_cm3 = exact("reference_concentration_cm3", calculated.C0);
        const auto length_m = length * exact("characteristic_length_um", calculated.L0);
        const auto mobility = 1e-4L * exact("reference_mobility_cm2_V_s", calculated.mu0);
        const auto concentration_m3 = unitConcentration_m3 * concentration_cm3;
        normalization = ExtendedPopulationNormalization{concentration_cm3,
            geometry.permittivityReference_F_per_m * vt /
                (constants::arithmeticConstants<long double>().q * concentration_m3 * length_m * length_m),
            1.L / (concentration_m3 * mobility * vt)};
        for (auto scale : {concentration_m3, normalization->poissonRowFactor, normalization->continuityRowFactor})
            if (!(scale > 0) || !std::isfinite(scale))
                throw std::invalid_argument("Invalid extended SLOT normalization scales");
    }
    const auto concentrationScale_m3 = unitConcentration_m3 * concentration_cm3;
    const auto* equationScales = normalization ? &*normalization : nullptr;
    ArithmeticPoissonAssembler poisson(problem.mesh, problem.matdb, problem.doping, {}, {}, {}, po, vt);
    nlohmann::json out{
        {"arithmetic_precision", "extended80"}, {"scope", publish ? extendedSlotStageProfile : "initial_stage_check"},
        {"accepted_points", 0}, {"converged", false}, {"coupled_converged", false},
        {"potential_scale_V", vt}, {"nodes", problem.mesh.numNodes()}, {"poisson_geometry", "arithmetic_averagebox"},
        {"decimal_permittivities", permittivities.size()},
        {"arithmetic_transport_couples", std::ranges::count_if(externalCouples, [](auto v) { return v >= 0; })}};
    if (normalization)
        out["equation_normalization"] = {{"potential_scale_V", extendedStageDecimal(vt)},
            {"reference_concentration_cm3", extendedStageDecimal(concentration_cm3)},
            {"poisson_row_factor", extendedStageDecimal(normalization->poissonRowFactor)},
            {"continuity_row_factor", extendedStageDecimal(normalization->continuityRowFactor)}};
    ExtendedPopulationState state;
    // Stage iteration controls, output paths and applied drain bias may change.
    // The model and mesh identity must survive every accepted-state handoff.
    nlohmann::json identityInput{{"files", inputHashes}, {"port", currentContact},
                                 {"contacts", deck.at("contacts")}};
    for (auto& contact : identityInput.at("contacts"))
        contact.erase("bias");
    for (const auto key : {"scaling", "physical_constants_profile", "discretization", "mesh_geometry", "doping"})
        if (deck.contains(key)) identityInput[key] = deck.at(key);
    for (const auto key : {"mobility", "impact_ionization", "bandgap_narrowing", "temperature_K",
                           "recombination", "srh_doping_dependence", "auger_density_dependence",
                           "recombination_materials", "auger_with_generation", "auger_excess_product",
                           "contact_boundary_reconstruction"})
        if (input.contains(key)) identityInput["solver"][key] = input.at(key);
    identityInput["decimal_sources"] = deck.at("decimal_stage_controls").at("solver").at("slot_source_coefficients");
    identityInput["thermal_voltage_V"] = extendedStageDecimal(vt);
    if (normalization)
        identityInput["equation_normalization"] = out.at("equation_normalization");
    const ExtendedPopulationArchiveIdentity identity{extendedStateMeshIdentity(problem.mesh, length),
        stateSha256(identityInput.dump()), ExtendedPopulationArchiveRole::Accepted};
    if (coldStage) {
        auto cold = solveExtendedBoltzmannColdStart(poisson, problem.doping, material, targets,
            vt, unitConcentration_m3, sweep.at("initialization").at("poisson_max_iter").get<int>(), rawInput.value("abstol", 1e-10L), &geometry);
        out.update({{"poisson_converged", cold.converged}, {"poisson_updates", cold.updates},
                    {"poisson_residual", cold.residual}});
        if (!cold.converged)
            return out;
        state = std::move(cold.state);
        if (normalization) {
            state.coordinates.tail(2 * problem.mesh.numNodes()) /= concentration_cm3;
            state.concentrationScale_m3 = concentrationScale_m3;
        }
    } else {
        poisson.assemble(geometry);
        state = readExtendedPopulationArchive(resolvePath(directory, sweep.at("initial_state_file").get<std::string>()),
                                              problem.mesh.numNodes(), identity);
        if (state.supplyBias_V != start || state.potentialScale_V != vt ||
            state.concentrationScale_m3 != concentrationScale_m3)
            throw std::invalid_argument("Extended SLOT parent control or coordinate scale mismatch");
    }

    prepareIalTransportGeometry(config.mobility, problem.mesh);
    std::map<std::string, std::array<long double, 2>> mobilityByMaterial;
    for (const auto& [name, values] : config.mobility.ialmob->constantMaterialMobility) {
        const auto prefix = "ialmob.constant_material_mobility_cm2_per_Vs." + name + ".";
        mobilityByMaterial.emplace(name, std::array<long double, 2>{
            coefficient(prefix + "electron", values[0]), coefficient(prefix + "hole", values[1])});
    }
    for (const auto& axis : {config.mobility.ialmob->crystalX, config.mobility.ialmob->crystalY})
        if (std::ranges::any_of(axis, [](Real value) { return value != 0 && std::abs(value) != 1; }))
            throw std::invalid_argument("Extended SLOT requires axis-aligned crystal coordinates");
    auto measures = readExtendedIalTransportMeasures(config.mobility, problem.mesh, length, mobilityByMaterial);
    out["arithmetic_ial_surface_cells"] = measures.surface.size();
    out["arithmetic_ial_vertex_cells"] = measures.vertex.size();
    out["arithmetic_ial_avalanche_cells"] = measures.avalanche.size();
    auto sourceConfig = recombinationModelConfig(config.recombination, config.taun, config.taup,
                                                 config.srhDopingDependence);
    sourceConfig.activeMaterials = config.recombinationMaterials;
    sourceConfig.volumePolicy = config.recombinationVolumePolicy;
    sourceConfig.augerCn = config.augerCn;
    sourceConfig.augerCp = config.augerCp;
    sourceConfig.augerDensityDependence = config.augerDensityDependence;
    sourceConfig.augerWithGeneration = config.augerWithGeneration;
    sourceConfig.augerExcessProduct = config.augerExcessProduct;
    sourceConfig.bandToBand = config.bandToBand;
    const RecombinationModel sources(sourceConfig);
    if (sourceConfig.srhDopingDependence.temperatureDependence)
        throw std::invalid_argument("Extended SLOT SRH temperature law is not bound to decimal input");
    const auto lifetime = [&](const std::string& carrier, const SRHLifetimeParameters& p) {
        const auto prefix = "srh_doping_dependence." + carrier + ".";
        return ArithmeticBoltzmannSourceCoefficients::Lifetime{
            coefficient(prefix + "tau_min_s", p.tauMin), coefficient(prefix + "tau_max_s", p.tauMax),
            coefficient(prefix + "reference_doping_m3", p.referenceDoping), coefficient(prefix + "gamma", p.gamma)};
    };
    const auto& srh = sourceConfig.srhDopingDependence;
    const auto& auger = sourceConfig.augerDensityDependence;
    const ArithmeticBoltzmannSourceCoefficients sourceCoefficients{
        coefficient("taun", sourceConfig.taun), coefficient("taup", sourceConfig.taup),
        {lifetime("electron", srh.electron), lifetime("hole", srh.hole)},
        {coefficient("auger_cn_m6_per_s", sourceConfig.augerCn),
         coefficient("auger_cp_m6_per_s", sourceConfig.augerCp)},
        {coefficient("auger_density_dependence.electron.enhancement", auger.electronEnhancement),
         coefficient("auger_density_dependence.hole.enhancement", auger.holeEnhancement)},
        {coefficient("auger_density_dependence.electron.reference_density_m3", auger.electronReferenceDensity),
         coefficient("auger_density_dependence.hole.reference_density_m3", auger.holeReferenceDensity)}};
    out["arithmetic_recombination_coefficients"] = 16;
    const auto& impact = config.impactIonization;
    if (impact.parameterSet != "default" || impact.aScale != 1 || impact.bScale != 1)
        throw std::invalid_argument("Extended SLOT requires unmodified explicit Van Overstraeten coefficients");
    const auto impactCoefficient = [&](const std::string& key, Real value) {
        return coefficient("impact_ionization." + key, value);
    };
    const auto inverseLength = static_cast<long double>(units.inverseLengthMInvPerInternal());
    const auto field = static_cast<long double>(units.electricFieldVPerMPerInternal());
    const auto switchField = field * impactCoefficient("switch_field_V_m", impact.switchField);
    const auto curve = [&](std::string carrier, Real lowA, Real lowB, Real highA, Real highB) {
        return ArithmeticElementIonizationCurve{
            inverseLength * impactCoefficient(carrier + "_a_low_m_inv", lowA),
            field * impactCoefficient(carrier + "_b_low_V_m", lowB),
            inverseLength * impactCoefficient(carrier + "_a_high_m_inv", highA),
            field * impactCoefficient(carrier + "_b_high_V_m", highB), switchField};
    };
    measures.ionization = ArithmeticIalIonizationCoefficients{
        curve("electron", impact.electronALow, impact.electronBLow, impact.electronAHigh, impact.electronBHigh),
        curve("hole", impact.holeALow, impact.holeBLow, impact.holeAHigh, impact.holeBHigh),
        field * impactCoefficient("minimum_field_V_m", impact.minimumField)};
    out["arithmetic_ionization_coefficients"] = 10;
    const auto cellMaterials = detail::buildCellMaterials(problem.mesh, problem.matdb, config.temperature_K);
    const auto included = detail::recombinationNodeMask(problem.mesh, cellMaterials, sourceConfig.activeMaterials);
    std::vector<long double> areas(problem.mesh.numNodes());
    for (Index node = 0; node < problem.mesh.numNodes(); ++node)
        if (included[node])
            areas[node] = geometry.chargeVolumes[node] * length * length;
    const auto avalanche = detail::ialElementAvalancheOptionsFromConfig(config.impactIonization);
    const auto withUnitPopulations = [&](const ExtendedPopulationState& point, const auto& evaluate) {
        if (!normalization)
            return evaluate(point);
        auto physical = point;
        physical.coordinates.tail(2 * problem.mesh.numNodes()) *= concentration_cm3;
        physical.concentrationScale_m3 = unitConcentration_m3;
        return evaluate(physical);
    };
    const auto continuity = [&](const ExtendedPopulationState& point, ExtendedContactCurrents* allContacts = nullptr) {
        return assembleExtendedIalPopulationContinuity(config.mobility, problem.mesh, problem.doping, point,
            material.electronIntrinsic, material.holeIntrinsic, length, unitConcentration_m3, vt, temperature,
            config.bandgapNarrowing.model != "none", config.mobility.ialmobTransportFeedbackDerivatives,
            avalanche, port->node_ids, ExtendedRecombinationSupport{sources, areas, &sourceCoefficients}, externalCouples, &measures, allContacts);
    };
    const auto assemble = [&](const ExtendedPopulationState& point) {
        // At fixed bias all potential constraints are independent.
        return withUnitPopulations(point, [&](const auto& physical) {
            return assembleExtendedPopulationDevice(poisson, physical, material, targets, {},
                continuity(physical), &geometry, equationScales);
        });
    };
    const auto& compatibility = config.populationUpdateConvergence;
    const auto rawMetric = rawInput.value("population_update_convergence", WideJson::object());
    ExtendedFixedBiasUpdateReferences metric{compatibility.digits,
        rawMetric.value("potential_reference_V", 0.025852L),
        rawMetric.value("electron_reference_cm3", 1e10L), rawMetric.value("hole_reference_cm3", 1e10L)};
    if (const auto found = solverText.find("population_update_convergence"); found != solverText.end()) {
        const auto read = [&](const char* key, double view, long double raw) {
            return found->contains(key) ? decimalStageControl(found->at(key), key, view) : raw;
        };
        metric.potential_V = read("potential_reference_V", compatibility.potentialReference_V, metric.potential_V);
        metric.electron_cm3 = read("electron_reference_cm3", compatibility.electronReference_cm3, metric.electron_cm3);
        metric.hole_cm3 = read("hole_reference_cm3", compatibility.holeReference_cm3, metric.hole_cm3);
    }
    ExtendedFixedBiasNewtonConfig corrector{
        config.maxIter, CoupledLoadLineConfig{}.maxLineSearchSteps, config.linearRefinementIterations,
        rawInput.value("damping_factor", 1.L), 1e-10L, 1e-8L, metric,
        config.populationUpdate, config.lineSearch, {}};
    const auto rawProtection = rawInput.value("population_update", WideJson::object());
    corrector.exactProtection = BasicPopulationTrialLimits<long double>{
        rawProtection.value("local_retained_fraction", 0.01L), rawProtection.value("fraction_tau", 0.99L),
        rawProtection.value("minimum_density_cm3", rawProtection.value("minimum_density_m3", 0.L)) / concentration_cm3};
    ExtendedFixedBiasNewtonResult result;
    ExtendedContactCurrents currents;
    using Row = std::vector<std::string>;
    std::vector<Row> rows;
    const auto outputArchive = sweep.value("write_state_file", "");
    const auto outputCsv = deck.value("output_csv", "");
    if (publish && (outputArchive.empty() || outputCsv.empty() ||
        resolvePath(directory, outputArchive) == resolvePath(directory, outputCsv)))
        throw std::invalid_argument("Extended SLOT requires distinct CSV and accepted archive outputs");
    const auto& print = extendedStageDecimal;
    const auto& sweepText = deck.at("decimal_stage_controls").at("sweep");
    std::vector<long double> checkpoints;
    if (const auto exact = sweepText.find("vtk_bias_points"); exact != sweepText.end()) {
        const auto& view = sweep.at("vtk_bias_points");
        if (view.size() != exact->size())
            throw std::invalid_argument("Extended SLOT checkpoint number mismatch");
        for (std::size_t k = 0; k < exact->size(); ++k)
            checkpoints.push_back(decimalStageControl(exact->at(k), "vtk_bias_points", view.at(k).get<double>()));
    }
    const auto writePlot = [&] {
        if (!publish || !sweep.value("write_vtk", false))
            return;
        if (!checkpoints.empty() && state.supplyBias_V != stop &&
            std::ranges::find(checkpoints, state.supplyBias_V) == checkpoints.end())
            return;
        VtuWriter plot(problem.mesh, sweep.value("vtk_binary", false) ? VtuWriter::Encoding::Binary
                                                                      : VtuWriter::Encoding::Ascii);
        const auto nodes = problem.mesh.numNodes();
        TdrFieldSelection fields(sweep.value("vtk_fields", std::vector<std::string>{}));
        const auto add = [&](const char* name, const char* unit, const auto& value) {
            if (!fields.wantsScalar(name))
                return;
            std::vector<Real> display(nodes);
            for (Index i = 0; i < nodes; ++i)
                display[i] = static_cast<Real>(value(i));
            plot.addPointScalar(name, unit, display);
        };
        add("ElectrostaticPotential", "V", [&](Index i) { return vt * state.coordinates[i]; });
        add("eDensity", "cm^-3", [&](Index i) { return state.coordinates[nodes + i] * concentration_cm3; });
        add("hDensity", "cm^-3", [&](Index i) { return state.coordinates[2 * nodes + i] * concentration_cm3; });
        add("DonorConcentration", "cm^-3", [&](Index i) { return problem.doping.dopingOperandsForAssembly(i)[0]; });
        add("AcceptorConcentration", "cm^-3", [&](Index i) { return problem.doping.dopingOperandsForAssembly(i)[1]; });
        add("DopingConcentration", "cm^-3", [&](Index i) { return problem.doping.netDopingForAssembly(i); });
        add("EffectiveIntrinsicDensity", "cm^-3", [&](Index i) { return material.intrinsic[i]; });
        add("eQuasiFermiPotential", "V", [&](Index i) {
            const auto n = state.coordinates[nodes + i] * concentration_cm3;
            return n > 0 && material.electronIntrinsic[i] > 0
                ? vt * (state.coordinates[i] - std::log(n / material.electronIntrinsic[i])) : 0.L;
        });
        add("hQuasiFermiPotential", "V", [&](Index i) {
            const auto p = state.coordinates[2 * nodes + i] * concentration_cm3;
            return p > 0 && material.holeIntrinsic[i] > 0
                ? vt * (state.coordinates[i] + std::log(p / material.holeIntrinsic[i])) : 0.L;
        });
        const auto arrays = plot.pointArrayNames();
        std::vector<std::string> skipped;
        for (const auto& field : sweep.value("vtk_fields", std::vector<std::string>{}))
            if (std::ranges::find(arrays, field) == arrays.end())
                skipped.push_back(field);
        out["vtk_fields_skipped"] = skipped;
        if (const auto metadata = sweep.find("vtk_field_data"); metadata != sweep.end()) {
            for (const auto& [name, value] : metadata->items())
                if (value.is_string()) plot.addFieldString(name, value.get<std::string>());
                else if (value.is_number()) plot.addFieldValue(name, value.get<Real>());
        }
        plot.addFieldValue("PointIndex", static_cast<Real>(rows.size() - 1));
        plot.addFieldValue("t", stop == start ? 0 : static_cast<Real>((state.supplyBias_V - start) / (stop - start)));
        plot.addFieldValue(currentContact + "_InnerVoltage_V", static_cast<Real>(state.deviceBias_V));
        plot.addFieldValue(currentContact + "_OuterVoltage_V", static_cast<Real>(state.supplyBias_V));
        std::ostringstream name;
        name << sweep.value("vtk_prefix", "stage00") << '_' << std::setw(6) << std::setfill('0') << rows.size() << ".vtu";
        plot.write(resolvePath(directory, name.str()));
    };
    long double loadResidual = 0;
    const auto recordAccepted = [&] {
        const auto value = withUnitPopulations(state, [&](const auto& physical) {
            return continuity(physical, &currents).current_A_per_um;
        });
        if (value != currents.at(currentContact))
            throw std::runtime_error("Inconsistent extended SLOT terminal current");
        state.drainCurrent_A_per_um = value[0] - value[1];
        Row row{print(state.supplyBias_V), print(state.supplyBias_V), print(state.deviceBias_V),
            currentContact, print(value[0]), print(-value[1]), print(state.drainCurrent_A_per_um),
            "1", std::to_string(result.updates), print(loadResidual)};
        for (const auto& [name, pair] : currents) {
            row.push_back(print(pair[0]));
            row.push_back(print(-pair[1]));
            row.push_back(print(pair[0] - pair[1]));
        }
        rows.push_back(std::move(row));
        writePlot();
    };
    if (coldStage) {
        result = solveExtendedFixedBiasNewton(state, corrector, assemble);
        if (result.converged) {
            state = std::move(result.state);
            recordAccepted();
        }
    } else if (loadOnly) {
        const auto value = withUnitPopulations(state, [&](const auto& physical) {
            return continuity(physical).current_A_per_um;
        });
        if (value[0] - value[1] != state.drainCurrent_A_per_um)
            throw std::runtime_error("Extended SLOT Save/Load did not reproduce the saved terminal current");
        result.converged = true;
        recordAccepted();
    } else {
        const auto& circuit = sweep.at("external_circuit");
        if (circuit.value("mode", "") != "series_resistor" || circuit.value("solver", "") != "coupled_newton" ||
            circuit.value("coupled_linear_solver", "") != "direct_bordered" ||
            circuit.value("coupled_line_search_mode", "") != "full_step" ||
            circuit.value("coupled_step_control", "") != "quasistationary" ||
            circuit.value("coupled_state_predictor", "") != "constant" ||
            circuit.value("physical_acceptance_mode", "") != "report" ||
            circuit.value("coupled_diagnostic_min_updates", 0) != 0)
            throw std::invalid_argument("Unsupported extended SLOT circuit controls");
        corrector.lineSearch = false;
        corrector.maxUpdates = circuit.at("max_iterations").get<int>();
        // Correct against the original system using the same factors. The
        // unchanged linear gate ends refinement early once it is satisfied.
        corrector.refinementCorrections = std::max(corrector.refinementCorrections, 3);
        const auto& rawCircuit = wideDeck.at("sweep").at("external_circuit");
        if (const auto growth = rawCircuit.find("coupled_rhs_growth"); growth != rawCircuit.end()) {
            if (!growth->is_object() || growth->size() != 2 ||
                !growth->contains("mode") || !growth->contains("factor"))
                throw std::invalid_argument("coupled_rhs_growth must contain mode and factor");
            const auto growthMode = growth->at("mode").get<std::string>();
            const auto& factorValue = growth->at("factor");
            if ((growthMode != "report" && growthMode != "enforce") || !factorValue.is_number())
                throw std::invalid_argument("coupled_rhs_growth requires mode report/enforce and a numeric factor");
            const auto factor = factorValue.get<long double>();
            if (!std::isfinite(factor) || factor < 1.L)
                throw std::invalid_argument("coupled_rhs_growth.factor must be finite and at least 1");
            corrector.rhsGrowth = ExtendedRhsGrowthControl{factor, growthMode == "enforce"};
        }
        const auto& update = rawCircuit.at("coupled_update_convergence");
        corrector.updateError = {update.at("digits").get<int>(), update.at("potential_reference_V").get<long double>(),
            update.at("electron_reference_cm3").get<long double>(), update.at("hole_reference_cm3").get<long double>()};
        const auto& loadMode = circuit.at("load_line_convergence_mode");
        if (loadMode != "report" && loadMode != "enforce")
            throw std::invalid_argument("Unsupported extended SLOT load-line convergence mode");
        const auto rawDirection = circuit.value("current_direction", 1.);
        if (rawDirection != 1. && rawDirection != -1.)
            throw std::invalid_argument("Invalid extended SLOT current direction");
        const ExtendedSeriesResistorNewtonConfig series{
            stageControl(deck, "sweep", "resistance_ohm_um", circuit.at("resistance_ohm_um").get<double>()),
            rawCircuit.at("residual_tolerance_V").get<long double>(), update.at("contact_reference_V").get<long double>(),
            port->node_ids.size(), static_cast<int>(rawDirection), loadMode == "enforce"};
        const auto assembleCircuit = [&](const ExtendedPopulationState& point) {
            return withUnitPopulations(point, [&](const auto& physical) {
                const auto portState = continuity(physical);
                ExtendedVector portDerivative = portState.currentDerivative[0] - portState.currentDerivative[1];
                if (normalization)
                    portDerivative.tail(2 * problem.mesh.numNodes()) *= concentration_cm3;
                return ExtendedSeriesResistorFrame{
                    assembleExtendedPopulationDevice(poisson, physical, material, targets,
                        port->node_ids, portState, &geometry, equationScales),
                    std::move(portDerivative), portState.current_A_per_um[0] - portState.current_A_per_um[1]};
            });
        };
        const auto& breakCriteria = sweep.value("breakdown", nlohmann::json::object());
        const auto breakContact = breakCriteria.value("current_contact", currentContact);
        if (breakContact != currentContact)
            throw std::invalid_argument("Unsupported extended SLOT current-break contact");
        const auto currentLimit = wideDeck.at("sweep").value("breakdown", WideJson::object()).value(
            "current_abs_A_per_um", std::numeric_limits<long double>::infinity());
        bool reachedCurrent = false;
        detail::BasicDCSweepStepControlConfig<long double> control;
        control.start = start;
        control.stop = stop;
        const auto step = [&](const char* key) {
            return stageControl(deck, "sweep", key, sweep.at(key).get<double>());
        };
        control.step = step("step");
        control.initialStep = step("initial_step");
        control.minStep = step("min_step");
        control.maxStep = step("max_step");
        control.growthFactor = step("growth_factor");
        control.shrinkFactor = step("shrink_factor");
        control.growthMode = "newton_iterations";
        control.maxRetries = circuit.value("coupled_max_step_retries", std::numeric_limits<int>::max());
        control.newtonIterationLimit = corrector.maxUpdates;
        control.stopBelowMinimum = true;
        control.preserveStepAfterOutputClip = true;
        control.stopRequested = [&] { return reachedCurrent; };
        control.mandatoryPoints = checkpoints;
        detail::runTypedDCSweepStepControl<long double>(control,
            [&](long double supply, long double stepSize, int retry) {
                auto seed = state;
                seed.supplyBias_V = supply;
                try {
                    result = solveExtendedSeriesResistorNewton(seed, corrector, series, assembleCircuit);
                } catch (const std::runtime_error& error) {
                    result = {};
                    result.failure = error.what();
                }
                if (corrector.rhsGrowth) {
                    nlohmann::json report{{"supply_V", print(supply)}, {"updates", result.updates},
                        {"parent_supply_V", print(state.supplyBias_V)}, {"step_V", print(stepSize)},
                        {"retry", retry}, {"converged", result.converged}, {"failure_reason", result.failure},
                        {"mode", corrector.rhsGrowth->enforce ? "enforce" : "report"},
                        {"factor", print(corrector.rhsGrowth->factor)},
                        {"initial_norm", nullptr}, {"peak_norm", nullptr}, {"peak_ratio", nullptr},
                        {"first_exceeded_update", nullptr}};
                    if (const auto& growth = result.rhsGrowth) {
                        report["initial_norm"] = print(growth->initialNorm);
                        report["peak_norm"] = print(growth->peakNorm);
                        if (growth->initialNorm > 0)
                            report["peak_ratio"] = print(growth->peakNorm / growth->initialNorm);
                        if (growth->firstExceededUpdate)
                            report["first_exceeded_update"] = *growth->firstExceededUpdate;
                    }
                    std::cout << "extended80 rhs_growth " << report.dump() << '\n';
                }
                return detail::DCSweepStepAttemptResult{result.converged, result.updates};
            }, [&](const auto& event) {
                if (!event.converged) {
                    std::cout << "extended80 rejected " << print(event.voltage) << " V: " << result.failure
                              << " backward=" << print(result.linearBackwardError)
                              << " scaled_defect=" << print(result.scaledLinearDefect) << '\n';
                    return;
                }
                state = std::move(result.state);
                loadResidual = ArithmeticLoadLineResidual::evaluate(1.L, state.deviceBias_V,
                    series.direction * series.resistance_ohm_um, state.drainCurrent_A_per_um, state.supplyBias_V).voltage_V;
                recordAccepted();
                reachedCurrent = std::abs(state.drainCurrent_A_per_um) >= currentLimit;
                std::cout << "extended80 accepted " << print(state.supplyBias_V)
                          << " V updates=" << result.updates << " load_V=" << print(loadResidual) << '\n';
            });
        result.converged = !rows.empty() && (state.supplyBias_V == stop || reachedCurrent);
    }

    out.update({{"converged", result.converged}, {"coupled_converged", result.converged},
        {"coupled_updates", result.updates}, {"scaled_device_residual", result.scaledResidual},
        {"linear_backward_error", result.linearBackwardError}, {"scaled_linear_defect", result.scaledLinearDefect},
        {"accepted_alpha", result.acceptedAlpha}, {"factored_dimension", result.factoredDimension},
        {"update_error", result.updateError.valid ? nlohmann::json(result.updateError.rms) : nlohmann::json()},
        {"potential_scale_exact_V", print(vt)},
        {"update_references", {{"digits", corrector.updateError.digits}, {"potential_V", print(corrector.updateError.potential_V)},
            {"electron_cm3", print(corrector.updateError.electron_cm3)}, {"hole_cm3", print(corrector.updateError.hole_cm3)}}},
        {"failure_reason", result.failure}});
    if (!rows.empty()) {
        const auto& drained = currents.at(currentContact);
        out["electron_current_A_per_um"] = drained[0];
        out["hole_current_A_per_um"] = -drained[1];
        out["total_current_A_per_um"] = state.drainCurrent_A_per_um;
    }
    if (!publish || rows.empty())
        return out;
    if (out.contains("vtk_fields_skipped") && !out["vtk_fields_skipped"].empty())
        std::cout << "vtk_fields_skipped " << out["vtk_fields_skipped"].dump() << '\n';
    if (inputHashes != stateInputProvenance(deck, directory) || deckHash != extendedStageFileHash(configFile))
        throw std::runtime_error("Extended SLOT inputs changed during the stage");
    const auto archivePath = resolvePath(directory, outputArchive);
    const auto csvPath = resolvePath(directory, outputCsv);
    const auto pendingCsv = csvPath + ".extended80.pending";
    {
        std::ofstream csv(pendingCsv, std::ios::binary | std::ios::trunc);
        csv.exceptions(std::ios::badbit | std::ios::failbit);
        csv << "bias_V,outer_voltage_V,inner_voltage_V,current_contact,current_electron_A_per_um,"
               "current_hole_A_per_um,current_total_A_per_um,converged,iterations,load_line_residual_V";
        for (const auto& [name, pair] : currents)
            csv << ",contact_" << name << "_electron_A_per_um,contact_" << name
                << "_hole_A_per_um,contact_" << name << "_total_A_per_um";
        csv << '\n';
        for (const auto& row : rows) {
            for (std::size_t i = 0; i < row.size(); ++i)
                csv << (i ? "," : "") << row[i];
            csv << '\n';
        }
        csv.close();
    }
    try {
        writeExtendedPopulationArchive(archivePath, state, identity);
        std::filesystem::rename(pendingCsv, csvPath);
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(pendingCsv, ignored);
        throw;
    }
    nlohmann::json voltages = nlohmann::json::object(), totals = nlohmann::json::object();
    for (const auto& contact : contacts) {
        const auto& pair = currents.at(contact.name);
        voltages[contact.name] = print(contact.name == currentContact ? state.deviceBias_V : 0.L);
        totals[contact.name] = print(pair[0] - pair[1]);
    }
    out["accepted_points"] = rows.size();
    out["points"] = rows.size();
    out["accepted_extended80_state"] = {
        {"schema", "vela.extended80-stage-state.v1"}, {"role", "accepted"}, {"archive", outputArchive},
        {"archive_sha256", extendedStageFileHash(archivePath)}, {"csv_sha256", extendedStageFileHash(csvPath)},
        {"deck_sha256", deckHash}, {"accepted_points", rows.size()},
        {"node_count", problem.mesh.numNodes()}, {"mesh_sha256", identity.meshSha256},
        {"input_sha256", identity.inputSha256}, {"voltages_V", voltages}, {"currents_A_per_um", totals},
        {"outer_voltage_V", print(state.supplyBias_V)}};
    return out;
}

nlohmann::json checkExtendedInitialStage(const std::string& configFile)
{
    std::ifstream in(configFile);
    if (!in)
        throw std::invalid_argument("Cannot open the extended80 stage configuration");
    return runExtendedSlotStage(configFile, vela::canonicalizeDeck(nlohmann::json::parse(in)), false);
}
