#include "vela/solver/GummelSolver.h"
#include "vela/boundary/BoundaryCondition.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/DDAssembler.h"
#include "vela/solver/LinearSolver.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/io/VTKWriter.h"
#include "vela/post/ElectricFieldDiagnostics.h"
#include "vela/post/PathIonizationIntegral.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <optional>

namespace vela {



// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

namespace {

void parseImpactIonizationDrivingForceInterpolation(
    const nlohmann::json& value,
    const UnitScalingConfig& scaling,
    ImpactIonizationModelConfig& config,
    const char* context)
{
    if (!value.contains("driving_force_interpolation"))
        return;

    const auto& interpolation = value.at("driving_force_interpolation");
    if (interpolation.is_string()) {
        config.drivingForceInterpolation = interpolation.get<std::string>();
        return;
    }
    if (!interpolation.is_object()) {
        throw std::invalid_argument(
            std::string(context) +
            ": impact_ionization.driving_force_interpolation must be a string or object.");
    }

    config.drivingForceInterpolation = interpolation.value(
        "mode", config.drivingForceInterpolation);
    if (interpolation.contains("electron_ref_density_m3")) {
        config.electronDrivingForceRefDensity = scaling.concentrationToInternal(
            interpolation.at("electron_ref_density_m3").get<Real>());
    }
    if (interpolation.contains("hole_ref_density_m3")) {
        config.holeDrivingForceRefDensity = scaling.concentrationToInternal(
            interpolation.at("hole_ref_density_m3").get<Real>());
    }
}

/// Thermal voltage at temperature T [K].
inline double thermalVoltage(double T)
{
    if (T <= 0.0)
        throw std::invalid_argument("thermalVoltage: temperature_K must be positive.");
    return constants::kb * T / constants::q;
}

Real ohmicContactNetDoping(const DopingModel& doping,
                           const Contact& contact,
                           Index nodeId)
{
    const Real local = doping.netDoping(nodeId);
    if (contact.node_ids.empty())
        return local;

    Real sum = 0.0;
    for (Index nid : contact.node_ids)
        sum += doping.netDoping(nid);
    const Real mean = sum / static_cast<Real>(contact.node_ids.size());

    // Imported compensated/tie nodes can have a flipped node-owned sign at one
    // contact endpoint. Preserve local values when polarity matches, but align
    // opposite-sign outliers with the contact-local mean for Ohmic BC setup.
    if (mean == 0.0 || local == 0.0)
        return mean != 0.0 ? mean : local;
    if ((local > 0.0 && mean > 0.0) || (local < 0.0 && mean < 0.0))
        return local;
    return mean;
}

Real maxRelativePermittivityAcrossRegions(const DeviceMesh& mesh,
                                          const MaterialDatabase& matdb,
                                          Real temperature_K)
{
    Real maxEpsr = 0.0;
    for (const Region& region : mesh.regions()) {
        const Material& material = matdb.getMaterial(region.material, temperature_K);
        maxEpsr = std::max(maxEpsr, material.eps_r);
    }
    return std::max(maxEpsr, 1.0);
}

} // anonymous namespace


GummelConfig gummelConfigFromJson(const nlohmann::json& json, UnitScalingConfig scaling)
{
    GummelConfig cfg;
    cfg.inputScaling = scaling;
    cfg.maxIter = json.value("max_iter", cfg.maxIter);
    cfg.reltol = json.value("reltol", cfg.reltol);
    cfg.abstol = json.value("abstol", cfg.abstol);
    cfg.temperature_K = json.value("temperature_K", cfg.temperature_K);
    if (json.contains("carrier_statistics")) {
        const auto& statistics = json.at("carrier_statistics");
        if (statistics.is_string())
            cfg.carrierStatistics.model = statistics.get<std::string>();
        else if (statistics.is_object())
            cfg.carrierStatistics.model = statistics.value("model", cfg.carrierStatistics.model);
        else
            throw std::invalid_argument(
                "gummelConfigFromJson: carrier_statistics must be a string or object.");
    }
    if (json.contains("statistics"))
        cfg.carrierStatistics.model = json.at("statistics").get<std::string>();
    cfg.dampingPsi = json.value("damping_psi", cfg.dampingPsi);
    cfg.taun = json.value("taun", cfg.taun);
    cfg.taup = json.value("taup", cfg.taup);
    if (json.contains("srh_doping_dependence")) {
        cfg.srhDopingDependence = srhDopingDependenceConfigFromJson(
            json.at("srh_doping_dependence"), scaling);
    }
    cfg.srhDopingDependence.temperature_K = cfg.temperature_K;
    // The compiled Auger defaults are SI literals; express them in the active
    // internal unit before applying any deck override, which is already
    // written in internal units.
    cfg.augerCn = scaling.unitSystem().m6PerSToInternalAugerCoefficient(cfg.augerCn);
    cfg.augerCp = scaling.unitSystem().m6PerSToInternalAugerCoefficient(cfg.augerCp);
    cfg.augerCn = scaling.augerCoefficientToInternal(
        json.value("auger_cn_m6_per_s", cfg.augerCn));
    cfg.augerCp = scaling.augerCoefficientToInternal(
        json.value("auger_cp_m6_per_s", cfg.augerCp));
    cfg.carrierFloor = json.value("carrier_floor_m3", cfg.carrierFloor);
    if (json.contains("mobility"))
        cfg.mobility = mobilityModelConfigFromJson(json.at("mobility"), scaling);
    if (json.contains("bandgap_narrowing")) {
        const auto& value = json.at("bandgap_narrowing");
        if (value.is_string()) {
            cfg.bandgapNarrowing =
                bandgapNarrowingConfig(value.get<std::string>(), scaling);
        } else if (value.is_object()) {
            cfg.bandgapNarrowing = bandgapNarrowingConfig(
                value.value("model", cfg.bandgapNarrowing.model), scaling);
            if (value.contains("reference_doping_m3")) {
                cfg.bandgapNarrowing.referenceDoping = scaling.concentrationToInternal(
                    value.at("reference_doping_m3").get<Real>());
            }
            cfg.bandgapNarrowing.coefficient = value.value(
                "coefficient_eV", cfg.bandgapNarrowing.coefficient);
            cfg.bandgapNarrowing.smoothing = value.value(
                "smoothing", cfg.bandgapNarrowing.smoothing);
            cfg.bandgapNarrowing.offset = value.value(
                "offset_eV", cfg.bandgapNarrowing.offset);
            cfg.bandgapNarrowing.fermiStatisticsCorrection = value.value(
                "fermi_statistics_correction",
                cfg.bandgapNarrowing.fermiStatisticsCorrection);
        } else {
            throw std::invalid_argument(
                "gummelConfigFromJson: bandgap_narrowing must be a string or object.");
        }
    }

    if (json.contains("recombination")) {
        const auto& value = json.at("recombination");
        if (value.is_array())
            cfg.recombination = value.get<std::vector<std::string>>();
        else if (value.is_string())
            cfg.recombination = {value.get<std::string>()};
        else
            throw std::invalid_argument(
                "gummelConfigFromJson: recombination must be a string or string array.");
    }
    if (json.contains("band_to_band")) {
        cfg.bandToBand = bandToBandTunnelingConfigFromJson(
            json.at("band_to_band"), "gummelConfigFromJson");
    }


    if (json.contains("impact_ionization")) {
        const auto& value = json.at("impact_ionization");
        if (value.is_string()) {
            cfg.impactIonization = impactIonizationModelConfig(
                value.get<std::string>(), scaling);
        } else if (value.is_object()) {
            cfg.impactIonization = impactIonizationModelConfig(
                value.value("model", cfg.impactIonization.model), scaling);
            cfg.impactIonization.couplingMode = value.value(
                "coupling_mode", cfg.impactIonization.couplingMode);
            cfg.impactIonization.parameterSet = value.value(
                "parameter_set", cfg.impactIonization.parameterSet);
            cfg.impactIonization.drivingForce = value.value(
                "driving_force", cfg.impactIonization.drivingForce);
            cfg.impactIonization.generation = value.value(
                "generation", cfg.impactIonization.generation);
            cfg.impactIonization.currentApproximation = value.value(
                "current_approximation", cfg.impactIonization.currentApproximation);
            cfg.impactIonization.currentMagnitudeMode = value.value(
                "current_magnitude_mode", cfg.impactIonization.currentMagnitudeMode);
            cfg.impactIonization.eparallelFieldRecovery = value.value(
                "eparallel_field_recovery", cfg.impactIonization.eparallelFieldRecovery);
            cfg.impactIonization.cellReconstructedMidpointDensity = value.value(
                "cell_reconstructed_midpoint_density",
                cfg.impactIonization.cellReconstructedMidpointDensity);
            cfg.impactIonization.quasiFermiGradientDiscretization = value.value(
                "quasi_fermi_gradient_discretization",
                cfg.impactIonization.quasiFermiGradientDiscretization);
            cfg.impactIonization.contactElectricFieldFallback = value.value(
                "contact_electric_field_fallback",
                cfg.impactIonization.contactElectricFieldFallback);
            cfg.impactIonization.contactElectricFieldFallbackScope = value.value(
                "contact_electric_field_fallback_scope",
                cfg.impactIonization.contactElectricFieldFallbackScope);
            cfg.impactIonization.contactElectricFieldFallbackMode = value.value(
                "contact_electric_field_fallback_mode",
                cfg.impactIonization.contactElectricFieldFallbackMode);
            parseImpactIonizationDrivingForceInterpolation(
                value, scaling, cfg.impactIonization, "gummelConfigFromJson");
            cfg.impactIonization.sourceGeometryScale = value.value(
                "source_geometry_scale", cfg.impactIonization.sourceGeometryScale);
            cfg.impactIonization.sourceVolumePolicy = value.value(
                "source_volume_policy", cfg.impactIonization.sourceVolumePolicy);
            cfg.impactIonization.sourceVolumeFactor = value.value(
                "source_volume_factor", cfg.impactIonization.sourceVolumeFactor);
            cfg.impactIonization.sourceMappingMode = value.value(
                "source_mapping_mode", cfg.impactIonization.sourceMappingMode);
            cfg.impactIonization.edgeSourcePartition = value.value(
                "edge_source_partition", cfg.impactIonization.edgeSourcePartition);
            cfg.impactIonization.quasiFermiCarrierTruncation = value.value(
                "quasi_fermi_carrier_truncation",
                cfg.impactIonization.quasiFermiCarrierTruncation);
            cfg.impactIonization.quasiFermiCarrierTruncation = value.value(
                "quasi_fermi_carrier_trucation",
                cfg.impactIonization.quasiFermiCarrierTruncation);
            cfg.impactIonization.minimumField = scaling.electricFieldToInternal(value.value(
                "minimum_field_V_m", cfg.impactIonization.minimumField));
            cfg.impactIonization.debugRawVanOverstraeten = value.value(
                "debug_raw_vanoverstraeten",
                cfg.impactIonization.debugRawVanOverstraeten);
            cfg.impactIonization.aScale = value.value(
                "A_scale", cfg.impactIonization.aScale);
            cfg.impactIonization.bScale = value.value(
                "B_scale", cfg.impactIonization.bScale);
            if (value.contains("electron_A_m_inv")) {
                cfg.impactIonization.electronA = scaling.inverseLengthToInternal(
                    value.at("electron_A_m_inv").get<Real>());
            }
            if (value.contains("electron_B_V_m")) {
                cfg.impactIonization.electronB = scaling.electricFieldToInternal(
                    value.at("electron_B_V_m").get<Real>());
            }
            if (value.contains("hole_A_m_inv")) {
                cfg.impactIonization.holeA = scaling.inverseLengthToInternal(
                    value.at("hole_A_m_inv").get<Real>());
            }
            if (value.contains("hole_B_V_m")) {
                cfg.impactIonization.holeB = scaling.electricFieldToInternal(
                    value.at("hole_B_V_m").get<Real>());
            }
            if (value.contains("electron_a_low_m_inv")) {
                cfg.impactIonization.electronALow = scaling.inverseLengthToInternal(
                    value.at("electron_a_low_m_inv").get<Real>());
            }
            if (value.contains("electron_a_high_m_inv")) {
                cfg.impactIonization.electronAHigh = scaling.inverseLengthToInternal(
                    value.at("electron_a_high_m_inv").get<Real>());
            }
            if (value.contains("electron_b_low_V_m")) {
                cfg.impactIonization.electronBLow = scaling.electricFieldToInternal(
                    value.at("electron_b_low_V_m").get<Real>());
            }
            if (value.contains("electron_b_high_V_m")) {
                cfg.impactIonization.electronBHigh = scaling.electricFieldToInternal(
                    value.at("electron_b_high_V_m").get<Real>());
            }
            if (value.contains("hole_a_low_m_inv")) {
                cfg.impactIonization.holeALow = scaling.inverseLengthToInternal(
                    value.at("hole_a_low_m_inv").get<Real>());
            }
            if (value.contains("hole_a_high_m_inv")) {
                cfg.impactIonization.holeAHigh = scaling.inverseLengthToInternal(
                    value.at("hole_a_high_m_inv").get<Real>());
            }
            if (value.contains("hole_b_low_V_m")) {
                cfg.impactIonization.holeBLow = scaling.electricFieldToInternal(
                    value.at("hole_b_low_V_m").get<Real>());
            }
            if (value.contains("hole_b_high_V_m")) {
                cfg.impactIonization.holeBHigh = scaling.electricFieldToInternal(
                    value.at("hole_b_high_V_m").get<Real>());
            }
            if (value.contains("switch_field_V_m")) {
                cfg.impactIonization.switchField = scaling.electricFieldToInternal(
                    value.at("switch_field_V_m").get<Real>());
            }
            cfg.impactIonization.phononEnergy = value.value(
                "phonon_energy_eV", cfg.impactIonization.phononEnergy);
            cfg.impactIonization.referenceTemperature_K = value.value(
                "reference_temperature_K", cfg.impactIonization.referenceTemperature_K);
            cfg.impactIonization.temperature_K = value.value(
                "temperature_K", cfg.impactIonization.temperature_K);
            cfg.impactIonization.carrierVelocity = value.value(
                "carrier_velocity_m_s", cfg.impactIonization.carrierVelocity);
        } else {
            throw std::invalid_argument(
                "gummelConfigFromJson: impact_ionization must be a string or object.");
        }
    }
    detail::validateImpactIonizationDrivingForce(cfg.impactIonization, "gummelConfigFromJson");
    (void)BandToBandTunnelingModel(cfg.bandToBand);

    if (cfg.temperature_K <= 0.0)
        throw std::invalid_argument("gummelConfigFromJson: temperature_K must be positive.");
    (void)usesFermiDirac(cfg.carrierStatistics);
    if (cfg.carrierFloor < 0.0 || !std::isfinite(cfg.carrierFloor))
        throw std::invalid_argument(
            "gummelConfigFromJson: carrier_floor_m3 must be non-negative and finite.");

    return cfg;
}

// ---------------------------------------------------------------------------
// runGummel
// ---------------------------------------------------------------------------

namespace {

DDSolution runGummelImpl(const DeviceMesh&                          mesh,
                         const MaterialDatabase&                     matdb,
                         const DopingModel&                          doping,
                         const std::unordered_map<std::string, Real>& contactBiases,
                         const ContactSpecsMap&                       contactSpecs,
                         const GummelConfig&                          cfg,
                         const DDSolution*                           initialGuess,
                         std::vector<RegionFixedChargeSpec>           fixedCharges,
                         std::vector<InterfaceSheetChargeSpec>        sheetCharges)
{
    const Index  N   = mesh.numNodes();
    const double Vt  = thermalVoltage(cfg.temperature_K);

    // Per-node effective ni, including optional bandgap narrowing.
    std::vector<double> ni_v = detail::buildValidatedEffectiveNodeNi(
        "runGummel", mesh, matdb, doping, cfg.bandgapNarrowing, Vt);
    const std::vector<Real> Nc = detail::buildNodeDensityOfStates(
        mesh, matdb, cfg.temperature_K, true);
    const std::vector<Real> Nv = detail::buildNodeDensityOfStates(
        mesh, matdb, cfg.temperature_K, false);
    if (usesFermiDirac(cfg.carrierStatistics)) {
        for (Index node = 0; node < N; ++node) {
            if (ni_v[node] > 0.0 && (!(Nc[node] > 0.0) || !(Nv[node] > 0.0))) {
                throw std::invalid_argument(
                    "runGummel: Fermi-Dirac statistics require Nc_m3 and Nv_m3 "
                    "for every transport node.");
            }
        }
    }

    const bool useScaledUnknowns = cfg.inputScaling.isUnitScaling();
    DDScalingSpec ddScaling;
    if (useScaledUnknowns) {
        const Real epsRef = constants::eps0 *
            maxRelativePermittivityAcrossRegions(mesh, matdb, cfg.temperature_K);
        const Real niFloor = std::max(*std::max_element(ni_v.begin(), ni_v.end()), 1.0);
        const UnitScalingSystem::AutoInputs autoInputs = UnitScalingSystem::autoInputsFrom(
            mesh, doping, matdb, niFloor);
        const UnitScalingSystem scalingSystem = UnitScalingSystem::fromInputs(
            cfg.temperature_K, epsRef, autoInputs, cfg.unitScalingRefs, cfg.inputScaling.unitSystem());

        ddScaling.enabled = true;
        ddScaling.V0 = scalingSystem.V0();
        ddScaling.C0 = scalingSystem.C0();
        ddScaling.mu0 = scalingSystem.mu0();
        ddScaling.D0 = scalingSystem.D0();
        ddScaling.L0 = scalingSystem.L0();
        ddScaling.permittivityReference_F_per_m = epsRef;
        ddScaling.unitSystem = cfg.inputScaling.unitSystem();
        ddScaling.chargeAreaFactor = cfg.inputScaling.unitSystem().chargeAreaFactor();
        ddScaling.chargeLineFactor = cfg.inputScaling.unitSystem().chargeLineFactor();
        ddScaling.fieldFromCoordinateDeltaFactor = cfg.inputScaling.unitSystem().fieldFromCoordinateDeltaFactor();
        ddScaling.currentDensityLineIntegralFactor =
            cfg.inputScaling.unitSystem().currentDensityAM2PerInternal() *
            cfg.inputScaling.unitSystem().lengthMPerInternal();
    }

    // Look up a contact-region material for Schottky barrier helpers.
    auto contactMaterial = [&](const Contact& contact) -> std::optional<Material> {
        if (contact.region_id < mesh.numRegions()) {
            const Region& region = mesh.getRegion(contact.region_id);
            if (matdb.hasMaterial(region.material))
                return matdb.getMaterial(region.material, cfg.temperature_K);
        }
        // Fallback: pick any cell that contains a contact node.
        for (Index c = 0; c < mesh.numCells(); ++c) {
            const Cell& cell = mesh.getCell(c);
            for (Index nid : cell.node_ids) {
                for (Index cnid : contact.node_ids) {
                    if (nid == cnid) {
                        const Region& region = mesh.getRegion(cell.region_id);
                        if (matdb.hasMaterial(region.material))
                            return matdb.getMaterial(region.material, cfg.temperature_K);
                    }
                }
            }
        }
        return std::nullopt;
    };

    // ------------------------------------------------------------------
    // Build contact Dirichlet BCs
    // ------------------------------------------------------------------
    //   Ohmic (default):
    //     psi_contact  = V_bias + Vt * ln(n_eq / ni)
    //     n_contact    = n_eq
    //     p_contact    = ni^2 / n_eq
    //     phin = phip  = V_bias
    //   Metal gate:
    //     psi_contact = effective electrostatic gate potential
    //     no carrier-density or quasi-Fermi Dirichlet rows
    //   Schottky (prototype Dirichlet barrier):
    //     ContactState from computeSchottkyContactState()
    // ------------------------------------------------------------------
    std::unordered_map<Index, Real> psiBC;
    std::unordered_map<Index, Real> nBC;
    std::unordered_map<Index, Real> pBC;
    std::unordered_map<Index, Real> phinBC;
    std::unordered_map<Index, Real> phipBC;

    for (Index c = 0; c < mesh.numContacts(); ++c) {
        const Contact& contact = mesh.getContact(c);
        auto it = contactBiases.find(contact.name);
        if (it == contactBiases.end()) continue;
        const double Vbias = it->second;

        auto specIt = contactSpecs.find(contact.name);
        const ContactBoundarySpec* spec =
            (specIt != contactSpecs.end()) ? &specIt->second : nullptr;
        const bool isSchottky =
            spec != nullptr && spec->type == ContactType::Schottky;
        const bool isMetalGate =
            spec != nullptr && spec->type == ContactType::MetalGate;

        if (isMetalGate) {
            for (Index nid : contact.node_ids)
                psiBC[nid] = Vbias;
            continue;
        }

        if (isSchottky) {
            ContactBoundarySpec effSpec = *spec;
            effSpec.bias = Vbias; // honour the swept bias point

            const auto materialOpt = contactMaterial(contact);
            const Real bandgap_eV = materialOpt && materialOpt->bandgap_eV
                ? *materialOpt->bandgap_eV
                : std::numeric_limits<Real>::quiet_NaN();
            const Real affinity_eV = materialOpt && materialOpt->electron_affinity_eV
                ? *materialOpt->electron_affinity_eV
                : std::numeric_limits<Real>::quiet_NaN();
            const Real NcVal = materialOpt && materialOpt->Nc_m3
                ? *materialOpt->Nc_m3 : 0.0;
            const Real NvVal = materialOpt && materialOpt->Nv_m3
                ? *materialOpt->Nv_m3 : 0.0;

            for (Index nid : contact.node_ids) {
                const double ni_node = ni_v[nid];
                const double Ndop = doping.netDoping(nid);
                const ContactState state = computeSchottkyContactState(
                    effSpec, ni_node, NcVal, NvVal,
                    bandgap_eV, affinity_eV, Ndop, cfg.temperature_K);
                psiBC [nid] = state.psi;
                nBC   [nid] = state.n;
                pBC   [nid] = state.p;
                phinBC[nid] = state.phin;
                phipBC[nid] = state.phip;
            }
            continue;
        }

        for (Index nid : contact.node_ids) {
            const double ni_node  = ni_v[nid];
            const double Ndop     = ohmicContactNetDoping(doping, contact, nid);
            const EquilibriumCarrierState equilibrium = equilibriumCarrierState(
                Ndop, ni_node, Nc[nid], Nv[nid], Vt, cfg.carrierStatistics);

            psiBC [nid] = Vbias + equilibrium.potential;
            nBC   [nid] = equilibrium.n;
            pBC   [nid] = equilibrium.p;
            phinBC[nid] = Vbias;
            phipBC[nid] = Vbias;
        }
    }

    const auto scalePotential = [&](Real value) {
        return useScaledUnknowns ? (value / ddScaling.V0) : value;
    };
    const auto scaleCarrier = [&](Real value) {
        return useScaledUnknowns ? (value / ddScaling.C0) : value;
    };

    std::unordered_map<Index, Real> psiBCSolve;
    std::unordered_map<Index, Real> nBCSolve;
    std::unordered_map<Index, Real> pBCSolve;
    std::unordered_map<Index, Real> phinBCSolve;
    std::unordered_map<Index, Real> phipBCSolve;
    psiBCSolve.reserve(psiBC.size());
    nBCSolve.reserve(nBC.size());
    pBCSolve.reserve(pBC.size());
    phinBCSolve.reserve(phinBC.size());
    phipBCSolve.reserve(phipBC.size());
    for (const auto& [nid, val] : psiBC) psiBCSolve[nid] = scalePotential(val);
    for (const auto& [nid, val] : nBC) nBCSolve[nid] = scaleCarrier(val);
    for (const auto& [nid, val] : pBC) pBCSolve[nid] = scaleCarrier(val);
    for (const auto& [nid, val] : phinBC) phinBCSolve[nid] = scalePotential(val);
    for (const auto& [nid, val] : phipBC) phipBCSolve[nid] = scalePotential(val);


    // ------------------------------------------------------------------
    // Initial guess: solve linear Poisson (no carriers) for psi
    // ------------------------------------------------------------------
    const MobilityModelConfig mobilityConfig = cfg.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(
            cfg.recombination, cfg.taun, cfg.taup, cfg.srhDopingDependence);
    recombinationConfig.augerCn = cfg.augerCn;
    recombinationConfig.augerCp = cfg.augerCp;
    recombinationConfig.bandToBand = cfg.bandToBand;
    DDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg.bandgapNarrowing,
        cfg.impactIonization,
        fixedCharges,
        sheetCharges,
        ddScaling,
        cfg.carrierStatistics);

    VectorXd psi_init  = VectorXd::Zero(static_cast<int>(N));
    VectorXd n_init    = VectorXd::Zero(static_cast<int>(N));
    VectorXd p_init    = VectorXd::Zero(static_cast<int>(N));
    VectorXd phin      = VectorXd::Zero(static_cast<int>(N));
    VectorXd phip      = VectorXd::Zero(static_cast<int>(N));

    const bool useInitialGuess = initialGuess != nullptr &&
        initialGuess->psi.size() == static_cast<int>(N) &&
        initialGuess->phin.size() == static_cast<int>(N) &&
        initialGuess->phip.size() == static_cast<int>(N) &&
        initialGuess->n.size() == static_cast<int>(N) &&
        initialGuess->p.size() == static_cast<int>(N);

    if (useInitialGuess) {
        if (useScaledUnknowns) {
            psi_init = initialGuess->psi / ddScaling.V0;
            phin = initialGuess->phin / ddScaling.V0;
            phip = initialGuess->phip / ddScaling.V0;
            n_init = initialGuess->n / ddScaling.C0;
            p_init = initialGuess->p / ddScaling.C0;
        } else {
            psi_init = initialGuess->psi;
            phin = initialGuess->phin;
            phip = initialGuess->phip;
            n_init = initialGuess->n;
            p_init = initialGuess->p;
        }
    } else {
        // Override contact nodes before solving an initial Poisson problem.
        for (const auto& [nid, val] : psiBCSolve)
            psi_init(static_cast<int>(nid)) = val;

        // Solve initial Poisson (no free carriers)
        {
            VectorXd n_zero = VectorXd::Zero(static_cast<int>(N));
            VectorXd p_zero = VectorXd::Zero(static_cast<int>(N));
            assembler.assemblePoissonWithCarriers(n_zero, p_zero, psi_init);
            assembler.applyDirichlet(psiBCSolve);
            LinearSolver ls;
            psi_init = ls.solve(assembler.matrix(), assembler.rhs());
        }

        for (const auto& [nid, val] : phinBCSolve) phin(static_cast<int>(nid)) = val;
        for (const auto& [nid, val] : phipBCSolve) phip(static_cast<int>(nid)) = val;

        for (Index i = 0; i < N; ++i) {
            const int    ii     = static_cast<int>(i);
            const double ni_i   = ni_v[i];
            const double psiSi = useScaledUnknowns ? psi_init(ii) * ddScaling.V0 : psi_init(ii);
            const double phinSi = useScaledUnknowns ? phin(ii) * ddScaling.V0 : phin(ii);
            const double phipSi = useScaledUnknowns ? phip(ii) * ddScaling.V0 : phip(ii);
            const double nSi = electronDensity(
                ni_i, Nc[i], psiSi, phinSi, Vt, cfg.carrierStatistics);
            const double pSi = holeDensity(
                ni_i, Nv[i], psiSi, phipSi, Vt, cfg.carrierStatistics);
            n_init(ii) = useScaledUnknowns ? (nSi / ddScaling.C0) : nSi;
            p_init(ii) = useScaledUnknowns ? (pSi / ddScaling.C0) : pSi;
        }
    }

    // Enforce the new bias point's Ohmic contact values even when the
    // previous bias point is used as the initial guess.
    for (const auto& [nid, val] : psiBCSolve)  psi_init(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : phinBCSolve) phin(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : phipBCSolve) phip(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : nBCSolve)    n_init(static_cast<int>(nid)) = val;
    for (const auto& [nid, val] : pBCSolve)    p_init(static_cast<int>(nid)) = val;

    // ------------------------------------------------------------------
    // Working copies
    // ------------------------------------------------------------------
    VectorXd psi = psi_init;
    VectorXd n   = n_init;
    VectorXd p   = p_init;

    LinearSolver ls;
    int iters = 0;
    bool converged = false;
    const double carrierFloorSolve = useScaledUnknowns
        ? cfg.carrierFloor / ddScaling.C0
        : cfg.carrierFloor;

    // ------------------------------------------------------------------
    // Gummel iteration
    // ------------------------------------------------------------------
    for (int iter = 0; iter < cfg.maxIter; ++iter) {
        ++iters;

        // ---- a. Solve Poisson with current carriers ----
        assembler.assemblePoissonWithCarriers(n, p, psi);
        assembler.applyDirichlet(psiBCSolve);
        VectorXd psi_new = ls.solve(assembler.matrix(), assembler.rhs());

        // Damped update: compute full correction then apply with damping factor
        const VectorXd dpsi         = psi_new - psi;
        const VectorXd dpsi_applied = cfg.dampingPsi * dpsi;
        psi += dpsi_applied;

        // ---- b. Solve electron continuity for n ----
        assembler.assembleElectronContinuity(psi, n, p);
        assembler.applyDirichlet(nBCSolve);
        VectorXd n_new = ls.solve(assembler.matrix(), assembler.rhs());

        // Enforce positivity and keep quasi-Fermi reconstruction well-defined.
        for (Index i = 0; i < N; ++i) {
            const int ii = static_cast<int>(i);
            if (nBCSolve.find(i) == nBCSolve.end() && n_new(ii) < carrierFloorSolve)
                n_new(ii) = carrierFloorSolve;
        }

        // ---- c. Solve hole continuity for p ----
        assembler.assembleHoleContinuity(psi, n_new, p);
        assembler.applyDirichlet(pBCSolve);
        VectorXd p_new = ls.solve(assembler.matrix(), assembler.rhs());

        for (Index i = 0; i < N; ++i) {
            const int ii = static_cast<int>(i);
            if (pBCSolve.find(i) == pBCSolve.end() && p_new(ii) < carrierFloorSolve)
                p_new(ii) = carrierFloorSolve;
        }

        // ---- d. Update quasi-Fermi potentials ----
        for (Index i = 0; i < N; ++i) {
            const int    ii   = static_cast<int>(i);
            const double ni_i = ni_v[i];
            const double psiSi = useScaledUnknowns ? psi(ii) * ddScaling.V0 : psi(ii);
            const double nSi = useScaledUnknowns ? n_new(ii) * ddScaling.C0 : n_new(ii);
            const double pSi = useScaledUnknowns ? p_new(ii) * ddScaling.C0 : p_new(ii);
            if (ni_i > 0.0 && nSi > 0.0) {
                const double phinSi = electronQuasiFermiPotential(
                    ni_i, Nc[i], psiSi, nSi, Vt, cfg.carrierStatistics);
                phin(ii) = useScaledUnknowns ? (phinSi / ddScaling.V0) : phinSi;
            }
            if (ni_i > 0.0 && pSi > 0.0) {
                const double phipSi = holeQuasiFermiPotential(
                    ni_i, Nv[i], psiSi, pSi, Vt, cfg.carrierStatistics);
                phip(ii) = useScaledUnknowns ? (phipSi / ddScaling.V0) : phipSi;
            }
        }
        // Restore Dirichlet quasi-Fermi values at contacts
        for (const auto& [nid, val] : phinBCSolve) phin(static_cast<int>(nid)) = val;
        for (const auto& [nid, val] : phipBCSolve) phip(static_cast<int>(nid)) = val;

        // ---- e. Convergence check ----
        // Base residuals on the actual applied updates so that damping < 1
        // does not prevent detection of convergence.
        const double dpsi_norm = dpsi_applied.norm();
        const double psi_norm  = psi.norm();
        const double rel_err   = (psi_norm > 1.0e-30)
                                     ? dpsi_norm / psi_norm
                                     : dpsi_norm;

        const VectorXd dn = n_new - n;
        const double dn_norm = dn.norm();
        const double n_norm  = n.norm();
        const double rel_n   = (n_norm > 1.0e-30) ? dn_norm / n_norm : dn_norm;

        const VectorXd dp = p_new - p;
        const double dp_norm = dp.norm();
        const double p_norm  = p.norm();
        const double rel_p   = (p_norm > 1.0e-30) ? dp_norm / p_norm : dp_norm;

        n = n_new;
        p = p_new;

        const bool relativeConverged =
            rel_err < cfg.reltol && rel_n < cfg.reltol && rel_p < cfg.reltol;
        const bool absoluteConverged = cfg.abstol > 0.0 &&
            dpsi_norm <= cfg.abstol && dn_norm <= cfg.abstol && dp_norm <= cfg.abstol;
        if (relativeConverged || absoluteConverged) {
            converged = true;
            break;
        }
    }

    DDSolution sol;
    if (useScaledUnknowns) {
        sol.psi = psi * ddScaling.V0;
        sol.phin = phin * ddScaling.V0;
        sol.phip = phip * ddScaling.V0;
        sol.n = n * ddScaling.C0;
        sol.p = p * ddScaling.C0;
    } else {
        sol.psi   = psi;
        sol.phin  = phin;
        sol.phip  = phip;
        sol.n     = n;
        sol.p     = p;
    }
    sol.iters = iters;
    sol.converged = converged;
    return sol;
}

} // anonymous namespace

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                          cfg)
{
    ContactSpecsMap emptySpecs;
    return runGummelImpl(mesh, matdb, doping, contactBiases, emptySpecs, cfg, nullptr, {}, {});
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess)
{
    ContactSpecsMap emptySpecs;
    return runGummelImpl(mesh, matdb, doping, contactBiases, emptySpecs, cfg, &initialGuess, {}, {});
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg)
{
    return runGummelImpl(mesh, matdb, doping, contactBiases, contactSpecs, cfg, nullptr, {}, {});
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess)
{
    return runGummelImpl(mesh, matdb, doping, contactBiases, contactSpecs, cfg, &initialGuess, {}, {});
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     std::vector<RegionFixedChargeSpec>           fixedCharges,
                     std::vector<InterfaceSheetChargeSpec>        sheetCharges)
{
    return runGummelImpl(
        mesh, matdb, doping, contactBiases, contactSpecs, cfg, nullptr,
        std::move(fixedCharges), std::move(sheetCharges));
}

DDSolution runGummel(const DeviceMesh&                          mesh,
                     const MaterialDatabase&                     matdb,
                     const DopingModel&                          doping,
                     const std::unordered_map<std::string, Real>& contactBiases,
                     const ContactSpecsMap&                       contactSpecs,
                     const GummelConfig&                          cfg,
                     const DDSolution&                           initialGuess,
                     std::vector<RegionFixedChargeSpec>           fixedCharges,
                     std::vector<InterfaceSheetChargeSpec>        sheetCharges)
{
    return runGummelImpl(
        mesh, matdb, doping, contactBiases, contactSpecs, cfg, &initialGuess,
        std::move(fixedCharges), std::move(sheetCharges));
}


namespace {

std::vector<Point3> cellFieldVectorsVcm(const std::vector<CellField2>& fields, Real scaleToVcm)
{
    std::vector<Point3> out(fields.size(), Point3::Zero());
    for (std::size_t i = 0; i < fields.size(); ++i)
        out[i] = Point3{fields[i].vector.x() * scaleToVcm, fields[i].vector.y() * scaleToVcm, 0.0};
    return out;
}

std::vector<Real> cellFieldMagnitudesVcm(const std::vector<CellField2>& fields, Real scaleToVcm)
{
    std::vector<Real> out(fields.size(), 0.0);
    for (std::size_t i = 0; i < fields.size(); ++i)
        out[i] = fields[i].magnitude * scaleToVcm;
    return out;
}

std::vector<Point3> nodeFieldVectorsVcm(const std::vector<NodeField2>& fields, Real scaleToVcm)
{
    std::vector<Point3> out(fields.size(), Point3::Zero());
    for (std::size_t i = 0; i < fields.size(); ++i)
        out[i] = Point3{fields[i].vector.x() * scaleToVcm, fields[i].vector.y() * scaleToVcm, 0.0};
    return out;
}

std::vector<Real> nodeFieldMagnitudesVcm(const std::vector<NodeField2>& fields, Real scaleToVcm)
{
    std::vector<Real> out(fields.size(), 0.0);
    for (std::size_t i = 0; i < fields.size(); ++i)
        out[i] = fields[i].magnitude * scaleToVcm;
    return out;
}

void writeRecoveredElectricFields(VTKWriter& writer,
                                  const DeviceMesh& mesh,
                                  const DDSolution& sol,
                                  const UnitScalingConfig& scaling)
{
    const Real scaleToVcm = scaling.unitSystem().fieldFromCoordinateDeltaFactor() *
        scaling.unitSystem().electricFieldVPerMPerInternal() / 100.0;
    const auto cellElectric = computeCellElectricField(mesh, sol.psi);
    const auto cellElectronQf = computeCellGradElectronQuasiFermi(mesh, sol.phin);
    const auto cellHoleQf = computeCellGradHoleQuasiFermi(mesh, sol.phip);
    writer.addCellVector("CellElectricField", cellFieldVectorsVcm(cellElectric, scaleToVcm));
    writer.addCellScalar("CellElectricFieldMagnitude", cellFieldMagnitudesVcm(cellElectric, scaleToVcm));
    writer.addCellVector("CellGradElectronQuasiFermi", cellFieldVectorsVcm(cellElectronQf, scaleToVcm));
    writer.addCellScalar("CellGradElectronQuasiFermiMagnitude", cellFieldMagnitudesVcm(cellElectronQf, scaleToVcm));
    writer.addCellVector("CellGradHoleQuasiFermi", cellFieldVectorsVcm(cellHoleQf, scaleToVcm));
    writer.addCellScalar("CellGradHoleQuasiFermiMagnitude", cellFieldMagnitudesVcm(cellHoleQf, scaleToVcm));

    const auto area = computeNodeElectricFieldAreaAverage(mesh, sol.psi);
    const auto ls1d = computeNodeElectricFieldLeastSquares(
        mesh, sol.psi, ElectricFieldLeastSquaresWeight::InverseDistance);
    const auto ls1d2 = computeNodeElectricFieldLeastSquares(
        mesh, sol.psi, ElectricFieldLeastSquaresWeight::InverseDistanceSquared);
    const auto circum1d = computeNodeElectricFieldCircumcenterRecovery(
        mesh, sol.psi, ElectricFieldCircumcenterWeight::InverseDistance);
    const auto circumArea1d = computeNodeElectricFieldCircumcenterRecovery(
        mesh, sol.psi, ElectricFieldCircumcenterWeight::AreaOverDistance);
    const auto spr = computeNodeElectricFieldSPR(mesh, sol.psi);

    writer.addNodeScalar("NodeElectricField_AreaAverage", nodeFieldMagnitudesVcm(area, scaleToVcm));
    writer.addNodeVector("NodeElectricField_AreaAverageVector", nodeFieldVectorsVcm(area, scaleToVcm));
    writer.addNodeScalar("NodeElectricField_LS_1overD", nodeFieldMagnitudesVcm(ls1d, scaleToVcm));
    writer.addNodeVector("NodeElectricField_LS_1overDVector", nodeFieldVectorsVcm(ls1d, scaleToVcm));
    writer.addNodeScalar("NodeElectricField_LS_1overD2", nodeFieldMagnitudesVcm(ls1d2, scaleToVcm));
    writer.addNodeVector("NodeElectricField_LS_1overD2Vector", nodeFieldVectorsVcm(ls1d2, scaleToVcm));
    writer.addNodeScalar("NodeElectricField_Circumcenter1overD", nodeFieldMagnitudesVcm(circum1d, scaleToVcm));
    writer.addNodeVector("NodeElectricField_Circumcenter1overDVector", nodeFieldVectorsVcm(circum1d, scaleToVcm));
    writer.addNodeScalar("NodeElectricField_CircumcenterAreaOverD", nodeFieldMagnitudesVcm(circumArea1d, scaleToVcm));
    writer.addNodeVector("NodeElectricField_CircumcenterAreaOverDVector", nodeFieldVectorsVcm(circumArea1d, scaleToVcm));
    writer.addNodeScalar("NodeElectricField_SPR", nodeFieldMagnitudesVcm(spr, scaleToVcm));
    writer.addNodeVector("NodeElectricField_SPRVector", nodeFieldVectorsVcm(spr, scaleToVcm));
}

} // namespace
// ---------------------------------------------------------------------------
// writeDDSolutionVTK
// ---------------------------------------------------------------------------

void writeDDSolutionVTK(const std::string& filename,
                        const DeviceMesh&  mesh,
                        const DopingModel& doping,
                        const DDSolution&  sol,
                        UnitScalingConfig scaling)
{
    const Index N = mesh.numNodes();

    VTKWriter writer(filename, mesh);
    writer.write();
    writeRecoveredElectricFields(writer, mesh, sol, scaling);

    auto toVec = [&](const VectorXd& v) {
        std::vector<Real> out(N);
        for (Index i = 0; i < N; ++i)
            out[i] = v(static_cast<int>(i));
        return out;
    };

    writer.addNodeScalar("Potential",            toVec(sol.psi));
    writer.addNodeScalar("ElectronQuasiFermi",   toVec(sol.phin));
    writer.addNodeScalar("HoleQuasiFermi",       toVec(sol.phip));
    writer.addNodeScalar("Electrons",            toVec(sol.n));
    writer.addNodeScalar("Holes",                toVec(sol.p));

    std::vector<Real> netDop(N);
    for (Index i = 0; i < N; ++i)
        netDop[i] = doping.netDoping(i);
    writer.addNodeScalar("NetDoping", netDop);
}

void writeDDSolutionVTK(const std::string& filename,
                        const DeviceMesh& mesh,
                        const MaterialDatabase& matdb,
                        const DopingModel& doping,
                        const DDSolution& sol,
                        const MobilityModelConfig& mobilityConfig,
                        const RecombinationModelConfig& recombinationConfig,
                        const ImpactIonizationModelConfig& impactIonizationConfig,
                        const BandgapNarrowingConfig& bandgapNarrowingConfig,
                        Real temperature_K,
                        UnitScalingConfig scaling,
                        const CarrierStatisticsConfig& carrierStatistics)
{
    const Index N = mesh.numNodes();
    VTKWriter writer(filename, mesh);
    writer.write();
    writeRecoveredElectricFields(writer, mesh, sol, scaling);

    auto toVec = [&](const VectorXd& v) {
        std::vector<Real> out(N);
        for (Index i = 0; i < N; ++i)
            out[i] = v(static_cast<int>(i));
        return out;
    };

    writer.addNodeScalar("Potential",            toVec(sol.psi));
    writer.addNodeScalar("ElectronQuasiFermi",   toVec(sol.phin));
    writer.addNodeScalar("HoleQuasiFermi",       toVec(sol.phip));
    writer.addNodeScalar("Electrons",            toVec(sol.n));
    writer.addNodeScalar("Holes",                toVec(sol.p));

    std::vector<Real> netDop(N);
    for (Index i = 0; i < N; ++i)
        netDop[i] = doping.netDoping(i);
    writer.addNodeScalar("NetDoping", netDop);

    const auto nodeCells = detail::buildNodeCellMap(mesh);
    const Real fieldFactor = scaling.unitSystem().fieldFromCoordinateDeltaFactor();
    const Real scaleToVcm = scaling.unitSystem().electricFieldVPerMPerInternal() / 100.0;
    const std::vector<Point2> electricFieldGradient_V_m =
        detail::computeNodeWeightedLeastSquaresGradients(
            mesh, nodeCells, [&](Index node) { return sol.psi(static_cast<int>(node)); });
    const bool hasReferencedElectronQf =
        sol.phinIncrement.size() == sol.phin.size();
    const bool hasReferencedHoleQf =
        sol.phipIncrement.size() == sol.phip.size();
    const VectorXd& electronQfForDifferences =
        hasReferencedElectronQf ? sol.phinIncrement : sol.phin;
    const VectorXd& holeQfForDifferences =
        hasReferencedHoleQf ? sol.phipIncrement : sol.phip;
    const std::vector<Point2> electronQfGradientVector_V_m =
        detail::computeNodeWeightedLeastSquaresGradients(
            mesh, nodeCells, [&](Index node) {
                return electronQfForDifferences(static_cast<int>(node));
            });
    const std::vector<Point2> holeQfGradientVector_V_m =
        detail::computeNodeWeightedLeastSquaresGradients(
            mesh, nodeCells, [&](Index node) {
                return holeQfForDifferences(static_cast<int>(node));
            });
    std::vector<Real> electricField_V_m(N, 0.0);
    std::vector<Real> electronQfGradient_V_m(N, 0.0);
    std::vector<Real> holeQfGradient_V_m(N, 0.0);
    for (Index i = 0; i < N; ++i) {
        electricField_V_m[i] = electricFieldGradient_V_m[i].norm() * fieldFactor;
        electronQfGradient_V_m[i] = electronQfGradientVector_V_m[i].norm() * fieldFactor;
        holeQfGradient_V_m[i] = holeQfGradientVector_V_m[i].norm() * fieldFactor;
    }
    const bool qfMobility =
        mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
    const std::vector<Real>& electronMobilityDrive_V_m = qfMobility
        ? electronQfGradient_V_m
        : electricField_V_m;
    const std::vector<Real>& holeMobilityDrive_V_m = qfMobility
        ? holeQfGradient_V_m
        : electricField_V_m;

    std::vector<Real> electricField_V_cm(N, 0.0);
    std::vector<Point3> electricFieldVector_V_cm(N, Point3::Zero());
    for (Index i = 0; i < N; ++i) {
        electricField_V_cm[i] = electricField_V_m[i] * scaleToVcm;
        electricFieldVector_V_cm[i] = Point3{
            -electricFieldGradient_V_m[i].x() * fieldFactor * scaleToVcm,
            -electricFieldGradient_V_m[i].y() * fieldFactor * scaleToVcm,
            0.0};
    }
    writer.addNodeScalar("ElectricField", electricField_V_cm);
    writer.addNodeVector("ElectricFieldVector", electricFieldVector_V_cm);

    const Real Vt = constants::kb * temperature_K / constants::q;
    const std::unique_ptr<BandgapNarrowing> bgn =
        makeBandgapNarrowingModel(bandgapNarrowingConfig);
    const RecombinationModel recombination(recombinationConfig);
    const std::unique_ptr<ImpactIonizationModel> impact =
        makeImpactIonizationModel(impactIonizationConfig);
    const std::unique_ptr<MobilityModel> mobility = makeMobilityModel(mobilityConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh, matdb, temperature_K);
    const std::vector<Real> effectiveNi = detail::buildValidatedEffectiveNodeNi(
        "writeDDSolutionVTK",
        mesh,
        matdb,
        doping,
        bandgapNarrowingConfig,
        Vt);
    const std::vector<Real> densityOfStatesNc =
        detail::buildNodeDensityOfStates(mesh, matdb, temperature_K, true);
    const std::vector<Real> densityOfStatesNv =
        detail::buildNodeDensityOfStates(mesh, matdb, temperature_K, false);
    writer.addNodeScalar("EffectiveIntrinsicDensity", effectiveNi);
    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(impactIonizationConfig);
    const std::vector<Real> electronImpactQfGradient_V_m = qfImpact
        ? detail::computeElectronAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig, mesh, nodeCells, sol.psi, sol.phin, sol.n,
            effectiveNi, Vt)
        : std::vector<Real>{};
    const std::vector<Real> holeImpactQfGradient_V_m = qfImpact
        ? detail::computeHoleAvalancheNodeQuasiFermiDrivingFields(
            impactIonizationConfig, mesh, nodeCells, sol.psi, sol.phip, sol.p,
            effectiveNi, Vt)
        : std::vector<Real>{};
    const std::vector<Real>& electronDrivingField_V_m = qfImpact
        ? electronImpactQfGradient_V_m
        : electricField_V_m;
    const std::vector<Real>& holeDrivingField_V_m = qfImpact
        ? holeImpactQfGradient_V_m
        : electricField_V_m;

    std::vector<Material> nodeMaterials(N);
    std::vector<bool> seen(N, false);
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        const Region& region = mesh.getRegion(cell.region_id);
        if (!matdb.hasMaterial(region.material))
            continue;
        const Material material = matdb.getMaterial(region.material, temperature_K);
        for (Index nodeId : cell.node_ids) {
            if (!seen[nodeId]) {
                nodeMaterials[nodeId] = material;
                seen[nodeId] = true;
            }
        }
    }

    std::vector<Real> srh(N, 0.0);
    std::vector<Real> bandToBandGeneration(N, 0.0);
    std::vector<Real> avalanche(N, 0.0);
    std::vector<Real> electronMobility(N, 0.0);
    std::vector<Real> holeMobility(N, 0.0);
    std::vector<Real> electronLowFieldMobility(N, 0.0);
    std::vector<Real> holeLowFieldMobility(N, 0.0);
    std::vector<Real> electronHighFieldDrive_V_cm(N, 0.0);
    std::vector<Real> holeHighFieldDrive_V_cm(N, 0.0);
    std::vector<Real> electronMobilityLimiter(N, 0.0);
    std::vector<Real> holeMobilityLimiter(N, 0.0);
    std::vector<Point3> electronDriftCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Point3> electronDiffusionCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Point3> electronCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Point3> holeDriftCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Point3> holeDiffusionCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Point3> holeCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Point3> totalCurrentDensity_A_cm2(N, Point3::Zero());
    std::vector<Real> electronVelocity(N, 0.0);
    std::vector<Real> holeVelocity(N, 0.0);
    std::vector<Real> electronAlphaAvalanche(N, 0.0);
    std::vector<Real> holeAlphaAvalanche(N, 0.0);
    std::vector<Real> electronImpactIonizationDrive(N, 0.0);
    std::vector<Real> holeImpactIonizationDrive(N, 0.0);
    std::vector<Real> electronIonIntegral(N, 0.0);
    std::vector<Real> holeIonIntegral(N, 0.0);
    std::vector<Real> meanIonIntegral(N, 0.0);
    const std::vector<Real> bandToBandSourceIntegrals =
        detail::bandToBandGenerationNodeSourceIntegrals(
            recombination.bandToBand(), scaling.unitSystem(), mesh,
            cellMaterials, sol.psi, fieldFactor);
    const std::vector<Real> transportNodeAreas =
        detail::transportNodeLumpedAreas(mesh, cellMaterials);
    const bool triangleGssAvalanche =
        detail::usesTriangleGssAvalancheSource(impactIonizationConfig);
    const std::vector<detail::TriangleGssAvalancheSourceRecord>
        triangleGssAvalancheRecords = triangleGssAvalanche
        ? detail::triangleGssAvalancheSourceRecords(
            impactIonizationConfig,
            *impact,
            mobilityConfig,
            *mobility,
            edgeCells,
            mesh,
            doping,
            cellMaterials,
            sol.psi,
            sol.phin,
            sol.phip,
            sol.n,
            sol.p,
            effectiveNi,
            Vt,
            fieldFactor)
        : std::vector<detail::TriangleGssAvalancheSourceRecord>{};
    const std::vector<detail::SgEdgeCurrentAvalancheSourceRecord> sgAvalancheRecords =
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig)
            && !triangleGssAvalanche
        ? detail::sgEdgeCurrentAvalancheSourceRecords(
            impactIonizationConfig,
            *impact,
            mobilityConfig,
            *mobility,
            edgeCells,
            mesh,
            doping,
            cellMaterials,
            sol.psi,
            sol.phin,
            sol.phip,
            sol.n,
            sol.p,
            effectiveNi,
            Vt,
            fieldFactor,
            false,
            carrierStatistics,
            densityOfStatesNc,
            densityOfStatesNv)
        : std::vector<detail::SgEdgeCurrentAvalancheSourceRecord>{};
    const std::vector<Real> currentDensityAvalancheSourceIntegrals =
        detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig)
        ? detail::currentDensityAvalancheSourceIntegrals(
            impactIonizationConfig,
            *impact,
            mobilityConfig,
            *mobility,
            edgeCells,
            mesh,
            doping,
            cellMaterials,
            sol.psi,
            sol.phin,
            sol.phip,
            sol.n,
            sol.p,
            effectiveNi,
            Vt,
            fieldFactor,
            carrierStatistics,
            densityOfStatesNc,
            densityOfStatesNv)
        : std::vector<Real>{};
    for (Index i = 0; i < N; ++i) {
        const int row = static_cast<int>(i);
        const Real n = sol.n(row);
        const Real p = sol.p(row);
        const Real deltaEg = bgn->deltaEg(doping.totalImpurity(i), n, p);
        const Real ni = effectiveIntrinsicDensity(nodeMaterials[i].ni, Vt, deltaEg);
        // Match the cancellation-free source used by CoupledDDAssembler.
        // Re-forming n*p-ni^2 from rounded output densities can create a
        // completely artificial SRH rate in deep depletion.  Referenced
        // quasi-Fermi increments retain the sub-ulp carrier imbalance even
        // when the absolute contact bias is tens of volts.
        const Real electronQf = electronQfForDifferences(row)
            + (hasReferencedElectronQf ? sol.electronQfReference_V : 0.0);
        const Real holeQf = holeQfForDifferences(row)
            + (hasReferencedHoleQf ? sol.holeQfReference_V : 0.0);
        const Real excessProduct = ni * ni *
            std::expm1((holeQf - electronQf) / Vt);
        srh[i] = recombination.srhRateFromExcessProduct(
            excessProduct, n, p, ni,
            recombination.srhDopingConcentration(
                doping.donors(i), doping.acceptors(i)));
        bandToBandGeneration[i] = transportNodeAreas[i] > 0.0
            ? bandToBandSourceIntegrals[i] / transportNodeAreas[i]
            : 0.0;
        if (detail::usesEdgeCurrentAvalancheSource(impactIonizationConfig)) {
            avalanche[i] = mesh.getNode(i).volume > 0.0
                ? currentDensityAvalancheSourceIntegrals[i] / mesh.getNode(i).volume
                : 0.0;
        } else {
            avalanche[i] = detail::impactIonizationGenerationRate(
                impactIonizationConfig,
                *impact,
                mobilityConfig,
                *mobility,
                nodeCells,
                mesh,
                doping,
                cellMaterials,
                i,
                electricField_V_m[i],
                electronDrivingField_V_m[i],
                holeDrivingField_V_m[i],
                n,
                p);
        }
        const Real electronImpactField = detail::electronAvalancheDrivingField(
            impactIonizationConfig, electronDrivingField_V_m[i], electricField_V_m[i], n);
        const Real holeImpactField = detail::holeAvalancheDrivingField(
            impactIonizationConfig, holeDrivingField_V_m[i], electricField_V_m[i], p);
        electronImpactIonizationDrive[i] = std::abs(electronImpactField);
        holeImpactIonizationDrive[i] = std::abs(holeImpactField);
        electronAlphaAvalanche[i] = impact->electronCoefficient(electronImpactField);
        holeAlphaAvalanche[i] = impact->holeCoefficient(holeImpactField);
        const Real electronMobilityField = electronMobilityDrive_V_m[i];
        const Real holeMobilityField = holeMobilityDrive_V_m[i];
        electronHighFieldDrive_V_cm[i] = electronMobilityField * scaleToVcm;
        holeHighFieldDrive_V_cm[i] = holeMobilityField * scaleToVcm;
        electronLowFieldMobility[i] = detail::nodeMobility(
            nodeCells, mesh, doping, *mobility, cellMaterials, i,
            CarrierType::Electron, 0.0, &mobilityConfig);
        holeLowFieldMobility[i] = detail::nodeMobility(
            nodeCells, mesh, doping, *mobility, cellMaterials, i,
            CarrierType::Hole, 0.0, &mobilityConfig);
        electronMobility[i] = detail::nodeMobility(
            nodeCells, mesh, doping, *mobility, cellMaterials, i,
            CarrierType::Electron, electronMobilityField, &mobilityConfig);
        holeMobility[i] = detail::nodeMobility(
            nodeCells, mesh, doping, *mobility, cellMaterials, i,
            CarrierType::Hole, holeMobilityField, &mobilityConfig);
        if (electronLowFieldMobility[i] > 0.0)
            electronMobilityLimiter[i] = electronMobility[i] / electronLowFieldMobility[i];
        if (holeLowFieldMobility[i] > 0.0)
            holeMobilityLimiter[i] = holeMobility[i] / holeLowFieldMobility[i];
        const Real electronDriftScale = constants::q * electronMobility[i] * n / 1.0e4;
        const Real holeDriftScale = constants::q * holeMobility[i] * p / 1.0e4;
        electronDriftCurrentDensity_A_cm2[i] = Point3{
            electronDriftScale * electricFieldGradient_V_m[i].x(),
            electronDriftScale * electricFieldGradient_V_m[i].y(),
            0.0};
        electronCurrentDensity_A_cm2[i] = Point3{
            electronDriftScale * electronQfGradientVector_V_m[i].x(),
            electronDriftScale * electronQfGradientVector_V_m[i].y(),
            0.0};
        electronDiffusionCurrentDensity_A_cm2[i] =
            electronCurrentDensity_A_cm2[i] - electronDriftCurrentDensity_A_cm2[i];
        holeDriftCurrentDensity_A_cm2[i] = Point3{
            holeDriftScale * electricFieldGradient_V_m[i].x(),
            holeDriftScale * electricFieldGradient_V_m[i].y(),
            0.0};
        holeCurrentDensity_A_cm2[i] = Point3{
            holeDriftScale * holeQfGradientVector_V_m[i].x(),
            holeDriftScale * holeQfGradientVector_V_m[i].y(),
            0.0};
        holeDiffusionCurrentDensity_A_cm2[i] =
            holeCurrentDensity_A_cm2[i] - holeDriftCurrentDensity_A_cm2[i];
        totalCurrentDensity_A_cm2[i] = electronCurrentDensity_A_cm2[i] + holeCurrentDensity_A_cm2[i];
        electronVelocity[i] = electronMobility[i] * std::abs(electronImpactField);
        holeVelocity[i] = holeMobility[i] * std::abs(holeImpactField);
    }

    const Real alphaLengthUnitFactor =
        scaling.unitSystem().inverseLengthMInvPerInternal() *
        scaling.unitSystem().lengthMPerInternal();
    if (!triangleGssAvalancheRecords.empty()) {
        for (const auto& record : triangleGssAvalancheRecords) {
            const Real halfLength =
                0.5 * record.edgeLength * alphaLengthUnitFactor;
            electronIonIntegral[record.node0] += record.electronAlpha * halfLength;
            electronIonIntegral[record.node1] += record.electronAlpha * halfLength;
            holeIonIntegral[record.node0] += record.holeAlpha * halfLength;
            holeIonIntegral[record.node1] += record.holeAlpha * halfLength;
        }
    } else if (!sgAvalancheRecords.empty()) {
        for (const auto& record : sgAvalancheRecords) {
            const Real halfLength =
                0.5 * record.edgeLength * alphaLengthUnitFactor;
            electronIonIntegral[record.node0] += record.electronAlpha * halfLength;
            electronIonIntegral[record.node1] += record.electronAlpha * halfLength;
            holeIonIntegral[record.node0] += record.holeAlpha * halfLength;
            holeIonIntegral[record.node1] += record.holeAlpha * halfLength;
        }
    } else {
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            const Real halfLength =
                0.5 * edge.length * alphaLengthUnitFactor;
            const Real electronAlpha =
                0.5 * (electronAlphaAvalanche[edge.n0] + electronAlphaAvalanche[edge.n1]);
            const Real holeAlpha =
                0.5 * (holeAlphaAvalanche[edge.n0] + holeAlphaAvalanche[edge.n1]);
            electronIonIntegral[edge.n0] += electronAlpha * halfLength;
            electronIonIntegral[edge.n1] += electronAlpha * halfLength;
            holeIonIntegral[edge.n0] += holeAlpha * halfLength;
            holeIonIntegral[edge.n1] += holeAlpha * halfLength;
        }
    }
    for (Index i = 0; i < N; ++i) {
        meanIonIntegral[i] = 0.5 * (electronIonIntegral[i] + holeIonIntegral[i]);
    }
    std::vector<PathIonizationSegment> pathIonizationSamples;
    if (!sgAvalancheRecords.empty()) {
        pathIonizationSamples.reserve(sgAvalancheRecords.size());
        for (const auto& record : sgAvalancheRecords) {
            pathIonizationSamples.push_back({
                record.edgeId,
                record.node0,
                record.node1,
                scaling.unitSystem().internalLengthToMeters(record.edgeLength),
                scaling.unitSystem().internalInverseLengthToMInv(record.electronAlpha),
                scaling.unitSystem().internalInverseLengthToMInv(record.holeAlpha),
                scaling.unitSystem().internalElectricFieldToVPerM(
                    std::abs(record.electricField))});
        }
    } else if (!triangleGssAvalancheRecords.empty()) {
        pathIonizationSamples.reserve(triangleGssAvalancheRecords.size());
        for (const auto& record : triangleGssAvalancheRecords) {
            pathIonizationSamples.push_back({
                record.edgeId,
                record.node0,
                record.node1,
                scaling.unitSystem().internalLengthToMeters(record.edgeLength),
                scaling.unitSystem().internalInverseLengthToMInv(record.electronAlpha),
                scaling.unitSystem().internalInverseLengthToMInv(record.holeAlpha),
                scaling.unitSystem().internalElectricFieldToVPerM(std::max(
                    std::abs(record.electronImpactField),
                    std::abs(record.holeImpactField)))});
        }
    } else {
        pathIonizationSamples.reserve(mesh.numEdges());
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            if (!(edge.length > 0.0))
                continue;
            pathIonizationSamples.push_back({
                edgeId,
                edge.n0,
                edge.n1,
                scaling.unitSystem().internalLengthToMeters(edge.length),
                scaling.unitSystem().internalInverseLengthToMInv(
                    0.5 * (electronAlphaAvalanche[edge.n0] +
                           electronAlphaAvalanche[edge.n1])),
                scaling.unitSystem().internalInverseLengthToMInv(
                    0.5 * (holeAlphaAvalanche[edge.n0] +
                           holeAlphaAvalanche[edge.n1])),
                scaling.unitSystem().internalElectricFieldToVPerM(
                    std::abs(sol.psi(static_cast<int>(edge.n1)) -
                             sol.psi(static_cast<int>(edge.n0))) / edge.length)});
        }
    }
    const PathIonizationAnalysis pathIonization = analyzeIonizationPaths(
        static_cast<std::size_t>(N), sol.psi, pathIonizationSamples);
    writer.addNodeScalar("SRHRecombination", srh);
    writer.addNodeScalar("Band2BandGeneration", bandToBandGeneration);
    writer.addNodeScalar("AvalancheGeneration", avalanche);
    writer.addNodeVector("J_n_drift", electronDriftCurrentDensity_A_cm2);
    writer.addNodeVector("J_n_diffusion", electronDiffusionCurrentDensity_A_cm2);
    writer.addNodeVector("J_n_total", electronCurrentDensity_A_cm2);
    writer.addNodeVector("J_p_drift", holeDriftCurrentDensity_A_cm2);
    writer.addNodeVector("J_p_diffusion", holeDiffusionCurrentDensity_A_cm2);
    writer.addNodeVector("J_p_total", holeCurrentDensity_A_cm2);
    writer.addNodeVector("ElectronCurrentDensityVector", electronCurrentDensity_A_cm2);
    writer.addNodeVector("HoleCurrentDensityVector", holeCurrentDensity_A_cm2);
    writer.addNodeVector("TotalCurrentDensityVector", totalCurrentDensity_A_cm2);
    writer.addNodeVector(
        "NodeReconstructedElectronCurrentDensityVector",
        electronCurrentDensity_A_cm2);
    writer.addNodeVector(
        "NodeReconstructedHoleCurrentDensityVector",
        holeCurrentDensity_A_cm2);
    writer.addNodeVector(
        "NodeReconstructedTotalCurrentDensityVector",
        totalCurrentDensity_A_cm2);
    writer.addNodeScalar("ElectronMobility", electronMobility);
    writer.addNodeScalar("HoleMobility", holeMobility);
    writer.addNodeScalar(
        "NodeReconstructedElectronMobility", electronMobility);
    writer.addNodeScalar(
        "NodeReconstructedHoleMobility", holeMobility);
    writer.addNodeScalar("ElectronVelocity", electronVelocity);
    writer.addNodeScalar("HoleVelocity", holeVelocity);
    writer.addNodeScalar("ElectronAlphaAvalanche", electronAlphaAvalanche);
    writer.addNodeScalar("HoleAlphaAvalanche", holeAlphaAvalanche);
    writer.addNodeScalar("ElectronImpactIonizationDrive", electronImpactIonizationDrive);
    writer.addNodeScalar("HoleImpactIonizationDrive", holeImpactIonizationDrive);
    writer.addNodeScalar("ElectronIonIntegral", electronIonIntegral);
    writer.addNodeScalar("HoleIonIntegral", holeIonIntegral);
    writer.addNodeScalar("MeanIonIntegral", meanIonIntegral);
    writer.addNodeScalar(
        "ElectronPathIonIntegral", pathIonization.electronNodeValues);
    writer.addNodeScalar(
        "HolePathIonIntegral", pathIonization.holeNodeValues);
    writer.addNodeScalar(
        "MeanPathIonIntegral", pathIonization.meanNodeValues);
    writer.addNodeScalar(
        "LocalElectronAlphaLengthProxy", electronIonIntegral);
    writer.addNodeScalar(
        "LocalHoleAlphaLengthProxy", holeIonIntegral);
    writer.addNodeScalar(
        "LocalMeanAlphaLengthProxy", meanIonIntegral);
    writer.addNodeScalar("ElectronLowFieldMobility", electronLowFieldMobility);
    writer.addNodeScalar("HoleLowFieldMobility", holeLowFieldMobility);
    writer.addNodeScalar("ElectronHighFieldDrive", electronHighFieldDrive_V_cm);
    writer.addNodeScalar("HoleHighFieldDrive", holeHighFieldDrive_V_cm);
    writer.addNodeScalar("ElectronMobilityLimiter", electronMobilityLimiter);
    writer.addNodeScalar("HoleMobilityLimiter", holeMobilityLimiter);
}

} // namespace vela
