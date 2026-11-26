#pragma once

#include "vela/equation/IalTransport.h"
#include "vela/physics/RecombinationModel.h"
#include "vela/solver/ExtendedPopulationState.h"
#include <span>

namespace vela {

struct ExtendedIalPopulationContinuity {
    // Unconstrained electron/hole rows in physical particles per metre of
    // device depth per second. Columns are the state's internal [psi,n,p].
    ExtendedVector residual;
    ExtendedSparseMatrix jacobian;
    std::array<long double, 2> current_A_per_um{};
    std::array<ExtendedVector, 2> currentDerivative;
};

using ExtendedContactCurrents = std::map<std::string, std::array<long double, 2>>;

struct ExtendedRecombinationSupport {
    // Coefficients and concentrations use the state and doping's internal
    // concentration units. A zero area excludes that node from recombination.
    const RecombinationModel& model;
    std::span<const long double> nodeArea_m2;
    const ArithmeticBoltzmannSourceCoefficients* coefficients = nullptr;
};

// Assemble SG transport and the selected IALMob avalanche contribution from
// one authoritative population state. Optional SRH/Auger uses the same
// populations and the supplied material-local areas. Boundary replacements
// and Poisson are separate operators. Current uses transport edges only.
ExtendedIalPopulationContinuity assembleExtendedIalPopulationContinuity(
    const MobilityModelConfig& config, const DeviceMesh& mesh, const DopingModel& doping,
    const ExtendedPopulationState& state,
    const Eigen::Ref<const ExtendedVector>& electronIntrinsic,
    const Eigen::Ref<const ExtendedVector>& holeIntrinsic,
    long double lengthScale_m, long double dopingScale_m3, long double thermalVoltage_V, long double temperature_K,
    bool bandgapNarrowing, bool includeMobilityDerivatives,
    std::optional<IalElementAvalancheOptions> avalanche, std::span<const Index> contactNodes,
    std::optional<ExtendedRecombinationSupport> recombination = {},
    std::span<const long double> externalCouples = {},
    const ExtendedIalTransportMeasures* measures = nullptr,
    ExtendedContactCurrents* allContacts = nullptr, bool omitAvalancheDerivatives = false);

} // namespace vela
