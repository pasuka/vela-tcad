#pragma once

#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/ExtendedSparseLu.h"

namespace vela {

struct ExtendedEquilibriumMaterials {
    CarrierStatisticsModel statistics = CarrierStatisticsModel::Boltzmann;
    ExtendedVector intrinsic, electronIntrinsic, holeIntrinsic;
    ExtendedVector intrinsicPotential_V, electronDOS, holeDOS;
};

// Use the bound PAR at the operating temperature for each semiconductor.
// Concentrations are in the same internal units as the input doping model.
ExtendedEquilibriumMaterials buildExtendedEquilibriumMaterials(
    const DeviceMesh& mesh, const MaterialDatabase& database, const DopingModel& doping,
    const BandgapNarrowingConfig& narrowing, long double temperature_K, long double thermalVoltage_V);

struct ExtendedBandgapCoefficients {
    long double referenceDoping, coefficient, smoothing, offset;
};
ExtendedEquilibriumMaterials buildExtendedEquilibriumMaterials(
    const DeviceMesh& mesh, const MaterialDatabase& database, const DopingModel& doping,
    const BandgapNarrowingConfig& narrowing, long double temperature_K, long double thermalVoltage_V,
    const ExtendedBandgapCoefficients& coefficients,
    CarrierStatisticsModel statistics = CarrierStatisticsModel::Boltzmann);

struct ExtendedOhmicPopulation {
    long double potential, electrons, holes;
};

// The equilibrium populations are independent of the common contact bias.
// The potential target is expressed in the caller's state coordinate.
ExtendedOhmicPopulation extendedBoltzmannOhmicTarget(
    long double netDoping, long double intrinsic, long double intrinsicPotential_V,
    long double thermalVoltage_V, long double appliedBias_V,
    long double potentialScale_V, long double potentialOrigin_V);

ExtendedOhmicPopulation extendedEquilibriumOhmicTarget(
    long double netDoping, const ExtendedEquilibriumMaterials& materials, Index node,
    long double thermalVoltage_V);

} // namespace vela
