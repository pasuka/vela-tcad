#pragma once

#include "vela/equation/ExtendedPopulationDevice.h"
#include "vela/physics/ElementLombardi.h"

namespace vela {

// All coefficients are in SI units. Binding to a deck is separate from the
// operator so residual and Jacobian use identical authoritative operands.
struct ArithmeticPhuMobTransportCoefficients {
    ExtendedPhuMobParameters bulk;
    struct Lombardi {
        ArithmeticReal B, C, N0, N2, lambda, k, delta, A, alpha, eta;
        ArithmeticReal criticalLength, acousticFactor, roughnessFactor;
    };
    struct HighField {
        ArithmeticReal saturationVelocity, beta;
    };
    std::array<Lombardi, 2> surface;
    std::array<HighField, 2> highField;
    ArithmeticIalIonizationCoefficients ionization;
};

struct ArithmeticPhuMobTriangle {
    using Scalar = ArithmeticReal;
    // Coordinates are in metres. Couplings are face length divided by edge
    // length in cyclic vertex order. Measures are material-local areas in m².
    std::array<std::array<Scalar, 2>, 3> coordinates;
    std::array<Scalar, 3> edgeCouple, vertexMeasure, interfaceDistance;
    std::array<Scalar, 2> surfaceNormal;
    std::array<Scalar, 4> partialLayerGradient;
    std::size_t sourceVertex;
    bool touchesContact;
};

struct ArithmeticFermiPopulationVertex {
    ArithmeticReal potential, electrons, holes, electronIntrinsic, holeIntrinsic;
    ArithmeticReal electronDOS, holeDOS, donors, acceptors;
};

struct ArithmeticPhuMobTriangleResult {
    using Dual = detail::BasicTri3LocalForwardDual<ArithmeticReal>;
    // Oriented edge flux and paired avalanche source, in particles/(m s).
    std::array<std::array<Dual, 3>, 2> edgeFlux;
    std::array<Dual, 3> avalanche;
};

// Optional inexact-Newton and source-continuation controls. The tangent
// options preserve the residual at a fixed avalanche strength.
struct ExtendedPhuMobTangentOptions {
    ArithmeticReal avalancheRoundoffFactor = 0;
    bool freezeHighFieldDrive = false;
    ArithmeticReal avalancheStrength = 1;
};

// Outputs carry derivatives in caller [psi,n,p] per vertex, never QF columns.
ArithmeticPhuMobTriangleResult evaluateArithmeticPhuMobTriangle(
    const ArithmeticPhuMobTriangle& geometry, const std::array<ArithmeticFermiPopulationVertex, 3>& state,
    const ArithmeticPhuMobTransportCoefficients& coefficients,
    ArithmeticReal thermalVoltage_V, ArithmeticReal temperature_K,
    bool mobilityContactField, bool ionizationContactField, ArithmeticReal potentialScale_V = 1,
    CarrierStatisticsModel model = CarrierStatisticsModel::FermiDirac,
    const ExtendedPhuMobTangentOptions& tangent = {});


// Support for the SLOT Si/PolySilicon element-box operator. Geometry and
// material role are fixed across Newton evaluations.
struct ExtendedPhuMobMeshSupport {
    struct Cell {
        std::array<Index, 3> nodes;
        ArithmeticPhuMobTriangle geometry;
        bool silicon;
    };
    std::vector<Cell> cells;
    std::vector<ArithmeticReal> siliconArea_m2;
    std::vector<int> electrodeOwner;
    std::vector<std::string> electrodes;
};

ExtendedPhuMobMeshSupport prepareExtendedPhuMobMesh(
    const DeviceMesh& mesh, ArithmeticReal lengthScale_m);

ExtendedIalPopulationContinuity assembleExtendedPhuMobPopulationContinuity(
    const ExtendedPhuMobMeshSupport& support, const DopingModel& doping,
    const ExtendedPopulationState& state, const ExtendedEquilibriumMaterials& material,
    const ArithmeticPhuMobTransportCoefficients& coefficients,
    const std::array<ArithmeticReal, 2>& polyMobility_m2_per_Vs,
    ArithmeticReal thermalVoltage_V, ArithmeticReal temperature_K,
    bool mobilityContactField, bool ionizationContactField,
    const RecombinationModel& recombination, const ArithmeticBoltzmannSourceCoefficients& source,
    std::span<const Index> portNodes = {}, ExtendedContactCurrents* allContacts = nullptr,
    // One flag per prepared cell. Only the avalanche tangent is omitted.
    std::span<const unsigned char> omitAvalancheJacobian = {},
    const ExtendedPhuMobTangentOptions& tangent = {});

} // namespace vela
