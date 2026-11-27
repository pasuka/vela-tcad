#include "vela/equation/PoissonAssembler.h"
#include "vela/equation/ExtendedAveragePoissonGeometry.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/physics/CarrierStatistics.h"
#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace vela {

template<class Scalar>
BasicPoissonAssembler<Scalar>::BasicPoissonAssembler(
    const DeviceMesh&      mesh,
    const MaterialDatabase& matdb,
    const DopingModel&      doping,
    std::vector<RegionFixedChargeSpec> fixedCharges,
    std::vector<InterfaceSheetChargeSpec> sheetCharges,
    std::vector<PoissonNeumannBoundarySpec> neumannBoundaries,
    PoissonScalingSpec scaling, std::optional<Scalar> authoritativePotentialScale_V)
    : mesh_(mesh)
    , matdb_(matdb)
    , doping_(doping)
    , fixedCharges_(std::move(fixedCharges))
    , sheetCharges_(std::move(sheetCharges))
    , neumannBoundaries_(std::move(neumannBoundaries))
    , scaling_(scaling)
    , potentialScale_V_(authoritativePotentialScale_V.value_or(scaling.potentialScale_V))
    , A_(static_cast<int>(mesh.numNodes()),
         static_cast<int>(mesh.numNodes()))
    , b_(Vector::Zero(static_cast<int>(mesh.numNodes())))
{
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument(
            "PoissonAssembler: doping model size does not match mesh node count.");

    if (scaling_.enabled) {
        if (potentialScale_V_ <= 0 || !std::isfinite(potentialScale_V_)) {
            throw std::invalid_argument(
                "PoissonAssembler: scaling.potentialScale_V must be positive and finite.");
        }
        if (scaling_.permittivityReference_F_per_m <= 0.0 ||
            !std::isfinite(scaling_.permittivityReference_F_per_m)) {
            throw std::invalid_argument(
                "PoissonAssembler: scaling.permittivityReference_F_per_m must be positive and finite.");
        }
    }
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

template<class Scalar>
void BasicPoissonAssembler<Scalar>::assemble()
{
    assemble(nullptr);
}

template<class Scalar>
void BasicPoissonAssembler<Scalar>::assemble(const ExtendedAveragePoissonGeometry& geometry)
{
    assemble(&geometry);
}

template<class Scalar>
void BasicPoissonAssembler<Scalar>::assemble(const ExtendedAveragePoissonGeometry* geometry)
{
    const Index N = mesh_.numNodes();
    const auto& edgeCells = detail::buildEdgeCellMap(mesh_);
    if (geometry) {
        if (std::is_same_v<Scalar, Real> ||
            scaling_.poissonChargeVolumePolicy != PoissonChargeVolumePolicy::MaterialLocalBox ||
            geometry->chargeVolumes.size() != N || geometry->edgeStiffness.size() != mesh_.numEdges() ||
            !(geometry->permittivityReference_F_per_m > 0) ||
            !std::isfinite(geometry->permittivityReference_F_per_m) ||
            !fixedCharges_.empty() || !sheetCharges_.empty() || !neumannBoundaries_.empty())
            throw std::invalid_argument("Inconsistent extended AverageBox Poisson assembly");
        chargeVolumes_.clear();
    } else {
        const auto cellMaterials = detail::buildCellMaterials(mesh_, matdb_, constants::T0);
        chargeVolumes_ = detail::computePoissonChargeVolumes(mesh_, cellMaterials, scaling_.poissonChargeVolumePolicy);
    }

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(mesh_.numEdges() * 4 + N);

    b_ = Vector::Zero(static_cast<int>(N));

    // ---- Off-diagonal terms from edge fluxes ----
    for (Index e = 0; e < mesh_.numEdges(); ++e) {
        const Edge& edge = mesh_.getEdge(e);
        const Real  h    = edge.length;
        if (!geometry && h < 1e-30) continue; // degenerate edge guard

        const Scalar G = geometry ? static_cast<Scalar>(geometry->edgeStiffness[e])
            : static_cast<Scalar>(detail::poissonEdgeCoefficient(mesh_, matdb_, edgeCells, e));

        auto i = static_cast<int>(edge.n0);
        auto j = static_cast<int>(edge.n1);

        triplets.emplace_back(i, i,  G);
        triplets.emplace_back(j, j,  G);
        triplets.emplace_back(i, j, -G);
        triplets.emplace_back(j, i, -G);
    }

    A_.setFromTriplets(triplets.begin(), triplets.end());

    // ---- Mobile-free source term: rhs_i = +q * netDoping_i * vol_i ----
    const Real chargeAreaFactor = scaling_.enabled ? scaling_.chargeAreaFactor : 1.0;
    const Real chargeLineFactor = scaling_.enabled ? scaling_.chargeLineFactor : 1.0;
    const Real boundaryLengthFactor = scaling_.enabled ? scaling_.boundaryLengthFactor : 1.0;
    for (Index i = 0; i < N; ++i) {
        if constexpr (std::is_same_v<Scalar, Real>) {
            b_(static_cast<int>(i)) = doping_.hasDecimalNodeDoping()
                ? static_cast<Real>(constants::arithmeticConstants<ArithmeticReal>().q *
                                    doping_.netDopingForAssembly(i) * chargeVolumes_[i] * chargeAreaFactor)
                : constants::q * doping_.netDoping(i) * chargeVolumes_[i] * chargeAreaFactor;
        } else {
            const auto area = geometry ? geometry->chargeVolumes[i] : chargeVolumes_[i];
            b_(static_cast<int>(i)) = constants::arithmeticConstants<Scalar>().q *
                                      doping_.netDopingForAssembly(i) * area * chargeAreaFactor;
        }
    }

    // ---- Region fixed charge: q * fixed_charge_m3 * cell_area / 3 ----
    b_ += detail::computeFixedAndInterfaceChargeRhs(
        mesh_, edgeCells, fixedCharges_, sheetCharges_, "PoissonAssembler",
        chargeAreaFactor, chargeLineFactor).template cast<Scalar>();

    // ---- Neumann boundary conditions ----
    // For each boundary segment defined by a polyline of node IDs, compute the
    // RHS contribution from the normal displacement D.n [C/m^2].
    // For each edge in the polyline: rhs += D_n * edge_length / 2 to each endpoint.
    for (const auto& neumannSpec : neumannBoundaries_) {
        if (neumannSpec.node_ids.size() < 2) continue;

        for (size_t i = 0; i + 1 < neumannSpec.node_ids.size(); ++i) {
            const Index n0 = neumannSpec.node_ids[i];
            const Index n1 = neumannSpec.node_ids[i + 1];

            if (n0 >= N || n1 >= N) {
                throw std::out_of_range(
                    "PoissonAssembler: Neumann boundary node ID out of range.");
            }

            const Node& node0 = mesh_.getNode(n0);
            const Node& node1 = mesh_.getNode(n1);
            const Real dx = node1.x - node0.x;
            const Real dy = node1.y - node0.y;
            const Real edgeLength = std::sqrt(dx * dx + dy * dy);

            if (edgeLength < 1e-30) continue; // Skip degenerate edges

            const Real endpointContribution =
                neumannSpec.normalDisplacement * edgeLength * boundaryLengthFactor * 0.5;
            b_(static_cast<int>(n0)) += endpointContribution;
            b_(static_cast<int>(n1)) += endpointContribution;
        }
    }

    if (scaling_.enabled) {
        const Scalar eps = geometry ? static_cast<Scalar>(geometry->permittivityReference_F_per_m)
                                    : scaling_.permittivityReference_F_per_m;
        const Scalar matrixScale = Scalar(1) / eps;
        const Scalar rhsScale = matrixScale / potentialScale_V_;
        A_ *= matrixScale;
        b_ *= rhsScale;
    }
}

// ---------------------------------------------------------------------------
// Dirichlet boundary conditions
// ---------------------------------------------------------------------------

template<class Scalar>
void BasicPoissonAssembler<Scalar>::applyDirichlet(
    const std::unordered_map<Index, Scalar>& bcs)
{
    if (!scaling_.enabled) {
        detail::applyDirichletBC(A_, b_, bcs);
        return;
    }

    std::unordered_map<Index, Scalar> scaledBcs;
    scaledBcs.reserve(bcs.size());
    for (const auto& [node, value] : bcs)
        scaledBcs.emplace(node, value / potentialScale_V_);
    detail::applyDirichletBC(A_, b_, scaledBcs);
}

template<class Scalar>
Scalar BasicPoissonAssembler<Scalar>::mobileChargeScale(const ExtendedAveragePoissonGeometry* geometry) const
{
    auto charge = constants::arithmeticConstants<Scalar>().q;
    if (scaling_.enabled) {
        const Scalar eps = geometry ? static_cast<Scalar>(geometry->permittivityReference_F_per_m)
                                    : scaling_.permittivityReference_F_per_m;
        charge *= scaling_.chargeAreaFactor /
                  (potentialScale_V_ * eps);
    }
    return charge;
}

template<class Scalar>
Scalar BasicPoissonAssembler<Scalar>::mobileChargeWeight(
    Index node, const ExtendedAveragePoissonGeometry& geometry) const
{
    return mobileChargeScale(&geometry) * static_cast<Scalar>(geometry.chargeVolumes.at(node));
}

template<class Scalar>
typename BasicPoissonAssembler<Scalar>::EquilibriumLinearization
BasicPoissonAssembler<Scalar>::linearizeBoltzmannEquilibrium(
    const Vector& potential, const Vector& electronQf, const Vector& holeQf,
    const Vector& electronIntrinsic, const Vector& holeIntrinsic,
    Scalar thermalVoltage_V, const std::unordered_map<Index, Scalar>& dirichlet,
    DirichletCoordinate coordinate, const ExtendedAveragePoissonGeometry* geometry)
{
    const auto count = static_cast<Eigen::Index>(mesh_.numNodes());
    for (const auto* values : {&electronIntrinsic, &holeIntrinsic})
        if (values->size() != count || !values->allFinite() || (values->array() < 0).any())
            throw std::invalid_argument("Invalid Boltzmann Poisson statistics");
    const auto population = [&](Eigen::Index node, bool hole, Scalar exponent) {
        if (!std::isfinite(exponent))
            throw std::domain_error("Nonfinite Poisson equilibrium exponent");
        const auto maxExponent = Scalar(500);
        const auto intrinsic = hole ? holeIntrinsic[node] : electronIntrinsic[node];
        const auto value = intrinsic * std::exp(std::clamp(exponent, -maxExponent, maxExponent));
        const auto slope = std::abs(exponent) < maxExponent ? value : Scalar(0);
        return std::pair{value, slope};
    };
    return linearizeEquilibrium(potential, electronQf, holeQf, thermalVoltage_V,
                               dirichlet, coordinate, geometry, population);
}

template<class Scalar>
typename BasicPoissonAssembler<Scalar>::EquilibriumLinearization
BasicPoissonAssembler<Scalar>::linearizeFermiEquilibrium(
    const Vector& potential, const Vector& electronQf, const Vector& holeQf,
    const Vector& electronIntrinsic, const Vector& holeIntrinsic,
    const Vector& electronDOS, const Vector& holeDOS,
    Scalar thermalVoltage_V, const std::unordered_map<Index, Scalar>& dirichlet,
    DirichletCoordinate coordinate, const ExtendedAveragePoissonGeometry* geometry,
    CarrierStatisticsModel model)
{
    const auto count = static_cast<Eigen::Index>(mesh_.numNodes());
    for (const auto* values : {&electronIntrinsic, &holeIntrinsic, &electronDOS, &holeDOS})
        if (values->size() != count || !values->allFinite() || (values->array() < 0).any())
            throw std::invalid_argument("Invalid Fermi Poisson statistics");
    const auto population = [&](Eigen::Index node, bool hole, Scalar exponent) {
        const auto intrinsic = hole ? holeIntrinsic[node] : electronIntrinsic[node];
        const auto dos = hole ? holeDOS[node] : electronDOS[node];
        if (!std::isfinite(exponent) || ((intrinsic == 0) != (dos == 0)))
            throw std::domain_error("Invalid Fermi Poisson population arguments");
        if (dos == 0)
            return std::pair{Scalar(0), Scalar(0)};
        if constexpr (std::is_same_v<Scalar, long double>) {
            return extendedFermiPopulation(dos, intrinsic, exponent, model);
        } else {
            const auto eta = exponent + std::log(intrinsic / dos);
            const auto density = scaledFermiDiracHalf(dos, eta, model);
            const auto slope = eta < -700 ? density : dos * fermiDiracHalfDerivative(eta, model);
            return std::pair{density, slope};
        }
    };
    return linearizeEquilibrium(potential, electronQf, holeQf, thermalVoltage_V,
                               dirichlet, coordinate, geometry, population);
}

template<class Scalar>
template<class Population>
typename BasicPoissonAssembler<Scalar>::EquilibriumLinearization
BasicPoissonAssembler<Scalar>::linearizeEquilibrium(
    const Vector& potential, const Vector& electronQf, const Vector& holeQf,
    Scalar thermalVoltage_V, const std::unordered_map<Index, Scalar>& dirichlet,
    DirichletCoordinate coordinate, const ExtendedAveragePoissonGeometry* geometry,
    Population population)
{
    const auto count = static_cast<Eigen::Index>(mesh_.numNodes());
    for (const auto* values : {&potential, &electronQf, &holeQf})
        if (values->size() != count || !values->allFinite())
            throw std::invalid_argument("Invalid Poisson equilibrium state");
    if (!(thermalVoltage_V > 0) || !std::isfinite(thermalVoltage_V))
        throw std::invalid_argument("Invalid Poisson equilibrium thermal voltage");
    const auto potentialScale = scaling_.enabled ? potentialScale_V_ : Scalar(1);
    const auto energyFactor = potentialScale / thermalVoltage_V;
    const auto charge = mobileChargeScale(geometry);

    assemble(geometry);
    EquilibriumLinearization result{A_, b_, Vector(count), Vector(count)};
    // Accumulate physical fluxes from differences in the solver coordinate.
    // Multiplying a large common potential by the full matrix would mix the
    // diagonal summation rounding into a near-equilibrium Poisson residual.
    for (Eigen::Index col = 0; col < A_.outerSize(); ++col)
        for (typename Matrix::InnerIterator entry(A_, col); entry; ++entry)
            if (entry.row() != col)
                result.rhs[entry.row()] += entry.value() * (potential[entry.row()] - potential[col]);
    for (Eigen::Index i = 0; i < count; ++i) {
        const auto [n, dn] = population(i, false, (potential[i] - electronQf[i]) * energyFactor);
        const auto [p, dp] = population(i, true, (holeQf[i] - potential[i]) * energyFactor);
        const Scalar area = geometry ? static_cast<Scalar>(geometry->chargeVolumes[i]) : chargeVolumes_[i];
        const auto weight = charge * area;
        const auto diagonal = weight * (dn + dp) * energyFactor;
        result.electrons[i] = n;
        result.holes[i] = p;
        result.rhs[i] += weight * (p - n);
        if (!std::isfinite(diagonal))
            throw std::domain_error("Nonfinite Poisson equilibrium derivative");
        result.jacobian.coeffRef(i, i) += diagonal;
    }
    if (!result.electrons.allFinite() || !result.holes.allFinite() || !result.rhs.allFinite())
        throw std::domain_error("Nonfinite Poisson equilibrium residual");
    std::unordered_map<Index, Scalar> correction;
    correction.reserve(dirichlet.size());
    for (const auto& [node, value] : dirichlet) {
        if (node >= mesh_.numNodes() || !std::isfinite(value))
            throw std::invalid_argument("Invalid Poisson equilibrium Dirichlet value");
        const auto target = coordinate == DirichletCoordinate::PhysicalVoltage
            ? value / potentialScale : value;
        correction.emplace(node, target - potential[node]);
    }
    detail::applyDirichletBC(result.jacobian, result.rhs, correction);
    if (!result.rhs.allFinite())
        throw std::domain_error("Nonfinite constrained Poisson equilibrium residual");
    return result;
}

template class BasicPoissonAssembler<Real>;
#if defined(VELA_ARITHMETIC_EXTENDED80)
template class BasicPoissonAssembler<ArithmeticReal>;
#endif

} // namespace vela
