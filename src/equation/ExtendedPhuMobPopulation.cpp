#include "vela/equation/ExtendedPhuMobPopulation.h"

#include "vela/discretization/ScharfetterGummel.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/mesh/AverageBox.h"
#include "vela/physics/CarrierStatistics.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {
namespace {

using Scalar = ArithmeticReal;
using Dual = ArithmeticPhuMobTriangleResult::Dual;
using Vector = std::array<Dual, 2>;
using TriangleValues = std::array<Dual, 3>;

Dual magnitude(const Vector& vector)
{
    Dual result(std::hypot(vector[0].value, vector[1].value));
    if (result.value > 0)
        for (std::size_t k = 0; k < result.derivative.size(); ++k)
            result.derivative[k] =
                (vector[0].value / result.value) * vector[0].derivative[k] +
                (vector[1].value / result.value) * vector[1].derivative[k];
    return result;
}

Vector gradient(const std::array<Scalar, 4>& gradients, const Dual& first, const Dual& second)
{
    return {first * gradients[0] + second * gradients[2],
            first * gradients[1] + second * gradients[3]};
}

Dual surfaceMobility(const Dual& bulk, const Dual& field, Scalar distance, Scalar impurity,
                     Scalar temperature, const ArithmeticPhuMobTransportCoefficients::Lombardi& p)
{
    if (!(bulk.value > 0) || !std::isfinite(distance) || field.value == 0)
        return bulk;
    const auto F = field.value;
    const auto inverse = elementLombardiInverse(F, distance, impurity, temperature, p);
    const auto mu = 1 / (1 / bulk.value + inverse);
    const auto z = p.C * std::pow(temperature / 300, -p.k) *
        std::pow((impurity + p.N2) / p.N0, p.lambda) * std::pow(F, Scalar(2) / 3);
    const auto den = p.B + z;
    const auto acoustic = (p.B + z / 3) / (den * den);
    const auto roughness = p.A * std::pow(F / 100, p.A) / (p.delta * F) + 3 * F * F / p.eta;
    const auto damping = p.criticalLength > 1 ? Scalar(1) : std::exp(-distance / p.criticalLength);
    const auto bulkSlope = mu / bulk.value;
    const auto fieldSlope = -mu * mu * damping *
        (p.acousticFactor * acoustic + p.roughnessFactor * roughness);
    return Dual(mu) + (bulk - bulk.value) * (bulkSlope * bulkSlope) + (field - F) * fieldSlope;
}

Dual saturatedMobility(const Dual& low, const Dual& field,
                       const ArithmeticPhuMobTransportCoefficients::HighField& p)
{
    constexpr Scalar cutoff = 100;
    const auto mu = elementCanali(low.value, field.value, p, cutoff);
    if (field.value < cutoff || low.value <= 0)
        return low;
    const auto z = low.value * field.value / p.saturationVelocity;
    const auto slope = std::pow(1 + std::pow(z, p.beta), -1 - 1 / p.beta);
    const auto fieldSlope = -low.value * low.value / p.saturationVelocity *
        std::pow(z, p.beta - 1) * slope;
    return Dual(mu) + (low - low.value) * slope + (field - field.value) * fieldSlope;
}

Dual ionization(const Dual& field, const ArithmeticElementIonizationCurve& p, Scalar minimum)
{
    if (field.value <= 0 || field.value < minimum)
        return {};
    const bool high = field.value >= p.switchField_V_per_m;
    const auto a = high ? p.highA_per_m : p.lowA_per_m;
    const auto b = high ? p.highB_V_per_m : p.lowB_V_per_m;
    return Dual(a) * detail::dualExp(-Dual(b) / field);
}

} // namespace

ArithmeticPhuMobTriangleResult evaluateArithmeticPhuMobTriangle(
    const ArithmeticPhuMobTriangle& geometry, const std::array<ArithmeticFermiPopulationVertex, 3>& state,
    const ArithmeticPhuMobTransportCoefficients& coefficients,
    Scalar vt, Scalar temperature, bool mobilityContactField, bool ionizationContactField, Scalar potentialScale_V,
    CarrierStatisticsModel model, const ExtendedPhuMobTangentOptions& tangent)
{
    if (!(vt > 0 && temperature > 0 && potentialScale_V > 0) || geometry.sourceVertex >= 3)
        throw std::invalid_argument("Invalid extended PhuMob triangle");
    if (!(tangent.avalancheRoundoffFactor >= 0) || !std::isfinite(tangent.avalancheRoundoffFactor))
        throw std::invalid_argument("Invalid avalanche roundoff tangent scale");
    if (!(tangent.avalancheStrength >= 0 && tangent.avalancheStrength <= 1))
        throw std::invalid_argument("Invalid avalanche continuation strength");
    const auto x1 = geometry.coordinates[1][0] - geometry.coordinates[0][0];
    const auto y1 = geometry.coordinates[1][1] - geometry.coordinates[0][1];
    const auto x2 = geometry.coordinates[2][0] - geometry.coordinates[0][0];
    const auto y2 = geometry.coordinates[2][1] - geometry.coordinates[0][1];
    const auto determinant = x1 * y2 - x2 * y1;
    if (determinant == 0 || !std::isfinite(determinant))
        throw std::invalid_argument("Degenerate extended PhuMob triangle");
    const std::array fullGradient{y2 / determinant, -x2 / determinant, -y1 / determinant, x1 / determinant};
    TriangleValues psi, electrons, holes;
    for (std::size_t k = 0; k < 3; ++k) {
        psi[k] = Dual::variable(state[k].potential, 3 * k);
        electrons[k] = Dual::variable(state[k].electrons, 3 * k + 1);
        holes[k] = Dual::variable(state[k].holes, 3 * k + 2);
    }
    const auto electric = gradient(fullGradient, (psi[1] - psi[0]) * potentialScale_V,
                                   (psi[2] - psi[0]) * potentialScale_V);
    const auto electricMagnitude = magnitude(electric);
    const auto normal = detail::dualAbs(electric[0] * geometry.surfaceNormal[0] +
                                        electric[1] * geometry.surfaceNormal[1]);
    const auto cellMeasure = geometry.vertexMeasure[0] + geometry.vertexMeasure[1] + geometry.vertexMeasure[2];
    if (cellMeasure == 0 && std::ranges::all_of(geometry.vertexMeasure, [](auto value) { return value == 0; }))
        return {};
    if (!(cellMeasure > 0) || !std::isfinite(cellMeasure))
        throw std::invalid_argument("Invalid extended PhuMob control volumes");

    ArithmeticPhuMobTriangleResult out;
    for (std::size_t carrier = 0; carrier < 2; ++carrier) {
        const bool electron = carrier == 0;
        const auto species = electron ? CarrierType::Electron : CarrierType::Hole;
        const auto& density = electron ? electrons : holes;
        const auto sign = electron ? Scalar(1) : Scalar(-1);
        std::array<Scalar, 3> intrinsic, dos;
        TriangleValues correction;
        for (std::size_t k = 0; k < 3; ++k) {
            intrinsic[k] = electron ? state[k].electronIntrinsic : state[k].holeIntrinsic;
            dos[k] = electron ? state[k].electronDOS : state[k].holeDOS;
            if (!(density[k].value > 0 && intrinsic[k] > 0 && dos[k] > 0))
                throw std::domain_error("Extended PhuMob requires positive populations and DOS");
            const auto factor = extendedFermiDegeneracy(
                static_cast<long double>(density[k].value) / dos[k], model);
            correction[k] = Dual(static_cast<Scalar>(factor.logCorrection)) +
                (density[k] - density[k].value) * static_cast<Scalar>(factor.logSlope / density[k].value);
        }
        const auto drop = [&](std::size_t i, std::size_t j) {
            Dual logRatio(arithmeticPopulationLogRatio(density[i].value, density[j].value));
            logRatio.derivative[3 * i + carrier + 1] += 1 / density[i].value;
            logRatio.derivative[3 * j + carrier + 1] -= 1 / density[j].value;
            const auto classical = (psi[j] - psi[i]) * (sign * potentialScale_V / vt) +
                arithmeticPopulationLogRatio(intrinsic[j], intrinsic[i]);
            return (classical + logRatio + (correction[i] - correction[j])) * vt;
        };
        const auto drive = magnitude(gradient(geometry.partialLayerGradient, drop(0, 1), drop(0, 2)));
        const auto hfsField = geometry.touchesContact && mobilityContactField ? electricMagnitude : drive;
        const auto ionizationField = geometry.touchesContact && ionizationContactField ? electricMagnitude : drive;
        const auto flux = [&](std::size_t i, std::size_t j) {
            const auto value = sgFermiPopulationEdge(psi[i].value, psi[j].value,
                density[i].value, density[j].value, intrinsic[i], intrinsic[j], dos[i], dos[j],
                sign * potentialScale_V, vt, vt, model);
            Dual result(value.flux);
            result.derivative[3 * i] = -value.potential1;
            result.derivative[3 * j] = value.potential1;
            result.derivative[3 * i + carrier + 1] = value.population0;
            result.derivative[3 * j + carrier + 1] = value.population1;
            return result;
        };

        TriangleValues mobility;
        Dual mean;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto& vertex = state[k];
            const ExtendedPhuMobScalarState local{
                vertex.donors, vertex.acceptors, vertex.electrons, vertex.holes, temperature};
            const auto value = evaluateExtendedPhuMobScalar(species, local, coefficients.bulk).mobility;
            const auto slopes = evaluateExtendedPhuMobLogDensityDerivatives(species, local, coefficients.bulk);
            Dual bulk(static_cast<Scalar>(value));
            bulk.derivative[3 * k + 1] = static_cast<Scalar>(slopes.electrons / vertex.electrons);
            bulk.derivative[3 * k + 2] = static_cast<Scalar>(slopes.holes / vertex.holes);
            mobility[k] = saturatedMobility(surfaceMobility(bulk, normal, geometry.interfaceDistance[k],
                vertex.donors + vertex.acceptors, temperature, coefficients.surface[carrier]),
                tangent.freezeHighFieldDrive ? Dual(hfsField.value) : hfsField, coefficients.highField[carrier]);
            mean = mean + mobility[k] * (geometry.vertexMeasure[k] / cellMeasure);
        }
        TriangleValues edge;
        for (std::size_t k = 0; k < 3; ++k) {
            edge[k] = flux(k, (k + 1) % 3);
            out.edgeFlux[carrier][k] = edge[k] * mean * geometry.edgeCouple[k];
        }
        if (tangent.avalancheStrength == 0)
            continue;
        const auto best = geometry.sourceVertex;
        TriangleValues fromSource;
        fromSource[(best + 1) % 3] = edge[best];
        fromSource[(best + 2) % 3] = -edge[(best + 2) % 3];
        const auto& stencil = geometry.touchesContact ? fullGradient : geometry.partialLayerGradient;
        const auto j = gradient(stencil, fromSource[1] - fromSource[0], fromSource[2] - fromSource[0]);
        const auto alpha = ionization(ionizationField,
            electron ? coefficients.ionization.electron : coefficients.ionization.hole,
            coefficients.ionization.minimumField_V_per_m);
        auto current = magnitude(j);
        // The current direction can be determined by the subtraction of two
        // almost equal SG terms. Attenuate only the norm tangent at that scale.
        if (tangent.avalancheRoundoffFactor > 0 && current.value > 0) {
            const auto edgeRoundoff = [&](std::size_t i, std::size_t k) {
                const auto drift = (psi[k].value - psi[i].value) * (sign * potentialScale_V / vt) +
                    arithmeticPopulationLogRatio(intrinsic[k], intrinsic[i]) +
                    (correction[i].value - correction[k].value);
                return (std::numeric_limits<Scalar>::epsilon() * vt) *
                    std::max(density[i].value, density[k].value) * (1 + std::abs(drift));
            };
            const auto next = (best + 1) % 3;
            const auto prev = (best + 2) % 3;
            std::array<Scalar, 3> nodalRadius{};
            nodalRadius[next] = edgeRoundoff(best, next);
            nodalRadius[prev] = edgeRoundoff(best, prev);
            const auto radius = [&](std::size_t axis) {
                return std::abs(stencil[axis] + stencil[axis + 2]) * nodalRadius[0] +
                    std::abs(stencil[axis]) * nodalRadius[1] +
                    std::abs(stencil[axis + 2]) * nodalRadius[2];
            };
            const auto uncertainty = tangent.avalancheRoundoffFactor * std::hypot(radius(0), radius(1));
            const auto weight = current.value / std::hypot(current.value, uncertainty);
            for (auto& partial : current.derivative)
                partial *= weight;
        }
        const auto source = alpha * mobility[best] * current * tangent.avalancheStrength;
        for (std::size_t k = 0; k < 3; ++k)
            out.avalanche[k] = out.avalanche[k] + source * geometry.vertexMeasure[k];
    }
    return out;
}


ExtendedPhuMobMeshSupport prepareExtendedPhuMobMesh(
    const DeviceMesh& mesh, Scalar lengthScale_m)
{
    if (!(lengthScale_m > 0) || !std::isfinite(lengthScale_m) ||
        mesh.poissonPermittivityPolicy() != BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial)
        throw std::invalid_argument("Extended PhuMob requires cell-material AverageBox geometry");

    ExtendedPhuMobMeshSupport result;
    result.siliconArea_m2.resize(mesh.numNodes());
    result.electrodeOwner.assign(mesh.numNodes(), -1);
    for (const auto& contact : mesh.contacts()) {
        for (auto node : contact.node_ids) {
            auto& owner = result.electrodeOwner.at(node);
            if (owner >= 0)
                throw std::invalid_argument("Overlapping extended PhuMob contacts");
            owner = static_cast<int>(result.electrodes.size());
        }
        result.electrodes.push_back(contact.name);
    }

    using Point = std::array<Scalar, 2>;
    const auto coordinate = [&](Index node) {
        const auto [x, y] = mesh.arithmeticCoordinates(node);
        return Point{x, y};
    };
    const auto delta = [](const Point& a, const Point& b) { return Point{b[0] - a[0], b[1] - a[1]}; };
    const auto dot = [](const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1]; };
    const auto adjacency = detail::buildEdgeCellMap(mesh);
    const auto edgeIds = detail::buildCellEdgeMap(adjacency, mesh);
    const auto faces = buildArithmeticAverageCellBox(mesh);
    std::vector<int> role(mesh.numCells());
    for (const auto& cell : mesh.cells()) {
        const auto& name = mesh.getRegion(cell.region_id).material;
        role[cell.id] = name == "Si" ? 1 : (name == "PolySi" || name == "PolySilicon") ? 2 : 0;
    }
    std::vector<std::array<Point, 2>> interfaces;
    std::vector<bool> boundary(mesh.numEdges());
    for (const auto& edge : mesh.edges()) {
        const auto& adjacent = adjacency.at(edge.id);
        const auto transport = std::ranges::count_if(adjacent, [&](auto cell) { return role[cell] != 0; });
        boundary[edge.id] = transport == 1;
        if (adjacent.size() == 2 && transport == 1)
            interfaces.push_back({coordinate(edge.n0), coordinate(edge.n1)});
    }
    std::vector<Scalar> distance(mesh.numNodes(), std::numeric_limits<Scalar>::infinity());
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        const auto point = coordinate(node);
        for (const auto& segment : interfaces) {
            const auto direction = delta(segment[0], segment[1]);
            const auto position = delta(segment[0], point);
            const auto projection = std::clamp(dot(position, direction) / dot(direction, direction), Scalar(0), Scalar(1));
            distance[node] = std::min(distance[node], std::hypot(position[0] - projection * direction[0],
                                                                  position[1] - projection * direction[1]));
        }
    }
    for (const auto& cell : mesh.cells()) {
        if (!role[cell.id])
            continue;
        if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
            throw std::invalid_argument("Extended PhuMob requires Tri3 cells");
        ExtendedPhuMobMeshSupport::Cell entry{};
        entry.silicon = role[cell.id] == 1;
        auto& geom = entry.geometry;
        geom.edgeCouple = faces.at(cell.id);
        for (std::size_t k = 0; k < 3; ++k) {
            const auto i = cell.node_ids[k], j = cell.node_ids[(k + 1) % 3];
            entry.nodes[k] = i;
            const auto position = coordinate(i);
            geom.coordinates[k] = {(position[0] - coordinate(cell.node_ids[0])[0]) * lengthScale_m,
                (position[1] - coordinate(cell.node_ids[0])[1]) * lengthScale_m};
            geom.interfaceDistance[k] = distance[i] * lengthScale_m;
            const auto edge = delta(position, coordinate(j));
            const auto area = geom.edgeCouple[k] * dot(edge, edge) * lengthScale_m * lengthScale_m / 4;
            geom.vertexMeasure[k] += area;
            geom.vertexMeasure[(k + 1) % 3] += area;
            geom.touchesContact |= result.electrodeOwner[i] >= 0;
        }
        for (std::size_t k = 0; k < 3; ++k) {
            if (entry.silicon)
                result.siliconArea_m2[entry.nodes[k]] += geom.vertexMeasure[k];
            const auto best = geom.sourceVertex;
            if (geom.vertexMeasure[k] > geom.vertexMeasure[best] ||
                (geom.vertexMeasure[k] == geom.vertexMeasure[best] && entry.nodes[k] < entry.nodes[best]))
                geom.sourceVertex = k;
        }

        const auto a = delta(coordinate(entry.nodes[0]), coordinate(entry.nodes[1]));
        const auto b = delta(coordinate(entry.nodes[0]), coordinate(entry.nodes[2]));
        const auto determinant = (a[0] * b[1] - a[1] * b[0]) * lengthScale_m;
        if (determinant == 0 || !std::isfinite(determinant))
            throw std::invalid_argument("Degenerate extended PhuMob geometry");
        Point first{b[1] / determinant, -b[0] / determinant};
        Point second{-a[1] / determinant, a[0] / determinant};
        if (std::ranges::all_of(geom.interfaceDistance, [](auto x) { return std::isfinite(x); })) {
            const auto d1 = (distance[entry.nodes[1]] - distance[entry.nodes[0]]) * lengthScale_m;
            const auto d2 = (distance[entry.nodes[2]] - distance[entry.nodes[0]]) * lengthScale_m;
            geom.surfaceNormal = {first[0] * d1 + second[0] * d2, first[1] * d1 + second[1] * d2};
        }
        std::optional<Point> tangent;
        bool corner = false;
        if (!geom.touchesContact)
            for (auto edgeId : edgeIds.at(cell.id)) {
                if (!boundary[edgeId]) continue;
                const auto& edge = mesh.getEdge(edgeId);
                auto direction = delta(coordinate(edge.n0), coordinate(edge.n1));
                const auto length = std::hypot(direction[0], direction[1]);
                direction = {direction[0] / length, direction[1] / length};
                if (tangent && std::abs((*tangent)[0] * direction[1] - (*tangent)[1] * direction[0]) > 1e-12L)
                    corner = true;
                tangent = direction;
            }
        if (corner) {
            first = second = {};
        } else if (tangent) {
            const auto alongFirst = dot(first, *tangent), alongSecond = dot(second, *tangent);
            first = {alongFirst * (*tangent)[0], alongFirst * (*tangent)[1]};
            second = {alongSecond * (*tangent)[0], alongSecond * (*tangent)[1]};
        }
        geom.partialLayerGradient = {first[0], first[1], second[0], second[1]};
        result.cells.push_back(entry);
    }
    return result;
}

ExtendedIalPopulationContinuity assembleExtendedPhuMobPopulationContinuity(
    const ExtendedPhuMobMeshSupport& support, const DopingModel& doping,
    const ExtendedPopulationState& state, const ExtendedEquilibriumMaterials& material,
    const ArithmeticPhuMobTransportCoefficients& coefficients, const std::array<Scalar, 2>& polyMobility,
    Scalar vt, Scalar temperature, bool mobilityContactField, bool ionizationContactField,
    const RecombinationModel& recombination, const ArithmeticBoltzmannSourceCoefficients& source,
    std::span<const Index> portNodes, ExtendedContactCurrents* allContacts,
    std::span<const unsigned char> omitAvalancheJacobian,
    const ExtendedPhuMobTangentOptions& tangent)
{
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::logic_error("Extended PhuMob assembly requires the extended80 runner");
    const auto size = static_cast<int>(doping.numNodes());
    if (size <= 0 || doping.numNodes() > std::numeric_limits<int>::max() / 3 ||
        std::ssize(support.siliconArea_m2) != size || state.coordinates.size() != 3 * size ||
        !usesFermiDirac(material.statistics) ||
        material.electronDOS.size() != size || material.holeDOS.size() != size ||
        material.electronIntrinsic.size() != size || material.holeIntrinsic.size() != size)
        throw std::invalid_argument("Invalid extended PhuMob population state");
    if (!omitAvalancheJacobian.empty() && omitAvalancheJacobian.size() != support.cells.size())
        throw std::invalid_argument("Avalanche Jacobian mask must match extended PhuMob cells");
    // The shared device assembler supplies TCAD populations and Poisson rows.
    // Convert only the local transport operands and convert tangents back.
    constexpr Scalar tcadToSI = 1000000;
    if (state.concentrationScale_m3 != tcadToSI)
        throw std::invalid_argument("Extended PhuMob uses the TCAD population view");
    std::vector<bool> port(size);
    for (auto node : portNodes)
        port.at(node) = true;
    ExtendedIalPopulationContinuity out;
    out.residual = ExtendedVector::Zero(2 * size);
    for (auto& derivative : out.currentDerivative)
        derivative = ExtendedVector::Zero(3 * size);
    std::vector<Eigen::Triplet<long double>> entries;
    entries.reserve(support.cells.size() * 50);
    const auto add = [&](int row, int col, long double value) {
        if (value != 0)
            entries.emplace_back(row, col, value);
    };

    const long double currentScale = -constants::arithmeticConstants<long double>().q * 1e-6L;
    std::array<long double, 2> portCorrection{};
    std::vector<std::array<long double, 2>> current(support.electrodes.size()), currentCorrection(current.size());
    const auto accumulate = [](long double term, long double& sum, long double& correction) {
        const auto next = sum + term;
        correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
        sum = next;
    };
    for (std::size_t cellIndex = 0; cellIndex < support.cells.size(); ++cellIndex) {
        const auto& cell = support.cells[cellIndex];
        const bool omitAvalanche = !omitAvalancheJacobian.empty() && omitAvalancheJacobian[cellIndex];
        std::array<ArithmeticFermiPopulationVertex, 3> point;
        for (int k = 0; k < 3; ++k) {
            const auto i = cell.nodes[k];
            const auto impurity = doping.dopingOperandsForAssembly(i);
            point[k] = {static_cast<Scalar>(state.coordinates[i]),
                static_cast<Scalar>(state.coordinates[i + size] * tcadToSI),
                static_cast<Scalar>(state.coordinates[i + 2 * size] * tcadToSI),
                static_cast<Scalar>(material.electronIntrinsic[i] * tcadToSI), static_cast<Scalar>(material.holeIntrinsic[i] * tcadToSI),
                static_cast<Scalar>(material.electronDOS[i] * tcadToSI), static_cast<Scalar>(material.holeDOS[i] * tcadToSI),
                impurity[0] * tcadToSI, impurity[1] * tcadToSI};
        }
        ArithmeticPhuMobTriangleResult local{};
        if (cell.silicon) {
            local = evaluateArithmeticPhuMobTriangle(cell.geometry, point, coefficients, vt, temperature,
                mobilityContactField, ionizationContactField,
                static_cast<Scalar>(state.potentialScale_V), material.statistics, tangent);
        } else {
            for (int c = 0; c < 2; ++c)
                for (int k = 0; k < 3; ++k) {
                    const int j = (k + 1) % 3;
                    const auto& a = point[k];
                    const auto& b = point[j];
                    const auto value = c == 0
                        ? sgFermiPopulationEdge(a.potential, b.potential, a.electrons, b.electrons,
                            a.electronIntrinsic, b.electronIntrinsic, a.electronDOS, b.electronDOS,
                            static_cast<Scalar>(state.potentialScale_V), vt,
                            vt * polyMobility[c] * cell.geometry.edgeCouple[k], material.statistics)
                        : sgFermiPopulationEdge(a.potential, b.potential, a.holes, b.holes,
                            a.holeIntrinsic, b.holeIntrinsic, a.holeDOS, b.holeDOS,
                            -static_cast<Scalar>(state.potentialScale_V), vt,
                            vt * polyMobility[c] * cell.geometry.edgeCouple[k], material.statistics);
                    auto& flux = local.edgeFlux[c][k];
                    flux.value = value.flux;
                    flux.derivative[3 * k] = -value.potential1;
                    flux.derivative[3 * j] = value.potential1;
                    flux.derivative[3 * k + c + 1] = value.population0;
                    flux.derivative[3 * j + c + 1] = value.population1;
                }
        }
        for (int c = 0; c < 2; ++c)
            for (int k = 0; k < 3; ++k) {
                const auto i = cell.nodes[k], j = cell.nodes[(k + 1) % 3];
                const auto row = c * size + static_cast<int>(i);
                const auto opposite = c * size + static_cast<int>(j);
                const auto& edge = local.edgeFlux[c][k];
                const auto& avalanche = local.avalanche[k];
                out.residual[row] += edge.value - avalanche.value;
                out.residual[opposite] -= edge.value;
                const auto mask = static_cast<int>(port[i]) - static_cast<int>(port[j]);
                for (int d = 0; d < 9; ++d) {
                    const auto col = static_cast<int>(cell.nodes[d / 3]) + (d % 3) * size;
                    const auto scale = d % 3 == 0 ? Scalar(1) : tcadToSI;
                    const auto partial = edge.derivative[d] * scale;
                    add(row, col, partial - (omitAvalanche ? 0 : avalanche.derivative[d] * scale));
                    add(opposite, col, -partial);
                    if (mask)
                        out.currentDerivative[c][col] += currentScale * mask * partial;
                }
                if (mask)
                    accumulate(currentScale * mask * edge.value, out.current_A_per_um[c], portCorrection[c]);
                if (allContacts) {
                    const auto first = support.electrodeOwner[i], second = support.electrodeOwner[j];
                    if (first == second)
                        continue;
                    if (first >= 0)
                        accumulate(currentScale * edge.value, current[first][c], currentCorrection[first][c]);
                    if (second >= 0)
                        accumulate(-currentScale * edge.value, current[second][c], currentCorrection[second][c]);
                }
            }
    }
    for (int i = 0; i < size; ++i) {
        const auto area = support.siliconArea_m2[i];
        if (area == 0)
            continue;
        const auto impurity = doping.dopingOperandsForAssembly(i);
        const auto rates = recombination.fermiPopulationLinearization(
            static_cast<Scalar>(state.coordinates[size + i]), static_cast<Scalar>(state.coordinates[2 * size + i]),
            material.intrinsic[i], material.electronDOS[i], material.holeDOS[i],
            impurity[0], impurity[1], source, material.statistics);
        for (const auto row : {i, size + i}) {
            out.residual[row] += area * tcadToSI * (rates.srh + rates.auger);
            add(row, size + i, area * tcadToSI * (rates.dSrhDn + rates.dAugerDn));
            add(row, 2 * size + i, area * tcadToSI * (rates.dSrhDp + rates.dAugerDp));
        }
    }
    for (int c = 0; c < 2; ++c)
        out.current_A_per_um[c] += portCorrection[c];
    if (allContacts) {
        allContacts->clear();
        for (std::size_t k = 0; k < current.size(); ++k) {
            auto values = current[k];
            for (int c = 0; c < 2; ++c)
                values[c] += currentCorrection[k][c];
            allContacts->emplace(support.electrodes[k], values);
        }
    }
    out.jacobian.resize(2 * size, 3 * size);
    out.jacobian.setFromTriplets(entries.begin(), entries.end());
    if (!out.residual.allFinite() || !out.currentDerivative[0].allFinite() || !out.currentDerivative[1].allFinite() ||
        !std::ranges::all_of(out.current_A_per_um, [](auto v) { return std::isfinite(v); }) ||
        !std::all_of(out.jacobian.valuePtr(), out.jacobian.valuePtr() + out.jacobian.nonZeros(),
                     [](auto v) { return std::isfinite(v); }))
        throw std::runtime_error("Nonfinite extended PhuMob assembly");
    return out;
}

} // namespace vela
