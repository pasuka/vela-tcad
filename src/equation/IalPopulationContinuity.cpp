#include "vela/equation/IalPopulationContinuity.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace vela {

ExtendedIalPopulationContinuity assembleExtendedIalPopulationContinuity(
    const MobilityModelConfig& config, const DeviceMesh& mesh, const DopingModel& doping,
    const ExtendedPopulationState& state,
    const Eigen::Ref<const ExtendedVector>& electronIntrinsic,
    const Eigen::Ref<const ExtendedVector>& holeIntrinsic,
    long double lengthScale_m, long double dopingScale_m3, long double thermalVoltage_V,
    long double temperature_K, bool bandgapNarrowing, bool includeMobilityDerivatives,
    std::optional<IalElementAvalancheOptions> avalanche, std::span<const Index> contactNodes,
    std::optional<ExtendedRecombinationSupport> recombination,
    std::span<const long double> externalCouples,
    const ExtendedIalTransportMeasures* measures,
    ExtendedContactCurrents* allContacts)
{
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::logic_error("Extended IALMob continuity requires the extended80 runner");
    if (!config.ialmobGeometry || config.model != "ialmob" ||
        config.carrierCurrentDiscretization != "scharfetter_gummel_edge" ||
        !mesh.numNodes() || mesh.numNodes() > std::numeric_limits<int>::max() / 3 ||
        doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument("Extended continuity requires the bound SLOT IALMob SG operator");
    const int size = static_cast<int>(mesh.numNodes());
    if (!externalCouples.empty() && externalCouples.size() != mesh.numEdges())
        throw std::invalid_argument("Extended external couples must be indexed by mesh edge");
    if (state.coordinates.size() != 3 * size || electronIntrinsic.size() != size ||
        holeIntrinsic.size() != size || !state.coordinates.allFinite() ||
        !electronIntrinsic.allFinite() || !holeIntrinsic.allFinite() ||
        (state.coordinates.tail(2 * size).array() < 0).any() ||
        (electronIntrinsic.array() < 0).any() || (holeIntrinsic.array() < 0).any())
        throw std::invalid_argument("Invalid authoritative populations for extended continuity");
    for (auto scale : {state.potentialScale_V, state.concentrationScale_m3,
                       lengthScale_m, dopingScale_m3, thermalVoltage_V, temperature_K,
                       static_cast<long double>(config.internalMobilityToM2PerVS)})
        if (!(scale > 0) || !std::isfinite(scale))
            throw std::invalid_argument("Invalid coordinate scale for extended continuity");
    if (recombination) {
        const auto validArea = [](long double area) { return area >= 0 && std::isfinite(area); };
        if (recombination->nodeArea_m2.size() != mesh.numNodes() ||
            recombination->model.bandToBandEnabled() ||
            !std::ranges::all_of(recombination->nodeArea_m2, validArea))
            throw std::invalid_argument("Invalid SRH/Auger support for extended continuity");
    }
    const auto& geometry = *config.ialmobGeometry;
    std::vector<bool> contact(size);
    for (auto node : contactNodes) {
        if (node >= mesh.numNodes() || contact[node])
            throw std::invalid_argument("Invalid port node in extended continuity");
        contact[node] = true;
    }
    std::vector<int> contactOwner(size, -1);
    std::vector<std::array<long double, 2>*> contactSums;
    std::vector<std::array<long double, 2>> contactCorrections;
    if (allContacts) {
        allContacts->clear();
        for (const auto& electrode : mesh.contacts()) {
            const auto id = static_cast<int>(contactSums.size());
            contactSums.push_back(&allContacts->try_emplace(electrode.name).first->second);
            for (auto node : electrode.node_ids) {
                if (contactOwner.at(node) >= 0)
                    throw std::invalid_argument("Overlapping contacts in extended port assembly");
                contactOwner[node] = id;
            }
        }
        contactCorrections.resize(contactSums.size());
    }
    const auto concentration = state.concentrationScale_m3;
    const auto& primary = state.coordinates;
    const auto vertex = [&](Index i) {
        return IalPopulationVertex{
            static_cast<ArithmeticReal>(primary[i]),
            static_cast<ArithmeticReal>(primary[i + size] * concentration),
            static_cast<ArithmeticReal>(primary[i + 2 * size] * concentration),
            static_cast<ArithmeticReal>(electronIntrinsic[i] * concentration),
            static_cast<ArithmeticReal>(holeIntrinsic[i] * concentration)};
    };
    const auto column = [&](const Cell& cell, int component) {
        return static_cast<int>(cell.node_ids[component / 3]) + (component % 3) * size;
    };
    const auto direct = [&](long double differential, int dof) {
        if (dof < size || differential == 0)
            return differential;
        return differential / primary[dof];
    };

    ExtendedIalPopulationContinuity out;
    out.residual = ExtendedVector::Zero(2 * size);
    for (auto& derivative : out.currentDerivative)
        derivative = ExtendedVector::Zero(3 * size);
    std::vector<Eigen::Triplet<long double>> entries;
    entries.reserve(24 * mesh.numEdges());
    const auto add = [&](int row, int col, long double value) {
        if (value != 0)
            entries.emplace_back(row, col, value);
    };
    std::vector<ExtendedIalTransportCell> cells;
    cells.reserve(geometry.cells.size());
    for (const auto& support : geometry.cells) {
        const auto& cell = mesh.getCell(support.cellId);
        if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
            throw std::invalid_argument("Extended IALMob continuity requires Tri3 cells");
        const std::array<IalPopulationVertex, 3> vertices{
            vertex(cell.node_ids[0]), vertex(cell.node_ids[1]), vertex(cell.node_ids[2])};
        cells.push_back(evaluateExtendedIalTransportCell(config, mesh, doping, cell.id, vertices,
            lengthScale_m, dopingScale_m3, state.potentialScale_V, thermalVoltage_V, temperature_K,
            avalanche, measures));
    }

    const auto amperes = -constants::arithmeticConstants<long double>().q * 1e-6L;
    std::array<long double, 2> sumCorrection{};
    const auto addCurrent = [](long double term, long double& sum, long double& correction) {
        const auto next = sum + term;
        correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
        sum = next;
    };
    for (const auto& edge : mesh.edges()) {
        const auto length = mesh.arithmeticEdgeLength(edge.id);
        const auto external = externalCouples.empty() ? -1.L : externalCouples[edge.id];
        if (!std::isfinite(external) || (external < 0 && external != -1.L))
            throw std::invalid_argument("Invalid extended external couple");
        const long double measure = external >= 0 ? external :
            edge.transport_couple >= 0 ? edge.transport_couple : edge.couple;
        if (length <= 0 || measure == 0)
            continue;
        const auto i = static_cast<int>(edge.n0), j = static_cast<int>(edge.n1);
        const auto& support = geometry.edges.at(edge.id);
        const auto geometryFactor = thermalVoltage_V * measure / length;
        const auto portFactor = amperes * (static_cast<int>(contact[i]) - static_cast<int>(contact[j]));
        std::array<long double, 2> mu{
            static_cast<long double>(geometry.constantElectronEdges.at(edge.id)) * config.internalMobilityToM2PerVS,
            static_cast<long double>(geometry.constantHoleEdges.at(edge.id)) * config.internalMobilityToM2PerVS};
        if (measures && !measures->constantMobility_m2_per_Vs.empty())
            mu = measures->constantMobility_m2_per_Vs.at(edge.id);
        const auto edgeWeight = [&](std::size_t k) -> long double {
            return measures ? measures->edgeWeights.at(edge.id).at(k) : support[k].weight;
        };
        for (std::size_t k = 0; k < support.size(); ++k)
            for (int carrier = 0; carrier < 2; ++carrier)
                mu[carrier] += edgeWeight(k) * cells.at(support[k].support).mobility_m2_per_Vs[carrier].value;
        for (int carrier = 0; carrier < 2; ++carrier) {
            if (!(mu[carrier] > 0))
                continue;
            const int block = carrier * size;
            const int population = (carrier + 1) * size;
            const auto& intrinsic = carrier == 0 ? electronIntrinsic : holeIntrinsic;
            const auto signedPotentialScale = carrier == 0 ? state.potentialScale_V : -state.potentialScale_V;
            const auto local = sgBoltzmannPopulationEdge(primary[i], primary[j],
                primary[population + i] * concentration, primary[population + j] * concentration,
                intrinsic[i] * concentration, intrinsic[j] * concentration, signedPotentialScale,
                thermalVoltage_V, mu[carrier] * geometryFactor, bandgapNarrowing);
            out.residual[block + i] += local.flux;
            out.residual[block + j] -= local.flux;
            const auto scatter = [&](int dof, long double derivative) {
                add(block + i, dof, derivative);
                add(block + j, dof, -derivative);
                if (portFactor != 0)
                    out.currentDerivative[carrier][dof] += portFactor * derivative;
            };
            scatter(i, -local.potential1);
            scatter(j, local.potential1);
            scatter(population + i, local.population0 * concentration);
            scatter(population + j, local.population1 * concentration);
            if (includeMobilityDerivatives) {
                const auto mobilitySensitivity = static_cast<long double>(local.flux) / mu[carrier];
                for (std::size_t supportIndex = 0; supportIndex < support.size(); ++supportIndex) {
                    const auto index = support[supportIndex].support;
                    const auto& cell = mesh.getCell(geometry.cells.at(index).cellId);
                    const auto& local = cells.at(index).mobility_m2_per_Vs[carrier];
                    for (int k = 0; k < 9; ++k) {
                        const int col = column(cell, k);
                        scatter(col, direct(mobilitySensitivity * edgeWeight(supportIndex) *
                                            local.derivative[k], col));
                    }
                }
            }
            if (portFactor != 0)
                addCurrent(portFactor * local.flux, out.current_A_per_um[carrier], sumCorrection[carrier]);
            if (allContacts && contactOwner[i] != contactOwner[j]) {
                const auto term = amperes * local.flux;
                const auto first = contactOwner[i], second = contactOwner[j];
                if (first >= 0)
                    addCurrent(term, (*contactSums[first])[carrier], contactCorrections[first][carrier]);
                if (second >= 0)
                    addCurrent(-term, (*contactSums[second])[carrier], contactCorrections[second][carrier]);
            }
        }
    }
    if (avalanche) {
        for (std::size_t k = 0; k < cells.size(); ++k) {
            const auto& cell = mesh.getCell(geometry.cells[k].cellId);
            const auto& source = *cells[k].avalanche;
            for (int j = 0; j < 3; ++j) {
                const auto row = static_cast<int>(cell.node_ids[j]);
                const auto local = source.electron.source_per_m_s[j] + source.hole.source_per_m_s[j];
                for (int carrier = 0; carrier < 2; ++carrier) {
                    const int equation = carrier * size + row;
                    out.residual[equation] -= local.value;
                    for (int i = 0; i < 9; ++i) {
                        const int col = column(cell, i);
                        add(equation, col, -direct(local.derivative[i], col));
                    }
                }
            }
        }
    }
    if (recombination) {
        for (int i = 0; i < size; ++i) {
            const auto area = recombination->nodeArea_m2[i];
            if (area == 0 || electronIntrinsic[i] == 0 || holeIntrinsic[i] == 0)
                continue;
            const auto intrinsic = std::sqrt(electronIntrinsic[i] * holeIntrinsic[i]);
            const auto impurities = doping.dopingOperandsForAssembly(i);
            const auto n = primary[size + i], p = primary[2 * size + i];
            const auto source = recombination->coefficients
                ? recombination->model.boltzmannPopulationLinearization(
                    n, p, intrinsic, impurities[0], impurities[1], *recombination->coefficients)
                : recombination->model.boltzmannPopulationLinearization(n, p, intrinsic, impurities[0], impurities[1]);
            const auto factor = area * concentration;
            const auto value = factor * (source.srh + source.auger);
            const auto dN = factor * (source.dSrhDn + source.dAugerDn);
            const auto dP = factor * (source.dSrhDp + source.dAugerDp);
            for (const auto row : {i, size + i}) {
                out.residual[row] += value;
                add(row, size + i, dN);
                add(row, 2 * size + i, dP);
            }
        }
    }
    out.jacobian.resize(2 * size, 3 * size);
    out.jacobian.setFromTriplets(entries.begin(), entries.end());
    const auto valid = [](long double value) { return std::isfinite(value); };
    for (int carrier = 0; carrier < 2; ++carrier)
        out.current_A_per_um[carrier] += sumCorrection[carrier];
    for (std::size_t i = 0; i < contactSums.size(); ++i)
        for (int carrier = 0; carrier < 2; ++carrier) {
            auto& value = (*contactSums[i])[carrier];
            value += contactCorrections[i][carrier];
            if (!std::isfinite(value))
                throw std::runtime_error("Nonfinite extended electrode current");
        }
    if (!out.residual.allFinite() || !out.currentDerivative[0].allFinite() ||
        !out.currentDerivative[1].allFinite() ||
        !std::isfinite(out.current_A_per_um[0]) || !std::isfinite(out.current_A_per_um[1]) ||
        !std::all_of(out.jacobian.valuePtr(), out.jacobian.valuePtr() + out.jacobian.nonZeros(), valid))
        throw std::runtime_error("Nonfinite extended continuity assembly");
    return out;
}

} // namespace vela
