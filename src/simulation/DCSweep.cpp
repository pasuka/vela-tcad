#include "vela/simulation/DCSweep.h"
#include "vela/boundary/BoundaryCondition.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/simulation/DCSweepPredictor.h"
#include "vela/simulation/DCSweepStepControl.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/io/CSVWriter.h"
#include "vela/io/CsvUtils.h"
#include "vela/io/DDSolutionCsv.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/RecombinationModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/post/ElectricFieldDiagnostics.h"
#include "vela/post/TerminalCharge.h"
#include "vela/post/StoredCharge.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/solver/SolutionValidation.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace vela {

namespace {

enum class SolverMethod {
    Gummel,
    Newton,
    GummelNewton,
};

struct HybridHandoffConfig {
    bool fallbackToGummelOnNewtonFailure = false;
    bool requireGummelConvergence = true;
    int gummelMaxIter = -1;
    int newtonMaxIter = -1;
};

struct SweepRecombinationDiagnostics {
    Real maxAbsRate_m3_per_s = 0.0;
    Real meanAbsRate_m3_per_s = 0.0;
    Real maxCarrierProductRatio = 0.0;
};

struct SweepTransportDiagnostics {
    Real meanElectronMobility_m2_V_s = 0.0;
    Real meanHoleMobility_m2_V_s = 0.0;
    Real minElectronMobility_m2_V_s = 0.0;
    Real minHoleMobility_m2_V_s = 0.0;
    Real maxElectricField_V_per_cm = 0.0;
    Real meanElectronQfGradient_V_per_cm = 0.0;
    Real meanHoleQfGradient_V_per_cm = 0.0;
    Real meanElectronHighFieldDrive_V_per_cm = 0.0;
    Real meanHoleHighFieldDrive_V_per_cm = 0.0;
    Real minElectronMobilityLimiter = 0.0;
    Real minHoleMobilityLimiter = 0.0;
    Real meanElectronMobilityLimiter = 0.0;
    Real meanHoleMobilityLimiter = 0.0;
};

struct ContinuityBalanceDiagnosticRow {
    std::string contact;
    std::string carrier;
    Index contactNode = 0;
    Index interiorNode = 0;
    Index contactEdgeId = 0;
    Real contactEdgeFlux = 0.0;
    Real neighborEdgeFlux = 0.0;
    Real recombinationTerm = 0.0;
    Real continuityResidual = 0.0;
    Real interiorVolume_m2 = 0.0;
    Real qfContact_V = 0.0;
    Real qfInterior_V = 0.0;
    Real carrierDensityInterior_m3 = 0.0;
};

std::string classifySgAvalancheEdge(
    const DeviceMesh& mesh,
    const std::vector<std::vector<Index>>& edgeCells,
    const detail::SgEdgeCurrentAvalancheSourceRecord& record)
{
    bool node0Contact = false;
    bool node1Contact = false;
    for (Index contactId = 0; contactId < mesh.numContacts(); ++contactId) {
        const Contact& contact = mesh.getContact(contactId);
        node0Contact = node0Contact ||
            std::find(contact.node_ids.begin(), contact.node_ids.end(), record.node0) != contact.node_ids.end();
        node1Contact = node1Contact ||
            std::find(contact.node_ids.begin(), contact.node_ids.end(), record.node1) != contact.node_ids.end();
    }
    if (node0Contact || node1Contact)
        return "contact_edge";
    if (record.edgeId < edgeCells.size() && edgeCells[record.edgeId].size() < 2)
        return "boundary";
    return "interior_bulk";
}

std::string formatReal(Real value)
{
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

std::string biasToken(Real bias)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << std::abs(bias);
    std::string token = out.str();
    std::replace(token.begin(), token.end(), '.', 'p');
    return (bias < 0.0 ? "m" : "") + token;
}

std::string formatIndexOrMinusOne(Index value)
{
    if (value == std::numeric_limits<Index>::max())
        return "-1";
    return std::to_string(value);
}


std::string sanitizedColumnToken(std::string value)
{
    for (char& ch : value) {
        const unsigned char uch = static_cast<unsigned char>(ch);
        if (!std::isalnum(uch))
            ch = '_';
        else
            ch = static_cast<char>(std::tolower(uch));
    }
    if (value.empty())
        value = "terminal";
    return value;
}

std::string terminalChargeName(const TerminalChargeConfig& cfg, std::size_t index)
{
    if (!cfg.name.empty())
        return sanitizedColumnToken(cfg.name);
    if (!cfg.contact.empty())
        return sanitizedColumnToken(cfg.contact);
    return "terminal" + std::to_string(index + 1);
}

std::string capacitanceMnemonic(const std::string& sweepContact, const std::string& terminalName)
{
    return "C" + sanitizedColumnToken(sweepContact) + "_" + sanitizedColumnToken(terminalName);
}

Real perMeterToPerMicron(Real value)
{
    return value / 1.0e6;
}
Real voltsPerMeterToVoltsPerCm(Real value)
{
    return value / 100.0;
}

Real terminalCurrentConsistencyRatio(const ContactCurrentResult& current)
{
    constexpr Real floor = 1.0e-300;
    const Real componentMagnitude =
        std::abs(current.electronDriftCurrent) +
        std::abs(current.electronDiffusionCurrent) +
        std::abs(current.holeDriftCurrent) +
        std::abs(current.holeDiffusionCurrent);
    if (componentMagnitude <= floor)
        return 1.0;
    return std::abs(current.totalCurrent) / componentMagnitude;
}

std::vector<Real> buildEffectiveIntrinsicDensityVector(const DeviceMesh& mesh,
                                                       const MaterialDatabase& matdb,
                                                       const DopingModel& doping,
                                                       Real temperature_K,
                                                       const BandgapNarrowingConfig& bgnCfg)
{
    const Index nodeCount = mesh.numNodes();
    std::vector<Real> ni(nodeCount, 0.0);
    std::vector<bool> seen(nodeCount, false);

    for (Index c = 0; c < mesh.numCells(); ++c) {
        const Cell& cell = mesh.getCell(c);
        const Region& region = mesh.getRegion(cell.region_id);
        Real niMaterial = 0.0;
        if (matdb.hasMaterial(region.material))
            niMaterial = matdb.getMaterial(region.material, temperature_K).ni;
        for (Index nodeId : cell.node_ids) {
            if (!seen[nodeId]) {
                ni[nodeId] = niMaterial;
                seen[nodeId] = true;
            }
        }
    }

    const Real Vt = constants::kb * temperature_K / constants::q;
    const std::unique_ptr<BandgapNarrowing> bgn = makeBandgapNarrowingModel(bgnCfg);
    for (Index i = 0; i < nodeCount; ++i) {
        const Real deltaEg = bgn->deltaEg(doping.totalImpurity(i), 0.0, 0.0);
        ni[i] = effectiveIntrinsicDensity(ni[i], Vt, deltaEg);
    }
    return ni;
}

SweepRecombinationDiagnostics computeSweepRecombinationDiagnostics(
    const DDSolution& sol,
    const std::vector<Real>& effectiveNi,
    const RecombinationModelConfig& recombinationCfg)
{
    SweepRecombinationDiagnostics diagnostics;
    const Index nodeCount = static_cast<Index>(effectiveNi.size());
    if (nodeCount == 0)
        return diagnostics;

    const RecombinationModel recombination(recombinationCfg);
    Real absSum = 0.0;
    for (Index i = 0; i < nodeCount; ++i) {
        const int row = static_cast<int>(i);
        const Real n = sol.n(row);
        const Real p = sol.p(row);
        const Real ni = effectiveNi[i];
        const Real rate = recombination.totalRate(n, p, ni);
        const Real absRate = std::abs(rate);
        if (std::isfinite(absRate)) {
            diagnostics.maxAbsRate_m3_per_s = std::max(diagnostics.maxAbsRate_m3_per_s, absRate);
            absSum += absRate;
        }

        const Real ni2 = ni * ni;
        if (ni2 > 0.0) {
            const Real ratio = std::abs(n * p) / ni2;
            if (std::isfinite(ratio))
                diagnostics.maxCarrierProductRatio = std::max(diagnostics.maxCarrierProductRatio, ratio);
        }
    }

    diagnostics.meanAbsRate_m3_per_s = absSum / static_cast<Real>(nodeCount);
    return diagnostics;
}

SweepTransportDiagnostics computeSweepTransportDiagnostics(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const MobilityModelConfig& mobilityConfig,
    Real temperature_K,
    const DDSolution& sol)
{
    SweepTransportDiagnostics diagnostics;
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh, matdb, temperature_K);
    const std::unique_ptr<MobilityModel> mobility = makeMobilityModel(mobilityConfig);

    Real electronMobilitySum = 0.0;
    Real holeMobilitySum = 0.0;
    Real electronQfGradientSum = 0.0;
    Real holeQfGradientSum = 0.0;
    Real electronHighFieldDriveSum = 0.0;
    Real holeHighFieldDriveSum = 0.0;
    Real electronLimiterSum = 0.0;
    Real holeLimiterSum = 0.0;
    Index electronMobilityCount = 0;
    Index holeMobilityCount = 0;
    Index qfGradientCount = 0;
    Index highFieldDriveCount = 0;
    Index electronLimiterCount = 0;
    Index holeLimiterCount = 0;
    bool hasElectronMobility = false;
    bool hasHoleMobility = false;
    bool hasElectronLimiter = false;
    bool hasHoleLimiter = false;

    for (Index e = 0; e < mesh.numEdges(); ++e) {
        const Edge& edge = mesh.getEdge(e);
        const Node& n0 = mesh.getNode(edge.n0);
        const Node& n1 = mesh.getNode(edge.n1);
        const Real dx = n1.x - n0.x;
        const Real dy = n1.y - n0.y;
        const Real length = std::sqrt(dx * dx + dy * dy);
        if (length <= 0.0)
            continue;

        const int i0 = static_cast<int>(edge.n0);
        const int i1 = static_cast<int>(edge.n1);
        const Real electricField = std::abs(sol.psi(i1) - sol.psi(i0)) / length;
        diagnostics.maxElectricField_V_per_cm =
            std::max(diagnostics.maxElectricField_V_per_cm,
                     voltsPerMeterToVoltsPerCm(electricField));
        const Real electronMobilityField =
            mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient"
            ? std::abs(sol.phin(i1) - sol.phin(i0)) / length
            : electricField;
        const Real holeMobilityField =
            mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient"
            ? std::abs(sol.phip(i1) - sol.phip(i0)) / length
            : electricField;
        electronHighFieldDriveSum += voltsPerMeterToVoltsPerCm(electronMobilityField);
        holeHighFieldDriveSum += voltsPerMeterToVoltsPerCm(holeMobilityField);
        ++highFieldDriveCount;

        const Real electronLowFieldMobility = detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, e, CarrierType::Electron,
            0.0, &mobilityConfig, nullptr);
        const Real holeLowFieldMobility = detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, e, CarrierType::Hole,
            0.0, &mobilityConfig, nullptr);

        const Real electronMobility = detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, e, CarrierType::Electron,
            electronMobilityField, &mobilityConfig, &sol.psi);
        if (electronMobility > 0.0 && std::isfinite(electronMobility)) {
            electronMobilitySum += electronMobility;
            diagnostics.minElectronMobility_m2_V_s = hasElectronMobility
                ? std::min(diagnostics.minElectronMobility_m2_V_s, electronMobility)
                : electronMobility;
            hasElectronMobility = true;
            ++electronMobilityCount;
            if (electronLowFieldMobility > 0.0 && std::isfinite(electronLowFieldMobility)) {
                const Real limiter = electronMobility / electronLowFieldMobility;
                if (std::isfinite(limiter)) {
                    electronLimiterSum += limiter;
                    diagnostics.minElectronMobilityLimiter = hasElectronLimiter
                        ? std::min(diagnostics.minElectronMobilityLimiter, limiter)
                        : limiter;
                    hasElectronLimiter = true;
                    ++electronLimiterCount;
                }
            }
        }

        const Real holeMobility = detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, e, CarrierType::Hole,
            holeMobilityField, &mobilityConfig, &sol.psi);
        if (holeMobility > 0.0 && std::isfinite(holeMobility)) {
            holeMobilitySum += holeMobility;
            diagnostics.minHoleMobility_m2_V_s = hasHoleMobility
                ? std::min(diagnostics.minHoleMobility_m2_V_s, holeMobility)
                : holeMobility;
            hasHoleMobility = true;
            ++holeMobilityCount;
            if (holeLowFieldMobility > 0.0 && std::isfinite(holeLowFieldMobility)) {
                const Real limiter = holeMobility / holeLowFieldMobility;
                if (std::isfinite(limiter)) {
                    holeLimiterSum += limiter;
                    diagnostics.minHoleMobilityLimiter = hasHoleLimiter
                        ? std::min(diagnostics.minHoleMobilityLimiter, limiter)
                        : limiter;
                    hasHoleLimiter = true;
                    ++holeLimiterCount;
                }
            }
        }

        electronQfGradientSum += voltsPerMeterToVoltsPerCm(
            std::abs(sol.phin(i1) - sol.phin(i0)) / length);
        holeQfGradientSum += voltsPerMeterToVoltsPerCm(
            std::abs(sol.phip(i1) - sol.phip(i0)) / length);
        ++qfGradientCount;
    }

    if (electronMobilityCount > 0)
        diagnostics.meanElectronMobility_m2_V_s =
            electronMobilitySum / static_cast<Real>(electronMobilityCount);
    if (holeMobilityCount > 0)
        diagnostics.meanHoleMobility_m2_V_s =
            holeMobilitySum / static_cast<Real>(holeMobilityCount);
    if (qfGradientCount > 0) {
        diagnostics.meanElectronQfGradient_V_per_cm =
            electronQfGradientSum / static_cast<Real>(qfGradientCount);
        diagnostics.meanHoleQfGradient_V_per_cm =
            holeQfGradientSum / static_cast<Real>(qfGradientCount);
    }
    if (highFieldDriveCount > 0) {
        diagnostics.meanElectronHighFieldDrive_V_per_cm =
            electronHighFieldDriveSum / static_cast<Real>(highFieldDriveCount);
        diagnostics.meanHoleHighFieldDrive_V_per_cm =
            holeHighFieldDriveSum / static_cast<Real>(highFieldDriveCount);
    }
    if (electronLimiterCount > 0)
        diagnostics.meanElectronMobilityLimiter =
            electronLimiterSum / static_cast<Real>(electronLimiterCount);
    if (holeLimiterCount > 0)
        diagnostics.meanHoleMobilityLimiter =
            holeLimiterSum / static_cast<Real>(holeLimiterCount);
    return diagnostics;
}

std::vector<ContinuityBalanceDiagnosticRow> computeContinuityBalanceDiagnostics(
    const DeviceMesh& mesh,
    const MaterialDatabase& matdb,
    const DopingModel& doping,
    const MobilityModelConfig& mobilityConfig,
    Real temperature_K,
    const DDSolution& sol,
    const std::vector<Real>& effectiveNi,
    const RecombinationModelConfig& recombinationCfg,
    const std::vector<std::string>& contacts)
{
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh, matdb, temperature_K);
    const std::unique_ptr<MobilityModel> mobility = makeMobilityModel(mobilityConfig);
    const RecombinationModel recombination(recombinationCfg);
    const Real Vt = constants::kb * temperature_K / constants::q;

    std::unordered_set<std::string> requested(contacts.begin(), contacts.end());
    std::vector<std::vector<Index>> nodeEdges(mesh.numNodes());
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        nodeEdges[edge.n0].push_back(edgeId);
        nodeEdges[edge.n1].push_back(edgeId);
    }

    auto edgeFluxForCarrier = [&](Index edgeId, CarrierType carrier) {
        const Edge& edge = mesh.getEdge(edgeId);
        const Real h = edge.length;
        if (h <= 1.0e-30)
            return Real{0.0};
        const int i = static_cast<int>(edge.n0);
        const int j = static_cast<int>(edge.n1);
        const Real electricField = std::abs(sol.psi(j) - sol.psi(i)) / h;
        const Real drivingField = mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient"
            ? ((carrier == CarrierType::Electron)
                ? std::abs(sol.phin(j) - sol.phin(i)) / h
                : std::abs(sol.phip(j) - sol.phip(i)) / h)
            : electricField;
        const Real mu = detail::edgeMobility(
            edgeCells, mesh, doping, *mobility, cellMaterials, edgeId, carrier,
            drivingField, &mobilityConfig, &sol.psi);
        if (mu <= 0.0)
            return Real{0.0};
        const Real coef = mu * Vt * edge.couple / h;
        if (carrier == CarrierType::Electron) {
            return sgElectronContinuityFluxFromQuasiFermiVariableNi(
                effectiveNi[edge.n0],
                effectiveNi[edge.n1],
                sol.psi(i),
                sol.psi(j),
                sol.phin(i),
                sol.phin(j),
                Vt,
                coef);
        }
        return sgHoleContinuityFluxFromQuasiFermiVariableNi(
            effectiveNi[edge.n0],
            effectiveNi[edge.n1],
            sol.psi(i),
            sol.psi(j),
            sol.phip(i),
            sol.phip(j),
            Vt,
            coef);
    };

    auto nodeContribution = [&](Index edgeId, Index node, CarrierType carrier) {
        const Edge& edge = mesh.getEdge(edgeId);
        const Real flux = edgeFluxForCarrier(edgeId, carrier);
        if (edge.n0 == node)
            return flux;
        if (edge.n1 == node)
            return -flux;
        return Real{0.0};
    };

    auto recombinationTerm = [&](Index node) {
        const int row = static_cast<int>(node);
        const Real ni = effectiveNi[node];
        if (ni <= 0.0)
            return Real{0.0};
        return recombination.totalRate(sol.n(row), sol.p(row), ni) *
            mesh.getNode(node).volume;
    };

    std::vector<ContinuityBalanceDiagnosticRow> rows;
    for (const Contact& contact : mesh.contacts()) {
        if (!requested.empty() && !requested.contains(contact.name))
            continue;
        const std::unordered_set<Index> contactNodes(
            contact.node_ids.begin(), contact.node_ids.end());
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            const bool n0Contact = contactNodes.contains(edge.n0);
            const bool n1Contact = contactNodes.contains(edge.n1);
            if (n0Contact == n1Contact)
                continue;
            const Index contactNode = n0Contact ? edge.n0 : edge.n1;
            const Index interiorNode = n0Contact ? edge.n1 : edge.n0;

            for (CarrierType carrier : {CarrierType::Electron, CarrierType::Hole}) {
                Real contactFlux = 0.0;
                Real neighborFlux = 0.0;
                for (Index adjacentEdge : nodeEdges[interiorNode]) {
                    const Real contribution =
                        nodeContribution(adjacentEdge, interiorNode, carrier);
                    if (adjacentEdge == edgeId)
                        contactFlux += contribution;
                    else
                        neighborFlux += contribution;
                }
                const Real recombination = recombinationTerm(interiorNode);
                const int contactRow = static_cast<int>(contactNode);
                const int interiorRow = static_cast<int>(interiorNode);
                ContinuityBalanceDiagnosticRow row;
                row.contact = contact.name;
                row.carrier = carrier == CarrierType::Electron ? "electron" : "hole";
                row.contactNode = contactNode;
                row.interiorNode = interiorNode;
                row.contactEdgeId = edgeId;
                row.contactEdgeFlux = contactFlux;
                row.neighborEdgeFlux = neighborFlux;
                row.recombinationTerm = recombination;
                row.continuityResidual = contactFlux + neighborFlux + recombination;
                row.interiorVolume_m2 = mesh.getNode(interiorNode).volume;
                row.qfContact_V = carrier == CarrierType::Electron
                    ? sol.phin(contactRow)
                    : sol.phip(contactRow);
                row.qfInterior_V = carrier == CarrierType::Electron
                    ? sol.phin(interiorRow)
                    : sol.phip(interiorRow);
                row.carrierDensityInterior_m3 = carrier == CarrierType::Electron
                    ? sol.n(interiorRow)
                    : sol.p(interiorRow);
                rows.push_back(std::move(row));
            }
        }
    }
    return rows;
}

const Contact& requireContactByName(const DeviceMesh& mesh, const std::string& contactName)
{
    for (const Contact& contact : mesh.contacts()) {
        if (contact.name == contactName)
            return contact;
    }
    throw std::runtime_error(
        "DCSweep: contact_current_qf_floor references unknown contact '" +
        contactName + "'.");
}

ContactCurrentEdgeOverrides buildContactCurrentQfFloorOverrides(
    const DeviceMesh& mesh,
    const DDSolution& initial,
    const std::vector<std::string>& contacts)
{
    ContactCurrentEdgeOverrides overrides;
    if (initial.phip.size() != static_cast<int>(mesh.numNodes()))
        return overrides;

    for (const std::string& contactName : contacts) {
        const Contact& contact = requireContactByName(mesh, contactName);
        const std::unordered_set<Index> contactNodes(
            contact.node_ids.begin(), contact.node_ids.end());
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            const bool n0Contact = contactNodes.contains(edge.n0);
            const bool n1Contact = contactNodes.contains(edge.n1);
            if (n0Contact == n1Contact)
                continue;

            const Real drop =
                initial.phip(static_cast<int>(edge.n1)) -
                initial.phip(static_cast<int>(edge.n0));
            if (std::isfinite(drop) && drop != 0.0)
                overrides.holeQuasiFermiDropByEdge[edgeId] = drop;
        }
    }
    return overrides;
}

std::filesystem::path resolveConfigPath(const std::filesystem::path& cfgDir,
                                        const std::string& path)
{
    std::filesystem::path fp(path);
    if (fp.is_relative())
        fp = cfgDir / fp;
    return fp;
}

long long parseNodeDopingNodeId(const std::string& nodeIdText)
{
    std::size_t consumed = 0;
    long long parsedNodeId = 0;
    try {
        parsedNodeId = std::stoll(nodeIdText, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "DCSweep: node_doping_file has invalid node id '" + nodeIdText + "'");
    }
    if (consumed != nodeIdText.size()) {
        throw std::runtime_error(
            "DCSweep: node_doping_file has invalid node id '" + nodeIdText + "'");
    }
    return parsedNodeId;
}

Real parseNodeDopingConcentration(const std::string& value,
                                  const std::string& column,
                                  long long nodeId,
                                  UnitScalingConfig scaling)
{
    std::size_t consumed = 0;
    Real parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "DCSweep: node_doping_file has invalid " + column + " '" + value +
            "' for node id " + std::to_string(nodeId));
    }
    if (consumed != value.size()) {
        throw std::runtime_error(
            "DCSweep: node_doping_file has invalid " + column + " '" + value +
            "' for node id " + std::to_string(nodeId));
    }
    if (!std::isfinite(parsed)) {
        throw std::runtime_error(
            "DCSweep: node_doping_file has non-finite " + column + " '" + value +
            "' for node id " + std::to_string(nodeId));
    }
    return scaling.concentrationToSI(parsed);
}

TerminalChargeConfig terminalChargeConfigFromJson(const nlohmann::json& chargeCfg,
                                                  const DCSweepConfig& sweep,
                                                  std::size_t index,
                                                  UnitScalingConfig scaling)
{
    TerminalChargeConfig config;
    config.name = chargeCfg.value("name", std::string{});
    config.contact = chargeCfg.value("contact", sweep.chargeContact.empty() ? sweep.contact : sweep.chargeContact);
    config.regions = chargeCfg.value("regions", sweep.chargeRegions);
    config.contactRadius = chargeCfg.contains("contact_radius")
        ? scaling.lengthToSI(chargeCfg.at("contact_radius").get<Real>())
        : sweep.chargeContactRadius;
    config.includeMobileCharge = chargeCfg.value("include_mobile_charge", config.includeMobileCharge);
    config.includeIonizedDopants = chargeCfg.value("include_ionized_dopants", config.includeIonizedDopants);
    config.perMeter = chargeCfg.value("per_meter", sweep.chargePerMeter);
    config.depth_m = chargeCfg.contains("depth_m")
        ? scaling.lengthToSI(chargeCfg.at("depth_m").get<Real>())
        : sweep.chargeDepth_m;
    if (config.name.empty())
        config.name = terminalChargeName(config, index);
    if (!config.perMeter && config.depth_m <= 0.0)
        throw std::invalid_argument("DCSweep: sweep terminal charge depth_m must be positive.");
    return config;
}

DDSolutionValidationOptions validationOptionsFromJson(const nlohmann::json& cfg)
{
    DDSolutionValidationOptions options;
    const auto validation = cfg.value("validation", nlohmann::json::object());
    options.carrierFloor = validation.value("carrier_floor", options.carrierFloor);
    options.enforceMinimumCarrierDensity =
        validation.value("enforce_minimum_carrier_density", options.enforceMinimumCarrierDensity);
    options.minimumCarrierDensity =
        validation.value("minimum_carrier_density", options.minimumCarrierDensity);
    options.checkContactQuasiFermiBias =
        validation.value("check_contact_quasi_fermi_bias", options.checkContactQuasiFermiBias);
    options.contactPotentialAbsTolerance =
        validation.value("contact_potential_abs_tolerance", options.contactPotentialAbsTolerance);
    options.contactPotentialRelTolerance =
        validation.value("contact_potential_rel_tolerance", options.contactPotentialRelTolerance);
    if (options.carrierFloor < 0.0)
        throw std::invalid_argument("DCSweep: validation.carrier_floor must be non-negative.");
    if (options.enforceMinimumCarrierDensity && options.minimumCarrierDensity < 0.0)
        throw std::invalid_argument(
            "DCSweep: validation.minimum_carrier_density must be non-negative.");
    if (options.contactPotentialAbsTolerance < 0.0 ||
        options.contactPotentialRelTolerance < 0.0) {
        throw std::invalid_argument(
            "DCSweep: validation contact tolerances must be non-negative.");
    }
    return options;
}


std::string normalizedSolverMethod(const nlohmann::json& cfg)
{
    std::string method = "gummel";
    if (cfg.contains("solver")) {
        const auto& solver = cfg.at("solver");
        if (solver.contains("method"))
            method = solver.at("method").get<std::string>();
        else if (solver.contains("type"))
            method = solver.at("type").get<std::string>();
    }

    std::transform(method.begin(), method.end(), method.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return method;
}

SolverMethod solverMethodFromJson(const nlohmann::json& cfg)
{
    const std::string method = normalizedSolverMethod(cfg);
    if (method == "gummel")
        return SolverMethod::Gummel;
    if (method == "newton")
        return SolverMethod::Newton;
    if (method == "gummel_newton" || method == "hybrid")
        return SolverMethod::GummelNewton;
    throw std::invalid_argument(
        "DCSweep: solver.method/type must be 'gummel', 'newton', or 'gummel_newton'.");
}

HybridHandoffConfig hybridHandoffConfigFromJson(const nlohmann::json& solverJson)
{
    HybridHandoffConfig hybrid;
    if (!solverJson.contains("handoff"))
        return hybrid;

    const auto& handoff = solverJson.at("handoff");
    if (!handoff.is_object())
        throw std::invalid_argument("DCSweep: solver.handoff must be an object.");

    const std::string fallback = handoff.value("fallback", std::string("none"));
    if (fallback == "none") {
        hybrid.fallbackToGummelOnNewtonFailure = false;
    } else if (fallback == "gummel_on_newton_failure") {
        hybrid.fallbackToGummelOnNewtonFailure = true;
    } else {
        throw std::invalid_argument(
            "DCSweep: solver.handoff.fallback must be 'none' or "
            "'gummel_on_newton_failure'.");
    }

    hybrid.requireGummelConvergence =
        handoff.value("require_gummel_convergence", hybrid.requireGummelConvergence);
    if (handoff.contains("gummel_max_iter")) {
        hybrid.gummelMaxIter = handoff.at("gummel_max_iter").get<int>();
        if (hybrid.gummelMaxIter < 0)
            throw std::invalid_argument(
                "DCSweep: solver.handoff.gummel_max_iter must be non-negative.");
    }
    if (handoff.contains("newton_max_iter")) {
        hybrid.newtonMaxIter = handoff.at("newton_max_iter").get<int>();
        if (hybrid.newtonMaxIter < 0)
            throw std::invalid_argument(
                "DCSweep: solver.handoff.newton_max_iter must be non-negative.");
    }
    return hybrid;
}

bool isAllowedSweepPredictorMode(const std::string& mode)
{
    return mode == "none" || mode == "constant" || mode == "linear" || mode == "secant";
}

bool isAllowedSweepPredictorField(const std::string& field)
{
    return field == "psi" || field == "phin" || field == "phip";
}

void parseSweepContinuationConfig(const nlohmann::json& sweepJson,
                                  DCSweepConfig&        sweep)
{
    if (!sweepJson.contains("continuation"))
        return;
    const auto& continuation = sweepJson.at("continuation");
    if (!continuation.is_object())
        throw std::invalid_argument("DCSweep: sweep.continuation must be an object.");

    if (continuation.contains("predictor")) {
        const auto& predictor = continuation.at("predictor");
        if (!predictor.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.predictor must be an object.");
        sweep.continuation.predictor.mode =
            predictor.value("mode", sweep.continuation.predictor.mode);
        if (!isAllowedSweepPredictorMode(sweep.continuation.predictor.mode)) {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.predictor.mode must be one of "
                "'none', 'constant', 'linear', or 'secant'.");
        }
        if (predictor.contains("fields")) {
            const auto& fields = predictor.at("fields");
            if (!fields.is_array())
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.predictor.fields must be an array.");
            sweep.continuation.predictor.fields.clear();
            for (const auto& entry : fields) {
                const std::string field = entry.get<std::string>();
                if (!isAllowedSweepPredictorField(field)) {
                    throw std::invalid_argument(
                        "DCSweep: sweep.continuation.predictor.fields entries must be "
                        "'psi', 'phin', or 'phip'.");
                }
                sweep.continuation.predictor.fields.push_back(field);
            }
        }
        sweep.continuation.predictor.maxExtrapolationRatio = predictor.value(
            "max_extrapolation_ratio",
            sweep.continuation.predictor.maxExtrapolationRatio);
        if (!std::isfinite(sweep.continuation.predictor.maxExtrapolationRatio) ||
            sweep.continuation.predictor.maxExtrapolationRatio < 1.0) {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.predictor.max_extrapolation_ratio "
                "must be finite and at least 1.");
        }
        if (sweep.continuation.predictor.mode != "none" &&
            sweep.continuation.predictor.fields.empty()) {
            sweep.continuation.predictor.fields = {"psi", "phin", "phip"};
        }
    }

    if (continuation.contains("branch_acceptance")) {
        const auto& branchAcceptance = continuation.at("branch_acceptance");
        if (!branchAcceptance.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.branch_acceptance must be an object.");
        sweep.continuation.branchAcceptance.terminalCurrentConsistency =
            branchAcceptance.value(
                "terminal_current_consistency",
                sweep.continuation.branchAcceptance.terminalCurrentConsistency);
        sweep.continuation.branchAcceptance.minTerminalCurrentRatio =
            branchAcceptance.value(
                "min_terminal_current_ratio",
                sweep.continuation.branchAcceptance.minTerminalCurrentRatio);
        if (!std::isfinite(sweep.continuation.branchAcceptance.minTerminalCurrentRatio) ||
            sweep.continuation.branchAcceptance.minTerminalCurrentRatio < 0.0) {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.branch_acceptance.min_terminal_current_ratio "
                "must be finite and non-negative.");
        }
        sweep.continuation.branchAcceptance.psiPhinJump =
            branchAcceptance.value(
                "psi_phin_jump",
                sweep.continuation.branchAcceptance.psiPhinJump);
        if (branchAcceptance.contains("max_psi_phin_jump_V")) {
            sweep.continuation.branchAcceptance.maxPsiPhinJump_V =
                branchAcceptance.at("max_psi_phin_jump_V").get<Real>();
        }
        if (sweep.continuation.branchAcceptance.psiPhinJump &&
            (!std::isfinite(sweep.continuation.branchAcceptance.maxPsiPhinJump_V) ||
             sweep.continuation.branchAcceptance.maxPsiPhinJump_V < 0.0)) {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.branch_acceptance.max_psi_phin_jump_V "
                "must be finite and non-negative when psi_phin_jump is enabled.");
        }
        sweep.continuation.branchAcceptance.carrierDensityJump =
            branchAcceptance.value(
                "carrier_density_jump",
                sweep.continuation.branchAcceptance.carrierDensityJump);
        if (branchAcceptance.contains("max_electron_density_jump_dex")) {
            sweep.continuation.branchAcceptance.maxElectronDensityJumpDex =
                branchAcceptance.at("max_electron_density_jump_dex").get<Real>();
        }
        if (branchAcceptance.contains("max_electron_density_jump_p95_abs_dex")) {
            sweep.continuation.branchAcceptance.maxElectronDensityJumpP95AbsDex =
                branchAcceptance.at("max_electron_density_jump_p95_abs_dex").get<Real>();
        }
        if (sweep.continuation.branchAcceptance.carrierDensityJump &&
            (!std::isfinite(sweep.continuation.branchAcceptance.maxElectronDensityJumpDex) ||
             sweep.continuation.branchAcceptance.maxElectronDensityJumpDex < 0.0)) {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.branch_acceptance.max_electron_density_jump_dex "
                "must be finite and non-negative when carrier_density_jump is enabled.");
        }
        if (sweep.continuation.branchAcceptance.carrierDensityJump &&
            branchAcceptance.contains("max_electron_density_jump_p95_abs_dex") &&
            (!std::isfinite(sweep.continuation.branchAcceptance.maxElectronDensityJumpP95AbsDex) ||
             sweep.continuation.branchAcceptance.maxElectronDensityJumpP95AbsDex < 0.0)) {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.branch_acceptance."
                "max_electron_density_jump_p95_abs_dex "
                "must be finite and non-negative when carrier_density_jump is enabled.");
        }
    }

    if (continuation.contains("arclength")) {
        const auto& arclength = continuation.at("arclength");
        if (!arclength.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.arclength must be an object.");
        SweepArclengthConfig& cfg = sweep.continuation.arclength;
        cfg.enabled = arclength.value("enabled", cfg.enabled);
        cfg.predictor = arclength.value("predictor", cfg.predictor);
        if (cfg.predictor != "tangent") {
            throw std::invalid_argument(
                "DCSweep: sweep.continuation.arclength.predictor must be 'tangent'.");
        }
        cfg.core.enabled = cfg.enabled;
        cfg.core.initialStep = arclength.value("initial_step", cfg.core.initialStep);
        cfg.core.minStep = arclength.value("min_step", cfg.core.minStep);
        cfg.core.maxStep = arclength.value("max_step", cfg.core.maxStep);
        cfg.core.growthFactor = arclength.value("growth_factor", cfg.core.growthFactor);
        cfg.core.shrinkFactor = arclength.value("shrink_factor", cfg.core.shrinkFactor);
        cfg.core.maxCorrectorIterations =
            arclength.value("max_corrector_iterations", cfg.core.maxCorrectorIterations);
        cfg.core.correctorTolerance =
            arclength.value("corrector_tolerance", cfg.core.correctorTolerance);
        cfg.core.maxStepRetries = arclength.value("max_step_retries", cfg.core.maxStepRetries);
        cfg.core.parameterScale = arclength.value("parameter_scale", cfg.core.parameterScale);
        cfg.core.stateWeight = arclength.value("state_weight", cfg.core.stateWeight);
        cfg.core.dampingFactor = arclength.value("damping_factor", cfg.core.dampingFactor);
        cfg.core.maxLineSearchSteps =
            arclength.value("max_line_search_steps", cfg.core.maxLineSearchSteps);
        cfg.biasFiniteDifferenceStep_V =
            arclength.value("bias_finite_difference_step_V", cfg.biasFiniteDifferenceStep_V);

        if (cfg.enabled) {
            auto requirePositive = [&](Real value, const char* key) {
                if (!std::isfinite(value) || value <= 0.0) {
                    throw std::invalid_argument(
                        std::string("DCSweep: sweep.continuation.arclength.") + key +
                        " must be finite and positive when arclength continuation is enabled.");
                }
            };
            requirePositive(cfg.core.initialStep, "initial_step");
            requirePositive(cfg.core.minStep, "min_step");
            requirePositive(cfg.core.maxStep, "max_step");
            requirePositive(cfg.core.parameterScale, "parameter_scale");
            requirePositive(cfg.biasFiniteDifferenceStep_V, "bias_finite_difference_step_V");
            if (!std::isfinite(cfg.core.stateWeight) || cfg.core.stateWeight < 0.0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.state_weight must be finite "
                    "and non-negative.");
            }
            if (!std::isfinite(cfg.core.dampingFactor) ||
                cfg.core.dampingFactor <= 0.0 || cfg.core.dampingFactor > 1.0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.damping_factor must be in (0, 1].");
            }
            if (cfg.core.maxLineSearchSteps < 0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.max_line_search_steps "
                    "must be non-negative.");
            }
            if (cfg.core.minStep > cfg.core.maxStep) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.min_step must not exceed max_step.");
            }
            if (cfg.core.initialStep < cfg.core.minStep ||
                cfg.core.initialStep > cfg.core.maxStep) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.initial_step must lie within "
                    "[min_step, max_step].");
            }
            if (!std::isfinite(cfg.core.growthFactor) || cfg.core.growthFactor < 1.0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.growth_factor must be finite "
                    "and at least 1.");
            }
            if (!std::isfinite(cfg.core.shrinkFactor) ||
                cfg.core.shrinkFactor <= 0.0 || cfg.core.shrinkFactor >= 1.0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.shrink_factor must be in (0, 1).");
            }
            if (cfg.core.maxCorrectorIterations <= 0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.max_corrector_iterations "
                    "must be positive.");
            }
            if (!std::isfinite(cfg.core.correctorTolerance) ||
                cfg.core.correctorTolerance <= 0.0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.corrector_tolerance "
                    "must be finite and positive.");
            }
            if (cfg.core.maxStepRetries < 0) {
                throw std::invalid_argument(
                    "DCSweep: sweep.continuation.arclength.max_step_retries "
                    "must be non-negative.");
            }
        }
    }
}

DCSweepConfig dcSweepConfigFromJson(const nlohmann::json& cfg,
                                    const std::filesystem::path& cfgDir,
                                    UnitScalingConfig scaling)
{
    const auto& j = cfg.at("sweep");
    DCSweepConfig sweep;
    sweep.scaling = scaling;
    sweep.mode = curveSweepModeFromString(j.value("mode", std::string("iv")));
    sweep.contact = j.at("contact").get<std::string>();
    sweep.start = j.at("start").get<Real>();
    sweep.stop = j.at("stop").get<Real>();
    sweep.step = j.at("step").get<Real>();
    if (j.contains("bias_points")) {
        const auto& biasPoints = j.at("bias_points");
        if (!biasPoints.is_array())
            throw std::invalid_argument("DCSweep: sweep.bias_points must be an array.");
        if (biasPoints.empty())
            throw std::invalid_argument("DCSweep: sweep.bias_points must not be empty.");
        for (const auto& entry : biasPoints) {
            const Real bias = entry.get<Real>();
            if (!std::isfinite(bias))
                throw std::invalid_argument("DCSweep: sweep.bias_points entries must be finite.");
            sweep.biasPoints.push_back(bias);
        }
    }
    const Real nominalStep = std::abs(sweep.step);
    sweep.shrinkFactor = j.value("shrink_factor", sweep.shrinkFactor);
    sweep.growthFactor = j.value("growth_factor", sweep.growthFactor);
    sweep.maxRetries = j.value("max_retries", sweep.maxRetries);
    sweep.minStep = j.value("min_step", nominalStep * std::pow(sweep.shrinkFactor, sweep.maxRetries));
    sweep.maxStep = j.value("max_step", nominalStep);
    sweep.stopOnFailure = j.value("stop_on_failure", sweep.stopOnFailure);
    sweep.currentContact = j.value("current_contact", sweep.contact);
    sweep.writeVtk = j.value("write_vtk", cfg.value("write_vtk", false));
    sweep.csvFile = j.value("csv_file", cfg.value("output_csv", sweep.csvFile));
    sweep.vtkPrefix = j.value("vtk_prefix", cfg.value("output_vtk_prefix", std::string("dc_sweep")));
    sweep.initialStateFile = j.value("initial_state_file", std::string{});
    sweep.writeStateFile = j.value("write_state_file", std::string{});
    sweep.writeStateEveryPointPrefix =
        j.value("write_state_every_point_prefix", std::string{});
    parseSweepContinuationConfig(j, sweep);

    const auto chargeCfg = j.value("terminal_charge", nlohmann::json::object());
    sweep.chargeContact = chargeCfg.value("contact", j.value("charge_contact", sweep.contact));
    sweep.chargeRegions = chargeCfg.value("regions", j.value("charge_regions", std::vector<std::string>{}));
    sweep.chargeContactRadius = scaling.lengthToSI(
        chargeCfg.value("contact_radius", j.value("charge_contact_radius", 0.0)));
    sweep.chargePerMeter = chargeCfg.value("per_meter", j.value("charge_per_meter", true));
    sweep.chargeDepth_m = scaling.lengthToSI(
        chargeCfg.value("depth_m", j.value("charge_depth_m", 1.0)));

    sweep.storedChargeEnabled = j.contains("stored_charge");
    const auto storedCfg = j.value("stored_charge", nlohmann::json::object());
    sweep.storedCharge.regions = storedCfg.value("regions", std::vector<std::string>{});
    sweep.storedCharge.perMeter = storedCfg.value("per_meter", true);
    sweep.storedCharge.depth_m = scaling.lengthToSI(storedCfg.value("depth_m", 1.0));

    const auto diagnosticsCfg = j.value("diagnostics", nlohmann::json::object());
    if (diagnosticsCfg.contains("terminal_balance")) {
        const auto& terminalBalanceCfg = diagnosticsCfg.at("terminal_balance");
        if (!terminalBalanceCfg.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.diagnostics.terminal_balance must be an object.");
        sweep.diagnostics.terminalBalance.enabled =
            terminalBalanceCfg.value("enabled", sweep.diagnostics.terminalBalance.enabled);
        sweep.diagnostics.terminalBalance.contacts =
            terminalBalanceCfg.value("contacts", std::vector<std::string>{});
        sweep.diagnostics.terminalBalance.csvFile =
            terminalBalanceCfg.value("csv_file", std::string{});
    }
    if (diagnosticsCfg.contains("contact_edge")) {
        const auto& contactEdgeCfg = diagnosticsCfg.at("contact_edge");
        if (!contactEdgeCfg.is_object())
            throw std::invalid_argument("DCSweep: sweep.diagnostics.contact_edge must be an object.");
        sweep.diagnostics.contactEdge.enabled =
            contactEdgeCfg.value("enabled", sweep.diagnostics.contactEdge.enabled);
        sweep.diagnostics.contactEdge.contacts =
            contactEdgeCfg.value("contacts", std::vector<std::string>{});
        sweep.diagnostics.contactEdge.csvFile =
            contactEdgeCfg.value("csv_file", std::string{});
    }
    if (diagnosticsCfg.contains("transport")) {
        const auto& transportCfg = diagnosticsCfg.at("transport");
        if (!transportCfg.is_object())
            throw std::invalid_argument("DCSweep: sweep.diagnostics.transport must be an object.");
        sweep.diagnostics.transport.enabled =
            transportCfg.value("enabled", sweep.diagnostics.transport.enabled);
    }
    if (diagnosticsCfg.contains("continuity_balance")) {
        const auto& continuityCfg = diagnosticsCfg.at("continuity_balance");
        if (!continuityCfg.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.diagnostics.continuity_balance must be an object.");
        sweep.diagnostics.continuityBalance.enabled =
            continuityCfg.value("enabled", sweep.diagnostics.continuityBalance.enabled);
        sweep.diagnostics.continuityBalance.contacts =
            continuityCfg.value("contacts", std::vector<std::string>{});
        sweep.diagnostics.continuityBalance.csvFile =
            continuityCfg.value("csv_file", std::string{});
    }
    if (diagnosticsCfg.contains("sg_avalanche_edges")) {
        const auto& sgAvalancheCfg = diagnosticsCfg.at("sg_avalanche_edges");
        if (!sgAvalancheCfg.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.diagnostics.sg_avalanche_edges must be an object.");
        sweep.diagnostics.sgAvalancheEdges.enabled =
            sgAvalancheCfg.value("enabled", sweep.diagnostics.sgAvalancheEdges.enabled);
        sweep.diagnostics.sgAvalancheEdges.csvFile =
            sgAvalancheCfg.value("csv_file", std::string{});
    }
    if (diagnosticsCfg.contains("newton_history")) {
        const auto& newtonHistoryCfg = diagnosticsCfg.at("newton_history");
        if (!newtonHistoryCfg.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.diagnostics.newton_history must be an object.");
        sweep.diagnostics.newtonHistory.enabled =
            newtonHistoryCfg.value("enabled", sweep.diagnostics.newtonHistory.enabled);
        sweep.diagnostics.newtonHistory.csvFile =
            newtonHistoryCfg.value("csv_file", std::string{});
    }
    if (diagnosticsCfg.contains("contact_current_qf_floor")) {
        const auto& qfFloorCfg = diagnosticsCfg.at("contact_current_qf_floor");
        if (!qfFloorCfg.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.diagnostics.contact_current_qf_floor must be an object.");
        sweep.diagnostics.contactCurrentQfFloor.enabled =
            qfFloorCfg.value("enabled", sweep.diagnostics.contactCurrentQfFloor.enabled);
        sweep.diagnostics.contactCurrentQfFloor.contacts =
            qfFloorCfg.value("contacts", std::vector<std::string>{});
    }
    const auto contactCurrentReportingCfg =
        j.value("contact_current_reporting", nlohmann::json::object());
    if (!contactCurrentReportingCfg.is_object())
        throw std::invalid_argument(
            "DCSweep: sweep.contact_current_reporting must be an object.");
    if (contactCurrentReportingCfg.contains("endpoint_qf_floor")) {
        const auto& qfFloorCfg = contactCurrentReportingCfg.at("endpoint_qf_floor");
        if (!qfFloorCfg.is_object())
            throw std::invalid_argument(
                "DCSweep: sweep.contact_current_reporting.endpoint_qf_floor "
                "must be an object.");
        sweep.diagnostics.contactCurrentQfFloor.enabled =
            qfFloorCfg.value("enabled", sweep.diagnostics.contactCurrentQfFloor.enabled);
        sweep.diagnostics.contactCurrentQfFloor.contacts =
            qfFloorCfg.value("contacts", sweep.diagnostics.contactCurrentQfFloor.contacts);
    }

    if (j.contains("terminal_charges")) {
        const auto& charges = j.at("terminal_charges");
        if (!charges.is_array())
            throw std::invalid_argument("DCSweep: sweep.terminal_charges must be an array.");
        if (charges.empty())
            throw std::invalid_argument("DCSweep: sweep.terminal_charges must not be empty.");
        std::unordered_set<std::string> names;
        for (std::size_t i = 0; i < charges.size(); ++i) {
            TerminalChargeConfig config =
                terminalChargeConfigFromJson(charges.at(i), sweep, i, scaling);
            const std::string name = terminalChargeName(config, i);
            if (!names.insert(name).second)
                throw std::invalid_argument("DCSweep: duplicate terminal_charges name '" + name + "'.");
            config.name = name;
            sweep.terminalCharges.push_back(std::move(config));
        }
    } else {
        sweep.terminalCharges.push_back(
            terminalChargeConfigFromJson(chargeCfg, sweep, 0, scaling));
    }

    const auto bvCfg = j.value("breakdown", nlohmann::json::object());
    sweep.breakdown.maxElectricField_V_per_m = scaling.electricFieldToSI(
        bvCfg.value("max_electric_field_V_per_m",
                    j.value("breakdown_max_electric_field_V_per_m", 0.0)));
    sweep.breakdown.currentJumpRatio = bvCfg.value("current_jump_ratio", j.value("breakdown_current_jump_ratio", 0.0));
    sweep.breakdown.nonConvergenceBreakdown = bvCfg.value("non_convergence", j.value("breakdown_on_non_convergence", true));

    auto resolve = [&](std::string path) {
        std::filesystem::path fp(path);
        if (fp.is_relative())
            fp = cfgDir / fp;
        return fp.string();
    };
    sweep.csvFile = resolve(sweep.csvFile);
    sweep.vtkPrefix = resolve(sweep.vtkPrefix);
    if (!sweep.initialStateFile.empty())
        sweep.initialStateFile = resolve(sweep.initialStateFile);
    if (!sweep.writeStateFile.empty())
        sweep.writeStateFile = resolve(sweep.writeStateFile);
    if (!sweep.writeStateEveryPointPrefix.empty())
        sweep.writeStateEveryPointPrefix = resolve(sweep.writeStateEveryPointPrefix);
    if (sweep.diagnostics.terminalBalance.enabled) {
        if (sweep.diagnostics.terminalBalance.csvFile.empty()) {
            const std::filesystem::path csvPath(sweep.csvFile);
            sweep.diagnostics.terminalBalance.csvFile =
                (csvPath.parent_path() / (csvPath.stem().string() + "_terminal_balance.csv")).string();
        } else {
            sweep.diagnostics.terminalBalance.csvFile = resolve(sweep.diagnostics.terminalBalance.csvFile);
        }
    }
    if (sweep.diagnostics.contactEdge.enabled) {
        if (sweep.diagnostics.contactEdge.csvFile.empty()) {
            const std::filesystem::path csvPath(sweep.csvFile);
            sweep.diagnostics.contactEdge.csvFile =
                (csvPath.parent_path() / (csvPath.stem().string() + "_contact_edges.csv")).string();
        } else {
            sweep.diagnostics.contactEdge.csvFile = resolve(sweep.diagnostics.contactEdge.csvFile);
        }
    }
    if (sweep.diagnostics.continuityBalance.enabled) {
        if (sweep.diagnostics.continuityBalance.csvFile.empty()) {
            const std::filesystem::path csvPath(sweep.csvFile);
            sweep.diagnostics.continuityBalance.csvFile =
                (csvPath.parent_path() / (csvPath.stem().string() + "_continuity_balance.csv")).string();
        } else {
            sweep.diagnostics.continuityBalance.csvFile =
                resolve(sweep.diagnostics.continuityBalance.csvFile);
        }
    }
    if (sweep.diagnostics.sgAvalancheEdges.enabled) {
        if (sweep.diagnostics.sgAvalancheEdges.csvFile.empty()) {
            const std::filesystem::path csvPath(sweep.csvFile);
            sweep.diagnostics.sgAvalancheEdges.csvFile =
                (csvPath.parent_path() / (csvPath.stem().string() + "_sg_avalanche_edges.csv")).string();
        } else {
            sweep.diagnostics.sgAvalancheEdges.csvFile =
                resolve(sweep.diagnostics.sgAvalancheEdges.csvFile);
        }
    }
    if (sweep.diagnostics.newtonHistory.enabled) {
        if (sweep.diagnostics.newtonHistory.csvFile.empty()) {
            const std::filesystem::path csvPath(sweep.csvFile);
            sweep.diagnostics.newtonHistory.csvFile =
                (csvPath.parent_path() / (csvPath.stem().string() + "_newton_history.csv")).string();
        } else {
            sweep.diagnostics.newtonHistory.csvFile =
                resolve(sweep.diagnostics.newtonHistory.csvFile);
        }
    }

    if (sweep.step == 0.0)
        throw std::invalid_argument("DCSweep: sweep.step must be non-zero.");
    if ((sweep.stop - sweep.start) * sweep.step < 0.0)
        throw std::invalid_argument("DCSweep: sweep.step sign must move start toward stop.");
    if (sweep.minStep <= 0.0)
        throw std::invalid_argument("DCSweep: sweep.min_step must be positive.");
    if (sweep.maxStep <= 0.0)
        throw std::invalid_argument("DCSweep: sweep.max_step must be positive.");
    if (sweep.minStep > sweep.maxStep)
        throw std::invalid_argument("DCSweep: sweep.min_step must not exceed sweep.max_step.");
    if (sweep.growthFactor < 1.0)
        throw std::invalid_argument("DCSweep: sweep.growth_factor must be at least 1.");
    if (sweep.shrinkFactor <= 0.0 || sweep.shrinkFactor >= 1.0)
        throw std::invalid_argument("DCSweep: sweep.shrink_factor must be greater than 0 and less than 1.");
    if (sweep.maxRetries < 0)
        throw std::invalid_argument("DCSweep: sweep.max_retries must be non-negative.");
    if (sweep.mode == CurveSweepMode::BVReverse && sweep.start * sweep.stop < 0.0)
        throw std::invalid_argument("DCSweep: bv_reverse sweeps must stay on one reverse-bias polarity side.");
    if (!sweep.chargePerMeter && sweep.chargeDepth_m <= 0.0)
        throw std::invalid_argument("DCSweep: sweep terminal charge depth_m must be positive.");
    if (sweep.storedChargeEnabled && !sweep.storedCharge.perMeter && sweep.storedCharge.depth_m <= 0.0)
        throw std::invalid_argument("DCSweep: sweep stored_charge depth_m must be positive.");
    return sweep;
}

DopingModel dopingFromJson(const DeviceMesh& mesh,
                           const nlohmann::json& cfg,
                           const std::filesystem::path& cfgDir,
                           UnitScalingConfig scaling)
{
    if (cfg.contains("node_doping_file")) {
        const std::filesystem::path path =
            resolveConfigPath(cfgDir, cfg.at("node_doping_file").get<std::string>());
        std::ifstream input(path);
        if (!input.is_open())
            throw std::runtime_error("DCSweep: cannot open node_doping_file: " + path.string());

        std::string line;
        if (!std::getline(input, line))
            throw std::runtime_error("DCSweep: node_doping_file is empty: " + path.string());

        const std::vector<std::string> header = splitCsvLine(
            line,
            "DCSweep: node_doping_file does not support quoted fields.");
        std::unordered_map<std::string, std::size_t> columns;
        for (std::size_t i = 0; i < header.size(); ++i)
            columns[header.at(i)] = i;
        for (const std::string& required : {"node_id", "donors_cm3", "acceptors_cm3"}) {
            if (!columns.contains(required))
                throw std::runtime_error(
                    "DCSweep: node_doping_file missing required column '" + required + "'.");
        }

        DopingModel model(mesh.numNodes());
        std::vector<bool> seen(static_cast<std::size_t>(mesh.numNodes()), false);
        while (std::getline(input, line)) {
            if (trimCsvToken(line).empty())
                continue;
            const std::vector<std::string> row = splitCsvLine(
                line,
                "DCSweep: node_doping_file does not support quoted fields.");
            const auto requireColumn = [&](const std::string& name) -> const std::string& {
                const std::size_t index = columns.at(name);
                if (index >= row.size())
                    throw std::runtime_error(
                        "DCSweep: node_doping_file row missing column '" + name + "'.");
                return row.at(index);
            };

            const std::string& nodeIdText = requireColumn("node_id");
            const long long parsedNodeId = parseNodeDopingNodeId(nodeIdText);
            if (parsedNodeId < 0 ||
                static_cast<std::size_t>(parsedNodeId) >= static_cast<std::size_t>(mesh.numNodes())) {
                throw std::runtime_error(
                    "DCSweep: node_doping_file references missing node id " +
                    std::to_string(parsedNodeId));
            }
            const Index nodeId = static_cast<Index>(parsedNodeId);
            const std::size_t nodeIndex = static_cast<std::size_t>(parsedNodeId);
            if (seen.at(nodeIndex)) {
                throw std::runtime_error(
                    "DCSweep: node_doping_file has duplicate row for node id " +
                    std::to_string(parsedNodeId));
            }
            seen.at(nodeIndex) = true;
            const Real donors = parseNodeDopingConcentration(
                requireColumn("donors_cm3"), "donors_cm3", parsedNodeId, scaling);
            const Real acceptors = parseNodeDopingConcentration(
                requireColumn("acceptors_cm3"), "acceptors_cm3", parsedNodeId, scaling);
            model.setNodeDoping(nodeId, donors, acceptors);
        }
        for (std::size_t nodeId = 0; nodeId < seen.size(); ++nodeId) {
            if (!seen.at(nodeId)) {
                throw std::runtime_error(
                    "DCSweep: node_doping_file missing row for node id " +
                    std::to_string(nodeId));
            }
        }
        return model;
    }

    const std::vector<RegionDopingSpec> specs = parseDopingSpecs(cfg, scaling);
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

struct ContactConfig {
    std::unordered_map<std::string, Real> biases;
    ContactSpecsMap specs;
};

ContactConfig contactConfigFromJson(const nlohmann::json& cfg)
{
    // Route legacy ``contacts[]`` parsing through the unified boundary parser
    // so the optional ``type`` field is recognised and normalised.  For the
    // DD/Gummel paths the bias map preserves the historical semantics:
    // ohmic contacts apply the raw bias, while Dirichlet/MetalGate use the
    // same effective potential expression as the Poisson driver.  Schottky
    // contacts route through a barrier-aware ``ContactSpecsMap`` and use the
    // raw bias so the swept value reaches the metal Fermi level.  Floating
    // contacts have no DD physics yet and are rejected so a misconfigured
    // deck fails loudly instead of silently downgrading.
    ContactConfig out;
    for (const auto& spec : parseContactBoundarySpecs(cfg)) {
        switch (spec.type) {
            case ContactType::Ohmic:
                out.biases[spec.name] = spec.bias;
                break;
            case ContactType::Dirichlet:
            case ContactType::MetalGate:
                out.biases[spec.name] = effectivePoissonDirichletPotential(spec);
                break;
            case ContactType::Schottky:
                out.biases[spec.name] = spec.bias;
                out.specs[spec.name] = spec;
                break;
            case ContactType::Floating:
                throw std::runtime_error(
                    "DCSweep: contact '" + spec.name + "' has type '" +
                    toString(spec.type) + "' which is not yet implemented "
                    "for drift-diffusion sweeps. Use ohmic or schottky for now.");
        }
    }
    return out;
}


std::string vtkFilename(const std::string& prefix, int index, Real voltage)
{
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(4) << std::setfill('0') << index
        << "_" << std::setprecision(6) << std::defaultfloat << voltage << "V.vtk";
    return oss.str();
}


std::string stepDiagnostics(const DCSweepPoint& point)
{
    return "attempted_step=" + formatReal(point.attemptedStep) +
           ";accepted_step=" + formatReal(point.acceptedStep) +
           ";retry_count=" + std::to_string(point.retryCount);
}

nlohmann::json residualBlockJson(const NewtonBlockResidualInfo& blocks)
{
    return {
        {"psi", blocks.psi},
        {"phin", blocks.phin},
        {"phip", blocks.phip},
        {"combined", blocks.combined},
    };
}

nlohmann::json carrierDiagnosticsJson(const NewtonCarrierDiagnostics& diagnostics)
{
    return {
        {"positive_finite", diagnostics.positiveFinite},
        {"min_electron_density_m3", diagnostics.minElectronDensity},
        {"min_hole_density_m3", diagnostics.minHoleDensity},
        {"nonfinite_electron_count", diagnostics.nonfiniteElectronCount},
        {"nonfinite_hole_count", diagnostics.nonfiniteHoleCount},
        {"nonpositive_electron_count", diagnostics.nonpositiveElectronCount},
        {"nonpositive_hole_count", diagnostics.nonpositiveHoleCount},
    };
}

nlohmann::json lineSearchHistoryJson(const std::vector<LineSearchIterationInfo>& history)
{
    nlohmann::json out = nlohmann::json::array();
    for (const LineSearchIterationInfo& item : history) {
        out.push_back({
            {"attempt", item.attempt},
            {"damping", item.damping},
            {"residual_norm", item.residualNorm},
            {"target_residual_norm", item.targetResidualNorm},
            {"finite", item.finite},
            {"carrier_positive_finite", item.acceptedByCaller},
            {"sufficient_decrease", item.sufficientDecrease},
            {"accepted", item.accepted},
            {"rejection_reason", item.rejectionReason},
        });
    }
    return out;
}

nlohmann::json topResidualNodesJson(const std::vector<NewtonTopResidualNode>& nodes)
{
    nlohmann::json out = nlohmann::json::array();
    for (const NewtonTopResidualNode& node : nodes) {
        out.push_back({
            {"node_id", node.nodeId},
            {"x_m", node.x},
            {"y_m", node.y},
            {"x_um", node.x * 1.0e6},
            {"y_um", node.y * 1.0e6},
            {"poisson_residual", node.poissonResidual},
            {"abs_poisson_residual", node.absPoissonResidual},
            {"donors_m3", node.donors},
            {"acceptors_m3", node.acceptors},
            {"net_doping_m3", node.netDoping},
            {"donors_cm3", node.donors / 1.0e6},
            {"acceptors_cm3", node.acceptors / 1.0e6},
            {"net_doping_cm3", node.netDoping / 1.0e6},
            {"ni_eff_m3", node.effectiveIntrinsicDensity},
            {"ni_eff_cm3", node.effectiveIntrinsicDensity / 1.0e6},
        });
    }
    return out;
}

nlohmann::json newtonFailureDiagnosticsJson(const NewtonFailureDiagnostics& diagnostics)
{
    return {
        {"failure_reason", diagnostics.failureReason},
        {"failed_iteration", diagnostics.failedIteration},
        {"residual_norm", diagnostics.residualNorm},
        {"step_norm", diagnostics.stepNorm},
        {"damping_factor", diagnostics.dampingFactor},
        {"line_search_attempts", diagnostics.lineSearchAttempts},
        {"line_search_failure_reason", diagnostics.lineSearchFailureReason},
        {"block_residuals", residualBlockJson(diagnostics.blockResiduals)},
        {"carrier_diagnostics", carrierDiagnosticsJson(diagnostics.carrierDiagnostics)},
        {"line_search_history", lineSearchHistoryJson(diagnostics.lineSearchHistory)},
        {"top_poisson_residual_nodes", topResidualNodesJson(diagnostics.topPoissonResidualNodes)},
    };
}

std::filesystem::path failureDiagnosticsJsonPath(const DCSweepConfig& sweep)
{
    const std::filesystem::path csvPath(sweep.csvFile);
    return csvPath.parent_path() /
        (csvPath.stem().string() + "_newton_failure_diagnostics.json");
}


} // namespace

std::vector<DCSweepPoint> DCSweep::run(const std::string& configFile) const
{
    return runWithResult(configFile).points;
}

DCSweepResult DCSweep::runWithResult(const std::string& configFile) const
{
    std::ifstream ifs(configFile);
    if (!ifs.is_open())
        throw std::runtime_error("DCSweep: cannot open config file: " + configFile);

    nlohmann::json cfg;
    ifs >> cfg;
    const UnitScalingConfig scaling = parseUnitScalingConfig(cfg);
    const UnitScalingReferenceConfig scalingRefs = parseUnitScalingReferenceConfig(cfg);

    const std::filesystem::path cfgDir = std::filesystem::path(configFile).parent_path();
    auto resolve = [&](const std::string& path) {
        std::filesystem::path fp(path);
        if (fp.is_relative())
            fp = cfgDir / fp;
        return fp.string();
    };

    JsonMeshReader reader;
    DeviceMesh mesh = reader.read(resolve(cfg.at("mesh_file").get<std::string>()), scaling);
    mesh.buildBoxGeometry(parseBoxGeometryOptions(cfg));
    MaterialDatabase matdb;
    if (cfg.contains("materials_file"))
        matdb.loadJson(resolve(cfg.at("materials_file").get<std::string>()), scaling);
    DopingModel doping = dopingFromJson(mesh, cfg, cfgDir, scaling);
    std::vector<RegionFixedChargeSpec> fixedChargeSpecs =
        parseRegionFixedChargeSpecs(cfg, scaling);
    std::vector<InterfaceSheetChargeSpec> sheetChargeSpecs =
        parseInterfaceSheetChargeSpecs(cfg, scaling);
    ContactConfig contactConfig = contactConfigFromJson(cfg);
    std::unordered_map<std::string, Real>& baseBiases = contactConfig.biases;
    ContactSpecsMap& contactSpecs = contactConfig.specs;
    DCSweepConfig sweep = dcSweepConfigFromJson(cfg, cfgDir, scaling);

    const nlohmann::json solverCfg = cfg.value("solver", nlohmann::json::object());
    const SolverMethod solverMethod = solverMethodFromJson(cfg);
    const HybridHandoffConfig hybrid = solverMethod == SolverMethod::GummelNewton
        ? hybridHandoffConfigFromJson(solverCfg)
        : HybridHandoffConfig{};
    const DDSolutionValidationOptions validationOptions = validationOptionsFromJson(cfg);
    GummelConfig gummel;
    NewtonConfig newton;
    MobilityModelConfig mobilityConfig;
    if (solverMethod == SolverMethod::Newton) {
        newton = newtonConfigFromJson(solverCfg, scaling);
        newton.unitScalingRefs = scalingRefs;
        mobilityConfig = newton.mobility;
    } else if (solverMethod == SolverMethod::GummelNewton) {
        gummel = gummelConfigFromJson(solverCfg, scaling);
        gummel.unitScalingRefs = scalingRefs;
        newton = newtonConfigFromJson(solverCfg, scaling);
        newton.unitScalingRefs = scalingRefs;
        mobilityConfig = newton.mobility;
    } else {
        gummel = gummelConfigFromJson(solverCfg, scaling);
        gummel.unitScalingRefs = scalingRefs;
        mobilityConfig = gummel.mobility;
    }
    const Real temperature_K = (solverMethod == SolverMethod::Newton ||
                                solverMethod == SolverMethod::GummelNewton)
        ? newton.temperature_K
        : gummel.temperature_K;
    const bool recombinationDiagnosticsEnabled =
        (solverMethod == SolverMethod::Newton || solverMethod == SolverMethod::GummelNewton) &&
        newton.diagnostics;

    RecombinationModelConfig sweepRecombinationConfig;
    BandgapNarrowingConfig sweepBgnConfig;
    ImpactIonizationModelConfig sweepImpactIonizationConfig;
    if (solverMethod == SolverMethod::Newton || solverMethod == SolverMethod::GummelNewton) {
        sweepRecombinationConfig = recombinationModelConfig(newton.recombination, newton.taun, newton.taup);
        sweepRecombinationConfig.augerCn = newton.augerCn;
        sweepRecombinationConfig.augerCp = newton.augerCp;
        sweepBgnConfig = newton.bandgapNarrowing;
        sweepImpactIonizationConfig = newton.impactIonization;
    } else {
        sweepRecombinationConfig = recombinationModelConfig(gummel.recombination, gummel.taun, gummel.taup);
        sweepRecombinationConfig.augerCn = gummel.augerCn;
        sweepRecombinationConfig.augerCp = gummel.augerCp;
        sweepBgnConfig = gummel.bandgapNarrowing;
        sweepImpactIonizationConfig = gummel.impactIonization;
    }
    const std::vector<Real> effectiveNi = buildEffectiveIntrinsicDensityVector(
        mesh, matdb, doping, temperature_K, sweepBgnConfig);
    const auto sweepEdgeCells = detail::buildEdgeCellMap(mesh);
    const auto sweepCellMaterials = detail::buildCellMaterials(mesh, matdb, temperature_K);
    const auto sweepMobility = makeMobilityModel(mobilityConfig);
    const auto sweepImpactIonization = makeImpactIonizationModel(sweepImpactIonizationConfig);
    if (sweep.diagnostics.sgAvalancheEdges.enabled &&
        !detail::usesEdgeCurrentAvalancheSource(sweepImpactIonizationConfig)) {
        throw std::invalid_argument(
            "DCSweep: sweep.diagnostics.sg_avalanche_edges requires "
            "impact_ionization.generation='current_density' and "
            "impact_ionization.current_approximation='density_gradient' or 'grad_qf'.");
    }
    // Build DDScalingSpec for contact current post-processing.
    DDScalingSpec ddScaling;
    if (sweep.scaling.isUnitScaling()) {
        // Derive DDScalingSpec from the unit scaling system.
        UnitScalingSystem sc = UnitScalingSystem::fromInputs(
            temperature_K,
            (11.7 * vela::constants::eps0),
            UnitScalingSystem::autoInputsFrom(mesh, doping, matdb, 1e10),
            UnitScalingReferenceConfig{}
        );
        ddScaling.enabled = true;
        ddScaling.V0 = sc.V0();
        ddScaling.C0 = sc.C0();
        ddScaling.mu0 = sc.mu0();
        ddScaling.D0 = sc.D0();
        ddScaling.L0 = sc.L0();
        ddScaling.permittivityReference_F_per_m = (11.7 * vela::constants::eps0);
    }
    ContactCurrent contactCurrent(mesh, matdb, doping, mobilityConfig, temperature_K, ddScaling, sweepBgnConfig);
    TerminalCharge terminalCharge(mesh, doping);
    StoredCharge storedCharge(mesh);
    const bool hasMultiTerminalCharges = cfg.at("sweep").contains("terminal_charges");
    const TerminalChargeConfig& legacyChargeConfig = sweep.terminalCharges.front();
    const bool continuationDiagnosticsEnabled =
        sweep.continuation.predictor.mode != "none" ||
        sweep.continuation.branchAcceptance.terminalCurrentConsistency ||
        sweep.continuation.branchAcceptance.psiPhinJump ||
        sweep.continuation.branchAcceptance.carrierDensityJump ||
        sweep.continuation.arclength.enabled;

    CSVWriter csv(sweep.csvFile);
    std::vector<std::string> header = {"mode", "bias_contact", "bias_V",
        "current_contact", "current_electron", "current_electron_drift",
        "current_electron_diffusion", "current_hole", "current_hole_drift",
        "current_hole_diffusion", "current_total",
        "converged", "iterations", "solver_method", "gummel_iterations",
        "newton_iterations", "handoff_stage", "step_diagnostics",
        "validation_diagnostics", "failure_reason", "newton_failure_class",
        "newton_failure_diagnostics_json"};
    const bool writeUnitScaledColumns = sweep.scaling.isUnitScaling();
    if (writeUnitScaledColumns) {
        header.push_back("current_total_A_per_um");
        header.push_back("current_electron_A_per_um");
        header.push_back("current_electron_drift_A_per_um");
        header.push_back("current_electron_diffusion_A_per_um");
        header.push_back("current_hole_A_per_um");
        header.push_back("current_hole_drift_A_per_um");
        header.push_back("current_hole_diffusion_A_per_um");
    }
    std::vector<std::pair<std::string, std::string>> chargeColumns;
    std::vector<std::pair<std::string, std::string>> capacitanceColumns;
    if (sweep.storedChargeEnabled)
        header.push_back(sweep.storedCharge.perMeter ? "stored_charge_C_per_m" : "stored_charge_C");
    if (sweep.mode == CurveSweepMode::CVQuasistatic) {
        header.push_back(legacyChargeConfig.perMeter ? "charge_C_per_m" : "charge_C");
        header.push_back(legacyChargeConfig.perMeter ? "capacitance_F_per_m" : "capacitance_F");
        if (writeUnitScaledColumns && legacyChargeConfig.perMeter) {
            header.push_back("charge_C_per_um");
            header.push_back("capacitance_F_per_um");
        }
        if (hasMultiTerminalCharges) {
            for (std::size_t i = 0; i < sweep.terminalCharges.size(); ++i) {
                const TerminalChargeConfig& config = sweep.terminalCharges[i];
                const std::string name = terminalChargeName(config, i);
                const std::string chargeColumn = "charge_" + name + (config.perMeter ? "_C_per_m" : "_C");
                const std::string capColumn = "capacitance_" + capacitanceMnemonic(sweep.contact, name) +
                    (config.perMeter ? "_F_per_m" : "_F");
                chargeColumns.emplace_back(name, chargeColumn);
                capacitanceColumns.emplace_back(name, capColumn);
                header.push_back(chargeColumn);
                header.push_back(capColumn);
            }
        }
    }
    if (sweep.mode == CurveSweepMode::BVReverse) {
        header.push_back("max_electric_field_V_per_m");
        if (writeUnitScaledColumns)
            header.push_back("max_electric_field_V_per_cm");
        header.push_back("current_jump_ratio");
        header.push_back("breakdown_detected");
        header.push_back("breakdown_voltage");
        header.push_back("criterion");
        header.push_back("last_stable_bias");
        header.push_back("failed_bias");
        header.push_back("breakdown_failure_reason");
    }
    if (recombinationDiagnosticsEnabled) {
        header.push_back("recombination_max_abs_rate_m3_per_s");
        header.push_back("recombination_mean_abs_rate_m3_per_s");
        header.push_back("carrier_product_max_np_over_ni2");
    }
    if (sweep.diagnostics.transport.enabled) {
        header.push_back("mean_electron_mobility_m2_V_s");
        header.push_back("mean_hole_mobility_m2_V_s");
        header.push_back("min_electron_mobility_m2_V_s");
        header.push_back("min_hole_mobility_m2_V_s");
        header.push_back("max_electric_field_V_per_cm");
        header.push_back("mean_electron_qf_gradient_V_per_cm");
        header.push_back("mean_hole_qf_gradient_V_per_cm");
        header.push_back("mean_electron_high_field_drive_V_per_cm");
        header.push_back("mean_hole_high_field_drive_V_per_cm");
        header.push_back("min_electron_mobility_limiter");
        header.push_back("min_hole_mobility_limiter");
        header.push_back("mean_electron_mobility_limiter");
        header.push_back("mean_hole_mobility_limiter");
    }
    if (continuationDiagnosticsEnabled) {
        header.push_back("predictor_mode");
        header.push_back("predicted_initial_state");
        header.push_back("branch_acceptance_status");
        header.push_back("branch_acceptance_reason");
        header.push_back("terminal_current_consistency_ratio");
        header.push_back("psi_phin_max_jump_V");
        header.push_back("electron_density_jump_median_dex");
        header.push_back("electron_density_jump_p95_abs_dex");
        header.push_back("electron_density_jump_max_abs_dex");
        header.push_back("electron_density_jump_max_node");
    }
    csv.writeHeader(header);

    auto diagnosticContacts = [](const std::vector<std::string>& configured,
                                 const std::string& fallback) {
        if (!configured.empty())
            return configured;
        return std::vector<std::string>{fallback};
    };
    const std::vector<std::string> terminalBalanceContacts =
        diagnosticContacts(sweep.diagnostics.terminalBalance.contacts, sweep.currentContact);
    const std::vector<std::string> contactEdgeContacts =
        diagnosticContacts(sweep.diagnostics.contactEdge.contacts, sweep.currentContact);
    const std::vector<std::string> continuityBalanceContacts =
        diagnosticContacts(sweep.diagnostics.continuityBalance.contacts, sweep.currentContact);
    const std::vector<std::string> contactCurrentQfFloorContacts =
        diagnosticContacts(sweep.diagnostics.contactCurrentQfFloor.contacts, sweep.currentContact);

    std::unique_ptr<CSVWriter> terminalBalanceCsv;
    if (sweep.diagnostics.terminalBalance.enabled) {
        const std::filesystem::path diagPath(sweep.diagnostics.terminalBalance.csvFile);
        if (!diagPath.parent_path().empty())
            std::filesystem::create_directories(diagPath.parent_path());
        terminalBalanceCsv = std::make_unique<CSVWriter>(diagPath.string());
        std::vector<std::string> diagHeader = {
            "point_index",
            "bias_V",
            "bias_contact",
            "contact",
            "current_electron",
            "current_hole",
            "electron_minus_hole",
            "current_total",
            "electron_plus_hole"};
        if (writeUnitScaledColumns) {
            diagHeader.push_back("current_electron_A_per_um");
            diagHeader.push_back("current_hole_A_per_um");
            diagHeader.push_back("electron_minus_hole_A_per_um");
            diagHeader.push_back("current_total_A_per_um");
            diagHeader.push_back("electron_plus_hole_A_per_um");
        }
        diagHeader.push_back("converged");
        diagHeader.push_back("solver_method");
        diagHeader.push_back("gummel_iterations");
        diagHeader.push_back("newton_iterations");
        diagHeader.push_back("handoff_stage");
        terminalBalanceCsv->writeHeader(diagHeader);
    }

    std::unique_ptr<CSVWriter> contactEdgeCsv;
    if (sweep.diagnostics.contactEdge.enabled) {
        const std::filesystem::path diagPath(sweep.diagnostics.contactEdge.csvFile);
        if (!diagPath.parent_path().empty())
            std::filesystem::create_directories(diagPath.parent_path());
        contactEdgeCsv = std::make_unique<CSVWriter>(diagPath.string());
        std::vector<std::string> diagHeader = {
            "point_index",
            "bias_V",
            "current_contact",
            "edge_id",
            "node0",
            "node1",
            "edge_length_m",
            "edge_couple_m",
            "outward_sign",
            "bernoulli_u",
            "bernoulli_bplus",
            "bernoulli_bminus",
            "electron_branch",
            "hole_branch",
            "psi0",
            "psi1",
            "phin0",
            "phin1",
            "phip0",
            "phip1",
            "n0",
            "n1",
            "p0",
            "p1",
            "ni0",
            "ni1",
            "mun",
            "mup",
            "electron_continuity_flux",
            "hole_continuity_flux",
            "current_electron",
            "current_electron_drift",
            "current_electron_diffusion",
            "current_hole",
            "current_hole_drift",
            "current_hole_diffusion",
            "hole_qf_drop_override_applied",
            "current_total"};
        if (writeUnitScaledColumns)
            diagHeader.push_back("current_total_A_per_um");
        contactEdgeCsv->writeHeader(diagHeader);
    }

    std::unique_ptr<CSVWriter> continuityBalanceCsv;
    if (sweep.diagnostics.continuityBalance.enabled) {
        const std::filesystem::path diagPath(sweep.diagnostics.continuityBalance.csvFile);
        if (!diagPath.parent_path().empty())
            std::filesystem::create_directories(diagPath.parent_path());
        continuityBalanceCsv = std::make_unique<CSVWriter>(diagPath.string());
        continuityBalanceCsv->writeHeader({
            "point_index",
            "bias_V",
            "contact",
            "carrier",
            "contact_node",
            "interior_node",
            "contact_edge_id",
            "contact_edge_flux",
            "neighbor_edge_flux",
            "recombination_term",
            "continuity_residual",
            "interior_volume_m2",
            "qf_contact_V",
            "qf_interior_V",
            "qf_drop_V",
            "carrier_density_interior_m3"});
    }

    std::unique_ptr<CSVWriter> sgAvalancheEdgesCsv;
    if (sweep.diagnostics.sgAvalancheEdges.enabled) {
        const std::filesystem::path diagPath(sweep.diagnostics.sgAvalancheEdges.csvFile);
        if (!diagPath.parent_path().empty())
            std::filesystem::create_directories(diagPath.parent_path());
        sgAvalancheEdgesCsv = std::make_unique<CSVWriter>(diagPath.string());
        sgAvalancheEdgesCsv->writeHeader({
            "point_index",
            "bias_V",
            "edge_id",
            "node0",
            "node1",
            "x0_um",
            "y0_um",
            "x1_um",
            "y1_um",
            "edge_length_m",
            "edge_couple_m",
            "edge_area_proxy_m2",
            "electric_field_V_per_m",
            "electron_impact_field_V_per_m",
            "hole_impact_field_V_per_m",
            "electron_alpha_m_inv",
            "hole_alpha_m_inv",
            "electron_mobility_m2_V_s",
            "hole_mobility_m2_V_s",
            "electron_flux_proxy",
            "hole_flux_proxy",
            "electron_source_integral",
            "hole_source_integral",
            "edge_source_integral",
            "node0_source_integral",
            "node1_source_integral",
            "edge_class"});
    }

    std::unique_ptr<CSVWriter> newtonHistoryCsv;
    if (sweep.diagnostics.newtonHistory.enabled) {
        const std::filesystem::path diagPath(sweep.diagnostics.newtonHistory.csvFile);
        if (!diagPath.parent_path().empty())
            std::filesystem::create_directories(diagPath.parent_path());
        newtonHistoryCsv = std::make_unique<CSVWriter>(diagPath.string());
        newtonHistoryCsv->writeHeader({
            "point_index",
            "bias_V",
            "bias_contact",
            "solver_method",
            "handoff_stage",
            "iteration",
            "residual_norm",
            "relative_residual_norm",
            "raw_step_norm",
            "applied_step_norm",
            "damping_factor",
            "line_search_attempts",
            "line_search_accepted",
            "block_psi",
            "block_phin",
            "block_phip",
            "block_combined"});
    }

    std::vector<DCSweepPoint> points;
    DDSolution previousSolution;
    DDSolution predictorPreviousSolution;
    Real predictorPreviousBias = 0.0;
    Real currentSolutionBias = 0.0;
    bool hasPredictorPreviousSolution = false;
    bool hasCurrentSolutionBias = false;
    int vtkIndex = 0;

    struct SolvePointAttempt {
        bool ok = false;
        DDSolution solution;
        std::string failureReason;
        std::string validationDiagnostics;
        std::string solverMethod;
        int gummelIterations = 0;
        int newtonIterations = 0;
        std::string handoffStage;
        NewtonFailureDiagnostics newtonFailureDiagnostics;
        std::vector<NewtonIterationInfo> newtonHistory;
        bool predictedInitialState = false;
        std::string branchAcceptanceStatus;
        std::string branchAcceptanceReason;
        Real terminalCurrentConsistencyRatio = 1.0;
        Real psiPhinMaxJump_V = 0.0;
        Real electronDensityJumpMedianDex = 0.0;
        Real electronDensityJumpP95AbsDex = 0.0;
        Real electronDensityJumpMaxAbsDex = 0.0;
        Index electronDensityJumpMaxNode = -1;
        ContactCurrentEdgeOverrides contactCurrentOverrides;
    };

    if ((solverMethod == SolverMethod::Newton ||
         solverMethod == SolverMethod::GummelNewton) &&
        !contactSpecs.empty()) {
        // The Newton coupled-DD path does not yet construct Schottky-aware
        // boundary conditions.  Fail loudly so a user mis-configuring the
        // deck does not silently get the legacy Ohmic interpretation.
        std::string firstSchottky;
        for (const auto& [name, _] : contactSpecs) {
            firstSchottky = name;
            break;
        }
        throw std::runtime_error(
            "DCSweep: solver.method='newton' or 'gummel_newton' does not yet support "
            "Schottky contacts (contact '" + firstSchottky + "'). Use solver.method='gummel' for the "
            "Schottky prototype, or switch the contact to type='ohmic'.");
    }

    auto solvePoint = [&](Real voltage,
                          const DDSolution* initial,
                          bool allowContactCurrentQfFloorCapture) -> SolvePointAttempt {
        auto biases = baseBiases;
        biases[sweep.contact] = voltage;
        try {
            bool solverConverged = false;
            DDSolution sol;
            SolvePointAttempt attempt;
            if (sweep.diagnostics.contactCurrentQfFloor.enabled &&
                allowContactCurrentQfFloorCapture &&
                initial != nullptr) {
                attempt.contactCurrentOverrides = buildContactCurrentQfFloorOverrides(
                    mesh, *initial, contactCurrentQfFloorContacts);
            }
            if (solverMethod == SolverMethod::Newton) {
                NewtonResult result = initial != nullptr
                    ? runNewton(mesh, matdb, doping, biases, *initial, newton, fixedChargeSpecs, sheetChargeSpecs)
                    : runNewton(mesh, matdb, doping, biases, newton, fixedChargeSpecs, sheetChargeSpecs);
                solverConverged = result.converged;
                attempt.solverMethod = "newton";
                attempt.gummelIterations = 0;
                attempt.newtonIterations = result.iters;
                attempt.handoffStage = solverConverged ? "newton" : "newton_failed";
                attempt.newtonFailureDiagnostics = result.failureDiagnostics;
                attempt.newtonHistory = result.history;
                if (!solverConverged) {
                    attempt.failureReason = result.failureDiagnostics.failureReason.empty()
                        ? "newton_non_convergence"
                        : result.failureDiagnostics.failureReason;
                }
                sol = std::move(result.solution);
            } else if (solverMethod == SolverMethod::GummelNewton) {
                attempt.solverMethod = "gummel_newton";

                NewtonConfig handoffNewton = newton;
                handoffNewton.warmStart = true;
                if (hybrid.newtonMaxIter >= 0)
                    handoffNewton.maxIter = hybrid.newtonMaxIter;

                DDSolution gummelInitial;
                DDSolutionValidationResult gummelValidation;
                NewtonResult result;
                if (hybrid.gummelMaxIter == 0) {
                    attempt.gummelIterations = 0;
                    result = initial != nullptr
                        ? runNewton(mesh, matdb, doping, biases, *initial,
                                    handoffNewton, fixedChargeSpecs, sheetChargeSpecs)
                        : runNewton(mesh, matdb, doping, biases,
                                    handoffNewton, fixedChargeSpecs, sheetChargeSpecs);
                    attempt.newtonFailureDiagnostics = result.failureDiagnostics;
                    attempt.newtonHistory = result.history;
                } else {
                    GummelConfig initializerGummel = gummel;
                    if (hybrid.gummelMaxIter >= 0)
                        initializerGummel.maxIter = hybrid.gummelMaxIter;
                    gummelInitial = initial != nullptr
                        ? runGummel(mesh, matdb, doping, biases, contactSpecs, initializerGummel, *initial,
                                    fixedChargeSpecs, sheetChargeSpecs)
                        : runGummel(mesh, matdb, doping, biases, contactSpecs, initializerGummel,
                                    fixedChargeSpecs, sheetChargeSpecs);
                    attempt.gummelIterations = gummelInitial.iters;

                    if (!gummelInitial.converged && hybrid.requireGummelConvergence) {
                        attempt.ok = false;
                        attempt.solution = std::move(gummelInitial);
                        attempt.handoffStage = "gummel_failed";
                        attempt.failureReason = "gummel_non_convergence";
                        return attempt;
                    }

                    gummelValidation = validateDDSolution(gummelInitial, mesh, biases, validationOptions);
                    if (!gummelValidation.valid) {
                        attempt.ok = false;
                        attempt.solution = std::move(gummelInitial);
                        attempt.handoffStage = "gummel_validation_failed";
                        attempt.failureReason = "gummel_validation_failed";
                        attempt.validationDiagnostics = gummelValidation.diagnosticsString();
                        return attempt;
                    }

                    result = runNewton(mesh, matdb, doping, biases, gummelInitial,
                                       handoffNewton, fixedChargeSpecs, sheetChargeSpecs);
                    attempt.newtonFailureDiagnostics = result.failureDiagnostics;
                    attempt.newtonHistory = result.history;
                    const bool acceptedNewton =
                        result.converged && !(hybrid.newtonMaxIter == 0 && result.iters == 0);
                    if (!acceptedNewton && hybrid.fallbackToGummelOnNewtonFailure) {
                        attempt.ok = true;
                        attempt.solution = std::move(gummelInitial);
                        attempt.newtonIterations = result.iters;
                        attempt.handoffStage = "gummel_fallback";
                        attempt.validationDiagnostics = gummelValidation.diagnosticsString();
                        return attempt;
                    }
                }
                solverConverged =
                    result.converged && !(hybrid.newtonMaxIter == 0 && result.iters == 0);
                attempt.newtonIterations = result.iters;
                attempt.handoffStage = solverConverged ? "newton" : "newton_failed";
                if (!solverConverged) {
                    attempt.failureReason = result.failureDiagnostics.failureReason.empty()
                        ? "newton_non_convergence"
                        : result.failureDiagnostics.failureReason;
                }
                sol = std::move(result.solution);
            } else {
                sol = initial != nullptr
                    ? runGummel(mesh, matdb, doping, biases, contactSpecs, gummel, *initial, fixedChargeSpecs, sheetChargeSpecs)
                    : runGummel(mesh, matdb, doping, biases, contactSpecs, gummel, fixedChargeSpecs, sheetChargeSpecs);
                solverConverged = sol.converged;
                attempt.solverMethod = "gummel";
                attempt.gummelIterations = sol.iters;
                attempt.newtonIterations = 0;
                attempt.handoffStage = solverConverged ? "gummel" : "gummel_failed";
            }


            const DDSolutionValidationResult validation =
                validateDDSolution(sol, mesh, biases, validationOptions);
            attempt.ok = solverConverged && validation.valid;
            attempt.solution = std::move(sol);
            attempt.validationDiagnostics = validation.diagnosticsString();
            if (!solverConverged && attempt.failureReason.empty())
                attempt.failureReason = "non_convergence";
            else if (!validation.valid)
                attempt.failureReason = "validation_failed";
            return attempt;
        } catch (const std::exception& ex) {
            std::throw_with_nested(std::runtime_error(
                "DCSweep: solver threw at voltage " + formatReal(voltage) +
                " V: " + ex.what()));
        } catch (...) {
            std::throw_with_nested(std::runtime_error(
                "DCSweep: solver threw an unknown exception at voltage " +
                formatReal(voltage) + " V."));
        }
    };

    auto applyBranchAcceptance = [&](SolvePointAttempt& attempt) {
        const bool checkTerminalCurrent =
            sweep.continuation.branchAcceptance.terminalCurrentConsistency;
        const bool checkPsiPhinJump =
            sweep.continuation.branchAcceptance.psiPhinJump;
        const bool checkCarrierDensityJump =
            sweep.continuation.branchAcceptance.carrierDensityJump;
        if (!checkTerminalCurrent && !checkPsiPhinJump && !checkCarrierDensityJump)
            return;
        attempt.branchAcceptanceStatus = "not_checked";
        attempt.terminalCurrentConsistencyRatio = 1.0;
        attempt.psiPhinMaxJump_V = 0.0;
        attempt.electronDensityJumpMedianDex = 0.0;
        attempt.electronDensityJumpP95AbsDex = 0.0;
        attempt.electronDensityJumpMaxAbsDex = 0.0;
        attempt.electronDensityJumpMaxNode = -1;
        if (!attempt.ok)
            return;

        bool checked = false;
        if (checkTerminalCurrent) {
            checked = true;
            const ContactCurrentResult branchCurrent =
                contactCurrent.compute(attempt.solution, sweep.currentContact);
            attempt.terminalCurrentConsistencyRatio =
                terminalCurrentConsistencyRatio(branchCurrent);
            if (!std::isfinite(attempt.terminalCurrentConsistencyRatio) ||
                attempt.terminalCurrentConsistencyRatio <
                    sweep.continuation.branchAcceptance.minTerminalCurrentRatio) {
                attempt.ok = false;
                attempt.failureReason = "branch_acceptance_failed";
                attempt.branchAcceptanceStatus = "rejected";
                attempt.branchAcceptanceReason = "terminal_current_inconsistent";
                return;
            }
        }
        if (checkPsiPhinJump && hasCurrentSolutionBias) {
            checked = true;
            attempt.psiPhinMaxJump_V =
                detail::maxPsiPhinJump(previousSolution, attempt.solution);
            if (!std::isfinite(attempt.psiPhinMaxJump_V) ||
                attempt.psiPhinMaxJump_V >
                    sweep.continuation.branchAcceptance.maxPsiPhinJump_V) {
                attempt.ok = false;
                attempt.failureReason = "branch_acceptance_failed";
                attempt.branchAcceptanceStatus = "rejected";
                attempt.branchAcceptanceReason = "psi_phin_jump_exceeded";
                return;
            }
        }
        if (checkCarrierDensityJump && hasCurrentSolutionBias) {
            checked = true;
            const auto stats =
                detail::electronDensityJumpStats(previousSolution, attempt.solution);
            attempt.electronDensityJumpMedianDex = stats.medianDex;
            attempt.electronDensityJumpP95AbsDex = stats.p95AbsDex;
            attempt.electronDensityJumpMaxAbsDex = stats.maxAbsDex;
            attempt.electronDensityJumpMaxNode = stats.maxNode;
            const std::string densityJumpFailure =
                detail::electronDensityJumpAcceptanceFailure(
                    sweep.continuation.branchAcceptance, stats);
            if (!densityJumpFailure.empty()) {
                attempt.ok = false;
                attempt.failureReason = "branch_acceptance_failed";
                attempt.branchAcceptanceStatus = "rejected";
                attempt.branchAcceptanceReason = densityJumpFailure;
                return;
            }
        }
        if (!checked) {
            attempt.branchAcceptanceStatus = "not_checked";
            attempt.branchAcceptanceReason = "no_previous_solution";
            return;
        }
        attempt.branchAcceptanceStatus = "accepted";
        attempt.branchAcceptanceReason.clear();
    };

    auto solvePointWithContinuation = [&](Real voltage,
                                          const DDSolution* initial,
                                          bool allowContactCurrentQfFloorCapture,
                                          int retryCount = 0) -> SolvePointAttempt {
        if (initial == nullptr ||
            sweep.continuation.predictor.mode == "none" ||
            retryCount > 0) {
            SolvePointAttempt attempt = solvePoint(
                voltage, initial, allowContactCurrentQfFloorCapture);
            applyBranchAcceptance(attempt);
            return attempt;
        }

        DDSolution predicted = detail::predictDCSweepInitialState(
            sweep.continuation.predictor,
            hasPredictorPreviousSolution ? &predictorPreviousSolution : nullptr,
            *initial,
            predictorPreviousBias,
            currentSolutionBias,
            voltage,
            retryCount);
        SolvePointAttempt attempt = solvePoint(voltage, &predicted, false);
        if (sweep.diagnostics.contactCurrentQfFloor.enabled &&
            allowContactCurrentQfFloorCapture) {
            attempt.contactCurrentOverrides =
                buildContactCurrentQfFloorOverrides(
                    mesh, *initial, contactCurrentQfFloorContacts);
        }
        attempt.predictedInitialState = true;
        applyBranchAcceptance(attempt);
        return attempt;
    };

    auto acceptPredictorHistory = [&](const DDSolution& previous,
                                      Real              previousBias,
                                      Real              acceptedBias) {
        if (hasCurrentSolutionBias) {
            predictorPreviousSolution = previous;
            predictorPreviousBias = previousBias;
            hasPredictorPreviousSolution = true;
        }
        currentSolutionBias = acceptedBias;
        hasCurrentSolutionBias = true;
    };

    auto lastConvergedPoint = [&]() -> const DCSweepPoint* {
        for (auto it = points.rbegin(); it != points.rend(); ++it) {
            if (it->converged)
                return &(*it);
        }
        return nullptr;
    };

    const std::filesystem::path newtonFailureJsonPath = failureDiagnosticsJsonPath(sweep);
    nlohmann::json newtonFailureReports = nlohmann::json::array();
    auto writeNewtonFailureDiagnostics = [&](std::size_t pointIndex,
                                             const DCSweepPoint& point) {
        if (point.newtonFailureDiagnostics.failureReason.empty())
            return;
        nlohmann::json entry = newtonFailureDiagnosticsJson(point.newtonFailureDiagnostics);
        entry["point_index"] = pointIndex;
        entry["bias_V"] = point.bias;
        entry["bias_contact"] = sweep.contact;
        entry["solver_method"] = point.solverMethod;
        entry["handoff_stage"] = point.handoffStage;
        entry["failure_reason"] = point.failureReason;
        newtonFailureReports.push_back(std::move(entry));
        if (!newtonFailureJsonPath.parent_path().empty())
            std::filesystem::create_directories(newtonFailureJsonPath.parent_path());
        std::ofstream output(newtonFailureJsonPath);
        output << std::setw(2) << newtonFailureReports << '\n';
    };

    auto recordPoint = [&](Real voltage, const SolvePointAttempt& attempt, bool converged,
                           Real attemptedStep, Real acceptedStep, int retryCount,
                           const std::string& failureReason = std::string(),
                           const std::string& validationDiagnostics = std::string()) {
        ContactCurrentResult current{};
        ContactCurrentDetailedResult currentDetailed{};
        const DDSolution& sol = attempt.solution;
        std::unordered_map<std::string, ContactCurrentDetailedResult> detailedByContact;
        auto detailedForContact = [&](const std::string& contactName) -> const ContactCurrentDetailedResult& {
            auto it = detailedByContact.find(contactName);
            if (it == detailedByContact.end()) {
                auto inserted = detailedByContact.emplace(
                    contactName,
                    sweep.diagnostics.contactCurrentQfFloor.enabled
                        ? contactCurrent.computeDetailed(
                              sol, contactName, attempt.contactCurrentOverrides)
                        : contactCurrent.computeDetailed(sol, contactName));
                it = inserted.first;
            }
            return it->second;
        };
        if (converged) {
            if (sweep.diagnostics.contactEdge.enabled ||
                sweep.diagnostics.terminalBalance.enabled) {
                currentDetailed = detailedForContact(sweep.currentContact);
                current = currentDetailed.totals;
            } else {
                current = sweep.diagnostics.contactCurrentQfFloor.enabled
                    ? contactCurrent.computeDetailed(
                          sol, sweep.currentContact, attempt.contactCurrentOverrides).totals
                    : contactCurrent.compute(sol, sweep.currentContact);
            }
        }

        DCSweepPoint point;
        point.voltage = voltage;
        point.bias = voltage;
        point.outputCsv = sweep.csvFile;
        point.electronCurrent = current.electronCurrent;
        point.electronDriftCurrent = current.electronDriftCurrent;
        point.electronDiffusionCurrent = current.electronDiffusionCurrent;
        point.holeCurrent = current.holeCurrent;
        point.holeDriftCurrent = current.holeDriftCurrent;
        point.holeDiffusionCurrent = current.holeDiffusionCurrent;
        point.totalCurrent = current.totalCurrent;
        point.converged = converged;
        point.iterations = sol.iters;
        point.solverMethod = attempt.solverMethod;
        point.gummelIterations = attempt.gummelIterations;
        point.newtonIterations = attempt.newtonIterations;
        point.handoffStage = attempt.handoffStage;
        point.attemptedStep = attemptedStep;
        point.acceptedStep = acceptedStep;
        point.retryCount = retryCount;
        point.validationDiagnostics = validationDiagnostics;
        point.newtonFailureDiagnostics = attempt.newtonFailureDiagnostics;
        point.predictorMode = continuationDiagnosticsEnabled
            ? sweep.continuation.predictor.mode
            : std::string();
        point.predictedInitialState = attempt.predictedInitialState;
        point.branchAcceptanceStatus = attempt.branchAcceptanceStatus;
        point.branchAcceptanceReason = attempt.branchAcceptanceReason;
        point.terminalCurrentConsistencyRatio =
            attempt.terminalCurrentConsistencyRatio;
        point.psiPhinMaxJump_V = attempt.psiPhinMaxJump_V;
        point.electronDensityJumpMedianDex = attempt.electronDensityJumpMedianDex;
        point.electronDensityJumpP95AbsDex = attempt.electronDensityJumpP95AbsDex;
        point.electronDensityJumpMaxAbsDex = attempt.electronDensityJumpMaxAbsDex;
        point.electronDensityJumpMaxNode = attempt.electronDensityJumpMaxNode;
        if (!converged && !point.newtonFailureDiagnostics.failureReason.empty()) {
            point.newtonFailureClass = point.newtonFailureDiagnostics.failureReason;
            point.failureDiagnosticsJson = newtonFailureJsonPath.string();
        }
        if (converged && sweep.storedChargeEnabled) {
            point.extraFields.emplace_back(
                sweep.storedCharge.perMeter ? "stored_charge_C_per_m" : "stored_charge_C",
                storedCharge.compute(sol, sweep.storedCharge).charge);
        }
        if (!converged && sweep.mode != CurveSweepMode::BVReverse)
            point.failureReason = failureReason;
        if (converged && sweep.mode == CurveSweepMode::CVQuasistatic) {
            for (std::size_t i = 0; i < sweep.terminalCharges.size(); ++i) {
                const TerminalChargeConfig& config = sweep.terminalCharges[i];
                const std::string name = terminalChargeName(config, i);
                const Real charge = terminalCharge.compute(sol, config).charge;
                point.terminalChargeValues.emplace_back(name, charge);
                Real capacitance = 0.0;
                if (!points.empty()) {
                    const DCSweepPoint& prev = points.back();
                    const Real dV = point.bias - prev.bias;
                    auto prevIt = std::find_if(prev.terminalChargeValues.begin(), prev.terminalChargeValues.end(),
                        [&](const auto& entry) { return entry.first == name; });
                    if (dV != 0.0 && prevIt != prev.terminalChargeValues.end())
                        capacitance = (charge - prevIt->second) / dV;
                }
                point.terminalCapacitanceValues.emplace_back(name, capacitance);
            }
            if (!point.terminalChargeValues.empty())
                point.terminalCharge = point.terminalChargeValues.front().second;
            if (!point.terminalCapacitanceValues.empty())
                point.capacitance = point.terminalCapacitanceValues.front().second;
            if (hasMultiTerminalCharges) {
                for (std::size_t i = 0; i < point.terminalChargeValues.size(); ++i) {
                    point.extraFields.emplace_back(chargeColumns.at(i).second, point.terminalChargeValues.at(i).second);
                    point.extraFields.emplace_back(capacitanceColumns.at(i).second, point.terminalCapacitanceValues.at(i).second);
                }
            }
        }
        if (converged && sweep.mode == CurveSweepMode::BVReverse) {
            point.maxElectricField = maxEdgeElectricFieldMagnitude(mesh, sol.psi);
            if (!points.empty()) {
                const Real previous = std::abs(points.back().totalCurrent);
                const Real currentAbs = std::abs(point.totalCurrent);
                if (previous > 0.0)
                    point.currentJumpRatio = currentAbs / previous;
                else if (currentAbs > 0.0)
                    point.currentJumpRatio = std::numeric_limits<Real>::infinity();
            }
            if (sweep.breakdown.maxElectricField_V_per_m > 0.0 &&
                point.maxElectricField >= sweep.breakdown.maxElectricField_V_per_m) {
                point.breakdownDetected = true;
                point.breakdownVoltage = point.bias;
                point.breakdownCriterion = "max_electric_field";
            } else if (sweep.breakdown.currentJumpRatio > 0.0 &&
                       point.currentJumpRatio >= sweep.breakdown.currentJumpRatio) {
                point.breakdownDetected = true;
                point.breakdownVoltage = point.bias;
                point.breakdownCriterion = "current_jump";
            }
        } else if (!converged && sweep.mode == CurveSweepMode::BVReverse &&
                   sweep.breakdown.nonConvergenceBreakdown) {
            point.failed = true;
            point.failedBias = voltage;
            point.failureReason = failureReason.empty() ? "non_convergence" : failureReason;
            if (const DCSweepPoint* stable = lastConvergedPoint()) {
                point.lastStableBias = stable->bias;
                point.breakdownDetected = true;
                point.breakdownVoltage = point.lastStableBias;
                point.breakdownCriterion = "last_stable_before_nonconvergence";
            }
        }
        if (converged && recombinationDiagnosticsEnabled) {
            const SweepRecombinationDiagnostics diagnostics =
                computeSweepRecombinationDiagnostics(sol, effectiveNi, sweepRecombinationConfig);
            point.extraFields.emplace_back(
                "recombination_max_abs_rate_m3_per_s", diagnostics.maxAbsRate_m3_per_s);
            point.extraFields.emplace_back(
                "recombination_mean_abs_rate_m3_per_s", diagnostics.meanAbsRate_m3_per_s);
            point.extraFields.emplace_back(
                "carrier_product_max_np_over_ni2", diagnostics.maxCarrierProductRatio);
        }
        if (converged && sweep.diagnostics.transport.enabled) {
            const SweepTransportDiagnostics diagnostics =
                computeSweepTransportDiagnostics(
                    mesh, matdb, doping, mobilityConfig, temperature_K, sol);
            point.extraFields.emplace_back(
                "mean_electron_mobility_m2_V_s", diagnostics.meanElectronMobility_m2_V_s);
            point.extraFields.emplace_back(
                "mean_hole_mobility_m2_V_s", diagnostics.meanHoleMobility_m2_V_s);
            point.extraFields.emplace_back(
                "min_electron_mobility_m2_V_s", diagnostics.minElectronMobility_m2_V_s);
            point.extraFields.emplace_back(
                "min_hole_mobility_m2_V_s", diagnostics.minHoleMobility_m2_V_s);
            point.extraFields.emplace_back(
                "max_electric_field_V_per_cm", diagnostics.maxElectricField_V_per_cm);
            point.extraFields.emplace_back(
                "mean_electron_qf_gradient_V_per_cm", diagnostics.meanElectronQfGradient_V_per_cm);
            point.extraFields.emplace_back(
                "mean_hole_qf_gradient_V_per_cm", diagnostics.meanHoleQfGradient_V_per_cm);
            point.extraFields.emplace_back(
                "mean_electron_high_field_drive_V_per_cm",
                diagnostics.meanElectronHighFieldDrive_V_per_cm);
            point.extraFields.emplace_back(
                "mean_hole_high_field_drive_V_per_cm",
                diagnostics.meanHoleHighFieldDrive_V_per_cm);
            point.extraFields.emplace_back(
                "min_electron_mobility_limiter", diagnostics.minElectronMobilityLimiter);
            point.extraFields.emplace_back(
                "min_hole_mobility_limiter", diagnostics.minHoleMobilityLimiter);
            point.extraFields.emplace_back(
                "mean_electron_mobility_limiter", diagnostics.meanElectronMobilityLimiter);
            point.extraFields.emplace_back(
                "mean_hole_mobility_limiter", diagnostics.meanHoleMobilityLimiter);
        }
        if (writeUnitScaledColumns) {
            point.extraFields.emplace_back(
                "current_total_A_per_um", perMeterToPerMicron(point.totalCurrent));
            point.extraFields.emplace_back(
                "current_electron_A_per_um", perMeterToPerMicron(point.electronCurrent));
            point.extraFields.emplace_back(
                "current_electron_drift_A_per_um",
                perMeterToPerMicron(point.electronDriftCurrent));
            point.extraFields.emplace_back(
                "current_electron_diffusion_A_per_um",
                perMeterToPerMicron(point.electronDiffusionCurrent));
            point.extraFields.emplace_back(
                "current_hole_A_per_um", perMeterToPerMicron(point.holeCurrent));
            point.extraFields.emplace_back(
                "current_hole_drift_A_per_um",
                perMeterToPerMicron(point.holeDriftCurrent));
            point.extraFields.emplace_back(
                "current_hole_diffusion_A_per_um",
                perMeterToPerMicron(point.holeDiffusionCurrent));
            if (sweep.mode == CurveSweepMode::CVQuasistatic && legacyChargeConfig.perMeter) {
                point.extraFields.emplace_back(
                    "charge_C_per_um", perMeterToPerMicron(point.terminalCharge));
                point.extraFields.emplace_back(
                    "capacitance_F_per_um", perMeterToPerMicron(point.capacitance));
            }
            if (sweep.mode == CurveSweepMode::BVReverse) {
                point.extraFields.emplace_back(
                    "max_electric_field_V_per_cm",
                    voltsPerMeterToVoltsPerCm(point.maxElectricField));
            }
        }
        std::vector<std::string> row = {
            toString(sweep.mode),
            sweep.contact,
            formatReal(point.bias),
            sweep.currentContact,
            formatReal(point.electronCurrent),
            formatReal(point.electronDriftCurrent),
            formatReal(point.electronDiffusionCurrent),
            formatReal(point.holeCurrent),
            formatReal(point.holeDriftCurrent),
            formatReal(point.holeDiffusionCurrent),
            formatReal(point.totalCurrent),
            point.converged ? "1" : "0",
            std::to_string(point.iterations),
            point.solverMethod,
            std::to_string(point.gummelIterations),
            std::to_string(point.newtonIterations),
            point.handoffStage,
            stepDiagnostics(point),
            point.validationDiagnostics,
            point.failureReason,
            point.newtonFailureClass,
            point.failureDiagnosticsJson};
        if (writeUnitScaledColumns) {
            row.push_back(formatReal(perMeterToPerMicron(point.totalCurrent)));
            row.push_back(formatReal(perMeterToPerMicron(point.electronCurrent)));
            row.push_back(formatReal(perMeterToPerMicron(point.electronDriftCurrent)));
            row.push_back(formatReal(perMeterToPerMicron(point.electronDiffusionCurrent)));
            row.push_back(formatReal(perMeterToPerMicron(point.holeCurrent)));
            row.push_back(formatReal(perMeterToPerMicron(point.holeDriftCurrent)));
            row.push_back(formatReal(perMeterToPerMicron(point.holeDiffusionCurrent)));
        }
        if (sweep.storedChargeEnabled) {
            const char* storedColumn = sweep.storedCharge.perMeter ? "stored_charge_C_per_m" : "stored_charge_C";
            Real storedValue = 0.0;
            for (const auto& [name, value] : point.extraFields) {
                if (name == storedColumn) {
                    storedValue = value;
                    break;
                }
            }
            row.push_back(formatReal(storedValue));
        }
        if (sweep.mode == CurveSweepMode::CVQuasistatic) {
            row.push_back(formatReal(point.terminalCharge));
            row.push_back(formatReal(point.capacitance));
            if (writeUnitScaledColumns && legacyChargeConfig.perMeter) {
                row.push_back(formatReal(perMeterToPerMicron(point.terminalCharge)));
                row.push_back(formatReal(perMeterToPerMicron(point.capacitance)));
            }
            if (hasMultiTerminalCharges) {
                std::unordered_map<std::string, Real> extraFieldValues;
                for (const auto& [name, value] : point.extraFields)
                    extraFieldValues[name] = value;
                for (std::size_t i = 0; i < chargeColumns.size(); ++i) {
                    const auto chargeIt = extraFieldValues.find(chargeColumns.at(i).second);
                    row.push_back(formatReal(chargeIt != extraFieldValues.end() ? chargeIt->second : 0.0));
                    const auto capacitanceIt = extraFieldValues.find(capacitanceColumns.at(i).second);
                    row.push_back(formatReal(capacitanceIt != extraFieldValues.end() ? capacitanceIt->second : 0.0));
                }
            }
        }
        if (sweep.mode == CurveSweepMode::BVReverse) {
            row.push_back(formatReal(point.maxElectricField));
            if (writeUnitScaledColumns)
                row.push_back(formatReal(voltsPerMeterToVoltsPerCm(point.maxElectricField)));
            row.push_back(formatReal(point.currentJumpRatio));
            row.push_back(point.breakdownDetected ? "1" : "0");
            row.push_back(formatReal(point.breakdownVoltage));
            row.push_back(point.breakdownCriterion);
            row.push_back(formatReal(point.lastStableBias));
            row.push_back(formatReal(point.failedBias));
            row.push_back(point.failureReason);
        }
        if (recombinationDiagnosticsEnabled) {
            std::unordered_map<std::string, Real> extraFieldValues;
            for (const auto& [name, value] : point.extraFields)
                extraFieldValues[name] = value;
            for (const std::string& name : {
                     std::string("recombination_max_abs_rate_m3_per_s"),
                     std::string("recombination_mean_abs_rate_m3_per_s"),
                     std::string("carrier_product_max_np_over_ni2")}) {
                const auto it = extraFieldValues.find(name);
                row.push_back(formatReal(it != extraFieldValues.end() ? it->second : 0.0));
            }
        }
        if (sweep.diagnostics.transport.enabled) {
            std::unordered_map<std::string, Real> extraFieldValues;
            for (const auto& [name, value] : point.extraFields)
                extraFieldValues[name] = value;
            for (const std::string& name : {
                     std::string("mean_electron_mobility_m2_V_s"),
                     std::string("mean_hole_mobility_m2_V_s"),
                     std::string("min_electron_mobility_m2_V_s"),
                     std::string("min_hole_mobility_m2_V_s"),
                     std::string("max_electric_field_V_per_cm"),
                     std::string("mean_electron_qf_gradient_V_per_cm"),
                     std::string("mean_hole_qf_gradient_V_per_cm"),
                     std::string("mean_electron_high_field_drive_V_per_cm"),
                     std::string("mean_hole_high_field_drive_V_per_cm"),
                     std::string("min_electron_mobility_limiter"),
                     std::string("min_hole_mobility_limiter"),
                     std::string("mean_electron_mobility_limiter"),
                     std::string("mean_hole_mobility_limiter")}) {
                const auto it = extraFieldValues.find(name);
                row.push_back(formatReal(it != extraFieldValues.end() ? it->second : 0.0));
            }
        }
        if (continuationDiagnosticsEnabled) {
            row.push_back(point.predictorMode);
            row.push_back(point.predictedInitialState ? "1" : "0");
            row.push_back(point.branchAcceptanceStatus);
            row.push_back(point.branchAcceptanceReason);
            row.push_back(formatReal(point.terminalCurrentConsistencyRatio));
            row.push_back(formatReal(point.psiPhinMaxJump_V));
            row.push_back(formatReal(point.electronDensityJumpMedianDex));
            row.push_back(formatReal(point.electronDensityJumpP95AbsDex));
            row.push_back(formatReal(point.electronDensityJumpMaxAbsDex));
            row.push_back(formatIndexOrMinusOne(point.electronDensityJumpMaxNode));
        }
        csv.writeRow(row);
        if (!converged && !point.newtonFailureDiagnostics.failureReason.empty())
            writeNewtonFailureDiagnostics(points.size(), point);

        if (converged && newtonHistoryCsv != nullptr) {
            const std::size_t pointIndex = points.size();
            for (const NewtonIterationInfo& info : attempt.newtonHistory) {
                newtonHistoryCsv->writeRow({
                    std::to_string(pointIndex),
                    formatReal(point.bias),
                    sweep.contact,
                    point.solverMethod,
                    point.handoffStage,
                    std::to_string(info.iter),
                    formatReal(info.residualNorm),
                    formatReal(info.relativeResidualNorm),
                    formatReal(info.rawStepNorm),
                    formatReal(info.stepNorm),
                    formatReal(info.dampingFactor),
                    std::to_string(info.lineSearchAttempts),
                    info.lineSearchAccepted ? "1" : "0",
                    formatReal(info.blockResiduals.psi),
                    formatReal(info.blockResiduals.phin),
                    formatReal(info.blockResiduals.phip),
                    formatReal(info.blockResiduals.combined)});
            }
        }

        if (converged && terminalBalanceCsv != nullptr) {
            const std::size_t pointIndex = points.size();
            for (const std::string& contactName : terminalBalanceContacts) {
                const ContactCurrentResult& terminalCurrent =
                    detailedForContact(contactName).totals;
                const Real electronMinusHole =
                    terminalCurrent.electronCurrent - terminalCurrent.holeCurrent;
                const Real electronPlusHole =
                    terminalCurrent.electronCurrent + terminalCurrent.holeCurrent;
                std::vector<std::string> diagRow = {
                    std::to_string(pointIndex),
                    formatReal(point.bias),
                    sweep.contact,
                    contactName,
                    formatReal(terminalCurrent.electronCurrent),
                    formatReal(terminalCurrent.holeCurrent),
                    formatReal(electronMinusHole),
                    formatReal(terminalCurrent.totalCurrent),
                    formatReal(electronPlusHole)};
                if (writeUnitScaledColumns) {
                    diagRow.push_back(formatReal(perMeterToPerMicron(terminalCurrent.electronCurrent)));
                    diagRow.push_back(formatReal(perMeterToPerMicron(terminalCurrent.holeCurrent)));
                    diagRow.push_back(formatReal(perMeterToPerMicron(electronMinusHole)));
                    diagRow.push_back(formatReal(perMeterToPerMicron(terminalCurrent.totalCurrent)));
                    diagRow.push_back(formatReal(perMeterToPerMicron(electronPlusHole)));
                }
                diagRow.push_back(point.converged ? "1" : "0");
                diagRow.push_back(point.solverMethod);
                diagRow.push_back(std::to_string(point.gummelIterations));
                diagRow.push_back(std::to_string(point.newtonIterations));
                diagRow.push_back(point.handoffStage);
                terminalBalanceCsv->writeRow(diagRow);
            }
        }

        if (converged && contactEdgeCsv != nullptr) {
            const std::size_t pointIndex = points.size();
            for (const std::string& contactName : contactEdgeContacts) {
                const ContactCurrentDetailedResult& detailed =
                    detailedForContact(contactName);
                for (const ContactCurrentEdgeDiagnostic& edgeDiag : detailed.edges) {
                    std::vector<std::string> diagRow = {
                        std::to_string(pointIndex),
                        formatReal(point.bias),
                        contactName,
                        std::to_string(edgeDiag.edgeId),
                        std::to_string(edgeDiag.node0),
                        std::to_string(edgeDiag.node1),
                        formatReal(edgeDiag.edgeLength_m),
                        formatReal(edgeDiag.edgeCouple_m),
                        formatReal(edgeDiag.outwardSign),
                        formatReal(edgeDiag.bernoulliU),
                        formatReal(edgeDiag.bernoulliBplus),
                        formatReal(edgeDiag.bernoulliBminus),
                        edgeDiag.electronUsedQuasiFermi ? "quasi_fermi" : "density",
                        edgeDiag.holeUsedQuasiFermi ? "quasi_fermi" : "density",
                        formatReal(edgeDiag.psi0),
                        formatReal(edgeDiag.psi1),
                        formatReal(edgeDiag.phin0),
                        formatReal(edgeDiag.phin1),
                        formatReal(edgeDiag.phip0),
                        formatReal(edgeDiag.phip1),
                        formatReal(edgeDiag.n0),
                        formatReal(edgeDiag.n1),
                        formatReal(edgeDiag.p0),
                        formatReal(edgeDiag.p1),
                        formatReal(edgeDiag.ni0),
                        formatReal(edgeDiag.ni1),
                        formatReal(edgeDiag.mun),
                        formatReal(edgeDiag.mup),
                        formatReal(edgeDiag.electronContinuityFlux),
                        formatReal(edgeDiag.holeContinuityFlux),
                        formatReal(edgeDiag.electronCurrent),
                        formatReal(edgeDiag.electronDriftCurrent),
                        formatReal(edgeDiag.electronDiffusionCurrent),
                        formatReal(edgeDiag.holeCurrent),
                        formatReal(edgeDiag.holeDriftCurrent),
                        formatReal(edgeDiag.holeDiffusionCurrent),
                        edgeDiag.holeQfDropOverrideApplied ? "1" : "0",
                        formatReal(edgeDiag.totalCurrent)};
                    if (writeUnitScaledColumns)
                        diagRow.push_back(formatReal(perMeterToPerMicron(edgeDiag.totalCurrent)));
                    contactEdgeCsv->writeRow(diagRow);
                }
            }
        }

        if (converged && continuityBalanceCsv != nullptr) {
            const std::size_t pointIndex = points.size();
            const std::vector<ContinuityBalanceDiagnosticRow> balanceRows =
                computeContinuityBalanceDiagnostics(
                    mesh,
                    matdb,
                    doping,
                    mobilityConfig,
                    temperature_K,
                    sol,
                    effectiveNi,
                    sweepRecombinationConfig,
                    continuityBalanceContacts);
            for (const ContinuityBalanceDiagnosticRow& balance : balanceRows) {
                continuityBalanceCsv->writeRow({
                    std::to_string(pointIndex),
                    formatReal(point.bias),
                    balance.contact,
                    balance.carrier,
                    std::to_string(balance.contactNode),
                    std::to_string(balance.interiorNode),
                    std::to_string(balance.contactEdgeId),
                    formatReal(balance.contactEdgeFlux),
                    formatReal(balance.neighborEdgeFlux),
                    formatReal(balance.recombinationTerm),
                    formatReal(balance.continuityResidual),
                    formatReal(balance.interiorVolume_m2),
                    formatReal(balance.qfContact_V),
                    formatReal(balance.qfInterior_V),
                    formatReal(std::abs(balance.qfInterior_V - balance.qfContact_V)),
                    formatReal(balance.carrierDensityInterior_m3)});
            }
        }

        if (converged && sgAvalancheEdgesCsv != nullptr) {
            const std::size_t pointIndex = points.size();
            const std::vector<detail::SgEdgeCurrentAvalancheSourceRecord> records =
                detail::sgEdgeCurrentAvalancheSourceRecords(
                    sweepImpactIonizationConfig,
                    *sweepImpactIonization,
                    mobilityConfig,
                    *sweepMobility,
                    sweepEdgeCells,
                    mesh,
                    doping,
                    sweepCellMaterials,
                    sol.psi,
                    sol.phin,
                    sol.phip,
                    sol.n,
                    sol.p,
                    effectiveNi,
                    constants::kb * temperature_K / constants::q);
            for (const detail::SgEdgeCurrentAvalancheSourceRecord& record : records) {
                const Node& node0 = mesh.getNode(record.node0);
                const Node& node1 = mesh.getNode(record.node1);
                sgAvalancheEdgesCsv->writeRow({
                    std::to_string(pointIndex),
                    formatReal(point.bias),
                    std::to_string(record.edgeId),
                    std::to_string(record.node0),
                    std::to_string(record.node1),
                    formatReal(node0.x * 1.0e6),
                    formatReal(node0.y * 1.0e6),
                    formatReal(node1.x * 1.0e6),
                    formatReal(node1.y * 1.0e6),
                    formatReal(record.edgeLength),
                    formatReal(record.edgeCouple),
                    formatReal(record.edgeAreaProxy),
                    formatReal(record.electricField),
                    formatReal(record.electronImpactField),
                    formatReal(record.holeImpactField),
                    formatReal(record.electronAlpha),
                    formatReal(record.holeAlpha),
                    formatReal(record.electronMobility),
                    formatReal(record.holeMobility),
                    formatReal(record.electronFluxProxy),
                    formatReal(record.holeFluxProxy),
                    formatReal(record.electronSourceIntegral),
                    formatReal(record.holeSourceIntegral),
                    formatReal(record.edgeSourceIntegral),
                    formatReal(record.node0SourceIntegral),
                    formatReal(record.node1SourceIntegral),
                    classifySgAvalancheEdge(mesh, sweepEdgeCells, record)});
            }
        }

        if (converged && sweep.writeVtk) {
            point.outputVtk = vtkFilename(sweep.vtkPrefix, vtkIndex++, voltage);
            writeDDSolutionVTK(point.outputVtk,
                               mesh,
                               matdb,
                               doping,
                               sol,
                               mobilityConfig,
                               sweepRecombinationConfig,
                               sweepImpactIonizationConfig,
                               sweepBgnConfig,
                               temperature_K);
        }
        if (converged && !sweep.writeStateFile.empty())
            writeDDSolutionStateCsv(sweep.writeStateFile, sol);
        if (converged && !sweep.writeStateEveryPointPrefix.empty()) {
            const std::filesystem::path prefix(sweep.writeStateEveryPointPrefix);
            const std::filesystem::path path =
                prefix.parent_path() /
                (prefix.filename().string() + "_bias_" + biasToken(voltage) + ".csv");
            writeDDSolutionStateCsv(path, sol);
        }

        points.push_back(std::move(point));
    };

    std::unique_ptr<DDSolution> initialState;
    if (!sweep.initialStateFile.empty()) {
        initialState = std::make_unique<DDSolution>(
            readDDSolutionStateCsv(sweep.initialStateFile, mesh.numNodes()));
    }

    if (!sweep.biasPoints.empty()) {
        const DDSolution* initial = initialState.get();
        Real previousBias = 0.0;
        bool havePreviousBias = false;
        for (Real bias : sweep.biasPoints) {
            SolvePointAttempt attempt;
            bool ok = false;
            std::string failureReason;
            std::string validationDiagnostics;
            Real recordedVoltage = bias;
            Real attemptedStep = havePreviousBias ? bias - previousBias : 0.0;
            Real acceptedStep = 0.0;
            int retryCount = 0;
            try {
                if (!havePreviousBias) {
                    attempt = solvePointWithContinuation(
                        bias, initial, initial != nullptr && initial == initialState.get());
                    ok = attempt.ok;
                    acceptedStep = ok ? attemptedStep : 0.0;
                    failureReason = attempt.failureReason;
                    validationDiagnostics = attempt.validationDiagnostics;
                    if (!ok && failureReason.empty())
                        failureReason = "non_convergence";
                } else {
                    const Real direction = (bias >= previousBias) ? 1.0 : -1.0;
                    detail::DCSweepStepControlConfig pointStepControl;
                    pointStepControl.start = previousBias;
                    pointStepControl.stop = bias;
                    pointStepControl.step = direction * std::min(
                        std::abs(sweep.step), std::abs(bias - previousBias));
                    pointStepControl.minStep = sweep.minStep;
                    pointStepControl.maxStep = sweep.maxStep;
                    pointStepControl.growthFactor = sweep.growthFactor;
                    pointStepControl.shrinkFactor = sweep.shrinkFactor;
                    pointStepControl.maxRetries = sweep.maxRetries;
                    pointStepControl.stopOnFailure = sweep.stopOnFailure;

                    DDSolution localPreviousSolution = previousSolution;
                    SolvePointAttempt lastPointAttempt;
                    std::string lastPointFailureReason;
                    std::string lastPointValidationDiagnostics;

                    detail::runDCSweepStepControl(
                        pointStepControl,
                        [&](Real voltage, Real, int stepRetryCount) {
                            try {
                                SolvePointAttempt pointAttempt =
                                    solvePointWithContinuation(
                                        voltage, &localPreviousSolution, false, stepRetryCount);
                                const bool pointOk = pointAttempt.ok;
                                lastPointAttempt = std::move(pointAttempt);
                                lastPointFailureReason = pointOk
                                    ? std::string()
                                    : lastPointAttempt.failureReason;
                                lastPointValidationDiagnostics =
                                    lastPointAttempt.validationDiagnostics;
                                return pointOk;
                            } catch (const std::exception&) {
                                if (sweep.mode == CurveSweepMode::BVReverse &&
                                    sweep.breakdown.nonConvergenceBreakdown) {
                                    lastPointAttempt = SolvePointAttempt{};
                                    lastPointFailureReason = "solver_exception";
                                    lastPointValidationDiagnostics.clear();
                                    return false;
                                }
                                throw;
                            }
                        },
                        [&](const detail::DCSweepStepControlEvent& event) {
                            recordedVoltage = event.voltage;
                            attemptedStep = event.attemptedStep;
                            acceptedStep = event.acceptedStep;
                            retryCount = event.retryCount;
                            ok = event.converged;
                            attempt = std::move(lastPointAttempt);
                            validationDiagnostics = lastPointValidationDiagnostics;
                            if (!event.converged) {
                                failureReason = !lastPointFailureReason.empty()
                                    ? lastPointFailureReason
                                    : event.failureReason;
                                return;
                            }
                            failureReason.clear();
                            localPreviousSolution = attempt.solution;
                        });

                    if (ok) {
                        acceptPredictorHistory(previousSolution, currentSolutionBias, recordedVoltage);
                        previousSolution = std::move(localPreviousSolution);
                        initial = &previousSolution;
                    } else if (failureReason.empty()) {
                        failureReason = "non_convergence";
                    }
                }
            } catch (const std::exception&) {
                if (sweep.mode == CurveSweepMode::BVReverse &&
                    sweep.breakdown.nonConvergenceBreakdown) {
                    failureReason = "solver_exception";
                    validationDiagnostics.clear();
                } else {
                    throw;
                }
            }

            recordPoint(recordedVoltage, attempt, ok, attemptedStep, acceptedStep, retryCount,
                        failureReason, validationDiagnostics);
            if (!ok) {
                if (sweep.stopOnFailure)
                    return DCSweepResult{std::move(mesh), std::move(points)};
                previousBias = recordedVoltage;
                havePreviousBias = true;
                continue;
            }
            if (!havePreviousBias) {
                previousSolution = std::move(attempt.solution);
                currentSolutionBias = bias;
                hasCurrentSolutionBias = true;
                initial = &previousSolution;
            }
            previousBias = bias;
            havePreviousBias = true;
        }
        return DCSweepResult{std::move(mesh), std::move(points)};
    }

    bool startOk = false;
    SolvePointAttempt startAttempt;
    std::string startFailureReason;
    std::string startValidationDiagnostics;
    try {
        startAttempt = solvePointWithContinuation(
            sweep.start, initialState.get(), initialState != nullptr);
        startOk = startAttempt.ok;
        startFailureReason = startAttempt.failureReason;
        startValidationDiagnostics = startAttempt.validationDiagnostics;
        if (!startOk && startFailureReason.empty())
            startFailureReason = "non_convergence";
    } catch (const std::exception&) {
        if (sweep.mode == CurveSweepMode::BVReverse && sweep.breakdown.nonConvergenceBreakdown)
            startFailureReason = "solver_exception";
        else
            throw;
    }
    recordPoint(sweep.start, startAttempt, startOk, 0.0, 0.0, 0, startFailureReason,
                startValidationDiagnostics);
    if (!startOk)
        return DCSweepResult{std::move(mesh), std::move(points)};
    previousSolution = std::move(startAttempt.solution);
    currentSolutionBias = sweep.start;
    hasCurrentSolutionBias = true;

    if (sweep.mode == CurveSweepMode::BVReverse &&
        sweep.continuation.arclength.enabled) {
        if (solverMethod != SolverMethod::Newton &&
            solverMethod != SolverMethod::GummelNewton) {
            throw std::invalid_argument(
                "DCSweep: bv_reverse arclength continuation requires "
                "solver.method='newton' or 'gummel_newton'.");
        }

        const Real directionSign = sweep.stop >= sweep.start ? 1.0 : -1.0;
        auto reachedStop = [&](Real lambda) {
            return directionSign > 0.0 ? lambda >= sweep.stop : lambda <= sweep.stop;
        };
        if (reachedStop(sweep.start))
            return DCSweepResult{std::move(mesh), std::move(points)};

        NewtonSolver arclengthSolver(
            mesh,
            matdb,
            doping,
            baseBiases,
            newton,
            fixedChargeSpecs,
            sheetChargeSpecs);
        ArclengthSystem arclengthSystem = arclengthSolver.makeArclengthSystem(
            sweep.contact,
            sweep.continuation.arclength.biasFiniteDifferenceStep_V);
        PseudoArclengthContinuation continuation(
            arclengthSystem,
            sweep.continuation.arclength.core);

        ArclengthState anchor;
        anchor.x = arclengthSolver.packArclengthState(previousSolution);
        anchor.lambda = sweep.start;
        Real deltaS = sweep.continuation.arclength.core.initialStep;
        ArclengthTangent previousTangent;
        bool havePreviousTangent = false;

        constexpr int maxArclengthPoints = 10000;
        for (int pointCount = 0; pointCount < maxArclengthPoints; ++pointCount) {
            const ArclengthTangent tangent = continuation.computeTangent(
                anchor,
                directionSign,
                havePreviousTangent ? &previousTangent : nullptr);
            const ArclengthStepResult stepResult = continuation.step(
                anchor,
                tangent,
                deltaS);

            auto recordArclengthFailure = [&](Real failedBias,
                                              const std::string& reason,
                                              int retryCount) {
                SolvePointAttempt failedAttempt;
                failedAttempt.ok = false;
                failedAttempt.solution = previousSolution;
                failedAttempt.failureReason = reason.empty()
                    ? std::string("arclength_non_convergence")
                    : reason;
                failedAttempt.solverMethod = "arclength";
                failedAttempt.handoffStage = "arclength_failed";
                failedAttempt.newtonIterations = stepResult.correctorIterations;
                failedAttempt.branchAcceptanceStatus = "not_checked";
                failedAttempt.branchAcceptanceReason.clear();
                recordPoint(
                    failedBias,
                    failedAttempt,
                    false,
                    deltaS,
                    0.0,
                    retryCount,
                    failedAttempt.failureReason,
                    failedAttempt.validationDiagnostics);
            };

            if (!stepResult.converged) {
                const Real failedBias = anchor.lambda + stepResult.arclengthStep * tangent.lambdaDot;
                recordArclengthFailure(
                    failedBias,
                    stepResult.failureReason,
                    stepResult.retries);
                return DCSweepResult{std::move(mesh), std::move(points)};
            }

            SolvePointAttempt attempt;
            attempt.ok = true;
            attempt.solution = arclengthSolver.unpackArclengthState(stepResult.state.x);
            attempt.solution.converged = true;
            attempt.solution.iters = stepResult.correctorIterations;
            attempt.solverMethod = "arclength";
            attempt.gummelIterations = 0;
            attempt.newtonIterations = stepResult.correctorIterations;
            attempt.handoffStage = "arclength";
            attempt.predictedInitialState = true;

            auto arclengthBiases = baseBiases;
            arclengthBiases[sweep.contact] = stepResult.state.lambda;
            const DDSolutionValidationResult validation = validateDDSolution(
                attempt.solution,
                mesh,
                arclengthBiases,
                validationOptions);
            attempt.validationDiagnostics = validation.diagnosticsString();
            if (!validation.valid) {
                attempt.ok = false;
                attempt.failureReason = "validation_failed";
            }

            applyBranchAcceptance(attempt);
            if (!attempt.ok) {
                const Real shrunkStep = stepResult.arclengthStep *
                    sweep.continuation.arclength.core.shrinkFactor;
                if (shrunkStep >= sweep.continuation.arclength.core.minStep) {
                    deltaS = shrunkStep;
                    continue;
                }
                recordPoint(
                    stepResult.state.lambda,
                    attempt,
                    false,
                    stepResult.arclengthStep,
                    0.0,
                    stepResult.retries,
                    attempt.failureReason,
                    attempt.validationDiagnostics);
                return DCSweepResult{std::move(mesh), std::move(points)};
            }

            recordPoint(
                stepResult.state.lambda,
                attempt,
                true,
                stepResult.arclengthStep,
                stepResult.arclengthStep,
                stepResult.retries,
                std::string(),
                attempt.validationDiagnostics);
            acceptPredictorHistory(previousSolution, currentSolutionBias, stepResult.state.lambda);
            previousSolution = std::move(attempt.solution);
            currentSolutionBias = stepResult.state.lambda;
            anchor = stepResult.state;
            previousTangent = tangent;
            havePreviousTangent = true;
            deltaS = continuation.nextStep(stepResult);

            if (reachedStop(anchor.lambda))
                return DCSweepResult{std::move(mesh), std::move(points)};
        }

        throw std::runtime_error(
            "DCSweep: bv_reverse arclength continuation exceeded point budget.");
    }

    SolvePointAttempt lastStepAttempt;
    std::string lastStepFailureReason;
    std::string lastStepValidationDiagnostics;
    detail::DCSweepStepControlConfig stepControl;
    stepControl.start = sweep.start;
    stepControl.stop = sweep.stop;
    stepControl.step = sweep.step;
    stepControl.minStep = sweep.minStep;
    stepControl.maxStep = sweep.maxStep;
    stepControl.growthFactor = sweep.growthFactor;
    stepControl.shrinkFactor = sweep.shrinkFactor;
    stepControl.maxRetries = sweep.maxRetries;
    stepControl.stopOnFailure = sweep.stopOnFailure;

    detail::runDCSweepStepControl(
        stepControl,
        [&](Real voltage, Real, int stepRetryCount) {
            try {
                SolvePointAttempt attempt = solvePointWithContinuation(
                    voltage, &previousSolution, false, stepRetryCount);
                lastStepAttempt = std::move(attempt);
                lastStepFailureReason = lastStepAttempt.ok ? std::string() : lastStepAttempt.failureReason;
                lastStepValidationDiagnostics = lastStepAttempt.validationDiagnostics;
                return lastStepAttempt.ok;
            } catch (const std::exception&) {
                if (sweep.mode == CurveSweepMode::BVReverse &&
                    sweep.breakdown.nonConvergenceBreakdown) {
                    lastStepAttempt = SolvePointAttempt{};
                    lastStepFailureReason = "solver_exception";
                    lastStepValidationDiagnostics.clear();
                    return false;
                }
                throw;
            }
        },
        [&](const detail::DCSweepStepControlEvent& event) {
            std::string failureReason;
            if (!event.converged)
                failureReason = !lastStepFailureReason.empty()
                    ? lastStepFailureReason
                    : event.failureReason;
            recordPoint(event.voltage, lastStepAttempt, event.converged,
                        event.attemptedStep, event.acceptedStep, event.retryCount, failureReason,
                        lastStepValidationDiagnostics);
            if (event.converged) {
                acceptPredictorHistory(previousSolution, currentSolutionBias, event.voltage);
                previousSolution = std::move(lastStepAttempt.solution);
            }
        });

    return DCSweepResult{std::move(mesh), std::move(points)};
}

namespace detail {

namespace {

void validateDCSweepStepControlConfig(const DCSweepStepControlConfig& cfg)
{
    const auto requireFinite = [](Real value, const char* name) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(std::string("DCSweep step control: ") + name +
                                        " must be finite.");
        }
    };

    requireFinite(cfg.start, "start");
    requireFinite(cfg.stop, "stop");
    requireFinite(cfg.step, "step");
    requireFinite(cfg.minStep, "minStep");
    requireFinite(cfg.maxStep, "maxStep");
    requireFinite(cfg.growthFactor, "growthFactor");
    requireFinite(cfg.shrinkFactor, "shrinkFactor");

    if (cfg.step == 0.0)
        throw std::invalid_argument("DCSweep step control: step must be non-zero.");
    if ((cfg.stop - cfg.start) * cfg.step < 0.0) {
        throw std::invalid_argument(
            "DCSweep step control: step sign must move start toward stop.");
    }
    if (cfg.minStep <= 0.0)
        throw std::invalid_argument("DCSweep step control: minStep must be positive.");
    if (cfg.maxStep <= 0.0)
        throw std::invalid_argument("DCSweep step control: maxStep must be positive.");
    if (cfg.minStep > cfg.maxStep) {
        throw std::invalid_argument(
            "DCSweep step control: minStep must not exceed maxStep.");
    }
    if (cfg.growthFactor < 1.0) {
        throw std::invalid_argument(
            "DCSweep step control: growthFactor must be at least 1.");
    }
    if (cfg.shrinkFactor <= 0.0 || cfg.shrinkFactor >= 1.0) {
        throw std::invalid_argument(
            "DCSweep step control: shrinkFactor must be greater than 0 and less than 1.");
    }
    if (cfg.maxRetries < 0) {
        throw std::invalid_argument(
            "DCSweep step control: maxRetries must be non-negative.");
    }
}

} // namespace

void runDCSweepStepControl(const DCSweepStepControlConfig& cfg,
                           const DCSweepStepAttempt& attempt,
                           const DCSweepStepRecorder& record)
{
    validateDCSweepStepControlConfig(cfg);

    Real previousVoltage = cfg.start;
    const Real direction = (cfg.step > 0.0) ? 1.0 : -1.0;
    const Real tolerance = 1.0e-12;
    Real adaptiveStep = std::min(std::abs(cfg.step), cfg.maxStep);

    auto limitedTarget = [&](Real target, Real stepMagnitude) {
        const Real remaining = direction * (target - previousVoltage);
        const Real limited = previousVoltage + direction * std::min(stepMagnitude, remaining);
        return limited;
    };

    auto advanceToward = [&](Real target) -> bool {
        int retryCount = 0;
        Real trialStep = std::min(adaptiveStep, cfg.maxStep);
        Real lastAttempted = 0.0;
        Real lastCandidate = previousVoltage;

        while (true) {
            const Real remaining = direction * (target - previousVoltage);
            if (remaining <= tolerance)
                return true;

            const Real stepMagnitude = std::min(trialStep, remaining);
            const Real candidate = limitedTarget(target, stepMagnitude);
            const Real attemptedStep = candidate - previousVoltage;
            const bool ok = attempt(candidate, attemptedStep, retryCount);
            lastAttempted = attemptedStep;
            lastCandidate = candidate;

            if (ok) {
                record({candidate, true, attemptedStep, attemptedStep, retryCount});
                previousVoltage = candidate;
                adaptiveStep = std::min(cfg.maxStep, stepMagnitude * cfg.growthFactor);
                return true;
            }

            std::string failureReason;
            if (retryCount >= cfg.maxRetries) {
                failureReason = "non_convergence";
                record({lastCandidate, false, lastAttempted, 0.0, retryCount, failureReason});
                adaptiveStep = std::max(cfg.minStep,
                                        std::min(cfg.maxStep, std::abs(lastAttempted) * cfg.shrinkFactor));
                return false;
            }

            const Real shrunken = stepMagnitude * cfg.shrinkFactor;
            if (shrunken < cfg.minStep - std::numeric_limits<Real>::epsilon()) {
                failureReason = "min_step_exhausted";
                record({lastCandidate, false, lastAttempted, 0.0, retryCount, failureReason});
                adaptiveStep = std::max(cfg.minStep,
                                        std::min(cfg.maxStep, std::abs(lastAttempted) * cfg.shrinkFactor));
                return false;
            }

            trialStep = shrunken;
            ++retryCount;
        }
    };

    bool blockedByFailedStep = false;
    Real nominalTarget = cfg.start + cfg.step;
    while (!blockedByFailedStep && direction * (nominalTarget - cfg.stop) <= tolerance) {
        while (direction * (previousVoltage - nominalTarget) < -tolerance) {
            if (!advanceToward(nominalTarget)) {
                if (cfg.stopOnFailure)
                    return;
                blockedByFailedStep = true;
                break;
            }
        }
        nominalTarget += cfg.step;
    }

    while (!blockedByFailedStep && direction * (previousVoltage - cfg.stop) < -tolerance) {
        if (!advanceToward(cfg.stop)) {
            if (cfg.stopOnFailure)
                return;
            break;
        }
    }
}

} // namespace detail

} // namespace vela

