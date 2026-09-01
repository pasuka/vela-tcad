#include "vela/equation/FixedStateOperatorAudit.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"

#include <algorithm>
#include <array>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace vela {
namespace {

void requireFiniteVector(const VectorXd& value, const char* name, Index expectedSize)
{
    if (value.size() != static_cast<Eigen::Index>(expectedSize)) {
        throw std::invalid_argument(
            std::string("fixed-state audit state vector size mismatch for ") + name);
    }
    if (!value.allFinite()) {
        throw std::invalid_argument(
            std::string("fixed-state audit requires finite state field ") + name);
    }
}

std::vector<unsigned char> vectorBytes(const VectorXd& value)
{
    const auto* first = reinterpret_cast<const unsigned char*>(value.data());
    return {first, first + static_cast<std::size_t>(value.size()) * sizeof(Real)};
}

void requireUnchanged(const DDSolution& state,
                      const std::array<std::vector<unsigned char>, 5>& before,
                      int iters,
                      bool converged)
{
    const std::array<const VectorXd*, 5> fields = {
        &state.psi, &state.phin, &state.phip, &state.n, &state.p};
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (vectorBytes(*fields[i]) != before[i])
            throw std::logic_error("fixed-state audit mutated supplied DDSolution bytes");
    }
    if (state.iters != iters || state.converged != converged)
        throw std::logic_error("fixed-state audit mutated supplied DDSolution metadata");
}

DopingModel makeDopingModel(const VectorXd& netDoping)
{
    DopingModel model(static_cast<Index>(netDoping.size()));
    for (Index nodeId = 0; nodeId < static_cast<Index>(netDoping.size()); ++nodeId) {
        const Real value = netDoping(static_cast<Eigen::Index>(nodeId));
        model.setNodeDoping(nodeId, std::max(value, 0.0), std::max(-value, 0.0));
    }
    return model;
}

std::vector<Real> buildIntrinsicDensity(
    const DeviceMesh& mesh,
    const std::vector<Material>& cellMaterials)
{
    std::vector<Real> ni(static_cast<std::size_t>(mesh.numNodes()), 0.0);
    std::vector<bool> assigned(static_cast<std::size_t>(mesh.numNodes()), false);
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        for (Index nodeId : mesh.getCell(cellId).node_ids) {
            if (!assigned[nodeId]) {
                ni[nodeId] = cellMaterials[cellId].ni;
                assigned[nodeId] = true;
            }
        }
    }
    return ni;
}

Real signedDoubleArea(const DeviceMesh& mesh, const Cell& cell)
{
    const Node& a = mesh.getNode(cell.node_ids[0]);
    const Node& b = mesh.getNode(cell.node_ids[1]);
    const Node& c = mesh.getNode(cell.node_ids[2]);
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

} // namespace

FixedStateOperatorAuditResult evaluateFixedStateOperators(
    const DeviceMesh& mesh,
    const VectorXd& doping,
    const DDSolution& state,
    const GummelConfig& config)
{
    const MaterialDatabase materials(config.inputScaling);
    return evaluateFixedStateOperators(mesh, doping, state, config, materials);
}

FixedStateOperatorAuditResult evaluateFixedStateOperators(
    const DeviceMesh& mesh,
    const VectorXd& doping,
    const DDSolution& state,
    const GummelConfig& config,
    const MaterialDatabase& materials)
{
    return evaluateFixedStateOperators(
        mesh, makeDopingModel(doping), state, config, materials, {});
}

FixedStateOperatorAuditResult evaluateFixedStateOperators(
    const DeviceMesh& mesh,
    const DopingModel& dopingModel,
    const DDSolution& state,
    const GummelConfig& config,
    const MaterialDatabase& materials,
    FixedStateOperatorAuditOptions options)
{
    if (!options.allowGeneralTri3) {
        if (mesh.numNodes() != 6)
            throw std::invalid_argument("fixed-state audit requires exactly 6 nodes");
        if (mesh.numCells() != 4)
            throw std::invalid_argument("fixed-state audit requires exactly 4 cells");
        if (mesh.numEdges() != 9)
            throw std::invalid_argument("fixed-state audit requires exactly 9 canonical edges");
    }
    const std::set<std::array<Index, 3>> sketchTriangles = {
        std::array<Index, 3>{0, 1, 4},
        std::array<Index, 3>{1, 4, 5},
        std::array<Index, 3>{1, 3, 5},
        std::array<Index, 3>{1, 2, 3},
    };
    const std::set<std::array<Index, 3>> mirrorTriangles = {
        std::array<Index, 3>{0, 4, 5},
        std::array<Index, 3>{0, 1, 5},
        std::array<Index, 3>{1, 2, 5},
        std::array<Index, 3>{2, 3, 5},
    };
    std::set<std::array<Index, 3>> triangleKeys;

    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        const Node& node = mesh.getNode(nodeId);
        if (!std::isfinite(node.x) || !std::isfinite(node.y))
            throw std::invalid_argument("fixed-state audit requires finite node coordinates");
    }
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3)
            throw std::invalid_argument("fixed-state audit requires every cell to be Tri3");
        const Real area2 = signedDoubleArea(mesh, cell);
        if (!(area2 > 0.0) || !std::isfinite(area2)) {
            throw std::invalid_argument(
                "fixed-state audit requires counter-clockwise nondegenerate triangles");
        }
        std::array<Index, 3> key{
            cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
        std::sort(key.begin(), key.end());
        if (!triangleKeys.insert(key).second)
            throw std::invalid_argument("fixed-state audit rejects duplicate triangle keys");
    }
    if (!options.allowGeneralTri3 &&
        triangleKeys != sketchTriangles && triangleKeys != mirrorTriangles) {
        throw std::invalid_argument(
            "fixed-state audit rejects unexpected minimal6 connectivity");
    }
    if (dopingModel.numNodes() != mesh.numNodes())
        throw std::invalid_argument("fixed-state audit doping size mismatch");
    for (Index nodeId = 0; nodeId < dopingModel.numNodes(); ++nodeId) {
        if (!std::isfinite(dopingModel.donors(nodeId)) ||
            !std::isfinite(dopingModel.acceptors(nodeId))) {
            throw std::invalid_argument("fixed-state audit requires finite doping");
        }
    }

    requireFiniteVector(state.psi, "psi", mesh.numNodes());
    requireFiniteVector(state.phin, "phin", mesh.numNodes());
    requireFiniteVector(state.phip, "phip", mesh.numNodes());
    requireFiniteVector(state.n, "n", mesh.numNodes());
    requireFiniteVector(state.p, "p", mesh.numNodes());

    const std::array<std::vector<unsigned char>, 5> before = {
        vectorBytes(state.psi), vectorBytes(state.phin), vectorBytes(state.phip),
        vectorBytes(state.n), vectorBytes(state.p)};
    const int beforeIters = state.iters;
    const bool beforeConverged = state.converged;

    FixedStateOperatorAuditResult result;
    result.nodes.reserve(mesh.numNodes());

    for (Index nodeId = 0; nodeId < mesh.numNodes(); ++nodeId) {
        const Eigen::Index i = static_cast<Eigen::Index>(nodeId);
        result.nodes.push_back(FixedStateNodeRecord{
            nodeId, state.psi(i), state.phin(i), state.phip(i), state.n(i), state.p(i)});
    }

    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh, materials, config.temperature_K);
    const Real thermalVoltage =
        constants::kb * config.temperature_K / constants::q;
    const std::vector<Real> ni = detail::buildValidatedEffectiveNodeNi(
        "FixedStateOperatorAudit", mesh, materials, dopingModel,
        config.bandgapNarrowing, thermalVoltage);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    MobilityModelConfig mobilityConfig = config.mobility;
    if (isSurfaceMobilityModel(mobilityConfig)) {
        detail::updateSurfaceMobilityCellGeometry(
            mobilityConfig, mesh, edgeCells, state.psi,
            mobilityConfig.surface.coordinateFieldFactor, &cellMaterials);
    }
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(config.impactIonization);
    const Real fieldFactor =
        config.inputScaling.unitSystem().fieldFromCoordinateDeltaFactor();

    const bool qfImpact =
        detail::usesQuasiFermiAvalancheDrivingForce(config.impactIonization);
    const bool currentAlignedImpact =
        !config.impactIonization.debugRawVanOverstraeten &&
        detail::usesCurrentAlignedAvalancheDrivingForce(config.impactIonization);
    const bool qfMobility =
        config.mobility.highFieldDrivingForce == "quasi_fermi_gradient";
    const std::vector<bool> contactNodes = detail::contactNodeMask(mesh);
    const auto electronQfGradients = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) {
            const Eigen::Index i = static_cast<Eigen::Index>(node);
            return detail::electronQfForAvalancheGradient(
                state.psi(i), state.phin(i), state.n(i), ni[node],
                thermalVoltage, config.impactIonization);
        });
    const auto holeQfGradients = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) {
            const Eigen::Index i = static_cast<Eigen::Index>(node);
            return detail::holeQfForAvalancheGradient(
                state.psi(i), state.phip(i), state.p(i), ni[node],
                thermalVoltage, config.impactIonization);
        });

    std::vector<Real> electronMobilityByEdge(mesh.numEdges(), 0.0);
    std::vector<Real> holeMobilityByEdge(mesh.numEdges(), 0.0);
    std::vector<Real> electricFieldByEdge(mesh.numEdges(), 0.0);
    std::vector<Real> rawElectronFlux(mesh.numEdges(), 0.0);
    std::vector<Real> rawHoleFlux(mesh.numEdges(), 0.0);
    std::vector<Real> rawSignedElectronFlux(mesh.numEdges(), 0.0);
    std::vector<Real> rawSignedHoleFlux(mesh.numEdges(), 0.0);

    result.edges.reserve(mesh.numEdges());
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        if (edge.n0 >= edge.n1)
            throw std::invalid_argument("fixed-state audit requires canonical edge endpoints");
        if (!(edge.length > 1.0e-30) || !std::isfinite(edge.length))
            throw std::invalid_argument("fixed-state audit requires finite nonzero edge lengths");

        const Eigen::Index i = static_cast<Eigen::Index>(edge.n0);
        const Eigen::Index j = static_cast<Eigen::Index>(edge.n1);
        const Real electronQf0 = detail::electronQfForAvalancheGradient(
            state.psi(i), state.phin(i), state.n(i), ni[edge.n0],
            thermalVoltage, config.impactIonization);
        const Real electronQf1 = detail::electronQfForAvalancheGradient(
            state.psi(j), state.phin(j), state.n(j), ni[edge.n1],
            thermalVoltage, config.impactIonization);
        const Real holeQf0 = detail::holeQfForAvalancheGradient(
            state.psi(i), state.phip(i), state.p(i), ni[edge.n0],
            thermalVoltage, config.impactIonization);
        const Real holeQf1 = detail::holeQfForAvalancheGradient(
            state.psi(j), state.phip(j), state.p(j), ni[edge.n1],
            thermalVoltage, config.impactIonization);
        const Real electricField =
            std::abs((state.psi(j) - state.psi(i)) / edge.length) * fieldFactor;
        const Real electronQfField =
            std::abs((electronQf1 - electronQf0) / edge.length) * fieldFactor;
        const Real holeQfField =
            std::abs((holeQf1 - holeQf0) / edge.length) * fieldFactor;
        const Real electronCoefficientField = qfImpact
            ? detail::edgeQuasiFermiCoefficientField(
                  config.impactIonization, electronQfField, electricField,
                  edgeCells, mesh, edgeId, contactNodes, electronQfGradients,
                  fieldFactor)
            : electricField;
        const Real holeCoefficientField = qfImpact
            ? detail::edgeQuasiFermiCoefficientField(
                  config.impactIonization, holeQfField, electricField,
                  edgeCells, mesh, edgeId, contactNodes, holeQfGradients,
                  fieldFactor)
            : electricField;
        const Real electronMobilityField =
            qfMobility ? electronQfField : electricField;
        const Real holeMobilityField = qfMobility ? holeQfField : electricField;
        const detail::EdgeMobilityCarrierState mobilityState{
            state.n(i), state.n(j), state.p(i), state.p(j)};
        const Real electronMobility = detail::edgeMobility(
            edgeCells, mesh, dopingModel, *mobility, cellMaterials, edgeId,
            CarrierType::Electron, electronMobilityField, &mobilityConfig,
            &state.psi, &mobilityState);
        const Real holeMobility = detail::edgeMobility(
            edgeCells, mesh, dopingModel, *mobility, cellMaterials, edgeId,
            CarrierType::Hole, holeMobilityField, &mobilityConfig, &state.psi,
            &mobilityState);
        const SGBoltzmannFluxPolicy boltzmannFluxPolicy{
            config.bandgapNarrowing.model != "none",
            config.bandgapNarrowing.equalNiFluxEvaluation ==
                "compensated_log_expm1"};
        const Real electronRawSignedFlux = electronMobility > 0.0
            ? sgElectronBoltzmannContinuityFlux(
                  ni[edge.n0], ni[edge.n1], state.psi(i), state.psi(j),
                  state.phin(i), state.phin(j), thermalVoltage,
                  electronMobility * thermalVoltage * fieldFactor / edge.length,
                  boltzmannFluxPolicy)
            : 0.0;
        const Real holeRawSignedFlux = holeMobility > 0.0
            ? sgHoleBoltzmannContinuityFlux(
                  ni[edge.n0], ni[edge.n1], state.psi(i), state.psi(j),
                  state.phip(i), state.phip(j), thermalVoltage,
                  holeMobility * thermalVoltage * fieldFactor / edge.length,
                  boltzmannFluxPolicy)
            : 0.0;
        electronMobilityByEdge[edgeId] = electronMobility;
        holeMobilityByEdge[edgeId] = holeMobility;
        electricFieldByEdge[edgeId] = electricField;
        rawSignedElectronFlux[edgeId] = electronRawSignedFlux;
        rawSignedHoleFlux[edgeId] = holeRawSignedFlux;
        rawElectronFlux[edgeId] = std::abs(electronRawSignedFlux);
        rawHoleFlux[edgeId] = std::abs(holeRawSignedFlux);
        const Real electronMidpoint =
            detail::cellReconstructedAvalancheMidpointDensity(
                config.impactIonization, state.n(i), state.n(j),
                state.psi(i), state.psi(j), thermalVoltage);
        const Real holeMidpoint =
            detail::cellReconstructedAvalancheMidpointDensity(
                config.impactIonization, state.p(i), state.p(j),
                state.psi(j), state.psi(i), thermalVoltage);
        const Real signedElectricField =
            -(state.psi(j) - state.psi(i)) / edge.length * fieldFactor;
        Real electronImpactField = 0.0;
        Real holeImpactField = 0.0;
        Real electronAlpha = 0.0;
        Real holeAlpha = 0.0;
        if (electronMobility > 0.0) {
            electronImpactField = currentAlignedImpact
                ? detail::parallelCurrentAvalancheDrivingField(
                      signedElectricField, electronRawSignedFlux)
                : detail::electronAvalancheDrivingField(
                      config.impactIonization, electronCoefficientField,
                      electricField, 0.5 * (state.n(i) + state.n(j)));
            electronAlpha = impact->electronCoefficient(electronImpactField);
        }
        if (holeMobility > 0.0) {
            holeImpactField = currentAlignedImpact
                ? detail::parallelCurrentAvalancheDrivingField(
                      signedElectricField, holeRawSignedFlux)
                : detail::holeAvalancheDrivingField(
                      config.impactIonization, holeCoefficientField,
                      electricField, 0.5 * (state.p(i) + state.p(j)));
            holeAlpha = impact->holeCoefficient(holeImpactField);
        }

        result.edges.push_back(FixedStateEdgeRecord{
            edgeId, edge.n0, edge.n1, edge.length,
            electronRawSignedFlux, holeRawSignedFlux,
            electronMidpoint, holeMidpoint,
            electronImpactField, holeImpactField,
            electronAlpha, holeAlpha,
            detail::avalancheSourceEdgeArea(
                config.impactIonization, edgeCells, mesh, edgeId,
                &cellMaterials),
        });
    }

    const auto productionEdges = detail::sgEdgeCurrentAvalancheSourceRecords(
        config.impactIonization, *impact, config.mobility, *mobility,
        edgeCells, mesh, dopingModel, cellMaterials,
        state.psi, state.phin, state.phip, state.n, state.p, ni,
        thermalVoltage, fieldFactor);
    std::vector<bool> hasProductionRecord(mesh.numEdges(), false);
    for (const auto& production : productionEdges) {
        if (production.edgeId >= result.edges.size() ||
            hasProductionRecord[production.edgeId]) {
            throw std::runtime_error(
                "fixed-state audit received invalid production edge records");
        }
        hasProductionRecord[production.edgeId] = true;
        auto& edge = result.edges[production.edgeId];
        edge.electronRawSignedFlux = production.electronRawSignedFluxProxy;
        edge.holeRawSignedFlux = production.holeRawSignedFluxProxy;
        edge.electronImpactField = production.electronImpactField;
        edge.holeImpactField = production.holeImpactField;
        edge.electronAlpha = production.electronAlpha;
        edge.holeAlpha = production.holeAlpha;
        edge.edgeArea = production.edgeAreaProxy;
        rawSignedElectronFlux[production.edgeId] =
            production.electronRawSignedFluxProxy;
        rawSignedHoleFlux[production.edgeId] =
            production.holeRawSignedFluxProxy;
        rawElectronFlux[production.edgeId] =
            std::abs(production.electronRawSignedFluxProxy);
        rawHoleFlux[production.edgeId] =
            std::abs(production.holeRawSignedFluxProxy);
    }
    if (!options.allowGeneralTri3 && productionEdges.size() != 7) {
        throw std::runtime_error(
            "fixed-state audit expected 7 positive-couple production edges");
    }

    const auto cellEdges = detail::buildCellEdgeMap(edgeCells, mesh);
    std::vector<Real> reconstructedElectronFlux(mesh.numEdges(), 0.0);
    std::vector<Real> reconstructedHoleFlux(mesh.numEdges(), 0.0);
    const bool dualFace =
        config.impactIonization.currentMagnitudeMode == "dual_face_vector_mag";
    const bool cellVector = detail::usesCellVectorCurrentReconstructedAvalancheCurrent(
        config.impactIonization);
    const bool cellSmoothed = detail::usesCellCurrentReconstructedAvalancheCurrent(
        config.impactIonization);
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        reconstructedElectronFlux[edgeId] = dualFace
            ? detail::medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                  edgeId, rawSignedElectronFlux, edgeCells, cellEdges, mesh)
            : (cellVector
                ? detail::cellVectorReconstructedEdgeFluxMagnitude(
                      edgeId, rawSignedElectronFlux, edgeCells, cellEdges, mesh)
                : (cellSmoothed
                    ? detail::cellSmoothedEdgeFluxMagnitude(
                          edgeId, rawElectronFlux, edgeCells, cellEdges)
                    : rawElectronFlux[edgeId]));
        reconstructedHoleFlux[edgeId] = dualFace
            ? detail::medianDualFaceVectorReconstructedEdgeFluxMagnitude(
                  edgeId, rawSignedHoleFlux, edgeCells, cellEdges, mesh)
            : (cellVector
                ? detail::cellVectorReconstructedEdgeFluxMagnitude(
                      edgeId, rawSignedHoleFlux, edgeCells, cellEdges, mesh)
                : (cellSmoothed
                    ? detail::cellSmoothedEdgeFluxMagnitude(
                          edgeId, rawHoleFlux, edgeCells, cellEdges)
                    : rawHoleFlux[edgeId]));
    }

    Index supplementalCount = 0;
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        if (hasProductionRecord[edgeId])
            continue;
        ++supplementalCount;
        const auto& edge = result.edges[edgeId];
        const Real conservedFlux = detail::conservedTotalCurrentFluxMagnitude(
            edge.electronRawSignedFlux, edge.holeRawSignedFlux);
        const Real electronSelectedFlux = detail::selectAvalancheCurrentFluxProxy(
            config.impactIonization, std::abs(edge.electronRawSignedFlux),
            reconstructedElectronFlux[edgeId], electronMobilityByEdge[edgeId],
            edge.electronMidpointDensity, edge.electronImpactField,
            electricFieldByEdge[edgeId], conservedFlux);
        const Real holeSelectedFlux = detail::selectAvalancheCurrentFluxProxy(
            config.impactIonization, std::abs(edge.holeRawSignedFlux),
            reconstructedHoleFlux[edgeId], holeMobilityByEdge[edgeId],
            edge.holeMidpointDensity, edge.holeImpactField,
            electricFieldByEdge[edgeId], conservedFlux);
        if (!std::isfinite(electronSelectedFlux) ||
            !std::isfinite(holeSelectedFlux)) {
            throw std::runtime_error(
                "fixed-state audit supplemental current proxy is non-finite");
        }
    }
    if (!options.allowGeneralTri3 && supplementalCount != 2)
        throw std::runtime_error("fixed-state audit expected 2 zero-couple edges");

    const auto gradPsi = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) { return state.psi(static_cast<Eigen::Index>(node)); });
    const auto gradPhin = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) { return state.phin(static_cast<Eigen::Index>(node)); });
    const auto gradPhip = detail::computeCellScalarGradientCache(
        mesh, [&](Index node) { return state.phip(static_cast<Eigen::Index>(node)); });
    const auto triangleEdges = detail::triangleGssAvalancheSourceRecords(
        config.impactIonization, *impact, config.mobility, *mobility,
        edgeCells, mesh, dopingModel, cellMaterials,
        state.psi, state.phin, state.phip, state.n, state.p, ni,
        thermalVoltage, fieldFactor);

    result.triangles.reserve(mesh.numCells());
    for (Index cellId = 0; cellId < mesh.numCells(); ++cellId) {
        const Cell& cell = mesh.getCell(cellId);
        FixedStateTriangleRecord record{
            cellId,
            {cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]},
            signedDoubleArea(mesh, cell),
            gradPsi.gradients[cellId],
            gradPhin.gradients[cellId],
            gradPhip.gradients[cellId],
            {},
        };
        for (const auto& localEdge : triangleEdges) {
            if (localEdge.cellId == cellId)
                record.localEdges.push_back(localEdge);
        }
        std::sort(record.localEdges.begin(), record.localEdges.end(),
                  [](const auto& a, const auto& b) {
                      return a.localEdge < b.localEdge;
                  });
        if (record.localEdges.size() != 3) {
            throw std::runtime_error(
                "fixed-state audit production triangle evaluation did not return 3 local edges");
        }
        result.triangles.push_back(std::move(record));
    }

    if (detail::usesElementEdgeGssLauxAvalancheSource(
            config.impactIonization)) {
        result.elementEdgeGssLauxTriangles =
            detail::elementEdgeGssLauxAvalancheSourceRecords(
                config.impactIonization, *impact, config.mobility, *mobility,
                edgeCells, mesh, dopingModel, cellMaterials, state.psi,
                state.phin, state.phip, state.n, state.p, ni,
                thermalVoltage, fieldFactor);
    }

    result.processProbe = evaluateBVProcessProbe(
        mesh,
        dopingModel,
        state,
        config.mobility,
        config.impactIonization,
        config.bandgapNarrowing,
        materials,
        config.temperature_K,
        fieldFactor);

    requireUnchanged(state, before, beforeIters, beforeConverged);
    return result;
}

} // namespace vela
