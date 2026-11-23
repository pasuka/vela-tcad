#include "vela/physics/DensityGradientQuantumPotential.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/solver/LinearSolver.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace vela {
namespace {

Real safeExponential(Real argument)
{
    constexpr Real limit = 700.0;
    return std::exp(std::clamp(argument, -limit, limit));
}

Real safeExponentialMinusOne(Real argument)
{
    constexpr Real limit = 700.0;
    return std::expm1(std::clamp(argument, -limit, limit));
}

// GSS Eq. 9.128 in Vela's w=-Phi convention. For h=(w_j-w_i)/2,
// the decaying direction retains the exponential while the growing
// direction uses Wettstein's second-order branch to avoid overflow.
Real gssFittedJump(Real h)
{
    return h < 0.0 ? safeExponentialMinusOne(h) : h + 0.5 * h * h;
}

Real gssFittedJumpDerivative(Real h)
{
    return h < 0.0 ? safeExponential(h) : 1.0 + h;
}

std::array<Real, 3> averageBoxTriangleMeasures(
    const Real x[3], const Real y[3], Real area)
{
    std::array<Real, 3> measures{};
    // Sentaurus AverageBoxMethod uses the signed circumcentric dual measure.
    // In particular, it does not replace an obtuse triangle by the positive
    // Meyer mixed-area 1/2,1/4,1/4 fallback.  MeasureCoefficients.debug from
    // the SingleDevice mesh confirms the signed rule at the 100.73-degree
    // cell adjacent to node 1816.  Keeping the signed local measures is also
    // necessary for the element contributions to cancel the fitted stiffness
    // at the converged Sentaurus state.
    const auto squaredDistance = [&](int a, int b) {
        const Real dx = x[b] - x[a];
        const Real dy = y[b] - y[a];
        return dx * dx + dy * dy;
    };
    const auto cotangent = [&](int vertex, int a, int b) {
        const Real uax = x[a] - x[vertex];
        const Real uay = y[a] - y[vertex];
        const Real ubx = x[b] - x[vertex];
        const Real uby = y[b] - y[vertex];
        const Real cross = std::abs(uax * uby - uay * ubx);
        return (uax * ubx + uay * uby) / cross;
    };
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const int k = (i + 2) % 3;
        measures[static_cast<std::size_t>(i)] = 0.125 * (
            squaredDistance(i, k) * cotangent(j, i, k) +
            squaredDistance(i, j) * cotangent(k, i, j));
    }
    const Real sum = measures[0] + measures[1] + measures[2];
    if (sum > 0.0) {
        for (Real& measure : measures)
            measure *= area / sum;
    }
    return measures;
}

struct Eq231CellResidualRecord {
    Index cellId = 0;
    Index regionId = 0;
    bool isTransport = true;
    std::array<Index, 3> nodes{};
    std::array<Real, 3> lambda_V{};
    std::array<Real, 3> auxiliary{};
    std::array<Real, 3> stiffness{};
    std::array<Real, 3> gradientSquared{};
    std::array<Real, 3> reaction{};
    Real area_m2 = 0.0;
    Real gradientSquaredPerM2 = 0.0;
};

struct Eq231FittedEdgeRecord {
    Index cellId = 0;
    Index regionId = 0;
    Index rowNode = 0;
    Index columnNode = 0;
    Real rowW = 0.0;
    Real columnW = 0.0;
    Real halfJump = 0.0;
    bool exponentialBranch = false;
    Real stiffness = 0.0;
    Real fluxContribution = 0.0;
    Real jacobianContribution = 0.0;
};

struct Eq231ResidualDiagnostic {
    std::vector<Eq231CellResidualRecord> cells;
    std::vector<Eq231FittedEdgeRecord> fittedEdges;
    VectorXd stiffness;
    VectorXd gradientSquared;
    VectorXd reaction;
    VectorXd interfaceBoundary;
};

bool isSiliconDioxide(const std::string& material)
{
    std::string lower;
    lower.reserve(material.size());
    for (const unsigned char ch : material)
        lower.push_back(static_cast<char>(std::tolower(ch)));
    return lower == "sio2" || lower == "oxide" ||
           lower == "silicondioxide";
}

bool isPolySilicon(const std::string& material)
{
    std::string lower;
    lower.reserve(material.size());
    for (const unsigned char ch : material)
        lower.push_back(static_cast<char>(std::tolower(ch)));
    return lower == "polysilicon" || lower == "poly" ||
           lower == "poly_si" || lower == "poly-si";
}

std::string csvField(const std::string& value)
{
    if (value.find_first_of(",\"\r\n") == std::string::npos)
        return value;
    std::string quoted = "\"";
    for (const char ch : value) {
        quoted += ch;
        if (ch == '"')
            quoted += '"';
    }
    quoted += '"';
    return quoted;
}

void writeEq231ResidualDiagnostic(
    const std::string& prefixValue,
    const DeviceMesh& mesh,
    const Eq231ResidualDiagnostic& diagnostic,
    const VectorXd& state,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    bool stateIsLambda)
{
    namespace fs = std::filesystem;
    const fs::path prefix(prefixValue);
    if (!prefix.parent_path().empty())
        fs::create_directories(prefix.parent_path());
    const auto open = [&](const std::string& suffix) {
        std::ofstream stream(prefix.string() + suffix);
        if (!stream)
            throw std::runtime_error(
                "failed to create Eq. 231 residual diagnostic: " +
                prefix.string() + suffix);
        stream << std::setprecision(17);
        return stream;
    };

    using EdgeKey = std::pair<Index, Index>;
    std::map<EdgeKey, std::vector<Index>> edgeCells;
    for (const auto& record : diagnostic.cells) {
        for (int local = 0; local < 3; ++local) {
            Index a = record.nodes[local];
            Index b = record.nodes[(local + 1) % 3];
            if (a > b)
                std::swap(a, b);
            edgeCells[{a, b}].push_back(record.cellId);
        }
    }
    std::map<Index, std::set<std::pair<Index, Index>>> interfacePairsByCell;
    for (const auto& [edge, cells] : edgeCells) {
        (void)edge;
        if (cells.size() != 2)
            continue;
        const Index region0 = mesh.getCell(cells[0]).region_id;
        const Index region1 = mesh.getCell(cells[1]).region_id;
        if (region0 == region1)
            continue;
        const auto pair = std::minmax(region0, region1);
        interfacePairsByCell[cells[0]].insert(pair);
        interfacePairsByCell[cells[1]].insert(pair);
    }
    const auto pairText = [&](Index cellId) {
        std::ostringstream value;
        const auto found = interfacePairsByCell.find(cellId);
        if (found == interfacePairsByCell.end())
            return value.str();
        bool first = true;
        for (const auto& [region0, region1] : found->second) {
            if (!first)
                value << ';';
            first = false;
            value << mesh.getRegion(region0).name << '|'
                  << mesh.getRegion(region1).name;
        }
        return value.str();
    };

    auto cells = open("_cells.csv");
    cells << "cell_id,region_id,region_name,material,is_interface_cell,"
             "interface_pairs,is_transport,local_node,node_id,x_internal,y_internal,"
             "is_active,is_dirichlet,area_m2,lambda_V,w,grad_w_squared_per_m2,"
             "stiffness,gradient_squared,reaction,total\n";
    for (const auto& record : diagnostic.cells) {
        const Region& region = mesh.getRegion(record.regionId);
        const bool interfaceCell = interfacePairsByCell.contains(record.cellId);
        const std::string pairs = pairText(record.cellId);
        for (int local = 0; local < 3; ++local) {
            const Index node = record.nodes[local];
            const Node& point = mesh.getNode(node);
            const Real total = record.stiffness[local] +
                record.gradientSquared[local] + record.reaction[local];
            cells << record.cellId << ',' << record.regionId << ','
                  << csvField(region.name) << ',' << csvField(region.material) << ','
                  << (interfaceCell ? 1 : 0) << ',' << csvField(pairs) << ','
                  << (record.isTransport ? 1 : 0) << ','
                  << local << ',' << node << ',' << point.x << ',' << point.y << ','
                  << (activeNodes[node] ? 1 : 0) << ','
                  << (dirichletLambda_V.contains(node) ? 1 : 0) << ','
                  << record.area_m2 << ',' << record.lambda_V[local] << ','
                  << record.auxiliary[local] << ','
                  << record.gradientSquaredPerM2 << ','
                  << record.stiffness[local] << ','
                  << record.gradientSquared[local] << ','
                  << record.reaction[local] << ',' << total << '\n';
        }
    }

    auto nodes = open("_nodes.csv");
    nodes << "node_id,x_internal,y_internal,is_active,is_dirichlet,"
             "potential_like_V,output_lambda_V,stiffness,gradient_squared,"
             "reaction,raw_total,solver_total,interface_boundary\n";
    Index maxFreeNode = 0;
    Real maxFreeResidual = -1.0;
    for (int row = 0; row < state.size(); ++row) {
        const Index node = static_cast<Index>(row);
        const Node& point = mesh.getNode(node);
        const bool fixed = dirichletLambda_V.contains(node);
        const Real rawTotal = diagnostic.stiffness(row) +
            diagnostic.gradientSquared(row) + diagnostic.reaction(row) +
            diagnostic.interfaceBoundary(row);
        Real solverTotal = rawTotal;
        if (!activeNodes[node])
            solverTotal = state(row);
        if (fixed)
            solverTotal = state(row) - (stateIsLambda
                ? dirichletLambda_V.at(node)
                : dirichletLambda_V.at(node) - nodeOutputBandDrive_V(row));
        if (activeNodes[node] && !fixed && std::abs(rawTotal) > maxFreeResidual) {
            maxFreeResidual = std::abs(rawTotal);
            maxFreeNode = node;
        }
        nodes << node << ',' << point.x << ',' << point.y << ','
              << (activeNodes[node] ? 1 : 0) << ',' << (fixed ? 1 : 0) << ','
              << (stateIsLambda
                    ? state(row) - nodeOutputBandDrive_V(row)
                    : state(row))
              << ',' << (stateIsLambda
                    ? state(row)
                    : state(row) + nodeOutputBandDrive_V(row)) << ','
              << diagnostic.stiffness(row) << ','
              << diagnostic.gradientSquared(row) << ','
              << diagnostic.reaction(row) << ',' << rawTotal << ','
              << solverTotal << ',' << diagnostic.interfaceBoundary(row) << '\n';
    }

    if (!diagnostic.fittedEdges.empty()) {
        auto edges = open("_fitted_edges.csv");
        edges << "cell_id,region_id,region_name,material,row_node,column_node,"
                 "row_w,column_w,half_jump,branch,stiffness,flux_contribution,"
                 "jacobian_contribution\n";
        for (const auto& edge : diagnostic.fittedEdges) {
            const Region& region = mesh.getRegion(edge.regionId);
            edges << edge.cellId << ',' << edge.regionId << ','
                  << csvField(region.name) << ',' << csvField(region.material)
                  << ',' << edge.rowNode << ',' << edge.columnNode << ','
                  << edge.rowW << ',' << edge.columnW << ',' << edge.halfJump
                  << ',' << (edge.exponentialBranch ? "exponential" : "quadratic")
                  << ',' << edge.stiffness << ',' << edge.fluxContribution
                  << ',' << edge.jacobianContribution << '\n';
        }
    }

    struct RegionAggregate {
        std::set<Index> cells;
        std::set<Index> interfaceCells;
        Real stiffnessL1Free = 0.0;
        Real gradientL1Free = 0.0;
        Real reactionL1Free = 0.0;
        Real totalL1Free = 0.0;
        Real totalL1All = 0.0;
        Real interfaceTotalL1Free = 0.0;
        Real maxCellResidual = -1.0;
        Index maxCellId = 0;
    };
    std::map<Index, RegionAggregate> aggregates;
    Real globalStiffnessL1Free = 0.0;
    Real globalGradientL1Free = 0.0;
    Real globalReactionL1Free = 0.0;
    Real globalInterfaceBoundaryL1Free = 0.0;
    Real globalTotalL1Free = 0.0;
    for (const auto& record : diagnostic.cells) {
        auto& aggregate = aggregates[record.regionId];
        aggregate.cells.insert(record.cellId);
        const bool interfaceCell = interfacePairsByCell.contains(record.cellId);
        if (interfaceCell)
            aggregate.interfaceCells.insert(record.cellId);
        Real cellMaxFree = 0.0;
        for (int local = 0; local < 3; ++local) {
            const Index node = record.nodes[local];
            const Real total = record.stiffness[local] +
                record.gradientSquared[local] + record.reaction[local];
            aggregate.totalL1All += std::abs(total);
            if (!activeNodes[node] || dirichletLambda_V.contains(node))
                continue;
            aggregate.stiffnessL1Free += std::abs(record.stiffness[local]);
            aggregate.gradientL1Free += std::abs(record.gradientSquared[local]);
            aggregate.reactionL1Free += std::abs(record.reaction[local]);
            aggregate.totalL1Free += std::abs(total);
            globalStiffnessL1Free += std::abs(record.stiffness[local]);
            globalGradientL1Free += std::abs(record.gradientSquared[local]);
            globalReactionL1Free += std::abs(record.reaction[local]);
            globalTotalL1Free += std::abs(total);
            if (interfaceCell)
                aggregate.interfaceTotalL1Free += std::abs(total);
            cellMaxFree = std::max(cellMaxFree, std::abs(total));
        }
        if (cellMaxFree > aggregate.maxCellResidual) {
            aggregate.maxCellResidual = cellMaxFree;
            aggregate.maxCellId = record.cellId;
        }
    }
    for (int row = 0; row < diagnostic.interfaceBoundary.size(); ++row) {
        const Index node = static_cast<Index>(row);
        if (activeNodes[node] && !dirichletLambda_V.contains(node))
            globalInterfaceBoundaryL1Free +=
                std::abs(diagnostic.interfaceBoundary(row));
    }
    auto regions = open("_regions.csv");
    regions << "region_id,region_name,material,cell_count,interface_cell_count,"
               "stiffness_l1_free,gradient_squared_l1_free,reaction_l1_free,"
               "total_l1_free,total_l1_all,interface_total_l1_free,"
               "max_cell_residual_free,max_cell_id\n";
    for (const auto& [regionId, aggregate] : aggregates) {
        const Region& region = mesh.getRegion(regionId);
        regions << regionId << ',' << csvField(region.name) << ','
                << csvField(region.material) << ',' << aggregate.cells.size() << ','
                << aggregate.interfaceCells.size() << ','
                << aggregate.stiffnessL1Free << ','
                << aggregate.gradientL1Free << ','
                << aggregate.reactionL1Free << ','
                << aggregate.totalL1Free << ',' << aggregate.totalL1All << ','
                << aggregate.interfaceTotalL1Free << ','
                << aggregate.maxCellResidual << ',' << aggregate.maxCellId << '\n';
    }

    auto summary = open("_summary.txt");
    summary << "max_free_node=" << maxFreeNode << '\n'
            << "max_free_residual=" << maxFreeResidual << '\n'
            << "stiffness_l1_free=" << globalStiffnessL1Free << '\n'
            << "gradient_squared_l1_free=" << globalGradientL1Free << '\n'
            << "reaction_l1_free=" << globalReactionL1Free << '\n'
            << "interface_boundary_l1_free="
            << globalInterfaceBoundaryL1Free << '\n'
            << "cell_total_l1_free=" << globalTotalL1Free << '\n';
}

} // namespace

Real densityGradientStepFunction(Real x)
{
    if (!std::isfinite(x))
        throw std::invalid_argument("density-gradient step argument must be finite.");
    const Real ax = std::abs(x);
    if (ax < 1.0e-4) {
        // 1/2 + x/6 + x^2/24 + x^3/120 + x^4/720
        return 0.5 + x * (1.0 / 6.0 + x * (1.0 / 24.0 +
            x * (1.0 / 120.0 + x / 720.0)));
    }
    if (x > 700.0)
        return std::numeric_limits<Real>::max();
    return (std::expm1(x) - x) / (x * x);
}

Real densityGradientStepFunctionDerivative(Real x)
{
    if (!std::isfinite(x))
        throw std::invalid_argument("density-gradient step argument must be finite.");
    const Real ax = std::abs(x);
    if (ax < 1.0e-4) {
        // 1/6 + x/12 + x^2/40 + x^3/180 + x^4/1008
        return 1.0 / 6.0 + x * (1.0 / 12.0 + x * (1.0 / 40.0 +
            x * (1.0 / 180.0 + x / 1008.0)));
    }
    if (x > 700.0)
        return std::numeric_limits<Real>::max();
    return (std::exp(x) * (x - 2.0) + x + 2.0) / (x * x * x);
}

Real densityGradientCoefficientVm2(
    const DensityGradientQuantumPotentialConfig& config)
{
    return densityGradientCoefficientVm2(
        config.gamma, config.effectiveMassRatio);
}

Real densityGradientCoefficientVm2(Real gamma, Real effectiveMassRatio)
{
    if (!(gamma > 0.0) || !std::isfinite(gamma) ||
        !(effectiveMassRatio > 0.0) ||
        !std::isfinite(effectiveMassRatio)) {
        throw std::invalid_argument(
            "density-gradient gamma and effective mass must be finite and positive.");
    }
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    return gamma * hbar * hbar /
        (6.0 * effectiveMassRatio * constants::m0 * constants::q);
}

Real densityGradientOxidePenetrationDepthM(
    Real barrierEffectiveMassRatio,
    Real barrierHeight_V)
{
    if (!(barrierEffectiveMassRatio > 0.0) ||
        !std::isfinite(barrierEffectiveMassRatio) ||
        !(barrierHeight_V > 0.0) || !std::isfinite(barrierHeight_V)) {
        throw std::invalid_argument(
            "density-gradient oxide WKB mass and barrier must be positive.");
    }
    const Real hbar = constants::h / (2.0 * std::acos(-1.0));
    return hbar / std::sqrt(
        2.0 * barrierEffectiveMassRatio * constants::m0 *
        constants::q * barrierHeight_V);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    const VectorXd coefficient = VectorXd::Constant(
        static_cast<int>(mesh.numNodes()), densityGradientCoefficientVm2(config));
    return solveElectronDensityGradientPotential(
        mesh, classicalDensity, coefficient, activeNodes, dirichletPotential_V, {},
        thermalVoltage_V, units, std::move(config), initialPotential_V);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const VectorXd& coefficientVm2,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    return solveElectronDensityGradientPotential(
        mesh, classicalDensity, coefficientVm2, activeNodes,
        dirichletPotential_V, {}, thermalVoltage_V, units,
        std::move(config), initialPotential_V);
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotential(
    const DeviceMesh& mesh,
    const VectorXd& classicalDensity,
    const VectorXd& coefficientVm2,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletPotential_V,
    const std::vector<DensityGradientStepBoundary>& stepBoundaries,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialPotential_V)
{
    const int nodeCount = static_cast<int>(mesh.numNodes());
    if (classicalDensity.size() != nodeCount ||
        coefficientVm2.size() != nodeCount ||
        activeNodes.size() != static_cast<std::size_t>(nodeCount)) {
        throw std::invalid_argument(
            "density-gradient input sizes must match the mesh node count.");
    }
    if (!(thermalVoltage_V > 0.0) || !std::isfinite(thermalVoltage_V))
        throw std::invalid_argument("density-gradient thermal voltage must be positive.");
    if (config.maxIterations <= 0 || !(config.relativeTolerance >= 0.0) ||
        !(config.absoluteTolerance_V >= 0.0) || !(config.damping > 0.0) ||
        config.damping > 1.0 || !std::isfinite(config.damping) ||
        !(config.maxUpdate_V >= 0.0) || !std::isfinite(config.maxUpdate_V)) {
        throw std::invalid_argument("invalid density-gradient iteration configuration.");
    }

    VectorXd potential = initialPotential_V.size() == nodeCount
        ? initialPotential_V : VectorXd::Zero(nodeCount);
    for (int i = 0; i < nodeCount; ++i) {
        if (!std::isfinite(classicalDensity(i)) || classicalDensity(i) < 0.0)
            throw std::invalid_argument("density-gradient density must be finite and non-negative.");
        if (activeNodes[static_cast<std::size_t>(i)] &&
            (!(coefficientVm2(i) > 0.0) || !std::isfinite(coefficientVm2(i)))) {
            throw std::invalid_argument(
                "density-gradient material coefficient must be finite and positive.");
        }
        if (!std::isfinite(potential(i)))
            potential(i) = 0.0;
        if (!activeNodes[static_cast<std::size_t>(i)])
            potential(i) = 0.0;
    }
    for (const auto& [node, value] : dirichletPotential_V) {
        if (node >= mesh.numNodes() || !std::isfinite(value))
            throw std::invalid_argument("invalid density-gradient Dirichlet boundary.");
        potential(static_cast<int>(node)) = value;
    }
    for (const auto& boundary : stepBoundaries) {
        if (boundary.edgeId >= mesh.numEdges() ||
            !std::isfinite(boundary.barrierHeightN0_V) ||
            !std::isfinite(boundary.barrierHeightN1_V) ||
            !(boundary.barrierEffectiveMassRatio > 0.0) ||
            !std::isfinite(boundary.barrierEffectiveMassRatio) ||
            !(boundary.solvedGamma > 0.0) ||
            !std::isfinite(boundary.solvedGamma) ||
            !(boundary.theta > 0.0) || !std::isfinite(boundary.theta) ||
            !(boundary.alphaDeterminantCubeRoot > 0.0) ||
            !std::isfinite(boundary.alphaDeterminantCubeRoot) ||
            !std::isfinite(boundary.normalDrivingPotentialGradient_V_per_m)) {
            throw std::invalid_argument("invalid density-gradient step boundary.");
        }
        const Edge& edge = mesh.getEdge(boundary.edgeId);
        if (!activeNodes[edge.n0] || !activeNodes[edge.n1])
            throw std::invalid_argument(
                "density-gradient step boundary endpoints must be active.");
    }

    const Real lengthScale_m = units.lengthMPerInternal();
    const Real areaScale_m2 = lengthScale_m * lengthScale_m;
    const Real referenceDensity = std::max<Real>(classicalDensity.maxCoeff(), 1.0);
    const Real densityFloor = std::numeric_limits<Real>::min();
    LinearSolver solver;
    DensityGradientQuantumPotentialResult result;

    for (int iteration = 1; iteration <= config.maxIterations; ++iteration) {
        VectorXd rootDensity(nodeCount);
        VectorXd derivative(nodeCount);
        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)]) {
                rootDensity(i) = 0.0;
                derivative(i) = 0.0;
                continue;
            }
            const Real normalized = std::max(classicalDensity(i), densityFloor) /
                referenceDensity;
            rootDensity(i) = std::sqrt(normalized) *
                safeExponential(-potential(i) / (2.0 * thermalVoltage_V));
            derivative(i) = -rootDensity(i) / (2.0 * thermalVoltage_V);
        }

        VectorXd residual = VectorXd::Zero(nodeCount);
        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(static_cast<std::size_t>(nodeCount) * 8);
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            const int i = static_cast<int>(edge.n0);
            const int j = static_cast<int>(edge.n1);
            if (!activeNodes[edge.n0] || !activeNodes[edge.n1] ||
                edge.length <= 1.0e-30 || edge.couple <= 0.0)
                continue;
            const Real geometryWeight = edge.couple / edge.length;
            const Real laplaceWeightI = coefficientVm2(i) * geometryWeight;
            const Real laplaceWeightJ = coefficientVm2(j) * geometryWeight;
            residual(i) += laplaceWeightI * (rootDensity(j) - rootDensity(i));
            residual(j) += laplaceWeightJ * (rootDensity(i) - rootDensity(j));
            triplets.emplace_back(i, i, -laplaceWeightI * derivative(i));
            triplets.emplace_back(i, j,  laplaceWeightI * derivative(j));
            triplets.emplace_back(j, i,  laplaceWeightJ * derivative(i));
            triplets.emplace_back(j, j, -laplaceWeightJ * derivative(j));
        }

        // Sentaurus Eq. 233: analytic 1-D step condition on an interface
        // between a solved semiconductor and an unsolved nonmetal region.
        // In the transformed u equation, C*dn(u) =
        // C*u*(dn(drive)-dn(Lambda))/(2*Vt).  Linear edge shape functions
        // give each endpoint half of the physical interface-segment length.
        const Real hbar = constants::h / (2.0 * std::acos(-1.0));
        const Real kT_J = thermalVoltage_V * constants::q;
        for (const auto& boundary : stepBoundaries) {
            const Edge& edge = mesh.getEdge(boundary.edgeId);
            const Real endpointLength_m = 0.5 * edge.length * lengthScale_m;
            const Real kappa = std::sqrt(
                24.0 * boundary.barrierEffectiveMassRatio * constants::m0 * kT_J /
                (hbar * hbar * boundary.solvedGamma *
                 boundary.alphaDeterminantCubeRoot));
            const int rows[2] = {
                static_cast<int>(edge.n0), static_cast<int>(edge.n1)};
            const Real barriers[2] = {
                boundary.barrierHeightN0_V, boundary.barrierHeightN1_V};
            for (int endpoint = 0; endpoint < 2; ++endpoint) {
                const int row = rows[endpoint];
                const Real lambda = potential(row);
                const Real barrier = barriers[endpoint];
                const Real argument =
                    2.0 * boundary.theta * (lambda - barrier) / thermalVoltage_V;
                const Real step = densityGradientStepFunction(argument);
                const Real stepDerivative =
                    densityGradientStepFunctionDerivative(argument);
                const Real normalLambdaGradient =
                    kappa * (barrier - lambda) * step;
                const Real normalLambdaGradientDerivative = -kappa * (
                    step + (lambda - barrier) *
                        (2.0 * boundary.theta / thermalVoltage_V) *
                        stepDerivative);
                const Real scale = coefficientVm2(row) * endpointLength_m /
                    (2.0 * thermalVoltage_V);
                const Real transformedNormalGradient =
                    boundary.normalDrivingPotentialGradient_V_per_m -
                    normalLambdaGradient;
                residual(row) += scale * rootDensity(row) * transformedNormalGradient;
                triplets.emplace_back(
                    row, row,
                    scale * (derivative(row) * transformedNormalGradient -
                        rootDensity(row) * normalLambdaGradientDerivative));
            }
        }

        for (int i = 0; i < nodeCount; ++i) {
            if (!activeNodes[static_cast<std::size_t>(i)]) {
                residual(i) = potential(i);
                triplets.emplace_back(i, i, 1.0);
                continue;
            }
            const Real volume_m2 = mesh.getNode(static_cast<Index>(i)).volume * areaScale_m2;
            residual(i) += potential(i) * rootDensity(i) * volume_m2;
            triplets.emplace_back(
                i, i, (rootDensity(i) + potential(i) * derivative(i)) * volume_m2);
        }

        for (const auto& [node, value] : dirichletPotential_V) {
            const int row = static_cast<int>(node);
            residual(row) = potential(row) - value;
            triplets.emplace_back(row, row, 1.0);
        }

        SparseMatrixd jacobian(nodeCount, nodeCount);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
        // Replace constrained rows after summation, including edge entries.
        for (const auto& [node, value] : dirichletPotential_V) {
            (void)value;
            const int row = static_cast<int>(node);
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        for (int row = 0; row < nodeCount; ++row) {
            if (activeNodes[static_cast<std::size_t>(row)])
                continue;
            for (int outer = 0; outer < jacobian.outerSize(); ++outer) {
                for (SparseMatrixd::InnerIterator entry(jacobian, outer); entry; ++entry) {
                    if (entry.row() == row)
                        entry.valueRef() = entry.col() == row ? 1.0 : 0.0;
                }
            }
        }
        jacobian.prune(0.0);

        result.residualInfinityNorm = residual.lpNorm<Eigen::Infinity>();
        VectorXd update = solver.solve(jacobian, -residual);
        Real updateNorm = update.lpNorm<Eigen::Infinity>();
        Eigen::Index maxUpdateIndex = 0;
        update.cwiseAbs().maxCoeff(&maxUpdateIndex);
        if (config.maxUpdate_V > 0.0 && updateNorm > config.maxUpdate_V) {
            update *= config.maxUpdate_V / updateNorm;
            updateNorm = config.maxUpdate_V;
        }
        potential += config.damping * update;
        for (const auto& [node, value] : dirichletPotential_V)
            potential(static_cast<int>(node)) = value;
        result.iterations = iteration;
        result.lastUpdateInfinityNorm_V = config.damping * updateNorm;
        result.maxUpdateNode = static_cast<Index>(maxUpdateIndex);
        result.maxUpdateNodeValue_V = config.damping * update(maxUpdateIndex);
        const Real scale = std::max<Real>(potential.lpNorm<Eigen::Infinity>(), 1.0);
        if (updateNorm <= config.absoluteTolerance_V + config.relativeTolerance * scale) {
            result.converged = true;
            break;
        }
    }
    result.potentialInfinityNorm_V = potential.lpNorm<Eigen::Infinity>();
    result.potential_V = std::move(potential);
    return result;
}

static DensityGradientQuantumPotentialResult
solveElectronDensityGradientPotentialLikeGlobalImpl(
    const DeviceMesh& mesh,
    const std::vector<DensityGradientCellMaterial>& cellMaterials,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialLambda_V,
    const VectorXd& initialPotentialLike_V,
    DensityGradientBoxLinearization* linearization,
    bool linearizeBox)
{
    const int nodeCount = static_cast<int>(mesh.numNodes());
    if (nodeOutputBandDrive_V.size() != nodeCount ||
        activeNodes.size() != static_cast<std::size_t>(nodeCount) ||
        !(thermalVoltage_V > 0.0) || !std::isfinite(thermalVoltage_V)) {
        throw std::invalid_argument("invalid global density-gradient input.");
    }
    if (config.globalDiscretization != "p1_direct" &&
        config.globalDiscretization != "exponential_fitted" &&
        config.globalDiscretization != "cvfem_full" &&
        config.globalDiscretization != "p1_lambda_direct" &&
        config.globalDiscretization != "gss_potentiallike_fitted" &&
        config.globalDiscretization != "sentaurus_box" &&
        config.globalDiscretization != "conservative_sqrt_fitted" &&
        config.globalDiscretization != "gss_density_fitted") {
        throw std::invalid_argument(
            "unsupported global density-gradient discretization.");
    }
    if (config.oxideBoundary != "none" &&
        config.oxideBoundary != "devsim_wkb") {
        throw std::invalid_argument(
            "unsupported global density-gradient oxide boundary.");
    }
    if (config.oxideBoundary == "devsim_wkb" &&
        (!(config.oxideQuantumMassRatio > 0.0) ||
         !(config.oxideBarrierMassRatio > 0.0) ||
         !(config.oxideBarrierHeight_V > 0.0) ||
         !std::isfinite(config.oxideQuantumMassRatio) ||
         !std::isfinite(config.oxideBarrierMassRatio) ||
         !std::isfinite(config.oxideBarrierHeight_V))) {
        throw std::invalid_argument(
            "invalid global density-gradient oxide boundary parameters.");
    }
    const bool stateIsLambda =
        config.globalDiscretization == "gss_density_fitted" ||
        config.globalDiscretization == "p1_lambda_direct";
    const bool useGssFittedFlux =
        config.globalDiscretization == "gss_density_fitted" ||
        config.globalDiscretization == "gss_potentiallike_fitted" ||
        config.globalDiscretization == "sentaurus_box";
    const bool useSentaurusBox =
        config.globalDiscretization == "sentaurus_box";
    if (linearization != nullptr &&
        (!useSentaurusBox || std::abs(config.theta - 0.5) >= 1.0e-14 ||
         initialPotentialLike_V.size() != nodeCount || !initialPotentialLike_V.allFinite())) {
        throw std::invalid_argument(
            "DG box linearization requires a finite potential-like state and fitted sentaurus_box discretization.");
    }
    const bool useConservativeSqrtFlux =
        config.globalDiscretization == "conservative_sqrt_fitted";
    if ((useGssFittedFlux || useConservativeSqrtFlux) &&
        config.oxideBoundary != "none") {
        throw std::invalid_argument(
            "explicit-material fitted density-gradient modes cannot use an "
            "oxide WKB boundary.");
    }
    for (const auto& data : cellMaterials) {
        if (data.cellId >= mesh.numCells() || !(data.coefficientVm2 > 0.0) ||
            !std::isfinite(data.coefficientVm2) ||
            mesh.getCell(data.cellId).node_ids.size() != 3) {
            throw std::invalid_argument("invalid global density-gradient cell material.");
        }
        for (int local = 0; local < 3; ++local) {
            if (!std::isfinite(data.materialBandDrive_V[local]) ||
                !std::isfinite(data.dynamicDrivingPotential_V[local]) ||
                !std::isfinite(data.initialLambda_V[local]))
                throw std::invalid_argument("nonfinite global density-gradient cell drive.");
        }
    }
    struct OxideInterfaceSegment {
        Index edgeId = 0;
        Index oxideCellId = 0;
    };
    std::vector<OxideInterfaceSegment> oxideInterfaceSegments;
    const bool useDevsimOxideBoundary = config.oxideBoundary == "devsim_wkb";
    const Real devsimOxideCoefficientVm2 = useDevsimOxideBoundary
        ? densityGradientCoefficientVm2(1.0, config.oxideQuantumMassRatio)
        : 0.0;
    const Real oxidePenetrationDepth_m = useDevsimOxideBoundary
        ? densityGradientOxidePenetrationDepthM(
            config.oxideBarrierMassRatio, config.oxideBarrierHeight_V)
        : 0.0;
    if (useDevsimOxideBoundary) {
        using EdgeKey = std::pair<Index, Index>;
        struct EdgeSide {
            Index cellId = 0;
            bool isTransport = false;
            bool isOxide = false;
        };
        std::map<EdgeKey, std::vector<EdgeSide>> sides;
        for (const auto& data : cellMaterials) {
            const Cell& cell = mesh.getCell(data.cellId);
            const Region& region = mesh.getRegion(cell.region_id);
            for (int local = 0; local < 3; ++local) {
                Index a = cell.node_ids[local];
                Index b = cell.node_ids[(local + 1) % 3];
                if (a > b)
                    std::swap(a, b);
                sides[{a, b}].push_back({
                    data.cellId, data.isTransport,
                    isSiliconDioxide(region.material)});
            }
        }
        std::map<EdgeKey, Index> edgeIds;
        for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
            const Edge& edge = mesh.getEdge(edgeId);
            edgeIds[std::minmax(edge.n0, edge.n1)] = edgeId;
        }
        for (const auto& [key, adjacent] : sides) {
            if (adjacent.size() != 2)
                continue;
            const EdgeSide* oxide = nullptr;
            const EdgeSide* transport = nullptr;
            for (const auto& side : adjacent) {
                if (!side.isTransport && side.isOxide)
                    oxide = &side;
                if (side.isTransport)
                    transport = &side;
            }
            if (oxide != nullptr && transport != nullptr)
                oxideInterfaceSegments.push_back({edgeIds.at(key), oxide->cellId});
        }
    }

    const Real lengthScale_m = units.lengthMPerInternal();
    const Real areaScale_m2 = lengthScale_m * lengthScale_m;
    VectorXd maximumMaterialBandDrive = VectorXd::Constant(
        nodeCount, -std::numeric_limits<Real>::infinity());
    std::vector<bool> materialInterfaceNodes(
        static_cast<std::size_t>(nodeCount), false);
    std::vector<bool> interfaceHasSilicon(
        static_cast<std::size_t>(nodeCount), false);
    std::vector<bool> interfaceHasPolysilicon(
        static_cast<std::size_t>(nodeCount), false);
    std::vector<bool> insulatorReentrantCornerNodes(
        static_cast<std::size_t>(nodeCount), false);
    if (useSentaurusBox) {
        std::vector<bool> hasTransport(
            static_cast<std::size_t>(nodeCount), false);
        std::vector<bool> hasNonTransport(
            static_cast<std::size_t>(nodeCount), false);
        std::vector<std::map<std::string, int>> nonTransportCellCounts(
            static_cast<std::size_t>(nodeCount));
        for (const auto& data : cellMaterials) {
            const Cell& cell = mesh.getCell(data.cellId);
            const Region& region = mesh.getRegion(cell.region_id);
            for (int local = 0; local < 3; ++local) {
                const int node = static_cast<int>(cell.node_ids[local]);
                maximumMaterialBandDrive(node) = std::max(
                    maximumMaterialBandDrive(node),
                    data.materialBandDrive_V[local]);
                if (data.isTransport)
                    hasTransport[static_cast<std::size_t>(node)] = true;
                else {
                    hasNonTransport[static_cast<std::size_t>(node)] = true;
                    ++nonTransportCellCounts[static_cast<std::size_t>(node)]
                        [region.material];
                }
                if (data.isTransport) {
                    if (isPolySilicon(region.material))
                        interfaceHasPolysilicon[static_cast<std::size_t>(node)] = true;
                    else
                        interfaceHasSilicon[static_cast<std::size_t>(node)] = true;
                }
            }
        }
        for (int node = 0; node < nodeCount; ++node) {
            materialInterfaceNodes[static_cast<std::size_t>(node)] =
                hasTransport[static_cast<std::size_t>(node)] &&
                hasNonTransport[static_cast<std::size_t>(node)];
            const auto& counts = nonTransportCellCounts[
                static_cast<std::size_t>(node)];
            if (!hasTransport[static_cast<std::size_t>(node)] &&
                counts.size() == 2) {
                const int first = counts.begin()->second;
                const int second = std::next(counts.begin())->second;
                insulatorReentrantCornerNodes[static_cast<std::size_t>(node)] =
                    first + second == 6 &&
                    std::min(first, second) == 2 &&
                    std::max(first, second) == 4;
            }
        }
    }
    VectorXd stateBandDrive_V = nodeOutputBandDrive_V;
    if (useSentaurusBox) {
        for (int node = 0; node < nodeCount; ++node) {
            if (std::isfinite(maximumMaterialBandDrive(node)))
                stateBandDrive_V(node) = maximumMaterialBandDrive(node);
        }
    }
    VectorXd continuousPotential = VectorXd::Zero(nodeCount);
    if (!stateIsLambda)
        continuousPotential = -stateBandDrive_V;
    if (stateIsLambda) {
        if (initialLambda_V.size() == nodeCount)
            continuousPotential = initialLambda_V;
        else if (initialPotentialLike_V.size() == nodeCount)
            continuousPotential = initialPotentialLike_V + stateBandDrive_V;
    } else if (initialPotentialLike_V.size() == nodeCount) {
        continuousPotential = initialPotentialLike_V;
    } else if (initialLambda_V.size() == nodeCount) {
        continuousPotential = initialLambda_V - stateBandDrive_V;
    }
    std::vector<bool> constrainedRows(activeNodes.size(), false);
    for (int i = 0; i < nodeCount; ++i) {
        if (!std::isfinite(continuousPotential(i)))
            continuousPotential(i) = stateIsLambda
                ? 0.0 : -stateBandDrive_V(i);
        if (!activeNodes[static_cast<std::size_t>(i)]) {
            constrainedRows[static_cast<std::size_t>(i)] = true;
            if (linearization == nullptr)
                continuousPotential(i) = 0.0;
        }
    }
    for (const auto& [node, lambda] : dirichletLambda_V) {
        if (node >= mesh.numNodes() || !std::isfinite(lambda))
            throw std::invalid_argument(
                "invalid global density-gradient Dirichlet boundary.");
        constrainedRows[node] = true;
        if (linearization == nullptr)
            continuousPotential(static_cast<int>(node)) = stateIsLambda
                ? lambda
                : lambda - stateBandDrive_V(static_cast<int>(node));
    }
    VectorXd conservativeRowLogReference = VectorXd::Zero(nodeCount);
    if (useConservativeSqrtFlux) {
        conservativeRowLogReference = VectorXd::Constant(
            nodeCount, -std::numeric_limits<Real>::infinity());
        for (const auto& data : cellMaterials) {
            const Cell& cell = mesh.getCell(data.cellId);
            for (int local = 0; local < 3; ++local) {
                const int node = static_cast<int>(cell.node_ids[local]);
                const Real initialW =
                    (data.dynamicDrivingPotential_V[local] -
                     continuousPotential(node)) / thermalVoltage_V;
                conservativeRowLogReference(node) = std::max(
                    conservativeRowLogReference(node), 0.5 * initialW);
            }
        }
        for (int node = 0; node < nodeCount; ++node) {
            if (!std::isfinite(conservativeRowLogReference(node)))
                conservativeRowLogReference(node) = 0.0;
        }
    }
    LinearSolver solver;
    DensityGradientQuantumPotentialResult result;
    Real initialResidualNorm = -1.0;
    std::vector<Eigen::Triplet<Real>> transportDrive;
    std::vector<Eigen::Triplet<Real>> insulatorDrive;
    if (linearizeBox)
        linearization->materialDriveDiagonal = VectorXd::Zero(nodeCount);
    auto assemble = [&](const VectorXd& state,
                        VectorXd& residual,
                        std::vector<Eigen::Triplet<Real>>* triplets,
                        Eq231ResidualDiagnostic* diagnostic) {
        residual = VectorXd::Zero(nodeCount);
        if (triplets != nullptr) {
            triplets->clear();
            triplets->reserve(cellMaterials.size() * 18 + nodeCount);
        }
        if (diagnostic != nullptr) {
            diagnostic->cells.clear();
            diagnostic->cells.reserve(cellMaterials.size());
            diagnostic->fittedEdges.clear();
            diagnostic->fittedEdges.reserve(cellMaterials.size() * 6);
            diagnostic->stiffness = VectorXd::Zero(nodeCount);
            diagnostic->gradientSquared = VectorXd::Zero(nodeCount);
            diagnostic->reaction = VectorXd::Zero(nodeCount);
            diagnostic->interfaceBoundary = VectorXd::Zero(nodeCount);
        }
        auto addJacobian = [&](int row, int col, Real value) {
            if (triplets != nullptr && !constrainedRows[static_cast<std::size_t>(row)])
                triplets->emplace_back(row, col, value);
        };
        for (const auto& data : cellMaterials) {
            const Cell& cell = mesh.getCell(data.cellId);
            const Index nodes[3] = {
                cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
            const Node& p0 = mesh.getNode(nodes[0]);
            const Node& p1 = mesh.getNode(nodes[1]);
            const Node& p2 = mesh.getNode(nodes[2]);
            const Real twiceAreaInternal = std::abs(
                (p1.x - p0.x) * (p2.y - p0.y) -
                (p1.y - p0.y) * (p2.x - p0.x));
            if (!(twiceAreaInternal > 0.0))
                continue;
            const Real cellArea_m2 = 0.5 * twiceAreaInternal * areaScale_m2;
            const Region& cellRegion = mesh.getRegion(cell.region_id);
            const Real cellCoefficientVm2 =
                useDevsimOxideBoundary &&
                    isSiliconDioxide(cellRegion.material)
                ? devsimOxideCoefficientVm2
                : data.coefficientVm2;
            Real lambda[3]{};
            Real w[3]{};
            for (int local = 0; local < 3; ++local) {
                const int node = static_cast<int>(nodes[local]);
                if (stateIsLambda) {
                    lambda[local] = state(node);
                    // The existing w convention is minus the GSS Phi:
                    // w=-Phi=-(Eqc-EFn)/kT. Lambda is continuous, whereas
                    // the band/DOS contribution and therefore w remain
                    // material-side traces at an internal interface.
                    w[local] = (data.dynamicDrivingPotential_V[local] +
                        data.materialBandDrive_V[local] - state(node)) /
                        thermalVoltage_V;
                } else {
                    lambda[local] =
                        state(node) + data.materialBandDrive_V[local];
                    w[local] = (data.dynamicDrivingPotential_V[local] -
                        state(node)) / thermalVoltage_V;
                    // Direct NewtonPlot and adjacent-row probes recover
                    // separate oxide, Silicon, and PolySilicon traces.
                    // Apply the side trace in w so the residual and analytic
                    // Jacobian use the identical shifted state.  All offsets
                    // are explicit and neutral by default.
                    if (useSentaurusBox && materialInterfaceNodes[
                            static_cast<std::size_t>(node)]) {
                        Real sideOffset =
                            config.sentaurusInterfaceInsulatorHalfJumpOffset;
                        if (data.isTransport) {
                            const Region& region = mesh.getRegion(cell.region_id);
                            sideOffset = isPolySilicon(region.material)
                                ? config.
                                    sentaurusInterfacePolysiliconHalfJumpOffset
                                : config.sentaurusInterfaceSiliconHalfJumpOffset;
                        }
                        w[local] += 2.0 * sideOffset;
                    }
                }
            }
            const Real x[3] = {p0.x * lengthScale_m,
                               p1.x * lengthScale_m,
                               p2.x * lengthScale_m};
            const Real y[3] = {p0.y * lengthScale_m,
                               p1.y * lengthScale_m,
                               p2.y * lengthScale_m};
            const Real b[3] = {y[1] - y[2], y[2] - y[0], y[0] - y[1]};
            const Real c[3] = {x[2] - x[1], x[0] - x[2], x[1] - x[0]};
            const Real fourArea2 = 4.0 * cellArea_m2 * cellArea_m2;
            const std::array<Real, 3> sentaurusBoxMeasures =
                useSentaurusBox
                ? averageBoxTriangleMeasures(x, y, cellArea_m2)
                : std::array<Real, 3>{};
            Real cvfemFlux[3][3]{};
            if (config.globalDiscretization == "cvfem_full") {
                const Real centroidX = (x[0] + x[1] + x[2]) / 3.0;
                const Real centroidY = (y[0] + y[1] + y[2]) / 3.0;
                constexpr int pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto& pair : pairs) {
                    const int a = pair[0];
                    const int neighbour = pair[1];
                    const Real midpointX = 0.5 * (x[a] + x[neighbour]);
                    const Real midpointY = 0.5 * (y[a] + y[neighbour]);
                    const Real tangentX = centroidX - midpointX;
                    const Real tangentY = centroidY - midpointY;
                    Real normalLengthX = tangentY;
                    Real normalLengthY = -tangentX;
                    const Real towardNeighbourX = x[neighbour] - x[a];
                    const Real towardNeighbourY = y[neighbour] - y[a];
                    if (normalLengthX * towardNeighbourX +
                            normalLengthY * towardNeighbourY < 0.0) {
                        normalLengthX = -normalLengthX;
                        normalLengthY = -normalLengthY;
                    }
                    for (int d = 0; d < 3; ++d) {
                        const Real coefficient =
                            (b[d] * normalLengthX + c[d] * normalLengthY) /
                            (2.0 * cellArea_m2);
                        cvfemFlux[a][d] += coefficient;
                        cvfemFlux[neighbour][d] -= coefficient;
                    }
                }
            }
            Real gradientSquared = 0.0;
            for (int d = 0; d < 3; ++d) {
                for (int e = 0; e < 3; ++e) {
                    gradientSquared += w[d] * w[e] *
                        (b[d] * b[e] + c[d] * c[e]) / fourArea2;
                }
            }
            Eq231CellResidualRecord record;
            if (diagnostic != nullptr) {
                record.cellId = data.cellId;
                record.regionId = cell.region_id;
                record.isTransport = data.isTransport;
                record.nodes = {nodes[0], nodes[1], nodes[2]};
                record.area_m2 = cellArea_m2;
                record.gradientSquaredPerM2 = gradientSquared;
                for (int local = 0; local < 3; ++local) {
                    record.lambda_V[local] = lambda[local];
                    record.auxiliary[local] = w[local];
                }
            }
            for (int a = 0; a < 3; ++a) {
                const int row = static_cast<int>(nodes[a]);
                if (!activeNodes[static_cast<std::size_t>(row)])
                    continue;
                const Real lumpedVolume = useSentaurusBox
                    ? sentaurusBoxMeasures[static_cast<std::size_t>(a)]
                    : cellArea_m2 / 3.0;
                // Eq. 231 is written for the dimensionless beta-weighted
                // driving potential w.  With C=gamma*hbar^2/(6*m*q), its
                // voltage form is div(grad(w)) + theta*|grad(w)|^2
                // + 2*Lambda/C = 0.  Lambda is already stored in volts, so
                // no additional thermal-voltage factor belongs here.
                const Real equationScale = 2.0 / cellCoefficientVm2;
                Real gradientContribution = 0.0;
                Real reactionWeight = 1.0;
                Real reactionWeightDerivative = 0.0;
                if (useGssFittedFlux &&
                    std::abs(config.theta - 0.5) < 1.0e-14) {
                    // GSS Eq. 9.125/9.127 written with u=sqrt(n/Nref)
                    // and w=2*log(u).  Rewriting the P1 Laplace row as
                    // edge differences gives
                    //   -2*K_ad*expm1((w_d-w_a)/2), d != a.
                    // Eq. 9.128 keeps the decaying exponential but replaces
                    // the growing direction by its second-order branch;
                    // expm1 also avoids cancellation near zero. Each cell
                    // evaluates its own w trace, so
                    // band/DOS steps do not force u to be continuous even
                    // though the shared state Lambda is continuous.
                    for (int d = 0; d < 3; ++d) {
                        if (d == a)
                            continue;
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        const Real halfJump = 0.5 * (w[d] - w[a]);
                        const Real contribution = -2.0 * stiffness *
                            gssFittedJump(halfJump);
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                        const Real ratio =
                            gssFittedJumpDerivative(halfJump);
                        const Real derivative =
                            stiffness * ratio / thermalVoltage_V;
                        addJacobian(
                            row, static_cast<int>(nodes[d]), derivative);
                        addJacobian(row, row, -derivative);
                        if (linearizeBox && !constrainedRows[row]) {
                            // w = (drive - Phi/q)/VT, so only the flux
                            // derivative reverses sign for a drive column.
                            auto& drive = data.isTransport ? transportDrive : insulatorDrive;
                            drive.emplace_back(row, static_cast<int>(nodes[d]), -derivative);
                            drive.emplace_back(row, row, derivative);
                        }
                        if (diagnostic != nullptr) {
                            diagnostic->fittedEdges.push_back({
                                data.cellId,
                                cell.region_id,
                                nodes[a],
                                nodes[d],
                                w[a],
                                w[d],
                                halfJump,
                                halfJump < 0.0,
                                stiffness,
                                contribution,
                                derivative,
                            });
                        }
                    }
                } else if (useConservativeSqrtFlux &&
                    std::abs(config.theta - 0.5) < 1.0e-14) {
                    // Exact conservative theta=1/2 weak form after
                    // u=exp(w/2): -2*K*u + (2/C)*Lambda*u*V = 0.
                    // A fixed positive scale is applied to each assembled
                    // row. It preserves the nonlinear root and keeps the
                    // imported oxide exponent range finite without adding a
                    // state-dependent Jacobian term.
                    const Real rowReference =
                        conservativeRowLogReference(row);
                    Real scaledU[3]{};
                    for (int d = 0; d < 3; ++d) {
                        scaledU[d] = safeExponential(
                            0.5 * w[d] - rowReference);
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        const Real contribution =
                            -2.0 * stiffness * scaledU[d];
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                        addJacobian(
                            row, static_cast<int>(nodes[d]),
                            stiffness * scaledU[d] / thermalVoltage_V);
                    }
                    reactionWeight = scaledU[a];
                    reactionWeightDerivative =
                        -0.5 * reactionWeight / thermalVoltage_V;
                } else if (config.globalDiscretization == "exponential_fitted" &&
                    std::abs(config.theta - 0.5) < 1.0e-14) {
                    // For theta=1/2, div(grad(w))+|grad(w)|^2/2 =
                    // 2*laplace(u)/u with u=exp(w/2). Normalize u per cell;
                    // the arbitrary Nref then cancels without overflow.
                    const Real maximumW = std::max({w[0], w[1], w[2]});
                    Real u[3]{};
                    for (int d = 0; d < 3; ++d)
                        u[d] = safeExponential(0.5 * (w[d] - maximumW));
                    const Real inverseUa = 1.0 / u[a];
                    for (int d = 0; d < 3; ++d) {
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        const Real contribution =
                            -2.0 * stiffness * u[d] * inverseUa;
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                    }
                    Real rowWeightedU = 0.0;
                    for (int d = 0; d < 3; ++d) {
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        rowWeightedU += stiffness * u[d] * inverseUa;
                    }
                    for (int d = 0; d < 3; ++d) {
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        Real derivative = stiffness * u[d] * inverseUa /
                            thermalVoltage_V;
                        if (d == a)
                            derivative -= rowWeightedU / thermalVoltage_V;
                        addJacobian(
                            row, static_cast<int>(nodes[d]), derivative);
                    }
                } else if (config.globalDiscretization == "cvfem_full") {
                    // Control-volume finite element integration over the
                    // median-dual subcell.  The Laplacian is evaluated as
                    // flux through the two centroid-to-edge-midpoint faces.
                    // The nonlinear field source is constant in a Tri3 and
                    // therefore receives exactly one third of the cell area.
                    for (int d = 0; d < 3; ++d) {
                        const int col = static_cast<int>(nodes[d]);
                        const Real contribution = cvfemFlux[a][d] * w[d];
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                        addJacobian(
                            row, col,
                            -cvfemFlux[a][d] / thermalVoltage_V);
                    }
                    gradientContribution =
                        config.theta * gradientSquared * lumpedVolume;
                } else {
                    for (int d = 0; d < 3; ++d) {
                        const int col = static_cast<int>(nodes[d]);
                        const Real stiffness = cellArea_m2 *
                            (b[a] * b[d] + c[a] * c[d]) / fourArea2;
                        const Real contribution = -stiffness * w[d];
                        residual(row) += contribution;
                        if (diagnostic != nullptr) {
                            record.stiffness[a] += contribution;
                            diagnostic->stiffness(row) += contribution;
                        }
                        addJacobian(row, col, stiffness / thermalVoltage_V);
                    }
                    gradientContribution =
                        config.theta * gradientSquared * lumpedVolume;
                }
                const Real reactionLambda = useSentaurusBox
                    ? state(row) + maximumMaterialBandDrive(row)
                    : lambda[a];
                Real interfaceReactionMultiplier = 1.0;
                Real interfaceReactionOffset_V = 0.0;
                if (useSentaurusBox && materialInterfaceNodes[
                        static_cast<std::size_t>(row)]) {
                    if (data.isTransport) {
                        const bool poly = isPolySilicon(cellRegion.material);
                        interfaceReactionMultiplier = poly
                            ? config.sentaurusInterfacePolysiliconReactionWeight
                            : config.sentaurusInterfaceSiliconReactionWeight;
                        interfaceReactionOffset_V = poly
                            ? config.sentaurusInterfacePolysiliconReactionOffset_V
                            : config.sentaurusInterfaceSiliconReactionOffset_V;
                    } else {
                        const bool polyOnly = interfaceHasPolysilicon[
                            static_cast<std::size_t>(row)] &&
                            !interfaceHasSilicon[static_cast<std::size_t>(row)];
                        interfaceReactionMultiplier = polyOnly
                            ? config.
                                sentaurusInterfaceInsulatorAtPolysiliconReactionWeight
                            : config.
                                sentaurusInterfaceInsulatorAtSiliconReactionWeight;
                        interfaceReactionOffset_V = polyOnly
                            ? config.
                                sentaurusInterfaceInsulatorAtPolysiliconReactionOffset_V
                            : config.
                                sentaurusInterfaceInsulatorAtSiliconReactionOffset_V;
                    }
                } else if (useSentaurusBox && !data.isTransport &&
                           insulatorReentrantCornerNodes[
                               static_cast<std::size_t>(row)]) {
                    interfaceReactionMultiplier = config.
                        sentaurusInsulatorReentrantCornerReactionWeight;
                }
                const Real reactionContribution = equationScale *
                    (reactionLambda * interfaceReactionMultiplier +
                     interfaceReactionOffset_V) * lumpedVolume * reactionWeight;
                residual(row) += gradientContribution + reactionContribution;
                if (linearizeBox && !constrainedRows[row])
                    linearization->materialDriveDiagonal(row) +=
                        equationScale * lumpedVolume * interfaceReactionMultiplier;
                if (diagnostic != nullptr) {
                    record.gradientSquared[a] = gradientContribution;
                    record.reaction[a] = reactionContribution;
                    diagnostic->gradientSquared(row) += gradientContribution;
                    diagnostic->reaction(row) += reactionContribution;
                }
                addJacobian(
                    row, row, equationScale * lumpedVolume *
                        (reactionWeight +
                         reactionLambda * reactionWeightDerivative) *
                        interfaceReactionMultiplier);
                const bool fittedThetaHalf =
                    (config.globalDiscretization == "exponential_fitted" ||
                     useGssFittedFlux || useConservativeSqrtFlux) &&
                    std::abs(config.theta - 0.5) < 1.0e-14;
                if (!fittedThetaHalf) {
                    for (int d = 0; d < 3; ++d) {
                        Real gradientDerivative = 0.0;
                        for (int e = 0; e < 3; ++e) {
                            gradientDerivative += -2.0 * w[e] *
                                (b[d] * b[e] + c[d] * c[e]) /
                                (fourArea2 * thermalVoltage_V);
                        }
                        addJacobian(
                            row, static_cast<int>(nodes[d]),
                            config.theta * gradientDerivative * lumpedVolume);
                    }
                }
            }
            if (diagnostic != nullptr)
                diagnostic->cells.push_back(std::move(record));
        }
        // DEVSIM/Garcia-Asenov oxide closure.  In the oxide integrated
        // equation Lambda*V - b_n,ox*S/x_n = 0.  Multiplication by the
        // expanded Eq. 231 row scale 2/b_n,ox gives the constant face source
        // below.  Its Jacobian is exactly zero; the Lambda reaction Jacobian
        // is already assembled in the oxide cells.
        for (const auto& segment : oxideInterfaceSegments) {
            const Edge& edge = mesh.getEdge(segment.edgeId);
            const Real endpointLength_m =
                0.5 * edge.length * lengthScale_m;
            const Real boundaryContribution =
                -2.0 * endpointLength_m / oxidePenetrationDepth_m;
            const int rows[2] = {
                static_cast<int>(edge.n0), static_cast<int>(edge.n1)};
            for (const int row : rows) {
                if (!activeNodes[static_cast<std::size_t>(row)] ||
                    dirichletLambda_V.contains(static_cast<Index>(row)))
                    continue;
                residual(row) += boundaryContribution;
                if (diagnostic != nullptr)
                    diagnostic->interfaceBoundary(row) += boundaryContribution;
            }
        }
        for (int row = 0; row < nodeCount; ++row) {
            if (activeNodes[static_cast<std::size_t>(row)])
                continue;
            residual(row) = state(row);
        }
        for (const auto& [node, lambda] : dirichletLambda_V) {
            const int row = static_cast<int>(node);
            residual(row) = state(row) - (stateIsLambda
                ? lambda : lambda - stateBandDrive_V(row));
            if (linearizeBox)
                linearization->materialDriveDiagonal(row) = 1.0;
        }
        // Constrained equations have a single unit diagonal. Keep their
        // column contributions in every free row.
        if (triplets != nullptr) {
            for (int row = 0; row < nodeCount; ++row) {
                if (constrainedRows[static_cast<std::size_t>(row)])
                    triplets->emplace_back(row, row, 1.0);
            }
        }
    };
    if (linearization != nullptr) {
        std::vector<Eigen::Triplet<Real>> triplets;
        if (!linearizeBox) {
            assemble(continuousPotential, linearization->residual, nullptr, nullptr);
            return result;
        }
        assemble(continuousPotential, linearization->residual, &triplets, nullptr);
        linearization->potentialLike.resize(nodeCount, nodeCount);
        linearization->potentialLike.setFromTriplets(triplets.begin(), triplets.end());
        linearization->transportDrive.resize(nodeCount, nodeCount);
        linearization->transportDrive.setFromTriplets(transportDrive.begin(), transportDrive.end());
        linearization->insulatorDrive.resize(nodeCount, nodeCount);
        linearization->insulatorDrive.setFromTriplets(insulatorDrive.begin(), insulatorDrive.end());
        return result;
    }
    if (!config.residualDiagnosticPrefix.empty()) {
        VectorXd diagnosticResidual;
        Eq231ResidualDiagnostic diagnostic;
        assemble(continuousPotential, diagnosticResidual, nullptr, &diagnostic);
        writeEq231ResidualDiagnostic(
            config.residualDiagnosticPrefix, mesh, diagnostic,
            continuousPotential, stateBandDrive_V, activeNodes,
            dirichletLambda_V, stateIsLambda);
    }
    VectorXd residual(nodeCount);
    std::vector<Eigen::Triplet<Real>> triplets;
    for (int iteration = 1; iteration <= config.maxIterations; ++iteration) {
        assemble(continuousPotential, residual, &triplets, nullptr);
        result.residualInfinityNorm = residual.lpNorm<Eigen::Infinity>();
        if (initialResidualNorm < 0.0)
            initialResidualNorm = result.residualInfinityNorm;
        const Real residualTolerance = config.absoluteTolerance_V +
            config.relativeTolerance * std::max<Real>(initialResidualNorm, 1.0);
        if (result.residualInfinityNorm <= residualTolerance) {
            result.iterations = iteration;
            result.converged = true;
            break;
        }
        SparseMatrixd jacobian(nodeCount, nodeCount);
        jacobian.setFromTriplets(triplets.begin(), triplets.end());
        jacobian.prune(0.0);
        VectorXd update = solver.solve(jacobian, -residual);
        Real updateNorm = update.lpNorm<Eigen::Infinity>();
        Eigen::Index maxUpdateIndex = 0;
        update.cwiseAbs().maxCoeff(&maxUpdateIndex);
        if (config.maxUpdate_V > 0.0 && updateNorm > config.maxUpdate_V) {
            update *= config.maxUpdate_V / updateNorm;
            updateNorm = config.maxUpdate_V;
        }
        Real lineScale = config.damping;
        VectorXd trialResidual;
        VectorXd trial = continuousPotential + lineScale * update;
        assemble(trial, trialResidual, nullptr, nullptr);
        while (trialResidual.lpNorm<Eigen::Infinity>() >=
                   result.residualInfinityNorm && lineScale > 1.0e-6) {
            lineScale *= 0.5;
            trial = continuousPotential + lineScale * update;
            assemble(trial, trialResidual, nullptr, nullptr);
        }
        const Real trialResidualNorm =
            trialResidual.lpNorm<Eigen::Infinity>();
        // A vanishing line-search scale is a stagnation/failure condition,
        // not Newton convergence.  Previously the scaled update alone could
        // accept a state whose residual was still many orders of magnitude
        // larger than the incoming residual.
        if (!(trialResidualNorm < result.residualInfinityNorm)) {
            result.lastUpdateInfinityNorm_V = 0.0;
            result.maxUpdateNode = static_cast<Index>(maxUpdateIndex);
            result.maxUpdateNodeValue_V = 0.0;
            break;
        }
        continuousPotential = std::move(trial);
        for (const auto& [node, lambda] : dirichletLambda_V) {
            continuousPotential(static_cast<int>(node)) = stateIsLambda
                ? lambda
                : lambda - stateBandDrive_V(static_cast<int>(node));
        }
        result.iterations = iteration;
        result.lastUpdateInfinityNorm_V = lineScale * updateNorm;
        result.maxUpdateNode = static_cast<Index>(maxUpdateIndex);
        result.maxUpdateNodeValue_V = lineScale * update(maxUpdateIndex);
        result.residualInfinityNorm = trialResidualNorm;
        const Real scale = std::max<Real>(continuousPotential.lpNorm<Eigen::Infinity>(), 1.0);
        if (lineScale * updateNorm <= config.absoluteTolerance_V +
                config.relativeTolerance * scale &&
            trialResidualNorm <= residualTolerance) {
            result.converged = true;
            break;
        }
    }
    if (stateIsLambda) {
        result.potentialLike_V.resize(0);
        result.potential_V = continuousPotential;
    } else {
        result.potentialLike_V = continuousPotential;
        result.potential_V = continuousPotential + stateBandDrive_V;
    }
    result.potentialInfinityNorm_V =
        result.potential_V.lpNorm<Eigen::Infinity>();
    return result;
}

DensityGradientQuantumPotentialResult solveElectronDensityGradientPotentialLikeGlobal(
    const DeviceMesh& mesh,
    const std::vector<DensityGradientCellMaterial>& cellMaterials,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    DensityGradientQuantumPotentialConfig config,
    const VectorXd& initialLambda_V,
    const VectorXd& initialPotentialLike_V)
{
    return solveElectronDensityGradientPotentialLikeGlobalImpl(
        mesh, cellMaterials, nodeOutputBandDrive_V, activeNodes, dirichletLambda_V,
        thermalVoltage_V, units, std::move(config), initialLambda_V,
        initialPotentialLike_V, nullptr, false);
}

DensityGradientBoxLinearization linearizeElectronDensityGradientBox(
    const DeviceMesh& mesh,
    const std::vector<DensityGradientCellMaterial>& cellMaterials,
    const VectorXd& nodeOutputBandDrive_V,
    const std::vector<bool>& activeNodes,
    const std::unordered_map<Index, Real>& dirichletLambda_V,
    Real thermalVoltage_V,
    PhysicalUnitSystem units,
    const DensityGradientQuantumPotentialConfig& config,
    const VectorXd& potentialLike_V,
    bool assembleJacobian)
{
    DensityGradientBoxLinearization linearization;
    solveElectronDensityGradientPotentialLikeGlobalImpl(
        mesh, cellMaterials, nodeOutputBandDrive_V, activeNodes, dirichletLambda_V,
        thermalVoltage_V, units, config, {}, potentialLike_V, &linearization,
        assembleJacobian);
    return linearization;
}

} // namespace vela
