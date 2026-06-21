#include "vela/core/PhysicalConstants.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/BandgapNarrowing.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/ImpactIonizationModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/RecombinationModel.h"

#include <Eigen/Dense>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace vela;

namespace {

struct Args {
    std::filesystem::path output;
    std::vector<Real> biases;
};

DeviceMesh makeAuditMesh(Real L)
{
    DeviceMesh mesh;

    Node n0; n0.id = 0; n0.x = 0.0;     n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;       n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = L;       n2.y = L;       mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0;     n3.y = L;       mesh.addNode(n3);
    Node n4; n4.id = 4; n4.x = 0.5 * L; n4.y = 0.5 * L; mesh.addNode(n4);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 4}; mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 0; c1.node_ids = {1, 2, 4}; mesh.addCell(c1);
    Cell c2; c2.id = 2; c2.type = CellType::Tri3; c2.region_id = 1; c2.node_ids = {2, 3, 4}; mesh.addCell(c2);
    Cell c3; c3.id = 3; c3.type = CellType::Tri3; c3.region_id = 1; c3.node_ids = {3, 0, 4}; mesh.addCell(c3);

    Region r0; r0.id = 0; r0.name = "n_region"; r0.material = "Si"; r0.cell_ids = {0, 1}; mesh.addRegion(r0);
    Region r1; r1.id = 1; r1.name = "p_region"; r1.material = "Si"; r1.cell_ids = {2, 3}; mesh.addRegion(r1);

    Contact cathode; cathode.id = 0; cathode.name = "cathode"; cathode.region_id = 0; cathode.node_ids = {1, 2}; mesh.addContact(cathode);
    Contact anode; anode.id = 1; anode.name = "anode"; anode.region_id = 1; anode.node_ids = {0, 3}; mesh.addContact(anode);

    mesh.buildEdges();
    return mesh;
}

DopingModel makeAuditDoping(const DeviceMesh& mesh)
{
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

CoupledDDState makeAuditState(int nodeCount, Real bias)
{
    CoupledDDState state;
    state.psi.resize(nodeCount);
    state.phin.resize(nodeCount);
    state.phip.resize(nodeCount);

    const Real center = 0.01 * bias;
    state.psi = Eigen::VectorXd::LinSpaced(nodeCount, center - 0.08, center + 0.08);
    state.phin = Eigen::VectorXd::LinSpaced(nodeCount, center + 0.32, center - 0.28);
    state.phip = Eigen::VectorXd::LinSpaced(nodeCount, center - 0.27, center + 0.31);
    return state;
}

ImpactIonizationModelConfig sgAvalancheConfig()
{
    ImpactIonizationModelConfig impact;
    impact.model = "van_overstraeten";
    impact.drivingForce = "quasi_fermi_gradient";
    impact.generation = "current_density";
    impact.currentApproximation = "density_gradient";
    impact.drivingForceInterpolation = "quasi_fermi_to_electric_field";
    impact.electronDrivingForceRefDensity = 1.0e20;
    impact.holeDrivingForceRefDensity = 1.0e20;
    return impact;
}

CoupledDDBoundaryConditions auditBoundaryConditions(const CoupledDDState& state)
{
    CoupledDDBoundaryConditions bcs;
    for (Index node : {Index{0}, Index{2}}) {
        const int ii = static_cast<int>(node);
        bcs.psi[node] = state.psi(ii);
        bcs.phin[node] = state.phin(ii);
        bcs.phip[node] = state.phip(ii);
    }
    return bcs;
}

Eigen::MatrixXd denseJacobian(const CoupledDDAssembler& assembler,
                              const Eigen::VectorXd& x,
                              const CoupledDDBoundaryConditions& bcs,
                              bool finiteDifference)
{
    const SparseMatrixd jacobian = finiteDifference
        ? assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7)
        : assembler.assembleJacobian(x, bcs);
    return Eigen::MatrixXd(jacobian);
}

struct MatrixPair {
    Eigen::MatrixXd analytic;
    Eigen::MatrixXd fd;
};

MatrixPair matrixPair(const DeviceMesh& mesh,
                      const MaterialDatabase& matdb,
                      const DopingModel& doping,
                      const CoupledDDState& state,
                      const MobilityModelConfig& mobility,
                      const RecombinationModelConfig& recombination,
                      const ImpactIonizationModelConfig& impact,
                      const CoupledDDBoundaryConditions& bcs = {})
{
    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        mobility,
        recombination,
        BandgapNarrowingConfig{},
        impact);
    const Eigen::VectorXd x = assembler.pack(state);
    return {
        denseJacobian(assembler, x, bcs, false),
        denseJacobian(assembler, x, bcs, true),
    };
}

std::vector<int> rowsForBlock(const std::string& block, int nodeCount)
{
    std::vector<int> rows;
    if (block == "poisson") {
        for (int i = 0; i < nodeCount; ++i)
            rows.push_back(i);
    } else {
        for (int i = nodeCount; i < 3 * nodeCount; ++i)
            rows.push_back(i);
    }
    return rows;
}

Real restrictedNorm(const Eigen::MatrixXd& matrix, const std::vector<int>& rows)
{
    Real sum = 0.0;
    for (int row : rows) {
        for (int col = 0; col < matrix.cols(); ++col) {
            const Real value = matrix(row, col);
            sum += value * value;
        }
    }
    return std::sqrt(sum);
}

void writeRow(std::ofstream& out,
              Real bias,
              const std::string& block,
              const Eigen::MatrixXd& analytic,
              const Eigen::MatrixXd& fd,
              const std::vector<int>& rows)
{
    const Real analyticNorm = restrictedNorm(analytic, rows);
    const Real fdNorm = restrictedNorm(fd, rows);
    const Real diffNorm = restrictedNorm(analytic - fd, rows);
    const Real relDiff = diffNorm / std::max<Real>(1.0, fdNorm);
    out << bias << ',' << block << ','
        << analyticNorm << ',' << fdNorm << ',' << diffNorm << ',' << relDiff << '\n';
}

Args parseArgs(int argc, char** argv)
{
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc) {
            args.output = argv[++i];
        } else if (arg == "--bias" && i + 1 < argc) {
            args.biases.push_back(std::stod(argv[++i]));
        } else {
            throw std::invalid_argument("unknown or incomplete argument: " + arg);
        }
    }
    if (args.output.empty())
        throw std::invalid_argument("--output is required");
    if (args.biases.empty())
        throw std::invalid_argument("at least one --bias is required");
    return args;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Args args = parseArgs(argc, argv);
        std::filesystem::create_directories(args.output.parent_path());

        const DeviceMesh strongFieldMesh = makeAuditMesh(1.0e-8);
        const DeviceMesh mildFieldMesh = makeAuditMesh(1.0e-6);
        MaterialDatabase matdb;
        const DopingModel strongFieldDoping = makeAuditDoping(strongFieldMesh);
        const DopingModel mildFieldDoping = makeAuditDoping(mildFieldMesh);
        const int nodeCount = static_cast<int>(strongFieldMesh.numNodes());
        const MobilityModelConfig mobility = mobilityModelConfig("constant");
        const RecombinationModelConfig noRecomb = recombinationModelConfig({"none"});
        const RecombinationModelConfig recomb = recombinationModelConfig({"srh", "auger"});
        const ImpactIonizationModelConfig noImpact{};
        const ImpactIonizationModelConfig impact = sgAvalancheConfig();

        std::ofstream out(args.output);
        if (!out)
            throw std::runtime_error("failed to open output CSV");
        out << "bias_V,block,analytic_norm,fd_norm,diff_norm,rel_diff\n";

        for (const Real bias : args.biases) {
            const CoupledDDState strongState = makeAuditState(nodeCount, bias);
            const CoupledDDState mildState = makeAuditState(nodeCount, bias);
            const MatrixPair base =
                matrixPair(strongFieldMesh, matdb, strongFieldDoping, strongState, mobility, noRecomb, noImpact);
            const MatrixPair withImpact =
                matrixPair(strongFieldMesh, matdb, strongFieldDoping, strongState, mobility, noRecomb, impact);
            const MatrixPair mildBase =
                matrixPair(mildFieldMesh, matdb, mildFieldDoping, mildState, mobility, noRecomb, noImpact);
            const MatrixPair withRecomb =
                matrixPair(mildFieldMesh, matdb, mildFieldDoping, mildState, mobility, recomb, noImpact);
            const CoupledDDBoundaryConditions bcs = auditBoundaryConditions(strongState);
            const MatrixPair constrained =
                matrixPair(strongFieldMesh, matdb, strongFieldDoping, strongState, mobility, recomb, impact, bcs);

            writeRow(out, bias, "poisson", base.analytic, base.fd, rowsForBlock("poisson", nodeCount));
            writeRow(out, bias, "transport", base.analytic, base.fd, rowsForBlock("transport", nodeCount));
            writeRow(
                out,
                bias,
                "srh_auger",
                withRecomb.analytic - mildBase.analytic,
                withRecomb.fd - mildBase.fd,
                rowsForBlock("srh_auger", nodeCount));
            writeRow(
                out,
                bias,
                "sg_avalanche",
                withImpact.analytic - base.analytic,
                withImpact.fd - base.fd,
                rowsForBlock("sg_avalanche", nodeCount));
            writeRow(
                out,
                bias,
                "dirichlet_or_gauge",
                constrained.analytic,
                constrained.fd,
                std::vector<int>{0, 2, nodeCount, nodeCount + 2, 2 * nodeCount, 2 * nodeCount + 2});
        }
    } catch (const std::exception& ex) {
        std::cerr << "pn2d_jacobian_block_audit: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
