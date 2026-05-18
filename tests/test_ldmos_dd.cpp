#include <catch2/catch_test_macros.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/DDAssembler.h"
#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/GummelSolver.h"

#include <cmath>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

using namespace vela;

namespace {

DeviceMesh loadLdmosMesh()
{
    const std::filesystem::path meshPath =
        std::filesystem::path(VELA_SOURCE_DIR) / "examples" / "ldmos2d" / "mesh.json";
    return JsonMeshReader().read(meshPath.string());
}

DopingModel makeLdmosDoping(const DeviceMesh& mesh)
{
    return DopingModel::fromMeshAndRegions(mesh, {
        {"p_body_drift", 0.0, 1.0e21},
        {"n_source", 1.0e22, 0.0},
        {"n_drain_drift", 2.0e21, 0.0},
        {"field_oxide", 0.0, 0.0},
    });
}

void requireFiniteVector(const VectorXd& values)
{
    for (Eigen::Index i = 0; i < values.size(); ++i)
        REQUIRE(std::isfinite(values(i)));
}

void requireFiniteSparseMatrix(const SparseMatrixd& matrix)
{
    REQUIRE(matrix.rows() > 0);
    REQUIRE(matrix.cols() > 0);
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrixd::InnerIterator it(matrix, outer); it; ++it)
            REQUIRE(std::isfinite(it.value()));
    }
}

void requireFiniteScalarSystem(const DDAssembler& assembler, Index expectedNodes)
{
    REQUIRE(assembler.matrix().rows() == static_cast<int>(expectedNodes));
    REQUIRE(assembler.matrix().cols() == static_cast<int>(expectedNodes));
    REQUIRE(assembler.rhs().size() == static_cast<int>(expectedNodes));
    requireFiniteSparseMatrix(assembler.matrix());
    requireFiniteVector(assembler.rhs());
}

bool allFinite(const VectorXd& values)
{
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values(i)))
            return false;
    }
    return true;
}

} // namespace

TEST_CASE("LDMOS mixed-material drift-diffusion assembly stays finite", "[ldmos][dd]")
{
    const DeviceMesh mesh = loadLdmosMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = makeLdmosDoping(mesh);

    const int nNodes = static_cast<int>(mesh.numNodes());
    const VectorXd psi = VectorXd::Zero(nNodes);
    VectorXd n = VectorXd::Zero(nNodes);
    VectorXd p = VectorXd::Zero(nNodes);
    for (Index node = 0; node < mesh.numNodes(); ++node) {
        n(static_cast<int>(node)) = std::max(1.0e10, doping.donors(node));
        p(static_cast<int>(node)) = std::max(1.0e10, doping.acceptors(node));
    }

    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};
    DDAssembler assembler(mesh,
                          matdb,
                          doping,
                          constants::Vt_300,
                          MobilityModelConfig{},
                          recombination,
                          BandgapNarrowingConfig{},
                          ImpactIonizationModelConfig{},
                          {RegionFixedChargeSpec{"field_oxide", 5.0e20}},
                          {});

    REQUIRE_NOTHROW(assembler.assemblePoissonWithCarriers(n, p, psi));
    requireFiniteScalarSystem(assembler, mesh.numNodes());

    REQUIRE_NOTHROW(assembler.assembleElectronContinuity(psi, n, p));
    requireFiniteScalarSystem(assembler, mesh.numNodes());

    REQUIRE_NOTHROW(assembler.assembleHoleContinuity(psi, n, p));
    requireFiniteScalarSystem(assembler, mesh.numNodes());
}

TEST_CASE("LDMOS low-bias multi-terminal Gummel solve converges", "[ldmos][gummel]")
{
    const DeviceMesh mesh = loadLdmosMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = makeLdmosDoping(mesh);

    const std::unordered_map<std::string, Real> biases = {
        {"body", 0.0},
        {"source", 0.0},
        {"gate", 0.1},
        {"drain", 0.05},
    };

    GummelConfig cfg;
    cfg.maxIter = 120;
    cfg.reltol = 1.0e-5;
    cfg.abstol = 1.0e18;
    cfg.dampingPsi = 0.5;

    const DDSolution solution = runGummel(
        mesh,
        matdb,
        doping,
        biases,
        ContactSpecsMap{},
        cfg,
        {RegionFixedChargeSpec{"field_oxide", 5.0e20}},
        {});

    REQUIRE(solution.converged);
    REQUIRE(solution.iters > 0);
    REQUIRE(allFinite(solution.psi));
    REQUIRE(allFinite(solution.phin));
    REQUIRE(allFinite(solution.phip));
    REQUIRE(allFinite(solution.n));
    REQUIRE(allFinite(solution.p));
}
