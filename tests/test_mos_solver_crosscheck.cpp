#include <catch2/catch_test_macros.hpp>

#include "vela/io/MeshReader.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/post/ContactCurrent.h"
#include "vela/solver/GummelSolver.h"
#include "vela/solver/NewtonSolver.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

using namespace vela;

namespace {

struct MosCase {
    std::string exampleName;
    Real gateBias;
    Real drainBias;
};

nlohmann::json readJson(const std::filesystem::path& path)
{
    std::ifstream input(path);
    REQUIRE(input.is_open());
    nlohmann::json json;
    input >> json;
    return json;
}

DopingModel dopingFromDeck(const DeviceMesh& mesh, const nlohmann::json& cfg)
{
    std::vector<RegionDopingSpec> specs;
    for (const auto& entry : cfg.at("doping")) {
        specs.push_back({
            entry.at("region").get<std::string>(),
            entry.at("donors").get<Real>(),
            entry.at("acceptors").get<Real>(),
        });
    }
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

std::unordered_map<std::string, Real> biasesFromDeck(const nlohmann::json& cfg,
                                                     Real gateBias,
                                                     Real drainBias)
{
    std::unordered_map<std::string, Real> biases;
    for (const auto& contact : cfg.at("contacts"))
        biases[contact.at("name").get<std::string>()] = contact.at("bias").get<Real>();
    biases["gate"] = gateBias;
    biases["drain"] = drainBias;
    return biases;
}

void requireFinitePhysicalSolution(const DDSolution& sol,
                                   Index numNodes,
                                   const std::string& label,
                                   bool requireConverged)
{
    INFO(label);
    if (requireConverged)
        REQUIRE(sol.converged);
    REQUIRE(sol.psi.size() == static_cast<int>(numNodes));
    REQUIRE(sol.phin.size() == static_cast<int>(numNodes));
    REQUIRE(sol.phip.size() == static_cast<int>(numNodes));
    REQUIRE(sol.n.size() == static_cast<int>(numNodes));
    REQUIRE(sol.p.size() == static_cast<int>(numNodes));

    for (int i = 0; i < static_cast<int>(numNodes); ++i) {
        REQUIRE(std::isfinite(sol.psi(i)));
        REQUIRE(std::isfinite(sol.phin(i)));
        REQUIRE(std::isfinite(sol.phip(i)));
        REQUIRE(std::isfinite(sol.n(i)));
        REQUIRE(std::isfinite(sol.p(i)));
        REQUIRE(sol.n(i) >= 0.0);
        REQUIRE(sol.p(i) >= 0.0);
    }
}

void checkMosCase(const MosCase& mos)
{
    const std::filesystem::path exampleDir =
        std::filesystem::path(VELA_SOURCE_DIR) / "examples" / mos.exampleName;
    const nlohmann::json cfg = readJson(exampleDir / "simulation_iv.json");

    JsonMeshReader reader;
    DeviceMesh mesh = reader.read((exampleDir / "mesh.json").string());
    MaterialDatabase matdb;
    DopingModel doping = dopingFromDeck(mesh, cfg);
    const auto biases = biasesFromDeck(cfg, mos.gateBias, mos.drainBias);
    const Real expectedSign =
        cfg.at("regression").at("mos").value("drain_current_sign", 1.0);

    GummelConfig gummelCfg = gummelConfigFromJson(cfg.at("solver"));
    gummelCfg.maxIter = 100;
    gummelCfg.reltol = 1.0e-5;
    DDSolution gummel = runGummel(mesh, matdb, doping, biases, gummelCfg);
    requireFinitePhysicalSolution(gummel,
                                  mesh.numNodes(),
                                  mos.exampleName + " Gummel",
                                  true);

    NewtonConfig newtonCfg;
    newtonCfg.maxIter = 25;
    // MOS DD Newton coverage is a smoke/stability check: require an accepted
    // residual-reducing Newton update without imposing PN-grade tolerances.
    newtonCfg.reltol = 9.9e-1;
    newtonCfg.abstol = 1.0e-18;
    newtonCfg.dampingFactor = 1.0;
    newtonCfg.lineSearch = true;
    newtonCfg.verbose = false;
    newtonCfg.warmStart = true;
    NewtonResult newton = runNewton(mesh, matdb, doping, biases, gummel, newtonCfg);
    INFO(mos.exampleName << " Newton initial=" << newton.initialResidualNorm
                         << " final=" << newton.finalResidualNorm
                         << " iters=" << newton.iters);
    REQUIRE(newton.iters > 0);
    REQUIRE(newton.finalResidualNorm < newton.initialResidualNorm);
    requireFinitePhysicalSolution(newton.solution,
                                  mesh.numNodes(),
                                  mos.exampleName + " Newton",
                                  false);


    // 保持legacy路径，scaling默认关闭
    ContactCurrent current(mesh, matdb, doping, {}, 300.0, {});
    const Real gummelDrainCurrent = current.compute(gummel, "drain").totalCurrent;
    const Real newtonDrainCurrent = current.compute(newton.solution, "drain").totalCurrent;

    REQUIRE(std::isfinite(gummelDrainCurrent));
    REQUIRE(std::isfinite(newtonDrainCurrent));
    REQUIRE(expectedSign * gummelDrainCurrent > 0.0);
    REQUIRE(expectedSign * newtonDrainCurrent > 0.0);

    const Real reference = std::max({std::abs(gummelDrainCurrent), 1.0e-30});
    const Real relativeDifference =
        std::abs(newtonDrainCurrent - gummelDrainCurrent) / reference;
    REQUIRE(relativeDifference < 1.0e2);
}

} // namespace

TEST_CASE("NMOS and PMOS single-bias Newton/Gummel cross-check", "[mos_solver_crosscheck][newton][gummel][mos]")
{
    checkMosCase({"nmos2d_dd", 0.01, 0.01});
    checkMosCase({"pmos2d_dd", -0.01, -0.01});
}
