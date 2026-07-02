#include <catch2/catch_test_macros.hpp>

#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/solver/NewtonSolver.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using vela::Cell;
using vela::CellType;
using vela::Contact;
using vela::DeviceMesh;
using vela::DopingModel;
using vela::Index;
using vela::MaterialDatabase;
using vela::NewtonConfig;
using vela::NewtonResult;
using vela::NewtonSolver;
using vela::Node;
using vela::Real;
using vela::Region;
using vela::RegionDopingSpec;

DeviceMesh makePNMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

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

DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"n_region", 1.0e21, 0.0},
        {"p_region", 0.0, 1.0e21},
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

NewtonConfig quietNewtonConfig()
{
    NewtonConfig cfg;
    cfg.maxIter = 12;
    cfg.reltol = 1.0e-8;
    cfg.abstol = 1.0e-18;
    cfg.dampingFactor = 1.0;
    cfg.lineSearch = true;
    cfg.warmStart = true;
    cfg.verbose = false;
    return cfg;
}

std::filesystem::path makeTempDir(const std::string& name)
{
    const auto base = std::filesystem::temp_directory_path() /
        ("vela_external_state_" + name + "_" + std::to_string(std::rand()));
    std::filesystem::create_directories(base);
    return base;
}

void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    REQUIRE(out.is_open());
    out << text;
}

void writeMeshJson(const std::filesystem::path& path)
{
    writeText(path, R"json({
  "nodes": [
    {"id": 0, "x": 0.0, "y": 0.0},
    {"id": 1, "x": 1.0e-6, "y": 0.0},
    {"id": 2, "x": 1.0e-6, "y": 1.0e-6},
    {"id": 3, "x": 0.0, "y": 1.0e-6},
    {"id": 4, "x": 5.0e-7, "y": 5.0e-7}
  ],
  "triangles": [
    {"id": 0, "region_id": 0, "node_ids": [0, 1, 4]},
    {"id": 1, "region_id": 0, "node_ids": [1, 2, 4]},
    {"id": 2, "region_id": 1, "node_ids": [2, 3, 4]},
    {"id": 3, "region_id": 1, "node_ids": [3, 0, 4]}
  ],
  "regions": [
    {"id": 0, "name": "n_region", "material": "Si", "cell_ids": [0, 1]},
    {"id": 1, "name": "p_region", "material": "Si", "cell_ids": [2, 3]}
  ],
  "contacts": [
    {"id": 0, "name": "cathode", "region_id": 0, "node_ids": [1, 2]},
    {"id": 1, "name": "anode", "region_id": 1, "node_ids": [0, 3]}
  ]
})json");
}

void writeScalarField(const std::filesystem::path& path,
                      const std::vector<Real>& values,
                      const std::string& header = "node_id,component0\n",
                      const std::vector<Index>& nodeOrder = {0, 1, 2, 3, 4})
{
    std::ofstream out(path);
    REQUIRE(out.is_open());
    out << header;
    for (const Index node : nodeOrder)
        out << node << ',' << values.at(static_cast<std::size_t>(node)) << '\n';
}

struct RunnerCase {
    std::filesystem::path dir;
    std::filesystem::path fieldsDir;
    std::filesystem::path configPath;
    NewtonResult equilibrium;
};

RunnerCase makeRunnerCase(const std::string& name)
{
    RunnerCase c;
    c.dir = makeTempDir(name);
    c.fieldsDir = c.dir / "fields";
    std::filesystem::create_directories(c.fieldsDir);
    writeMeshJson(c.dir / "mesh.json");

    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    NewtonConfig cfg = quietNewtonConfig();
    c.equilibrium = NewtonSolver(
        mesh,
        matdb,
        doping,
        {{"anode", 0.0}, {"cathode", 0.0}},
        cfg).solve();
    REQUIRE(c.equilibrium.converged);

    std::vector<Real> psi(5);
    std::vector<Real> phin(5);
    std::vector<Real> phip(5);
    for (int i = 0; i < 5; ++i) {
        psi[static_cast<std::size_t>(i)] = c.equilibrium.solution.psi(i);
        phin[static_cast<std::size_t>(i)] = c.equilibrium.solution.phin(i);
        phip[static_cast<std::size_t>(i)] = c.equilibrium.solution.phip(i);
    }
    writeScalarField(c.fieldsDir / "ElectrostaticPotential_region0.csv", psi);
    writeScalarField(c.fieldsDir / "eQuasiFermiPotential_region0.csv", phin);
    writeScalarField(c.fieldsDir / "hQuasiFermiPotential_region0.csv", phip);

    c.configPath = c.dir / "config.json";
    const nlohmann::json config = {
        {"simulation_type", "newton_solve_from_state"},
        {"mesh_file", "mesh.json"},
        {"state_fields_dir", "fields"},
        {"output_state_file", "out/state.csv"},
        {"contacts", nlohmann::json::array({
            {{"name", "anode"}, {"bias", 0.0}},
            {{"name", "cathode"}, {"bias", 0.0}},
        })},
        {"doping", nlohmann::json::array({
            {{"region", "n_region"}, {"donors", 1.0e21}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e21}},
        })},
        {"solver", {
            {"method", "newton"},
            {"max_iter", 0},
            {"abstol", 1.0e-6},
            {"reltol", 1.0e-8},
            {"line_search", true},
            {"warm_start", true},
            {"verbose", false},
        }},
    };
    writeText(c.configPath, config.dump(2));
    return c;
}

struct RunnerOutput {
    int exitCode = 0;
    std::string out;
    std::string err;
};

std::string readFile(const std::filesystem::path& path)
{
    std::ifstream in(path);
    REQUIRE(in.is_open());
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

RunnerOutput runCase(const std::filesystem::path& config)
{
    const std::filesystem::path dir = config.parent_path();
    const std::filesystem::path out = dir / "runner.stdout";
    const std::filesystem::path err = dir / "runner.stderr";
    const std::string command = std::string("cd /D \"") + dir.string() +
        "\" && \"" + VELA_EXAMPLE_RUNNER_EXE + "\" --config \"" +
        config.filename().string() + "\" > \"" + out.filename().string() +
        "\" 2> \"" + err.filename().string() + "\"";
    RunnerOutput result;
    result.exitCode = std::system(command.c_str());
    result.out = readFile(out);
    result.err = readFile(err);
    return result;
}

} // namespace

TEST_CASE("newton_solve_from_state accepts a converged external field state",
          "[external_state][newton]")
{
    RunnerCase c = makeRunnerCase("success");

    const RunnerOutput run = runCase(c.configPath);

    REQUIRE(run.exitCode == 0);
    const auto status = nlohmann::json::parse(run.out);
    REQUIRE(status.at("simulation_type") == "newton_solve_from_state");
    REQUIRE(status.at("converged").get<bool>());
    REQUIRE(status.at("iterations").get<int>() == 0);
    REQUIRE(status.at("initial_residual").get<double>() <= 1.0e-6);
    REQUIRE(std::filesystem::exists(c.dir / "out" / "state.csv"));
}

TEST_CASE("newton_solve_from_state returns nonzero when Newton does not converge",
          "[external_state][newton]")
{
    RunnerCase c = makeRunnerCase("nonconverged");
    writeScalarField(
        c.fieldsDir / "ElectrostaticPotential_region0.csv",
        {0.0, 0.0, 0.0, 0.0, 10.0});
    nlohmann::json config = nlohmann::json::parse(readFile(c.configPath));
    config["solver"]["abstol"] = 1.0e-30;
    config["solver"]["max_iter"] = 0;
    writeText(c.configPath, config.dump(2));

    const RunnerOutput run = runCase(c.configPath);

    REQUIRE(run.exitCode != 0);
    const auto status = nlohmann::json::parse(run.out);
    REQUIRE_FALSE(status.at("converged").get<bool>());
}

TEST_CASE("external state field reader rejects malformed scalar CSVs",
          "[external_state][csv]")
{
    struct MalformedCase {
        std::string name;
        std::string file;
        std::string expected;
        void (*mutate)(RunnerCase&);
    };

    const std::vector<MalformedCase> cases = {
        {
            "missing_e_quasi_fermi_file",
            "eQuasiFermiPotential_region0.csv",
            "Cannot open scalar field CSV",
            [](RunnerCase& c) {
                std::filesystem::remove(c.fieldsDir / "eQuasiFermiPotential_region0.csv");
            },
        },
        {
            "missing_node_id_column",
            "ElectrostaticPotential_region0.csv",
            "node_id",
            [](RunnerCase& c) {
                writeScalarField(
                    c.fieldsDir / "ElectrostaticPotential_region0.csv",
                    {0, 0, 0, 0, 0},
                    "id,component0\n");
            },
        },
        {
            "missing_component0_column",
            "ElectrostaticPotential_region0.csv",
            "component0",
            [](RunnerCase& c) {
                writeScalarField(
                    c.fieldsDir / "ElectrostaticPotential_region0.csv",
                    {0, 0, 0, 0, 0},
                    "node_id,value\n");
            },
        },
        {
            "duplicate_node_id",
            "ElectrostaticPotential_region0.csv",
            "duplicate node_id",
            [](RunnerCase& c) {
                writeScalarField(
                    c.fieldsDir / "ElectrostaticPotential_region0.csv",
                    {0, 0, 0, 0, 0},
                    "node_id,component0\n",
                    {0, 1, 2, 3, 4, 0});
            },
        },
        {
            "out_of_range_node_id",
            "ElectrostaticPotential_region0.csv",
            "node_id out of range",
            [](RunnerCase& c) {
                writeScalarField(
                    c.fieldsDir / "ElectrostaticPotential_region0.csv",
                    {0, 0, 0, 0, 0, 0},
                    "node_id,component0\n",
                    {0, 1, 2, 3, 5});
            },
        },
        {
            "missing_node_row",
            "ElectrostaticPotential_region0.csv",
            "missing a node row",
            [](RunnerCase& c) {
                writeScalarField(
                    c.fieldsDir / "ElectrostaticPotential_region0.csv",
                    {0, 0, 0, 0, 0},
                    "node_id,component0\n",
                    {0, 1, 2, 3});
            },
        },
    };

    for (const MalformedCase& malformed : cases) {
        INFO(malformed.name);
        RunnerCase c = makeRunnerCase(malformed.name);
        malformed.mutate(c);

        const RunnerOutput run = runCase(c.configPath);

        REQUIRE(run.exitCode != 0);
        REQUIRE(run.err.find(malformed.file) != std::string::npos);
        REQUIRE(run.err.find(malformed.expected) != std::string::npos);
    }
}