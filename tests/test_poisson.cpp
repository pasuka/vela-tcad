#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/solver/LinearSolver.h"
#include "vela/io/VTKWriter.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/simulation/PoissonSimulation.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <nlohmann/json.hpp>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace vela;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/**
 * Build a simple 2-D p-n junction mesh:
 *
 *   3 -------- 2
 *   |  p-reg  /|
 *   |  (T1)  / |
 *   |       /  |
 *   |      /   |
 *   |     /    |
 *   |    /  T0 |
 *   |   / n-reg|
 *   |  /       |
 *   | /        |
 *   0 -------- 1
 *
 *  Nodes:  0=(0,0), 1=(1e-6,0), 2=(1e-6,1e-6), 3=(0,1e-6)  [1 um square]
 *  Cells:  T0={0,1,2} region 0 (n-Si),  T1={0,2,3} region 1 (p-Si)
 *  Contacts:
 *    cathode (n): nodes 1, 2   bias = 0 V
 *    anode   (p): nodes 0, 3   bias = 0 V
 */
static DeviceMesh makePNMesh(const std::string& material = "Si")
{
    DeviceMesh mesh;

    const double L = 1.0e-6; // 1 um

    Node n0; n0.id=0; n0.x=0;  n0.y=0;  mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L;  n1.y=0;  mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=L;  n2.y=L;  mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=0;  n3.y=L;  mesh.addNode(n3);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=1;
    c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region r0; r0.id=0; r0.name="n_region"; r0.material=material; r0.cell_ids={0};
    mesh.addRegion(r0);

    Region r1; r1.id=1; r1.name="p_region"; r1.material=material; r1.cell_ids={1};
    mesh.addRegion(r1);

    Contact anode;   anode.id=0;   anode.name="anode";   anode.region_id=1;
    anode.node_ids = {0, 3};
    mesh.addContact(anode);

    Contact cathode; cathode.id=1; cathode.name="cathode"; cathode.region_id=0;
    cathode.node_ids = {1, 2};
    mesh.addContact(cathode);

    mesh.buildEdges();
    return mesh;
}

// Build a DopingModel for the p-n mesh:
//   n_region: Nd = 1e23 m^-3, Na = 0
//   p_region: Nd = 0,         Na = 1e23 m^-3
static DopingModel makePNDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        { "n_region", 1e23, 0.0 },
        { "p_region", 0.0,  1e23 }
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static DeviceMesh makeMOSCapChargeMesh();
static DopingModel makeZeroMOSCapDoping(const DeviceMesh& mesh);

static void writePNMeshJson(const std::filesystem::path& meshPath,
                            const std::string& material = "Si")
{
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0}, {"y", 0.0}},
            {{"id", 1}, {"x", 1.0e-6}, {"y", 0.0}},
            {{"id", 2}, {"x", 1.0e-6}, {"y", 1.0e-6}},
            {{"id", 3}, {"x", 0.0}, {"y", 1.0e-6}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 2}}},
            {{"id", 1}, {"region_id", 1}, {"node_ids", {0, 2, 3}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "n_region"}, {"material", material}, {"cell_ids", {0}}},
            {{"id", 1}, {"name", "p_region"}, {"material", material}, {"cell_ids", {1}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "anode"}, {"region_id", 1}, {"node_ids", {0, 3}}},
            {{"id", 1}, {"name", "cathode"}, {"region_id", 0}, {"node_ids", {1, 2}}}
        }}
    };
    std::ofstream(meshPath) << mesh.dump(2);
}

static std::filesystem::path writePNMeshJsonToDir(const std::filesystem::path& dir,
                                                  const std::string& material = "Si")
{
    const auto meshPath = dir / "mesh.json";
    writePNMeshJson(meshPath, material);
    return meshPath;
}

static std::filesystem::path writeMOSCapMeshJsonToDir(const std::filesystem::path& dir,
                                                      Real lengthScale = 1.0e-6)
{
    const auto meshPath = dir / "moscap_mesh.json";
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0}, {"y", 0.0}},
            {{"id", 1}, {"x", 1.0 * lengthScale}, {"y", 0.0}},
            {{"id", 2}, {"x", 0.0}, {"y", 1.0 * lengthScale}},
            {{"id", 3}, {"x", 1.0 * lengthScale}, {"y", 1.0 * lengthScale}},
            {{"id", 4}, {"x", 0.0}, {"y", 2.0 * lengthScale}},
            {{"id", 5}, {"x", 1.0 * lengthScale}, {"y", 2.0 * lengthScale}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 3}}},
            {{"id", 1}, {"region_id", 0}, {"node_ids", {0, 3, 2}}},
            {{"id", 2}, {"region_id", 1}, {"node_ids", {2, 3, 5}}},
            {{"id", 3}, {"region_id", 1}, {"node_ids", {2, 5, 4}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "silicon"}, {"material", "Si"}, {"cell_ids", {0, 1}}},
            {{"id", 1}, {"name", "oxide"}, {"material", "SiO2"}, {"cell_ids", {2, 3}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "body"}, {"region_id", 0}, {"node_ids", {0, 1}}},
            {{"id", 1}, {"name", "gate"}, {"region_id", 1}, {"node_ids", {4, 5}}}
        }}
    };
    std::ofstream(meshPath) << mesh.dump(2);
    return meshPath;
}

struct ScopedPoissonTempDir {
    std::filesystem::path path;

    explicit ScopedPoissonTempDir(std::filesystem::path dir)
        : path(std::move(dir))
    {
    }

    ScopedPoissonTempDir(const ScopedPoissonTempDir&) = delete;
    ScopedPoissonTempDir& operator=(const ScopedPoissonTempDir&) = delete;
    ScopedPoissonTempDir(ScopedPoissonTempDir&& other) noexcept
        : path(std::move(other.path))
    {
        other.path.clear();
    }

    ScopedPoissonTempDir& operator=(ScopedPoissonTempDir&& other) noexcept
    {
        if (this != &other) {
            path = std::move(other.path);
            other.path.clear();
        }
        return *this;
    }

    ~ScopedPoissonTempDir()
    {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
};

static ScopedPoissonTempDir makePoissonTempDir(const std::string& name)
{
    constexpr int maxAttempts = 8;
    const auto base = std::filesystem::temp_directory_path();
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto tid = std::hash<std::thread::id>{}(std::this_thread::get_id());
    std::mt19937_64 rng(static_cast<std::mt19937_64::result_type>(stamp ^ tid));
    std::uniform_int_distribution<unsigned long long> dist;

    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        const auto dir = base /
            (name + "_" + std::to_string(stamp) + "_" + std::to_string(dist(rng)));
        std::error_code ec;
        if (std::filesystem::create_directory(dir, ec))
            return ScopedPoissonTempDir(dir);
    }

    throw std::runtime_error("Failed to create a unique temp directory for Poisson test.");
}

static void writePoissonConfigJson(const std::filesystem::path& cfgPath,
                                   const std::filesystem::path& meshPath,
                                   const std::filesystem::path& outputPath,
                                   const std::filesystem::path& materialsPath = {})
{
    nlohmann::json cfg = {
        {"mesh_file", meshPath.filename().string()},
        {"output_vtk", outputPath.filename().string()},
        {"doping", {
            {{"region", "n_region"}, {"donors", 1.0e23}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e23}}
        }},
        {"contacts", {
            {{"name", "anode"}, {"bias", 0.0}},
            {{"name", "cathode"}, {"bias", 0.0}}
        }}
    };
    if (!materialsPath.empty())
        cfg["materials_file"] = materialsPath.filename().string();
    std::ofstream(cfgPath) << cfg.dump(2);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("DopingModel: net doping has correct sign per region", "[doping]")
{
    DeviceMesh mesh = makePNMesh();
    DopingModel doping = makePNDoping(mesh);

    REQUIRE(doping.numNodes() == mesh.numNodes());

    // Node 1 belongs only to n_region cell -> positive net doping
    REQUIRE(doping.netDoping(1) > 0.0);

    // Node 3 belongs only to p_region cell -> negative net doping
    REQUIRE(doping.netDoping(3) < 0.0);
}

TEST_CASE("DopingModel: total impurity is separate from signed net doping", "[doping]")
{
    DopingModel doping(1);
    doping.setNodeDoping(0, 1.0e23, 1.0e23);

    REQUIRE(doping.donors(0) == Catch::Approx(1.0e23));
    REQUIRE(doping.acceptors(0) == Catch::Approx(1.0e23));
    REQUIRE(doping.netDoping(0) == Catch::Approx(0.0));
    REQUIRE(doping.totalImpurity(0) == Catch::Approx(2.0e23));
}

TEST_CASE("ConfigParsing unit_scaling converts cm^-3 doping to m^-3",
          "[poisson][doping][scaling]")
{
    DeviceMesh mesh = makePNMesh();
    const nlohmann::json legacyCfg = {
        {"doping", {
            {{"region", "n_region"}, {"donors", 1.0e23}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e23}},
        }}
    };
    const nlohmann::json scaledCfg = {
        {"scaling", {{"mode", "unit_scaling"}}},
        {"doping", {
            {{"region", "n_region"}, {"donors", 1.0e17}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e17}},
        }}
    };

    const DopingModel legacy = DopingModel::fromMeshAndRegions(
        mesh, parseDopingSpecs(legacyCfg, parseUnitScalingConfig(legacyCfg)));
    const DopingModel scaled = DopingModel::fromMeshAndRegions(
        mesh, parseDopingSpecs(scaledCfg, parseUnitScalingConfig(scaledCfg)));

    REQUIRE(scaled.numNodes() == legacy.numNodes());
    for (Index i = 0; i < legacy.numNodes(); ++i) {
        REQUIRE(scaled.donors(i) == Catch::Approx(legacy.donors(i)));
        REQUIRE(scaled.acceptors(i) == Catch::Approx(legacy.acceptors(i)));
        REQUIRE(scaled.netDoping(i) == Catch::Approx(legacy.netDoping(i)));
    }
}

TEST_CASE("UnitScalingConfig rejects scaling.system aliases",
          "[poisson][config][scaling]")
{
    REQUIRE_THROWS_WITH(
        parseUnitScalingConfig(nlohmann::json{{"scaling", {{"system", "si"}}}}),
        Catch::Matchers::ContainsSubstring("scaling.system"));
    REQUIRE_THROWS_WITH(
        parseUnitScalingConfig(nlohmann::json{{"scaling", {{"system", "sentaurus"}}}}),
        Catch::Matchers::ContainsSubstring("scaling.system"));
    REQUIRE_THROWS_WITH(
        parseUnitScalingConfig(nlohmann::json{{"scaling", {{"mode", "si"}}}}),
        Catch::Matchers::ContainsSubstring("unit_scaling"));
}

TEST_CASE("PoissonAssembler: matrix dimensions match node count", "[poisson]")
{
    DeviceMesh      mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    const Index N = mesh.numNodes();
    REQUIRE(asm_.matrix().rows() == static_cast<int>(N));
    REQUIRE(asm_.matrix().cols() == static_cast<int>(N));
    REQUIRE(asm_.rhs().size()    == static_cast<int>(N));
}

TEST_CASE("PoissonAssembler: scaled RHS maps to physical RHS by eps_ref*V0",
          "[poisson][scaling][boundary][interface]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    std::vector<RegionFixedChargeSpec> fixedCharges = {
        RegionFixedChargeSpec{"oxide", 1.2e21}
    };
    std::vector<InterfaceSheetChargeSpec> sheetCharges = {
        InterfaceSheetChargeSpec{"silicon", "oxide", 2.5e15, -1.0e15, 4.0e15, 0.25}
    };
    std::vector<PoissonNeumannBoundarySpec> neumannSpecs = {
        PoissonNeumannBoundarySpec{{2, 3}, 3.0e-8}
    };

    PoissonAssembler legacy(mesh, matdb, doping, fixedCharges, sheetCharges, neumannSpecs);
    legacy.assemble();

    const Real epsRef = constants::eps0 * 11.7;
    const Real V0 = constants::Vt_300;
    PoissonScalingSpec scaledSpec;
    scaledSpec.enabled = true;
    scaledSpec.potentialScale_V = V0;
    scaledSpec.permittivityReference_F_per_m = epsRef;
    PoissonAssembler scaled(mesh,
                            matdb,
                            doping,
                            fixedCharges,
                            sheetCharges,
                            neumannSpecs,
                            scaledSpec);
    scaled.assemble();

    REQUIRE(scaled.rhs().size() == legacy.rhs().size());
    for (int i = 0; i < legacy.rhs().size(); ++i) {
        REQUIRE(scaled.rhs()(i) * epsRef * V0
                == Catch::Approx(legacy.rhs()(i)).epsilon(1e-10));
    }
}

TEST_CASE("PoissonAssembler: scaled solve recovers legacy physical potential",
          "[poisson][scaling]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    std::vector<RegionFixedChargeSpec> fixedCharges = {
        RegionFixedChargeSpec{"oxide", 5.0e20}
    };
    std::vector<InterfaceSheetChargeSpec> sheetCharges = {
        InterfaceSheetChargeSpec{"silicon", "oxide", 2.0e15}
    };

    std::unordered_map<Index, Real> bcs = {
        {0, -0.15}, {1, -0.15}, {4, 0.35}, {5, 0.35}
    };

    PoissonAssembler legacy(mesh, matdb, doping, fixedCharges, sheetCharges);
    legacy.assemble();
    legacy.applyDirichlet(bcs);
    LinearSolver solver;
    const VectorXd psiLegacy = solver.solve(legacy.matrix(), legacy.rhs());

    const Real epsRef = constants::eps0 * 11.7;
    const Real V0 = constants::Vt_300;
    PoissonScalingSpec scaledSpec;
    scaledSpec.enabled = true;
    scaledSpec.potentialScale_V = V0;
    scaledSpec.permittivityReference_F_per_m = epsRef;
    PoissonAssembler scaled(mesh,
                            matdb,
                            doping,
                            fixedCharges,
                            sheetCharges,
                            {},
                            scaledSpec);
    scaled.assemble();
    scaled.applyDirichlet(bcs);
    const VectorXd psiHat = solver.solve(scaled.matrix(), scaled.rhs());

    REQUIRE(psiHat.size() == psiLegacy.size());
    for (int i = 0; i < psiLegacy.size(); ++i) {
        REQUIRE(psiHat(i) * V0 == Catch::Approx(psiLegacy(i)).epsilon(1e-10));
    }
}

TEST_CASE("PoissonAssembler + LinearSolver: solve succeeds, no NaN", "[poisson]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    // Apply Dirichlet: both contacts at 0 V
    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c) {
        const Contact& ct = mesh.getContact(c);
        for (Index nid : ct.node_ids)
            bcs[nid] = 0.0;
    }
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(asm_.matrix(), asm_.rhs());

    // Solution length must match number of nodes
    REQUIRE(psi.size() == static_cast<int>(mesh.numNodes()));

    // No NaN values
    for (int i = 0; i < psi.size(); ++i)
        REQUIRE_FALSE(std::isnan(psi(i)));
}

TEST_CASE("PoissonAssembler: Dirichlet nodes take prescribed value", "[poisson]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    // Prescribe 0 V on all boundary nodes
    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c)
        for (Index nid : mesh.getContact(c).node_ids)
            bcs[nid] = 0.0;
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(asm_.matrix(), asm_.rhs());

    // All Dirichlet nodes must return exactly the prescribed value
    for (const auto& [nid, val] : bcs)
        REQUIRE(psi(static_cast<int>(nid)) == Catch::Approx(val).margin(1e-12));
}

TEST_CASE("VTKWriter: writes file with potential field", "[poisson][vtk]")
{
    DeviceMesh       mesh   = makePNMesh();
    MaterialDatabase matdb;
    DopingModel      doping = makePNDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping);
    asm_.assemble();

    std::unordered_map<Index, Real> bcs;
    for (Index c = 0; c < mesh.numContacts(); ++c)
        for (Index nid : mesh.getContact(c).node_ids)
            bcs[nid] = 0.0;
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    VectorXd psi = solver.solve(asm_.matrix(), asm_.rhs());

    // Write to a temporary file
    const std::string vtkPath =
        (std::filesystem::temp_directory_path() / "test_poisson_out.vtk").string();
    VTKWriter writer(vtkPath, mesh);
    writer.write();

    std::vector<Real> psiVec(mesh.numNodes());
    for (Index i = 0; i < mesh.numNodes(); ++i)
        psiVec[i] = psi(static_cast<int>(i));
    writer.addNodeScalar("potential_V", psiVec);

    // File must exist and be non-empty
    REQUIRE(std::filesystem::exists(vtkPath));
    REQUIRE(std::filesystem::file_size(vtkPath) > 0);

    // File must contain the field name
    std::ifstream ifs(vtkPath);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    REQUIRE(content.find("potential_V") != std::string::npos);
}



static DeviceMesh makeMOSCapChargeMesh()
{
    DeviceMesh mesh;

    const double L = 1.0e-6;
    const double H = 1.0e-6;

    Node n0; n0.id=0; n0.x=0; n0.y=0;     mesh.addNode(n0);
    Node n1; n1.id=1; n1.x=L; n1.y=0;     mesh.addNode(n1);
    Node n2; n2.id=2; n2.x=0; n2.y=H;     mesh.addNode(n2);
    Node n3; n3.id=3; n3.x=L; n3.y=H;     mesh.addNode(n3);
    Node n4; n4.id=4; n4.x=0; n4.y=2*H;   mesh.addNode(n4);
    Node n5; n5.id=5; n5.x=L; n5.y=2*H;   mesh.addNode(n5);

    Cell c0; c0.id=0; c0.type=CellType::Tri3; c0.region_id=0; c0.node_ids={0,1,3}; mesh.addCell(c0);
    Cell c1; c1.id=1; c1.type=CellType::Tri3; c1.region_id=0; c1.node_ids={0,3,2}; mesh.addCell(c1);
    Cell c2; c2.id=2; c2.type=CellType::Tri3; c2.region_id=1; c2.node_ids={2,3,5}; mesh.addCell(c2);
    Cell c3; c3.id=3; c3.type=CellType::Tri3; c3.region_id=1; c3.node_ids={2,5,4}; mesh.addCell(c3);

    Region silicon; silicon.id=0; silicon.name="silicon"; silicon.material="Si"; silicon.cell_ids={0,1};
    mesh.addRegion(silicon);

    Region oxide; oxide.id=1; oxide.name="oxide"; oxide.material="SiO2"; oxide.cell_ids={2,3};
    mesh.addRegion(oxide);

    Contact body; body.id=0; body.name="body"; body.region_id=0; body.node_ids={0,1};
    mesh.addContact(body);

    Contact gate; gate.id=1; gate.name="gate"; gate.region_id=1; gate.node_ids={4,5};
    mesh.addContact(gate);

    mesh.buildEdges();
    return mesh;
}

static DopingModel makeZeroMOSCapDoping(const DeviceMesh& mesh)
{
    std::vector<RegionDopingSpec> specs = {
        {"silicon", 0.0, 0.0},
        {"oxide", 0.0, 0.0}
    };
    return DopingModel::fromMeshAndRegions(mesh, specs);
}

static VectorXd solveMOSCapChargeCase(
    const std::vector<RegionFixedChargeSpec>& fixedCharges,
    const std::vector<InterfaceSheetChargeSpec>& sheetCharges = {})
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler asm_(mesh, matdb, doping, fixedCharges, sheetCharges);
    asm_.assemble();

    std::unordered_map<Index, Real> bcs = {
        {0, 0.0}, {1, 0.0}, {4, 0.0}, {5, 0.0}
    };
    asm_.applyDirichlet(bcs);

    LinearSolver solver;
    return solver.solve(asm_.matrix(), asm_.rhs());
}

TEST_CASE("PoissonAssembler: zero explicit charge matches legacy RHS", "[poisson][charge]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler legacy(mesh, matdb, doping);
    legacy.assemble();

    PoissonAssembler explicitZero(
        mesh,
        matdb,
        doping,
        {RegionFixedChargeSpec{"oxide", 0.0}},
        {InterfaceSheetChargeSpec{"silicon", "oxide", 0.0}});
    explicitZero.assemble();

    REQUIRE(explicitZero.rhs().size() == legacy.rhs().size());
    for (int i = 0; i < legacy.rhs().size(); ++i)
        REQUIRE(explicitZero.rhs()(i) == Catch::Approx(legacy.rhs()(i)).margin(1e-30));
}

TEST_CASE("PoissonSimulation: duplicate fixed charge config entries are rejected", "[poisson][charge]")
{
    const auto tempDir = makePoissonTempDir("vela_poisson_duplicate_fixed_charge_test");
    const auto& dir = tempDir.path;

    const auto meshPath = writePNMeshJsonToDir(dir);
    const auto configPath = dir / "poisson_duplicate_fixed_charge.json";
    const nlohmann::json cfg = {
        {"mesh_file", meshPath.string()},
        {"output_vtk", (dir / "out.vtk").string()},
        {"doping", {
            {{"region", "n_region"},
             {"donors", 1.0e23},
             {"acceptors", 0.0},
             {"fixed_charge_m3", 1.0e20}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e23}}
        }},
        {"regions", {
            {{"name", "n_region"}, {"fixed_charge_m3", 2.0e20}}
        }},
        {"contacts", {
            {{"name", "anode"}, {"bias", 0.0}},
            {{"name", "cathode"}, {"bias", 0.0}}
        }}
    };
    std::ofstream(configPath) << cfg.dump(2);

    PoissonSimulation sim;
    REQUIRE_THROWS_AS(sim.runWithResult(configPath.string()), std::runtime_error);
}


TEST_CASE("PoissonSimulation: trap occupancy outside unit interval is rejected", "[poisson][charge][traps]")
{
    const auto tempDir = makePoissonTempDir("vela_poisson_bad_trap_occupancy_test");
    const auto& dir = tempDir.path;

    const auto meshPath = writePNMeshJsonToDir(dir);
    const auto configPath = dir / "poisson_bad_trap_occupancy.json";
    nlohmann::json cfg = {
        {"mesh_file", meshPath.string()},
        {"output_vtk", (dir / "out.vtk").string()},
        {"doping", {
            {{"region", "n_region"}, {"donors", 1.0e23}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e23}}
        }},
        {"interfaces", {
            {{"regions", {"n_region", "p_region"}},
             {"trap_density_m2", 1.0e15},
             {"trap_occupancy", 1.25}}
        }},
        {"contacts", {
            {{"name", "anode"}, {"bias", 0.0}},
            {{"name", "cathode"}, {"bias", 0.0}}
        }}
    };
    std::ofstream(configPath) << cfg.dump(2);

    PoissonSimulation sim;
    REQUIRE_THROWS_WITH(
        sim.runWithResult(configPath.string()),
        Catch::Matchers::ContainsSubstring("trap_occupancy must be in [0, 1]"));
}

TEST_CASE("PoissonAssembler: duplicate fixed charge specs are rejected", "[poisson][charge]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler asm_(
        mesh,
        matdb,
        doping,
        {RegionFixedChargeSpec{"oxide", 1.0e21},
         RegionFixedChargeSpec{"oxide", 2.0e21}});

    REQUIRE_THROWS_AS(asm_.assemble(), std::invalid_argument);
}

TEST_CASE("PoissonAssembler: fixed charge sign shifts MOS capacitor potential", "[poisson][charge]")
{
    const VectorXd zero = solveMOSCapChargeCase({});
    const VectorXd positive = solveMOSCapChargeCase({RegionFixedChargeSpec{"oxide", 1.0e21}});
    const VectorXd negative = solveMOSCapChargeCase({RegionFixedChargeSpec{"oxide", -1.0e21}});

    const double zeroInterface = 0.5 * (zero(2) + zero(3));
    const double positiveInterface = 0.5 * (positive(2) + positive(3));
    const double negativeInterface = 0.5 * (negative(2) + negative(3));

    REQUIRE(zeroInterface == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(positiveInterface > zeroInterface);
    REQUIRE(negativeInterface < zeroInterface);
    REQUIRE(positiveInterface == Catch::Approx(-negativeInterface).epsilon(1e-12));
}

TEST_CASE("PoissonAssembler: sheet charge is split to shared interface endpoints", "[poisson][charge]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler asm_(
        mesh,
        matdb,
        doping,
        {},
        {InterfaceSheetChargeSpec{"silicon", "oxide", 1.25e15},
         InterfaceSheetChargeSpec{"oxide", "silicon", 0.75e15}});
    asm_.assemble();

    const Real expectedEndpointCharge = constants::q * 2.0e15 * 1.0e-6 * 0.5;
    REQUIRE(asm_.rhs()(2) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(asm_.rhs()(3) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(asm_.rhs()(0) == Catch::Approx(0.0).margin(1e-30));
    REQUIRE(asm_.rhs()(5) == Catch::Approx(0.0).margin(1e-30));
}


TEST_CASE("PoissonAssembler: interface fixed charge and traps add to sheet charge", "[poisson][charge][traps]")
{
    DeviceMesh mesh = makeMOSCapChargeMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroMOSCapDoping(mesh);

    PoissonAssembler asm_(
        mesh,
        matdb,
        doping,
        {},
        {InterfaceSheetChargeSpec{"silicon", "oxide", 1.0e15, 2.0e15, 4.0e15, 0.25}});
    asm_.assemble();

    const Real totalCharge = 1.0e15 + 2.0e15 + 4.0e15 * 0.25;
    const Real expectedEndpointCharge = constants::q * totalCharge * 1.0e-6 * 0.5;
    REQUIRE(asm_.rhs()(2) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(asm_.rhs()(3) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
}

TEST_CASE("MaterialDatabase: external file overrides built-ins", "[material]")
{
    const auto tempDir = makePoissonTempDir("vela_material_override_test");
    const auto& dir = tempDir.path;
    const auto materialsPath = dir / "materials.json";
    std::ofstream(materialsPath) << nlohmann::json{
        {"materials", {
            {{"name", "Si"}, {"eps_r", 12.5}, {"bandgap_eV", 1.11}, {"Nc_m3", 3.0e25}}
        }}
    }.dump(2);

    MaterialDatabase matdb;
    REQUIRE(matdb.getMaterial("Si").eps_r == Catch::Approx(11.7));
    matdb.loadJson(materialsPath.string());

    const Material& si = matdb.getMaterial("Si");
    REQUIRE(si.eps_r == Catch::Approx(12.5));
    REQUIRE(si.mun == Catch::Approx(0.135));
    REQUIRE(si.bandgap_eV.has_value());
    REQUIRE(*si.bandgap_eV == Catch::Approx(1.11));
    REQUIRE(si.Nc_m3.has_value());
    REQUIRE(*si.Nc_m3 == Catch::Approx(3.0e25));
    REQUIRE(matdb.hasMaterial("SiO2"));

}


TEST_CASE("MaterialDatabase: malformed materials array reports file path", "[material]")
{
    const auto tempDir = makePoissonTempDir("vela_material_malformed_test");
    const auto& dir = tempDir.path;
    const auto materialsPath = dir / "materials.json";
    std::ofstream(materialsPath) << nlohmann::json{{"materials", "not an array"}}.dump(2);

    MaterialDatabase matdb;
    try {
        matdb.loadJson(materialsPath.string());
        FAIL("Expected malformed materials JSON to throw");
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        REQUIRE(message.find(materialsPath.string()) != std::string::npos);
        REQUIRE(message.find("'materials' must be an array") != std::string::npos);
    }
}

TEST_CASE("PoissonAssembler: unknown material reports an error", "[poisson][material]")
{
    DeviceMesh mesh = makePNMesh("Unobtainium");

    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    PoissonAssembler asm_(mesh, matdb, doping);

    REQUIRE_THROWS_WITH(asm_.assemble(), Catch::Matchers::ContainsSubstring("unknown material"));
}

TEST_CASE("PoissonSimulation: external new material can be used for assembly", "[poisson][material]")
{
    const auto tempDir = makePoissonTempDir("vela_new_material_poisson_test");
    const auto& dir = tempDir.path;
    const auto meshPath = dir / "pn_mesh.json";
    const auto cfgPath = dir / "poisson.json";
    const auto outputPath = dir / "out.vtk";
    const auto materialsPath = dir / "materials.json";

    writePNMeshJson(meshPath, "GaAsLike");
    std::ofstream(materialsPath) << nlohmann::json{
        {"materials", {
            {{"name", "GaAsLike"}, {"eps_r", 12.9}, {"ni", 2.0e12},
             {"mun", 0.85}, {"mup", 0.04}, {"bandgap_eV", 1.42},
             {"electron_affinity_eV", 4.07}, {"Nc_m3", 4.7e23},
             {"Nv_m3", 7.0e24}, {"temperature_K", 300.0}}
        }}
    }.dump(2);
    writePoissonConfigJson(cfgPath, meshPath, outputPath, materialsPath);

    PoissonSimulation sim;
    PoissonResult result = sim.runWithResult(cfgPath.string());

    REQUIRE(result.potential.size() == static_cast<int>(result.mesh.numNodes()));
    REQUIRE(std::filesystem::exists(outputPath));
    for (int i = 0; i < result.potential.size(); ++i)
        REQUIRE(std::isfinite(result.potential(i)));

}

TEST_CASE("PoissonSimulation: unit_scaling Neumann displacement matches legacy physical result",
          "[poisson][boundary][scaling]")
{
    const auto tempLegacy = makePoissonTempDir("vela_poisson_legacy_neumann_scaling");
    const auto tempScaled = makePoissonTempDir("vela_poisson_unit_neumann_scaling");

    const auto meshLegacy = writeMOSCapMeshJsonToDir(tempLegacy.path, 1.0e-6);
    const auto meshScaled = writeMOSCapMeshJsonToDir(tempScaled.path, 1.0);
    const auto cfgLegacyPath = tempLegacy.path / "poisson_legacy.json";
    const auto cfgScaledPath = tempScaled.path / "poisson_unit_scaling.json";
    const auto outLegacy = tempLegacy.path / "out.vtk";
    const auto outScaled = tempScaled.path / "out.vtk";

    const nlohmann::json legacyCfg = {
        {"mesh_file", meshLegacy.filename().string()},
        {"output_vtk", outLegacy.filename().string()},
        {"doping", {
            {{"region", "silicon"}, {"donors", 0.0}, {"acceptors", 0.0}},
            {{"region", "oxide"}, {"donors", 0.0}, {"acceptors", 0.0}}
        }},
        {"contacts", {
            {{"name", "body"}, {"bias", 0.0}},
            {{"name", "gate"}, {"bias", 0.0}}
        }},
        {"boundaries", {
            {{"name", "midline_neumann"},
             {"type", "neumann"},
             {"node_ids", {2, 3}},
             {"normal_displacement_C_per_m2", 1.0e-8}}
        }}
    };

    const nlohmann::json scaledCfg = {
        {"scaling", {{"mode", "unit_scaling"}}},
        {"mesh_file", meshScaled.filename().string()},
        {"output_vtk", outScaled.filename().string()},
        {"doping", {
            {{"region", "silicon"}, {"donors", 0.0}, {"acceptors", 0.0}},
            {{"region", "oxide"}, {"donors", 0.0}, {"acceptors", 0.0}}
        }},
        {"contacts", {
            {{"name", "body"}, {"bias", 0.0}},
            {{"name", "gate"}, {"bias", 0.0}}
        }},
        {"boundaries", {
            {{"name", "midline_neumann"},
             {"type", "neumann"},
             {"node_ids", {2, 3}},
             {"normal_displacement_C_per_m2", 1.0e-12}}
        }}
    };

    std::ofstream(cfgLegacyPath) << legacyCfg.dump(2);
    std::ofstream(cfgScaledPath) << scaledCfg.dump(2);

    PoissonSimulation sim;
    const PoissonResult legacy = sim.runWithResult(cfgLegacyPath.string());
    const PoissonResult scaled = sim.runWithResult(cfgScaledPath.string());

    REQUIRE(scaled.potential.size() == legacy.potential.size());
    for (int i = 0; i < legacy.potential.size(); ++i) {
        REQUIRE(scaled.potential(i) == Catch::Approx(legacy.potential(i)).epsilon(1e-6));
    }
}
