#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>
#include <cmath>

#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScaling.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/simulation/ConfigParsing.h"

using namespace vela;

namespace {

DeviceMesh makeMOSInterfaceMesh()
{
    DeviceMesh mesh;
    const double L = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0; n0.y = 0.0; mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = L;   n1.y = 0.0; mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.0; n2.y = L;   mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = L;   n3.y = L;   mesh.addNode(n3);

    Cell c0; c0.id = 0; c0.type = CellType::Tri3; c0.region_id = 0; c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);
    Cell c1; c1.id = 1; c1.type = CellType::Tri3; c1.region_id = 1; c1.node_ids = {1, 3, 2};
    mesh.addCell(c1);

    Region silicon; silicon.id = 0; silicon.name = "silicon"; silicon.material = "Si"; silicon.cell_ids = {0};
    mesh.addRegion(silicon);
    Region oxide; oxide.id = 1; oxide.name = "oxide"; oxide.material = "SiO2"; oxide.cell_ids = {1};
    mesh.addRegion(oxide);

    mesh.buildEdges();
    return mesh;
}

DopingModel makeZeroDoping(const DeviceMesh& mesh)
{
    return DopingModel::fromMeshAndRegions(
        mesh,
        {RegionDopingSpec{"silicon", 0.0, 0.0}, RegionDopingSpec{"oxide", 0.0, 0.0}});
}

} // namespace

TEST_CASE("DDAssembler: Poisson substep includes MOS interface sheet charge", "[interface][gummel][charge]")
{
    DeviceMesh mesh = makeMOSInterfaceMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroDoping(mesh);

    const int N = static_cast<int>(mesh.numNodes());
    const VectorXd n = VectorXd::Zero(N);
    const VectorXd p = VectorXd::Zero(N);
    const VectorXd psi = VectorXd::Zero(N);

    DDAssembler neutral(mesh, matdb, doping, constants::Vt_300, 1.0e-7, 1.0e-7);
    neutral.assemblePoissonWithCarriers(n, p, psi);

    const Real sheet = 2.0e15;
    DDAssembler charged(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        1.0e-7,
        1.0e-7,
        {},
        {InterfaceSheetChargeSpec{"silicon", "oxide", sheet}});
    charged.assemblePoissonWithCarriers(n, p, psi);

    const Real expectedEndpointCharge = constants::q * sheet * std::sqrt(2.0) * 1.0e-6 * 0.5;
    REQUIRE(charged.rhs()(1) - neutral.rhs()(1) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(charged.rhs()(2) - neutral.rhs()(2) == Catch::Approx(expectedEndpointCharge).epsilon(1e-12));
    REQUIRE(charged.rhs()(0) - neutral.rhs()(0) == Catch::Approx(0.0).margin(1e-30));
    REQUIRE(charged.rhs()(3) - neutral.rhs()(3) == Catch::Approx(0.0).margin(1e-30));
}

TEST_CASE("CoupledDDAssembler: Newton Poisson residual consumes fixed trap occupancy", "[interface][newton][charge]")
{
    DeviceMesh mesh = makeMOSInterfaceMesh();
    MaterialDatabase matdb;
    DopingModel doping = makeZeroDoping(mesh);

    CoupledDDAssembler neutral(mesh, matdb, doping, constants::Vt_300, 1.0e-7, 1.0e-7);
    CoupledDDAssembler trapped(
        mesh,
        matdb,
        doping,
        constants::Vt_300,
        1.0e-7,
        1.0e-7,
        {},
        {InterfaceSheetChargeSpec{"silicon", "oxide", 0.0, 0.0, 4.0e15, 0.25}});

    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::Zero(N);
    state.phin = VectorXd::Zero(N);
    state.phip = VectorXd::Zero(N);
    CoupledDDBoundaryConditions bcs;

    const VectorXd r0 = neutral.residual(neutral.pack(state), bcs);
    const VectorXd r1 = trapped.residual(trapped.pack(state), bcs);

    const Real effectiveSheet = 1.0e15;
    const Real expectedEndpointResidualShift = -constants::q * effectiveSheet * std::sqrt(2.0) * 1.0e-6 * 0.5;
    REQUIRE(r1(1) - r0(1) == Catch::Approx(expectedEndpointResidualShift).epsilon(1e-12));
    REQUIRE(r1(2) - r0(2) == Catch::Approx(expectedEndpointResidualShift).epsilon(1e-12));
}

TEST_CASE("ConfigParsing: interface trap occupancy outside unit interval is rejected", "[interface][traps][config]")
{
    const nlohmann::json cfg = {
        {"interfaces", {{{"regions", {"silicon", "oxide"}},
                         {"trap_density_m2", 1.0e15},
                         {"trap_occupancy", -0.1}}}}
    };

    REQUIRE_THROWS_WITH(
        parseInterfaceSheetChargeSpecs(cfg),
        Catch::Matchers::ContainsSubstring("trap_occupancy must be in [0, 1]"));
}


TEST_CASE("ConfigParsing: interface region selectors validate preferred and legacy forms", "[interface][config]")
{
    const nlohmann::json preferred = {
        {"interfaces", {{{"regions", {"silicon", "oxide", "metal"}},
                         {"fixed_charge_m2", 1.0e14}}}}
    };
    REQUIRE_THROWS_WITH(
        parseInterfaceSheetChargeSpecs(preferred),
        Catch::Matchers::ContainsSubstring("interface regions must contain exactly two names"));

    const nlohmann::json legacy = {
        {"interfaces", {{{"region0", "silicon"},
                         {"region1", "oxide"},
                         {"sheet_charge_m2", 2.0e14}}}}
    };
    const auto specs = parseInterfaceSheetChargeSpecs(legacy);
    REQUIRE(specs.size() == 1);
    REQUIRE(specs[0].region0 == "silicon");
    REQUIRE(specs[0].region1 == "oxide");
    REQUIRE(specs[0].sheetCharge == Catch::Approx(2.0e14));
}

TEST_CASE("ConfigParsing unit_scaling converts cm^-2 interface charge to m^-2",
          "[interface][config][scaling]")
{
    const nlohmann::json cfg = {
        {"scaling", {{"mode", "unit_scaling"}}},
        {"interfaces", {{{"regions", {"silicon", "oxide"}},
                         {"sheet_charge_m2", 2.0e11},
                         {"fixed_charge_m2", -3.0e11},
                         {"trap_density_m2", 4.0e11},
                         {"trap_occupancy", 0.25}}}}
    };

    const auto specs = parseInterfaceSheetChargeSpecs(cfg, parseUnitScalingConfig(cfg));
    REQUIRE(specs.size() == 1);
    REQUIRE(specs[0].sheetCharge == Catch::Approx(2.0e15));
    REQUIRE(specs[0].fixedCharge == Catch::Approx(-3.0e15));
    REQUIRE(specs[0].trapDensity == Catch::Approx(4.0e15));
    REQUIRE(specs[0].trapOccupancy == Catch::Approx(0.25));
}

TEST_CASE("ConfigParsing: duplicate region fixed charge sources are rejected", "[interface][config]")
{
    const nlohmann::json cfg = {
        {"doping", {{{"region", "silicon"},
                     {"donors", 0.0},
                     {"acceptors", 0.0},
                     {"fixed_charge_m3", 1.0e20}}}},
        {"regions", {{{"name", "silicon"},
                      {"material", "Si"},
                      {"fixed_charge_m3", 2.0e20}}}}
    };

    REQUIRE_THROWS_WITH(
        parseRegionFixedChargeSpecs(cfg),
        Catch::Matchers::ContainsSubstring("duplicate fixed_charge_m3 for region 'silicon'"));
}

TEST_CASE("ConfigParsing: trap occupancy requires density and is always bounded", "[interface][traps][config]")
{
    const nlohmann::json missingDensity = {
        {"interfaces", {{{"regions", {"silicon", "oxide"}},
                         {"fixed_charge_m2", 1.0e14},
                         {"trap_occupancy", 0.5}}}}
    };
    REQUIRE_THROWS_WITH(
        parseInterfaceSheetChargeSpecs(missingDensity),
        Catch::Matchers::ContainsSubstring("trap_occupancy requires trap_density_m2"));

    const nlohmann::json outOfRangeWithoutCharge = {
        {"interfaces", {{{"regions", {"silicon", "oxide"}},
                         {"trap_density_m2", 1.0e15},
                         {"trap_occupancy", 1.1}}}}
    };
    REQUIRE_THROWS_WITH(
        parseInterfaceSheetChargeSpecs(outOfRangeWithoutCharge),
        Catch::Matchers::ContainsSubstring("trap_occupancy must be in [0, 1]"));
}
