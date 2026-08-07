#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>

#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingProfileEvaluator.h"

using namespace vela;

namespace {

DeviceMesh buildSingleRegionMesh()
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0.0, 0.0});
    mesh.addNode(Node{1, 1.0, 0.0});
    mesh.addNode(Node{2, 1.0, 1.0});
    mesh.addNode(Node{3, 0.0, 1.0});

    Cell c0;
    c0.id = 0;
    c0.region_id = 0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Cell c1;
    c1.id = 1;
    c1.region_id = 0;
    c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region region;
    region.id = 0;
    region.name = "R.Si";
    region.material = "Silicon";
    region.cell_ids = {0, 1};
    mesh.addRegion(region);
    return mesh;
}

} // namespace

TEST_CASE("doping evaluator applies constant region profile", "[doping][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {}});

    DopingProfileIr profile;
    profile.name = "ConstN";
    profile.kind = DopingProfileKind::Constant;
    profile.targetRegion = "R.Si";
    profile.donors_cm3 = 1.0e17;
    ir.dopingProfiles.push_back(profile);

    const DeviceMesh mesh = buildSingleRegionMesh();
    const DopingModel model = DopingProfileEvaluator::evaluate(mesh, ir);

    for (Index i = 0; i < mesh.numNodes(); ++i) {
        CHECK(model.donors(i) == Catch::Approx(1.0e17));
        CHECK(model.acceptors(i) == Catch::Approx(0.0));
    }
}

TEST_CASE("doping evaluator applies gaussian acceptor profile", "[doping][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {}});

    DopingProfileIr profile;
    profile.name = "GaussP";
    profile.kind = DopingProfileKind::Gaussian;
    profile.targetRegion = "R.Si";
    profile.gaussianPeak_cm3 = 5.0e18;
    profile.gaussianValueAtDepth_cm3 = 1.0e16;
    profile.gaussianPeakPosUm = Point2D{0.0, 0.0};
    profile.gaussianSigmaXUm = 0.2 / std::sqrt(-2.0 * std::log(1.0e16 / 5.0e18));
    profile.gaussianActsOnDonors = false;
    ir.dopingProfiles.push_back(profile);

    const DeviceMesh mesh = buildSingleRegionMesh();
    const DopingModel model = DopingProfileEvaluator::evaluate(mesh, ir);

    CHECK(model.acceptors(0) == Catch::Approx(5.0e18));
    CHECK(model.acceptors(1) == Catch::Approx(5.0e18 * std::exp(-0.5 * std::pow(1.0 / profile.gaussianSigmaXUm, 2.0))).epsilon(1e-6));
    CHECK(model.acceptors(3) == Catch::Approx(5.0e18));
    CHECK(model.acceptors(2) == Catch::Approx(5.0e18 * std::exp(-0.5 * std::pow(1.0 / profile.gaussianSigmaXUm, 2.0))).epsilon(1e-6));
}

TEST_CASE("doping evaluator gaussian matches peak depth symmetry and ignores y", "[doping][preprocess]")
{
    DeviceMesh mesh;
    mesh.addNode(Node{0, -0.5, 0.0});
    mesh.addNode(Node{1, -0.3, 0.0});
    mesh.addNode(Node{2, -0.1, 0.0});
    mesh.addNode(Node{3, -0.3, 0.8});

    Cell c0;
    c0.id = 0;
    c0.region_id = 0;
    c0.node_ids = {0, 1, 2};
    mesh.addCell(c0);

    Cell c1;
    c1.id = 1;
    c1.region_id = 0;
    c1.node_ids = {0, 2, 3};
    mesh.addCell(c1);

    Region region;
    region.id = 0;
    region.name = "R.Si";
    region.material = "Silicon";
    region.cell_ids = {0, 1};
    mesh.addRegion(region);

    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {}});

    DopingProfileIr profile;
    profile.name = "GaussSym";
    profile.kind = DopingProfileKind::Gaussian;
    profile.targetRegion = "R.Si";
    profile.gaussianPeak_cm3 = 5.0e18;
    profile.gaussianValueAtDepth_cm3 = 1.0e16;
    profile.gaussianPeakPosUm = Point2D{-0.3, 0.0};
    profile.gaussianSigmaXUm = 0.2 / std::sqrt(-2.0 * std::log(1.0e16 / 5.0e18));
    profile.gaussianActsOnDonors = true;
    ir.dopingProfiles.push_back(profile);

    const DopingModel model = DopingProfileEvaluator::evaluate(mesh, ir);
    CHECK(model.donors(1) == Catch::Approx(profile.gaussianPeak_cm3));
    CHECK(model.donors(0) == Catch::Approx(profile.gaussianValueAtDepth_cm3).epsilon(1e-6));
    CHECK(model.donors(2) == Catch::Approx(profile.gaussianValueAtDepth_cm3).epsilon(1e-6));
    CHECK(model.donors(3) == Catch::Approx(profile.gaussianPeak_cm3));
}

TEST_CASE("doping evaluator rejects invalid gaussian parameters", "[doping][preprocess]")
{
    DeviceIr2D ir;
    ir.regions.push_back(RegionIr{"R.Si", "Silicon", {}});

    DopingProfileIr profile;
    profile.name = "BadGauss";
    profile.kind = DopingProfileKind::Gaussian;
    profile.targetRegion = "R.Si";
    profile.gaussianPeak_cm3 = 1.0e16;
    profile.gaussianValueAtDepth_cm3 = 1.0e16;
    profile.gaussianPeakPosUm = Point2D{0.0, 0.0};
    profile.gaussianSigmaXUm = 0.2;
    profile.gaussianActsOnDonors = true;
    ir.dopingProfiles.push_back(profile);

    REQUIRE_THROWS_WITH(
        DopingProfileEvaluator::evaluate(buildSingleRegionMesh(), ir),
        Catch::Matchers::ContainsSubstring("invalid gaussian parameters"));
}
