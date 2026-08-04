#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

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
    profile.gaussianBackground_cm3 = 1.0e16;
    profile.gaussianCenterUm = Point2D{0.0, 0.0};
    profile.gaussianSigmaXUm = 0.2;
    profile.gaussianSigmaYUm = 0.2;
    profile.gaussianActsOnDonors = false;
    ir.dopingProfiles.push_back(profile);

    const DeviceMesh mesh = buildSingleRegionMesh();
    const DopingModel model = DopingProfileEvaluator::evaluate(mesh, ir);

    CHECK(model.acceptors(0) > model.acceptors(2));
    CHECK(model.acceptors(0) > 1.0e18);
    CHECK(model.acceptors(2) > 0.0);
}

