#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/equation/ExtendedPoissonResidual.h"
#include "vela/equation/DDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/physics/DopingModel.h"
#include "vela/simulation/ConfigParsing.h"

#include <limits>
#include <nlohmann/json.hpp>

using namespace vela;

namespace {

DeviceMesh makeHorizontalInterfaceMesh(bool skew = false, bool obtuse = false)
{
    DeviceMesh mesh;
    const Real length = 1.0e-6;

    Node n0; n0.id = 0; n0.x = 0.0;          n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = length;       n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.5 * length; n2.y = -length; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = (skew ? 0.7 : 0.5) * length; n3.y = (obtuse ? 0.2 : skew ? 2.0 : 1.0) * length; mesh.addNode(n3);

    Cell siliconCell;
    siliconCell.id = 0;
    siliconCell.type = CellType::Tri3;
    siliconCell.region_id = 0;
    siliconCell.node_ids = {0, 1, 2};
    mesh.addCell(siliconCell);

    Cell oxideCell;
    oxideCell.id = 1;
    oxideCell.type = CellType::Tri3;
    oxideCell.region_id = 1;
    oxideCell.node_ids = {0, 3, 1};
    mesh.addCell(oxideCell);

    Region silicon;
    silicon.id = 0;
    silicon.name = "channel";
    silicon.material = "Si";
    silicon.cell_ids = {0};
    mesh.addRegion(silicon);

    Region oxide;
    oxide.id = 1;
    oxide.name = "gate_oxide";
    oxide.material = "SiO2";
    oxide.cell_ids = {1};
    mesh.addRegion(oxide);

    mesh.buildEdges();
    return mesh;
}

DeviceMesh makeSiliconTriangleMesh(bool obtuse = false)
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0;    n0.y = 0.0;    mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0e-6; n1.y = 0.0;    mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = obtuse ? 0.5e-6 : 0.0; n2.y = obtuse ? 0.1e-6 : 1.0e-6; mesh.addNode(n2);
    Cell cell;
    cell.id = 0;
    cell.type = CellType::Tri3;
    cell.region_id = 0;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);
    Region silicon;
    silicon.id = 0;
    silicon.name = "silicon";
    silicon.material = "Si";
    silicon.cell_ids = {0};
    mesh.addRegion(silicon);
    mesh.buildEdges();
    return mesh;
}

} // namespace

TEST_CASE("Extended Poisson matches production material reference charge and leaves carrier rows unchanged",
          "[poisson][precision][material_reference]")
{
    const auto mesh=makeHorizontalInterfaceMesh();
    MaterialDatabase db;
    auto material=db.getMaterial("Si");
    material.intrinsic_reference_potential_V=GENERATE(-.008,0.,.007927767678211338);
    db.addMaterial(material);
    const auto doping=DopingModel::fromMeshAndRegions(mesh,{{"channel",1e17,0.}});
    DDScalingSpec scaling;
    scaling.poissonChargeVolumePolicy=PoissonChargeVolumePolicy::MaterialLocalCircumcentric;
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,{}, {}, {}, {}, {}, {},scaling);
    CoupledDDState state{VectorXd::Constant(4,.08),VectorXd::Constant(4,.01),VectorXd::Constant(4,-.02)};
    state.psi(1)=.09;state.psi(3)=.07;
    const auto x=assembler.pack(state);
    CoupledDDBoundaryConditions bc;bc.psi[3]=.06;bc.phin[2]=.013;bc.phip[2]=-.017;
    const auto original=assembler.residual(x,bc);
    const auto originalJ=assembler.assembleJacobian(x,bc);
    assembler.setExtendedPoissonResidual(true);
    const auto extended=assembler.residual(x,bc);
    const auto model=assembler.poissonArithmeticModel(x,bc);
    const auto exported=poisson_precision::evaluate(model);
    for(int i=0;i<4;++i)REQUIRE(static_cast<Real>(exported[i])==extended[i]);
    // A sub-ULP free-potential increment must remain visible to both dielectric
    // flux and charge, with the derivative of the same physical Jacobian.
    std::vector<double> low(4,0.);low[0]=1e-18;
    REQUIRE(x[0]+low[0]==x[0]);
    const auto shifted=poisson_precision::evaluate(model,low);
    for(int i=0;i<3;++i) {
        const auto derivative=(shifted[i]-exported[i])/poisson_precision::Q(low[0]);
        REQUIRE(static_cast<Real>(derivative)==Catch::Approx(originalJ.coeff(i,0)).epsilon(1e-11).margin(1e-25));
    }
    REQUIRE(shifted[3]==exported[3]);
    low[3]=1e-18;
    REQUIRE(shifted[3]!=poisson_precision::evaluate(model,low)[3]);
    REQUIRE((extended.tail(8).array()==original.tail(8).array()).all());
    REQUIRE(extended[3]==original[3]);
    REQUIRE((extended.head(3)-original.head(3)).norm()<1e-12*original.head(3).norm());
    // Arithmetic mode must not alter the physical analytic Jacobian.
    REQUIRE((assembler.assembleJacobian(x,bc)-originalJ).norm()==0.);
    // Verify the corrected charge response against the production derivative.
    // A 10 microvolt step resolves the weak charge derivative on top of the
    // dielectric residual; its central-difference truncation is O((h/Vt)^2).
    for(int j=0;j<12;++j) {
        auto xp=x,xm=x;xp[j]+=1e-5;xm[j]-=1e-5;
        const VectorXd fd=(assembler.residual(xp,bc)-assembler.residual(xm,bc))/(2e-5);
        for(int i=0;i<3;++i) {
            CAPTURE(i,j,material.intrinsic_reference_potential_V,fd[i],originalJ.coeff(i,j));
            REQUIRE(fd[i]==Catch::Approx(originalJ.coeff(i,j)).epsilon(2e-6).margin(1e-24));
        }
    }
}

TEST_CASE("per-cell dielectric diagnostic matches triangle energy and full Jacobian",
          "[poisson][permittivity][jacobian]")
{
    const auto mesh = makeHorizontalInterfaceMesh(true, GENERATE(false, true));
    MaterialDatabase db;
    const auto doping = DopingModel::fromMeshAndRegions(mesh,
        {{"channel", 1.0e20, 0.0}, {"gate_oxide", 0.0, 0.0}});
    DDScalingSpec spec;
    spec.poissonPerCellPermittivity = true;
    spec.poissonChargeVolumePolicy = GENERATE(PoissonChargeVolumePolicy::Global,
        PoissonChargeVolumePolicy::MaterialLocalCircumcentric);
    auto recombination = recombinationModelConfig({"srh","auger"},1e-6,1e-6);
    recombination.volumePolicy = GENERATE("legacy_node","material_local_circumcentric");
    CoupledDDAssembler assembler(mesh, db, doping, constants::Vt_300,
        MobilityModelConfig{}, recombination,
        BandgapNarrowingConfig{}, ImpactIonizationModelConfig{}, {}, {}, spec);
    CoupledDDState state{VectorXd::Constant(4,0.012),
        VectorXd::Constant(4,0.003),VectorXd::Constant(4,-0.002)};
    state.psi(1)=0.015;state.psi(2)=0.007;
    // Nonzero currents avoid a relative FD check of cancelling equilibrium SG terms.
    state.phin(1)=0.008;state.phin(2)=-0.003;
    state.phip(1)=0.004;state.phip(2)=-0.006;
    auto x=assembler.pack(state);
    CoupledDDBoundaryConditions bc;bc.psi[3]=0.012;
    const auto J=assembler.assembleJacobian(x,bc,false);
    const auto terms=assembler.poissonTermDiagnostics(x,bc);
    Eigen::MatrixXd energy=Eigen::MatrixXd::Zero(4,4);
    for (const auto& cell:mesh.cells()) {
        const auto& a=mesh.getNode(cell.node_ids[0]);
        const auto& b=mesh.getNode(cell.node_ids[1]);
        const auto& c=mesh.getNode(cell.node_ids[2]);
        const Real det=(b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y);
        Eigen::Matrix<double,3,2> grad;
        grad << b.y-c.y,c.x-b.x,c.y-a.y,a.x-c.x,a.y-b.y,b.x-a.x;
        grad/=det;
        const auto K=(db.getMaterial(mesh.getRegion(cell.region_id).material).eps_r*
            constants::eps0*std::abs(det)/2*(grad*grad.transpose())).eval();
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            energy(cell.node_ids[i],cell.node_ids[j])+=K(i,j);
    }
    // Independent dense products can round the two triangular halves differently.
    REQUIRE((energy-energy.transpose()).norm() <=
            16*std::numeric_limits<Real>::epsilon()*energy.norm());
    REQUIRE((energy*VectorXd::Ones(4)).norm()<1e-25);
    REQUIRE(state.psi.dot(energy*state.psi)>0);
    const auto flux=(energy*state.psi).eval();
    for(int i=0;i<4;++i) {
        REQUIRE(terms[i].dielectricFlux==Catch::Approx(flux(i)).margin(1e-25));
        REQUIRE(std::abs(terms[i].dielectricFlux)<=terms[i].dielectricOperandSum);
        REQUIRE(std::abs(terms[i].closureError)<1e-23);
        if(i==3)continue;
        for(int j=0;j<4;++j)if(i!=j)
            REQUIRE(J.coeff(i,j)==Catch::Approx(energy(i,j)).margin(1e-25));
    }
    for(int k=0;k<x.size();++k) {
        auto xp=x,xm=x;const Real h=1e-7;xp(k)+=h;xm(k)-=h;
        const auto fd=((assembler.residual(xp,bc)-assembler.residual(xm,bc))/(2*h)).eval();
        for(int i=0;i<x.size();++i) {
            CAPTURE(i,k,fd(i),J.coeff(i,k));
            REQUIRE(fd(i)==Catch::Approx(J.coeff(i,k)).epsilon(2e-6).margin(1e-18));
        }
    }
    REQUIRE(J.coeff(3,3)==1.0);
    REQUIRE_THROWS(DDAssembler(mesh,db,doping,constants::Vt_300,
        MobilityModelConfig{},recombinationModelConfig({"srh"},1e-6,1e-6),
        BandgapNarrowingConfig{},ImpactIonizationModelConfig{}, {}, {}, spec));
}

TEST_CASE("per-cell dielectric policy is opt-in and rejects unknown names",
          "[poisson][permittivity][config]")
{
    REQUIRE_FALSE(parsePoissonPerCellPermittivity(nlohmann::json::object()));
    REQUIRE_FALSE(parsePoissonPerCellPermittivity({{"discretization",{{"poisson_permittivity_coupling","edge_mean"}}}}));
    REQUIRE(parsePoissonPerCellPermittivity({{"discretization",{{"poisson_permittivity_coupling","per_cell"}}}}));
    REQUIRE_THROWS(parsePoissonPerCellPermittivity({{"discretization",{{"poisson_permittivity_coupling","unknown"}}}}));
}

TEST_CASE("circumcentric material volumes conserve area and reject nonpositive nodes",
          "[poisson][charge_volume][circumcentric]")
{
    const auto mesh=makeSiliconTriangleMesh();
    const MaterialDatabase db;
    const auto materials=detail::buildCellMaterials(mesh,db,constants::T0);
    const auto policy=PoissonChargeVolumePolicy::MaterialLocalCircumcentric;
    const auto v=detail::computePoissonChargeVolumes(mesh,materials,policy);
    REQUIRE(v[0]==Catch::Approx(0.25e-12).epsilon(1e-14));
    REQUIRE(v[1]==Catch::Approx(0.125e-12).epsilon(1e-14));
    REQUIRE(v[2]==Catch::Approx(0.125e-12).epsilon(1e-14));
    REQUIRE(v[0]!=mesh.getNode(0).volume); // No all-semiconductor legacy shortcut.
    REQUIRE(v[0]+v[1]+v[2]==Catch::Approx(0.5e-12).epsilon(1e-14));
    REQUIRE(detail::recombinationNodeVolumes(mesh,materials,"material_local_circumcentric")==v);
    const auto bad=makeSiliconTriangleMesh(true);
    const auto shares=BoxGeometryBuilder::circumcentricVertexVolumes(bad,bad.getCell(0));
    REQUIRE(*std::min_element(shares.begin(),shares.end())<0);
    REQUIRE_THROWS(detail::computePoissonChargeVolumes(bad,materials,policy));
    REQUIRE(parsePoissonChargeVolumePolicy({{"discretization",{{"poisson_charge_volume_policy","material_local_circumcentric"}}}})==policy);
}

TEST_CASE("material-local Poisson charge volume excludes oxide cell area",
          "[poisson][charge_volume]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const auto materials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto global = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::Global);
    const auto local = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::MaterialLocalBarycentric);

    REQUIRE(local.at(0) == Catch::Approx(0.5 * global.at(0)).epsilon(1.0e-14));
    REQUIRE(local.at(1) == Catch::Approx(0.5 * global.at(1)).epsilon(1.0e-14));
    REQUIRE(local.at(2) == global.at(2));
    REQUIRE(local.at(3) == 0.0);
}

TEST_CASE("material-local Poisson charge volume preserves silicon-only mesh",
          "[poisson][charge_volume]")
{
    const DeviceMesh mesh = makeSiliconTriangleMesh();
    const MaterialDatabase matdb;
    const auto materials = detail::buildCellMaterials(mesh, matdb, constants::T0);
    const auto global = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::Global);
    const auto local = detail::computePoissonChargeVolumes(
        mesh, materials, PoissonChargeVolumePolicy::MaterialLocalBarycentric);
    REQUIRE(local == global);
}

TEST_CASE("material-local Poisson charge volume changes residual and Jacobian together",
          "[poisson][charge_volume][jacobian]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});

    auto build = [&](PoissonChargeVolumePolicy policy) {
        DDScalingSpec scaling;
        scaling.poissonChargeVolumePolicy = policy;
        return CoupledDDAssembler(
            mesh, matdb, doping, constants::Vt_300,
            MobilityModelConfig{},
            recombinationModelConfig({"srh"}, 1.0e-6, 1.0e-6),
            BandgapNarrowingConfig{}, ImpactIonizationModelConfig{},
            {}, {}, scaling);
    };

    CoupledDDAssembler baseline = build(PoissonChargeVolumePolicy::Global);
    CoupledDDAssembler candidate = build(
        PoissonChargeVolumePolicy::MaterialLocalBarycentric);
    CoupledDDState state{
        VectorXd::Zero(4), VectorXd::Constant(4, 0.01),
        VectorXd::Constant(4, -0.02)};
    const VectorXd baselineState = baseline.pack(state);
    const VectorXd candidateState = candidate.pack(state);
    const CoupledDDBoundaryConditions bcs;

    const VectorXd baselineResidual = baseline.residual(baselineState, bcs);
    const VectorXd candidateResidual = candidate.residual(candidateState, bcs);
    REQUIRE((candidateResidual.segment(4, 8) -
             baselineResidual.segment(4, 8)).norm() == 0.0);
    const auto baselineTerms = baseline.poissonTermDiagnostics(baselineState, bcs);
    const auto candidateTerms = candidate.poissonTermDiagnostics(candidateState, bcs);
    REQUIRE(candidateTerms.at(0).dopingCharge ==
            Catch::Approx(0.5 * baselineTerms.at(0).dopingCharge).epsilon(1.0e-12));
    REQUIRE(candidateTerms.at(0).electronCharge ==
            Catch::Approx(0.5 * baselineTerms.at(0).electronCharge).epsilon(1.0e-12));

    const SparseMatrixd baselineJacobian = baseline.assembleJacobian(baselineState, bcs);
    const SparseMatrixd candidateJacobian = candidate.assembleJacobian(candidateState, bcs);
    REQUIRE(candidateJacobian.coeff(0, 4) ==
            Catch::Approx(0.5 * baselineJacobian.coeff(0, 4)).epsilon(1.0e-12));
    REQUIRE(candidateJacobian.coeff(0, 8) ==
            Catch::Approx(0.5 * baselineJacobian.coeff(0, 8)).epsilon(1.0e-12));
}

TEST_CASE("Poisson charge volume policy parsing is explicit and barycentric-only",
          "[poisson][charge_volume][config]")
{
    REQUIRE(parsePoissonChargeVolumePolicy(nlohmann::json::object()) ==
            PoissonChargeVolumePolicy::Global);
    REQUIRE(parsePoissonChargeVolumePolicy(nlohmann::json{
                {"discretization", {{"poisson_charge_volume_policy", "material_local"}}}}) ==
            PoissonChargeVolumePolicy::MaterialLocalBarycentric);
    REQUIRE_THROWS_WITH(
        parsePoissonChargeVolumePolicy(nlohmann::json{
            {"mesh_geometry", {{"node_volume_policy", "mixed_voronoi"}}},
            {"discretization", {{"poisson_charge_volume_policy", "material_local"}}}}),
        Catch::Matchers::ContainsSubstring("qualified only"));
}

TEST_CASE("Mixed material-local Poisson preserves volume and continuity with full charge derivatives",
          "[poisson][charge_volume][jacobian]")
{
    auto mesh=makeHorizontalInterfaceMesh();
    BoxGeometryBuilder::Options geometry;
    geometry.nodeVolumePolicy=BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
    mesh.buildBoxGeometry(geometry);
    MaterialDatabase materials;
    const auto cellMaterials=detail::buildCellMaterials(mesh,materials,300.);
    const auto local=detail::computePoissonChargeVolumes(mesh,cellMaterials,PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi);
    const auto global=detail::computeNodeVolumes(mesh);
    CHECK(local[0]==Catch::Approx(.5*global[0]).epsilon(1e-14));
    CHECK(local[1]==Catch::Approx(.5*global[1]).epsilon(1e-14));
    CHECK(local[2]==global[2]);CHECK(local[3]==0.);
    CHECK(local[0]+local[1]+local[2]==Catch::Approx(.5e-12).epsilon(1e-14));
    const auto doping=DopingModel::fromMeshAndRegions(mesh,{{"channel",1e21,0.}});
    DDScalingSpec scaling;scaling.poissonChargeVolumePolicy=PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi;
    CoupledDDAssembler candidate(mesh,materials,doping,constants::Vt_300,{}, {}, {}, {}, {}, {},scaling);
    scaling.poissonChargeVolumePolicy=PoissonChargeVolumePolicy::Global;
    CoupledDDAssembler baseline(mesh,materials,doping,constants::Vt_300,{}, {}, {}, {}, {}, {},scaling);
    CoupledDDState state{VectorXd::Constant(4,.03),VectorXd::Constant(4,.01),VectorXd::Constant(4,-.02)};
    const auto x=candidate.pack(state);CoupledDDBoundaryConditions bc;
    CHECK(candidate.residual(x,bc).tail(8)==baseline.residual(x,bc).tail(8));
    const auto jac=candidate.assembleJacobian(x,bc);
    for(int col=0;col<12;++col){
        auto a=x,b=x;a[col]+=1e-6;b[col]-=1e-6;
        const VectorXd fd=(candidate.residual(a,bc)-candidate.residual(b,bc))/(2e-6);
        const VectorXd analytic=VectorXd(jac.col(col));
        CHECK((analytic-fd).norm()/std::max({analytic.norm(),fd.norm(),1e-30})<1e-6);
    }
    const nlohmann::json cfg{{"mesh_geometry",{{"node_volume_policy","mixed_voronoi"}}},
        {"discretization",{{"poisson_charge_volume_policy","material_local_mixed_voronoi"}}}};
    CHECK(parsePoissonChargeVolumePolicy(cfg)==PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi);
    auto bad=cfg;bad["mesh_geometry"]["node_volume_policy"]="barycentric";
    CHECK_THROWS(parsePoissonChargeVolumePolicy(bad));
    auto single=makeSiliconTriangleMesh();single.buildBoxGeometry(geometry);
    CHECK(detail::computePoissonChargeVolumes(single,detail::buildCellMaterials(single,materials,300.),
        PoissonChargeVolumePolicy::MaterialLocalMixedVoronoi)==detail::computeNodeVolumes(single));
}

TEST_CASE("Poisson carrier assembly reports a non-finite carrier state",
          "[poisson][diagnostic]")
{
    const DeviceMesh mesh = makeHorizontalInterfaceMesh();
    const MaterialDatabase matdb;
    const DopingModel doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    const VectorXd psi = VectorXd::Zero(4);
    const VectorXd electrons = VectorXd::Constant(4, 1.0e16);
    VectorXd holes = VectorXd::Constant(4, 1.0e16);
    holes(0) = std::numeric_limits<Real>::infinity();

    RecombinationModelConfig recombination;
    recombination.mechanisms = {"none"};
    DDAssembler assembler(mesh, matdb, doping, constants::Vt_300,
                          MobilityModelConfig{}, recombination,
                          BandgapNarrowingConfig{},
                          ImpactIonizationModelConfig{}, {}, {});

    REQUIRE_THROWS_WITH(
        assembler.assemblePoissonWithCarriers(electrons, holes, psi),
        Catch::Matchers::ContainsSubstring("non-finite carrier diagonal at node 0"));
}

TEST_CASE("Poisson diagnostics preserve explicit interface geometry after policy integration",
          "[poisson][charge_volume][production_geometry]")
{
    DeviceMesh mesh = makeHorizontalInterfaceMesh();
    BoxGeometryBuilder::Options geometry;
    geometry.poissonPermittivityPolicy = BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
    mesh.buildBoxGeometry(geometry);
    const MaterialDatabase materials;
    const auto doping = DopingModel::fromMeshAndRegions(
        mesh, {{"channel", 1.0e21, 0.0}, {"gate_oxide", 0.0, 0.0}});
    for (const auto policy : {PoissonChargeVolumePolicy::Global,
                             PoissonChargeVolumePolicy::MaterialLocalBarycentric}) {
        for (const auto volume : {"inherit", "signed_transport"}) {
            DDScalingSpec scaling;
            scaling.poissonChargeVolumePolicy = policy;
            scaling.regionResolvedInterfaceAssembly.poissonChargeNodeVolume = volume;
            CoupledDDAssembler assembler(mesh, materials, doping, constants::Vt_300,
                MobilityModelConfig{}, recombinationModelConfig({"none"}),
                {}, {}, {}, {}, scaling);
            CoupledDDState state{VectorXd::LinSpaced(4, -0.03, 0.04),
                VectorXd::Constant(4, 0.01), VectorXd::Constant(4, -0.02)};
            const auto x = assembler.pack(state);
            const auto residual = assembler.residual(x, {});
            const auto terms = assembler.poissonTermDiagnostics(x, {});
            for (int i = 0; i < 4; ++i) {
                const auto& term = terms.at(i);
                const Real reconstructed = term.dielectricFlux + term.electronCharge +
                    term.holeCharge + term.dopingCharge + term.fixedInterfaceCharge;
                REQUIRE(reconstructed == Catch::Approx(residual(i)).epsilon(1e-12).margin(1e-30));
            }
        }
    }
}
