#include "vela/equation/LatticeHeatAssembler.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Eigen/SparseLU>
#include <limits>

using namespace vela;
using Catch::Approx;
namespace {
DeviceMesh rectangle(unsigned nx=4, unsigned ny=2, Real unit=1.) {
    DeviceMesh mesh;
    mesh.addRegion({0,"left","left",{}}); mesh.addRegion({1,"right","right",{}});
    for (unsigned y=0;y<=ny;++y) for (unsigned x=0;x<=nx;++x)
        mesh.addNode({mesh.numNodes(),2.*x/nx/unit,1.*y/ny/unit,0});
    for (unsigned y=0;y<ny;++y) for (unsigned x=0;x<nx;++x) {
        Index a=y*(nx+1)+x,b=a+1,c=a+nx+1,d=c+1,r=x<nx/2?0:1;
        mesh.addCell({mesh.numCells(),CellType::Tri3,r,{a,b,d}});
        mesh.addCell({mesh.numCells(),CellType::Tri3,r,{a,d,c}});
    }
    mesh.buildEdges(); return mesh;
}
std::vector<LatticeThermodeEdge> boundaries(unsigned nx=4,unsigned ny=2,Real left=310,Real right=300) {
    std::vector<LatticeThermodeEdge> out;
    for (unsigned y=0;y<ny;++y) {
        out.push_back({{y*(nx+1),(y+1)*(nx+1)},left,5.});
        out.push_back({{y*(nx+1)+nx,(y+1)*(nx+1)+nx},right,20.});
    }
    return out;
}
VectorXd step(const LatticeHeatAssembly& a) {
    Eigen::SparseLU<SparseMatrixd> lu;
    lu.compute(a.jacobian_W_per_m_K); REQUIRE(lu.info()==Eigen::Success);
    VectorXd result=lu.solve(-a.residual_W_per_m);
    REQUIRE(lu.info()==Eigen::Success); return result;
}
LatticeConductivity constant(Real k) { LatticeConductivity law;law.constant_W_per_m_K=k;return law; }
LatticeConductivity silicon() {
    LatticeConductivity law;
    law.model=LatticeConductivity::Model::InverseQuadratic;
    law.numerator=100.;law.denominator={-.0393,.00155,1.82e-6};return law;
}
}

TEST_CASE("Lattice heat preserves zero-source equilibrium and thermal isolation", "[lattice_heat]") {
    auto mesh=rectangle();auto t=VectorXd::Constant(mesh.numNodes(),300.);
    auto q=VectorXd::Zero(mesh.numCells());
    LatticeHeatAssembler closed(mesh,1.,{{0,silicon()},{1,constant(1.4)}},{});
    auto a=closed.assemble(t,q);
    REQUIRE(a.residual_W_per_m.norm()==0.);
    REQUIRE((a.jacobian_W_per_m_K*VectorXd::Ones(t.size())).norm()<1e-12);
    LatticeHeatAssembler robin(mesh,1.,{{0,silicon()},{1,constant(1.4)}},boundaries(4,2,300.,300.));
    a=robin.assemble(t,q);
    REQUIRE(a.residual_W_per_m.norm()==0.);
    REQUIRE(a.integrated_source_W_per_m==0.);
    REQUIRE(a.outward_boundary_heat_W_per_m==0.);
}

TEST_CASE("Lattice heat matches series thermal resistance and interface flux", "[lattice_heat]") {
    auto mesh=rectangle();
    LatticeHeatAssembler solver(mesh,1.,{{0,constant(10.)},{1,constant(2.)}},boundaries());
    VectorXd t=VectorXd::Constant(mesh.numNodes(),300.);
    const VectorXd zero=VectorXd::Zero(mesh.numCells());
    t+=step(solver.assemble(t,zero));
    const Real flux=10./(1./5.+1./10.+1./2.+1./20.);
    for (const auto& n:mesh.nodes()) {
        Real expected=310.-flux/5.-flux*(std::min(n.x,1.)/10.+std::max(n.x-1.,0.)/2.);
        REQUIRE(t[n.id]==Approx(expected).margin(1e-11));
    }
    const auto a=solver.assemble(t,zero);
    REQUIRE(a.residual_W_per_m.norm()<1e-10);
    REQUIRE(std::abs(a.outward_boundary_heat_W_per_m)<1e-10);
    // Independent gradients on either side of the material interface.
    REQUIRE(10.*(t[1]-t[2])/.5==Approx(flux).margin(1e-10));
    REQUIRE(2.*(t[2]-t[3])/.5==Approx(flux).margin(1e-10));
}

TEST_CASE("Lattice heat source balances boundary power in SI and micrometre meshes", "[lattice_heat]") {
    VectorXd reference;
    for (Real unit:{1.,1e-6}) {
        auto mesh=rectangle(8,4,unit);
        LatticeHeatAssembler solver(mesh,unit,{{0,constant(10.)},{1,constant(2.)}},boundaries(8,4,300.,300.));
        VectorXd t=VectorXd::Constant(mesh.numNodes(),300.);
        VectorXd q=VectorXd::Constant(mesh.numCells(),100.);
        t+=step(solver.assemble(t,q));
        const auto a=solver.assemble(t,q);
        REQUIRE(t.minCoeff()>300.);
        REQUIRE(a.integrated_source_W_per_m==Approx(200.).margin(1e-11));
        REQUIRE(a.outward_boundary_heat_W_per_m==Approx(200.).margin(1e-9));
        REQUIRE(a.residual_W_per_m.norm()<1e-9);
        REQUIRE(solver.nodalAreas_m2().sum()==Approx(2.));
        if (reference.size()) REQUIRE((reference-t).norm()<1e-9); else reference=t;
    }
}

TEST_CASE("Lattice temperature-dependent conductivity has consistent Newton derivatives", "[lattice_heat]") {
    auto mesh=rectangle();
    LatticeHeatAssembler solver(mesh,1.,{{0,silicon()},{1,constant(1.4)}},boundaries());
    VectorXd t(mesh.numNodes());for (const auto& n:mesh.nodes())t[n.id]=320.+10.*n.x+3.*n.y;
    const VectorXd q=VectorXd::Constant(mesh.numCells(),3.);
    const auto a=solver.assemble(t,q);
    for (Eigen::Index j=0;j<t.size();++j) {
        VectorXd plus=t,minus=t;plus[j]+=1e-3;minus[j]-=1e-3;
        const VectorXd fd=(solver.assemble(plus,q).residual_W_per_m-solver.assemble(minus,q).residual_W_per_m)/.002;
        const VectorXd analytic=a.jacobian_W_per_m_K.col(j);
        REQUIRE((fd-analytic).norm()/analytic.norm()<1e-8);
    }
    REQUIRE(a.residual_W_per_m.sum()==Approx(a.outward_boundary_heat_W_per_m-a.integrated_source_W_per_m).margin(1e-9));
    REQUIRE(silicon().valueAndDerivative(300.).first==Approx(169.63528413910093).epsilon(1e-12));
}

TEST_CASE("Lattice volumetric heating converges to a one dimensional analytic profile", "[lattice_heat]") {
    Real coarseError=0.;
    for (unsigned nx : {8u,32u}) {
        // Constant heat, insulated x=0/top/bottom, Robin cooling at x=2.
        // Exact T(x)=300 + Q*L/h + Q*(L^2-x^2)/(2*k).
        auto mesh=rectangle(nx,nx/2);
        auto robin=boundaries(nx,nx/2,300.,300.);
        std::vector<LatticeThermodeEdge> right;
        for (Index i=1;i<robin.size();i+=2) right.push_back(robin[i]);
        LatticeHeatAssembler solver(mesh,1.,{{0,constant(10.)},{1,constant(10.)}},right);
        VectorXd t=VectorXd::Constant(mesh.numNodes(),300.);
        const VectorXd q=VectorXd::Constant(mesh.numCells(),100.);
        t+=step(solver.assemble(t,q));
        Real error=0.;
        for (const auto& node:mesh.nodes()) {
            const Real exact=310.+5.*(4.-node.x*node.x);
            error=std::max(error,std::abs(t[node.id]-exact));
        }
        if (nx==8) coarseError=error;
        else {
            REQUIRE(error<coarseError/8.);
            REQUIRE(error<.02);
        }
    }
}

TEST_CASE("Lattice heat rejects ambiguous boundaries and invalid material inputs", "[lattice_heat]") {
    auto mesh=rectangle();auto b=boundaries();
    REQUIRE_THROWS_AS(LatticeHeatAssembler(mesh,1.,{{0,constant(1.)}},b),std::invalid_argument);
    REQUIRE_THROWS_AS(LatticeHeatAssembler(mesh,0.,{{0,constant(1.)},{1,constant(1.)}},b),std::invalid_argument);
    b.push_back(b[0]);
    REQUIRE_THROWS_AS(LatticeHeatAssembler(mesh,1.,{{0,constant(1.)},{1,constant(1.)}},b),std::invalid_argument);
    b={{{1,6},300.,5.}}; // Interior edge.
    REQUIRE_THROWS_AS(LatticeHeatAssembler(mesh,1.,{{0,constant(1.)},{1,constant(1.)}},b),std::invalid_argument);
    REQUIRE_THROWS_AS(silicon().valueAndDerivative(-1.),std::invalid_argument);
    REQUIRE_THROWS_AS(silicon().valueAndDerivative(std::numeric_limits<Real>::quiet_NaN()),std::invalid_argument);
    REQUIRE_THROWS_AS(constant(-1.).valueAndDerivative(300.),std::invalid_argument);
}
