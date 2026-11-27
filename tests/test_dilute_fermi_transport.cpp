#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/detail/MobilityChainProduct.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>
#include <limits>

using namespace vela;

TEST_CASE("Scaled Fermi endpoint populations preserve the dilute limit", "[dilute_transport]")
{
    using Wide = boost::multiprecision::cpp_dec_float_100;
    for (double dos : {3.1046289324748378e19, 2.8e25}) {
        for (double eta : {-700.001, -730.999507433708, -745.91221927522724, -750.}) {
            const double ref = static_cast<double>(Wide(dos) * exp(Wide(eta)));
            const double p = scaledFermiDiracHalf(dos, eta);
            CAPTURE(dos, eta, p, ref);
            REQUIRE(p > 0.);
            CHECK(std::abs(p/ref - 1.) < 3e-13);
            CHECK(p == holeDensity(dos,dos,0.,eta,1.,CarrierStatisticsModel::FermiDirac));
            CHECK(p == electronDensity(dos,dos,0.,-eta,1.,CarrierStatisticsModel::FermiDirac));
        }
        for (double eta : {-700., -30., -8., 0., 20.})
            CHECK(scaledFermiDiracHalf(dos,eta) == dos*fermiDiracHalf(eta));
    }
    CHECK(scaledFermiDiracHalf(3.1e19,-1000.) == 0.); // No artificial floor.
}

TEST_CASE("Effective band SG retains dilute flux and exposes invalid operands", "[dilute_transport][sg]")
{
    const double e0=-745.91221927522724,e1=-730.999507433708,dos=3.1046289324748378e19;
    const double p0=scaledFermiDiracHalf(dos,e0),p1=scaledFermiDiracHalf(dos,e1);
    const double vt=.025851995266484913,q0=-18.339136030651026,q1=-17.952782563196067;
    const double f=sgHoleFermiDiracContinuityFlux(p0,p1,e0,e1,0.,q0,q1,vt,1., FermiSgFormulation::EffectivePotential);
    const long double drop=-(static_cast<long double>(q1)-q0)/vt;
    const long double argument=std::log(static_cast<long double>(p1)/p0)+drop;
    const long double b=argument/std::expm1(argument);
    const long double ref=b*p1*std::expm1(drop);
    REQUIRE(std::isfinite(f)); REQUIRE(f<0.);
    CHECK(std::abs(static_cast<long double>(f)/ref-1)<2e-12L);
    const double reverse=sgHoleFermiDiracContinuityFlux(p1,p0,e1,e0,0.,q1,q0,vt,1., FermiSgFormulation::EffectivePotential);
    CHECK(std::abs(reverse/f+1.)<2e-12);
    CHECK(std::isnan(sgHoleFermiDiracContinuityFlux(0.,p1,e0,e1,0.,q0,q1,vt,1., FermiSgFormulation::EffectivePotential)));
    CHECK(std::isnan(sgElectronFermiDiracContinuityFlux(0.,p1,e0,e1,0.,q0,q1,vt,1., FermiSgFormulation::EffectivePotential)));
    CHECK(std::isnan(sgHoleFermiDiracContinuityFlux(p0,p1,e0,e1,0.,q0,q1,vt,
        std::numeric_limits<double>::quiet_NaN(), FermiSgFormulation::EffectivePotential)));
    // Legitimate finite-input overflow retains the existing saturation policy.
    CHECK(sgHoleFermiDiracContinuityFlux(1e100,1e100,0.,0.,0.,0.,1.,1.,1e300, FermiSgFormulation::EffectivePotential)
          == -std::numeric_limits<double>::max());
}

TEST_CASE("Zero mobility derivatives do not hide invalid flux", "[dilute_transport]")
{
    using detail::mobilityChainProduct;
    const double largest=std::numeric_limits<double>::max();
    CHECK(mobilityChainProduct(-largest,2.,0.,.02585)==0.);
    CHECK(mobilityChainProduct(-largest,2.,1.,.02585,0.)==0.);
    CHECK(std::isnan(mobilityChainProduct(std::numeric_limits<double>::quiet_NaN(),2.,0.)));
    CHECK(std::isnan(mobilityChainProduct(-std::numeric_limits<double>::infinity(),2.,0.)));
    CHECK(std::isnan(mobilityChainProduct(1.,0.,0.)));
    for(double derivative : {-2.,1e-8,3.})
        CHECK(mobilityChainProduct(1.2,.25,derivative,.026,.9)
              == 1.2/.25/.026*derivative*.9);
}

TEST_CASE("Coupled PhuMob dilute Fermi residual and Jacobian remain finite", "[dilute_transport][element_box]")
{
    DeviceMesh mesh;
    const std::array<std::array<double,2>,4> xy{{{0,0},{1,0},{.5,.3},{.5,-.9}}};
    for(Index i=0;i<4;++i){Node n;n.id=i;n.x=xy[i][0]*1e-6;n.y=xy[i][1]*1e-6;mesh.addNode(n);}
    for(Index i=0;i<2;++i){Cell c;c.id=i;c.type=CellType::Tri3;c.region_id=0;
        c.node_ids=i==0?std::vector<Index>{0,1,2}:std::vector<Index>{1,0,3};mesh.addCell(c);}
    Region region;region.id=0;region.name="Si";region.material="Si";region.cell_ids={0,1};mesh.addRegion(region);
    mesh.buildEdges();
    BoxGeometryBuilder::Options geometry;
    geometry.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
    geometry.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::DelaunayTransfer;
    mesh.buildBoxGeometry(geometry);
    MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e26,1e24);
    auto mobility=mobilityModelConfig("phumob");mobility.edgeAveraging="element_box_phumob";
    mobility.dopingConcentrationBasis="total_impurity";
    DDScalingSpec scaling;scaling.regionResolvedInterfaceAssembly.transportEdgeGeometry="element_box";
    const CarrierStatisticsConfig stats{.model = "fermi_dirac", .transport = "effective_band_sg"};
    CoupledDDAssembler assembler(mesh,db,doping,constants::Vt_300,mobility,
        recombinationModelConfig({"none"}),{},{},{},{},scaling,{},stats);
    const auto mat=db.getMaterial("Si",300.);
    for(bool cache : {false,true}) for(bool hole : {true,false}) {
        assembler.setFermiNodeCache(cache);
        CoupledDDState state{VectorXd(4),VectorXd(4),VectorXd(4)};
        for(int i=0;i<4;++i) {
            state.psi(i)=.02*i;
            const double diluteEta=i==0?-745.91221927522724:(i==1?-730.999507433708:-85.-i);
            const double en=hole?-2.:diluteEta,ep=hole?diluteEta:-2.;
            state.phin(i)=state.psi(i)-constants::Vt_300*(en-std::log(mat.ni/ *mat.Nc_m3));
            state.phip(i)=state.psi(i)+constants::Vt_300*(ep-std::log(mat.ni/ *mat.Nv_m3));
        }
        const auto x=assembler.pack(state);
        const auto density=hole?assembler.holeDensity(x):assembler.electronDensity(x);
        REQUIRE(density(0)>0.);
        for(double perturbation : {0.,-1e-7,1e-7}) {
            auto trial=x;trial((hole?8:4)+1)+=perturbation;
            const auto residual=assembler.residual(trial,{});
            REQUIRE(residual.allFinite());
            const auto jacobian=assembler.assembleJacobian(trial,{});
            for(int col=0;col<jacobian.outerSize();++col)
                for(SparseMatrixd::InnerIterator it(jacobian,col);it;++it)
                    REQUIRE(std::isfinite(it.value()));
        }
    }
}
