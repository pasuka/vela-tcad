#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/experimental/BvMethodsTri3Avalanche.h"
#include "vela/physics/experimental/BvMethodsTri3Transport.h"
#include "vela/physics/experimental/BvMethodsPrimitiveState.h"
#include <cmath>
#include "vela/physics/experimental/PopulationRecombination.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace vela::experimental;
using Catch::Approx;

TEST_CASE("Tri3 SG retains diffusion and both endpoint derivatives beyond exponential range",
          "[av_contract][tri3_ad]")
{
    AvJet left(1e20),right(1e-303);left.derivative[0]=left.value;right.derivative[1]=right.value;
    AvJet delta(std::log(left.value)-std::log(right.value));delta.derivative[0]=1.;delta.derivative[1]=-1.;
    const auto argument=avLogRatio(right,left)+delta;
    REQUIRE(argument.value==Approx(0.).margin(1e-12));
    const auto forward=avSgFlux(AvJet(1),left,right,delta,argument);
    const auto reverse=avSgFlux(AvJet(1),right,left,-delta,-argument);
    REQUIRE(forward.value==Approx(left.value));
    REQUIRE(forward.derivative[0]/left.value==Approx(1.));
    REQUIRE(forward.derivative[1]/right.value==Approx(-1.));
    REQUIRE(reverse.value==-forward.value);
    for(int k=0;k<9;++k){
        REQUIRE(std::isfinite(forward.derivative[k]));
        REQUIRE(reverse.derivative[k]==-forward.derivative[k]);
    }
}

TEST_CASE("Tri3 dual quotient preserves finite derivatives across density scales",
          "[av_contract][tri3_ad]")
{
    for(double scale:{1e-200,1.,1e200}){
        const auto a=AvJet::variable(5.*scale,0);
        const auto b=AvJet::variable(2.*scale,1);
        const auto ratio=a/b;
        REQUIRE(ratio.value==Approx(2.5));
        REQUIRE(ratio.derivative[0]*scale==Approx(.5));
        REQUIRE(ratio.derivative[1]*scale==Approx(-1.25));
        for(int k=2;k<9;++k)REQUIRE(ratio.derivative[k]==0.);
        const auto unity=a/a;
        REQUIRE(unity.value==1.);
        for(double d:unity.derivative)REQUIRE(d==0.);
    }
}

TEST_CASE("Tri3 SG transport and AV remain homogeneous at extreme positive populations",
          "[av_contract][tri3_ad]")
{
    const AvTri3Geometry geometry{{{{0.,0.},{2.,0.},{.5,1.}}},{.25,.4,.35},0};
    auto evaluate=[&](double scale,bool electron){
        std::array<AvJet,3> density,qf;
        for(int k=0;k<3;++k){
            qf[k]=AvJet::variable(.01*k,3+k);
            density[k]=AvJet(scale*(2.+k))*vela::detail::dualExp(qf[k]/AvJet(.02585));
        }
        const auto transport=evaluateBvMethodsTri3Transport({.25,.4,.35},{.2,.7,.3},
            density,qf,{AvJet(50.),AvJet(80.),AvJet(100.)},.02585,electron);
        const auto av=evaluateBvMethodsTri3Avalanche(geometry,density,density,qf,qf,
            AvJet(50.),AvJet(80.),AvJet(5e4),AvJet(6e4),.02585,
            {7e5,1.2e6,7e5,1.2e6,4e5},{1.6e6,2e6,1.6e6,2e6,4e5});
        return std::array<AvJet,4>{transport.nodeFlux_per_cm_s[0],transport.nodeFlux_per_cm_s[1],
            transport.nodeFlux_per_cm_s[2],av.generation_cm3_s};
    };
    for(bool electron:{true,false}){
        const auto normal=evaluate(1.,electron);
        for(double scale:{1e-200,1e200}){
            const auto scaled=evaluate(scale,electron);
            for(int k=0;k<4;++k){
                REQUIRE(std::isfinite(scaled[k].value));
                REQUIRE(scaled[k].value/scale==Approx(normal[k].value).epsilon(1e-12));
                for(int j=0;j<9;++j){
                    REQUIRE(std::isfinite(scaled[k].derivative[j]));
                    REQUIRE(scaled[k].derivative[j]/scale==Approx(normal[k].derivative[j]).epsilon(1e-11).margin(1e-25));
                }
            }
        }
    }
}

TEST_CASE("Experimental transport preserves uniform-density Ohm law and paired edge flux",
          "[av_contract]")
{
    const std::array<double,3> areas{.5,.25,.25},coefficients{2.,.3,.7};
    const std::array<AvJet,3> n{AvJet(3e16),AvJet(3e16),AvJet(3e16)};
    const std::array<AvJet,3> qf{AvJet(0.),AvJet(.013),AvJet(-.007)};
    const std::array<AvJet,3> mu{AvJet(2.),AvJet(4.),AvJet(8.)};
    for(bool electron:{true,false})for(double factor:{1.,3.}){
        const auto result=evaluateBvMethodsTri3Transport(areas,coefficients,n,qf,mu,.02585,electron,
            {AvJet(factor),AvJet(factor),AvJet(factor)});
        REQUIRE(result.elementMobility_cm2_V_s.value==Approx(4.));
        for(int k=0;k<3;++k){
            const double expected=(electron?1.:-1.)*coefficients[k]*4.*3e16
                *(qf[(k+1)%3].value-qf[k].value);
            REQUIRE(result.edgeFlux_per_cm_s[k].value==Approx(expected));
        }
        const double sum=result.nodeFlux_per_cm_s[0].value+result.nodeFlux_per_cm_s[1].value+result.nodeFlux_per_cm_s[2].value;
        REQUIRE(sum==Approx(0.).margin(1.)); // floating summation of O(1e15) rates
    }
}

TEST_CASE("Experimental transport derivative respects quasi-Fermi gauge and conservation",
          "[av_contract]")
{
    const std::array<AvJet,3> n{AvJet(1e16),AvJet(2e16),AvJet(4e16)};
    std::array<AvJet,3> qf{AvJet::variable(.001,3),AvJet::variable(.003,4),AvJet::variable(-.002,5)};
    const auto result=evaluateBvMethodsTri3Transport({.5,.3,.2},{.2,.7,.3},n,qf,
        {AvJet(50.),AvJet(80.),AvJet(100.)},.02585,true);
    for(int k=0;k<3;++k){
        const auto& f=result.nodeFlux_per_cm_s[k];
        const double norm=std::abs(f.derivative[3])+std::abs(f.derivative[4])+std::abs(f.derivative[5]);
        REQUIRE(std::abs(f.derivative[3]+f.derivative[4]+f.derivative[5])<1e-14*norm);
    }
    for(int j=3;j<6;++j){
        double sum=0.,norm=0.;
        for(const auto& f:result.nodeFlux_per_cm_s){sum+=f.derivative[j];norm+=std::abs(f.derivative[j]);}
        REQUIRE(std::abs(sum)<1e-14*norm);
    }
}

TEST_CASE("Experimental AV norm keeps transverse and zero-current directional rules",
          "[av_contract]")
{
    auto x=AvJet::variable(3.,0),y=AvJet::variable(4.,1);
    const auto a=avNorm({x,y});
    REQUIRE(a.value==Approx(5.));
    REQUIRE(a.derivative[0]==Approx(.6));
    REQUIRE(a.derivative[1]==Approx(.8));
    REQUIRE(avNormDirectional({3.,4.},{-4.,3.})==Approx(0.).margin(1e-15));
    REQUIRE(avNormDirectional({0.,0.},{3.,4.})==Approx(5.));
    REQUIRE(avNormDirectional({0.,0.},{-3.,-4.})==Approx(5.));
    REQUIRE(avNorm({AvJet::variable(0.,0),AvJet::variable(0.,1)}).derivative[0]==0.);
}

TEST_CASE("Experimental AV SG recovers exact diffusion on a nonright triangle",
          "[av_contract]")
{
    const double vt=.02585;
    AvTri3Geometry g{{{{0.,0.},{2.,0.},{.5,1.}}},{.25,.4,.35},0};
    std::array<AvJet,3> n{AvJet(2),AvJet(4),AvJet(5)},fn{};
    for(int i=0;i<3;++i)fn[i]=AvJet(-vt*std::log(n[i].value));
    for(std::size_t representative=0;representative<3;++representative){
        g.representative=representative;
        auto j=avParticleCurrentWithoutMobility(g,n,fn,vt,true);
        REQUIRE(j[0].value==Approx(vt));
        REQUIRE(j[1].value==Approx(2.5*vt));
        std::array<AvJet,3> fp{};
        for(int i=0;i<3;++i)fp[i]=-fn[i];
        const auto hole=avParticleCurrentWithoutMobility(g,n,fp,vt,false);
        REQUIRE(hole[0].value==Approx(-vt));
        REQUIRE(hole[1].value==Approx(-2.5*vt));
    }
}

TEST_CASE("Experimental AV full local derivative follows its finite source function",
          "[av_contract]")
{
    const AvTri3Geometry g{{{{0.,0.},{1e-7,0.},{.3e-7,1e-7}}},
                          {2e-15,1.5e-15,1.5e-15},1};
    const double vt=.02585;
    const std::array<double,9> point{.3,.307,.306,.01,-.005,.008,-.02,.004,.015};
    const AvCoefficient coefficient{7e5,1.2e6,6e5,1.1e6,8e4};
    // An explicit smooth state law exercises density, mobility and field chains
    // together; this is not a fit to a native matrix or its finite differences.
    auto source=[&](const std::array<double,9>& values){
        std::array<AvJet,3> psi,fn,fp,n,p;
        for(int k=0;k<3;++k){
            psi[k]=AvJet::variable(values[k],k);
            fn[k]=AvJet::variable(values[3+k],3+k);
            fp[k]=AvJet::variable(values[6+k],6+k);
            n[k]=AvJet(1e16)*vela::detail::dualExp((psi[k]-fn[k]-AvJet(.3))/AvJet(vt));
            p[k]=AvJet(1e16)*vela::detail::dualExp((-psi[k]+fp[k]+AvJet(.3))/AvJet(vt));
        }
        const auto ex=(psi[1]-psi[0])/AvJet(1e-7);
        const auto ey=(psi[2]-psi[0])/AvJet(1e-7)-AvJet(.3)*ex;
        const auto field=avNorm({ex,ey});
        const auto mun=AvJet(100)/(AvJet(1)+field/AvJet(1e5));
        const auto mup=AvJet(60)/(AvJet(1)+field/AvJet(8e4));
        return evaluateBvMethodsTri3Avalanche(g,n,p,fn,fp,mun,mup,
                    field,field,vt,coefficient,coefficient);
    };
    const auto center=source(point);
    REQUIRE_FALSE(center.electronZeroNorm);
    REQUIRE_FALSE(center.holeZeroNorm);
    for(int j=0;j<9;++j){
        auto plus=point,minus=point;
        const double h=1e-7;
        plus[j]+=h;minus[j]-=h;
        const auto a=source(plus),b=source(minus);
        for(int k=0;k<3;++k){
            const double fd=(a.nodeGeneration_per_cm_s[k].value-b.nodeGeneration_per_cm_s[k].value)/(2*h);
            REQUIRE(center.nodeGeneration_per_cm_s[k].derivative[j]==Approx(fd).epsilon(2e-6));
        }
    }
}

TEST_CASE("Experimental AV box source preserves total area and paired generation",
          "[av_contract]")
{
    AvTri3Geometry g{{{{0.,0.},{1.,0.},{0.,1.}}},{.25,.125,.125},0};
    std::array<AvJet,3> n{AvJet(1),AvJet(1),AvJet(1)},p=n;
    std::array<AvJet,3> fn{AvJet(0),AvJet::variable(.01,4),AvJet(0)},fp{};
    AvCoefficient alpha{2.,1.,2.,1.,4e5};
    auto a=evaluateBvMethodsTri3Avalanche(g,n,p,fn,fp,AvJet(3),AvJet(2),AvJet(100),AvJet(100),.02585,alpha,alpha);
    const auto sum=a.nodeGeneration_per_cm_s[0]+a.nodeGeneration_per_cm_s[1]+a.nodeGeneration_per_cm_s[2];
    REQUIRE(a.generation_cm3_s.value>0.);
    REQUIRE(sum.value==Approx(.5*a.generation_cm3_s.value));
    REQUIRE(sum.derivative[4]==Approx(.5*a.generation_cm3_s.derivative[4]));
    REQUIRE(a.nodeGeneration_per_cm_s[0].value==Approx(2*a.nodeGeneration_per_cm_s[1].value));
    REQUIRE_FALSE(a.electronZeroNorm);
    REQUIRE(a.holeZeroNorm);
}

TEST_CASE("Experimental AV coefficient has no artificial low-field floor",
          "[av_contract]")
{
    AvCoefficient p{703000.,1231000.,703000.,1231000.,400000.};
    REQUIRE(avCoefficient(AvJet::variable(1.,0),p).value==0.);
    REQUIRE(avCoefficient(AvJet::variable(1.,0),p).derivative[0]==0.);
    const double field=53000.;
    const auto a=avCoefficient(AvJet::variable(field,0),p);
    REQUIRE(a.derivative[0]==Approx(a.value*1231000./(field*field)));
}


#include "vela/core/PhysicalConstants.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/solver/NewtonSolver.h"
#include <nlohmann/json.hpp>
namespace {
struct ProfileFixture {
    vela::DeviceMesh mesh;
    vela::MaterialDatabase materials;
    vela::DopingModel doping{6};
    vela::MobilityModelConfig mobility=vela::mobilityModelConfig("masetti_field_lombardi");
    vela::ImpactIonizationModelConfig avalanche;
    vela::DDScalingSpec scaling;
    vela::CarrierStatisticsConfig statistics{"fermi_dirac_joyce_dixon"};
    ProfileFixture() {
        using namespace vela;
        for(const auto& n:std::vector<Node>{{0,0.,0.},{1,1e-8,0.},{2,0.,1e-8},{3,1e-8,1e-8},{4,2e-8,0.},{5,2e-8,1e-8}})mesh.addNode(n);
        mesh.addCell({0,CellType::Tri3,0,{0,1,2}});mesh.addCell({1,CellType::Tri3,0,{1,3,2}});
        mesh.addCell({2,CellType::Tri3,1,{1,4,3}});mesh.addCell({3,CellType::Tri3,1,{4,5,3}});
        mesh.addRegion({0,"silicon","Si",{0,1}});mesh.addRegion({1,"oxide","SiO2",{2,3}});
        mesh.addContact({0,"source",0,{0}});mesh.addContact({1,"drain",0,{2}});mesh.buildEdges();
        BoxGeometryBuilder::Options opts;opts.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
        opts.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::DelaunayTransfer;mesh.buildBoxGeometry(opts);
        for(int i=0;i<4;++i)doping.setNodeDoping(i,1e23,2e22);
        mobility.dopingConcentrationBasis="total_impurity";mobility.highFieldDrivingForce="quasi_fermi_gradient";
        avalanche=impactIonizationModelConfig("van_overstraeten");avalanche.sourceJacobianMode="local_ad";avalanche.drivingForce="quasi_fermi_gradient";avalanche.generation="current_density";
        scaling.coupledOperator="bvmethods_tri3_qf";
    }
    std::unique_ptr<vela::CoupledDDAssembler> assembler() {
        return std::make_unique<vela::CoupledDDAssembler>(mesh,materials,doping,vela::constants::Vt_300,mobility,
            vela::RecombinationModelConfig{},vela::BandgapNarrowingConfig{},avalanche,
            std::vector<vela::RegionFixedChargeSpec>{},std::vector<vela::InterfaceSheetChargeSpec>{},scaling,
            vela::CarrierDiagonalFloorRegularizationConfig{},statistics);
    }
    vela::VectorXd point() {
        vela::VectorXd x(18);x<<.05,.25,.07,.28,.25,.28,
            .001,.006,.002,.007,0.,0.,-.002,.003,-.001,.004,0.,0.;return x;
    }
};
}
TEST_CASE("BVMethods profile residual Jacobian and ledger share the same autonomous state", "[bvmethods_profile]") {
    using namespace vela;ProfileFixture f;auto a=f.assembler();auto x=f.point();auto F=a->residual(x,{});auto J=a->assembleJacobian(x,{});
    const auto rows=a->carrierContinuityEquationTermDiagnostics(x,{});
    for(int i=0;i<6;++i){
        REQUIRE(rows[i].electronResidual==Approx(F[6+i]).margin(1e-30));
        REQUIRE(rows[i].holeResidual==Approx(F[12+i]).margin(1e-30));
        REQUIRE(rows[i].electronResidual==Approx(rows[i].electronFlux+rows[i].electronRecombination+rows[i].electronImpact+rows[i].electronGauge).margin(1e-25));
        REQUIRE(rows[i].electronImpact==rows[i].holeImpact);
    }
    REQUIRE(a->impactIonizationSourceResidual(x,{}).norm()>0.);
    VectorXd v(18);for(int i=0;i<18;++i)v[i]=std::sin(i*.73)*.3;
    auto check=[&](const CoupledDDBoundaryConditions& bc){
        auto M=a->assembleJacobian(x,bc);VectorXd an=M*v;
        VectorXd fd=(a->residual(x+1e-7*v,bc)-a->residual(x-1e-7*v,bc))/(2e-7);
        for(int block=0;block<3;++block){const auto av=an.segment(6*block,6),fv=fd.segment(6*block,6);
            REQUIRE((av-fv).norm()/std::max(av.norm(),fv.norm())<2e-6);
        }
    };check({});CoupledDDBoundaryConditions bc;bc.psi[0]=.05;bc.phin[0]=.001;bc.phip[0]=-.002;check(bc);
    REQUIRE(a->residual(x,bc)[6]==Approx(0.).margin(1e-20));
    REQUIRE(a->assembleJacobian(x,bc).coeff(6,6)==1.);
    const auto bounded=a->carrierContinuityEquationTermDiagnostics(x,bc);REQUIRE_FALSE(bounded[0].electronContinuityActive);
}
TEST_CASE("BVMethods profile preserves equilibrium currents and rejects unsupported combinations", "[bvmethods_profile]") {
    using namespace vela;ProfileFixture f;auto a=f.assembler();auto x=f.point();x.segment(6,12).setZero();
    const auto rows=a->carrierContinuityEquationTermDiagnostics(x,{});
    for(int i=0;i<4;++i){REQUIRE(rows[i].electronFlux==0.);REQUIRE(rows[i].holeFlux==0.);REQUIRE(rows[i].electronImpact==0.);REQUIRE(rows[i].electronRecombination==0.);}
    f.statistics.model="fermi_dirac";auto unsupported=f.assembler();REQUIRE_THROWS_AS(unsupported->residual(x,{}),std::invalid_argument);
    f.statistics.model="fermi_dirac_joyce_dixon";f.mobility.dopingConcentrationBasis="net_doping";unsupported=f.assembler();REQUIRE_THROWS_AS(unsupported->residual(x,{}),std::invalid_argument);
    const auto legacy=newtonConfigFromJson(nlohmann::json::object());REQUIRE(legacy.coupledOperator=="legacy");
    REQUIRE_THROWS_AS(newtonConfigFromJson(nlohmann::json{{"coupled_operator","unknown"}}),std::invalid_argument);
}

#include "vela/post/ContactCurrent.h"
#include "vela/io/DDSolutionState.h"
#include "vela/simulation/CoupledLoadLine.h"
#include "vela/equation/BvMethodsTri3Operator.h"
#include "vela/solver/experimental/ExactValueBits.h"
TEST_CASE("BVMethods frozen mobility tangent preserves live transport and conservation", "[bvmethods_profile]") {
    using namespace vela;ProfileFixture f;auto a=f.assembler();const auto x=f.point();
    detail::BvMethodsQfState s{x.head(6),x.segment(6,6),x.tail(6),a->electronDensity(x),a->holeDensity(x)};
    std::vector<double> nc(6),nv(6);
    for(const auto& t:a->poissonTermDiagnostics(x,{})){nc[t.nodeId]=t.electronDensityOfStates;nv[t.nodeId]=t.holeDensityOfStates;}
    detail::BvMethodsTri3Context context{f.mesh,f.materials,f.doping,
        {f.mobility,f.avalanche,{},CarrierStatisticsModel::FermiDiracJoyceDixon}};
    detail::BvMethodsTri3Operator element(context,s,nc,nv,1.,1.,constants::Vt_300);
    const auto live=element.transport(s,true,true,true);
    const auto frozen=element.transport(s,true,true,true,true);
    REQUIRE((live.source.array()==frozen.source.array()).all());
    REQUIRE((live.absoluteSource.array()==frozen.absoluteSource.array()).all());
    const SparseMatrixd removed=live.jacobian-frozen.jacobian;
    REQUIRE(removed.norm()>0.);
    // Dropping mobility derivatives must retain edge-pair conservation and
    // invariance under a common electrostatic/QF reference shift.
    for(const auto* matrix:{&live.jacobian,&frozen.jacobian}) {
        for(int col=0;col<18;++col)for(int block:{1,2}){
            double sum=0.,magnitude=0.;for(int i=0;i<6;++i){double v=matrix->coeff(block*6+i,col);sum+=v;magnitude+=std::abs(v);}
            REQUIRE(std::abs(sum)<=1e-13*magnitude+1e-30);
        }
        const VectorXd gauge=(*matrix)*VectorXd::Ones(18);
        REQUIRE(gauge.norm()<1e-12*matrix->norm());
    }
    // A new state must recompute the live mobility even in the frozen-tangent
    // arm; this control freezes a derivative, never a value between iterations.
    s.phin[1]+=.002;s.n=a->electronDensity((VectorXd(18)<<s.psi,s.phin,s.phip).finished());
    const auto changed=element.transport(s,true,true,true,true);
    const auto changedLive=element.transport(s,true,true,true);
    REQUIRE((changed.source.array()==changedLive.source.array()).all());
    REQUIRE((changed.source-live.source).norm()>0.);
}
TEST_CASE("BVMethods profile terminal extraction follows its element residual", "[bvmethods_profile]") {
    using namespace vela;ProfileFixture f;auto a=f.assembler();auto x=f.point();DDSolution state;
    state.psi=x.head(6);state.phin=x.segment(6,6);state.phip=x.tail(6);state.n=a->electronDensity(x);state.p=a->holeDensity(x);
    ContactCurrent port(f.mesh,f.materials,f.doping,f.mobility,300.,f.scaling,{},f.statistics);
    const auto current=port.compute(state,"drain");const auto rows=a->carrierContinuityEquationTermDiagnostics(x,{});
    const double expected=constants::q*(rows[2].holeFlux-rows[2].electronFlux);
    REQUIRE(current.totalCurrent==Approx(expected).epsilon(2e-13));
    double sum=0.,absolute=0.;for(int i=0;i<4;++i){sum+=rows[i].holeFlux-rows[i].electronFlux;absolute+=std::abs(rows[i].holeFlux)+std::abs(rows[i].electronFlux);}
    REQUIRE(std::abs(sum)<1e-14*absolute);
    VectorXd invalid=x;invalid.segment(6,6).setConstant(1000.);REQUIRE_FALSE(a->residual(invalid,{}).allFinite());
    a->setQuasiFermiReferences(1.,1.);REQUIRE_THROWS_AS(a->residual(x,{}),std::invalid_argument);
}

TEST_CASE("BVMethods profile evaluation reuse preserves state and boundary sensitivity", "[bvmethods_profile]") {
    using namespace vela;using experimental::exactValueBits;
    ProfileFixture fixture;auto fresh=fixture.assembler(),reused=fixture.assembler();
    reused->setBvMethodsEvaluationReuse(true);
    auto x=fixture.point();CoupledDDBoundaryConditions bc;
    bc.psi[0]=.12;bc.phin[0]=.02;bc.phip[0]=-.03;
    for(int trial=0;trial<4;++trial) {
        if(trial==1)x[7]+=1e-13;
        if(trial==2)x[1]=std::nextafter(x[1],1.);
        if(trial==3){x.segment(6,12).setZero();bc.psi[0]=.2;}
        for(const auto& boundary:{CoupledDDBoundaryConditions{},bc}) {
            const auto expected=fresh->residual(x,boundary);
            REQUIRE(exactValueBits(reused->residual(x,boundary))==exactValueBits(expected));
            REQUIRE(exactValueBits(reused->residual(x,boundary))==exactValueBits(expected));
            REQUIRE(exactValueBits(reused->assembleJacobian(x,boundary))==exactValueBits(fresh->assembleJacobian(x,boundary)));
            REQUIRE(exactValueBits(reused->carrierContinuityEquationTermDiagnostics(x,boundary))==
                    exactValueBits(fresh->carrierContinuityEquationTermDiagnostics(x,boundary)));
            REQUIRE(exactValueBits(reused->bvMethodsTransportResidual(x))==exactValueBits(fresh->bvMethodsTransportResidual(x)));
            REQUIRE(exactValueBits(reused->bvMethodsTransportJacobian(x))==exactValueBits(fresh->bvMethodsTransportJacobian(x)));
        }
    }
    fresh->setQuasiFermiReferences(-0.,-0.);reused->setQuasiFermiReferences(-0.,-0.);
    REQUIRE(exactValueBits(reused->residual(x,bc))==exactValueBits(fresh->residual(x,bc)));
    fresh->setQuasiFermiReferences(0.,0.);reused->setQuasiFermiReferences(0.,0.);
    auto invalid=x;invalid.segment(6,6).setConstant(1000.);
    REQUIRE_FALSE(reused->residual(invalid,{}).allFinite());
    REQUIRE(exactValueBits(reused->residual(x,bc))==exactValueBits(fresh->residual(x,bc)));
    reused->setQuasiFermiReferences(1.,1.);
    REQUIRE_THROWS_AS(reused->residual(x,bc),std::invalid_argument);
    reused->setQuasiFermiReferences(0.,0.);reused->setBvMethodsEvaluationReuse(false);
    REQUIRE(exactValueBits(reused->residual(x,bc))==exactValueBits(fresh->residual(x,bc)));
}

TEST_CASE("BVMethods profile shared ports and boundary derivative preserve the full operator", "[bvmethods_profile]") {
    using namespace vela;using experimental::exactValueBits;
    ProfileFixture f;NewtonConfig cfg;cfg.coupledOperator="bvmethods_tri3_qf";
    cfg.carrierStatistics=f.statistics;cfg.mobility=f.mobility;cfg.impactIonization=f.avalanche;
    NewtonSolver baseline(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",0.}},cfg);
    cfg.diagnosticBvMethodsEvaluationReuse=true;cfg.diagnosticBvMethodsSharedPorts=true;
    NewtonSolver optimized(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",0.}},cfg);
    auto a=f.assembler();const auto physical=f.point();DDSolution state;
    state.psi=physical.head(6);state.phin=physical.segment(6,6);state.phip=physical.tail(6);
    state.n=a->electronDensity(physical);state.p=a->holeDensity(physical);
    auto x=baseline.packArclengthState(state);
    const auto original=baseline.makeArclengthSystem("drain");
    const auto reduced=optimized.makeArclengthSystem("drain",1e-4,true);
    auto source=optimized.makeArclengthContactCurrentFunctional("source",1.);
    auto drain=optimized.makeArclengthContactCurrentFunctional("drain",1.);
    auto oldSource=baseline.makeArclengthContactCurrentFunctional("source",1.);
    auto oldDrain=baseline.makeArclengthContactCurrentFunctional("drain",1.);
    for(double bias:{0.,1e-5,.2}) {
        REQUIRE(exactValueBits(reduced.residual(x,bias))==exactValueBits(original.residual(x,bias)));
        REQUIRE(exactValueBits(reduced.jacobian(x,bias))==exactValueBits(original.jacobian(x,bias)));
        REQUIRE(exactValueBits(reduced.parameterDerivative(x,bias))==exactValueBits(original.parameterDerivative(x,bias)));
        REQUIRE(exactValueBits(source.value(x,bias))==exactValueBits(oldSource.value(x,bias)));
        REQUIRE(exactValueBits(drain.value(x,bias))==exactValueBits(oldDrain.value(x,bias)));
        REQUIRE(exactValueBits(source.linearize(x,bias))==exactValueBits(oldSource.linearize(x,bias)));
        REQUIRE(exactValueBits(drain.linearize(x,bias))==exactValueBits(oldDrain.linearize(x,bias)));
        x[7]+=1e-13;
    }
    x.segment(6,6).setConstant(1e6);
    REQUIRE_FALSE(original.parameterDerivative(x,.1).allFinite());
    REQUIRE_FALSE(reduced.parameterDerivative(x,.1).allFinite());
}


#include "vela/solver/PopulationState.h"
TEST_CASE("Direct population BVMethods operator matches the physical QF residual and chain rule", "[population_core]") {
    using namespace vela;
    ProfileFixture f;auto q=f.assembler();const VectorXd x=f.point();
    const VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);
    f.scaling.stateVariables="potential_density";auto a=f.assembler();
    const VectorXd y=a->packPopulation(x.head(6),en,hp);
    const auto state=a->unpack(y);
    REQUIRE((state.n-en).norm()==0.);REQUIRE((state.p-hp).norm()==0.);
    REQUIRE((state.phin-x.segment(6,6)).head(4).norm()<2e-16);
    Eigen::MatrixXd T=Eigen::MatrixXd::Zero(18,18);
    for(int i=0;i<6;++i){T(i,i)=1.;if(i<4){
        T(6+i,i)=T(12+i,i)=1.;
        T(6+i,6+i)=-constants::Vt_300/experimental::JoyceDixonStatistics::densitySlopeFromPopulation(en[i],q->electronDensityOfStates()[i]);
        T(12+i,12+i)=constants::Vt_300/experimental::JoyceDixonStatistics::densitySlopeFromPopulation(hp[i],q->holeDensityOfStates()[i]);
    }else {T(6+i,6+i)=T(12+i,12+i)=1.;}}
    CoupledDDBoundaryConditions bc;bc.psi[0]=.04;bc.phin[0]=.003;bc.phip[0]=-.004;
    for(const auto& b:std::vector<CoupledDDBoundaryConditions>{{},bc}){
        const VectorXd old=q->residual(x,b),now=a->residual(y,b);
        for(int block=0;block<3;++block)
            REQUIRE((old-now).segment(block*6,6).norm()<1e-11*std::max(old.segment(block*6,6).norm(),1e-25));
        const Eigen::MatrixXd expected=Eigen::MatrixXd(q->assembleJacobian(x,b))*T;
        const Eigen::MatrixXd actual=a->assembleJacobian(y,b);
        const VectorXd beforeBootstrap=a->residual(y,b);
        const Eigen::MatrixXd bootstrap=a->populationBootstrapJacobian(y,b);
        REQUIRE((beforeBootstrap-a->residual(y,b)).norm()==0.);
        REQUIRE((bootstrap.topRows(6)-actual.topRows(6)).norm()==0.);
        REQUIRE(bootstrap.allFinite());
        VectorXd scale=VectorXd::Ones(18);scale.segment(6,6)=en;scale.tail(6)=hp;
        for(int i=4;i<6;++i)scale[6+i]=scale[12+i]=1.;
        const Eigen::MatrixXd es=expected*scale.asDiagonal(),as=actual*scale.asDiagonal();
        for(int block=0;block<3;++block)
            REQUIRE((es-as).middleRows(block*6,6).norm()<2e-11*std::max(es.middleRows(block*6,6).norm(),1e-25));
        VectorXd v=VectorXd::Zero(18);
        for(int i=0;i<6;++i)v[i]=.01*std::sin(i+.1);
        for(int i=0;i<4;++i){v[6+i]=en[i]*.03*std::sin(i+.4);v[12+i]=hp[i]*.03*std::cos(i+.8);}
        const VectorXd an=actual*v;
        const double h=1e-5;
        const VectorXd fd=(a->residual(y+h*v,b)-a->residual(y-h*v,b))/(2*h);
        for(int block=0;block<3;++block)
            REQUIRE((an-fd).segment(block*6,6).norm()<2e-6*std::max(an.segment(block*6,6).norm(),1e-25));
    }
    VectorXd broken=y;broken[6]=0.;REQUIRE_THROWS_AS(a->unpack(broken),std::domain_error);
    broken=y;broken[10]=1.;REQUIRE_THROWS_AS(a->unpack(broken),std::domain_error);
}
TEST_CASE("Core population protection repeats on each new candidate and separates global scaling", "[population_core]") {
    using namespace vela;
    VectorXd before(6),raw(6);before<<.2,.3,100.,200.,300.,400.;raw<<.01,.02,-120.,10.,20.,-500.;
    PopulationUpdateConfig c;c.localRetainedFraction=.5;
    const auto first=makePopulationTrial(before,raw,1.,{true,true},{true,true},c);
    REQUIRE(first.state[2]==50.);REQUIRE(first.state[3]==210.);
    REQUIRE(first.state[4]==320.);REQUIRE(first.state[5]==200.);
    REQUIRE(first.protectedElectrons==std::vector<int>{0});
    REQUIRE(first.protectedHoles==std::vector<int>{1});
    const auto second=makePopulationTrial(first.state,raw,1.,{true,true},{true,true},c);
    REQUIRE(second.state[2]==25.);REQUIRE(second.state[5]==100.);
    const auto shorter=makePopulationTrial(before,raw,.1,{true,true},{true,true},c);
    REQUIRE(shorter.protectedElectrons.empty());REQUIRE(shorter.protectedHoles.empty());
    c.protection="fraction_to_boundary";
    const auto global=makePopulationTrial(before,raw,1.,{true,true},{true,true},c);
    REQUIRE(global.alpha<1.);REQUIRE(global.protectedElectrons.empty());
    REQUIRE(global.state[0]-before[0]==Approx(global.alpha*raw[0]));
}

TEST_CASE("Population configuration rejects solver paths that would revert to QF coordinates", "[population_core]") {
    using namespace vela;
    nlohmann::json c={{"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}}};
    for(const auto* method:{"gummel","gummel_newton","frozen_state"}) {
        c["method"]=method;REQUIRE_THROWS_AS(newtonConfigFromJson(c),std::invalid_argument);
    }
    c["method"]="newton";
    REQUIRE(newtonConfigFromJson(c).stateVariables=="potential_density");
}

TEST_CASE("Population core callbacks preserve authoritative density and project each finite candidate", "[population_core]") {
    using namespace vela;
    ProfileFixture f;auto old=f.assembler();const auto x=f.point();
    DDSolution s;s.psi=x.head(6);s.phin=VectorXd::Constant(6,99.);s.phip=VectorXd::Constant(6,-99.);
    s.n=old->electronDensity(x);s.p=old->holeDensity(x);
    NewtonConfig cfg;cfg.coupledOperator="bvmethods_tri3_qf";cfg.stateVariables="potential_density";
    cfg.carrierStatistics=f.statistics;cfg.mobility=f.mobility;cfg.impactIonization=f.avalanche;
    cfg.populationUpdate.localRetainedFraction=.5;cfg.warmStart=true;
    NewtonSolver solver(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},cfg);
    const auto y=solver.packArclengthState(s);const auto restored=solver.unpackArclengthState(y);
    REQUIRE((restored.n-s.n).norm()==0.);REQUIRE((restored.p-s.p).norm()==0.);
    REQUIRE(restored.stateVariables=="potential_density");REQUIRE(restored.phinIncrement.size()==0);
    REQUIRE((restored.phin.head(4)-x.segment(6,4)).norm()<2e-16);
    const auto sys=solver.makeArclengthSystem("drain");
    REQUIRE(sys.prescribedRows(y,.01).size()==10);
    VectorXd raw=VectorXd::Zero(18);raw[7]=-2*y[7];raw[13]=-3*y[13];raw[1]=.001;
    const auto trial=sys.makeTrial({y,.01},raw,0.,1.);
    REQUIRE(trial.state.x[7]==.5*y[7]);REQUIRE(trial.state.x[13]==.5*y[13]);
    REQUIRE(trial.protectedElectrons==1);REQUIRE(trial.protectedHoles==1);
    const auto reduced=sys.makeTrial({y,.01},raw,0.,.1);
    REQUIRE(reduced.protectedElectrons==0);REQUIRE(reduced.protectedHoles==0);
    DensityUpdateConvergenceConfig errorConfig;
    const auto error=solver.makeArclengthUpdateErrorEvaluator("drain",errorConfig)(y,raw,trial.state.x,.01,0.);
    REQUIRE(error.valid);REQUIRE(error.counts[1]==2);REQUIRE(error.counts[2]==2);
    const auto port=solver.makeArclengthContactCurrentFunctional("drain",1.);
    const auto derivative=port.linearize(y,.01);
    // The postprocessing path must consume the same primary populations even
    // when a caller supplies stale derived QF vectors.
    auto stale=restored;stale.phin.setConstant(99.);stale.phip.setConstant(-99.);
    auto scaling=f.scaling;scaling.stateVariables="potential_density";
    ContactCurrent post(f.mesh,f.materials,f.doping,f.mobility,300.,scaling,{},f.statistics);
    const auto terminal=post.compute(stale,"drain");
    REQUIRE(terminal.totalCurrent==Approx(derivative.value).epsilon(1e-13));
    VectorXd v=VectorXd::Zero(18);v[1]=.003;v[7]=.01*y[7];v[13]=-.02*y[13];
    const double h=1e-5;
    const double fd=(port.value(y+h*v,.01)-port.value(y-h*v,.01))/(2*h);
    REQUIRE(derivative.stateDerivative.dot(v)==Approx(fd).epsilon(2e-6));
    auto reportConfig=cfg;reportConfig.globalContinuityClosure.mode="report";
    reportConfig.globalContinuityClosure.tolerance=1e-12;
    NewtonSolver reported(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},reportConfig);
    const auto reportGates=reported.makeArclengthConvergenceEvaluator("drain")(y,.01);
    REQUIRE_FALSE(reportGates.globalContinuity.satisfied);
    REQUIRE(reported.makeArclengthSystem("drain").convergenceSatisfied(y,.01));
    reportConfig.globalContinuityClosure.mode="enforce";
    NewtonSolver enforced(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},reportConfig);
    REQUIRE_FALSE(enforced.makeArclengthSystem("drain").convergenceSatisfied(y,.01));
}

TEST_CASE("Population continuation restart preserves the physical predictor and next coupled solution", "[population_core][state_archive]") {
    using namespace vela;
    ProfileFixture f;
    NewtonConfig cfg;cfg.coupledOperator="bvmethods_tri3_qf";cfg.stateVariables="potential_density";
    cfg.carrierStatistics=f.statistics;cfg.mobility=f.mobility;cfg.maxIter=100;cfg.verbose=false;
    // Normalize equations while retaining this fixture's physical SI mesh and
    // population units. Coordinate scaling is independent of input units.
    cfg.inputScaling.mode=UnitScalingMode::UnitScaling;
    NewtonSolver owner(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",0.}},cfg);
    const auto poisson=owner.solvePoissonOnly(owner.buildPoissonBlockInitialization().coldInitial);
    INFO(poisson.failureDiagnostics.failureReason);REQUIRE(poisson.converged);
    REQUIRE(poisson.solution.stateVariables=="potential_density");
    REQUIRE(poisson.solution.packedState.size()==0);
    const auto sys=owner.makeArclengthSystem("drain");
    const auto port=owner.makeArclengthContactCurrentFunctional("drain",1e-6);
    CoupledLoadLineConfig control;control.linearSolver="direct_bordered";control.lineSearchMode="full_step";
    control.updateConvergence.mode="enforce_joint";control.maxIterations=40;
    CoupledLoadLineNewton newton(sys,port.value,control,port,owner.makeArclengthUpdateErrorEvaluator("drain",control.updateConvergence));
    const ArclengthState previous{owner.packArclengthState(poisson.solution),0.};
    const auto first=newton.solve(previous,.001,1e3);
    INFO(first.failureReason);REQUIRE(first.converged);
    REQUIRE(sys.predictState(previous.x,previous.x,1.)==previous.x);
    REQUIRE(sys.predictState(previous.x,first.state.x,0.)==first.state.x);
    const auto metadata=nlohmann::json{{"mode","dd"},{"mesh_sha256",std::string(64,'a')},{"potential_origin_V",0.}};
    auto roundtrip=[&](const ArclengthState& state){
        const auto physical=owner.unpackArclengthState(state.x);
        const auto disk=restoreDDSolution(archiveDDSolution(physical,metadata,cfg.inputScaling),cfg.inputScaling);
        return ArclengthState{owner.packArclengthState(disk),state.lambda};
    };
    const auto originalPredictor=newton.secantPredict(previous,first.state,.002,1e3);
    const auto restoredPredictor=newton.secantPredict(roundtrip(previous),roundtrip(first.state),.002,1e3);
    REQUIRE(exactValueBits(originalPredictor.x)==exactValueBits(restoredPredictor.x));
    REQUIRE(originalPredictor.lambda==restoredPredictor.lambda);
    const auto uninterrupted=newton.solve(originalPredictor,.002,1e3);
    const auto restarted=newton.solve(restoredPredictor,.002,1e3);
    INFO(restarted.failureReason);REQUIRE(uninterrupted.converged);REQUIRE(restarted.converged);
    REQUIRE(exactValueBits(uninterrupted.state.x)==exactValueBits(restarted.state.x));
    REQUIRE(uninterrupted.state.lambda==restarted.state.lambda);
    REQUIRE(uninterrupted.iterations==restarted.iterations);
    REQUIRE(restarted.updateError.valid);REQUIRE(restarted.updateError.rms<1.);
}

TEST_CASE("Population log predictor preserves depletion ratios and rejects unrepresentable history", "[population_core]") {
    using namespace vela;
    VectorXd previous(6),current(6);previous<<0.,0.,100.,0.,1e-200,0.;current<<1.,2.,10.,0.,1e-100,0.;
    const auto next=predictPopulation(previous,current,1.);
    REQUIRE(next[0]==2.);REQUIRE(next[1]==4.);
    REQUIRE(next[2]==Catch::Approx(1.).epsilon(1e-13));REQUIRE(next[3]==0.);
    REQUIRE(next[4]==Catch::Approx(1.).epsilon(1e-12));REQUIRE(next[5]==0.);
    REQUIRE(previous[2]==100.);REQUIRE(current[2]==10.);
    current[2]=0.;REQUIRE_THROWS(predictPopulation(previous,current,1.));
    current[2]=1e-300;REQUIRE_THROWS(predictPopulation(previous,current,4.));
}


TEST_CASE("Population paired Bernoulli preserves the zero tangent and reversal", "[population_sg][tri3_ad]") {
    for(double value:{-1000.,-100.,-10.,-.5,-1e-8,0.,1e-8,.5,10.,100.,1000.}){
        const auto x=AvJet::variable(value,0);const auto pair=populationBernoulliPair(x);
        const auto expected=avBernoulli(x),opposite=avBernoulli(-x);
        REQUIRE(pair[0].value==Approx(expected.value).epsilon(2e-14).margin(1e-300));
        REQUIRE(pair[1].value==Approx(opposite.value).epsilon(2e-14).margin(1e-300));
        REQUIRE(pair[0].derivative[0]==Approx(expected.derivative[0]).epsilon(3e-13).margin(1e-300));
        REQUIRE(pair[1].derivative[0]==Approx(opposite.derivative[0]).epsilon(3e-13).margin(1e-300));
        const auto forward=populationSgFlux(AvJet(1),AvJet(2e19),AvJet(3e19),x);
        const auto reverse=populationSgFlux(AvJet(1),AvJet(3e19),AvJet(2e19),-x);
        REQUIRE(forward.value==-reverse.value);
        REQUIRE(forward.derivative[0]==-reverse.derivative[0]);
    }
}

TEST_CASE("Population band SG retains conservation and complete finite-state Jv", "[population_sg][population_core]") {
    using namespace vela;ProfileFixture f;auto q=f.assembler();const VectorXd x=f.point();
    const VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);
    f.scaling.stateVariables="potential_density";f.scaling.populationSgCurrent="band_population";auto a=f.assembler();
    const VectorXd y=a->packPopulation(x.head(6),en,hp);
    const VectorXd old=q->residual(x,{}),now=a->residual(y,{});
    for(int b=0;b<3;++b)REQUIRE((old-now).segment(b*6,6).norm()<2e-11*std::max(old.segment(b*6,6).norm(),1e-25));
    VectorXd v=VectorXd::Zero(18);for(int i=0;i<6;++i)v[i]=.01*std::sin(i+.1);
    for(int i=0;i<4;++i){v[6+i]=en[i]*.03*std::sin(i+.4);v[12+i]=hp[i]*.03*std::cos(i+.8);}
    const Eigen::MatrixXd J=a->assembleJacobian(y,{});const VectorXd an=J*v;
    const double h=1e-5;const VectorXd fd=(a->residual(y+h*v,{})-a->residual(y-h*v,{}))/(2*h);
    for(int b=0;b<3;++b)REQUIRE((an-fd).segment(b*6,6).norm()<2e-6*std::max(an.segment(b*6,6).norm(),1e-25));
    auto config=nlohmann::json{{"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},{"population_sg_current","band_population"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}}};
    REQUIRE(newtonConfigFromJson(config).populationSgCurrent=="band_population");
    config["state_variables"]="potential_quasi_fermi";REQUIRE_THROWS(newtonConfigFromJson(config));
}


TEST_CASE("Primitive population state preserves normalized authority and reference origin", "[population_primitive]") {
    BvMethodsPrimitiveReference r{.02585,1e10,4.07274,178.9,1.108};
    const double n=1.2345678901234567;const auto a=populationPrimitiveFromNormalized(AvJet(.5),AvJet(n),2e9,2.,r);
    REQUIRE(a.population_overN0.value==n);
    auto shifted=r;shifted.energyReference_overVT+=.5;
    const auto b=populationPrimitiveFromNormalized(AvJet(1.),AvJet(n),2e9,2.,shifted);
    REQUIRE(primitiveQuasiFermiPotential(a,r,true).value==primitiveQuasiFermiPotential(b,shifted,true).value);
    REQUIRE_THROWS(populationPrimitiveFromNormalized(AvJet(.5),AvJet(0),2e9,2.,r));
    auto invalid=r;invalid.concentration_cm3=0.;REQUIRE_THROWS(invalid.validate());
}
TEST_CASE("Primitive SG has complete finite-state population tangents and equivalent band drift", "[population_primitive]") {
    const BvMethodsPrimitiveReference r{.02585,1e10,4.07274,178.9,1.108};
    const std::array<double,4> initial{.3,.32,1e18,2e18};
    auto evaluate=[&](const std::array<double,4>& values,bool ad){
        std::array<AvJet,4> z;for(int k=0;k<4;++k)z[k]=ad?AvJet::variable(values[k],k):AvJet(values[k]);
        const auto a=populationPrimitiveFromPhysical(z[0],z[2],3e19,.02,r),b=populationPrimitiveFromPhysical(z[1],z[3],3e19,.04,r);
        const auto argument=primitiveSgArgument(a,b,true);
        const auto previous=populationSgArgument(z[2],z[3],z[0],z[1],1e10*std::exp(.02/(2*r.thermalVoltage_V)),1e10*std::exp(.04/(2*r.thermalVoltage_V)),3e19,3e19,r.thermalVoltage_V,true);
        REQUIRE(std::abs(argument.value-previous.value)<1e-13);
        for(int k=0;k<4;++k)REQUIRE(std::abs(argument.derivative[k]-previous.derivative[k])<=1e-13*std::max(1.,std::abs(previous.derivative[k])));
        return primitiveSgParticleFlux(a,b,r,1e-5,true);
    };
    const auto ad=evaluate(initial,true);const std::array<double,4> direction{.01,-.02,initial[2]*.03,-initial[3]*.04};
    double tangent=0.;for(int k=0;k<4;++k)tangent+=ad.derivative[k]*direction[k];
    auto lo=initial,hi=initial;const double h=1e-5;for(int k=0;k<4;++k){lo[k]-=h*direction[k];hi[k]+=h*direction[k];}
    const double finite=(evaluate(hi,false).value-evaluate(lo,false).value)/(2*h);
    REQUIRE(std::abs(finite-tangent)<2e-6*std::abs(tangent));
}

TEST_CASE("Primitive electron and hole QF share the material energy reference", "[population_primitive]") {
    using namespace vela;
    const double vt=.02585,nc=3e19,nv=2e19,ni0=1e10,psi=.3,bgn=.04;
    BvMethodsPrimitiveReference ref{vt,1e10,4.07274,4.07274/vt+std::log(nc/ni0),2*vt*std::log(std::sqrt(nc*nv)/ni0)};
    CarrierStatisticsConfig stats;stats.model="fermi_dirac_joyce_dixon";
    const double ni=ni0*std::exp(bgn/(2*vt));
    const auto en=populationPrimitiveFromPhysical(AvJet(psi),AvJet(1e18),nc,bgn,ref);
    const auto hp=populationPrimitiveFromPhysical(AvJet(psi),AvJet(2e17),nv,bgn,ref);
    REQUIRE(std::abs(primitiveQuasiFermiPotential(en,ref,true).value-electronQuasiFermiPotential(ni,nc,psi,1e18,vt,stats))<1e-13);
    REQUIRE(std::abs(primitiveQuasiFermiPotential(hp,ref,false).value-holeQuasiFermiPotential(ni,nv,psi,2e17,vt,stats))<1e-13);
}

TEST_CASE("Primitive current preserves canonical reversal and complete triangle tangent", "[population_primitive]") {
    BvMethodsPrimitiveReference ref{.02585,1e10,4.07274,178.9,1.108};
    AvTri3Geometry g;g.points_cm={{{0.,0.},{2e-5,0.},{0.,1e-5}}};g.nodeOrder={8,3,5};g.representative=1;
    g.currentProjection={.9,.1,-.1,1.1};
    std::array<double,6> z{.3,.32,.29,1e18,2e18,8e17};
    const std::array<double,6> d{.01,-.02,.03,3e16,-4e16,2e16};
    auto evaluate=[&](const std::array<double,6>& x,bool ad,bool electron){
        std::array<BvMethodsPrimitiveCarrier,3> nodes;
        for(int k=0;k<3;++k)nodes[k]=populationPrimitiveFromPhysical(ad?AvJet::variable(x[k],k):AvJet(x[k]),
            ad?AvJet::variable(x[3+k],3+k):AvJet(x[3+k]),3e19,.02+.01*k,ref);
        return primitiveAvCurrent(g,nodes,ref,electron);
    };
    for(bool electron:{true,false}){
        auto lo=z,hi=z;const double h=1e-5;for(int k=0;k<6;++k){lo[k]-=h*d[k];hi[k]+=h*d[k];}
        const auto ad=evaluate(z,true,electron),a=evaluate(lo,false,electron),b=evaluate(hi,false,electron);
        for(int axis=0;axis<2;++axis){double v=0.;for(int k=0;k<6;++k)v+=ad[axis].derivative[k]*d[k];
            REQUIRE(std::abs((b[axis].value-a[axis].value)/(2*h)-v)<2e-6*std::abs(v));}
        const auto original=evaluate(z,false,electron);g.nodeOrder={0,1,2};const auto reversed=evaluate(z,false,electron);
        for(int axis=0;axis<2;++axis)REQUIRE(original[axis].value==reversed[axis].value);
        g.nodeOrder={8,3,5};
    }
}

TEST_CASE("Primitive stable statistics preserves dilute gamma and its tangent", "[population_primitive]") {
    BvMethodsPrimitiveReference ref{.02585,1e10,4.05,178.,1.1};
    const auto n=AvJet::variable(1e-8,0);const double dos=2e9;
    const auto old=populationPrimitiveFromNormalized(AvJet(0),n,dos,.1,ref);
    REQUIRE(old.logGamma.value==populationLogGamma(n,dos).value);
    ref.stableStatistics=true;
    const auto stable=populationPrimitiveFromNormalized(AvJet(0),n,dos,.1,ref);
    REQUIRE(stable.logGamma.value==Approx(-.3535533905932737622*1e-8/dos).epsilon(2e-15));
    REQUIRE(stable.logGamma.derivative[0]==Approx(-.3535533905932737622/dos).epsilon(2e-15));
    REQUIRE(stable.logGamma.value<0.);
}


TEST_CASE("BVMethods residual values are independent of derivative request and refresh on changed populations", "[bvmethods_profile][population_core]") {
    using namespace vela;ProfileFixture f;auto a=f.assembler();const auto x=f.point();
    detail::BvMethodsQfState state{x.head(6),x.segment(6,6),x.tail(6),a->electronDensity(x),a->holeDensity(x)};
    const auto nc=a->electronDensityOfStates(),nv=a->holeDensityOfStates(),ni=a->intrinsicDensity();
    detail::BvMethodsTri3Context context{f.mesh,f.materials,f.doping,
        {f.mobility,f.avalanche,{},CarrierStatisticsModel::FermiDiracJoyceDixon}};
    for(bool population:{false,true})for(bool direct:{false,true}) {
        if(direct&&!population)continue;
        state.populationCoordinates=population;
        detail::BvMethodsTri3Operator element(context,state,nc,nv,1.,1.,constants::Vt_300,ni,direct);
        for(int trial=0;trial<2;++trial) {
            const auto value=element.transport(state,false,true,true),ad=element.transport(state,true,true,true);
            REQUIRE(exactValueBits(value.source)==exactValueBits(ad.source));
            REQUIRE(exactValueBits(value.absoluteSource)==exactValueBits(ad.absoluteSource));
            const auto av=element.evaluate(state,false),avAd=element.evaluate(state,true);
            REQUIRE(exactValueBits(av.source)==exactValueBits(avAd.source));
            REQUIRE(exactValueBits(av.electronGeneration)==exactValueBits(avAd.electronGeneration));
            REQUIRE(exactValueBits(av.holeGeneration)==exactValueBits(avAd.holeGeneration));
            auto changed=state;changed.n[1]*=1.01;changed.p[2]*=.99;
            const auto fresh=element.transport(changed,false,true,true);
            REQUIRE((fresh.source-value.source).norm()>0.);
            state=changed;
        }
    }
}


TEST_CASE("BVMethods fused diagnostics retain constrained and free physical equations", "[bvmethods_profile][population_core]") {
    using namespace vela;ProfileFixture f;auto qf=f.assembler();const auto x=f.point();
    for(bool population:{false,true}) {
        f.scaling.stateVariables=population?"potential_density":"potential_quasi_fermi";
        auto a=f.assembler();VectorXd y=population?a->packPopulation(x.head(6),qf->electronDensity(x),qf->holeDensity(x)):x;
        CoupledDDBoundaryConditions b;b.psi[0]=.1;b.phin[0]=.002;b.phip[0]=-.001;
        for(int trial=0;trial<3;++trial) {
            if(trial==1)y[7]=std::nextafter(y[7],1.);
            if(trial==2)b.psi[0]+=.01;
            const auto fused=a->bvMethodsPhysicalDiagnostics(y,b);
            REQUIRE(exactValueBits(fused.residual)==exactValueBits(a->residual(y,b)));
            REQUIRE(exactValueBits(fused.constrainedTerms)==exactValueBits(a->carrierContinuityEquationTermDiagnostics(y,b)));
            REQUIRE(exactValueBits(fused.unconstrainedTerms)==exactValueBits(a->carrierContinuityEquationTermDiagnostics(y,{})));
            REQUIRE_FALSE(fused.constrainedTerms[0].electronContinuityActive);
            REQUIRE(fused.unconstrainedTerms[0].electronContinuityActive);
        }
    }
}


TEST_CASE("Population numerical policy rejects inactive QF controls before assembly", "[population_core]") {
    using namespace vela;
    const nlohmann::json base={{"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}}};
    for(const auto& option:std::vector<nlohmann::json>{
        {{"quasi_fermi_update_limit_V",.1}},{{"quasi_fermi_update_limit_minority_V",.1}},{{"max_update",1.}},
        {{"quasi_fermi_update_limit_mode","uniform_trust_region"}},
        {{"quasi_fermi_trust_region_growth_factor",2.}},{{"quasi_fermi_recenter_on_stall",true}},
        {{"contact_majority_qf_branch_drop_limit_V",1e-9}},{{"stable_merit_comparison",true}},
        {{"carrier_row_qualified_stall_acceptance",true}},{{"linear_equilibration",{{"mode","l2_row_column"}}}},
        {{"jacobian","finite_difference"}}}) {
        auto c=base;c.update(option);CAPTURE(option.dump());
        REQUIRE_THROWS_AS(newtonConfigFromJson(c),std::invalid_argument);
    }
    ProfileFixture f;
    auto cfg=newtonConfigFromJson(base);cfg.mobility=f.mobility;cfg.impactIonization=f.avalanche;
    for(const auto& change:std::vector<std::function<void(NewtonConfig&)>>{
        [](auto& c){c.quasiFermiUpdateLimit_V=.1;},[](auto& c){c.maxUpdate=1.;},
        [](auto& c){c.quasiFermiRecenterOnInitialState=true;},
        [](auto& c){c.quasiFermiTrustRegionGrowthFactor=2.;},
        [](auto& c){c.contactMajorityQfBranchDropLimit_V=1e-9;},
        [](auto& c){c.stableMeritComparison=true;}}) {
        auto c=cfg;change(c);
        REQUIRE_THROWS_AS((NewtonSolver(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},c)),std::invalid_argument);
    }
    NewtonSolver valid(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},cfg);
    REQUIRE_FALSE(valid.makeArclengthSystem("drain").limitUpdate);
}

TEST_CASE("Population absolute block stopping replaces scalar tolerances without bypassing blocks", "[population_core]") {
    using namespace vela;ProfileFixture f;
    NewtonConfig cfg;cfg.coupledOperator="bvmethods_tri3_qf";cfg.stateVariables="potential_density";
    cfg.populationStoppingRule="legacy";cfg.primaryVariableCoordinates="physical";
    cfg.populationDerivativeRecovery="full";
    cfg.carrierStatistics=f.statistics;cfg.mobility=f.mobility;cfg.impactIonization=f.avalanche;
    cfg.maxIter=1;cfg.reltol=1e-40;cfg.abstol=1e-80;cfg.lineSearch=false;
    NewtonSolver owner(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},cfg);
    const auto qf=f.assembler();const auto x=f.point();DDSolution state;
    state.psi=x.head(6);state.n=qf->electronDensity(x);state.p=qf->holeDensity(x);
    state.phin=x.segment(6,6);state.phip=x.tail(6);
    const auto y=owner.packArclengthState(state);
    const auto rhs=owner.makeArclengthSystem("source").residual(y,0.);
    REQUIRE(rhs.allFinite());REQUIRE(rhs.norm()>0.);
    cfg.blockAbsoluteConvergence={"enforce",2*std::max(rhs.head(6).cwiseAbs().maxCoeff(),1e-100),
        2*std::max(rhs.segment(6,6).cwiseAbs().maxCoeff(),1e-100),
        2*std::max(rhs.tail(6).cwiseAbs().maxCoeff(),1e-100)};
    NewtonSolver blocks(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},cfg);
    const auto result=blocks.solve(state);
    INFO(result.failureDiagnostics.failureReason);REQUIRE(result.converged);REQUIRE(result.iters==0);
    REQUIRE(result.finalResidualNorm<=1.);
    cfg.blockAbsoluteConvergence.psiResidualCeiling=rhs.head(6).cwiseAbs().maxCoeff()/2;
    REQUIRE(cfg.blockAbsoluteConvergence.psiResidualCeiling>0.);
    NewtonSolver strict(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},cfg);
    REQUIRE_FALSE(strict.makeArclengthSystem("source").convergenceSatisfied(y,0.));
    const auto rejected=strict.solve(state);
    REQUIRE_FALSE((rejected.converged&&rejected.iters==0));
}

TEST_CASE("Log density ratio retains finite tangents across depleted populations", "[bvmethods_profile][population_core]") {
    using namespace vela::experimental;
    auto a=AvJet::variable(1e-200,0),b=AvJet::variable(1e-300,1);
    const auto g=avLogRatio(a,b);
    REQUIRE(std::isfinite(g.derivative[0]));REQUIRE(std::isfinite(g.derivative[1]));
    REQUIRE(g.derivative[0]==Catch::Approx(1e200));
    REQUIRE(g.derivative[1]==Catch::Approx(-1e300));
    const auto reverse=avLogRatio(b,a);
    REQUIRE(reverse.value==Catch::Approx(-g.value));
    REQUIRE(reverse.derivative[0]==Catch::Approx(-g.derivative[0]));
    REQUIRE(reverse.derivative[1]==Catch::Approx(-g.derivative[1]));
}

TEST_CASE("Population transport evaluates depleted tangents without overflowing intermediates", "[bvmethods_profile][population_core]") {
    using namespace vela;ProfileFixture f;auto q=f.assembler();const auto x=f.point();
    VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);
    hp[1]=1e-296;hp[2]=1e-194;en[1]=1e-232;
    f.scaling.stateVariables="potential_density";auto a=f.assembler();
    const VectorXd y=a->packPopulation(x.head(6),en,hp);
    REQUIRE(a->bvMethodsTransportResidual(y).allFinite());
    const auto j=a->bvMethodsTransportJacobian(y);
    for(int k=0;k<j.outerSize();++k)for(SparseMatrixd::InnerIterator v(j,k);v;++v)REQUIRE(std::isfinite(v.value()));
}

TEST_CASE("Depleted population sources keep complete bootstrap tangents finite", "[bvmethods_profile][population_core]") {
    using namespace vela;ProfileFixture f;auto q=f.assembler();const auto x=f.point();
    VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);hp[1]=1e-308;en[2]=1e-239;
    f.scaling.stateVariables="potential_density";auto a=f.assembler();
    const auto y=a->packPopulation(x.head(6),en,hp);
    REQUIRE(a->residual(y,{}).allFinite());
    const auto j=a->populationBootstrapJacobian(y,{});
    for(int k=0;k<j.outerSize();++k)for(SparseMatrixd::InnerIterator v(j,k);v;++v)REQUIRE(std::isfinite(v.value()));
}

TEST_CASE("Population minimum is distinct from positivity and excludes constrained slots", "[population_core]") {
    using namespace vela;VectorXd x(6);x<<.1,.2,1e-110,0.,1e-110,5e-120;
    VectorXd d=VectorXd::Zero(6);d[2]=-2e-110;d[4]=-.5e-110;
    PopulationUpdateConfig c;c.minimumDensity=1e-100;
    auto t=makePopulationTrial(x,d,1.,{true,false},{true,false},c);
    REQUIRE(t.state[2]==1e-100);REQUIRE(t.state[4]==1e-100);
    REQUIRE(t.state[3]==0.);REQUIRE(t.state[5]==5e-120);
    REQUIRE(t.protectedElectrons.size()==1);REQUIRE(t.protectedHoles.empty());
    REQUIRE(t.flooredElectrons.size()==1);REQUIRE(t.flooredHoles.size()==1);
    c.minimumDensity=0.;auto old=makePopulationTrial(x,d,1.,{true,false},{true,false},c);
    REQUIRE(old.state[2]==Catch::Approx(.01*x[2]));REQUIRE(old.state[4]==.5*x[4]);
}
TEST_CASE("Predictor low-density threshold retains density without clipping it", "[population_core]") {
    using namespace vela;VectorXd a(6),b(6);a<<0.,0.,1e-110,0.,1e-90,0.;b<<.1,.2,1e-120,0.,1e-80,0.;
    auto y=predictPopulation(a,b,1.,1e-100);
    REQUIRE(y[0]==.2);REQUIRE(y[2]==1e-120);REQUIRE(y[2]<1e-100);
    REQUIRE(y[4]==Catch::Approx(1e-70));REQUIRE(y[3]==0.);
    const auto baseline=predictPopulation(a,b,1.);REQUIRE(baseline[2]<b[2]);
    REQUIRE_THROWS(predictPopulation(a,b,1.,-1.));
}
TEST_CASE("Population density controls are independent JSON inputs with physical units", "[population_core]") {
    using namespace vela;
    nlohmann::json j={{"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}},
        {"population_update",{{"minimum_density_cm3",1e-100}}},
        {"population_predictor_low_density_limit_cm3",1e-100}};
    ProfileFixture f;auto cfg=newtonConfigFromJson(j);cfg.mobility=f.mobility;cfg.impactIonization=f.avalanche;
    REQUIRE(cfg.populationMinimumDensity_cm3==1e-100);REQUIRE(cfg.populationPredictorLowDensity_cm3==1e-100);
    auto a=f.assembler();const auto x=f.point();DDSolution state;
    state.psi=x.head(6);state.n=a->electronDensity(x);state.p=a->holeDensity(x);
    state.n[1]=state.p[1]=2e-110;
    NewtonSolver solver(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",.01}},cfg);
    const auto y=solver.packArclengthState(state);const auto sys=solver.makeArclengthSystem("drain");
    const auto t=sys.makeTrial({y,.01},VectorXd::Zero(18),0.,1.);
    REQUIRE(t.state.x[7]==Catch::Approx(1e-94).epsilon(1e-14));
    REQUIRE(t.state.x[13]==Catch::Approx(1e-94).epsilon(1e-14));
    VectorXd before=y;before[7]*=10.;before[13]*=10.;
    const auto predicted=sys.predictState(before,y,1.);
    REQUIRE(predicted[7]==y[7]);REQUIRE(predicted[13]==y[13]);
    auto bad=j;bad["population_update"]["minimum_density_cm3"]=-1.;REQUIRE_THROWS(newtonConfigFromJson(bad));
    bad=j;bad["population_predictor_low_density_limit_cm3"]=true;REQUIRE_THROWS(newtonConfigFromJson(bad));
}

TEST_CASE("Canonical SI density-limit aliases preserve physical concentration", "[population_core]") {
    using namespace vela;nlohmann::json j={{"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}},
        {"population_update",{{"minimum_density_m3",1e-94}}},{"population_predictor_low_density_limit_m3",1e-94}};
    const auto c=newtonConfigFromJson(j);REQUIRE(c.populationMinimumDensity_cm3==Catch::Approx(1e-100));
    REQUIRE(c.populationPredictorLowDensity_cm3==Catch::Approx(1e-100));
    j["population_update"]["minimum_density_cm3"]=1e-100;REQUIRE_THROWS(newtonConfigFromJson(j));
}

TEST_CASE("Format v2 density controls retain cm3 values through canonicalization", "[population_core]") {
    using namespace vela;
    nlohmann::json deck={{"format_version",2},{"solver",{
        {"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}},
        {"population_update",{{"minimum_density_cm3",1e-100}}},
        {"population_predictor_low_density_limit_cm3",1e-100}}}};
    const auto parsed=canonicalizeDeck(deck);
    const auto c=newtonConfigFromJson(parsed.at("solver"),parseUnitScalingConfig(deck));
    REQUIRE(c.populationMinimumDensity_cm3==Catch::Approx(1e-100));
    REQUIRE(c.populationPredictorLowDensity_cm3==Catch::Approx(1e-100));
}
TEST_CASE("Direct population columns agree with coordinate multiplication in the finite domain", "[population_core][relative_columns]") {
    using namespace vela;ProfileFixture f;auto q=f.assembler();const auto x=f.point();
    const VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);f.scaling.stateVariables="potential_density";auto a=f.assembler();
    const auto y=a->packPopulation(x.head(6),en,hp);VectorXd t=VectorXd::Ones(18);
    for(int i=6;i<18;++i)if(y[i]>0.)t[i]=y[i];
    CoupledDDBoundaryConditions b;b.psi[0]=x[0];b.phin[0]=x[6];b.phip[0]=x[12];
    const auto f0=a->residual(y,b);
    for(bool approximate:{false,true}){
        SparseMatrixd reference=approximate?a->populationBootstrapJacobian(y,b):a->assembleJacobian(y,b);
        for(int k=0;k<reference.outerSize();++k)for(SparseMatrixd::InnerIterator e(reference,k);e;++e)e.valueRef()*=t[e.col()];
        const auto scaled=a->populationJacobianInColumns(y,b,t,approximate);
        REQUIRE((scaled-reference).norm()/reference.norm()<1e-11);
    }
    REQUIRE((a->residual(y,b)-f0).isZero(0.));
    t[7]=0.;REQUIRE_THROWS(a->populationJacobianInColumns(y,b,t));
}
TEST_CASE("Extreme populations retain finite native column assembly", "[population_core][relative_columns]") {
    using namespace vela;ProfileFixture f;auto q=f.assembler();const auto x=f.point();VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);
    en[1]=1e-300;hp[1]=1e-316;
    f.scaling.stateVariables="potential_density";auto a=f.assembler();const auto y=a->packPopulation(x.head(6),en,hp);VectorXd t=VectorXd::Ones(18);t[7]=en[1];t[13]=hp[1];
    const auto rhs=a->residual(y,{});REQUIRE(rhs.allFinite());
    for(bool approximate:{false,true}){const auto j=a->populationJacobianInColumns(y,{},t,approximate);
        for(int k=0;k<j.outerSize();++k)for(SparseMatrixd::InnerIterator v(j,k);v;++v)REQUIRE(std::isfinite(v.value()));
    }
    REQUIRE((a->residual(y,{})-rhs).isZero(0.));
}

TEST_CASE("BVMethods dilute SRH population tangent retains its positive physical limit",
          "[bvmethods_profile][population_core][srh_stable]") {
    using namespace vela;
    ProfileFixture fixture;
    fixture.scaling.stateVariables="potential_density";
    fixture.avalanche.model="none";
    auto assembler=fixture.assembler();
    const RecombinationModelConfig srh;
    // Uniform dilute populations give zero flux and AV. Summing all carrier
    // rows removes the conservative transport tangent, leaving the local SRH
    // response. Its depletion limit is n*dR/dn = -R*n*taup/(ni*(taun+taup))
    // (and tau_n for holes); JD corrections here are O(ni/Nc), below 1e-8.
    for(double population:{1e-40,1e-84,1e-250}) {
        VectorXd n=VectorXd::Zero(6),p=n;n.head(4).setConstant(population);p.head(4).setConstant(population);
        const VectorXd x=assembler->packPopulation(VectorXd::Zero(6),n,p);
        const auto terms=assembler->carrierContinuityEquationTermDiagnostics(x,{});
        VectorXd columns=VectorXd::Ones(18);columns.segment(6,4)=n.head(4);columns.segment(12,4)=p.head(4);
        for(bool approximate:{false,true}) {
            const auto matrix=assembler->populationJacobianInColumns(x,{},columns,approximate);
            for(int species:{1,2}) {
                const int node=1,column=6*species+node;
                const double lifetime=species==1?srh.taup:srh.taun;
                const double expected=-terms[node].electronRecombination*
                    (population/assembler->intrinsicDensity()[node])*(lifetime/(srh.taun+srh.taup));
                REQUIRE(expected>0.);
                for(int equation:{1,2}) {
                    double actual=0.;for(int i=0;i<4;++i)actual+=matrix.coeff(equation*6+i,column);
                    REQUIRE(std::isfinite(actual));
                    REQUIRE(actual>0.);
                    REQUIRE(actual/expected==Approx(1.).epsilon(1e-7));
                }
            }
        }
    }
}

TEST_CASE("Population normalized packing preserves nonzero last-bit coordinates and converts changed scales", "[population_core][normalization]") {
    using namespace vela;ProfileFixture f;f.scaling.stateVariables="potential_density";
    f.scaling.enabled=true;f.scaling.V0=.025851999786435535;
    auto a=f.assembler();VectorXd x=VectorXd::Zero(18);
    x.head(6)<<21.200072448909626,-17.200072448909626,.1,.7,2.,3.;
    x.segment(6,4).setConstant(1e20);x.segment(12,4).setConstant(1e10);
    REQUIRE((x[0]*f.scaling.V0)/f.scaling.V0!=x[0]);
    REQUIRE(exactValueBits(a->pack(a->unpack(x)))==exactValueBits(x));
    NewtonConfig c;c.coupledOperator="bvmethods_tri3_qf";c.stateVariables="potential_density";
    c.carrierStatistics=f.statistics;c.mobility=f.mobility;c.impactIonization=f.avalanche;
    c.inputScaling.mode=UnitScalingMode::UnitScaling;
    NewtonSolver solver(f.mesh,f.materials,f.doping,{{"source",0.},{"drain",0.}},c);
    const auto physical=solver.unpackArclengthState(x);
    REQUIRE(physical.hasConsistentPopulationPotential());
    REQUIRE(exactValueBits(solver.packArclengthState(physical))==exactValueBits(x));
    auto changed=physical;changed.psi[1]+=.01;
    REQUIRE_FALSE(changed.hasConsistentPopulationPotential());
    auto z=solver.packArclengthState(changed);
    REQUIRE(z[1]==changed.psi[1]/physical.populationPotentialScale_V);
    // A genuine potential-scale change is an explicit conversion, not a stale reuse.
    changed=physical;changed.populationPotentialScale_V*=2.;changed.populationPsi/=2.;
    REQUIRE(changed.hasConsistentPopulationPotential());
    REQUIRE(solver.packArclengthState(changed)[0]==physical.psi[0]/physical.populationPotentialScale_V);
}

TEST_CASE("Local bootstrap SG feedback keeps values and potential tangents with conservative population coefficients",
          "[population_sg_bootstrap][tri3_ad]") {
    for(double drive:{.2,80.})for(bool electron:{false,true}) {
        const double sign=electron?1.:-1.,vt=.02585;
        std::array<AvJet,3> c{AvJet(1e-100),AvJet(2e-99),AvJet(5e-98)};
        std::array<AvJet,3> f{AvJet(0),AvJet(sign*drive*vt),AvJet(0)};
        std::array<AvJet,3> mu{AvJet(1),AvJet(1),AvJet(1)},g{AvJet(1),AvJet(1),AvJet(1)};
        for(int k=0;k<3;++k){c[k].derivative[3+k]=c[k].value;f[k].derivative[k]=1.;
            f[k].derivative[3+k]=sign*vt;g[k].derivative[3+k]=.07;}
        const std::array<double,3> areas{1,1,1},edges{1,0,0};
        const auto original=evaluateBvMethodsTri3Transport(areas,edges,c,f,mu,vt,electron,g);
        const auto frozen=evaluateBvMethodsTri3Transport(areas,edges,c,f,mu,vt,electron,g,nullptr,{0,1,2},{true,false,false});
        auto argument=avLogRatio(c[1],c[0])+AvJet(sign)*(f[1]-f[0])/(AvJet(vt)*g[0]);
        for(int z=3;z<9;++z)argument.derivative[z]=0.;
        auto coefficient=AvJet(vt)*g[0];for(int z=3;z<9;++z)coefficient.derivative[z]=0.;
        const auto expected=populationSgFlux(coefficient,c[0],c[1],argument);
        REQUIRE(original.edgeFlux_per_cm_s[0].value==frozen.edgeFlux_per_cm_s[0].value);
        for(int z=0;z<9;++z){
            INFO("drive="<<drive<<" electron="<<electron<<" slot="<<z);
            REQUIRE(std::isfinite(frozen.edgeFlux_per_cm_s[0].derivative[z]));
            if(z<3)REQUIRE(original.edgeFlux_per_cm_s[0].derivative[z]==frozen.edgeFlux_per_cm_s[0].derivative[z]);
            else REQUIRE(frozen.edgeFlux_per_cm_s[0].derivative[z]==Approx(expected.derivative[z]).epsilon(2e-12).margin(1e-300));
            REQUIRE(frozen.nodeFlux_per_cm_s[0].derivative[z]==-frozen.nodeFlux_per_cm_s[1].derivative[z]);
        }
        REQUIRE(frozen.edgeFlux_per_cm_s[0].derivative[3]>0.);
        REQUIRE(frozen.edgeFlux_per_cm_s[0].derivative[4]<0.);
        // An unselected edge remains byte-for-byte on its original AD path.
        for(int z=0;z<9;++z)REQUIRE(original.edgeFlux_per_cm_s[1].derivative[z]==frozen.edgeFlux_per_cm_s[1].derivative[z]);
    }
}
TEST_CASE("Statistical bootstrap is opt-in and restores full and above-floor matrices",
          "[population_sg_bootstrap][population_core]") {
    using namespace vela;ProfileFixture f;auto q=f.assembler();const VectorXd x=f.point();
    const VectorXd en=q->electronDensity(x),hp=q->holeDensity(x);
    f.scaling.stateVariables="potential_density";auto old=f.assembler();
    const VectorXd y=old->packPopulation(x.head(6),en,hp);
    const VectorXd rhs=old->residual(y,{});const Eigen::MatrixXd full=old->assembleJacobian(y,{});
    const Eigen::MatrixXd bootstrap=old->populationBootstrapJacobian(y,{});
    f.scaling.populationSgBootstrapFloor_cm3=1e30;auto controlled=f.assembler();
    REQUIRE((controlled->residual(y,{}).array()==rhs.array()).all());
    REQUIRE((Eigen::MatrixXd(controlled->assembleJacobian(y,{})).array()==full.array()).all());
    REQUIRE((Eigen::MatrixXd(controlled->populationBootstrapJacobian(y,{}))-bootstrap).norm()>0.);
    // The event policy omits mobility/AV tangents, never the separate
    // statistical floor-edge feedback. Its RHS is the same physical model.
    REQUIRE((Eigen::MatrixXd(controlled->populationProtectedJacobian(y,{})).array()==bootstrap.array()).all());
    REQUIRE((controlled->residual(y,{}).array()==rhs.array()).all());
    f.scaling.populationSgBootstrapFloor_cm3=1e-100;auto above=f.assembler();
    REQUIRE((Eigen::MatrixXd(above->populationBootstrapJacobian(y,{})).array()==bootstrap.array()).all());
    auto cfg=nlohmann::json{{"state_variables","potential_density"},{"coupled_operator","bvmethods_tri3_qf"},
        {"population_sg_bootstrap_linearization","floor_adjacent_edges"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}},
        {"population_update",{{"minimum_density_cm3",1e-100}}}};
    REQUIRE(newtonConfigFromJson(cfg).populationSgBootstrapLinearization=="floor_adjacent_edges");
    cfg["population_update"]["minimum_density_cm3"]=0.;REQUIRE_THROWS(newtonConfigFromJson(cfg));
    cfg["population_update"]["minimum_density_cm3"]=1e-100;
    cfg["population_sg_current"]="band_population";REQUIRE_THROWS(newtonConfigFromJson(cfg));
    cfg["population_sg_current"]="quasi_fermi";cfg["population_sg_bootstrap_linearization"]="invalid";REQUIRE_THROWS(newtonConfigFromJson(cfg));
}


TEST_CASE("AV exponential has finite directions at extremely small positive field",
          "[av_contract][tri3_ad][range]")
{
    const AvCoefficient coefficient{7.03e5,1.231e6,7.03e5,1.231e6,4e5};
    for(double field:{0.,1e-300,1e-200,1e-160,1e-140,1e-20,1.,1e4,4e5,1e6}) {
        const auto x=AvJet::variable(field,0);
        const auto value=avCoefficient(x,coefficient);
        REQUIRE(std::isfinite(value.value));
        for(double d:value.derivative)REQUIRE(std::isfinite(d));
        if(field<1.) {REQUIRE(value.value==0.);REQUIRE(value.derivative[0]==0.);}
        else {
            const long double exact=7.03e5L*std::exp(-1.231e6L/field);
            const long double slope=exact*1.231e6L/field/field;
            REQUIRE(value.value==Approx(double(exact)).epsilon(3e-14));
            REQUIRE(value.derivative[0]==Approx(double(slope)).epsilon(3e-14));
        }
    }
}

TEST_CASE("AD exponential scales finite directions before narrowing",
          "[av_contract][tri3_ad][range]")
{
    for(double x:{-800.,800.}) {
        AvJet input(x);input.derivative[0]=x>0.?1e-100:1e100;
        const double expected=double(std::exp((long double)x)*input.derivative[0]);
        for(const auto result:{vela::detail::dualExp(input),vela::detail::dualExpm1(input)}) {
            REQUIRE(std::isfinite(result.derivative[0]));
            REQUIRE(result.derivative[0]==Approx(expected).epsilon(3e-14));
            for(int k=1;k<9;++k)REQUIRE(result.derivative[k]==0.);
        }
    }
}

TEST_CASE("Population recombination preserves equilibrium cancellation and has no psi direction",
          "[av_contract][population_core][srh]")
{
    using HP=boost::multiprecision::cpp_dec_float_100;
    const double n0=1e10,p0=std::nextafter(n0,INFINITY);
    AvJet n=AvJet::variable(n0,3),p=AvJet::variable(p0,6);
    const auto excess=populationRecombinationExcess(n,p,AvJet(1),AvJet(1),n0);
    const double expected=double(HP(n0)*HP(p0)-HP(n0)*HP(n0));
    REQUIRE(excess.value!=0.);
    REQUIRE(excess.value==Approx(expected).epsilon(1e-15));
    REQUIRE(excess.derivative[0]==0.);
    REQUIRE(excess.derivative[3]==p0);
    REQUIRE(excess.derivative[6]==n0);
    REQUIRE(populationRecombinationExcess(n,n,AvJet(1),AvJet(1),n0).value==0.);
    // Both captured positive populations are finite even when the former
    // QF-splitting expm1 path overflowed (splitting > 709).
    for(const auto values:std::array<std::array<double,5>,2>{{
        {1.7738208624835845e23,1.1435920332085196e23,3.052789335414555e-174,2.1978091765134392e-122,8.337122390342656e10},
        {3.9192596591152703e23,1.4878424437485111e10,2.2163657908827803e-297,.9999999998305675,4.095286472506107e10}}}) {
        const auto e=populationRecombinationExcess(AvJet::variable(values[0],3),
            AvJet::variable(values[1],6),AvJet(values[2]),AvJet(values[3]),values[4]);
        const HP reference=HP(values[0])*HP(values[1])-HP(values[2])*HP(values[3])*HP(values[4])*HP(values[4]);
        REQUIRE(std::isfinite(e.value));
        REQUIRE(e.value==Approx(double(reference)).epsilon(2e-15));
        REQUIRE(e.derivative[0]==0.);
        REQUIRE(e.derivative[3]==values[1]);
        REQUIRE(e.derivative[6]==values[0]);
    }
}
