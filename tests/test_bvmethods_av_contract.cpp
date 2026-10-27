#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/physics/experimental/BvMethodsTri3Avalanche.h"
#include "vela/physics/experimental/BvMethodsTri3Transport.h"
#include <cmath>
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
