#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/solver/experimental/UpdateErrorReport.h"
#include "vela/solver/experimental/NewtonCapture.h"
#include "vela/solver/experimental/NewtonAudit.h"
#include "vela/solver/experimental/ExactValueBits.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/physics/CarrierStatistics.h"
#include <algorithm>

using namespace vela::experimental;
using Catch::Approx;

TEST_CASE("Exact result encoding preserves all exceptional floating bit patterns", "[update_error_report]") {
    using vela::experimental::exactValueBits;
    REQUIRE(exactValueBits(0.)!=exactValueBits(-0.));
    REQUIRE(exactValueBits(std::numeric_limits<double>::denorm_min())!=exactValueBits(0.));
    REQUIRE(exactValueBits(std::numeric_limits<double>::infinity())!=exactValueBits(-std::numeric_limits<double>::infinity()));
    const auto a=std::bit_cast<double>(std::uint64_t{0x7ff8000000000001});
    const auto b=std::bit_cast<double>(std::uint64_t{0x7ff8000000000002});
    REQUIRE(exactValueBits(a)!=exactValueBits(b));
    vela::NewtonResult result;
    const auto initial=exactValueBits(result);
    REQUIRE(initial==exactValueBits(result));
    auto changed=result;changed.solution.electronQuantumPotentialLike=vela::VectorXd::Constant(1,a);
    REQUIRE(initial!=exactValueBits(changed));
    changed=result;changed.failureDiagnostics.bestRejectedContactMajorityQfDrop=-0.;
    REQUIRE(initial!=exactValueBits(changed));
    changed=result;changed.trace.resize(1);const auto before=exactValueBits(changed);
    changed.trace[0].topHoleResidual.signedResidual=b;
    REQUIRE(before!=exactValueBits(changed));
    changed=result;changed.electronQuantumOuterHistory.resize(1);
    REQUIRE(initial!=exactValueBits(changed));
    // PFR visits declared fields automatically; this prevents a hand-written
    // subset from silently omitting new result members.
    boost::pfr::for_each_field(result,[&](const auto&,auto i){
        REQUIRE(initial.find(boost::pfr::get_name<decltype(i)::value,vela::NewtonResult>())!=std::string::npos);
    });
}

TEST_CASE("Independent witness owns its data and bounds diagnostic failure", "[update_error_report]") {
    NewtonAudit a,b;
    vela::VectorXd x=vela::VectorXd::Ones(3),r=x,w=x;
    vela::SparseMatrixd j(3,3);j.setIdentity();
    NewtonAuditScope::system(1,x,r,j,w);
    REQUIRE(a.iterations.empty());
    {
        NewtonAuditScope scope(&a);NewtonAuditScope::begin(true);
        NewtonAuditScope::system(1,x,r,j,w);
        x[0]=7.;j.coeffRef(0,0)=9.;
        REQUIRE(a.iterations[0].x[0]==1.);
        REQUIRE(a.iterations[0].jacobian.coeff(0,0)==1.);
        {NewtonAuditScope nested(&b);NewtonAuditScope::begin(true);
         NewtonAuditScope::record([](auto&){throw 42;});
         NewtonAuditScope::system(1,x,r,j,w);}
        REQUIRE(b.failed);REQUIRE(b.iterations.empty());
        NewtonAuditScope::direction(1,x);
        REQUIRE(a.iterations[0].linearDirection[0]==7.);
        a.limit=1;NewtonAuditScope::system(2,x,r,j,w);
        REQUIRE(a.truncated);REQUIRE(a.iterations.size()==1);
        NewtonAuditScope::direction(2,r);
        REQUIRE(a.iterations[0].linearDirection[0]==7.);
    }
}

TEST_CASE("Update error report is opt-in and preserves equation weighting", "[update_error_report]") {
    NewtonUpdateErrorReportConfig c;
    std::vector<UpdateErrorRow> rows{{UpdateErrorEquation::Poisson,1,true,0.,1.},
        {UpdateErrorEquation::Electron,2,true,0.,3.}, {UpdateErrorEquation::Electron,3,true,0.,3.},
        {UpdateErrorEquation::Contact,4,false,0.,1e100}};
    REQUIRE_FALSE(evaluateUpdateErrorReport(c,rows).enabled);
    c.mode=UpdateErrorReportMode::Report;c.digits=0;c.errRef.fill(1.);
    auto r=evaluateUpdateErrorReport(c,rows);
    REQUIRE(r.equations[0].count==1);REQUIRE(r.equations[1].count==2);
    REQUIRE(r.equations[3].count==0);
    REQUIRE(*r.rms==Approx(std::sqrt(19./3.)));
    REQUIRE_FALSE(r.completeEquationSet);
    rows.back().active=true;rows.back().delta=2.;
    auto full=evaluateUpdateErrorReport(c,rows,true);
    REQUIRE(full.equations[3].count==1);REQUIRE(*full.rms==Approx(std::sqrt(23./4.)));
    REQUIRE(full.completeEquationSet);
}

TEST_CASE("Only fixed-denominator tangent reports are homogeneous", "[update_error_report]") {
    NewtonUpdateErrorReportConfig c;c.mode=UpdateErrorReportMode::Report;
    std::vector<UpdateErrorRow> rows{{UpdateErrorEquation::Hole,918,true,1e-20,1e-5}};
    const auto before=rows;const double e=*evaluateUpdateErrorReport(c,rows).rms;
    rows[0].delta*=-7.;REQUIRE(*evaluateUpdateErrorReport(c,rows).rms==Approx(7.*e));
    REQUIRE(rows[0].denominatorState==before[0].denominatorState);
    const double p=1e-20,deta=-1e3;
    const double finite=p*std::expm1(deta),tangent=p*deta;
    REQUIRE(std::abs(finite)<=p);REQUIRE(std::abs(tangent)>p);
    rows[0].delta=finite;
    REQUIRE(*evaluateUpdateErrorReport(c,rows).rms<=p/(1e-5*(p+1e10))*(1.+1e-14));
}

TEST_CASE("Invalid or unaccepted update cannot become zero error", "[update_error_report]") {
    NewtonUpdateErrorReportConfig c;c.mode=UpdateErrorReportMode::Report;
    std::vector<UpdateErrorRow> rows{{UpdateErrorEquation::Hole,918,true,1e-20,0.,UpdateErrorValidity::NoAcceptedUpdate},
        {UpdateErrorEquation::Poisson,1,true,0.,0.}};
    auto r=evaluateUpdateErrorReport(c,rows);
    REQUIRE_FALSE(r.rms);REQUIRE_FALSE(r.equations[2].rms);REQUIRE(r.equations[2].invalidCount==1);
    REQUIRE(*r.equations[0].rms==0.);
    rows[0].validity=UpdateErrorValidity::Valid;rows[0].delta=std::numeric_limits<double>::infinity();
    REQUIRE_FALSE(evaluateUpdateErrorReport(c,rows).rms);
    rows[0].active=false;REQUIRE(*evaluateUpdateErrorReport(c,rows).rms==0.);
    c.errRef[0]=0.;REQUIRE(evaluateUpdateErrorReport(c,rows).invalidConfig);
}

TEST_CASE("Update RMS avoids overflow in finite operands and squared sums", "[update_error_report]") {
    NewtonUpdateErrorReportConfig c;c.mode=UpdateErrorReportMode::Report;c.digits=0;c.errRef.fill(1e308);
    std::vector<UpdateErrorRow> rows{{UpdateErrorEquation::Poisson,1,true,1e308,1e308}};
    REQUIRE(*evaluateUpdateErrorReport(c,rows).rms==Approx(.5));
    c.errRef.fill(1.);rows[0].denominatorState=0.;
    rows.push_back(rows[0]);REQUIRE(*evaluateUpdateErrorReport(c,rows).rms==Approx(1e308));
    rows.resize(1);rows[0].delta=1e-320;c.errRef.fill(1e7);c.digits=5;
    const auto underflow=evaluateUpdateErrorReport(c,rows);
    REQUIRE_FALSE(underflow.rms);
    REQUIRE(underflow.rows[0].row.validity==UpdateErrorValidity::UnderflowDegeneracy);
}

TEST_CASE("Fermi population tangent includes both potential and QF directions", "[update_error_report]") {
    const double ni=1e16,Nc=2.8e25,Nv=1.04e25,vt=.025852,psi=.2,fn=.05,fp=.1;
    const double dpsi=.013,dfn=-.007,dfp=.019;
    vela::CarrierStatisticsConfig statistics{"fermi_dirac"};
    const double dn=vela::electronDensityDerivativeEta(ni,Nc,psi,fn,vt,statistics)*(dpsi-dfn)/vt;
    const double dp=vela::holeDensityDerivativeEta(ni,Nv,psi,fp,vt,statistics)*(dfp-dpsi)/vt;
    for(double h:{1e-3,3e-4,1e-4}) {
        const double nd=(vela::electronDensity(ni,Nc,psi+h*dpsi,fn+h*dfn,vt,statistics)-vela::electronDensity(ni,Nc,psi-h*dpsi,fn-h*dfn,vt,statistics))/(2*h);
        const double pd=(vela::holeDensity(ni,Nv,psi+h*dpsi,fp+h*dfp,vt,statistics)-vela::holeDensity(ni,Nv,psi-h*dpsi,fp-h*dfp,vt,statistics))/(2*h);
        REQUIRE(nd==Approx(dn).epsilon(2e-7));REQUIRE(pd==Approx(dp).epsilon(2e-7));
    }
}

TEST_CASE("Metric aliases retain contributions without inventing independent unknowns", "[update_error_report]") {
    NewtonUpdateErrorReportConfig c;c.mode=UpdateErrorReportMode::Report;
    c.digits=0;c.errRef.fill(1.);
    std::vector<UpdateErrorRow> rows;
    rows.push_back({UpdateErrorEquation::Poisson,0,true,0.,1.,UpdateErrorValidity::Valid,0});
    for(std::size_t node=100;node<111;++node)
        rows.push_back({UpdateErrorEquation::Contact,node,true,0.,3.,UpdateErrorValidity::Valid,0});
    const auto r=evaluateUpdateErrorReport(c,rows,true);
    REQUIRE(r.metricEntryCount==12);
    REQUIRE(r.independentUnknownCount==2); // identity is equation-local
    REQUIRE(r.equations[3].count==11);
    REQUIRE(r.equations[3].independentUnknownCount==1);
    REQUIRE(*r.rms==Approx(std::sqrt(100./12.)));
    REQUIRE(*r.rms!=Approx(std::sqrt(10./2.))); // deduplicating changes the metric
    REQUIRE(r.rows.size()==12);
    REQUIRE(r.rows.back().row.id==110);
    rows.back().independentUnknownId.reset();
    const auto missing=evaluateUpdateErrorReport(c,rows,true);
    REQUIRE_FALSE(missing.independentUnknownCount);
    REQUIRE_FALSE(missing.equations[3].independentUnknownCount);
    REQUIRE(missing.metricEntryCount==12);
    REQUIRE(missing.rms==r.rms);
    rows.back().active=false;
    REQUIRE(evaluateUpdateErrorReport(c,rows).independentUnknownCount==2);
}

TEST_CASE("Raw numerator with applied denominator remains a distinct report combination", "[update_error_report]") {
    NewtonUpdateErrorReportConfig c;c.mode=UpdateErrorReportMode::Report;c.digits=0;c.errRef.fill(1.);
    const double old=10.,raw=-20.,applied=.01*old;
    std::vector<UpdateErrorRow> rows{{UpdateErrorEquation::Hole,918,true,applied,raw}};
    const UpdateErrorProvenance hybrid{UpdateErrorDeltaKind::RawLinearTangent,UpdateErrorDenominatorKind::AppliedAfter};
    const auto r=evaluateUpdateErrorReport(c,rows,false,hybrid);
    REQUIRE(r.provenance.delta==UpdateErrorDeltaKind::RawLinearTangent);
    REQUIRE(r.provenance.denominator==UpdateErrorDenominatorKind::AppliedAfter);
    REQUIRE(*r.rms==Approx(20./1.1));
    rows[0].delta=applied-old;
    const auto actual=evaluateUpdateErrorReport(c,rows,false,
        {UpdateErrorDeltaKind::AppliedDifference,UpdateErrorDenominatorKind::AppliedAfter});
    REQUIRE(*actual.rms==Approx(9.9/1.1));
    REQUIRE(actual.provenance.delta!=r.provenance.delta);
    REQUIRE(*actual.rms!=Approx(*r.rms));
    REQUIRE(evaluateUpdateErrorReport(c,rows).provenance.delta==UpdateErrorDeltaKind::Unspecified);
}
TEST_CASE("Passive Newton capture isolates copies and contains diagnostic failures", "[update_error_report]") {
    bool called = false;
    NewtonCaptureScope::record([&](auto&) { called = true; });
    REQUIRE_FALSE(called);
    NewtonCapture outer, inner;
    vela::VectorXd state = vela::VectorXd::Constant(3, 2.);
    {
        NewtonCaptureScope outerScope(&outer);
        NewtonCaptureScope::record([&](auto& sink) {
            NewtonCapturedIteration item; item.x = state;
            sink.iterations.push_back(std::move(item));
        });
        state[0] = 7.;
        REQUIRE(outer.iterations[0].x[0] == 2.);
        {
            NewtonCaptureScope innerScope(&inner);
            REQUIRE_NOTHROW(NewtonCaptureScope::record([](auto&) { throw 42; }));
            REQUIRE(inner.failed);
            NewtonCaptureScope::record([&](auto&) { called = true; });
            REQUIRE_FALSE(called);
        }
        NewtonCaptureScope::record([](auto& sink) { sink.truncated = true; });
        REQUIRE(outer.truncated);
        REQUIRE_FALSE(outer.failed);
    }
    NewtonCaptureScope::record([&](auto&) { called = true; });
    REQUIRE_FALSE(called);
}

TEST_CASE("Assembly population capture is scoped and copies existing operands", "[update_error_report]") {
    AssemblyPopulationCapture a,b;
    vela::VectorXd x=vela::VectorXd::Constant(3,.25),n=vela::VectorXd::Constant(1,2.),p=n;
    {
        AssemblyPopulationCaptureScope scope(&a);
        AssemblyPopulationCaptureScope::record(x,n,p);
        n[0]=7.;
        {AssemblyPopulationCaptureScope nested(&b);AssemblyPopulationCaptureScope::record(x,n,p);}
        REQUIRE(a.electrons[0]==2.);
        REQUIRE(b.electrons[0]==7.);
        {AssemblyPopulationCaptureScope disabled(nullptr);AssemblyPopulationCaptureScope::record(x,n,p);}
        REQUIRE(a.electrons[0]==2.);
        AssemblyPopulationCaptureScope::record(x,n,p);
        REQUIRE(a.electrons[0]==7.);
    }
    n[0]=9.;AssemblyPopulationCaptureScope::record(x,n,p);
    REQUIRE(a.electrons[0]==7.);
    REQUIRE(a.present);
    REQUIRE_FALSE(a.failed);
}
