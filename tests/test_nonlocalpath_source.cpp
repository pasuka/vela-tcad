#include "vela/physics/NonlocalPathSourceAssembly.h"
#include "vela/physics/NonlocalPathCellSource.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/core/PhysicalConstants.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include <cmath>
#include <numeric>

using namespace vela;
using namespace vela::experimental;
using Catch::Approx;
namespace {
constexpr Real q=1.602192e-19,F=1e8,L=1e-8,H=2e-9;
NonlocalPathBandToBandConfig config() {
    NonlocalPathBandToBandConfig c;c.regions={0};c.constants=TunnelingConstants::SentaurusT202203;
    c.ab=IndirectABParameters{4e20,1.9e9,q*F*L};c.phononEnergy_J=.037*q;c.maximumPathLength_m=4*L;return c;
}
TunnelingGeometryMesh mesh(Index nx=6) {
    TunnelingGeometryMesh m;m.revision=1;
    for (Index i=0;i<=nx;++i) for (Real y:{0.,H}) m.positions_m.emplace_back(3*L*i/nx,y);
    for (Index i=0;i<nx;++i) {
        m.cells.push_back({{2*i,2*i+2,2*i+1},0,0,true});
        m.cells.push_back({{2*i+2,2*i+3,2*i+1},0,0,true});
        m.boundaries.push_back({{2*i,2*i+2},TunnelingBoundaryKind::Reflecting});
        m.boundaries.push_back({{2*i+1,2*i+3},TunnelingBoundaryKind::Reflecting});
    }
    m.boundaries.push_back({{0,1},TunnelingBoundaryKind::Reflecting});
    m.boundaries.push_back({{2*nx,2*nx+1},TunnelingBoundaryKind::Reflecting});return m;
}
TunnelingBandState state(const TunnelingGeometryMesh& m) {
    TunnelingBandState s;s.identity.mesh=m.revision;s.identity.trial=7;
    for (const auto& p:m.positions_m) {
        TunnelingBandNode n;n.valence_J=-q*F*p.x();n.conduction_J=n.valence_J+q*F*L;
        n.electronFermi_J=-100*q;n.holeFermi_J=100*q;
        n.conductionDOS_per_m3=n.valenceDOS_per_m3=1e25;
        n.electrons_per_m3=n.holes_per_m3=1e16;s.nodes.push_back(n);
    }
    return s;
}
TunnelingWindowInput windows(const TunnelingBandState& s) {
    return {s.identity,std::vector<TunnelingSeedWindow>(s.nodes.size(),{100*q,100*q})};
}
long double sum(const std::vector<Real>& v) {return std::accumulate(v.begin(),v.end(),0.L);}
}
TEST_CASE("NLP N4 first crossing keeps uphill identities without duplicate collection", "[nonlocalpath][source][cell]") {
    const auto w=tunnelingFirstCrossingWeights({4.,2.,3.,0.},0.,4.);
    REQUIRE(w==std::vector<Real>{.25,.25,.125,.375});
    REQUIRE_THROWS(tunnelingFirstCrossingWeights({4.,2.},1.,1.));
}
TEST_CASE("NLP N4 fraction primitive has analytic constant gap limit", "[nonlocalpath][source][cell]") {
    // z=x/d, equal masses, gap=2: k=2/hbar. Each half integrates
    // sqrt(z) and 1/sqrt(z) exactly, including the integrable endpoint.
    const auto a=tunnelingFractionIntegrals({{0.,0.,2.,300.},{3.,-2.,0.,300.}},0.,1.,1.,1.);
    REQUIRE(a.valenceAction==Approx(4*std::pow(.5,1.5)).epsilon(1e-13));
    REQUIRE(a.conductionAction==Approx(a.valenceAction).epsilon(1e-13));
    REQUIRE(a.valenceInverse_m2==Approx(3*std::sqrt(.5)).epsilon(1e-13));
    REQUIRE(a.conductionInverse_m2==Approx(a.valenceInverse_m2).epsilon(1e-13));
}
TEST_CASE("NLP N4 cell source is conservative signed gauge invariant and coverage guarded", "[nonlocalpath][source][cell]") {
    auto m=mesh();auto s=state(m);auto c=config();
    m.boundaries.back().kind=TunnelingBoundaryKind::Contact;
    s.nodes[s.nodes.size()-1].contact=s.nodes[s.nodes.size()-2].contact=true;
    TunnelingCellIntegrationPolicy policy;policy.variationReference_per_m2_s=1e100;
    const auto a=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE(a.resolved);REQUIRE_NOTHROW(a.requireResolved(s));
    REQUIRE(a.records.size()==m.cells.size());REQUIRE(a.sources.balanceError()<1e-12);
    REQUIRE(sum(a.sources.electronIntegrals_per_m_s)>0);
    REQUIRE(a.sources.electronIntegrals_per_m_s!=a.sources.holeIntegrals_per_m_s);
    REQUIRE(a.contacts.contactElectron_per_m_s>0);
    auto partial=a;partial.records.pop_back();REQUIRE_THROWS(partial.requireResolved(s));
    for(auto& n:s.nodes) {n.conduction_J+=2*q;n.valence_J+=2*q;n.electronFermi_J+=2*q;n.holeFermi_J+=2*q;}
    const auto gauge=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE(static_cast<Real>(sum(gauge.sources.electronIntegrals_per_m_s)/sum(a.sources.electronIntegrals_per_m_s))==Approx(1).epsilon(1e-11));
    for(auto& n:s.nodes) std::swap(n.electronFermi_J,n.holeFermi_J);
    const auto negative=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE(static_cast<Real>(sum(negative.sources.electronIntegrals_per_m_s)/sum(a.sources.electronIntegrals_per_m_s))==Approx(-1).epsilon(1e-11));
    REQUIRE(negative.sources.balanceError()<1e-12);
    policy.maximumDepth=0;policy.minimumWidth_kT=-1;
    REQUIRE_THROWS(assembleNonlocalPathCellSources(m,s,c,policy));
}
TEST_CASE("NLP N4 independent uniform law integrates and transfers between different endpoints", "[nonlocalpath][source]") {
    for (const auto policy:{TunnelingBoxPolicy::Barycentric,TunnelingBoxPolicy::Delaunay}) {
        const auto m=mesh();const auto s=state(m);const auto c=config();
        const auto a=assembleNonlocalPathSources(m,s,c,windows(s),policy);
        REQUIRE(a.resolved);REQUIRE_NOTHROW(a.requireResolved(s));
        REQUIRE(a.sources.balanceError()<1e-12);
        const Real G=c.ab->A_per_m3_s*std::pow(F/100,2.5)*std::exp(-c.ab->B_V_per_m/F);
        long double expected=0,moment=0;
        for (Index i=0;i<s.nodes.size();++i) {
            if (m.positions_m[i].x()<=2*L*(1+1e-14)) expected+=a.boxes.seedAreas_m2[i]*G;
            moment+=static_cast<long double>(m.positions_m[i].x())*(a.sources.electronIntegrals_per_m_s[i]-a.sources.holeIntegrals_per_m_s[i]);
        }
        REQUIRE(static_cast<Real>(sum(a.sources.electronIntegrals_per_m_s)/expected)==Approx(1).epsilon(1e-8));
        REQUIRE(static_cast<Real>(moment/(expected*L))==Approx(1).epsilon(1e-8));
        REQUIRE(a.sources.electronIntegrals_per_m_s!=a.sources.holeIntegrals_per_m_s);
        REQUIRE(static_cast<Real>(sum(a.boxes.fieldAreas_m2))==Approx(3*L*H).epsilon(1e-13));
    }
}
TEST_CASE("NLP N4 cell diagnostics reject incomplete chords and energy resource exhaustion", "[nonlocalpath][source][cell]") {
    const auto m=mesh();auto s=state(m);auto c=config();TunnelingCellIntegrationPolicy policy;
    policy.variationReference_per_m2_s=1e100;
    c.maximumPathLength_m=L/100;
    auto narrow=s;for(auto& n:narrow.nodes) n.conduction_J=n.valence_J+q*F*L/100;
    const auto shortChord=assembleNonlocalPathCellSources(m,narrow,c,policy);
    REQUIRE_FALSE(shortChord.resolved);REQUIRE_THROWS(shortChord.requireResolved(s));
    // A linear barrier whose gap exceeds the entire chord's band drop has
    // no conduction turn even before the length cap: an exact zero proof.
    const auto emptyShort=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE(emptyShort.resolved);REQUIRE_NOTHROW(emptyShort.requireResolved(s));
    for(Real v:emptyShort.sources.electronIntegrals_per_m_s) REQUIRE(v==0);
    c=config();c.maximumIntegrationIntervals=2;
    const auto exhausted=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE_FALSE(exhausted.resolved);REQUIRE_THROWS(exhausted.requireResolved(s));
    c=config();for(auto& n:s.nodes) n.electronFermi_J=n.holeFermi_J=0;
    const auto equilibrium=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE(equilibrium.resolved);REQUIRE_NOTHROW(equilibrium.requireResolved(s));
    for(Real v:equilibrium.sources.electronIntegrals_per_m_s) REQUIRE(v==0);
    for(Real v:equilibrium.sources.holeIntegrals_per_m_s) REQUIRE(v==0);
}
TEST_CASE("NLP N4 spatial weighting retains representable subnormal events", "[nonlocalpath][source][cell]") {
    auto m=mesh();for(auto& p:m.positions_m) {p.x()*=100;p.y()*=1e-7;}
    auto s=state(m);auto c=config();c.maximumPathLength_m=1e-5;
    for(Index i=0;i<s.nodes.size();++i) {
        s.nodes[i].valence_J=-q*2.52e6*m.positions_m[i].x();
        s.nodes[i].conduction_J=s.nodes[i].valence_J+q*F*L;
    }
    TunnelingCellIntegrationPolicy policy;policy.variationReference_per_m2_s=1e100;
    const auto a=assembleNonlocalPathCellSources(m,s,c,policy);
    REQUIRE(a.resolved);REQUIRE_NOTHROW(a.requireResolved(s));
    REQUIRE(sum(a.sources.electronIntegrals_per_m_s)>0);
    REQUIRE(sum(a.sources.holeIntegrals_per_m_s)>0);
}
TEST_CASE("NLP N4 signed equilibrium and contact bookkeeping retain both carriers", "[nonlocalpath][source]") {
    const auto m=mesh();auto s=state(m);s.nodes[0].contact=true;s.nodes[1].contact=true;s.nodes[4].contact=true;
    auto a=assembleNonlocalPathSources(m,s,config(),windows(s),TunnelingBoxPolicy::Delaunay);
    REQUIRE(a.contacts.contactHole_per_m_s>0);REQUIRE(a.contacts.contactElectron_per_m_s>0);
    REQUIRE(a.contacts.freeElectron_per_m_s+a.contacts.contactElectron_per_m_s==Approx(static_cast<Real>(sum(a.sources.electronIntegrals_per_m_s))));
    REQUIRE(a.contacts.freeHole_per_m_s+a.contacts.contactHole_per_m_s==Approx(static_cast<Real>(sum(a.sources.holeIntegrals_per_m_s))));
    for (auto& n:s.nodes) std::swap(n.electronFermi_J,n.holeFermi_J);
    auto b=assembleNonlocalPathSources(m,s,config(),windows(s),TunnelingBoxPolicy::Delaunay);
    for (Index i=0;i<s.nodes.size();++i) {
        REQUIRE(b.sources.electronIntegrals_per_m_s[i]==-a.sources.electronIntegrals_per_m_s[i]);
        REQUIRE(b.sources.holeIntegrals_per_m_s[i]==-a.sources.holeIntegrals_per_m_s[i]);
        s.nodes[i].electronFermi_J=s.nodes[i].holeFermi_J=0;
    }
    const auto zero=assembleNonlocalPathSources(m,s,config(),windows(s),TunnelingBoxPolicy::Delaunay);
    REQUIRE(sum(zero.sources.electronIntegrals_per_m_s)==0);REQUIRE(sum(zero.sources.holeIntegrals_per_m_s)==0);
}
TEST_CASE("NLP N4 fields integrate back to source and carry no evaluated derivatives", "[nonlocalpath][source]") {
    const auto m=mesh();const auto s=state(m);
    const auto a=assembleNonlocalPathSources(m,s,config(),windows(s),TunnelingBoxPolicy::Delaunay);
    for (bool electron:{true,false}) {
        const auto rates=a.sources.volumeRates(electron,a.boxes.fieldAreas_m2);long double integral=0;
        for (Index i=0;i<s.nodes.size();++i) integral+=static_cast<long double>(rates[i])*a.boxes.fieldAreas_m2[i];
        REQUIRE(static_cast<Real>(integral)==Approx(static_cast<Real>(sum(electron?a.sources.electronIntegrals_per_m_s:a.sources.holeIntegrals_per_m_s))).epsilon(1e-13));
    }
    REQUIRE_FALSE(a.sources.evaluatedDerivativeMode.has_value());REQUIRE(a.sources.derivatives.empty());
    const auto audit=tunnelingSourceAudit(a,s);REQUIRE(audit["resolved"].get<bool>());
    REQUIRE_FALSE(audit["native_window_mapping_qualified"].get<bool>());
}
TEST_CASE("NLP N4 missing windows and incomplete contact paths fail closed", "[nonlocalpath][source]") {
    auto m=mesh();const auto s=state(m);auto w=windows(s);++w.identity.trial;
    REQUIRE_THROWS(assembleNonlocalPathSources(m,s,config(),w,TunnelingBoxPolicy::Barycentric));
    w=windows(s);w.nodes.pop_back();REQUIRE_THROWS(assembleNonlocalPathSources(m,s,config(),w,TunnelingBoxPolicy::Barycentric));
    w=windows(s);w.nodes[0].conductionAvailableEnergy_J=-q;
    REQUIRE_THROWS(assembleNonlocalPathSources(m,s,config(),w,TunnelingBoxPolicy::Barycentric));
    m.boundaries.back().kind=TunnelingBoundaryKind::Contact;
    auto a=assembleNonlocalPathSources(m,s,config(),windows(s),TunnelingBoxPolicy::Barycentric);
    REQUIRE_FALSE(a.resolved);REQUIRE_THROWS(a.requireResolved(s));
    REQUIRE(a.sources.balanceError()<1e-12); // Partial ledger is conservative, NOT complete.
    auto overlap=config();overlap.conductionOffset_J=-2*q*F*L;
    const auto unsupported=assembleNonlocalPathSources(m,s,overlap,windows(s),TunnelingBoxPolicy::Barycentric);
    REQUIRE_FALSE(unsupported.resolved);REQUIRE_THROWS(unsupported.requireResolved(s));
}
TEST_CASE("NLP N4 source box policy excludes oxide and respects selected regions", "[nonlocalpath][source]") {
    auto m=mesh();auto c=config();
    const auto base=tunnelingSourceBoxes(m,c,TunnelingBoxPolicy::Delaunay);
    for (Index i=0;i<m.cells.size();++i) m.cells[i].region=i<6?0:1;
    auto restricted=tunnelingSourceBoxes(m,c,TunnelingBoxPolicy::Delaunay);
    REQUIRE(sum(restricted.seedAreas_m2)==Approx(static_cast<Real>(sum(base.seedAreas_m2)/2)).epsilon(1e-13));
    REQUIRE(restricted.fieldAreas_m2==base.fieldAreas_m2);
    for (auto& cell:m.cells) if (cell.region==1) {cell.material=1;cell.semiconductor=false;}
    auto oxide=tunnelingSourceBoxes(m,c,TunnelingBoxPolicy::Delaunay);
    REQUIRE(sum(oxide.fieldAreas_m2)==Approx(static_cast<Real>(sum(base.fieldAreas_m2)/2)).epsilon(1e-13));
    REQUIRE(oxide.fieldAreas_m2.back()==0);
}
TEST_CASE("NLP N4 actual residual uses separate source rows SI scaling and contact replacement", "[nonlocalpath][source][coupled]") {
    for (bool tcad:{false,true}) {
        const auto m=mesh();auto s=state(m);auto c=config();c.ab->A_per_m3_s=1e35;
        UnitScalingConfig units(tcad?UnitScalingMode::UnitScaling:UnitScalingMode::LegacySI);
        const auto u=units.unitSystem();DeviceMesh dm;
        for (Index i=0;i<m.positions_m.size();++i) {Node n;n.id=i;n.x=u.metersToInternalLength(m.positions_m[i].x());n.y=u.metersToInternalLength(m.positions_m[i].y());dm.addNode(n);}
        for (Index i=0;i<m.cells.size();++i) {Cell cell;cell.id=i;cell.type=CellType::Tri3;cell.region_id=0;cell.node_ids.assign(m.cells[i].nodes.begin(),m.cells[i].nodes.end());dm.addCell(cell);}
        Region region;region.id=0;region.name="si";region.material="Si";for(Index i=0;i<m.cells.size();++i)region.cell_ids.push_back(i);dm.addRegion(region);dm.buildEdges();
        MaterialDatabase db(units);DopingModel doping(dm.numNodes());DDScalingSpec scale;
        if (tcad) {scale.enabled=true;scale.unitSystem=u;scale.V0=.025;scale.C0=1e10;scale.D0=1;scale.chargeAreaFactor=u.chargeAreaFactor();scale.fieldFromCoordinateDeltaFactor=u.fieldFromCoordinateDeltaFactor();}
        CoupledDDAssembler dd(dm,db,doping,constants::Vt_300,mobilityModelConfig("constant"),recombinationModelConfig({"none"}),{},{},{},{},scale);
        VectorXd x=VectorXd::Zero(3*dm.numNodes());
        for (Index i=0;i<dm.numNodes();++i) {x(dm.numNodes()+i)=.01/(tcad?scale.V0:1);x(2*dm.numNodes()+i)=-.01/(tcad?scale.V0:1);}
        const auto n=dd.electronDensity(x),p=dd.holeDensity(x);
        for (Index i=0;i<s.nodes.size();++i) {s.nodes[i].electronFermi_J=-q*.01;s.nodes[i].holeFermi_J=q*.01;s.nodes[i].electrons_per_m3=n(i)*u.concentrationM3PerInternal();s.nodes[i].holes_per_m3=p(i)*u.concentrationM3PerInternal();}
        TunnelingCellIntegrationPolicy integration;integration.variationReference_per_m2_s=1e100;
        const auto a=assembleNonlocalPathCellSources(m,s,c,integration);REQUIRE(a.resolved);
        CoupledDDBoundaryConditions bc;bc.phin[4]=x(dm.numNodes()+4);bc.phip[0]=x(2*dm.numNodes());
        const auto base=dd.residual(x,bc),actual=dd.experimentalNonlocalResidual(x,bc,a,s);
        // Independent SI particle-line-flux -> internal-current conversion.
        const Real factor=1/(u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*dd.continuityResidualScale());
        for (Index i=0;i<dm.numNodes();++i) {
            REQUIRE(actual(i)==base(i));
            for (bool electron:{true,false}) {
                const Index row=(electron?dm.numNodes():2*dm.numNodes())+i;
                const Real event=(electron?a.sources.electronIntegrals_per_m_s:a.sources.holeIntegrals_per_m_s)[i];
                const bool constrained=electron?bc.phin.contains(i):bc.phip.contains(i);
                if (constrained || event==0) REQUIRE(actual(row)==base(row));
                else REQUIRE((actual(row)-base(row))/(-event*factor)==Approx(1).epsilon(1e-8));
            }
        }
        auto stale=s;stale.nodes[0].electronFermi_J+=q*.001;REQUIRE_THROWS(dd.experimentalNonlocalResidual(x,bc,a,stale));
        CoupledDDFeedbackStateSubstitution frozen;
        frozen.replaceElectronDensity=frozen.replaceHoleDensity=true;
        frozen.electronDensity=n*1.02;frozen.holeDensity=p*1.02;
        auto fs=s;for(Index i=0;i<fs.nodes.size();++i) {
            fs.nodes[i].electrons_per_m3=frozen.electronDensity(i)*u.concentrationM3PerInternal();
            fs.nodes[i].holes_per_m3=frozen.holeDensity(i)*u.concentrationM3PerInternal();
        }
        const auto fa=assembleNonlocalPathCellSources(m,fs,c,integration);
        REQUIRE_THROWS(dd.experimentalNonlocalResidual(x,bc,fa,fs));
        const auto fb=dd.feedbackSubstitutionResidual(x,bc,frozen);
        const auto fr=dd.experimentalNonlocalFrozenResidual(x,bc,fa,fs,frozen);
        for(Index i=0;i<dm.numNodes();++i) {
            REQUIRE(fr(i)==fb(i));
            if(!bc.phin.contains(i) && fa.sources.electronIntegrals_per_m_s[i]!=0)
                REQUIRE((fr(dm.numNodes()+i)-fb(dm.numNodes()+i))/(-fa.sources.electronIntegrals_per_m_s[i]*factor)==Approx(1).epsilon(1e-8));
        }
        auto invalidFrozen=frozen;invalidFrozen.replaceHoleDensity=false;
        REQUIRE_THROWS(dd.experimentalNonlocalFrozenResidual(x,bc,fa,fs,invalidFrozen));
        invalidFrozen=frozen;invalidFrozen.electronDensity(0)*=2;
        REQUIRE_THROWS(dd.experimentalNonlocalFrozenResidual(x,bc,fa,fs,invalidFrozen));
        auto hot=fs;++hot.identity.temperature;
        for(auto& node:hot.nodes) node.temperature_K=400;
        const auto hotAssembly=assembleNonlocalPathCellSources(m,hot,c,integration);
        REQUIRE(hotAssembly.resolved);
        REQUIRE_THROWS(dd.experimentalNonlocalFrozenResidual(x,bc,hotAssembly,hot,frozen));
        ++stale.identity.trial;REQUIRE_THROWS(dd.experimentalNonlocalResidual(x,bc,a,stale));
        auto wrongMesh=a;wrongMesh.sourceMesh.positions_m[0].x()+=1e-12;
        REQUIRE_THROWS(dd.experimentalNonlocalResidual(x,bc,wrongMesh,s));
        wrongMesh=a;wrongMesh.sourceMesh.cells[0].region=1;
        REQUIRE_THROWS(dd.experimentalNonlocalResidual(x,bc,wrongMesh,s));
    }
}

TEST_CASE("NLP N4 repeated mesh-line crossings preserve vertex support under curved bands", "[nonlocalpath][source][geometry]") {
    auto m=mesh(32);auto s=state(m);auto c=config();
    for (Index i=0;i<s.nodes.size();++i) {
        const Real x=m.positions_m[i].x();
        s.nodes[i].valence_J=-q*F*(x+.2*x*x/L);
        s.nodes[i].conduction_J=s.nodes[i].valence_J+q*F*L;
    }
    const auto paths=NonlocalPathGeometry(m).traceFromNodes(s,c);
    REQUIRE(paths.size()==s.nodes.size());
    for (const auto& p:paths) {
        REQUIRE(p.path.termination!=TunnelingTermination::NumericalFailure);
        for (const auto& point:p.path.points) {
            Point2 rebuilt=Point2::Zero();for (const auto& [node,w]:point.support.weights)rebuilt+=w*m.positions_m[node];
            REQUIRE((rebuilt-point.position_m).norm()<128*std::numeric_limits<Real>::epsilon()*3*L);
        }
    }
}

TEST_CASE("Frozen SRH populations close the actual continuity ledger", "[nonlocalpath][source][coupled][srh]") {
    const auto tm=mesh();DeviceMesh dm;
    for(Index i=0;i<tm.positions_m.size();++i) {
        Node n;n.id=i;n.x=tm.positions_m[i].x();n.y=tm.positions_m[i].y();dm.addNode(n);
    }
    for(Index i=0;i<tm.cells.size();++i) {
        Cell c;c.id=i;c.type=CellType::Tri3;c.region_id=0;
        c.node_ids.assign(tm.cells[i].nodes.begin(),tm.cells[i].nodes.end());dm.addCell(c);
    }
    Region region;region.id=0;region.material="Si";region.name="si";
    for(Index i=0;i<dm.numCells();++i)region.cell_ids.push_back(i);
    dm.addRegion(region);dm.buildEdges();MaterialDatabase db;DopingModel doping(dm.numNodes());
    for(const std::string model:{"boltzmann","fermi_dirac"}) {
      for(const std::string coupling:{"quantum","sentaurus_default"}) {
        auto recombination=recombinationModelConfig({"srh"});
        recombination.srhDopingDependence.densityCoupling=coupling;
        CoupledDDAssembler dd(dm,db,doping,constants::Vt_300,mobilityModelConfig("constant"),
            recombination,{},{},{},{},{},{},CarrierStatisticsConfig{model});
        const Index N=dm.numNodes();VectorXd x=VectorXd::Zero(3*N);
        x.head(N).setConstant(.05);x.segment(N,N).setConstant(.01);x.tail(N).setConstant(-.03);
        CoupledDDFeedbackStateSubstitution frozen;frozen.replaceElectronDensity=frozen.replaceHoleDensity=true;
        frozen.electronDensity=7*dd.electronDensity(x);frozen.holeDensity=3*dd.holeDensity(x);
        CoupledDDBoundaryConditions bc;bc.phin[0]=.01;bc.phip[0]=-.03;
        const auto actual=dd.feedbackSubstitutionResidual(x,bc,frozen);
        const auto terms=dd.feedbackSubstitutionCarrierContinuityTermDiagnostics(x,bc,frozen);
        for(Index i=1;i<N;++i) {
            REQUIRE(terms[i].electronRecombination!=0);
            REQUIRE(actual(N+i)==Approx(terms[i].electronResidual).epsilon(2e-14));
            REQUIRE(actual(2*N+i)==Approx(terms[i].holeResidual).epsilon(2e-14));
        }
        REQUIRE(actual(N)==Approx(terms[0].electronResidual).margin(1e-18));
        REQUIRE(actual(2*N)==Approx(terms[0].holeResidual).margin(1e-18));
      }
    }
}

TEST_CASE("Frozen Poisson ledger retains absolute edge support and supplied charge", "[nonlocalpath][source][coupled][poisson]") {
    DeviceMesh dm;
    for(Index i=0;i<3;++i) {Node n;n.id=i;n.x=i==1?1e-6:0.;n.y=i==2?1e-6:0.;dm.addNode(n);}
    Cell cell;cell.id=0;cell.type=CellType::Tri3;cell.region_id=0;cell.node_ids={0,1,2};dm.addCell(cell);
    Region region;region.id=0;region.material="Si";region.name="si";region.cell_ids={0};
    dm.addRegion(region);dm.buildEdges();MaterialDatabase db;DopingModel doping(3);
    CoupledDDAssembler dd(dm,db,doping,constants::Vt_300,mobilityModelConfig("constant"),recombinationModelConfig({"none"}));
    VectorXd x=VectorXd::Zero(9);x(1)=.1;x(2)=-.1;
    CoupledDDFeedbackStateSubstitution frozen;frozen.replaceElectronDensity=frozen.replaceHoleDensity=true;
    frozen.electronDensity=dd.electronDensity(x)*4;frozen.holeDensity=dd.holeDensity(x)*2;
    CoupledDDBoundaryConditions bc;bc.psi[2]=-.1;
    const auto check=[&]<class Assembler>(Assembler& a) {
        if constexpr(requires{a.feedbackSubstitutionPoissonTermDiagnostics(x,bc,frozen);}) {
            const auto terms=a.feedbackSubstitutionPoissonTermDiagnostics(x,bc,frozen);
            const auto actual=a.feedbackSubstitutionResidual(x,bc,frozen);
            for(Index i=0;i<3;++i) {
                const auto& t=terms[i];
                REQUIRE(t.productionResidual==actual(i));
                REQUIRE(t.electronCharge==Approx(4*a.poissonTermDiagnostics(x,bc)[i].electronCharge).epsilon(2e-14));
                REQUIRE(t.holeCharge==Approx(2*a.poissonTermDiagnostics(x,bc)[i].holeCharge).epsilon(2e-14));
                REQUIRE(t.dielectricFluxAbsSum>=std::abs(t.dielectricFlux));
                REQUIRE(std::abs(t.closureError)<2e-14*(t.dielectricFluxAbsSum+std::abs(t.electronCharge)+std::abs(t.holeCharge)));
            }
            // Opposing dielectric edges cancel at node 0, while their absolute
            // support remains nonzero. A net flux cannot supply this denominator.
            REQUIRE(terms[0].dielectricFluxAbsSum>1e-13);
            REQUIRE(std::abs(terms[0].dielectricFlux)<1e-14*terms[0].dielectricFluxAbsSum);
            auto bad=frozen;bad.electronDensity(0)=-1;
            REQUIRE_THROWS(a.feedbackSubstitutionPoissonTermDiagnostics(x,bc,bad));
        } else FAIL("Poisson frozen-state and absolute-edge ledger is missing");
    };
    check(dd);
}

TEST_CASE("Scaled Poisson observer preserves the actual endpoint subtraction", "[nonlocalpath][source][coupled][poisson]") {
    DeviceMesh dm;
    for(Index i=0;i<3;++i) {Node n;n.id=i;n.x=i==1?1e-6:0.;n.y=i==2?1e-6:0.;dm.addNode(n);}
    Cell cell;cell.id=0;cell.type=CellType::Tri3;cell.region_id=0;cell.node_ids={0,1,2};dm.addCell(cell);
    Region region;region.id=0;region.material="Si";region.name="si";region.cell_ids={0};
    dm.addRegion(region);dm.buildEdges();MaterialDatabase db;DopingModel doping(3);
    DDScalingSpec scaling;scaling.enabled=true;scaling.V0=.025;
    CoupledDDAssembler dd(dm,db,doping,constants::Vt_300,mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),{},{},{},{},scaling);
    VectorXd x=VectorXd::Zero(9);x(0)=100.;x(1)=100.+1e-10;x(2)=100.-2e-10;
    CoupledDDFeedbackStateSubstitution frozen;frozen.replaceElectronDensity=frozen.replaceHoleDensity=true;
    frozen.electronDensity=VectorXd::Constant(3,1e-100);frozen.holeDensity=frozen.electronDensity;
    const auto actual=dd.feedbackSubstitutionResidual(x,{},frozen);
    const auto terms=dd.feedbackSubstitutionPoissonTermDiagnostics(x,{},frozen);
    for(Index i=0;i<3;++i) {
        REQUIRE(terms[i].dielectricFluxAbsSum>0);
        REQUIRE(std::abs(terms[i].dielectricFlux-actual(i))<5e-15*terms[i].dielectricFluxAbsSum);
        REQUIRE(std::abs(terms[i].closureError)<5e-15*terms[i].dielectricFluxAbsSum);
    }
}
