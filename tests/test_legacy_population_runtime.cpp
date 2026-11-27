#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/solver/PopulationState.h"
#include "vela/solver/SlotPopulationNumerics.h"
#include "vela/solver/PopulationWriteback.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/simulation/PseudoArclength.h"
#include <nlohmann/json.hpp>
using namespace vela;

TEST_CASE("Boltzmann writeback preserves common potential equilibrium and protection ownership", "[legacy_population][writeback]") {
    const double vt=.025852,psi=10.,n=1e20,p=1e10;
    VectorXd old(3);old<<psi,n,p;
    const double dp=.6*(std::nextafter(psi,INFINITY)-psi);
    VectorXd raw(3);raw<<dp,n*std::expm1(dp/vt),p*std::expm1(-dp/vt);
    PopulationUpdateConfig cfg;
    auto base=makePopulationTrial(old,raw,1.,{true},{true},cfg);
    auto fixed=compensateBoltzmannPsiWriteback(old,raw,base,1.,vt,{true});
    auto gap=[&](const VectorXd& x) {
        const long double dpsi=static_cast<long double>(x[0])-psi;
        const long double en=dpsi-vt*std::log(static_cast<long double>(x[1])/n);
        const long double ep=dpsi+vt*std::log(static_cast<long double>(x[2])/p);
        return std::max(std::abs(en),std::abs(ep));
    };
    CHECK(gap(fixed.state)<.02L*gap(base.state));
    CHECK(fixed.state[0]==base.state[0]);CHECK(fixed.alpha==base.alpha);
    CHECK(compensateBoltzmannPsiWriteback(old,raw,base,1.,vt,{false}).state==base.state);
    const auto zero=makePopulationTrial(old,VectorXd::Zero(3),1.,{true},{true},cfg);
    CHECK(compensateBoltzmannPsiWriteback(old,VectorXd::Zero(3),zero,1.,vt,{true}).state==old);
    raw[1]=-2*n;
    const auto projected=makePopulationTrial(old,raw,1.,{true},{true},cfg);
    const auto projectedFixed=compensateBoltzmannPsiWriteback(old,raw,projected,1.,vt,{true});
    CHECK(projectedFixed.state[1]==projected.state[1]);
    CHECK(projectedFixed.protectedElectrons==projected.protectedElectrons);
    cfg.protection="fraction_to_boundary";
    const auto bounded=makePopulationTrial(old,raw,1.,{true},{true},cfg);
    const auto boundedFixed=compensateBoltzmannPsiWriteback(old,raw,bounded,1.,vt,{true});
    CHECK(boundedFixed.alpha==bounded.alpha);CHECK(boundedFixed.state[1]>0.);
    CHECK_THROWS(compensateBoltzmannPsiWriteback(old,raw,bounded,0.,vt,{true}));
}

TEST_CASE("Psi writeback compensation has explicit model and coordinate guards", "[legacy_population][writeback]") {
    nlohmann::json deck={{"state_variables","potential_density"},{"coupled_operator","legacy"},
        {"carrier_statistics",{{"model","boltzmann"}}}};
    CHECK(newtonConfigFromJson(deck).populationWriteback=="independent");
    deck["population_writeback"]="boltzmann_psi_compensated";
    CHECK(newtonConfigFromJson(deck).populationWriteback=="boltzmann_psi_compensated");
    deck["carrier_statistics"]["model"]="fermi_dirac";CHECK_THROWS(newtonConfigFromJson(deck));
    deck["carrier_statistics"]["model"]="boltzmann";
    deck["state_variables"]="potential_quasi_fermi";CHECK_THROWS(newtonConfigFromJson(deck));
    deck["state_variables"]="potential_density";
    deck["population_writeback"]="typo";CHECK_THROWS(newtonConfigFromJson(deck));
}

namespace {
DeviceMesh populationMesh() {
    DeviceMesh m;
    for(Index i=0;i<4;++i) {Node n;n.id=i;n.x=(i==1?1.:.5)*1e-6;n.y=(i==2?.3:i==3?-.9:0.)*1e-6;if(i==0)n.x=0.;m.addNode(n);}
    for(Index i=0;i<2;++i) {Cell c;c.id=i;c.type=CellType::Tri3;c.region_id=0;c.node_ids=i?std::vector<Index>{1,0,3}:std::vector<Index>{0,1,2};m.addCell(c);}
    Region r;r.id=0;r.name=r.material="Si";r.cell_ids={0,1};m.addRegion(r);
    Contact s;s.name="source";s.node_ids={0};m.addContact(s);
    Contact d;d.id=1;d.name="drain";d.node_ids={1};m.addContact(d);
    m.buildEdges();BoxGeometryBuilder::Options o;
    o.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
    o.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::DelaunayTransfer;m.buildBoxGeometry(o);return m;
}
}

TEST_CASE("Prescribed population seed projects only boundary coordinates", "[legacy_population][contact_preprojection]") {
    const auto model = GENERATE(std::string("boltzmann"), std::string("fermi_dirac"));
    auto mesh = populationMesh();
    MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node)
        doping.setNodeDoping(node, 1e24, 0.);
    const Real vt = constants::Vt_300;
    DDScalingSpec scaling;
    scaling.enabled = true;
    scaling.V0 = vt;
    scaling.C0 = 1e23;
    scaling.D0 = .01;
    scaling.stateVariables = "potential_density";
    CoupledDDAssembler assembler(mesh, materials, doping, vt,
        mobilityModelConfig("constant"), recombinationModelConfig({}), bandgapNarrowingConfig("none"),
        {}, {}, {}, scaling, {}, {.model = model, .transport = model == "boltzmann" ? "generalized_einstein" : "effective_band_sg"});
    const int count = static_cast<int>(mesh.numNodes());
    VectorXd parent = assembler.packPopulation(VectorXd::Zero(count),
        VectorXd::Constant(count, 1e23), VectorXd::Constant(count, 1e10));
    CoupledDDBoundaryConditions boundary;
    boundary.psi[0] = .47 / vt;
    boundary.phin[0] = boundary.phip[0] = .1 / vt;
    boundary.psi[2] = .2 / vt;
    VectorXd seed = parent;
    assembler.projectPopulationBoundarySeed(seed, boundary);
    const auto residual = assembler.residual(seed, boundary);
    CHECK(seed[0] == boundary.psi[0]);
    CHECK(seed[2] == boundary.psi[2]);
    CHECK(std::abs(residual[count] * vt) < 5e-14);
    CHECK(std::abs(residual[2 * count] * vt) < 5e-14);
    for (const int node : {1, 2, 3}) {
        CHECK(seed[count + node] == parent[count + node]);
        CHECK(seed[2 * count + node] == parent[2 * count + node]);
    }
    CHECK(seed[1] == parent[1]);
    CHECK(seed[3] == parent[3]);
    CHECK(parent[0] == 0.);
    CHECK_FALSE(newtonConfigFromJson(nlohmann::json::object()).contactPreprojection);
    CHECK(newtonConfigFromJson({{"contact_preprojection", true}}).contactPreprojection);
}

TEST_CASE("Legacy population coordinates retain classical SLOT models and physical Jacobian", "[legacy_population]") {
    const auto model=GENERATE(std::string("boltzmann"),std::string("fermi_dirac"),
                              std::string("fermi_dirac_joyce_dixon"));
    const auto mobilityName=GENERATE(std::string("constant"),std::string("phumob"));
    const bool scaled=GENERATE(false,true);
    INFO(model);INFO(mobilityName);INFO(scaled);
    auto mesh=populationMesh();MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e24*(i+1),2e23);
    auto mobility=mobilityModelConfig(mobilityName);
    mobility.dopingConcentrationBasis="total_impurity";
    if(mobilityName=="phumob")mobility.edgeAveraging="element_box_phumob";
    auto rec=recombinationModelConfig({"srh","auger"});
    const bool isFermi=model!="boltzmann";
    auto bgn=bandgapNarrowingConfig(isFermi?"bennett_wilson":"old_slotboom");
    CarrierStatisticsConfig statistics{.model = model, .transport = isFermi ? "effective_band_sg" : "generalized_einstein"};
    DDScalingSpec q;q.enabled=scaled;q.V0=constants::Vt_300;q.C0=1e23;q.D0=.01;
    q.regionResolvedInterfaceAssembly.transportEdgeGeometry="element_box";
    auto p=q;p.stateVariables="potential_density";
    CoupledDDAssembler qa(mesh,db,doping,constants::Vt_300,mobility,rec,bgn,{},{},{},q,{},statistics);
    CoupledDDAssembler pa(mesh,db,doping,constants::Vt_300,mobility,rec,bgn,{},{},{},p,{},statistics);
    const double scale=scaled?q.V0:1.;
    CoupledDDState state;state.psi=VectorXd::LinSpaced(4,.01,.04)/scale;
    state.phin=VectorXd::LinSpaced(4,-.005,.008)/scale;state.phip=VectorXd::LinSpaced(4,.004,-.006)/scale;
    const auto x=qa.pack(state);const auto n=qa.electronDensity(x),h=qa.holeDensity(x);
    const auto y=pa.packPopulation(state.psi*scale,n,h);
    REQUIRE(pa.electronDensity(y)==n);REQUIRE(pa.holeDensity(y)==h);
    const auto restored=pa.unpack(y);REQUIRE(restored.n==n);REQUIRE(restored.p==h);
    const auto fields=pa.populationPhysicalFields(y);
    CHECK(fields.psi_V==restored.psi*scale);
    CHECK(fields.phin_V==restored.phin*scale);
    CHECK(fields.phip_V==restored.phip*scale);
    CHECK(fields.n==n);CHECK(fields.p==h);
    CoupledDDBoundaryConditions bc;bc.psi[0]=state.psi[0];bc.phin[0]=state.phin[0];bc.phip[0]=state.phip[0];
    const VectorXd a=qa.residual(x,bc),b=pa.residual(y,bc);
    for(int block=0;block<3;++block) {
        // A QF boundary identity has potential units, unlike a free PDE row.
        // Its inversion roundoff must not be normalized by a tiny PDE flux.
        const double denom=std::max(1e-30,a.segment(4*block+1,3).cwiseAbs().maxCoeff());
        CAPTURE(block,a.transpose(),b.transpose());
        CHECK((a-b).segment(4*block+1,3).cwiseAbs().maxCoeff()/denom<2e-10);
        const double bound=16*std::numeric_limits<Real>::epsilon()*std::max(.04,constants::Vt_300);
        CHECK(std::abs(a[4*block]-b[4*block])*scale<=bound);
    }
    const auto j=pa.assembleJacobian(y,bc);
    VectorXd v(12);v.head(4)=VectorXd::LinSpaced(4,-.003,.002)/scale;
    v.segment(4,4)=.03*n;v.tail(4)=-.02*h;
    const VectorXd tangent=j*v;
    for(double step:{1e-3,1e-4}) {
        const VectorXd fd=(pa.residual((y+step*v).eval(),bc)-pa.residual((y-step*v).eval(),bc))/(2*step);
        for(int block=0;block<3;++block) {
            const double denom=std::max(1e-30,tangent.segment(4*block+1,3).cwiseAbs().maxCoeff());
            CHECK((fd-tangent).segment(4*block+1,3).cwiseAbs().maxCoeff()/denom<3e-4);
            CHECK(std::abs(fd[4*block]-tangent[4*block])/std::max(1e-30,std::abs(tangent[4*block]))<3e-4);
        }
    }
    const auto terms=pa.carrierContinuityEquationTermDiagnostics(y,bc);
    for(int i=1;i<4;++i) {
        CHECK(terms[i].electronResidual==Catch::Approx(b[4+i]).margin(1e-30));
        CHECK(terms[i].holeResidual==Catch::Approx(b[8+i]).margin(1e-30));
    }
}

TEST_CASE("Population Poisson charge columns survive extremely dilute minority density", "[legacy_population]") {
    auto mesh=populationMesh();MaterialDatabase db;DopingModel doping(4);
    DDScalingSpec s;s.stateVariables="potential_density";
    CarrierStatisticsConfig stats;stats.model="fermi_dirac";
    CoupledDDAssembler a(mesh,db,doping,constants::Vt_300,mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),{},{},{},{},s,{},stats);
    const auto y=a.packPopulation(VectorXd::Zero(4),VectorXd::Constant(4,1e15),VectorXd::Constant(4,1e-300));
    const auto saved=a.unpack(y);REQUIRE(saved.p.minCoeff()==1e-300);
    REQUIRE(a.residual(y,{}).allFinite());
    const auto j=a.assembleJacobian(y,{});
    for(int i=0;i<4;++i) {REQUIRE(j.coeff(i,4+i)>0.);REQUIRE(j.coeff(i,8+i)<0.);}
    for(int c=0;c<j.outerSize();++c)for(SparseMatrixd::InnerIterator it(j,c);it;++it)REQUIRE(std::isfinite(it.value()));
}

TEST_CASE("Legacy population mode is explicit and retains incompatible-option guards", "[legacy_population]") {
    const auto c=newtonConfigFromJson({{"state_variables","potential_density"},{"coupled_operator","legacy"},
        {"carrier_statistics",{{"model","fermi_dirac"}}}});
    CHECK(c.stateVariables=="potential_density");CHECK(c.carrierStatistics.model=="fermi_dirac");
    for(const auto option:std::vector<nlohmann::json>{
        {{"split_dd_state",true}},{{"quasi_fermi_update_limit_V",.1}},
        {{"carrier_regularization_scale",1e-8}},{{"population_sg_current","band_population"}},
        {{"population_derivative_recovery","automatic"}}}) {
        nlohmann::json deck={{"state_variables","potential_density"},{"coupled_operator","legacy"},
            {"carrier_statistics",{{"model","fermi_dirac"}}}};
        deck.update(option);
        CHECK_THROWS(newtonConfigFromJson(deck));
    }
}

TEST_CASE("Production population trial uses opt-in compensation without changing its operator", "[legacy_population][writeback]") {
    auto mesh=populationMesh();MaterialDatabase db;DopingModel doping(4);
    auto cfg=newtonConfigFromJson({{"state_variables","potential_density"},{"coupled_operator","legacy"},
        {"carrier_statistics",{{"model","boltzmann"}}}});
    cfg.mobility=mobilityModelConfig("constant");
    NewtonSolver baseline(mesh,db,doping,{{"source",0.},{"drain",0.}},cfg);
    cfg.populationWriteback="boltzmann_psi_compensated";
    NewtonSolver corrected(mesh,db,doping,{{"source",0.},{"drain",0.}},cfg);
    DDSolution state;state.psi=VectorXd::Constant(4,.1);state.n=VectorXd::Constant(4,1e20);state.p=VectorXd::Constant(4,1e10);
    const auto y=baseline.packArclengthState(state);
    const auto a=baseline.makeArclengthSystem("drain"),b=corrected.makeArclengthSystem("drain");
    CHECK(a.residual(y,0.)==b.residual(y,0.));
    CHECK((a.jacobian(y,0.)-b.jacobian(y,0.)).norm()==0.);
    const VectorXd zero=VectorXd::Zero(12);
    const auto identity=b.makeTrial({y,0.},zero,0.,1.);
    CHECK(identity.state.x==y);CHECK(identity.state.lambda==0.);
    VectorXd raw=VectorXd::Zero(12);raw.head(4).setConstant(1e-17);
    raw.segment(4,4)=state.n*1e-16;raw.tail(4)=-state.p*1e-16;
    const auto old=a.makeTrial({y,0.},raw,1e-20,.5);
    const auto updated=b.makeTrial({y,0.},raw,1e-20,.5);
    std::vector<bool> active(4,true);
    PopulationTrial original;original.state=old.state.x;original.alpha=old.effectiveAlpha;
    const auto expected=compensateBoltzmannPsiWriteback(y,raw,original,1.,constants::Vt_300,active);
    CHECK(updated.state.x==expected.state);
    CHECK(updated.state.lambda==old.state.lambda);CHECK(updated.effectiveAlpha==old.effectiveAlpha);
}

TEST_CASE("Population-primary Fermi inverse controls relative rather than absolute density error", "[legacy_population]") {
    const double dos=2.8e25;
    for(double ratio:{1e-300,1e-20,5.94882e-13,1e-9,1e-7,.1,10.,100.}) {
        const double density=dos*ratio;
        const double eta=slot_population::populationEtaAccurate(density,dos);
        const double rebuilt=scaledFermiDiracHalf(dos,eta);
        CAPTURE(ratio,eta);
        CHECK(std::abs((static_cast<long double>(rebuilt)-density)/density)<3e-13L);
    }
}

TEST_CASE("Legacy population contact and load-line directions include one independent contact voltage", "[legacy_population]") {
    auto mesh=populationMesh();MaterialDatabase db;DopingModel doping(4);
    for(Index i=0;i<4;++i)doping.setNodeDoping(i,1e22,2e21);
    for(const auto model:{"boltzmann","fermi_dirac"}) {
        auto cfg=newtonConfigFromJson({{"state_variables","potential_density"},{"coupled_operator","legacy"},
            {"carrier_statistics",{{"model",model}}},{"reltol",0.},{"abstol",1e-10}});
        cfg.mobility=mobilityModelConfig("constant");
        NewtonSolver solver(mesh,db,doping,{{"source",0.},{"drain",.02}},cfg);
        const auto system=solver.makeArclengthSystem("drain");
        const auto port=solver.makeArclengthContactCurrentFunctional("drain",1e-6);
        CoupledDDAssembler q(mesh,db,doping,constants::Vt_300,cfg.mobility,
            recombinationModelConfig({"none"}),{},{},{},{},{},{},cfg.carrierStatistics);
        CoupledDDState s;s.psi=VectorXd::LinSpaced(4,.02,.04);
        s.phin=VectorXd::LinSpaced(4,.001,.008);s.phip=VectorXd::LinSpaced(4,.003,-.002);
        const auto x=q.pack(s);DDSolution state;state.psi=s.psi;state.phin=s.phin;state.phip=s.phip;
        state.n=q.electronDensity(x);state.p=q.holeDensity(x);
        const auto y=solver.packArclengthState(state);
        for(double tolerance:{0.,1e30}) {
            auto zeroCfg=cfg;zeroCfg.maxIter=0;zeroCfg.abstol=tolerance;
            NewtonSolver audit(mesh,db,doping,{{"source",0.},{"drain",.02}},zeroCfg);
            const auto result=audit.solve(state);
            CHECK(result.iters==0);CHECK(result.solution.psi==state.psi);
            CHECK(result.solution.n==state.n);CHECK(result.solution.p==state.p);
            CHECK(result.converged==(tolerance>0.));
        }
        VectorXd d(12);d.head(4)=VectorXd::LinSpaced(4,-.002,.003);
        d.segment(4,4)=.01*state.n;d.tail(4)=-.02*state.p;
        const Real dv=.003,bias=.02,resistance=1e12;
        const auto j=system.jacobian(y,bias);
        const auto fl=system.parameterDerivative(y,bias);
        const auto current=port.linearize(y,bias);
        const VectorXd action=j*d+fl*dv;
        REQUIRE(action.size()+1==13);
        const Real circuit=dv+resistance*(current.stateDerivative.dot(d)+current.parameterDerivative*dv);
        for(double h:{1e-3,1e-4}) {
            const auto yp=VectorXd(y+h*d),ym=VectorXd(y-h*d);
            const VectorXd fd=(system.residual(yp,bias+h*dv)-system.residual(ym,bias-h*dv))/(2*h);
            CHECK((fd-action).norm()/action.norm()<3e-4);
            const Real fdCircuit=dv+resistance*(port.value(yp,bias+h*dv)-port.value(ym,bias-h*dv))/(2*h);
            CHECK(std::abs(fdCircuit-circuit)/std::max(1e-30,std::abs(circuit))<3e-4);
        }
    }
}

TEST_CASE("Main population defaults remain scoped to BVMethods", "[legacy_population][config]") {
    const auto legacy=newtonConfigFromJson({{"coupled_operator","legacy"},
        {"state_variables","potential_density"},
        {"carrier_statistics",{{"model","boltzmann"}}}});
    CHECK(legacy.primaryVariableCoordinates=="physical");
    CHECK(legacy.populationStoppingRule=="legacy");
    CHECK(legacy.populationDerivativeRecovery=="full");
    CHECK(legacy.populationMinimumDensity_cm3==0.);
    CHECK(legacy.populationWriteback=="independent");
    const auto bv=newtonConfigFromJson({{"coupled_operator","bvmethods_tri3_qf"},
        {"carrier_statistics",{{"model","fermi_dirac_joyce_dixon"}}}});
    CHECK(bv.stateVariables=="potential_density");
    CHECK(bv.primaryVariableCoordinates=="unit");
    CHECK(bv.populationStoppingRule=="native_approx");
    CHECK(bv.populationDerivativeRecovery=="after_protection");
    CHECK(bv.populationMinimumDensity_cm3==1e-100);
}

TEST_CASE("Legacy Newton keeps its physical Jacobian for dilute populations",
          "[legacy_population][arclength]")
{
    auto mesh = populationMesh();
    MaterialDatabase database;
    DopingModel doping(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node)
        doping.setNodeDoping(node, 1e22, 0.);
    auto config = newtonConfigFromJson({
        {"coupled_operator", "legacy"},
        {"state_variables", "potential_density"},
        {"carrier_statistics", {{"model", "fermi_dirac"}, {"transport", "effective_band_sg"}}}
    });
    config.mobility = mobilityModelConfig("constant");
    config.maxIter = 1;
    NewtonSolver solver(mesh, database, doping, {{"source", 0.}, {"drain", .006}}, config);
    const auto system = solver.makeArclengthSystem("drain");

    DDSolution state;
    state.psi = VectorXd::Constant(4, .25);
    state.phin = VectorXd::Zero(4);
    state.phip = VectorXd::Zero(4);
    state.n = VectorXd::Constant(4, 1e22);
    state.p = VectorXd::Constant(4, 1e10);
    state.p[2] = 1e-160;
    const auto packed = solver.packArclengthState(state);
    const auto scales = system.columnScales({packed, .006});
    REQUIRE(scales[10] < std::sqrt(std::numeric_limits<Real>::min()));
    CHECK_FALSE(system.jacobianInColumns);
    CHECK_FALSE(system.bootstrapJacobianInColumns);
    CHECK(system.jacobian(packed, .006).coeff(10, 10) != 0.);
    CHECK_NOTHROW(solver.solve(state));
}
