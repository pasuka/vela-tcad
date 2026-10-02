from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=R/'tests/test_production_numerics.cpp'
s=(R/'tests/diagnostics/test_simplemos_linear_refinement.cpp').read_text().replace('"simplemos_linear_refinement.hpp"','"vela/solver/LinearRefinement.h"').replace('simplemos_diagnostic','vela::linear_refinement')
s+='\n'+(R/'tests/diagnostics/test_simplemos_stable_sg_derivative.cpp').read_text().replace('"simplemos_stable_sg_psi_derivative.hpp"','"vela/discretization/StableSGDerivative.h"').replace('simplemos_stable_sg','vela::stable_sg')
p.write_text(s,newline='\n')
p=R/'CMakeLists.txt';s=p.read_text().replace('add_executable(test_sg_flux tests/test_sg_flux.cpp)','add_executable(test_production_numerics tests/test_production_numerics.cpp)\ntarget_link_libraries(test_production_numerics PRIVATE vela_core Catch2::Catch2WithMain)\n\nadd_executable(test_sg_flux tests/test_sg_flux.cpp)').replace('vela_catch_discover_tests(test_sg_flux)','vela_catch_discover_tests(test_production_numerics)\nvela_catch_discover_tests(test_sg_flux)');p.write_text(s,newline='\n')
p=R/'tests/test_newton_solver.cpp'
s=p.read_text()+'\n'+r'''
TEST_CASE("Effective transport geometry reaches contact integration", "[newton][production_geometry]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = makePNDoping(mesh);
    const auto mobility = mobilityModelConfig("constant");
    const auto recombination = recombinationModelConfig({"none"});
    DDScalingSpec scaling;
    scaling.regionResolvedInterfaceAssembly.transportEdgeCouplingRatios.assign(mesh.edges().size(), 0.37);
    CoupledDDAssembler assembler(mesh, matdb, doping, constants::Vt_300,
        mobility, recombination, {}, {}, {}, {}, scaling);
    CoupledDDAssembler baseline(mesh, matdb, doping, constants::Vt_300, mobility, recombination);
    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state{VectorXd::LinSpaced(N,-.02,.03),
        VectorXd::LinSpaced(N,-.04,.01),VectorXd::LinSpaced(N,.02,-.03)};
    const auto x = assembler.pack(state);
    DDSolution solution;
    solution.psi=state.psi;solution.phin=state.phin;solution.phip=state.phip;
    solution.n=assembler.electronDensity(x);solution.p=assembler.holeDensity(x);
    ContactCurrent current(mesh,matdb,doping,mobility,constants::T0,scaling);
    ContactCurrent legacy(mesh,matdb,doping,mobility,constants::T0);
    const auto port=current.compute(solution,"anode");
    const auto residual=current.computeFromResidual(assembler,x,"anode");
    REQUIRE(port.totalCurrent == Catch::Approx(.37*legacy.compute(solution,"anode").totalCurrent).epsilon(1e-12));
    REQUIRE(port.totalCurrent == Catch::Approx(residual.totalCurrent).epsilon(1e-12));
    const auto rb=baseline.residual(x,{}), rc=assembler.residual(x,{});
    REQUIRE((rb.head(N)-rc.head(N)).norm()==0.);
    REQUIRE((rc.tail(2*N)-.37*rb.tail(2*N)).norm() < 1e-12*rb.tail(2*N).norm());
    scaling.regionResolvedInterfaceAssembly.transportEdgeCouplingRatios.pop_back();
    REQUIRE_THROWS_AS(ContactCurrent(mesh,matdb,doping,mobility,constants::T0,scaling),std::invalid_argument);
}

TEST_CASE("Charge-only signed volume leaves SRH and continuity Jacobian unchanged", "[newton][production_geometry]")
{
    const DeviceMesh mesh=makePNMesh();MaterialDatabase matdb;const auto doping=makePNDoping(mesh);
    const auto mobility=mobilityModelConfig("constant");const auto srh=recombinationModelConfig({"srh"});
    DDScalingSpec scaling;scaling.regionResolvedInterfaceAssembly.poissonChargeNodeVolume="signed_transport";
    CoupledDDAssembler baseline(mesh,matdb,doping,constants::Vt_300,mobility,srh);
    CoupledDDAssembler charge(mesh,matdb,doping,constants::Vt_300,mobility,srh,{}, {}, {}, {}, scaling);
    const int N=static_cast<int>(mesh.numNodes());
    CoupledDDState state{VectorXd::LinSpaced(N,-.02,.03),VectorXd::Constant(N,-.01),VectorXd::Constant(N,.01)};
    const auto x=baseline.pack(state);
    const auto r0=baseline.residual(x,{}),r1=charge.residual(x,{});
    REQUIRE((r0.tail(2*N)-r1.tail(2*N)).norm()==0.);
    const Eigen::MatrixXd j0=baseline.assembleJacobian(x,{}),j1=charge.assembleJacobian(x,{});
    REQUIRE((j0.bottomRows(2*N)-j1.bottomRows(2*N)).norm()==0.);
    REQUIRE((r0.head(N)-r1.head(N)).norm()>0.);
}

TEST_CASE("Production numerical policies reject ambiguous configuration", "[newton][production_geometry]")
{
    REQUIRE(newtonConfigFromJson(nlohmann::json::object()).linearRefinementIterations==0);
    REQUIRE(newtonConfigFromJson({{"linear_refinement_iterations",4}}).linearRefinementIterations==4);
    REQUIRE_THROWS_AS(newtonConfigFromJson({{"linear_refinement_iterations",-1}}),std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson({{"region_resolved_interface_assembly",{
        {"poisson_charge_node_volume","signed_transport"},{"poisson_electron_transport_node_volume",true}}}}),std::invalid_argument);
    REQUIRE_THROWS_AS(newtonConfigFromJson({{"region_resolved_interface_assembly",{
        {"transport_edge_coupling_ratios",{-1.}}}}}),std::invalid_argument);
}
'''
p.write_text(s,newline='\n')
