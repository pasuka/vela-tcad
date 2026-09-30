#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include "vela/core/UnitScaling.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/material/MaterialDatabase.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include "vela/physics/MobilityModel.h"
#include "vela/physics/RecombinationModel.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <array>

using namespace vela;

namespace {

DeviceMesh makePhuMobTriangle()
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0;    n0.y = 0.0;    mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0e-6; n1.y = 0.0;    mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.0;    n2.y = 1.0e-6; mesh.addNode(n2);
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

DeviceMesh makePhuMobMosTrianglePair()
{
    DeviceMesh mesh;
    Node n0; n0.id = 0; n0.x = 0.0;    n0.y = 0.0;     mesh.addNode(n0);
    Node n1; n1.id = 1; n1.x = 1.0e-6; n1.y = 0.0;     mesh.addNode(n1);
    Node n2; n2.id = 2; n2.x = 0.0;    n2.y = -1.0e-6; mesh.addNode(n2);
    Node n3; n3.id = 3; n3.x = 0.0;    n3.y = 1.0e-6;  mesh.addNode(n3);
    Cell siliconCell;
    siliconCell.id = 0;
    siliconCell.type = CellType::Tri3;
    siliconCell.region_id = 0;
    siliconCell.node_ids = {0, 2, 1};
    mesh.addCell(siliconCell);
    Cell oxideCell;
    oxideCell.id = 1;
    oxideCell.type = CellType::Tri3;
    oxideCell.region_id = 1;
    oxideCell.node_ids = {0, 1, 3};
    mesh.addCell(oxideCell);
    Region silicon;
    silicon.id = 0;
    silicon.name = "silicon";
    silicon.material = "Si";
    silicon.cell_ids = {0};
    mesh.addRegion(silicon);
    Region oxide;
    oxide.id = 1;
    oxide.name = "oxide";
    oxide.material = "SiO2";
    oxide.cell_ids = {1};
    mesh.addRegion(oxide);
    mesh.buildEdges();
    return mesh;
}

Index edgeByNodes(const DeviceMesh& mesh, Index a, Index b)
{
    for (Index edgeId = 0; edgeId < mesh.numEdges(); ++edgeId) {
        const Edge& edge = mesh.getEdge(edgeId);
        if ((edge.n0 == a && edge.n1 == b) ||
            (edge.n0 == b && edge.n1 == a)) {
            return edgeId;
        }
    }
    throw std::runtime_error("edge not found");
}

} // namespace

TEST_CASE("PhuMob chain rule resolves weak populations and the G floor",
          "[mobility][phumob][stable_chain]")
{
    // Independent 100-digit chain rule and 1e-20/1e-25 logarithmic-density
    // differences of the manual equations; references are in SI mobility units.
    struct Sample {
        PhuMobScalarState state;
        std::array<Real,4> expected;
        bool floor;
    };
    const std::array<Sample,3> samples{{
        {{1e23,2e21,8e22,3e20,300},
         {1.0385370458836108e-3,-5.5224972604801945e-5,
          -2.5310703066440914e-3,2.0004058349426506e-6},false},
        {{5e25,2e23,5e25,1e3,300},
         {4.8176497592391926e-3,7.043788402956757e-26,
          3.8496665629276428e-3,1.9462427440182106e-25},false},
        {{5e27,2e23,5e27,1e3,300},
         {9.507997427143333e-4,1.8252671486063443e-28,
          4.798362760818133e-3,1.230417748427984e-27},true}
    }};
    for (const auto& sample : samples) {
        for (int c=0;c<2;++c) {
            const auto carrier=c==0?CarrierType::Electron:CarrierType::Hole;
            const auto value=evaluatePhuMobLogDensityDerivatives(carrier,sample.state);
            CAPTURE(c,sample.state.donors);
            // Ratio checks have no absolute tolerance that could hide a zero.
            REQUIRE(std::abs(value.electrons/sample.expected[2*c]-1)<2e-12);
            REQUIRE(std::abs(value.holes/sample.expected[2*c+1]-1)<2e-12);
            if(sample.floor)
                REQUIRE(evaluatePhuMobScalar(carrier,sample.state).screeningGDerivative==0.0);
        }
    }
}

TEST_CASE("PhuMob log-density derivatives preserve units and temperature dependence",
          "[mobility][phumob][stable_chain]")
{
    for(Real temperature : {100.,300.,350.}) {
        PhuMobScalarState state{1e23,2e22,8e22,3e21,temperature};
        PhuMobParameters params;
        for(auto species : {PhuMobDonorSpecies::Arsenic,PhuMobDonorSpecies::Phosphorus}) {
            params.donorSpecies=species;
            for(auto carrier : {CarrierType::Electron,CarrierType::Hole}) {
                const auto d=evaluatePhuMobLogDensityDerivatives(carrier,state,params);
                for(int population=0;population<2;++population) {
                    auto plus=state,minus=state;
                    const Real step=1e-4;
                    if(population==0) { plus.electrons*=std::exp(step);minus.electrons*=std::exp(-step); }
                    else { plus.holes*=std::exp(step);minus.holes*=std::exp(-step); }
                    const Real fd=(evaluatePhuMobScalar(carrier,plus,params).mobility-
                                   evaluatePhuMobScalar(carrier,minus,params).mobility)/(2*step);
                    const Real derivative=population==0?d.electrons:d.holes;
                    CAPTURE(temperature,population);
                    REQUIRE(std::abs(fd/derivative-1)<2e-7);
                }
                auto converted=params;
                converted.internalConcentrationToCm3=1.;
                converted.internalMobilityToCm2PerVS=1.;
                for(auto cp : {&converted.electronArsenic,&converted.electronPhosphorus,&converted.holeBoron}) {
                    cp->muMax*=1e4;cp->muMin*=1e4;cp->nRef*=1e-6;
                }
                converted.donorClusterReference*=1e-6;converted.acceptorClusterReference*=1e-6;
                auto local=state;local.donors*=1e-6;local.acceptors*=1e-6;
                local.electrons*=1e-6;local.holes*=1e-6;
                const auto dc=evaluatePhuMobLogDensityDerivatives(carrier,local,converted);
                REQUIRE(std::abs(dc.electrons*1e-4/d.electrons-1)<2e-13);
                REQUIRE(std::abs(dc.holes*1e-4/d.holes-1)<2e-13);
            }
        }
    }
}

TEST_CASE("PhuMob assembly preserves a weak cross column and conservative signs",
          "[mobility][phumob][stable_chain][jacobian]")
{
    const auto mesh=makePhuMobTriangle();
    const MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    for(Index node=0;node<mesh.numNodes();++node)
        doping.setNodeDoping(node,5e27,2e23);
    CoupledDDAssembler assembler(mesh,materials,doping,constants::Vt_300,
        mobilityModelConfig("phumob"),recombinationModelConfig({"none"}));
    const Real ni=materials.getMaterial("Si",300.).ni;
    CoupledDDState state;
    state.psi=VectorXd(3);state.psi<<0.,.01,-.02;
    state.phin=state.psi.array()-constants::Vt_300*std::log(5e27/ni);
    state.phip=state.psi.array()+constants::Vt_300*std::log(1e3/ni);
    const auto x=assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const auto J=assembler.assembleJacobian(x,bcs);
    const auto edges=assembler.sgEdgeFluxDiagnostics(x,bcs);
    const auto edge=std::find_if(edges.begin(),edges.end(),[](const auto& e) {
        return (e.node0==0 && e.node1==1) || (e.node0==1 && e.node1==0);
    });
    REQUIRE(edge!=edges.end());
    // Independent 100-digit weak derivative at this uniform carrier state.
    const Real dmuDlogP=1.8252671486063443e-28;
    const Real expected=(edge->node0==0?1.:-1.)*edge->electronFlux/
        edge->electronMobility_m2_V_s*dmuDlogP/(2*constants::Vt_300);
    REQUIRE(expected!=0.);
    REQUIRE(std::abs(J.coeff(3,7)/expected-1)<1e-10);
    // Every unconstrained continuity edge contributes equal and opposite rows.
    for(int column=0;column<9;++column) {
        Real sum=0.,scale=0.;
        for(int row=3;row<6;++row) { sum+=J.coeff(row,column);scale+=std::abs(J.coeff(row,column)); }
        REQUIRE(std::abs(sum)<=4e-15*scale);
    }
}

TEST_CASE("PhuMob defaults reproduce T-2022.03 silicon parameters",
          "[mobility][phumob][parameters]")
{
    const MobilityModelConfig config = mobilityModelConfig("constant");
    const PhuMobParameters& p = config.phuMob;

    REQUIRE(p.donorSpecies == PhuMobDonorSpecies::Arsenic);
    REQUIRE(p.electronArsenic.muMax == Catch::Approx(1417.0e-4));
    REQUIRE(p.electronArsenic.muMin == Catch::Approx(52.2e-4));
    REQUIRE(p.electronArsenic.theta == Catch::Approx(2.285));
    REQUIRE(p.electronArsenic.nRef == Catch::Approx(9.68e22));
    REQUIRE(p.electronArsenic.alpha == Catch::Approx(0.68));
    REQUIRE(p.electronPhosphorus.muMax == Catch::Approx(1414.0e-4));
    REQUIRE(p.electronPhosphorus.muMin == Catch::Approx(68.5e-4));
    REQUIRE(p.electronPhosphorus.nRef == Catch::Approx(9.20e22));
    REQUIRE(p.electronPhosphorus.alpha == Catch::Approx(0.711));
    REQUIRE(p.holeBoron.muMax == Catch::Approx(470.5e-4));
    REQUIRE(p.holeBoron.muMin == Catch::Approx(44.9e-4));
    REQUIRE(p.holeBoron.theta == Catch::Approx(2.247));
    REQUIRE(p.holeBoron.nRef == Catch::Approx(2.23e23));
    REQUIRE(p.holeBoron.alpha == Catch::Approx(0.719));
    REQUIRE(p.donorClusterReference == Catch::Approx(4.0e26));
    REQUIRE(p.acceptorClusterReference == Catch::Approx(7.2e26));
    REQUIRE(p.gA == Catch::Approx(0.89233));
    REQUIRE(p.gAlphaPrime == Catch::Approx(0.72169));
    REQUIRE(p.gGamma == Catch::Approx(1.80618));
}

TEST_CASE("PhuMob scalar kernel matches T-2022.03 equations at 300 K",
          "[mobility][phumob][formula]")
{
    const PhuMobScalarState state{
        1.0e23, // donors: 1e17 cm^-3
        2.0e21, // acceptors: 2e15 cm^-3
        8.0e22, // electrons: 8e16 cm^-3
        3.0e20, // holes: 3e14 cm^-3
        300.0,
    };

    const PhuMobScalarResult electron = evaluatePhuMobScalar(
        CarrierType::Electron, state);
    const PhuMobScalarResult hole = evaluatePhuMobScalar(
        CarrierType::Hole, state);

    REQUIRE(electron.latticeMobility == Catch::Approx(0.1417).epsilon(1.0e-13));
    REQUIRE(electron.screeningParameter ==
            Catch::Approx(63.25703274767962).epsilon(1.0e-12));
    REQUIRE(electron.screeningF ==
            Catch::Approx(1.1441331997437554).epsilon(1.0e-12));
    REQUIRE(electron.screeningG ==
            Catch::Approx(0.6075913697623474).epsilon(1.0e-12));
    REQUIRE(electron.scatteringMobility ==
            Catch::Approx(0.14713091713819333).epsilon(1.0e-12));
    REQUIRE(electron.mobility ==
            Catch::Approx(0.07218219976259292).epsilon(1.0e-12));

    REQUIRE(hole.latticeMobility == Catch::Approx(0.04705).epsilon(1.0e-13));
    REQUIRE(hole.screeningParameter ==
            Catch::Approx(46.10665223671113).epsilon(1.0e-12));
    REQUIRE(hole.screeningF ==
            Catch::Approx(1.4851524381449164).epsilon(1.0e-12));
    REQUIRE(hole.screeningG ==
            Catch::Approx(0.5771285936866543).epsilon(1.0e-12));
    REQUIRE(hole.scatteringMobility ==
            Catch::Approx(0.0999653272864784).epsilon(1.0e-12));
    REQUIRE(hole.mobility ==
            Catch::Approx(0.0319923693375432).epsilon(1.0e-12));
}

TEST_CASE("PhuMob uses arsenic by default and supports phosphorus selection",
          "[mobility][phumob][species]")
{
    const PhuMobScalarState state{1.0e24, 0.0, 1.0e24, 1.0e10, 300.0};
    PhuMobParameters arsenic;
    PhuMobParameters phosphorus = arsenic;
    phosphorus.donorSpecies = PhuMobDonorSpecies::Phosphorus;

    const Real muAs = evaluatePhuMobScalar(
        CarrierType::Electron, state, arsenic).mobility;
    const Real muP = evaluatePhuMobScalar(
        CarrierType::Electron, state, phosphorus).mobility;

    REQUIRE(muAs > 0.0);
    REQUIRE(muP > 0.0);
    REQUIRE(muAs != Catch::Approx(muP));
}

TEST_CASE("PhuMob physical result is invariant across Vela unit systems",
          "[mobility][phumob][units]")
{
    const MobilityModelConfig siConfig = mobilityModelConfigFromJson(
        nlohmann::json{{"model", "constant"}});
    const MobilityModelConfig tcadConfig = mobilityModelConfigFromJson(
        nlohmann::json{{"model", "constant"}},
        UnitScalingConfig{UnitScalingMode::UnitScaling});
    const PhysicalUnitSystem tcad = PhysicalUnitSystem::tcadInternal();

    const PhuMobScalarState siState{
        2.5e23, 8.0e22, 1.2e23, 4.0e21, 350.0};
    const PhuMobScalarState tcadState{
        tcad.m3ToInternalConcentration(siState.donors),
        tcad.m3ToInternalConcentration(siState.acceptors),
        tcad.m3ToInternalConcentration(siState.electrons),
        tcad.m3ToInternalConcentration(siState.holes),
        siState.temperature_K,
    };

    for (const CarrierType carrier : {CarrierType::Electron, CarrierType::Hole}) {
        const PhuMobScalarResult si = evaluatePhuMobScalar(
            carrier, siState, siConfig.phuMob);
        const PhuMobScalarResult internal = evaluatePhuMobScalar(
            carrier, tcadState, tcadConfig.phuMob);
        REQUIRE(tcad.internalMobilityToM2PerVS(internal.mobility) ==
                Catch::Approx(si.mobility).epsilon(1.0e-12));
        REQUIRE(tcad.internalMobilityToM2PerVS(internal.latticeMobility) ==
                Catch::Approx(si.latticeMobility).epsilon(1.0e-12));
        REQUIRE(tcad.internalMobilityToM2PerVS(internal.scatteringMobility) ==
                Catch::Approx(si.scatteringMobility).epsilon(1.0e-12));
        REQUIRE(tcad.internalConcentrationToM3(
                    internal.effectiveScatteringConcentration) ==
                Catch::Approx(si.effectiveScatteringConcentration).epsilon(1.0e-12));
    }
}

TEST_CASE("PhuMob JSON parameters follow the active unit system",
          "[mobility][phumob][json][units]")
{
    const nlohmann::json json = {
        {"model", "constant"},
        {"phumob", {
            {"donor_species", "phosphorus"},
            {"electron_phosphorus", {
                {"mu_max_m2_V_s", 1500.0},
                {"mu_min_m2_V_s", 75.0},
                {"n_ref_m3", 1.1e17},
                {"theta", 2.1},
                {"alpha", 0.75},
            }},
            {"donor_cluster_reference_m3", 5.0e20},
            {"g_gamma", 1.9},
        }},
    };
    const MobilityModelConfig config = mobilityModelConfigFromJson(
        json, UnitScalingConfig{UnitScalingMode::UnitScaling});

    REQUIRE(config.phuMob.donorSpecies == PhuMobDonorSpecies::Phosphorus);
    REQUIRE(config.phuMob.electronPhosphorus.muMax == Catch::Approx(1500.0));
    REQUIRE(config.phuMob.electronPhosphorus.muMin == Catch::Approx(75.0));
    REQUIRE(config.phuMob.electronPhosphorus.nRef == Catch::Approx(1.1e17));
    REQUIRE(config.phuMob.electronPhosphorus.theta == Catch::Approx(2.1));
    REQUIRE(config.phuMob.electronPhosphorus.alpha == Catch::Approx(0.75));
    REQUIRE(config.phuMob.donorClusterReference == Catch::Approx(5.0e20));
    REQUIRE(config.phuMob.gGamma == Catch::Approx(1.9));
    REQUIRE(config.phuMob.internalConcentrationToCm3 == Catch::Approx(1.0));
    REQUIRE(config.phuMob.internalMobilityToCm2PerVS == Catch::Approx(1.0));
}

TEST_CASE("PhuMob format-version 2 accepts only TCAD unit key spellings",
          "[mobility][phumob][json][units][format-v2]")
{
    const nlohmann::json deck = {
        {"format_version", 2},
        {"solver", {
            {"mobility", {
                {"model", "constant"},
                {"phumob", {
                    {"electron_arsenic", {
                        {"mu_max_cm2_V_s", 1420.0},
                        {"mu_min_cm2_V_s", 53.0},
                        {"n_ref_cm3", 9.9e16},
                    }},
                    {"acceptor_cluster_reference_cm3", 7.3e20},
                }},
            }},
        }},
    };
    const nlohmann::json canonical = canonicalizeDeck(deck);
    const UnitScalingConfig scaling = parseUnitScalingConfig(deck);
    const MobilityModelConfig config = mobilityModelConfigFromJson(
        canonical.at("solver").at("mobility"), scaling);

    REQUIRE(config.phuMob.electronArsenic.muMax == Catch::Approx(1420.0));
    REQUIRE(config.phuMob.electronArsenic.muMin == Catch::Approx(53.0));
    REQUIRE(config.phuMob.electronArsenic.nRef == Catch::Approx(9.9e16));
    REQUIRE(config.phuMob.acceptorClusterReference == Catch::Approx(7.3e20));

    nlohmann::json legacyKeyDeck = deck;
    auto& arsenic = legacyKeyDeck["solver"]["mobility"]["phumob"]
        ["electron_arsenic"];
    arsenic["mu_max_m2_V_s"] = arsenic["mu_max_cm2_V_s"];
    arsenic.erase("mu_max_cm2_V_s");
    REQUIRE_THROWS_AS(canonicalizeDeck(legacyKeyDeck), std::invalid_argument);
}

TEST_CASE("PhuMob zero-scattering limit is the lattice mobility",
          "[mobility][phumob][limit]")
{
    const PhuMobScalarResult result = evaluatePhuMobScalar(
        CarrierType::Electron, PhuMobScalarState{});
    REQUIRE(result.mobility == Catch::Approx(0.1417));
    REQUIRE(result.latticeMobility == Catch::Approx(0.1417));
    REQUIRE(std::isinf(result.scatteringMobility));
    REQUIRE(std::isinf(result.screeningParameter));
}

TEST_CASE("PhuMob rejects nonphysical scalar states and parameters",
          "[mobility][phumob][validation]")
{
    PhuMobScalarState badState;
    badState.donors = -1.0;
    REQUIRE_THROWS_AS(
        evaluatePhuMobScalar(CarrierType::Electron, badState),
        std::invalid_argument);

    PhuMobParameters badParameters;
    badParameters.gGamma = 0.0;
    REQUIRE_THROWS_AS(
        evaluatePhuMobScalar(
            CarrierType::Electron, PhuMobScalarState{}, badParameters),
        std::invalid_argument);
}

TEST_CASE("PhuMob mobility model retains compensated impurities and carrier state",
          "[mobility][phumob][m6][dd]")
{
    const Material silicon = MaterialDatabase{}.getMaterial("Si", 300.0);
    const MobilityModelConfig config = mobilityModelConfig("phumob");
    const auto mobility = makeMobilityModel(config);
    const PhuMobScalarState state{
        1.0e23, 2.0e21, 8.0e22, 3.0e20, 300.0};

    const Real electron = mobility->electronMobilityWithIonizedImpurities(
        silicon, state.donors, state.acceptors, state.electrons, state.holes);
    const Real hole = mobility->holeMobilityWithIonizedImpurities(
        silicon, state.donors, state.acceptors, state.electrons, state.holes);
    REQUIRE(electron == Catch::Approx(
        evaluatePhuMobScalar(CarrierType::Electron, state).mobility)
        .epsilon(1.0e-13));
    REQUIRE(hole == Catch::Approx(
        evaluatePhuMobScalar(CarrierType::Hole, state).mobility)
        .epsilon(1.0e-13));

    const Real reconstructedNetOnly = mobility->electronMobility(
        silicon, state.donors - state.acceptors,
        state.electrons, state.holes);
    REQUIRE(reconstructedNetOnly != Catch::Approx(electron));
}

TEST_CASE("PhuMob field model composes low-field state with saturation",
          "[mobility][phumob][m6][field]")
{
    const Material silicon = MaterialDatabase{}.getMaterial("Si", 300.0);
    const auto lowField = makeMobilityModel(mobilityModelConfig("phumob"));
    const auto highField = makeMobilityModel(mobilityModelConfig("phumob_field"));
    const Real donors = 2.0e23;
    const Real acceptors = 1.0e22;
    const Real electrons = 1.2e23;
    const Real holes = 4.0e21;
    const Real baseline = lowField->electronMobilityWithIonizedImpurities(
        silicon, donors, acceptors, electrons, holes);
    const Real limited = highField->electronMobilityWithIonizedImpurities(
        silicon, donors, acceptors, electrons, holes, 2.0e7);
    REQUIRE(baseline > 0.0);
    REQUIRE(limited > 0.0);
    REQUIRE(limited < baseline);
}

TEST_CASE("PhuMob coupled DD Jacobian includes self-consistent carrier scattering",
          "[mobility][phumob][m6][jacobian]")
{
    const DeviceMesh mesh = makePhuMobTriangle();
    const MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(0, 1.0e23, 2.0e22);
    doping.setNodeDoping(1, 1.2e23, 1.0e22);
    doping.setNodeDoping(2, 8.0e22, 3.0e22);
    MobilityModelConfig mobility = mobilityModelConfig("phumob");
    const RecombinationModelConfig recombination =
        recombinationModelConfig({"none"});
    CoupledDDAssembler assembler(
        mesh, materials, doping, constants::Vt_300,
        mobility, recombination);

    CoupledDDState state;
    state.psi = VectorXd(3);
    state.phin = VectorXd(3);
    state.phip = VectorXd(3);
    state.psi << 0.02, 0.00, -0.01;
    state.phin << -0.35, -0.37, -0.34;
    state.phip << 0.35, 0.33, 0.36;
    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions boundaries;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, boundaries);
    const SparseMatrixd finiteDifference =
        assembler.finiteDifferenceJacobian(x, boundaries, 1.0e-7);
    const Real denominator = std::max<Real>(finiteDifference.norm(), 1.0e-30);
    const Real relativeError =
        (analytic - finiteDifference).norm() / denominator;
    CAPTURE(relativeError);
    REQUIRE(relativeError < 2.0e-4);

    const int nodeCount = static_cast<int>(mesh.numNodes());
    Real electronToHoleCrossNorm = 0.0;
    Real holeToElectronCrossNorm = 0.0;
    for (int row = 0; row < nodeCount; ++row) {
        for (int column = 0; column < nodeCount; ++column) {
            electronToHoleCrossNorm += std::abs(
                analytic.coeff(nodeCount + row, 2 * nodeCount + column));
            holeToElectronCrossNorm += std::abs(
                analytic.coeff(2 * nodeCount + row, nodeCount + column));
        }
    }
    REQUIRE(electronToHoleCrossNorm > 0.0);
    REQUIRE(holeToElectronCrossNorm > 0.0);
}

TEST_CASE("PhuMob Enormal defaults reproduce T-2022.03 Silicon parameters",
          "[mobility][phumob][m7][parameters]")
{
    const MobilityModelConfig config =
        mobilityModelConfig("phumob_field_lombardi");
    REQUIRE(isPhuMobModel(config));
    REQUIRE(isSurfaceMobilityModel(config));
    REQUIRE(config.electronLombardi.B == Catch::Approx(4.7500e5));
    REQUIRE(config.holeLombardi.B == Catch::Approx(9.9250e4));
    REQUIRE(config.electronLombardi.C == Catch::Approx(
        5.8000e2 * 4.641588833612778e-4));
    REQUIRE(config.holeLombardi.C == Catch::Approx(
        2.9470e3 * 4.641588833612778e-4));
    REQUIRE(config.electronLombardi.delta == Catch::Approx(5.8200e10));
    REQUIRE(config.holeLombardi.delta == Catch::Approx(2.0546e10));
    REQUIRE(config.electronLombardi.eta == Catch::Approx(5.8200e32));
    REQUIRE(config.holeLombardi.eta == Catch::Approx(2.0546e32));
    REQUIRE(config.electronLombardi.criticalLength == Catch::Approx(1.0e-8));
    REQUIRE(config.holeLombardi.criticalLength == Catch::Approx(1.0e-8));
}

TEST_CASE("PhuMob Enormal is low-field input to high-field saturation",
          "[mobility][phumob][m7][formula]")
{
    const Material silicon = MaterialDatabase{}.getMaterial("Si", 350.0);
    const MobilityModelConfig config =
        mobilityModelConfig("phumob_field_lombardi");
    const auto combined = makeMobilityModel(config);
    const PhuMobScalarState state{
        2.0e23, 5.0e22, 1.2e23, 4.0e21, 350.0};
    const Real parallelField = 2.0e7;
    const Real normalField = 1.0e7;
    const Real bulk = evaluatePhuMobScalar(
        CarrierType::Electron, state, config.phuMob).mobility;
    const auto& lp = config.electronLombardi;
    const Real totalImpurity = state.donors + state.acceptors;
    const Real muAc = lp.B / normalField +
        lp.C * std::pow(state.temperature_K / 300.0, -lp.k) *
        std::pow((totalImpurity + lp.N2) / lp.N0, lp.lambda) /
        std::cbrt(normalField);
    const Real exponent = lp.A + lp.alpha * state.electrons /
        std::pow(totalImpurity + lp.N1, lp.nu);
    const Real inverseSurface =
        std::pow(normalField / 100.0, exponent) / lp.delta +
        std::pow(normalField, 3.0) / lp.eta;
    const Real surfaceLowField = 1.0 /
        (1.0 / bulk + 1.0 / muAc + inverseSurface);
    const Real ratio = surfaceLowField * parallelField /
        config.electronField.saturationVelocity;
    const Real expected = surfaceLowField /
        std::pow(1.0 + std::pow(ratio, config.electronField.beta),
                 1.0 / config.electronField.beta);
    const Real actual = combined->electronMobilityWithIonizedImpurities(
        silicon, state.donors, state.acceptors, state.electrons, state.holes,
        parallelField, normalField, 0.0);

    REQUIRE(actual == Catch::Approx(expected).epsilon(1.0e-12));
    REQUIRE(actual < surfaceLowField);
    REQUIRE(surfaceLowField < bulk);
}

TEST_CASE("PhuMob Enormal re-evaluates the live spatial normal field",
          "[mobility][phumob][m7][spatial]")
{
    const DeviceMesh mesh = makePhuMobMosTrianglePair();
    const MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node)
        doping.setNodeDoping(node, 1.0e23, 2.0e22);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const auto cellMaterials = detail::buildCellMaterials(
        mesh, materials, 300.0);
    MobilityModelConfig config = mobilityModelConfig("phumob_lombardi");
    config.surface.surfaceInterface = {"silicon", "oxide"};
    VectorXd psi(4);
    psi << 0.0, 0.0, -0.5, 0.5;
    detail::updateSurfaceMobilityCellGeometry(
        config, mesh, edgeCells, psi, 1.0, &cellMaterials);
    const auto mobility = makeMobilityModel(config);
    const detail::EdgeMobilityCarrierState carriers{
        8.0e22, 8.0e22, 3.0e20, 3.0e20};
    const Index interfaceEdge = edgeByNodes(mesh, 0, 1);
    const Real first = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config, &psi, &carriers);
    const Real cachedField = config.surface.cellNormalFields.at(0);
    const Real firstLiveField = detail::liveSurfaceNormalFieldForCell(
        config, mesh, psi, 0);
    const Real firstExpected = mobility->electronMobilityWithIonizedImpurities(
        cellMaterials.at(0), 1.0e23, 2.0e22, 8.0e22, 3.0e20,
        0.0, 5.0e5, config.surface.cellDistances.at(0));
    psi(2) = -1.0;
    const Real second = detail::edgeMobility(
        edgeCells, mesh, doping, *mobility, cellMaterials, interfaceEdge,
        CarrierType::Electron, 0.0, &config, &psi, &carriers);
    const Real secondLiveField = detail::liveSurfaceNormalFieldForCell(
        config, mesh, psi, 0);
    const Real secondExpected = mobility->electronMobilityWithIonizedImpurities(
        cellMaterials.at(0), 1.0e23, 2.0e22, 8.0e22, 3.0e20,
        0.0, 1.0e6, config.surface.cellDistances.at(0));

    REQUIRE(cachedField == Catch::Approx(5.0e5));
    REQUIRE(firstLiveField == Catch::Approx(5.0e5));
    REQUIRE(config.surface.cellNormalFields.at(0) == cachedField);
    REQUIRE(secondLiveField == Catch::Approx(1.0e6));
    REQUIRE(first == Catch::Approx(firstExpected).epsilon(1.0e-12));
    REQUIRE(second == Catch::Approx(secondExpected).epsilon(1.0e-12));
    REQUIRE(second > 0.0);
    REQUIRE(second < first);
}

TEST_CASE("PhuMob Enormal coupled DD Jacobian includes spatial field closure",
          "[mobility][phumob][m7][jacobian]")
{
    const DeviceMesh mesh = makePhuMobMosTrianglePair();
    const MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    for (Index node = 0; node < mesh.numNodes(); ++node)
        doping.setNodeDoping(node, 1.0e23, 2.0e22);
    MobilityModelConfig mobility =
        mobilityModelConfig("phumob_lombardi");
    mobility.surface.surfaceInterface = {"silicon", "oxide"};
    const RecombinationModelConfig recombination =
        recombinationModelConfig({"none"});
    CoupledDDAssembler assembler(
        mesh, materials, doping, constants::Vt_300,
        mobility, recombination);

    CoupledDDState state;
    state.psi = VectorXd(4);
    state.phin = VectorXd(4);
    state.phip = VectorXd(4);
    state.psi << 0.02, 0.00, -0.25, 0.30;
    state.phin << -0.35, -0.37, -0.34, 0.0;
    state.phip << 0.35, 0.33, 0.36, 0.0;
    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions boundaries;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, boundaries);
    const SparseMatrixd finiteDifference =
        assembler.finiteDifferenceJacobian(x, boundaries, 1.0e-7);
    const Real denominator = std::max<Real>(finiteDifference.norm(), 1.0e-30);
    const Real relativeError =
        (analytic - finiteDifference).norm() / denominator;
    CAPTURE(relativeError);
    REQUIRE(relativeError < 4.0e-4);
}
