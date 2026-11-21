#include "vela/boundary/BoundaryCondition.h"
#include "vela/core/ArithmeticDecimal.h"
#include "vela/core/PhysicalConstants.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/discretization/ScharfetterGummel.h"
#include "vela/equation/Tri3LocalForwardAD.h"
#include "vela/equation/PoissonAssembler.h"
#include "vela/equation/IalElementAvalanche.h"
#include "vela/equation/IalPopulationContinuity.h"
#include "vela/equation/ExtendedBoltzmannMaterials.h"
#include "vela/equation/ExtendedAveragePoissonGeometry.h"
#include "vela/simulation/ConfigParsing.h"
#include "vela/equation/ExtendedPopulationDevice.h"
#include "vela/solver/ExtendedSparseLu.h"
#include "vela/solver/ExtendedBoltzmannColdStart.h"
#include "vela/solver/ExtendedFixedBiasNewton.h"
#include "vela/solver/DensityUpdateConvergence.h"
#include "vela/solver/LinearRefinement.h"
#include "vela/solver/PopulationState.h"
#include "vela/io/MeshReader.h"
#include "vela/io/StateIdentity.h"
#include "vela/physics/IalMobilityJson.h"
#include "vela/physics/IalInterfaceGeometry.h"
#include "vela/physics/RecombinationModel.h"
#include "vela/simulation/DCSweep.h"
#include "vela/simulation/DCSweepStepControl.h"
#include "vela/simulation/LoadLineResidual.h"
#include "vela/simulation/BorderedLoadLineSystem.h"
#include "vela/simulation/ExtendedBorderedLoadLine.h"
#include "vela/simulation/CoupledLoadLine.h"
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>

using namespace vela;
using nlohmann::json;

namespace {
constexpr auto low = 0x1p-60L;
constexpr auto nextToOneDecimal = "1.000000000000000000867361737988403547205962240695953369140625";
static_assert(compiledArithmeticPrecision == ArithmeticPrecision::Extended80);
}

TEST_CASE("Extended population extrapolation preserves positive densities", "[extended80]") {
    const auto old = (ExtendedVector(6) << 0.25L, 1.0L, 3.0L, 0.0L, 2.0L, 0.0L).finished();
    const auto current = (ExtendedVector(6) << 0.5L, 1.0L + low, 1.0L, 0.0L, 4.0L, 0.0L).finished();
    const auto predicted = predictPopulation<long double>(old, current, 2.0L);
    CHECK(predicted[0] == 1.0L);
    CHECK(predicted[1] == 1.0L + 3.0L * low);
    CHECK(std::abs(9.0L * predicted[2] - 1.0L) < 16 * std::numeric_limits<long double>::epsilon());
    CHECK(predicted[3] == 0.0L);
    CHECK(std::abs(predicted[4] / 16.0L - 1.0L) < 16 * std::numeric_limits<long double>::epsilon());
    CHECK(predicted[5] == 0.0L);

    auto invalid = old;
    invalid[2] = 0.0L;
    CHECK_THROWS_AS(predictPopulation<long double>(invalid, current, 2.0L), std::domain_error);
}

TEST_CASE("Prepared extended refinement defects preserve uncached arithmetic", "[extended80]") {
    ExtendedSparseMatrix matrix(3, 3);
    matrix.insert(0, 0) = 1.0L + low;
    matrix.insert(2, 0) = 0.0L;
    matrix.insert(0, 1) = -1.0L;
    matrix.insert(2, 1) = 0x1p+600L;
    matrix.insert(2, 2) = -0x1p+600L;
    const ExtendedVector residual = (ExtendedVector(3) << -low, 0.0L, 0x1p+540L).finished();
    ExtendedVector direction = ExtendedVector::Ones(3);
    const linear_refinement::PreparedMpDefect<long double> prepared(matrix);
    for (const long double shift : {0.0L, low, -low}) {
        direction(2) = 1.0L + shift;
        const auto expected = linear_refinement::computeHighPrecisionDefect(matrix, residual, direction);
        const auto actual = prepared.evaluate(residual, direction);
        CHECK((actual.values.array() == expected.values.array()).all());
        CHECK((actual.componentwiseBackwardErrors.array() == expected.componentwiseBackwardErrors.array()).all());
        CHECK(actual.backwardMax == expected.backwardMax);
    }
}

TEST_CASE("Load-line trial comparisons retain residual and decrease low bits", "[extended80]") {
    CoupledLoadLineConfig config;
    config.equationTolerance = 1;
    config.loadLineTolerance_V = 1;
    config.filterGamma = static_cast<Real>(low);
    CHECK(coupledNormalizedMerit(config, 1.L + low, .5L) == 1.L + low);
    CHECK(coupledNormalizedMerit(config, static_cast<Real>(1.L + low), .5) == 1.);
    CHECK_FALSE(coupledResidualFilterAccepts<long double>(
        config, 1.L + low, 2.L, 1.L + low, 2.L, 1.L));
    CHECK(coupledResidualFilterAccepts<long double>(
        config, 1.L + low, 2.L, 1.L, 2.L, 1.L));
}

TEST_CASE("SLOT contact and sweep consumers keep authoritative voltage input", "[extended80]") {
    const json deck = {
        {"contacts", json::array({{{"name", "drain"}, {"bias", 1.0}, {"flatband_voltage", 1.0}}})},
        {"sweep", {{"contact", "drain"}, {"start", 0.0}, {"stop", 2.0}, {"step", 1.0},
                   {"initial_step", 1.0}, {"min_step", 0.0625}, {"max_step", 1.0},
                   {"bias_points", {1.0, 2.0}}, {"vtk_bias_points", {1.0}}}},
        {"decimal_stage_controls",
         {{"contacts", {{"drain", nextToOneDecimal}}},
          {"sweep", {{"start", "0"}, {"stop", "2"}, {"step", nextToOneDecimal},
                     {"shrink_factor", "0.5"}, {"growth_factor", "1"},
                     {"initial_step", nextToOneDecimal}, {"min_step", "0.0625"},
                     {"max_step", nextToOneDecimal}, {"bias_points", {nextToOneDecimal, "2"}},
                     {"vtk_bias_points", {nextToOneDecimal}}}}}}
    };
    const auto contacts = parseContactBoundarySpecs(deck);
    REQUIRE(contacts.size() == 1);
    REQUIRE(contacts.front().decimalInitialBias_V.has_value());
    CHECK(effectivePoissonDirichletPotential(contacts.front()) == 0);
    CHECK(arithmeticPoissonDirichletPotential(contacts.front(), *contacts.front().decimalInitialBias_V) == low);
    CHECK(arithmeticPoissonDirichletPotential(contacts.front(), 1 + 2 * low) == 2 * low);
    const auto sweep = readArithmeticSweepVoltageControls(deck);
    REQUIRE(sweep.has_value());
    const auto& exact = *sweep;
    CHECK(static_cast<Real>(exact.initialStep) == 1.0);
    CHECK(exact.step - 1 == low);
    CHECK(exact.initialStep - 1 == low);
    CHECK(exact.maxStep - 1 == low);
    REQUIRE(exact.biasPoints.size() == 2);
    REQUIRE(exact.vtkBiasPoints.size() == 1);
    CHECK(exact.biasPoints[0] - 1 == low);
    CHECK(exact.vtkBiasPoints[0] - 1 == low);
}

TEST_CASE("Automatic normalization retains native input precision", "[extended80]") {
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0., 0.});
    mesh.addNode(Node{1, 1., 0.}, {1.L + low, 0.L});
    mesh.addRegion(Region{0, "body", "Si"});
    DopingModel doping(mesh.numNodes());
    doping.setNodeDopingFromDecimal(1, nextToOneDecimal, "0");
    const auto temperature = 300.L + 0x1p-53L;
    UnitScalingReferenceConfig refs;
    refs.decimalTemperature_K = UnitScalingReferenceConfig::TemperatureInput{temperature, 300.};
    const auto automatic = UnitScalingSystem::autoInputsFrom(mesh, doping, MaterialDatabase{}, .5);
    const auto scales = UnitScalingSystem::fromInputs(
        300., 1., automatic, refs, PhysicalUnitSystem::tcadInternal());
    CHECK(scales.C0() == 1.);
    CHECK(scales.L0() == 1.);
    CHECK(scales.arithmeticValues().C0 - 1.L == low);
    CHECK(scales.arithmeticValues().L0 - 1.L == low);
    CHECK(scales.arithmeticValues().V0 == constants::thermalVoltage(temperature));
    CHECK(arithmeticThermalVoltageForInputs(300., refs) == scales.arithmeticValues().V0);
    CHECK(scales.V0() == thermalVoltageForInputs(300., refs));
    CHECK(scales.arithmeticValues().V0 != static_cast<long double>(scales.V0()));
}

TEST_CASE("Imported decimal coordinate survives inside the mesh and its edge length", "[extended80]") {
    const auto path = std::filesystem::temp_directory_path() /
        ("vela_extended_mesh_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".json");
    struct RemoveFile {
        const std::filesystem::path& path;
        ~RemoveFile() { std::error_code ec; std::filesystem::remove(path, ec); }
    } cleanup{path};
    {
        std::ofstream out(path);
        out << json{
            {"nodes", json::array({
                {{"id", 0}, {"x", 0.0}, {"y", 0.0}},
                {{"id", 1}, {"x", 1.0}, {"y", 0.0}, {"coordinate_literals", {nextToOneDecimal, "0"}}},
                {{"id", 2}, {"x", 0.0}, {"y", 1.0}}
            })},
            {"triangles", json::array({{{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 2}}}})},
            {"regions", json::array({{{"id", 0}, {"material", "Si"}, {"name", "body"}, {"cell_ids", {0}}}})},
            {"contacts", json::array()}
        };
    }
    const auto mesh = JsonMeshReader{}.read(path.string());
    CHECK_FALSE(mesh.hasDecimalCoordinates(0));
    REQUIRE(mesh.hasDecimalCoordinates(1));
    CHECK(mesh.arithmeticCoordinates(1)[0] - 1 == low);
    bool found = false;
    for (const auto& edge : mesh.edges()) {
        if (edge.n0 != 0 || edge.n1 != 1) continue;
        found = true;
        CHECK(edge.length == 1.0);
        CHECK(mesh.arithmeticEdgeLength(edge.id) - 1 == low);
    }
    CHECK(found);
}

TEST_CASE("Wide continuity and port derivatives share the SG edge contributions", "[extended80]") {
    DeviceMesh mesh;
    mesh.addNode(Node{0, 0., 0.});
    mesh.addNode(Node{1, 1., 0.});
    mesh.addNode(Node{2, 0., 1.});
    mesh.addRegion(Region{0, "gate", "PolySi", {0}});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.buildEdges();
    for (const auto& edge : mesh.edges())
        mesh.setTransportCouple(edge.id, 1.);
    auto geometry = std::make_shared<IalTransportGeometry>();
    geometry->edges.resize(mesh.numEdges());
    geometry->constantElectronEdges.assign(mesh.numEdges(), .1);
    geometry->constantHoleEdges.assign(mesh.numEdges(), .05);
    MobilityModelConfig config;
    config.model = "ialmob";
    config.carrierCurrentDiscretization = "scharfetter_gummel_edge";
    config.ialmobGeometry = geometry;
    const DopingModel doping(mesh.numNodes());
    const std::array<Index, 1> contact{0};
    ExtendedPopulationState state;
    state.potentialScale_V = .025L;
    state.concentrationScale_m3 = 1e6L;
    state.coordinates.resize(9);
    state.coordinates << 1.L, 1.1L, .95L, 10.L, 12.L, 9.L, 20.L, 19.L, 24.L;
    const ExtendedVector intrinsic = ExtendedVector::Constant(3, 1.L);
    const auto assemble = [&](const ExtendedPopulationState& frame,
                              std::optional<ExtendedRecombinationSupport> source = {}) {
        return assembleExtendedIalPopulationContinuity(config, mesh, doping, frame,
            intrinsic, intrinsic, 1.L, 1e6L, .025L, 300.L, false, true, {}, contact, source);
    };
    const auto baseline = assemble(state);
    ExtendedVector direction(9);
    direction << .2L, -.1L, .05L, -.2L, .4L, .3L, .1L, .5L, -1.L;
    constexpr long double h = 1e-5L;
    auto plus = state, minus = state;
    plus.coordinates += h * direction;
    minus.coordinates -= h * direction;
    const auto upper = assemble(plus), lower = assemble(minus);
    const ExtendedVector difference = (upper.residual - lower.residual) / (2 * h);
    const ExtendedVector tangent = baseline.jacobian * direction;
    CHECK((difference - tangent).norm() <= 1e-10L * tangent.norm());
    const auto factor = -constants::arithmeticConstants<long double>().q * 1e-6L;
    for (int carrier = 0; carrier < 2; ++carrier) {
        const auto rows = baseline.residual.segment(3 * carrier, 3);
        CHECK(std::abs(rows.sum()) <= 1e-18L * rows.cwiseAbs().sum());
        const auto current = factor * rows[0];
        CHECK(std::abs(baseline.current_A_per_um[carrier] - current) <= 1e-18L * std::abs(current));
        const auto change = (upper.current_A_per_um[carrier] - lower.current_A_per_um[carrier]) / (2 * h);
        const auto predicted = baseline.currentDerivative[carrier].dot(direction);
        CHECK(std::abs(change - predicted) <= 1e-10L * std::abs(predicted));
    }
    auto sourceConfig = recombinationModelConfig({"srh", "auger"}, .2, .4);
    sourceConfig.augerCn = .001;
    sourceConfig.augerCp = .002;
    const RecombinationModel sourceModel(sourceConfig);
    const std::array<long double, 3> areas{0, .25L, .5L};
    const ExtendedRecombinationSupport source{sourceModel, areas};
    const auto combined = assemble(state, source);
    CHECK(combined.residual[0] == baseline.residual[0]);
    CHECK(combined.residual[3] == baseline.residual[3]);
    const auto atNode = sourceModel.boltzmannPopulationLinearization(12.L, 19.L, 1.L, 0.L, 0.L);
    const auto expected = .25L * state.concentrationScale_m3 * (atNode.srh + atNode.auger);
    for (int row : {1, 4})
        CHECK(std::abs(combined.residual[row] - baseline.residual[row] - expected) <= 1e-18L * expected);
    CHECK(combined.current_A_per_um == baseline.current_A_per_um);
    const ExtendedVector combinedDifference =
        (assemble(plus, source).residual - assemble(minus, source).residual) / (2 * h);
    const ExtendedVector combinedTangent = combined.jacobian * direction;
    CHECK((combinedDifference - combinedTangent).norm() <= 1e-10L * combinedTangent.norm());
}

TEST_CASE("Wide SG value and direct density Jacobian share edge operands", "[extended80]") {
    const auto small = sgBoltzmannPopulationEdge(1.0L, 1.0L + low, 1.0L, 1.0L,
                                                1.0L, 1.0L, 1.0L, 1.0L, 1.0L, true);
    CHECK(std::abs(small.flux - low) <= low * 0x1p-58L);
    CHECK(sgBoltzmannPopulationEdge(1, 1, 1, 1, 1, 1, 1, 1, 1, true).flux == 0);
    const auto evaluate = [](long double psi1, long double n0, long double n1) {
        return sgBoltzmannPopulationEdge(-0.2L, psi1, n0, n1, 0.6L, 0.75L,
                                         -0.021L, 0.026L, 1.3L, true);
    };
    const auto value = evaluate(0.15L, 2, 3);
    const long double h = 0x1p-20L;
    CHECK(std::abs(value.potential1 - (evaluate(0.15L + h, 2, 3).flux -
                                       evaluate(0.15L - h, 2, 3).flux) / (2 * h)) < 1e-11L);
    CHECK(std::abs(value.population0 - (evaluate(0.15L, 2 + h, 3).flux -
                                        evaluate(0.15L, 2 - h, 3).flux) / (2 * h)) < 1e-11L);
    CHECK(std::abs(value.population1 - (evaluate(0.15L, 2, 3 + h).flux -
                                        evaluate(0.15L, 2, 3 - h).flux) / (2 * h)) < 1e-11L);
}

TEST_CASE("Wide Tri3 local differentiation keeps near-equilibrium increments", "[extended80]") {
    using Dual = detail::BasicTri3LocalForwardDual<ArithmeticReal>;
    const auto delta = Dual::variable(1.0L + low, 1) - Dual::variable(1.0L, 0);
    const auto source = detail::dualExpm1(delta);
    CHECK(std::abs(source.value - low) <= low * low);
    CHECK(source.derivative[0] == -(1 + low));
    CHECK(source.derivative[1] == 1 + low);
    for (std::size_t i = 2; i < source.derivative.size(); ++i)
        CHECK(source.derivative[i] == 0);
}

TEST_CASE("Wide Poisson constraints preserve the parsed contact through assembly and solve",
          "[extended80]") {
    DeviceMesh mesh;
    for (const Node node : {Node{0, 0., 0.}, Node{1, 1., 0.}, Node{2, 0., 1.}})
        mesh.addNode(node);
    Cell cell;
    cell.id = 0;
    cell.region_id = 0;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);
    Region region;
    region.id = 0;
    region.material = "Si";
    region.name = "body";
    region.cell_ids = {0};
    mesh.addRegion(region);
    mesh.buildEdges();

    const auto contacts = parseContactBoundarySpecs({
        {"contacts", json::array({{{"name", "drain"}, {"bias", 1.0}, {"flatband_voltage", 1.0}}})},
        {"decimal_stage_controls", {{"contacts", {{"drain", nextToOneDecimal}}}}}
    });
    const auto target = arithmeticPoissonDirichletPotential(contacts[0], *contacts[0].decimalInitialBias_V);
    MaterialDatabase materials;
    DopingModel doping(mesh.numNodes());
    PoissonScalingSpec scaling;
    scaling.enabled = true;
    scaling.potentialScale_V = 1;
    ArithmeticPoissonAssembler assembler(mesh, materials, doping, {}, {}, {}, scaling, 1.0L + low);
    assembler.assemble();
    assembler.applyDirichlet({{0, target}, {1, target}});
    CHECK(assembler.rhs()[0] == target / (1.0L + low));
    const auto solution = solveExtendedSparseLu(assembler.matrix(), assembler.rhs());
    CHECK(solution[0] == assembler.rhs()[0]);
    CHECK(solution[1] == assembler.rhs()[1]);
    CHECK(std::abs(solution[2] / solution[0] - 1.0L) <= 0x1p-61L);

    const ExtendedVector shifted = ExtendedVector::Constant(3, 1 + low);
    const ExtendedVector qf = ExtendedVector::Ones(3);
    const ExtendedVector intrinsic = ExtendedVector::Ones(3);
    const auto equilibrium = assembler.linearizeBoltzmannEquilibrium(
        shifted, qf, qf, intrinsic, intrinsic, 1.0L + low, {});
    CHECK(equilibrium.electrons[0] == 1 + low);
    CHECK(equilibrium.holes[0] == 1 - low);
    for (const auto& edge : mesh.edges()) {
        const auto flux = sgBoltzmannPopulationEdge(
            shifted[edge.n0], shifted[edge.n1], equilibrium.electrons[edge.n0],
            equilibrium.electrons[edge.n1], 1, 1, 1 + low, 1 + low, 1, false);
        CHECK(flux.flux == 0);
    }
}

TEST_CASE("Population avalanche keeps the sign of a sub-binary64 quasi-Fermi gradient",
          "[extended80]") {
    using Dual = detail::BasicTri3LocalForwardDual<ArithmeticReal>;
    BasicIalElementGeometry<ArithmeticReal> geometry;
    geometry.coordinates_m = {{{0, 0}, {1, 0}, {0, 1}}};
    geometry.touchesEffectiveElectrode = true;
    IalElementAvalancheOptions options;
    options.sourceMeasure_m2 = {0.25, 0.25, 0};
    options.minimumField_V_per_m = 0;
    options.electron = {1, 0, 1, 0};
    options.hole = options.electron;
    const IalAvalanchePopulationVertex zero{1, 2, 3, 1, 1};
    std::array state{zero, zero, zero};
    const std::array mobility{Dual(1), Dual(1)};
    state[1].potential += low;
    auto plus = evaluateIalPopulationAvalancheFromMobility(geometry, state, 1, 1, mobility, options);
    state[1].potential -= 2 * low;
    auto minus = evaluateIalPopulationAvalancheFromMobility(geometry, state, 1, 1, mobility, options);
    CHECK(plus.electron.currentField_V_per_m.value == low);
    CHECK(minus.electron.currentField_V_per_m.value == low);
    CHECK(plus.electron.currentField_V_per_m.derivative[4] == -1);
    CHECK(minus.electron.currentField_V_per_m.derivative[4] == 1);
    CHECK(plus.electron.ionizationField_V_per_m.derivative[4] == 0);
    CHECK(plus.electron.generation_m3_per_s.value == 2 * low);
}

TEST_CASE("Wide cold start hands off authoritative potential and equilibrium populations",
          "[extended80]") {
    DeviceMesh mesh;
    for (const Node node : {Node{0, 0., 0.}, Node{1, 1., 0.}, Node{2, 0., 1.}})
        mesh.addNode(node);
    Cell cell;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);
    Region region;
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);
    Contact terminal;
    terminal.name = "source";
    terminal.node_ids = {0};
    mesh.addContact(terminal);
    mesh.buildEdges();
    const auto contacts = parseContactBoundarySpecs({
        {"contacts", json::array({{{"name", "source"}, {"type", "ohmic"}, {"bias", 0}}})},
        {"decimal_stage_controls", {{"contacts", {{"source", "0"}}}}}
    });

    DopingModel doping(mesh.numNodes());
    doping.setNodeDoping(1, 1, 0);
    constexpr auto vt = 0.025L + low;
    ExtendedBoltzmannMaterials material;
    material.intrinsic = ExtendedVector::Ones(3);
    material.intrinsicPotential_V = ExtendedVector::Constant(3, vt * (1 + low));
    material.electronIntrinsic = ExtendedVector::Constant(3, std::exp(-(1 + low)));
    material.holeIntrinsic = ExtendedVector::Constant(3, std::exp(1 + low));
    const auto targets = buildExtendedZeroBiasOhmicTargets(mesh, doping, material, contacts, vt, true);
    REQUIRE(targets.potential[0].has_value());
    CHECK(static_cast<double>(*targets.potential[0]) == 1.);
    CHECK(*targets.potential[0] > 1.L);
    CHECK_FALSE(targets.potential[1].has_value());
    CHECK(targets.electrons[0] == 1.L);
    CHECK(targets.holes[0] == 1.L);

    MaterialDatabase database;
    PoissonScalingSpec scaling;
    scaling.enabled = true;
    scaling.permittivityReference_F_per_m = constants::eps0;
    ArithmeticPoissonAssembler poisson(mesh, database, doping, {}, {}, {}, scaling, vt);
    const auto limited = solveExtendedBoltzmannColdStart(poisson, doping, material, targets, vt, 1, 0, 1e-15L);
    CHECK_FALSE(limited.converged);
    const auto cold = solveExtendedBoltzmannColdStart(poisson, doping, material, targets, vt, 2, 12, 1e-15L);
    REQUIRE(cold.converged);
    CHECK(cold.updates > 0);
    CHECK(cold.state.potentialScale_V == vt);
    CHECK(cold.state.concentrationScale_m3 == 2.L);
    CHECK(cold.state.coordinates[0] == *targets.potential[0]);
    for (int i = 0; i < 3; ++i) {
        const auto psi = cold.state.coordinates[i];
        CHECK(cold.state.coordinates[3 + i] == material.electronIntrinsic[i] * std::exp(psi));
        CHECK(cold.state.coordinates[6 + i] == material.holeIntrinsic[i] * std::exp(-psi));
    }
    ExtendedIalPopulationContinuity continuity;
    continuity.residual = ExtendedVector::Zero(6);
    continuity.jacobian.resize(6, 9);
    continuity.jacobian.insert(0, 2) = -5;
    continuity.jacobian.insert(1, 2) = 2;
    const std::array<Index, 1> driven{0};
    const auto system = assembleExtendedPopulationDevice(poisson, cold.state, material, targets, driven, continuity);
    CHECK(system.residual.head(3).cwiseAbs().maxCoeff() <= 1e-15L);
    CHECK(system.jacobian.coeff(3, 2) == 0.L);
    CHECK(system.jacobian.coeff(4, 2) == 2.L);
    CHECK(system.voltageDerivative[0] == -1 / vt);
    CHECK(system.prescribedPopulationRows == std::vector<int>{3, 6});

    auto perturbed = cold.state;
    constexpr auto dn = 0.125L;
    perturbed.coordinates[4] += dn;
    const auto trial = assembleExtendedPopulationDevice(poisson, perturbed, material, targets, driven, continuity);
    CHECK(std::abs((trial.residual[1] - system.residual[1]) - system.jacobian.coeff(1, 4) * dn) <= 1e-20L);

    perturbed.coordinates[3] += 0.25L;
    continuity.residual[1] = 0.5L;
    continuity.jacobian.insert(1, 4) = 3.L;
    const auto physical = assembleExtendedPopulationDevice(poisson, perturbed, material, targets, driven, continuity);
    const ExtendedPopulationNormalization scales{8.L, 0.125L, 0.25L};
    const auto normalized = assembleExtendedPopulationDevice(
        poisson, perturbed, material, targets, driven, continuity, nullptr, &scales);
    CHECK(normalized.voltageDerivative[0] == physical.voltageDerivative[0]);
    CHECK(normalized.residual[1] == physical.residual[1] * scales.poissonRowFactor);
    CHECK(normalized.residual[3] == physical.residual[3] / scales.concentration_cm3);
    CHECK(normalized.residual[4] == physical.residual[4] * scales.continuityRowFactor);
    CHECK(normalized.jacobian.coeff(1, 4) ==
          physical.jacobian.coeff(1, 4) * scales.poissonRowFactor * scales.concentration_cm3);
    CHECK(normalized.jacobian.coeff(3, 3) == 1.L);
    CHECK(normalized.jacobian.coeff(4, 2) == physical.jacobian.coeff(4, 2) * scales.continuityRowFactor);
    CHECK(normalized.jacobian.coeff(4, 4) ==
          physical.jacobian.coeff(4, 4) * scales.continuityRowFactor * scales.concentration_cm3);
}

TEST_CASE("Extended fixed-bias Newton retains the direction and rejects an exhausted budget", "[extended80]") {
    ExtendedPopulationState seed;
    seed.coordinates = ExtendedVector::Ones(3);
    const auto target = 1.L + low;
    const auto assemble = [&](const ExtendedPopulationState& state) {
        ExtendedPopulationDeviceSystem frame;
        frame.residual = state.coordinates - ExtendedVector::Constant(3, target);
        frame.jacobian.resize(3, 3);
        frame.jacobian.setIdentity();
        frame.electronFree = {true};
        frame.holeFree = {true};
        return frame;
    };
    ExtendedFixedBiasNewtonConfig config{1, 3, 1, 1.L, 1e-10L, 1e-8L, {}, {}, true};
    const auto step = solveExtendedFixedBiasNewton(seed, config, assemble);
    REQUIRE(step.converged);
    CHECK(step.updates == 1);
    CHECK(step.factoredDimension == 3);
    CHECK(step.linearBackwardError == 0);
    CHECK(step.state.deviceBias_V == seed.deviceBias_V);
    CHECK(step.state.coordinates[0] == target);
    CHECK(step.state.coordinates[1] == target);
    CHECK(step.state.coordinates[2] == target);
    CHECK(step.updateError.valid);
    CHECK(step.updateError.rms > 0);

    auto distant = seed;
    distant.coordinates.array() += 1.L;
    const auto incomplete = solveExtendedFixedBiasNewton(distant, config, assemble);
    CHECK_FALSE(incomplete.converged);
    CHECK(incomplete.updates == 1);
    CHECK(incomplete.failure.find("iteration budget") != std::string::npos);
}

TEST_CASE("Extended resistor corrector accepts a coupled wide update and preserves its parent",
          "[extended80]") {
    ExtendedPopulationState seed;
    seed.coordinates = ExtendedVector::Ones(3);
    seed.supplyBias_V = low;
    const auto assemble = [](const ExtendedPopulationState& state) {
        ExtendedSeriesResistorFrame frame;
        frame.device.residual = state.coordinates - ExtendedVector::Ones(3);
        frame.device.residual[0] -= state.deviceBias_V;
        frame.device.jacobian.resize(3, 3);
        frame.device.jacobian.setIdentity();
        frame.device.electronFree = {true};
        frame.device.holeFree = {true};
        frame.device.voltageDerivative = ExtendedVector::Zero(3);
        frame.device.voltageDerivative[0] = -1;
        frame.portDerivative = ExtendedVector::Zero(3);
        frame.portDerivative[0] = 1;
        frame.current_A_per_um = state.coordinates[0] - 1;
        return frame;
    };
    const ExtendedFixedBiasNewtonConfig config{1, 0, 1, 1.L, 1e-10L, 1e-8L, {}, {}, false};
    const ExtendedSeriesResistorNewtonConfig circuit{3.L, low, 0.025852L, 1, 1, true};
    const auto result = solveExtendedSeriesResistorNewton(seed, config, circuit, assemble);
    REQUIRE(result.converged);
    CHECK(result.state.deviceBias_V == low / 4);
    CHECK(result.state.drainCurrent_A_per_um == low / 4);
    CHECK(result.state.coordinates[0] == 1 + low / 4);
    CHECK(result.updateError.counts[3] == 1);
    CHECK(seed.deviceBias_V == 0);
    CHECK(seed.coordinates[0] == 1);
    auto distant = seed;
    distant.supplyBias_V = 0.1L;
    CHECK_FALSE(solveExtendedSeriesResistorNewton(distant, config, circuit, assemble).converged);
    CHECK(distant.coordinates[0] == 1);
}

TEST_CASE("Extended resistor RHS growth includes the load line and respects report mode", "[extended80]") {
    ExtendedPopulationState seed;
    seed.coordinates = ExtendedVector::Ones(3);
    seed.supplyBias_V = 0.25L;
    const auto assemble = [](const ExtendedPopulationState& state) {
        ExtendedSeriesResistorFrame frame;
        frame.device.residual = ExtendedVector::Zero(3);
        frame.device.jacobian.resize(3, 3);
        frame.device.jacobian.setIdentity();
        frame.device.electronFree = {true};
        frame.device.holeFree = {true};
        frame.device.voltageDerivative = ExtendedVector::Zero(3);
        frame.portDerivative = ExtendedVector::Zero(3);
        frame.current_A_per_um = state.deviceBias_V == 0 ? 0.L : 0.5L;
        return frame;
    };
    ExtendedFixedBiasNewtonConfig config{1, 0, 0, 1.L, 1e-10L, 1e-8L, {}, {}, false};
    const ExtendedSeriesResistorNewtonConfig circuit{1.L, 1e-6L, 1e10L, 1, 1, false};
    const auto baseline = solveExtendedSeriesResistorNewton(seed, config, circuit, assemble);
    REQUIRE(baseline.converged);
    CHECK_FALSE(baseline.rhsGrowth.has_value());

    config.rhsGrowth = ExtendedRhsGrowthControl{1.5L, false};
    const auto report = solveExtendedSeriesResistorNewton(seed, config, circuit, assemble);
    REQUIRE(report.converged);
    REQUIRE(report.rhsGrowth.has_value());
    CHECK(report.rhsGrowth->initialNorm == 0.25L);
    CHECK(report.rhsGrowth->peakNorm == 0.5L);
    CHECK(report.rhsGrowth->firstExceededUpdate == 1);
    CHECK(report.state.deviceBias_V == baseline.state.deviceBias_V);
    CHECK(report.updateError.rms == baseline.updateError.rms);

    config.rhsGrowth->enforce = true;
    const auto rejected = solveExtendedSeriesResistorNewton(seed, config, circuit, assemble);
    CHECK_FALSE(rejected.converged);
    CHECK(rejected.failure == "rhs_growth_limit");
    CHECK(rejected.updates == 1);
    CHECK(seed.deviceBias_V == 0.L);
    config.rhsGrowth->factor = 2.L;
    const auto atLimit = solveExtendedSeriesResistorNewton(seed, config, circuit, assemble);
    CHECK(atLimit.converged);
    REQUIRE(atLimit.rhsGrowth.has_value());
    CHECK_FALSE(atLimit.rhsGrowth->firstExceededUpdate.has_value());
    config.rhsGrowth->factor = 0.5L;
    CHECK_THROWS_AS(solveExtendedSeriesResistorNewton(seed, config, circuit, assemble), std::invalid_argument);
}

TEST_CASE("Extended AverageBox transfers coordinate low bits to Poisson and mobile charge", "[extended80]") {
    const auto triangle = [](long double height) {
        DeviceMesh mesh;
        mesh.addNode({0, 0., 0.});
        mesh.addNode({1, 1., 0.});
        mesh.addNode({2, 0., 1.}, {0.L, height});
        Cell cell;
        cell.node_ids = {0, 1, 2};
        mesh.addCell(cell);
        Region region;
        region.material = "Si";
        region.cell_ids = {0};
        mesh.addRegion(region);
        BoxGeometryBuilder::Options policy;
        policy.poissonPermittivityPolicy = BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
        policy.cellBoxPolicy = BoxGeometryBuilder::CellBoxPolicy::AverageIntersection;
        mesh.buildBoxGeometry(policy);
        return mesh;
    };
    const auto flat = triangle(1.L);
    const auto raised = triangle(1.L + low);
    REQUIRE(flat.poissonCellEdgeCoefficients(0) == raised.poissonCellEdgeCoefficients(0));

    MaterialDatabase material;
    const auto a = buildExtendedAveragePoissonGeometry(flat, material, 300.L);
    const auto b = buildExtendedAveragePoissonGeometry(raised, material, 300.L);
    CHECK(a.chargeVolumes[0] == 0.25L);
    CHECK(b.chargeVolumes[0] > a.chargeVolumes[0]);
    DopingModel doping(3);
    doping.setNodeDoping(0, 1., 0.);
    PoissonScalingSpec spec;
    spec.enabled = true;
    spec.poissonChargeVolumePolicy = PoissonChargeVolumePolicy::MaterialLocalBox;
    ArithmeticPoissonAssembler flatPoisson(flat, material, doping, {}, {}, {}, spec);
    ArithmeticPoissonAssembler raisedPoisson(raised, material, doping, {}, {}, {}, spec);
    flatPoisson.assemble(a);
    raisedPoisson.assemble(b);
    CHECK(raisedPoisson.matrix().coeff(0, 1) < flatPoisson.matrix().coeff(0, 1));
    CHECK(raisedPoisson.rhs()[0] > flatPoisson.rhs()[0]);
    CHECK(raisedPoisson.mobileChargeWeight(0, b) > flatPoisson.mobileChargeWeight(0, a));
}

TEST_CASE("Extended transport CSV and material bindings preserve sub-binary64 input changes",
          "[extended80]") {
    const auto folder = std::filesystem::temp_directory_path() /
        ("vela_extended_slot_weights_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    REQUIRE(std::filesystem::create_directory(folder));
    struct RemoveDirectory {
        const std::filesystem::path& path;
        ~RemoveDirectory() { std::error_code error; std::filesystem::remove_all(path, error); }
    } cleanup{folder};
    DeviceMesh mesh;
    for (const Node node : {Node{0, 0., 0.}, Node{1, 1., 0.}, Node{2, 0., 1.}})
        mesh.addNode(node);
    Cell cell;
    cell.node_ids = {0, 1, 2};
    mesh.addCell(cell);
    Region region;
    region.material = "Si";
    region.cell_ids = {0};
    mesh.addRegion(region);
    mesh.buildEdges();
    {
        std::ofstream csv(folder / "couples.csv");
        csv << "node0,node1,couple_m\n0,1,0.000001000000000000000001\n0,2,0\n";
    }
    const json deck{{"mesh_geometry", {
        {"carrier_transport_couple_profile", "slot_ldmos_external_averagebox"},
        {"external_averagebox_couples_file", "couples.csv"}, {"external_averagebox_expected_edges", 2}}}};
    std::vector<long double> wide;
    const auto report = applyCarrierTransportCoupleProfile(mesh, deck, folder, UnitScalingConfig{}, wide);
    REQUIRE(wide.size() == mesh.numEdges());
    CHECK(report.records == 2);
    CHECK(std::ranges::count(wide, 0.L) == 1);
    CHECK(std::ranges::count(wide, -1.L) == 1);
    for (const auto& edge : mesh.edges())
        if (edge.transport_couple > 0) {
            CHECK(wide[edge.id] > 1e-6L);
            CHECK(static_cast<Real>(wide[edge.id]) == edge.transport_couple);
        }

    MaterialDatabase materials;
    const auto materialFile = folder / "material.json";
    {
        std::ofstream out(materialFile);
        out << R"({"materials":[{"name":"Si","eps_r":11.7,"decimal_eps_r":"11.700000000000000001"}]})";
    }
    CHECK(readExtendedPermittivities(materialFile.string(), materials).at("Si") == 11.700000000000000001L);
    {
        std::ofstream out(materialFile);
        out << R"({"materials":[{"name":"Si","eps_r":11.7,"decimal_eps_r":"11.8"}]})";
    }
    CHECK_THROWS_AS(readExtendedPermittivities(materialFile.string(), materials), std::invalid_argument);
}

TEST_CASE("Extended IALMob interface distance preserves the mesh coordinate low bits", "[extended80]") {
    const auto axis = std::array<long double, 3>{1, 0, 0};
    const auto vertical = std::array<long double, 3>{0, 1, 0};
    const std::vector<std::array<long double, 2>> xy{{0,0}, {1,0}, {0,1 + low}};
    const std::vector<ExtendedIalInterfaceSegment> interface{{0, 1, {0, 1}}};
    const auto geometry = buildIalInterfaceGeometry(xy, interface, axis, vertical);
    CHECK(geometry[2].distance_m == 1 + low);
    CHECK(geometry[2].orientationFamily == 100);
    CHECK(static_cast<Real>(geometry[2].distance_m) == 1);
}

TEST_CASE("IALMob shared-edge weights retain geometry input below a double ULP", "[extended80]") {
    const auto folder = std::filesystem::temp_directory_path() /
        ("vela_extended_ial_geometry_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    REQUIRE(std::filesystem::create_directory(folder));
    struct RemoveDirectory {
        std::filesystem::path path;
        ~RemoveDirectory() { std::error_code error; std::filesystem::remove_all(path, error); }
    } cleanup{folder};
    DeviceMesh mesh;
    for (const auto node : {Node{0,0,0}, Node{1,1,0}, Node{2,0,1}, Node{3,1,1}})
        mesh.addNode(node);
    for (const auto nodes : {std::vector<Index>{0,1,2}, std::vector<Index>{2,1,3}}) {
        Cell cell;
        cell.node_ids = nodes;
        mesh.addCell(cell);
    }
    mesh.buildEdges();
    auto options = std::make_shared<IalTransportOptions>();
    options->geometryFile = (folder / "transport.json").string();
    auto bound = std::make_shared<IalTransportGeometry>();
    bound->edges.resize(mesh.numEdges());
    bound->cells.resize(2);
    bound->cells[0].cellId = 0;
    bound->cells[1].cellId = 1;
    for (auto& cell : bound->cells)
        cell.geometry.vertexMeasure_m2 = {1,1,1};
    const auto shared = std::ranges::find_if(mesh.edges(), [](const auto& edge) {
        return std::min(edge.n0, edge.n1) == 1 && std::max(edge.n0, edge.n1) == 2;
    });
    REQUIRE(shared != mesh.edges().end());
    bound->edges[shared->id] = {{0, .5}, {1, .5}};
    {
        std::ofstream out(options->geometryFile);
        out << R"({"cells":[
            {"cell_id":0,"node_ids":[0,1,2],"vertex_measure_m2":[1.000000000000000005,1,1],
             "edge_coefficients":[0,1.000000000000000005,0]},
            {"cell_id":1,"node_ids":[2,1,3],"vertex_measure_m2":[1,1,1],"edge_coefficients":[1,0,0]}
        ]})";
    }
    MobilityModelConfig config;
    config.ialmob = options;
    config.ialmobGeometry = bound;
    const auto exact = readExtendedIalTransportMeasures(config, mesh);
    CHECK(exact.vertex[0][0] == 1.000000000000000005L);
    REQUIRE(exact.edgeWeights[shared->id].size() == 2);
    CHECK(exact.edgeWeights[shared->id][0] > .5L);
    CHECK(exact.edgeWeights[shared->id][1] < .5L);
    CHECK(static_cast<double>(exact.edgeWeights[shared->id][0]) == .5);
    CHECK(static_cast<double>(exact.edgeWeights[shared->id][1]) == .5);
}

TEST_CASE("Material lookup preserves bound operating-temperature PAR values", "[extended80]") {
    const auto path = std::filesystem::temp_directory_path() /
        ("vela_extended_material_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".json");
    struct RemoveFile {
        const std::filesystem::path& path;
        ~RemoveFile() { std::error_code ec; std::filesystem::remove(path, ec); }
    } cleanup{path};
    const auto parameter = [](const char* value) {
        return json{{"literals", json::array({value})}};
    };
    {
        std::ofstream out(path);
        out << json{
            {"physical_constants_profile", constants::profileName},
            {"materials", json::array({{
                {"name", "PolySi"}, {"temperature_K", 300.0}, {"intrinsic_reference_potential_V", 0.125},
                {"decimal_semiconductor", {
                    {"temperature_K", "300"}, {"fold_old_slotboom", false},
                    {"sections", {
                        {"Bandgap", {{"Tpar", parameter("300")}, {"beta", parameter("1")},
                                     {"alpha", parameter("0")}, {"Eg0", parameter(nextToOneDecimal)},
                                     {"Chi0", parameter(nextToOneDecimal)}}},
                        {"eDOSMass", {{"Formula", parameter("1")}, {"a", parameter("0.1")},
                                      {"ml", parameter("1")}, {"mm", parameter("0")}}},
                        {"hDOSMass", {{"Formula", parameter("1")}, {"a", parameter("1")}}}
                    }}
                }}
            }})}
        };
    }
    MaterialDatabase db(path.string());
    const auto* exact = db.decimalSemiconductor("PolySi", 300);
    REQUIRE(exact);
    CHECK(exact->bandgap_eV - 1 == low);
    CHECK(exact->affinity_eV - 1 == low);
    CHECK(exact->boundReferencePotential_V == 0.125L);
    CHECK(exact->electronsDOS_internal > 0);
    CHECK(db.getMaterial("PolySi").bandgap_eV == 1.0);
    CHECK(db.decimalSemiconductor("PolySi", std::nextafter(300.0L, 301.0L)) == nullptr);

    DeviceMesh mesh;
    mesh.addNode(Node{0, 0., 0.});
    mesh.addNode(Node{1, 1., 0.});
    mesh.addNode(Node{2, 0., 1.});
    mesh.addNode(Node{3, 1., 1.});
    mesh.addRegion(Region{0, "oxide", "SiO2", {0}});
    mesh.addRegion(Region{1, "gate", "PolySi", {1}});
    mesh.addCell(Cell{0, CellType::Tri3, 0, {0, 1, 2}});
    mesh.addCell(Cell{1, CellType::Tri3, 1, {1, 3, 2}});
    DopingModel doping(mesh.numNodes());
    doping.setNodeDopingFromDecimal(2, "1e23", "1e22");
    const auto bgn = bandgapNarrowingConfig("old_slotboom");
    const auto wide = buildExtendedBoltzmannMaterials(mesh, db, doping, bgn, 300.L, .025L);
    CHECK(wide.intrinsic[0] == 0);
    CHECK(wide.intrinsic[2] > wide.intrinsic[1]);
    CHECK(wide.intrinsic[1] == exact->intrinsicDensity_internal);
    CHECK(wide.intrinsicPotential_V[2] == exact->boundReferencePotential_V);
    CHECK(wide.electronDOS[2] == exact->electronsDOS_internal);
    const auto intrinsicSquared = wide.intrinsic[2] * wide.intrinsic[2];
    CHECK(std::abs(wide.electronIntrinsic[2] * wide.holeIntrinsic[2] / intrinsicSquared - 1) < 1e-18L);
    for (const auto net : {50.L, -50.L}) {
        const auto target = extendedBoltzmannOhmicTarget(net, 2.L, .125L, .025L, 1.L, 1.L, 1.L);
        CHECK(std::abs(target.electrons - target.holes - net) < 1e-17L);
        CHECK(std::abs(target.electrons * target.holes - 4) < 1e-18L);
        const auto shifted = extendedBoltzmannOhmicTarget(net, 2.L, .125L, .025L, 1.L + low, 1.L, 1.L);
        CHECK(shifted.potential - target.potential == low);
        CHECK(shifted.electrons == target.electrons);
        CHECK(shifted.holes == target.holes);
    }
    const auto bennett = bandgapNarrowingConfig("bennett_wilson");
    const ExtendedBandgapCoefficients coefficients{3.162e24L, 6.84e-3L, 0, 0};
    // Strong compensation puts the total above the threshold while leaving zero net charge.
    doping.setNodeDopingFromDecimal(2, "5e24", "5e24");
    doping.setNodeDopingFromDecimal(3, "1.581e24", "1.581e24");
    const auto narrowed = buildExtendedBoltzmannMaterials(mesh, db, doping, bennett, 300.L, .025L, coefficients);
    const long double logRatio = std::log(1e25L / 3.162e24L);
    const long double expected = std::exp(6.84e-3L * logRatio * logRatio / .05L);
    CHECK(std::abs(narrowed.intrinsic[2] / exact->intrinsicDensity_internal - expected) < 1e-18L);
    CHECK(narrowed.intrinsic[3] == exact->intrinsicDensity_internal);
    const auto disabled = buildExtendedBoltzmannMaterials(
        mesh, db, doping, bandgapNarrowingConfig("none"), 300.L, .025L);
    CHECK(disabled.intrinsic[2] == exact->intrinsicDensity_internal);

    const auto replacement = db.getMaterial("PolySi");
    db.addMaterial(replacement);
    CHECK(db.decimalSemiconductor("PolySi", 300) == nullptr);
    CHECK_THROWS_AS(buildExtendedBoltzmannMaterials(mesh, db, doping, bgn, 300.L, .025L),
                    std::invalid_argument);
}

TEST_CASE("Wide voltage step controller retains low bits across shrink and output landings",
          "[extended80]") {
    detail::BasicDCSweepStepControlConfig<ArithmeticReal> config;
    config.start = 1.0L + low;
    config.stop = 2.0L + low;
    config.step = config.initialStep = config.maxStep = 0.5L;
    config.minStep = 0.125L;
    config.mandatoryPoints = {1.75L + low};
    std::vector<ArithmeticReal> accepted, attempted;
    detail::runTypedDCSweepStepControl(config,
        [&](ArithmeticReal bias, ArithmeticReal increment, int retry) {
            attempted.push_back(bias);
            if (attempted.size() == 1) {
                CHECK(increment == 0.5L);
                return false;
            }
            if (attempted.size() == 2) {
                CHECK(retry == 1);
                CHECK(increment == 0.25L);
            }
            return true;
        }, [&](const detail::BasicDCSweepStepControlEvent<ArithmeticReal>& event) {
            if (event.converged) accepted.push_back(event.voltage);
        });
    const std::vector<ArithmeticReal> expected{
        1.25L + low, 1.5L + low, 1.75L + low, 2.0L + low};
    REQUIRE(accepted == expected);
    CHECK(attempted.front() == 1.5L + low);
}

TEST_CASE("Extended sweep lands on distinct checkpoints below a double ULP", "[extended80]") {
    detail::BasicDCSweepStepControlConfig<long double> config;
    config.start = 1.L;
    config.stop = 1.L + 8 * low;
    config.step = config.initialStep = config.maxStep = 8 * low;
    config.minStep = low;
    config.mandatoryPoints = {1.L + low};
    std::vector<long double> accepted;
    detail::runTypedDCSweepStepControl<long double>(config,
        [](long double, long double, int) { return detail::DCSweepStepAttemptResult(true); },
        [&](const auto& event) {
            if (event.converged)
                accepted.push_back(event.voltage);
        });
    CHECK((accepted == std::vector<long double>{1.L + low, config.stop}));
}

TEST_CASE("Extended sweep evaluates the exact decimal endpoint after nominal step rounding",
          "[extended80]") {
    detail::BasicDCSweepStepControlConfig<long double> config;
    config.start = 0.L;
    config.stop = 0.10L;
    // A rounded nominal target just above stop must never be sent to Newton.
    config.step = std::nextafter(config.stop, 1.L);
    config.initialStep = config.maxStep = 0.06L;
    config.minStep = 0.03L;
    config.stopBelowMinimum = true;

    for (const bool rejectEndpoint : {false, true}) {
        CAPTURE(rejectEndpoint);
        std::vector<long double> attempted, accepted;
        std::string failureReason;
        detail::runTypedDCSweepStepControl<long double>(config,
            [&](long double bias, long double, int) {
                attempted.push_back(bias);
                return detail::DCSweepStepAttemptResult(!rejectEndpoint || bias <= config.initialStep);
            },
            [&](const auto& event) {
                if (event.converged)
                    accepted.push_back(event.voltage);
                else
                    failureReason = event.failureReason;
            });

        CHECK((attempted == std::vector<long double>{config.initialStep, config.stop}));
        if (rejectEndpoint) {
            CHECK(accepted == std::vector<long double>{config.initialStep});
            CHECK(failureReason == "min_step_exhausted");
        } else {
            CHECK(accepted == attempted);
            CHECK(failureReason.empty());
        }
    }
}

TEST_CASE("Extended mesh identity detects changes below a double ULP", "[extended80]") {
    DeviceMesh a, b;
    const Node compatible{0, 1.0, 0.0};
    a.addNode(compatible, {1.0L, 0.0L});
    b.addNode(compatible, {1.0L + low, 0.0L});
    CHECK(stateMeshIdentity(a, 1.0) == stateMeshIdentity(b, 1.0));
    CHECK(extendedStateMeshIdentity(a, 1.0L) != extendedStateMeshIdentity(b, 1.0L));
    CHECK(extendedStateMeshIdentity(a, 1.0L) != extendedStateMeshIdentity(a, 1.0L + low));
}

TEST_CASE("SLOT IALMob low-field kernel consumes the decimal PAR parameters", "[extended80]") {
    const json rounded{{"mumax", 1417.}};
    const auto standard = ial_json::extendedParameters(rounded, true, {{"mumax", "1417"}});
    const auto perturbed = ial_json::extendedParameters(
        rounded, true, {{"mumax", "1417.00000000000008"}});
    REQUIRE(static_cast<Real>(perturbed.muMax) == static_cast<Real>(standard.muMax));
    constexpr auto method = IalScreeningMethod::U202212SingleNewton;
    const BasicIalMobilityState<long double> empty;
    const auto first = evaluateExtendedSlotIalMobility(empty, standard, true, method);
    const auto second = evaluateExtendedSlotIalMobility(empty, perturbed, true, method);
    CHECK(second.result.mobility_m2_per_Vs > first.result.mobility_m2_per_Vs);

    BasicIalMobilityState<long double> state{1e23L, 4e21L, 2e23L, 4e14L, 3e6L, 2e-9L, 300.L};
    const IalMobility legacy(IalMobility::siliconDefaults(true), true, IalScreeningMethod::U202212SingleNewton);
    const auto legacyRate = legacy.evaluateWithDerivatives({1e23, 4e21, 2e23, 4e14, 3e6, 2e-9, 300.});
    const auto wide = evaluateExtendedSlotIalMobility(state, standard, true, method);
    CHECK(std::abs(wide.result.mobility_m2_per_Vs / legacyRate.result.mobility_m2_per_Vs - 1.L) < 1e-12L);
    constexpr auto delta = 1.L;
    state.normalField_V_per_m += delta;
    const auto plus = evaluateExtendedSlotIalMobility(state, standard, true, method);
    state.normalField_V_per_m -= 2.L * delta;
    const auto minus = evaluateExtendedSlotIalMobility(state, standard, true, method);
    const auto finiteDifference = (plus.result.mobility_m2_per_Vs - minus.result.mobility_m2_per_Vs) / (2.L * delta);
    CHECK(std::abs(finiteDifference / wide.derivative_SI[4] - 1.L) < 1e-9L);
}

TEST_CASE("Wide IALMob and avalanche share the direct density perturbation", "[extended80]") {
    BasicIalElementGeometry<ArithmeticReal> geometry;
    geometry.coordinates_m = {{{0,0},{1e-8L,0},{0,2e-8L}}};
    geometry.interfaceDistance_m = {2e-9L, 2e-9L, 3e-9L};
    geometry.touchesEffectiveElectrode = true;
    const auto electron = ialSiliconDefaults<long double>(true);
    const auto hole = ialSiliconDefaults<long double>(false);
    const std::array electrons{&electron, &electron, &electron};
    const std::array holes{&hole, &hole, &hole};
    const std::array<std::array<ArithmeticReal,2>,3> doping{{
        {1e23L, 1e21L}, {1e23L, 1e21L}, {1e23L, 1e21L}}};
    std::array<IalPopulationVertex,3> state{{
        {21.L, 1e23L, 1e12L, 1e16L, 1e16L},
        {23.L, 1.09e23L, 1.02e12L, 1e16L, 1e16L},
        {26.L, 1.02e23L, 1.01e12L, 1e16L, 1e16L}}};
    constexpr auto voltageScale = 0.02585L;
    const std::array<IalHighFieldDecimal,2> highField{{{1.07e5L,1.109L},{8.37e4L,1.213L}}};
    IalElementAvalancheOptions avalanche;
    avalanche.sourceMeasure_m2 = {2e-17,3e-17,4e-17};
    const auto evaluate = [&](const auto& x, const std::array<long double, 3>* shares = nullptr,
                              const ArithmeticIalIonizationCoefficients* coefficients = nullptr) {
        const auto fields = detail::buildIalPopulationFields(x, voltageScale, voltageScale);
        const auto mobility = evaluateIalPopulationMobility(geometry, fields, doping, electrons, holes,
            300.L, 0.L, highField, GradQfContactField::ElectricField, IalScreeningMethod::U202212SingleNewton);
        return evaluateIalPopulationAvalancheFromMobility(geometry, x, voltageScale, voltageScale,
            {mobility.electron[1],mobility.hole[1]}, avalanche, shares, coefficients);
    };
    const auto initial = evaluate(state);
    const auto asArithmetic = [](const auto& c) -> ArithmeticElementIonizationCurve {
        return {c.lowA_per_m, c.lowB_V_per_m, c.highA_per_m, c.highB_V_per_m, c.switchField_V_per_m};
    };
    ArithmeticIalIonizationCoefficients coefficients{asArithmetic(avalanche.electron),
        asArithmetic(avalanche.hole), avalanche.minimumField_V_per_m};
    coefficients.electron.lowA_per_m *= 1 + 32 * low;
    coefficients.electron.highA_per_m *= 1 + 32 * low;
    REQUIRE(static_cast<double>(coefficients.electron.lowA_per_m) == avalanche.electron.lowA_per_m);
    const auto altered = evaluate(state, nullptr, &coefficients);
    CHECK(altered.electron.alpha_per_m.value > initial.electron.alpha_per_m.value);
    CHECK(altered.hole.alpha_per_m.value == initial.hole.alpha_per_m.value);
    std::array<long double, 3> shares{avalanche.sourceMeasure_m2[0],
        avalanche.sourceMeasure_m2[1], avalanche.sourceMeasure_m2[2]};
    shares[1] *= 1.L + 32 * std::numeric_limits<long double>::epsilon();
    REQUIRE(static_cast<double>(shares[1]) == avalanche.sourceMeasure_m2[1]);
    auto exact = evaluate(state, &shares);
    CHECK(exact.electron.source_per_m_s[1].value > initial.electron.source_per_m_s[1].value);
    CHECK(exact.electron.source_per_m_s[0].value == initial.electron.source_per_m_s[0].value);
    CHECK(exact.electron.source_per_m_s[1].derivative[4] != initial.electron.source_per_m_s[1].derivative[4]);
    const auto previous = state[1].electrons_m3;
    constexpr auto step = 1e-8L;
    state[1].electrons_m3 = previous * (1 + step);
    const auto plus = evaluate(state);
    state[1].electrons_m3 = previous * (1 - step);
    const auto minus = evaluate(state);
    const auto fd = (plus.electron.source_per_m_s[1].value -
                     minus.electron.source_per_m_s[1].value) / (2 * step);
    const auto derivative = initial.electron.source_per_m_s[1].derivative[4];
    REQUIRE(derivative != 0.L);
    CHECK(std::abs(fd/derivative - 1.L) < 1e-6L);
    CHECK(initial.electron.ionizationField_V_per_m.derivative[4] == 0.L);
}

TEST_CASE("The physical load line retains low bits in the RHS and voltage column", "[extended80]") {
    const auto wide = ArithmeticLoadLineResidual::evaluate(1.L, 1.L+low, 1.L, -1.L, 0.L);
    CHECK(wide.voltage_V == low);
    CHECK(wide.physicalRhs() == -low);
    CHECK(wide.inVoltageUnits(0.5L) == 2.L*low);
    CHECK(ArithmeticLoadLineResidual::voltageDerivative(1.L, 1.L+low, -1.L, {}) == -low);
    CHECK(LoadLineResidual::evaluate(1., static_cast<Real>(1.L+low), 1., -1., 0.).voltage_V == 0.);
}

TEST_CASE("The scaled bordered matrix preserves a voltage increment below double precision", "[extended80]") {
    ExtendedSparseMatrix device(1,1);
    device.insert(0,0) = 1.L;
    const ExtendedVector residual = ExtendedVector::Constant(1,-1.L);
    const ExtendedVector potentialColumn = ExtendedVector::Ones(1);
    const ExtendedVector currentDerivative = ExtendedVector::Zero(1);
    std::vector<int> prescribed;
    std::string failure;
    const auto system = detail::buildBorderedLoadLineSystem(device, residual,
        potentialColumn, currentDerivative, 1.L, -low, ExtendedVector::Ones(2),
        prescribed, failure);
    REQUIRE(system);
    CHECK(prescribed == std::vector{1});
    CHECK(system->rawRhs[1] == low);
    const auto retained = solveExtendedBorderedSystem(*system, {}, prescribed, 2);
    REQUIRE(retained.factoredDimension == 2);
    CHECK(retained.update[0] == 1.L-low);
    CHECK(retained.update[1] == low);
    CHECK(retained.backwardError == 0.L);
    CHECK(retained.scaledResidualInfinityNorm == 0.L);
    const auto eliminated = solveExtendedBorderedSystem(*system, prescribed, {}, 2);
    REQUIRE(eliminated.factoredDimension == 1);
    CHECK(eliminated.update == retained.update);
}

TEST_CASE("The native update metric retains the side of the strict E threshold", "[extended80]") {
    DensityUpdateAccumulator<ArithmeticReal> wide(0);
    wide.add(0, 0.L, 1.L-low, 1.L);
    const auto result = wide.finish();
    REQUIRE(result.valid);
    CHECK(result.rms == 1.L-low);
    CHECK(result.rms < 1.L);

    DensityUpdateAccumulator<Real> compatible(0);
    compatible.add(0, 0., static_cast<Real>(1.L-low), 1.);
    REQUIRE(compatible.finish().valid);
    CHECK(compatible.finish().rms == 1.);
}

TEST_CASE("Separate dopant operands retain low bits for mobility and net charge", "[extended80]") {
    DopingModel doping(2);
    doping.setNodeDoping(1,8.,3.);
    doping.setNodeDopingFromDecimal(0,nextToOneDecimal,"1");
    CHECK(doping.donors(0) == doping.acceptors(0));
    CHECK(doping.dopingOperandsForAssembly(0)[0] == 1.L+low);
    CHECK(doping.dopingOperandsForAssembly(0)[1] == 1.L);
    CHECK(doping.netDopingForAssembly(0) == low);
    CHECK((doping.dopingOperandsForAssembly(1) == std::array<ArithmeticReal,2>{8.L,3.L}));
    doping.setNodeDoping(0,4.,2.);
    CHECK((doping.dopingOperandsForAssembly(0) == std::array<ArithmeticReal,2>{4.L,2.L}));
    CHECK(doping.netDopingForAssembly(0) == 2.L);
}

TEST_CASE("Direct concentration SRH and Auger retain the side of equilibrium", "[extended80]") {
    RecombinationModelConfig cfg;
    cfg.mechanisms = {"srh","auger"};
    cfg.taun = 3.;
    cfg.taup = 2.;
    cfg.augerCn = 1.;
    cfg.augerCp = 2.;
    cfg.augerExcessProduct = "classical_np";
    cfg.augerWithGeneration = false;
    const RecombinationModel model(cfg);
    const auto evaluate = [&](ArithmeticReal n, ArithmeticReal p) {
        return model.boltzmannPopulationLinearization(n,p,1.L,0.L,0.L);
    };
    const auto above = evaluate(1.L,1.L+low);
    const auto below = evaluate(1.L,1.L-low);
    CHECK(above.srh > 0.L);
    CHECK(above.auger > 0.L);
    CHECK(below.srh < 0.L);
    CHECK(below.auger == 0.L);
    CHECK(below.dAugerDn == 0.L);
    CHECK(below.dAugerDp == 0.L);

    constexpr auto h = 1e-6L;
    const auto at = evaluate(1.3L,1.4L);
    ArithmeticBoltzmannSourceCoefficients bound{cfg.taun, cfg.taup, {},
        {cfg.augerCn, cfg.augerCp}, {}, {1, 1}};
    bound.taup = 2.L + 32 * low;
    bound.augerCoefficient[0] = 1.L + 32 * low;
    REQUIRE(static_cast<double>(bound.taup) == cfg.taup);
    REQUIRE(static_cast<double>(bound.augerCoefficient[0]) == cfg.augerCn);
    const auto authored = model.boltzmannPopulationLinearization(1.3L, 1.4L, 1.L, 0.L, 0.L, bound);
    CHECK(authored.srh < at.srh);
    CHECK(authored.auger > at.auger);
    CHECK(authored.dAugerDp != at.dAugerDp);
    const auto plus = evaluate(1.3L,1.4L+h);
    const auto minus = evaluate(1.3L,1.4L-h);
    CHECK(std::abs((plus.srh-minus.srh)/(2.L*h)-at.dSrhDp) < 1e-12L);
    CHECK(std::abs((plus.auger-minus.auger)/(2.L*h)-at.dAugerDp) < 1e-10L);
}
