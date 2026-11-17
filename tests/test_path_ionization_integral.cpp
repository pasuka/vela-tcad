#include "vela/post/PathIonizationIntegral.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

using namespace vela;

TEST_CASE("continuous ABA peak window keeps every eligible cell",
          "[impact][path_ionization]")
{
    std::vector<CellIonizationSample> cells{
        {10, {0, 1, 2}, {Point2{0, 0}, Point2{1, 0}, Point2{0, 1}},
         {0, -1, 0}, Point2{10, 0}, 1, 1, 2, 1},
        {11, {1, 3, 2}, {Point2{1, 0}, Point2{1, 1}, Point2{0, 1}},
         {-1, -1, 0}, Point2{9.5, 0}, 1, 1, 2, 1},
        {12, {4, 5, 6}, {Point2{2, 0}, Point2{3, 0}, Point2{2, 1}},
         {0, -1, 0}, Point2{8, 0}, 1, 1, 2, 1},
    };
    const auto result = analyzeContinuousCellIonizationPaths(
        7, cells, 0, 0., 0.,
        ContinuousPathDirection::Bidirectional, {},
        ContinuousPathRetention::AllSeedTrajectories,
        ContinuousPathSeedMode::CellLocalMaxima,
        ContinuousPathTracingPolicy::CellElectricField, 0., 1.);
    REQUIRE(result.paths.size() == 2);
    for (const auto& path : result.paths) {
        REQUIRE((path.seedCellId == 10 || path.seedCellId == 11));
        REQUIRE(path.seedField >= 9.);
    }
}

TEST_CASE("Sentaurus path ionization integral preserves equal-alpha limit",
          "[impact][path-ionization]")
{
    const std::vector<PathIonizationSegment> path = {
        {0, 0, 1, 0.4, 2.0, 1.0, 1.0},
        {1, 1, 2, 0.6, 2.0, 1.0, 1.0},
    };
    const PathIonizationIntegral value = integrateIonizationPath(path);

    const Real expectedElectron = 2.0 * (1.0 - std::exp(-1.0));
    const Real expectedHole = std::exp(1.0) - 1.0;
    REQUIRE(value.electron == Catch::Approx(expectedElectron).epsilon(1.0e-13));
    REQUIRE(value.hole == Catch::Approx(expectedHole).epsilon(1.0e-13));
    REQUIRE(value.mean == Catch::Approx(0.5 * (expectedElectron + expectedHole)));
}

TEST_CASE("Sentaurus path integral uses stable zero-difference segment limit",
          "[impact][path-ionization]")
{
    const std::vector<PathIonizationSegment> path = {
        {0, 0, 1, 0.25, 3.0, 3.0, 1.0},
        {1, 1, 2, 0.75, 3.0, 3.0, 1.0},
    };
    const PathIonizationIntegral value = integrateIonizationPath(path);
    REQUIRE(value.electron == Catch::Approx(3.0).epsilon(1.0e-14));
    REQUIRE(value.hole == Catch::Approx(3.0).epsilon(1.0e-14));
    REQUIRE(value.mean == Catch::Approx(3.0).epsilon(1.0e-14));
}

TEST_CASE("Carrier path supports terminate independently around their driving-field maxima",
          "[impact][path-ionization][carrier-stop]")
{
    std::vector<PathIonizationSegment> path(3);
    path[0].length = 0.2;
    path[1].length = 0.3;
    path[2].length = 0.5;
    for (PathIonizationSegment& segment : path) {
        segment.electronAlpha = 2.0;
        segment.holeAlpha = 1.0;
    }
    path[0].electronDrivingField = 1.0;
    path[1].electronDrivingField = 5.0;
    path[2].electronDrivingField = 1.0;
    path[0].holeDrivingField = 4.0;
    path[1].holeDrivingField = 4.0;
    path[2].holeDrivingField = 0.5;

    PathIonizationIntegrationOptions options;
    options.electronMinimumDrivingField = 2.0;
    options.holeMinimumDrivingField = 2.0;
    options.meanDefinition =
        MeanIonizationDefinition::CarrierAlphaLengthArithmetic;
    const PathIonizationIntegral value = integrateIonizationPath(path, options);

    const PathIonizationIntegral electronReference =
        integrateIonizationPath({path[1]});
    const PathIonizationIntegral holeReference =
        integrateIonizationPath({path[0], path[1]});
    REQUIRE(value.electron ==
            Catch::Approx(electronReference.electron).epsilon(1.0e-13));
    REQUIRE(value.hole ==
            Catch::Approx(holeReference.hole).epsilon(1.0e-13));
    REQUIRE(value.electronSupportLength == Catch::Approx(0.3));
    REQUIRE(value.holeSupportLength == Catch::Approx(0.5));
    REQUIRE(value.mean == Catch::Approx(
        0.5 * 2.0 * 0.3 + 0.5 * 1.0 * (0.2 + 0.3)));

    options.electronMinimumDrivingField = 6.0;
    REQUIRE(integrateIonizationPath(path, options).electron == 0.0);
    options.electronMinimumDrivingField = -1.0;
    REQUIRE_THROWS_AS(
        integrateIonizationPath(path, options), std::invalid_argument);
}

TEST_CASE("BreakAtIonIntegral ordering flattens electron and hole path values",
          "[impact][path-ionization][break-ordering]")
{
    PathIonizationAnalysis analysis;
    analysis.paths.resize(2);
    analysis.paths[0].integral = {1.2, 0.9, 1.05};
    analysis.paths[1].integral = {1.1, 1.05, 1.075};

    REQUIRE(nthLargestCarrierIonizationIntegral(analysis, 1) ==
            Catch::Approx(1.2));
    REQUIRE(nthLargestCarrierIonizationIntegral(analysis, 2) ==
            Catch::Approx(1.1));
    REQUIRE(nthLargestCarrierIonizationIntegral(analysis, 3) ==
            Catch::Approx(1.05));
    REQUIRE(nthLargestCarrierIonizationIntegral(analysis, 5) == 0.0);
}

TEST_CASE("Path tracer crosses a local electric-field maximum",
          "[impact][path-ionization]")
{
    VectorXd potential(4);
    potential << 0.0, 1.0, 2.0, 3.0;
    const std::vector<PathIonizationSegment> edges = {
        {10, 0, 1, 1.0, 0.2, 0.1, 1.0},
        {11, 1, 2, 1.0, 0.2, 0.1, 3.0},
        {12, 2, 3, 1.0, 0.2, 0.1, 2.0},
    };
    const PathIonizationAnalysis analysis = analyzeIonizationPaths(
        4, potential, edges);

    REQUIRE(analysis.paths.size() == 1);
    REQUIRE(analysis.paths.front().seedEdgeId == 11);
    REQUIRE(analysis.paths.front().nodes == std::vector<Index>{3, 2, 1, 0});
    REQUIRE(analysis.paths.front().segments.size() == 3);
    REQUIRE(nthLargestMeanIonizationIntegral(analysis, 1) ==
            Catch::Approx(analysis.paths.front().integral.mean));
    REQUIRE(nthLargestMeanIonizationIntegral(analysis, 2) == 0.0);
    for (Real value : analysis.meanNodeValues)
        REQUIRE(value == Catch::Approx(analysis.paths.front().integral.mean));
}

TEST_CASE("Path tracer stops when the electric field falls below the configured threshold",
          "[impact][path-ionization][stop-field]")
{
    VectorXd potential(4);
    potential << 0.0, 1.0, 2.0, 3.0;
    const std::vector<PathIonizationSegment> edges = {
        {10, 0, 1, 1.0, 0.2, 0.1, 1.0},
        {11, 1, 2, 1.0, 0.2, 0.1, 3.0},
        {12, 2, 3, 1.0, 0.2, 0.1, 2.0},
    };

    const PathIonizationAnalysis analysis = analyzeIonizationPaths(
        4, potential, edges, 0, 1.5);
    REQUIRE(analysis.paths.size() == 1);
    REQUIRE(analysis.paths.front().nodes == std::vector<Index>{3, 2, 1});
    REQUIRE(analysis.paths.front().segments.size() == 2);
    REQUIRE(analysis.meanNodeValues[0] == 0.0);
    REQUIRE(analysis.meanNodeValues[1] > 0.0);
    REQUIRE_THROWS_AS(
        analyzeIonizationPaths(4, potential, edges, 0, -1.0),
        std::invalid_argument);
}

TEST_CASE("Path tracer follows the two-dimensional electric-field direction",
          "[impact][path-ionization][field-line]")
{
    VectorXd potential(5);
    potential << 0.0, 1.0, 2.0, 4.0, 3.0;
    const std::vector<PathIonizationSegment> edges = {
        {10, 0, 1, 1.0, 0.2, 0.1, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
        {11, 1, 2, 1.0, 0.2, 0.1, 3.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
        // Larger potential slope, but perpendicular to the local field line.
        {12, 2, 3, 1.0, 0.2, 0.1, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0},
        // Smaller potential slope and aligned with -E toward higher potential.
        {13, 2, 4, 1.0, 0.2, 0.1, 2.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
    };

    const PathIonizationAnalysis analysis = analyzeIonizationPaths(
        5, potential, edges, 1);
    REQUIRE(analysis.paths.size() == 1);
    REQUIRE(analysis.paths.front().nodes == std::vector<Index>{4, 2, 1, 0});
    REQUIRE(analysis.paths.front().segments.front().edgeId == 13);
}

TEST_CASE("Continuous field-line tracer crosses triangle interiors",
          "[impact][path-ionization][continuous-cell]")
{
    std::vector<CellIonizationSample> cells = {
        {10,
         {0, 1, 2},
         {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{0.0, 1.0}},
         {0.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
        {11,
         {1, 3, 2},
         {Point2{1.0, 0.0}, Point2{1.0, 1.0}, Point2{0.0, 1.0}},
         {-1.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
    };
    for (CellIonizationSample& cell : cells) {
        cell.vertexElectricFieldsVPerM = {
            Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}};
    }

    const PathIonizationAnalysis analysis =
        analyzeContinuousCellIonizationPaths(4, cells, 0, 0.0, 0.0);
    REQUIRE(analysis.paths.size() == 1);
    const IonizationPath& path = analysis.paths.front();
    REQUIRE(path.seedCellId == 10);
    REQUIRE(path.seedNodeId == 0);
    REQUIRE(path.segments.size() == 3);
    Real length = 0.0;
    for (const PathIonizationSegment& segment : path.segments) {
        REQUIRE(segment.explicitGeometry);
        length += segment.length;
    }
    REQUIRE(length == Catch::Approx(1.0).epsilon(1.0e-9));
    const Real expectedElectron = 2.0 * (1.0 - std::exp(-1.0));
    const Real expectedHole = std::exp(1.0) - 1.0;
    REQUIRE(path.integral.electron == Catch::Approx(expectedElectron).epsilon(1.0e-9));
    REQUIRE(path.integral.hole == Catch::Approx(expectedHole).epsilon(1.0e-9));
}

TEST_CASE("Continuous tracer honors an explicit carrier-current tangent",
          "[impact][path-ionization][continuous-cell][tracing-vector]")
{
    std::vector<CellIonizationSample> cells = {
        {10,
         {0, 1, 2},
         {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{0.0, 1.0}},
         {0.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
        {11,
         {1, 3, 2},
         {Point2{1.0, 0.0}, Point2{1.0, 1.0}, Point2{0.0, 1.0}},
         {-1.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
    };
    for (CellIonizationSample& cell : cells) {
        cell.vertexElectricFieldsVPerM = {
            Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}};
        cell.vertexTracingDirections = {
            Point2{0.0, 1.0}, Point2{0.0, 1.0}, Point2{0.0, 1.0}};
    }

    const PathIonizationAnalysis analysis =
        analyzeContinuousCellIonizationPaths(4, cells, 1, 0.0, 0.0);
    REQUIRE(analysis.paths.size() == 1);
    Real maximumX = 0.0;
    Real maximumY = 0.0;
    for (const PathIonizationSegment& segment : analysis.paths.front().segments) {
        maximumX = std::max(
            maximumX, std::max(segment.startPointM.x(), segment.endPointM.x()));
        maximumY = std::max(
            maximumY, std::max(segment.startPointM.y(), segment.endPointM.y()));
    }
    REQUIRE(maximumX < 1.0e-6);
    REQUIRE(maximumY == Catch::Approx(1.0).epsilon(1.0e-8));

    const PathIonizationAnalysis cellElectric =
        analyzeContinuousCellIonizationPaths(
            4, cells, 1, 0.0, 0.0,
            ContinuousPathDirection::Bidirectional, {},
            ContinuousPathRetention::DistinctLocalMaxima,
            ContinuousPathSeedMode::NodalLocalMaxima,
            ContinuousPathTracingPolicy::CellElectricField);
    REQUIRE(cellElectric.paths.size() == 1);
    Real cellMaximumX = 0.0;
    Real cellMaximumY = 0.0;
    for (const PathIonizationSegment& segment :
         cellElectric.paths.front().segments) {
        cellMaximumX = std::max(
            cellMaximumX,
            std::max(segment.startPointM.x(), segment.endPointM.x()));
        cellMaximumY = std::max(
            cellMaximumY,
            std::max(segment.startPointM.y(), segment.endPointM.y()));
    }
    REQUIRE(cellMaximumX == Catch::Approx(1.0).epsilon(1.0e-8));
    REQUIRE(cellMaximumY < 1.0e-6);

    const PathIonizationAnalysis rk4 = analyzeContinuousCellIonizationPaths(
        4, cells, 1, 0.0, 0.0,
        ContinuousPathDirection::Bidirectional, {},
        ContinuousPathRetention::DistinctLocalMaxima,
        ContinuousPathSeedMode::NodalLocalMaxima,
        ContinuousPathTracingPolicy::P1ElectricFieldRungeKutta);
    REQUIRE(rk4.paths.size() == 1);
    Real rk4MaximumX = 0.0;
    Real rk4MaximumY = 0.0;
    for (const PathIonizationSegment& segment : rk4.paths.front().segments) {
        rk4MaximumX = std::max(
            rk4MaximumX,
            std::max(segment.startPointM.x(), segment.endPointM.x()));
        rk4MaximumY = std::max(
            rk4MaximumY,
            std::max(segment.startPointM.y(), segment.endPointM.y()));
    }
    REQUIRE(rk4MaximumX < 1.0e-6);
    REQUIRE(rk4MaximumY == Catch::Approx(1.0).epsilon(1.0e-6));
}

TEST_CASE("Adaptive Eparallel tracing uses minority QF direction at an n-type boundary peak",
          "[impact][path-ionization][continuous-cell][adaptive-tracing]")
{
    std::vector<CellIonizationSample> cells = {
        {10,
         {0, 1, 2},
         {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{0.0, 1.0}},
         {0.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
        {11,
         {1, 3, 2},
         {Point2{1.0, 0.0}, Point2{1.0, 1.0}, Point2{0.0, 1.0}},
         {-1.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
    };
    for (CellIonizationSample& cell : cells) {
        cell.vertexElectricFieldsVPerM = {
            Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}};
        cell.vertexTracingDirections = {
            Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}};
        cell.vertexElectronCurrentDirections = cell.vertexTracingDirections;
        cell.vertexHoleQfDirections = {
            Point2{0.0, 1.0}, Point2{0.0, 1.0}, Point2{0.0, 1.0}};
        cell.vertexNetDoping = {1.0, 1.0, 1.0};
        cell.vertexTransportBoundary = {true, true, true};
    }

    const PathIonizationAnalysis analysis =
        analyzeContinuousCellIonizationPaths(
            4, cells, 1, 0.0, 0.0,
            ContinuousPathDirection::Bidirectional, {},
            ContinuousPathRetention::DistinctLocalMaxima,
            ContinuousPathSeedMode::NodalLocalMaxima,
            ContinuousPathTracingPolicy::SentaurusEparallelAdaptive);

    REQUIRE(analysis.paths.size() == 1);
    Real maximumX = 0.0;
    Real maximumY = 0.0;
    for (const PathIonizationSegment& segment : analysis.paths.front().segments) {
        maximumX = std::max(
            maximumX, std::max(segment.startPointM.x(), segment.endPointM.x()));
        maximumY = std::max(
            maximumY, std::max(segment.startPointM.y(), segment.endPointM.y()));
    }
    REQUIRE(maximumX < 1.0e-6);
    REQUIRE(maximumY == Catch::Approx(1.0).epsilon(1.0e-8));
}

TEST_CASE("Continuous tracer clusters two-ring nodal maxima into one rank",
          "[impact][path-ionization][continuous-cell][dedup]")
{
    std::vector<CellIonizationSample> cells = {
        {10,
         {0, 1, 2},
         {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{0.0, 1.0}},
         {0.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
        {11,
         {1, 3, 2},
         {Point2{1.0, 0.0}, Point2{1.0, 1.0}, Point2{0.0, 1.0}},
         {-1.0, -1.0, 0.0},
         Point2{1.0, 0.0},
         1.0, 1.0, 2.0, 1.0},
    };
    cells[0].vertexElectricFieldsVPerM = {
        Point2{3.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}};
    cells[1].vertexElectricFieldsVPerM = {
        Point2{1.0, 0.0}, Point2{3.0, 0.0}, Point2{1.0, 0.0}};

    const PathIonizationAnalysis all =
        analyzeContinuousCellIonizationPaths(
            4, cells, 0, 0.0, 0.0,
            ContinuousPathDirection::Bidirectional, {},
            ContinuousPathRetention::AllSeedTrajectories,
            ContinuousPathSeedMode::NodalLocalMaxima);
    const PathIonizationAnalysis distinct =
        analyzeContinuousCellIonizationPaths(
            4, cells, 0, 0.0, 0.0,
            ContinuousPathDirection::Bidirectional, {},
            ContinuousPathRetention::DistinctLocalMaxima,
            ContinuousPathSeedMode::NodalLocalMaxima);
    REQUIRE(all.paths.size() >= distinct.paths.size());
    REQUIRE(distinct.paths.size() == 2);
    const auto globalPeak = std::find_if(
        distinct.paths.begin(), distinct.paths.end(),
        [](const IonizationPath& path) { return path.seedNodeId == 0; });
    const auto secondaryPeak = std::find_if(
        distinct.paths.begin(), distinct.paths.end(),
        [](const IonizationPath& path) { return path.seedNodeId == 3; });
    REQUIRE(globalPeak != distinct.paths.end());
    REQUIRE(secondaryPeak != distinct.paths.end());
    REQUIRE(globalPeak->seedField == Catch::Approx(3.0));
    REQUIRE(globalPeak->saddleField == Catch::Approx(0.0));
    REQUIRE(globalPeak->peakProminenceRatio == Catch::Approx(1.0));
    REQUIRE(secondaryPeak->seedField == Catch::Approx(3.0));
    REQUIRE(secondaryPeak->saddleField == Catch::Approx(1.0));
    REQUIRE(secondaryPeak->peakProminenceRatio == Catch::Approx(2.0 / 3.0));
    REQUIRE(secondaryPeak->parentPeakNodeId == 0);

    const PathIonizationAnalysis numbered =
        analyzeContinuousCellIonizationPaths(
            4, cells, 0, 0.0, 0.0,
            ContinuousPathDirection::Bidirectional, {},
            ContinuousPathRetention::NumberedPeakGroups,
            ContinuousPathSeedMode::NodalLocalMaxima);
    REQUIRE(numbered.paths.size() == 2);
    REQUIRE(numbered.paths[0].seedNodeId != numbered.paths[1].seedNodeId);
    REQUIRE(numbered.paths[0].integral.mean ==
            Catch::Approx(numbered.paths[1].integral.mean));
    REQUIRE(numbered.paths[0].physicalPathGroupId ==
            numbered.paths[1].physicalPathGroupId);
    REQUIRE(nthLargestMeanIonizationIntegral(numbered, 1) ==
            Catch::Approx(numbered.paths[0].integral.mean));
    REQUIRE(nthLargestMeanIonizationIntegral(numbered, 2) == 0.0);

    const PathIonizationAnalysis analysis =
        analyzeContinuousCellIonizationPaths(
            4, cells, 0, 0.0, 0.0,
            ContinuousPathDirection::Bidirectional, {},
            ContinuousPathRetention::CorridorDeduplicated,
            ContinuousPathSeedMode::NodalLocalMaxima);
    REQUIRE(analysis.paths.size() == 1);
    REQUIRE((analysis.paths.front().seedNodeId == 0 ||
             analysis.paths.front().seedNodeId == 3));
}

TEST_CASE("Path ionization input validation rejects negative coefficients",
          "[impact][path-ionization]")
{
    REQUIRE_THROWS_AS(
        integrateIonizationPath({{0, 0, 1, 1.0, -1.0, 1.0, 1.0}}),
        std::invalid_argument);
}
