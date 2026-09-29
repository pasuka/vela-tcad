#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "vela/simulation/DCSweep.h"
#include "vela/simulation/DCSweepStepControl.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace vela;

namespace {

struct ScopedDirectoryCleanup {
    std::filesystem::path path;
    ~ScopedDirectoryCleanup()
    {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
};

std::filesystem::path makeUniqueDir()
{
    const auto base = std::filesystem::temp_directory_path();
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto tid = std::hash<std::thread::id>{}(std::this_thread::get_id());
    std::mt19937_64 rng(static_cast<std::mt19937_64::result_type>(stamp ^ tid));
    std::uniform_int_distribution<unsigned long long> dist;
    for (int attempt = 0; attempt < 8; ++attempt) {
        const auto dir = base / ("vela_dc_output_contract_" + std::to_string(stamp) + "_" + std::to_string(dist(rng)));
        if (std::filesystem::create_directory(dir))
            return dir;
    }
    throw std::runtime_error("Failed to create a unique temp directory.");
}

std::filesystem::path writePNMeshMicrometers(const std::filesystem::path& dir)
{
    nlohmann::json mesh = {
        {"nodes", {
            {{"id", 0}, {"x", 0.0}, {"y", 0.0}},
            {{"id", 1}, {"x", 1.0}, {"y", 0.0}},
            {{"id", 2}, {"x", 1.0}, {"y", 1.0}},
            {{"id", 3}, {"x", 0.0}, {"y", 1.0}}
        }},
        {"triangles", {
            {{"id", 0}, {"region_id", 0}, {"node_ids", {0, 1, 2}}},
            {{"id", 1}, {"region_id", 1}, {"node_ids", {0, 2, 3}}}
        }},
        {"regions", {
            {{"id", 0}, {"name", "n_region"}, {"material", "Si"}, {"cell_ids", {0}}},
            {{"id", 1}, {"name", "p_region"}, {"material", "Si"}, {"cell_ids", {1}}}
        }},
        {"contacts", {
            {{"id", 0}, {"name", "anode"}, {"region_id", 1}, {"node_ids", {0, 3}}},
            {{"id", 1}, {"name", "cathode"}, {"region_id", 0}, {"node_ids", {1, 2}}}
        }}
    };
    const auto meshPath = dir / "pn_mesh_um.json";
    std::ofstream(meshPath) << mesh.dump(2);
    return meshPath;
}

nlohmann::json baseConfig(const std::filesystem::path& dir, const std::filesystem::path& meshPath,
                          const std::filesystem::path& csvPath)
{
    return {
        {"mesh_file", meshPath.string()},
        {"output_csv", csvPath.string()},
        {"scaling", {{"mode", "unit_scaling"}}},
        {"doping", {
            {{"region", "n_region"}, {"donors", 1.0e17}, {"acceptors", 0.0}},
            {{"region", "p_region"}, {"donors", 0.0}, {"acceptors", 1.0e17}}
        }},
        {"contacts", {
            {{"name", "anode"}, {"bias", 0.0}},
            {{"name", "cathode"}, {"bias", 0.0}}
        }},
        {"solver", {
            {"max_iter", 80},
            {"reltol", 1.0e-5},
            {"damping_psi", 0.5}
        }},
        {"sweep", {
            {"contact", "anode"},
            {"start", 0.0},
            {"stop", 0.5},
            {"step", 0.25},
            {"current_contact", "anode"},
            {"write_vtk", false},
            {"vtk_prefix", (dir / "pn_out").string()}
        }}
    };
}

std::filesystem::path writeConfig(const std::filesystem::path& dir, const nlohmann::json& cfg,
                                  const std::string& name)
{
    const auto cfgPath = dir / name;
    std::ofstream(cfgPath) << cfg.dump(2);
    return cfgPath;
}

std::vector<std::vector<std::string>> readCsvRows(const std::filesystem::path& csvPath)
{
    std::vector<std::vector<std::string>> rows;
    std::ifstream in(csvPath);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        std::vector<std::string> row;
        std::string cell;
        std::istringstream stream(line);
        while (std::getline(stream, cell, ','))
            row.push_back(cell);
        rows.push_back(row);
    }
    return rows;
}

std::size_t column(const std::vector<std::string>& header, const std::string& name)
{
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name)
            return i;
    }
    throw std::runtime_error("missing CSV column " + name);
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

std::vector<std::string> pointDataNames(const std::string& vtu)
{
    const std::size_t start = vtu.find("<PointData>");
    const std::size_t end = vtu.find("</PointData>");
    REQUIRE(start != std::string::npos);
    REQUIRE(end != std::string::npos);
    const std::string section = vtu.substr(start, end - start);
    static const std::regex nameRe(R"(<DataArray[^>]*Name=\"([^\"]+)\")");
    std::vector<std::string> names;
    for (auto it = std::sregex_iterator(section.begin(), section.end(), nameRe); it != std::sregex_iterator(); ++it)
        names.push_back((*it)[1].str());
    return names;
}

std::vector<Real> visitedVoltages(const detail::DCSweepStepControlConfig& cfg,
                                  const std::function<bool(Real)>& converges)
{
    std::vector<Real> visited;
    detail::runDCSweepStepControl(
        cfg,
        [&](Real voltage, Real, int) -> detail::DCSweepStepAttemptResult { return converges(voltage); },
        [&](const detail::DCSweepStepControlEvent& event) {
            if (event.converged)
                visited.push_back(event.voltage);
        });
    return visited;
}

} // namespace

TEST_CASE("Step control lands exactly on mandatory points without changing the other targets",
          "[dc_sweep][output_contract][step_control]")
{
    detail::DCSweepStepControlConfig cfg;
    cfg.start = 0.0;
    cfg.stop = 1.0;
    cfg.step = 0.25;
    cfg.initialStep = 0.25;
    cfg.minStep = 0.01;
    cfg.maxStep = 0.25;
    cfg.growthFactor = 1.0;
    cfg.maxRetries = 3;
    const auto always = [](Real) { return true; };
    const std::vector<Real> baseline = visitedVoltages(cfg, always);
    REQUIRE(baseline.size() == 4);

    cfg.mandatoryPoints = {0.3, 0.7};
    const std::vector<Real> withPoints = visitedVoltages(cfg, always);
    REQUIRE(withPoints.size() == baseline.size() + 2);
    // The mandatory points are hit to 1e-12 V ...
    REQUIRE(std::abs(withPoints[1] - 0.3) <= 1.0e-12);
    REQUIRE(std::abs(withPoints[3] - 0.7) <= 1.0e-12);
    // ... and every other visited voltage is unchanged.
    std::vector<Real> others;
    for (std::size_t i = 0; i < withPoints.size(); ++i) {
        if (i != 1 && i != 3)
            others.push_back(withPoints[i]);
    }
    REQUIRE(others.size() == baseline.size());
    for (std::size_t i = 0; i < baseline.size(); ++i)
        REQUIRE(others[i] == Catch::Approx(baseline[i]).margin(1.0e-12));

    // Points outside the interval or duplicating a nominal target are ignored.
    cfg.mandatoryPoints = {-0.5, 0.5, 1.5, 0.75};
    REQUIRE(visitedVoltages(cfg, always) == baseline);

    // A mandatory point that needs several adaptive steps is still reached exactly.
    detail::DCSweepStepControlConfig fine = cfg;
    fine.mandatoryPoints = {0.3};
    fine.initialStep = 0.08;
    fine.minStep = 0.01;
    fine.maxStep = 0.08;
    const std::vector<Real> fineVisited = visitedVoltages(fine, always);
    REQUIRE(std::any_of(fineVisited.begin(), fineVisited.end(),
                        [](Real v) { return std::abs(v - 0.3) <= 1.0e-12; }));
    REQUIRE(std::abs(fineVisited.back() - 1.0) <= 1.0e-12);
}

TEST_CASE("DCSweep writes TDR-named .vtu files at mandatory bias points with the requested subset",
          "[dc_sweep][output_contract][vtu]")
{
    const auto dir = makeUniqueDir();
    const ScopedDirectoryCleanup cleanup{dir};
    const auto meshPath = writePNMeshMicrometers(dir);

    // Reference run: legacy stepping without mandatory points.
    const auto referenceCsv = dir / "reference.csv";
    nlohmann::json reference = baseConfig(dir, meshPath, referenceCsv);
    const auto referenceCfg = writeConfig(dir, reference, "reference.json");
    DCSweep referenceSweep;
    const DCSweepResult referenceResult = referenceSweep.runWithResult(referenceCfg.string());
    REQUIRE(referenceResult.points.size() == 3);

    const auto csvPath = dir / "vtu.csv";
    nlohmann::json cfg = baseConfig(dir, meshPath, csvPath);
    cfg["sweep"]["vtk_format"] = "vtu";
    cfg["sweep"]["vtk_bias_points"] = {0.1};
    cfg["sweep"]["vtk_fields"] = {"ElectrostaticPotential", "eDensity", "hDensity", "DopingConcentration",
                                  "ElectricField/Vector", "eCurrentDensity", "eMobility",
                                  "srhRecombination", "eTemperature"};
    cfg["sweep"]["csv_all_contacts"] = true;
    const auto cfgPath = writeConfig(dir, cfg, "vtu.json");

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    const auto& points = result.points;
    REQUIRE(points.size() == 4);
    REQUIRE(std::abs(points[1].bias - 0.1) <= 1.0e-12);
    for (std::size_t i = 0, j = 0; i < points.size(); ++i) {
        if (i == 1)
            continue;
        REQUIRE(points[i].bias == Catch::Approx(referenceResult.points[j].bias).margin(1.0e-12));
        REQUIRE(points[i].totalCurrent == Catch::Approx(referenceResult.points[j].totalCurrent).epsilon(1.0e-9));
        ++j;
    }

    // Only the mandatory point produced a file, and it is XML with TDR names.
    REQUIRE(points[0].outputVtk.empty());
    REQUIRE_FALSE(points[1].outputVtk.empty());
    REQUIRE(points[2].outputVtk.empty());
    REQUIRE(std::filesystem::path(points[1].outputVtk).extension() == ".vtu");
    REQUIRE(std::filesystem::exists(points[1].outputVtk));
    const std::string vtu = readText(points[1].outputVtk);
    REQUIRE(vtu.rfind("<?xml version=\"1.0\"?>", 0) == 0);
    REQUIRE(pointDataNames(vtu) == std::vector<std::string>{
        "ElectrostaticPotential", "eDensity", "hDensity", "DopingConcentration",
        "ElectricField", "eCurrentDensity", "eMobility", "srhRecombination"});
    REQUIRE(vtu.find("Name=\"ElectricField\" NumberOfComponents=\"3\"") != std::string::npos);
    REQUIRE(vtu.find("Name=\"eCurrentDensity\" NumberOfComponents") == std::string::npos);
    REQUIRE(vtu.find("Name=\"Unit:eDensity\"") != std::string::npos);
    REQUIRE(vtu.find("Name=\"Bias_V\"") != std::string::npos);
    REQUIRE(vtu.find("Name=\"anode_TotalCurrent_A_per_um\"") != std::string::npos);
    REQUIRE(vtu.find("Name=\"region_id\"") != std::string::npos);
    REQUIRE(vtu.find("Name=\"Potential\"") == std::string::npos);  // no legacy names
    REQUIRE(vtu.find("Name=\"Electrons\"") == std::string::npos);
    REQUIRE(points[1].skippedVtkFields == std::vector<std::string>{"eTemperature"});

    // csv_all_contacts appends per-contact columns consistent with the monitored current.
    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == points.size() + 1);
    const auto& header = rows.front();
    const std::size_t anodeTotal = column(header, "contact_anode_total_A_per_um");
    const std::size_t cathodeTotal = column(header, "contact_cathode_total_A_per_um");
    const std::size_t anodeElectron = column(header, "contact_anode_electron_A_per_um");
    const std::size_t anodeHole = column(header, "contact_anode_hole_A_per_um");
    const std::size_t monitored = column(header, "current_total_A_per_um");
    for (std::size_t r = 1; r < rows.size(); ++r) {
        const Real anode = std::stod(rows[r][anodeTotal]);
        const Real cathode = std::stod(rows[r][cathodeTotal]);
        REQUIRE(anode == Catch::Approx(std::stod(rows[r][monitored])).margin(1.0e-300));
        // ContactCurrent convention: total = electron - hole (see src/post/ContactCurrent.cpp).
        REQUIRE(std::stod(rows[r][anodeElectron]) - std::stod(rows[r][anodeHole]) ==
                Catch::Approx(anode).epsilon(1.0e-9).margin(1.0e-30));
        REQUIRE(std::isfinite(cathode));
        // Two-terminal device: the currents balance to solver precision.
        REQUIRE(std::abs(anode + cathode) <= 1.0e-6 * std::max(std::abs(anode), 1.0e-30) + 1.0e-20);
    }
}

TEST_CASE("DCSweep vtk_format legacy keeps the ASCII .vtk output and write_vtk semantics",
          "[dc_sweep][output_contract][vtu]")
{
    const auto dir = makeUniqueDir();
    const ScopedDirectoryCleanup cleanup{dir};
    const auto meshPath = writePNMeshMicrometers(dir);
    const auto csvPath = dir / "legacy.csv";
    nlohmann::json cfg = baseConfig(dir, meshPath, csvPath);
    cfg["sweep"]["write_vtk"] = true;
    cfg["sweep"]["stop"] = 0.25;
    const auto cfgPath = writeConfig(dir, cfg, "legacy.json");

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 2);
    for (const DCSweepPoint& point : result.points) {
        REQUIRE(std::filesystem::path(point.outputVtk).extension() == ".vtk");
        const std::string text = readText(point.outputVtk);
        REQUIRE(text.rfind("# vtk DataFile Version 3.0", 0) == 0);
        REQUIRE(text.find("SCALARS Potential double 1") != std::string::npos);
    }

    nlohmann::json bad = baseConfig(dir, meshPath, csvPath);
    bad["sweep"]["vtk_format"] = "xml";
    const auto badPath = writeConfig(dir, bad, "bad.json");
    DCSweep badSweep;
    REQUIRE_THROWS_AS(badSweep.runWithResult(badPath.string()), std::invalid_argument);
}

TEST_CASE("DCSweep terminates on terminal-current criteria and marks the CSV row",
          "[dc_sweep][output_contract][break_criteria]")
{
    const auto dir = makeUniqueDir();
    const ScopedDirectoryCleanup cleanup{dir};
    const auto meshPath = writePNMeshMicrometers(dir);

    // Baseline forward sweep to learn the monotone current sequence.
    const auto baselineCsv = dir / "baseline.csv";
    nlohmann::json baseline = baseConfig(dir, meshPath, baselineCsv);
    baseline["sweep"]["stop"] = 0.75;
    const auto baselinePath = writeConfig(dir, baseline, "baseline.json");
    DCSweep baselineSweep;
    const DCSweepResult baselineResult = baselineSweep.runWithResult(baselinePath.string());
    REQUIRE(baselineResult.points.size() == 4);
    const auto baselineRows = readCsvRows(baselineCsv);
    REQUIRE(baselineRows.size() == 5);
    const std::size_t baselineCurrent = column(baselineRows.front(), "current_total_A_per_um");
    std::vector<Real> currents_A_per_um;
    for (std::size_t r = 1; r < baselineRows.size(); ++r) {
        REQUIRE(baselineResult.points[r - 1].converged);
        currents_A_per_um.push_back(std::stod(baselineRows[r][baselineCurrent]));
    }
    for (std::size_t i = 1; i < currents_A_per_um.size(); ++i)
        REQUIRE(std::abs(currents_A_per_um[i]) > std::abs(currents_A_per_um[i - 1]));

    // Threshold between the 0.25 V and 0.5 V currents: the sweep must stop at 0.5 V.
    const Real threshold = 0.5 * (std::abs(currents_A_per_um[1]) + std::abs(currents_A_per_um[2]));
    const auto csvPath = dir / "break.csv";
    nlohmann::json cfg = baseConfig(dir, meshPath, csvPath);
    cfg["sweep"]["stop"] = 0.75;
    cfg["sweep"]["breakdown"] = {{"current_abs_A_per_um", threshold}};
    const auto cfgPath = writeConfig(dir, cfg, "break.json");
    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());
    REQUIRE(result.points.size() == 3);
    REQUIRE(result.points.back().bias == Catch::Approx(0.5).margin(1.0e-12));
    REQUIRE(result.points.back().breakdownDetected);
    REQUIRE(result.points.back().breakdownCriterion == "current_abs");
    REQUIRE_FALSE(result.points[1].breakdownDetected);

    const auto rows = readCsvRows(csvPath);
    REQUIRE(rows.size() == 4);
    const std::size_t criterion = column(rows.front(), "stop_criterion");
    const std::size_t stopCurrent = column(rows.front(), "stop_current_A_per_um");
    REQUIRE(rows[1][criterion].empty());
    REQUIRE(rows[2][criterion].empty());
    REQUIRE(rows[3][criterion] == "current_abs");
    REQUIRE(std::stod(rows[3][stopCurrent]) == Catch::Approx(currents_A_per_um[2]).epsilon(1.0e-9));

    // Signed min/max criteria on the same sequence (forward current is negative
    // at the anode in this fixture or positive; use the sign of the data).
    const Real signedTarget = currents_A_per_um[2];
    nlohmann::json signedCfg = baseConfig(dir, meshPath, dir / "signed.csv");
    signedCfg["sweep"]["stop"] = 0.75;
    if (signedTarget < 0.0)
        signedCfg["sweep"]["breakdown"] = {{"current_min_A_per_um", 0.5 * (currents_A_per_um[1] + currents_A_per_um[2])}};
    else
        signedCfg["sweep"]["breakdown"] = {{"current_max_A_per_um", 0.5 * (currents_A_per_um[1] + currents_A_per_um[2])}};
    const auto signedPath = writeConfig(dir, signedCfg, "signed.json");
    DCSweep signedSweep;
    const DCSweepResult signedResult = signedSweep.runWithResult(signedPath.string());
    REQUIRE(signedResult.points.size() == 3);
    REQUIRE(signedResult.points.back().breakdownCriterion == (signedTarget < 0.0 ? "current_min" : "current_max"));

    // Invalid criteria fail closed.
    nlohmann::json bad = baseConfig(dir, meshPath, dir / "bad.csv");
    bad["sweep"]["breakdown"] = {{"current_abs_A_per_um", -1.0}};
    const auto badPath = writeConfig(dir, bad, "bad_criteria.json");
    DCSweep badSweep;
    REQUIRE_THROWS_AS(badSweep.runWithResult(badPath.string()), std::invalid_argument);
}
