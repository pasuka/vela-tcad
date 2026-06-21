#include "vela/io/DDSolutionCsv.h"

#include "vela/io/CsvUtils.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace vela {
namespace {

std::string formatRestartReal(Real value)
{
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

Real parseRestartStateReal(const std::string& text,
                           const std::string& column,
                           Index nodeId)
{
    std::size_t consumed = 0;
    Real value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "DCSweep: initial_state_file has invalid " + column +
            " '" + text + "' for node id " + std::to_string(nodeId));
    }
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::runtime_error(
            "DCSweep: initial_state_file has invalid " + column +
            " '" + text + "' for node id " + std::to_string(nodeId));
    }
    return value;
}

long long parseRestartStateNodeId(const std::string& nodeIdText)
{
    std::size_t consumed = 0;
    long long parsedNodeId = 0;
    try {
        parsedNodeId = std::stoll(nodeIdText, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "DCSweep: initial_state_file has invalid node id '" + nodeIdText + "'");
    }
    if (consumed != nodeIdText.size()) {
        throw std::runtime_error(
            "DCSweep: initial_state_file has invalid node id '" + nodeIdText + "'");
    }
    return parsedNodeId;
}

} // namespace

DDSolution readDDSolutionStateCsv(const std::filesystem::path& path,
                                  Index expectedNodeCount)
{
    std::ifstream input(path);
    if (!input.is_open())
        throw std::runtime_error("DCSweep: cannot open initial_state_file: " + path.string());

    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("DCSweep: initial_state_file is empty: " + path.string());

    const std::vector<std::string> expectedHeader = {
        "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"};
    const std::vector<std::string> header = splitCsvLine(
        line,
        "DCSweep: initial_state_file does not support quoted fields.");
    if (header != expectedHeader)
        throw std::runtime_error(
            "DCSweep: initial_state_file header must be "
            "node_id,psi,phin,phip,electrons_m3,holes_m3");

    DDSolution solution;
    solution.psi = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.phin = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.phip = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.n = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.p = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.iters = 0;
    solution.converged = true;

    std::vector<bool> seen(expectedNodeCount, false);
    while (std::getline(input, line)) {
        if (trimCsvToken(line).empty())
            continue;
        const std::vector<std::string> row = splitCsvLine(
            line,
            "DCSweep: initial_state_file does not support quoted fields.");
        if (row.size() != expectedHeader.size())
            throw std::runtime_error(
                "DCSweep: initial_state_file rows must have 6 columns.");
        const long long parsedNodeId = parseRestartStateNodeId(row.at(0));
        if (parsedNodeId < 0 ||
            parsedNodeId >= static_cast<long long>(expectedNodeCount)) {
            throw std::runtime_error(
                "DCSweep: initial_state_file has out-of-range node id " +
                std::to_string(parsedNodeId));
        }
        const Index nodeId = static_cast<Index>(parsedNodeId);
        if (seen.at(nodeId)) {
            throw std::runtime_error(
                "DCSweep: initial_state_file has duplicate row for node id " +
                std::to_string(nodeId));
        }
        seen.at(nodeId) = true;
        const int rowIndex = static_cast<int>(nodeId);
        solution.psi(rowIndex) = parseRestartStateReal(row.at(1), "psi", nodeId);
        solution.phin(rowIndex) = parseRestartStateReal(row.at(2), "phin", nodeId);
        solution.phip(rowIndex) = parseRestartStateReal(row.at(3), "phip", nodeId);
        solution.n(rowIndex) = parseRestartStateReal(row.at(4), "electrons_m3", nodeId);
        solution.p(rowIndex) = parseRestartStateReal(row.at(5), "holes_m3", nodeId);
    }

    for (Index nodeId = 0; nodeId < expectedNodeCount; ++nodeId) {
        if (!seen.at(nodeId)) {
            throw std::runtime_error(
                "DCSweep: initial_state_file missing row for node id " +
                std::to_string(nodeId));
        }
    }
    return solution;
}

void writeDDSolutionStateCsv(const std::filesystem::path& path,
                             const DDSolution& solution)
{
    const auto fieldSize = solution.psi.size();
    if (solution.phin.size() != fieldSize ||
        solution.phip.size() != fieldSize ||
        solution.n.size() != fieldSize ||
        solution.p.size() != fieldSize) {
        throw std::runtime_error("DCSweep: cannot write restart state with inconsistent field sizes.");
    }

    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    if (!output.is_open())
        throw std::runtime_error("DCSweep: cannot open write_state_file: " + path.string());

    output << "node_id,psi,phin,phip,electrons_m3,holes_m3\n";
    for (int i = 0; i < fieldSize; ++i) {
        output << i << ','
               << formatRestartReal(solution.psi(i)) << ','
               << formatRestartReal(solution.phin(i)) << ','
               << formatRestartReal(solution.phip(i)) << ','
               << formatRestartReal(solution.n(i)) << ','
               << formatRestartReal(solution.p(i)) << '\n';
    }
}

} // namespace vela
