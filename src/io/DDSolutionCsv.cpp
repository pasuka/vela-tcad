#include "vela/io/DDSolutionCsv.h"

#include "vela/io/CsvUtils.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace vela {
namespace {

std::string formatRestartReal(Real value)
{
    if (!std::isfinite(value)) {
        throw std::runtime_error(
            "DCSweep: cannot write non-finite restart-state value.");
    }
    if (value != 0.0 &&
        std::abs(value) < std::numeric_limits<Real>::min()) {
        return "0";
    }
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

Real parseSplitTail(const std::string& text) {
    Real value=0.;auto result=std::from_chars(text.data(),text.data()+text.size(),value);
    if(result.ec!=std::errc{} || result.ptr!=text.data()+text.size() || !std::isfinite(value))throw std::runtime_error("Invalid exact split coordinate tail");
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

std::optional<std::size_t> columnIndex(
    const std::vector<std::string>& header,
    const std::string& name)
{
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name)
            return i;
    }
    return std::nullopt;
}

} // namespace

DDSolution readDDSolutionStateCsv(const std::filesystem::path& path,
                                  Index expectedNodeCount,
                                  UnitScalingConfig scaling)
{
    std::ifstream input(path);
    if (!input.is_open())
        throw std::runtime_error("DCSweep: cannot open initial_state_file: " + path.string());

    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("DCSweep: initial_state_file is empty: " + path.string());

    const std::vector<std::string> requiredHeader = {
        "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"};
    const std::vector<std::string> header = splitCsvLine(
        line,
        "DCSweep: initial_state_file does not support quoted fields.");
    if (header.size() < requiredHeader.size() ||
        !std::equal(requiredHeader.begin(), requiredHeader.end(), header.begin())) {
        throw std::runtime_error(
            "DCSweep: initial_state_file header must be "
            "node_id,psi,phin,phip,electrons_m3,holes_m3 with optional "
            "restart-state columns");
    }

    const std::vector<std::string> optionalColumns = {
        "electron_quantum_potential_V",
        "electron_quantum_potential_like_V",
        "electron_qf_increment_V",
        "hole_qf_increment_V",
        "electron_qf_reference_V",
        "hole_qf_reference_V",
        "packed_psi", "packed_electron_qf_increment", "packed_hole_qf_increment",
        "packed_potential_scale_V", "packed_psi_low", "packed_electron_qf_increment_low", "packed_hole_qf_increment_low", "split_state_schema", "split_mesh_fingerprint"};
    for (std::size_t i = requiredHeader.size(); i < header.size(); ++i) {
        if (std::find(optionalColumns.begin(), optionalColumns.end(), header[i]) ==
            optionalColumns.end()) {
            throw std::runtime_error(
                "DCSweep: initial_state_file has unsupported column '" +
                header[i] + "'");
        }
        if (std::count(header.begin(), header.end(), header[i]) != 1) {
            throw std::runtime_error(
                "DCSweep: initial_state_file has duplicate column '" +
                header[i] + "'");
        }
    }
    const auto quantumPotentialColumn =
        columnIndex(header, "electron_quantum_potential_V");
    const auto quantumPotentialLikeColumn =
        columnIndex(header, "electron_quantum_potential_like_V");
    const auto electronIncrementColumn =
        columnIndex(header, "electron_qf_increment_V");
    const auto holeIncrementColumn =
        columnIndex(header, "hole_qf_increment_V");
    const auto electronReferenceColumn =
        columnIndex(header, "electron_qf_reference_V");
    const auto holeReferenceColumn =
        columnIndex(header, "hole_qf_reference_V");
    const int qfCoordinateColumnCount =
        static_cast<int>(electronIncrementColumn.has_value()) +
        static_cast<int>(holeIncrementColumn.has_value()) +
        static_cast<int>(electronReferenceColumn.has_value()) +
        static_cast<int>(holeReferenceColumn.has_value());
    if (qfCoordinateColumnCount != 0 && qfCoordinateColumnCount != 4) {
        throw std::runtime_error(
            "DCSweep: initial_state_file must provide all four quasi-Fermi "
            "increment/reference columns together.");
    }
    if (quantumPotentialLikeColumn && !quantumPotentialColumn) {
        throw std::runtime_error(
            "DCSweep: electron_quantum_potential_like_V requires "
            "electron_quantum_potential_V.");
    }
    const bool hasReferencedQf = qfCoordinateColumnCount == 4;
    const auto packedPsiColumn=columnIndex(header,"packed_psi");
    const auto packedElectronColumn=columnIndex(header,"packed_electron_qf_increment");
    const auto packedHoleColumn=columnIndex(header,"packed_hole_qf_increment");
    const auto packedScaleColumn=columnIndex(header,"packed_potential_scale_V");
    const int packedColumnCount=int(packedPsiColumn.has_value())+int(packedElectronColumn.has_value())+
        int(packedHoleColumn.has_value())+int(packedScaleColumn.has_value());
    if(packedColumnCount!=0 && (packedColumnCount!=4 || !hasReferencedQf))
        throw std::runtime_error("DCSweep: packed restart requires all four packed columns and quasi-Fermi references/increments.");
    const bool hasPacked=packedColumnCount==4;
    const auto lowPsi=columnIndex(header,"packed_psi_low"),lowN=columnIndex(header,"packed_electron_qf_increment_low"),lowP=columnIndex(header,"packed_hole_qf_increment_low");
    const auto splitSchema=columnIndex(header,"split_state_schema"),splitMesh=columnIndex(header,"split_mesh_fingerprint");
    const int splitCount=int(lowPsi.has_value())+int(lowN.has_value())+int(lowP.has_value())+int(splitSchema.has_value())+int(splitMesh.has_value());
    if(splitCount!=0 && (splitCount!=5 || !hasPacked))throw std::runtime_error("Split restart requires all tails, schema, mesh and packed coordinates");
    const bool hasSplit=splitCount==5;


    DDSolution solution;
    solution.psi = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.phin = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    const PhysicalUnitSystem& units = scaling.unitSystem();

    solution.phip = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.n = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.p = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    solution.electronQuantumPotential =
        VectorXd::Zero(static_cast<int>(expectedNodeCount));
    if (quantumPotentialLikeColumn) {
        solution.electronQuantumPotentialLike =
            VectorXd::Zero(static_cast<int>(expectedNodeCount));
    }
    if (hasReferencedQf) {
        solution.phinIncrement = VectorXd::Zero(static_cast<int>(expectedNodeCount));
        solution.phipIncrement = VectorXd::Zero(static_cast<int>(expectedNodeCount));
        solution.electronQfReference = VectorXd::Zero(static_cast<int>(expectedNodeCount));
        solution.holeQfReference = VectorXd::Zero(static_cast<int>(expectedNodeCount));
    }
    solution.iters = 0;
    solution.converged = true;
    if(hasPacked)solution.packedState=VectorXd::Zero(3*static_cast<int>(expectedNodeCount));
    if(hasSplit)solution.packedLow=VectorXd::Zero(3*static_cast<int>(expectedNodeCount));

    std::vector<bool> seen(expectedNodeCount, false);
    while (std::getline(input, line)) {
        if (trimCsvToken(line).empty())
            continue;
        const std::vector<std::string> row = splitCsvLine(
            line,
            "DCSweep: initial_state_file does not support quoted fields.");
        if (row.size() != header.size())
            throw std::runtime_error(
                "DCSweep: initial_state_file rows do not match the header.");
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
        solution.n(rowIndex) = units.m3ToInternalConcentration(
            parseRestartStateReal(row.at(4), "electrons_m3", nodeId));
        solution.p(rowIndex) = units.m3ToInternalConcentration(
            parseRestartStateReal(row.at(5), "holes_m3", nodeId));
        if(hasPacked) {
            const int count=static_cast<int>(expectedNodeCount);
            solution.packedState(rowIndex)=parseRestartStateReal(row.at(*packedPsiColumn),"packed_psi",nodeId);
            solution.packedState(count+rowIndex)=parseRestartStateReal(row.at(*packedElectronColumn),"packed_electron_qf_increment",nodeId);
            solution.packedState(2*count+rowIndex)=parseRestartStateReal(row.at(*packedHoleColumn),"packed_hole_qf_increment",nodeId);
            const Real scale=parseRestartStateReal(row.at(*packedScaleColumn),"packed_potential_scale_V",nodeId);
            if(scale<=0. || (solution.packedPotentialScale_V!=0. && scale!=solution.packedPotentialScale_V))
                throw std::runtime_error("DCSweep: packed potential scale must be positive and identical on every node.");
            solution.packedPotentialScale_V=scale;
        }
        if(hasSplit) {
            const int count=expectedNodeCount;
            if(row.at(*splitSchema)!="vela.split-dd-state.v1" || row.at(*splitMesh).empty())throw std::runtime_error("Unsupported split restart schema/mesh");
            if(!solution.packedMeshFingerprint.empty() && solution.packedMeshFingerprint!=row.at(*splitMesh))throw std::runtime_error("Inconsistent split mesh metadata");
            solution.packedMeshFingerprint=row.at(*splitMesh);
            solution.packedLow(rowIndex)=parseSplitTail(row.at(*lowPsi));solution.packedLow(count+rowIndex)=parseSplitTail(row.at(*lowN));solution.packedLow(2*count+rowIndex)=parseSplitTail(row.at(*lowP));
        }
        if (quantumPotentialColumn) {
            solution.electronQuantumPotential(rowIndex) = parseRestartStateReal(
                row.at(*quantumPotentialColumn), "electron_quantum_potential_V", nodeId);
        }
        if (quantumPotentialLikeColumn) {
            solution.electronQuantumPotentialLike(rowIndex) = parseRestartStateReal(
                row.at(*quantumPotentialLikeColumn),
                "electron_quantum_potential_like_V", nodeId);
        }
        if (hasReferencedQf) {
            solution.phinIncrement(rowIndex) = parseRestartStateReal(
                row.at(*electronIncrementColumn), "electron_qf_increment_V", nodeId);
            solution.phipIncrement(rowIndex) = parseRestartStateReal(
                row.at(*holeIncrementColumn), "hole_qf_increment_V", nodeId);
            solution.electronQfReference(rowIndex) = parseRestartStateReal(
                row.at(*electronReferenceColumn), "electron_qf_reference_V", nodeId);
            solution.holeQfReference(rowIndex) = parseRestartStateReal(
                row.at(*holeReferenceColumn), "hole_qf_reference_V", nodeId);
        }
    }

    for (Index nodeId = 0; nodeId < expectedNodeCount; ++nodeId) {
        if (!seen.at(nodeId)) {
            throw std::runtime_error(
                "DCSweep: initial_state_file missing row for node id " +
                std::to_string(nodeId));
        }
    }
    if (hasReferencedQf && expectedNodeCount > 0) {
        solution.electronQfReference_V = solution.electronQfReference(0);
        solution.holeQfReference_V = solution.holeQfReference(0);
        for (int i = 0; i < solution.phin.size(); ++i) {
            const Real reconstructedElectronQf =
                solution.electronQfReference(i) + solution.phinIncrement(i);
            const Real reconstructedHoleQf =
                solution.holeQfReference(i) + solution.phipIncrement(i);
            const Real electronTolerance = 32.0 * std::numeric_limits<Real>::epsilon() *
                std::max({Real{1.0}, std::abs(solution.phin(i)),
                          std::abs(reconstructedElectronQf)});
            const Real holeTolerance = 32.0 * std::numeric_limits<Real>::epsilon() *
                std::max({Real{1.0}, std::abs(solution.phip(i)),
                          std::abs(reconstructedHoleQf)});
            if (std::abs(solution.phin(i) - reconstructedElectronQf) >
                    electronTolerance ||
                std::abs(solution.phip(i) - reconstructedHoleQf) >
                    holeTolerance) {
                throw std::runtime_error(
                    "DCSweep: initial_state_file physical and referenced "
                    "quasi-Fermi values are inconsistent at node id " +
                    std::to_string(i));
            }
        }
    }
    if(hasPacked && !solution.hasConsistentPackedState())
        throw std::runtime_error("DCSweep: packed restart coordinates do not match the physical fields/increments.");
    return solution;
}

void writeDDSolutionStateCsv(const std::filesystem::path& path,
                             const DDSolution& solution,
                             UnitScalingConfig scaling)
{
    if(solution.packedLow.size()!=0 && !solution.hasConsistentSplitPackedState())throw std::runtime_error("Cannot write stale or incomplete split checkpoint");
    const bool hasSplit=solution.packedLow.size()!=0;
    const auto fieldSize = solution.psi.size();
    if (solution.phin.size() != fieldSize ||
        solution.phip.size() != fieldSize ||
        solution.n.size() != fieldSize ||
        solution.p.size() != fieldSize) {
        throw std::runtime_error("DCSweep: cannot write restart state with inconsistent field sizes.");
    }

    const PhysicalUnitSystem& units = scaling.unitSystem();
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    if (!output.is_open())
        throw std::runtime_error("DCSweep: cannot open write_state_file: " + path.string());

    const bool hasQuantumPotential =
        solution.electronQuantumPotential.size() == fieldSize;
    const bool hasQuantumPotentialLike =
        hasQuantumPotential &&
        solution.electronQuantumPotentialLike.size() == fieldSize;
    const bool hasReferencedQf =
        solution.phinIncrement.size() == fieldSize &&
        solution.phipIncrement.size() == fieldSize;
    const bool hasPacked=solution.hasConsistentPackedState();
    const bool hasElectronReferenceField =
        solution.electronQfReference.size() == fieldSize;
    const bool hasHoleReferenceField =
        solution.holeQfReference.size() == fieldSize;
    if ((solution.phinIncrement.size() != 0 || solution.phipIncrement.size() != 0) &&
        !hasReferencedQf) {
        throw std::runtime_error(
            "DCSweep: cannot write restart state with partial quasi-Fermi increments.");
    }
    if ((solution.electronQfReference.size() != 0 && !hasElectronReferenceField) ||
        (solution.holeQfReference.size() != 0 && !hasHoleReferenceField)) {
        throw std::runtime_error(
            "DCSweep: cannot write restart state with inconsistent quasi-Fermi "
            "reference field sizes.");
    }
    output << "node_id,psi,phin,phip,electrons_m3,holes_m3";
    if (hasQuantumPotential)
        output << ",electron_quantum_potential_V";
    if (hasQuantumPotentialLike)
        output << ",electron_quantum_potential_like_V";
    if (hasReferencedQf) {
        output << ",electron_qf_increment_V,hole_qf_increment_V"
                  ",electron_qf_reference_V,hole_qf_reference_V";
    }
    if(hasPacked)
        output << ",packed_psi,packed_electron_qf_increment,packed_hole_qf_increment,packed_potential_scale_V";
    if(hasSplit)output << ",packed_psi_low,packed_electron_qf_increment_low,packed_hole_qf_increment_low,split_state_schema,split_mesh_fingerprint";
    output << '\n';
    for (int i = 0; i < fieldSize; ++i) {
        output << i << ','
               << formatRestartReal(solution.psi(i)) << ','
               << formatRestartReal(solution.phin(i)) << ','
               << formatRestartReal(solution.phip(i)) << ','
               << formatRestartReal(units.internalConcentrationToM3(solution.n(i))) << ','
               << formatRestartReal(units.internalConcentrationToM3(solution.p(i)));
        if (hasQuantumPotential)
            output << ',' << formatRestartReal(solution.electronQuantumPotential(i));
        if (hasQuantumPotentialLike)
            output << ',' << formatRestartReal(
                solution.electronQuantumPotentialLike(i));
        if (hasReferencedQf) {
            const Real electronReference = hasElectronReferenceField
                ? solution.electronQfReference(i)
                : solution.electronQfReference_V;
            const Real holeReference = hasHoleReferenceField
                ? solution.holeQfReference(i)
                : solution.holeQfReference_V;
            output << ',' << formatRestartReal(solution.phinIncrement(i))
                   << ',' << formatRestartReal(solution.phipIncrement(i))
                   << ',' << formatRestartReal(electronReference)
                   << ',' << formatRestartReal(holeReference);
        }
        if(hasPacked)
            output << ',' << formatRestartReal(solution.packedState(i))
                   << ',' << formatRestartReal(solution.packedState(fieldSize+i))
                   << ',' << formatRestartReal(solution.packedState(2*fieldSize+i))
                   << ',' << formatRestartReal(solution.packedPotentialScale_V);
        if(hasSplit)output << std::setprecision(17) << ',' << solution.packedLow(i) << ',' << solution.packedLow(fieldSize+i) << ',' << solution.packedLow(2*fieldSize+i) << ",vela.split-dd-state.v1," << solution.packedMeshFingerprint;
        output << '\n';
    }
}

} // namespace vela
