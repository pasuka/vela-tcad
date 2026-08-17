#include "vela/core/UnitScaling.h"
#include "vela/core/UnitScalingSystem.h"
#include "vela/equation/FixedStateOperatorAudit.h"
#include "vela/io/CsvUtils.h"
#include "vela/io/MeshReader.h"
#include "vela/solver/GummelSolver.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using vela::Index;
using vela::Real;

struct AuditConfiguration {
    vela::GummelConfig solver;
    vela::MaterialDatabase materials;
};

struct Arguments {
    std::filesystem::path mesh;
    std::filesystem::path doping;
    std::filesystem::path state;
    std::filesystem::path config;
    std::filesystem::path nodeOut;
    std::filesystem::path edgeOut;
    std::filesystem::path triangleOut;
    std::optional<std::filesystem::path> elementOut;
    std::optional<std::filesystem::path> processOut;
    bool generalTri3 = false;
};

std::string usage()
{
    return
        "pn2d_minimal6_operator_audit "
        "--mesh mesh.json --doping doping.csv --state state.csv "
        "--config audit.json --node-out vela_node_state.csv "
        "--edge-out vela_edge_audit.csv --triangle-out vela_triangle_audit.csv "
        "[--element-out vela_element_edge_gss_laux.csv] "
        "[--process-out vela_bv_process_probe.csv] "
        "[--scope minimal6|general_tri3]";
}

Arguments parseArguments(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << usage() << '\n';
        std::exit(0);
    }
    if (argc != 15 && argc != 17 && argc != 19 && argc != 21)
        throw std::invalid_argument("expected seven to ten option/value pairs; usage: " + usage());

    std::map<std::string, std::filesystem::path> values;
    for (int i = 1; i < argc; i += 2) {
        const std::string option = argv[i];
        if (option != "--mesh" && option != "--doping" && option != "--state" &&
            option != "--config" && option != "--node-out" &&
            option != "--edge-out" && option != "--triangle-out" &&
            option != "--element-out" && option != "--process-out" &&
            option != "--scope") {
            throw std::invalid_argument("unknown option '" + option + "'; usage: " + usage());
        }
        if (!values.emplace(option, argv[i + 1]).second)
            throw std::invalid_argument("duplicate option '" + option + "'");
    }
    const auto required = [&](const char* option) -> std::filesystem::path {
        const auto it = values.find(option);
        if (it == values.end() || it->second.empty())
            throw std::invalid_argument(std::string("missing required option ") + option);
        return it->second;
    };
    const auto optional = [&](const char* option)
        -> std::optional<std::filesystem::path> {
        const auto it = values.find(option);
        if (it == values.end())
            return std::nullopt;
        if (it->second.empty())
            throw std::invalid_argument(
                std::string("empty optional path ") + option);
        return it->second;
    };
    const std::string scope = optional("--scope").value_or("minimal6").string();
    if (scope != "minimal6" && scope != "general_tri3")
        throw std::invalid_argument("--scope must be minimal6 or general_tri3");
    return {
        required("--mesh"),
        required("--doping"),
        required("--state"),
        required("--config"),
        required("--node-out"),
        required("--edge-out"),
        required("--triangle-out"),
        optional("--element-out"),
        optional("--process-out"),
        scope == "general_tri3",
    };
}

std::ifstream openInput(const std::filesystem::path& path, const char* label)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error(std::string("cannot open ") + label + ": " + path.string());
    return input;
}

std::ofstream openOutput(const std::filesystem::path& path, const char* label)
{
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error(std::string("cannot open ") + label + ": " + path.string());
    output << std::setprecision(std::numeric_limits<Real>::max_digits10);
    return output;
}

Index parseNodeId(const std::string& token, Index nodeCount, const char* label)
{
    std::size_t parsed = 0;
    const unsigned long long raw = std::stoull(vela::trimCsvToken(token), &parsed);
    const std::string trimmed = vela::trimCsvToken(token);
    if (parsed != trimmed.size() || raw >= nodeCount)
        throw std::runtime_error(std::string(label) + " has invalid node_id '" + token + "'");
    return static_cast<Index>(raw);
}

Real parseFinite(const std::string& token, const char* column)
{
    const std::string trimmed = vela::trimCsvToken(token);
    std::size_t parsed = 0;
    const Real value = std::stod(trimmed, &parsed);
    if (parsed != trimmed.size() || !std::isfinite(value))
        throw std::runtime_error(std::string("invalid finite value in column ") + column);
    return value;
}

std::vector<std::string> readHeader(std::ifstream& input,
                                    const std::filesystem::path& path)
{
    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("empty CSV: " + path.string());
    return vela::splitCsvLine(line, "quoted CSV fields are not supported");
}

std::map<std::string, std::size_t> columnMap(
    const std::vector<std::string>& header)
{
    std::map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (!columns.emplace(header[i], i).second)
            throw std::runtime_error("duplicate CSV column '" + header[i] + "'");
    }
    return columns;
}

vela::DopingModel readDoping(const std::filesystem::path& path,
                          Index nodeCount,
                          const vela::UnitScalingConfig& scaling)
{
    std::ifstream input = openInput(path, "doping CSV");
    const auto columns = columnMap(readHeader(input, path));
    for (const char* name : {"node_id", "donors_cm3", "acceptors_cm3"}) {
        if (!columns.contains(name))
            throw std::runtime_error("doping CSV missing required column '" + std::string(name) + "'");
    }

    vela::DopingModel doping(nodeCount);
    std::vector<bool> seen(nodeCount, false);
    std::string line;
    while (std::getline(input, line)) {
        if (vela::trimCsvToken(line).empty())
            continue;
        const auto row = vela::splitCsvLine(line, "quoted CSV fields are not supported");
        if (row.size() != columns.size())
            throw std::runtime_error("doping CSV row width mismatch");
        const Index node = parseNodeId(row[columns.at("node_id")], nodeCount, "doping CSV");
        if (seen[node])
            throw std::runtime_error("doping CSV has duplicate node_id " + std::to_string(node));
        const Real donors = parseFinite(row[columns.at("donors_cm3")], "donors_cm3");
        const Real acceptors = parseFinite(row[columns.at("acceptors_cm3")], "acceptors_cm3");
        const Real donorsInternal = scaling.isUnitScaling()
            ? donors
            : donors * 1.0e6;
        const Real acceptorsInternal = scaling.isUnitScaling()
            ? acceptors
            : acceptors * 1.0e6;
        doping.setNodeDoping(node, donorsInternal, acceptorsInternal);
        seen[node] = true;
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end())
        throw std::runtime_error("doping CSV is missing one or more node rows");
    return doping;
}

vela::DDSolution readState(const std::filesystem::path& path,
                           Index nodeCount,
                           const vela::UnitScalingConfig& scaling)
{
    std::ifstream input = openInput(path, "state CSV");
    const std::vector<std::string> expected = {
        "node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"};
    const auto header = readHeader(input, path);
    if (header != expected) {
        throw std::runtime_error(
            "state CSV header must be exactly "
            "node_id,psi_V,phin_V,phip_V,n_m3,p_m3");
    }

    vela::DDSolution state;
    state.psi.resize(static_cast<Eigen::Index>(nodeCount));
    state.phin.resize(static_cast<Eigen::Index>(nodeCount));
    state.phip.resize(static_cast<Eigen::Index>(nodeCount));
    state.n.resize(static_cast<Eigen::Index>(nodeCount));
    state.p.resize(static_cast<Eigen::Index>(nodeCount));
    std::vector<bool> seen(nodeCount, false);
    std::string line;
    while (std::getline(input, line)) {
        if (vela::trimCsvToken(line).empty())
            continue;
        const auto row = vela::splitCsvLine(line, "quoted CSV fields are not supported");
        if (row.size() != expected.size())
            throw std::runtime_error("state CSV row width mismatch");
        const Index node = parseNodeId(row[0], nodeCount, "state CSV");
        if (seen[node])
            throw std::runtime_error("state CSV has duplicate node_id " + std::to_string(node));
        const Eigen::Index i = static_cast<Eigen::Index>(node);
        state.psi(i) = parseFinite(row[1], "psi_V");
        state.phin(i) = parseFinite(row[2], "phin_V");
        state.phip(i) = parseFinite(row[3], "phip_V");
        state.n(i) = scaling.unitSystem().m3ToInternalConcentration(
            parseFinite(row[4], "n_m3"));
        state.p(i) = scaling.unitSystem().m3ToInternalConcentration(
            parseFinite(row[5], "p_m3"));
        seen[node] = true;
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end())
        throw std::runtime_error("state CSV is missing one or more node rows");
    return state;
}

AuditConfiguration readConfig(const std::filesystem::path& path)
{
    std::ifstream input = openInput(path, "audit JSON");
    nlohmann::json json;
    input >> json;
    json = vela::canonicalizeDeck(json);
    const vela::UnitScalingConfig scaling = vela::parseUnitScalingConfig(json);
    vela::GummelConfig solver = json.contains("solver")
        ? vela::gummelConfigFromJson(json.at("solver"), scaling)
        : vela::gummelConfigFromJson(json, scaling);
    solver.unitScalingRefs = vela::parseUnitScalingReferenceConfig(json);
    vela::MaterialDatabase materials(scaling);
    if (json.contains("materials_file")) {
        std::filesystem::path materialsPath =
            json.at("materials_file").get<std::string>();
        if (materialsPath.is_relative())
            materialsPath = path.parent_path() / materialsPath;
        materials.loadJson(materialsPath.lexically_normal().string(), scaling);
    }
    return {solver, materials};
}

void writeNodes(const std::filesystem::path& path,
                const std::vector<vela::FixedStateNodeRecord>& nodes,
                const vela::PhysicalUnitSystem& units)
{
    std::ofstream out = openOutput(path, "node output");
    out << "node_id,psi_V,phin_V,phip_V,n_m3,p_m3\n";
    for (const auto& node : nodes) {
        out << node.nodeId << ',' << node.psi << ',' << node.phin << ','
            << node.phip << ','
            << units.internalConcentrationToM3(node.n) << ','
            << units.internalConcentrationToM3(node.p) << '\n';
    }
}

void writeEdges(const std::filesystem::path& path,
                const std::vector<vela::FixedStateEdgeRecord>& edges,
                const vela::PhysicalUnitSystem& units)
{
    std::ofstream out = openOutput(path, "edge output");
    out << "edge_id,node0,node1,length_m,"
           "electron_raw_signed_flux_per_m2_s,hole_raw_signed_flux_per_m2_s,"
           "electron_midpoint_density_m3,hole_midpoint_density_m3,"
           "electron_impact_field_V_per_m,hole_impact_field_V_per_m,"
           "electron_alpha_per_m,hole_alpha_per_m,edge_area_m2\n";
    std::map<std::pair<Index, Index>, const vela::FixedStateEdgeRecord*> ordered;
    for (const auto& edge : edges) {
        if (edge.node0 >= edge.node1)
            throw std::runtime_error("edge output requires canonical node-pair keys");
        if (!ordered.emplace(std::pair{edge.node0, edge.node1}, &edge).second)
            throw std::runtime_error("edge output has a duplicate canonical node-pair key");
    }
    for (const auto& [key, record] : ordered) {
        const auto& edge = *record;
        out << edge.edgeId << ',' << key.first << ',' << key.second << ','
            << units.internalLengthToMeters(edge.length) << ','
            << units.internalContinuityParticleFluxToPerM2PerS(
                   edge.electronRawSignedFlux) << ','
            << units.internalContinuityParticleFluxToPerM2PerS(
                   edge.holeRawSignedFlux) << ','
            << units.internalConcentrationToM3(edge.electronMidpointDensity) << ','
            << units.internalConcentrationToM3(edge.holeMidpointDensity) << ','
            << units.internalElectricFieldToVPerM(edge.electronImpactField) << ','
            << units.internalElectricFieldToVPerM(edge.holeImpactField) << ','
            << units.internalInverseLengthToMInv(edge.electronAlpha) << ','
            << units.internalInverseLengthToMInv(edge.holeAlpha) << ','
            << edge.edgeArea * units.areaM2PerInternal() << '\n';
    }
}

void writeTriangleLocalEdgeHeader(std::ofstream& out, int local)
{
    const std::string prefix = "local_edge" + std::to_string(local);
    out << ',' << prefix << "_edge_id"
        << ',' << prefix << "_node0"
        << ',' << prefix << "_node1"
        << ',' << prefix << "_truncated_partial_volume_m2"
        << ',' << prefix << "_electron_cell_qf_field_V_per_m"
        << ',' << prefix << "_hole_cell_qf_field_V_per_m"
        << ',' << prefix << "_electron_midpoint_density_m3"
        << ',' << prefix << "_hole_midpoint_density_m3"
        << ',' << prefix << "_electron_mobility_m2_per_V_s"
        << ',' << prefix << "_hole_mobility_m2_per_V_s"
        << ',' << prefix << "_electron_alpha_per_m"
        << ',' << prefix << "_hole_alpha_per_m"
        << ',' << prefix << "_electron_flux_proxy_per_m2_s"
        << ',' << prefix << "_hole_flux_proxy_per_m2_s"
        << ',' << prefix << "_electron_source_integral_per_m_s"
        << ',' << prefix << "_hole_source_integral_per_m_s";
}

void writeTriangleLocalEdge(std::ofstream& out,
                            const vela::TriangleGssAvalancheSourceRecord& edge,
                            const vela::PhysicalUnitSystem& units)
{
    const Real sourceIntegralToSi =
        units.inverseLengthMInvPerInternal() *
        units.currentDensityAM2PerInternal() *
        units.areaM2PerInternal();
    out << ',' << edge.edgeId
        << ',' << edge.node0
        << ',' << edge.node1
        << ',' << edge.truncatedPartialVolume * units.areaM2PerInternal()
        << ',' << units.internalElectricFieldToVPerM(edge.electronCellQfField)
        << ',' << units.internalElectricFieldToVPerM(edge.holeCellQfField)
        << ',' << units.internalConcentrationToM3(edge.electronMidpointDensity)
        << ',' << units.internalConcentrationToM3(edge.holeMidpointDensity)
        << ',' << units.internalMobilityToM2PerVS(edge.electronMobility)
        << ',' << units.internalMobilityToM2PerVS(edge.holeMobility)
        << ',' << units.internalInverseLengthToMInv(edge.electronAlpha)
        << ',' << units.internalInverseLengthToMInv(edge.holeAlpha)
        << ',' << units.internalContinuityParticleFluxToPerM2PerS(
                       edge.electronFluxProxy)
        << ',' << units.internalContinuityParticleFluxToPerM2PerS(
                       edge.holeFluxProxy)
        << ',' << edge.electronSourceIntegral * sourceIntegralToSi
        << ',' << edge.holeSourceIntegral * sourceIntegralToSi;
}

void writeTriangles(const std::filesystem::path& path,
                    const std::vector<vela::FixedStateTriangleRecord>& triangles,
                    const vela::PhysicalUnitSystem& units)
{
    std::ofstream out = openOutput(path, "triangle output");
    out << "cell_id,node0,node1,node2,signed_double_area_m2,"
           "grad_psi_x_V_per_m,grad_psi_y_V_per_m,"
           "grad_phin_x_V_per_m,grad_phin_y_V_per_m,"
           "grad_phip_x_V_per_m,grad_phip_y_V_per_m";
    for (int local = 0; local < 3; ++local)
        writeTriangleLocalEdgeHeader(out, local);
    out << '\n';

    for (const auto& triangle : triangles) {
        if (triangle.localEdges.size() != 3)
            throw std::runtime_error("triangle output requires exactly 3 local edges");
        const Real gradientToSi = 1.0 / units.lengthMPerInternal();
        out << triangle.cellId << ',' << triangle.nodes[0] << ','
            << triangle.nodes[1] << ',' << triangle.nodes[2] << ','
            << triangle.signedDoubleArea * units.areaM2PerInternal() << ','
            << triangle.gradPsi.x() * gradientToSi << ','
            << triangle.gradPsi.y() * gradientToSi << ','
            << triangle.gradPhin.x() * gradientToSi << ','
            << triangle.gradPhin.y() * gradientToSi << ','
            << triangle.gradPhip.x() * gradientToSi << ','
            << triangle.gradPhip.y() * gradientToSi;
        for (const auto& edge : triangle.localEdges)
            writeTriangleLocalEdge(out, edge, units);
        out << '\n';
    }
}

void writeElementEdgeGssLaux(
    const std::filesystem::path& path,
    const std::vector<vela::ElementEdgeGssLauxAvalancheSourceRecord>& records,
    const vela::DeviceMesh& mesh,
    const vela::PhysicalUnitSystem& units)
{
    std::ofstream out = openOutput(path, "element-edge GSS/Laux output");
    out << "cell_id,local_index,node_id,next_node_id,edge_id,"
           "edge_length_m,edge_partial_volume_m2,vertex_measure_m2,"
           "electron_mobility_m2_per_V_s,hole_mobility_m2_per_V_s,"
           "electron_signed_edge_flux_per_m2_s,"
           "hole_signed_edge_flux_per_m2_s,"
           "electron_current_x_per_m2_s,electron_current_y_per_m2_s,"
           "electron_current_magnitude_per_m2_s,"
           "hole_current_x_per_m2_s,hole_current_y_per_m2_s,"
           "hole_current_magnitude_per_m2_s,"
           "electron_impact_field_V_per_m,hole_impact_field_V_per_m,"
           "electron_alpha_per_m,hole_alpha_per_m,"
           "electron_source_integral_per_m_s,"
           "hole_source_integral_per_m_s,"
           "combined_source_integral_per_m_s\n";

    const Real lengthToSi = units.lengthMPerInternal();
    const Real areaToSi = units.areaM2PerInternal();
    const Real fluxToSi =
        units.internalContinuityParticleFluxToPerM2PerS(1.0);
    const Real sourceIntegralToSi =
        fluxToSi * units.inverseLengthMInvPerInternal() * areaToSi;
    const Real fieldToSi =
        units.internalElectricFieldToVPerM(1.0);

    for (const auto& record : records) {
        const auto& cell = mesh.getCell(record.cellId);
        if (cell.type != vela::CellType::Tri3 || cell.node_ids.size() != 3)
            throw std::runtime_error(
                "element-edge GSS/Laux output requires Tri3 cells");
        for (std::size_t local = 0; local < 3; ++local) {
            const Index node = cell.node_ids[local];
            const Index next = cell.node_ids[(local + 1) % 3];
            out << record.cellId << ',' << local << ',' << node << ','
                << next << ',' << record.edgeIds[local] << ','
                << record.edgeLengths[local] * lengthToSi << ','
                << record.edgePartialVolumes[local] * areaToSi << ','
                << record.vertexMeasures[local] * areaToSi << ','
                << units.internalMobilityToM2PerVS(
                       record.electronMobilities[local]) << ','
                << units.internalMobilityToM2PerVS(
                       record.holeMobilities[local]) << ','
                << record.electronSignedEdgeFlux[local] * fluxToSi << ','
                << record.holeSignedEdgeFlux[local] * fluxToSi << ','
                << record.electronCurrentVector.x() * fluxToSi << ','
                << record.electronCurrentVector.y() * fluxToSi << ','
                << record.electronCurrentVector.norm() * fluxToSi << ','
                << record.holeCurrentVector.x() * fluxToSi << ','
                << record.holeCurrentVector.y() * fluxToSi << ','
                << record.holeCurrentVector.norm() * fluxToSi << ','
                << record.electronImpactField * fieldToSi << ','
                << record.holeImpactField * fieldToSi << ','
                << units.internalInverseLengthToMInv(record.electronAlpha)
                << ','
                << units.internalInverseLengthToMInv(record.holeAlpha) << ','
                << record.electronSourceIntegrals[local] *
                       sourceIntegralToSi << ','
                << record.holeSourceIntegrals[local] *
                       sourceIntegralToSi << ','
                << record.combinedSourceIntegrals[local] *
                       sourceIntegralToSi << '\n';
        }
    }
}

void writeJoinedIndices(std::ofstream& out,
                        const vela::BVProcessProbeRecord& record)
{
    for (std::size_t i = 0; i < record.scatterCount; ++i) {
        if (i != 0)
            out << ';';
        out << record.scatterNodes[i];
    }
}

void writeJoinedValues(std::ofstream& out,
                       const vela::BVProcessProbeRecord& record,
                       const std::array<Real, 6>& values,
                       Real factor)
{
    for (std::size_t i = 0; i < record.scatterCount; ++i) {
        if (i != 0)
            out << ';';
        out << values[i] * factor;
    }
}

void writeProcessProbe(
    const std::filesystem::path& path,
    const vela::BVProcessProbeResult& probe,
    const vela::PhysicalUnitSystem& units)
{
    std::ofstream out = openOutput(path, "BV process probe output");
    out << "support_kind,carrier,cell_id,local_edge,edge_id,node0,node1,"
           "psi0_V,psi1_V,quasi_fermi0_V,quasi_fermi1_V,"
           "density0_m3,density1_m3,midpoint_density_m3,"
           "electric_field_x_V_per_m,electric_field_y_V_per_m,"
           "qf_gradient_x_V_per_m,qf_gradient_y_V_per_m,"
           "low_field_mobility_m2_per_V_s,high_field_drive_V_per_m,"
           "final_mobility_m2_per_V_s,mobility_limiter,"
           "directed_sg_flux_per_m2_s,selected_flux_magnitude_per_m2_s,"
           "current_vector_x_per_m2_s,current_vector_y_per_m2_s,"
           "current_vector_provenance,impact_field_V_per_m,alpha_per_m,"
           "source_measure_m2,generation_rate_per_m3_s,"
           "source_integral_per_m_s,qG_contribution_A_per_m,"
           "scatter_nodes,source_weights,"
           "electron_residual_contributions_per_m_s,"
           "hole_residual_contributions_per_m_s,"
           "solver_coupled,contact_adjacent,zero_measure,zero_mobility,"
           "zero_alpha,reconstructed_current,directional_partition,"
           "active_branches,configuration_fingerprint,"
           "active_branch_fingerprint\n";

    const Real concentrationToSi =
        units.internalConcentrationToM3(1.0);
    const Real fieldToSi =
        units.internalElectricFieldToVPerM(1.0);
    const Real mobilityToSi =
        units.internalMobilityToM2PerVS(1.0);
    const Real fluxToSi =
        units.internalContinuityParticleFluxToPerM2PerS(1.0);
    const Real inverseLengthToSi =
        units.inverseLengthMInvPerInternal();
    const Real areaToSi =
        units.areaM2PerInternal();
    const Real sourceIntegralToSi =
        fluxToSi * inverseLengthToSi * areaToSi;
    const Real generationRateToSi =
        fluxToSi * inverseLengthToSi;

    for (const auto& record : probe.records) {
        out << record.supportKind << ',' << record.carrier << ','
            << record.cellId << ',' << record.localEdge << ','
            << record.edgeId << ',' << record.node0 << ',' << record.node1
            << ',' << record.psi0 << ',' << record.psi1 << ','
            << record.quasiFermi0 << ',' << record.quasiFermi1 << ','
            << record.density0 * concentrationToSi << ','
            << record.density1 * concentrationToSi << ','
            << record.midpointDensity * concentrationToSi << ','
            << record.electricFieldVector.x() * fieldToSi << ','
            << record.electricFieldVector.y() * fieldToSi << ','
            << record.quasiFermiGradientVector.x() * fieldToSi << ','
            << record.quasiFermiGradientVector.y() * fieldToSi << ','
            << record.lowFieldMobility * mobilityToSi << ','
            << record.highFieldDrive * fieldToSi << ','
            << record.finalMobility * mobilityToSi << ','
            << record.mobilityLimiter << ','
            << record.directedSgFlux * fluxToSi << ','
            << record.selectedFluxMagnitude * fluxToSi << ','
            << record.currentVector.x() * fluxToSi << ','
            << record.currentVector.y() * fluxToSi << ','
            << record.currentVectorProvenance << ','
            << record.impactField * fieldToSi << ','
            << record.alpha * inverseLengthToSi << ','
            << record.sourceMeasure * areaToSi << ','
            << record.generationRate * generationRateToSi << ','
            << record.sourceIntegral * sourceIntegralToSi << ','
            << record.qGContribution * sourceIntegralToSi << ',';
        writeJoinedIndices(out, record);
        out << ',';
        writeJoinedValues(out, record, record.sourceWeights, 1.0);
        out << ',';
        writeJoinedValues(
            out, record, record.electronResidualContributions,
            sourceIntegralToSi);
        out << ',';
        writeJoinedValues(
            out, record, record.holeResidualContributions,
            sourceIntegralToSi);
        out << ',' << (record.solverCoupled ? 1 : 0)
            << ',' << (record.contactAdjacent ? 1 : 0)
            << ',' << (record.zeroMeasure ? 1 : 0)
            << ',' << (record.zeroMobility ? 1 : 0)
            << ',' << (record.zeroAlpha ? 1 : 0)
            << ',' << (record.reconstructedCurrent ? 1 : 0)
            << ',' << (record.directionalPartition ? 1 : 0)
            << ',' << record.activeBranches
            << ',' << record.configurationFingerprint
            << ',' << record.activeBranchFingerprint << '\n';
    }
}
} // namespace

int main(int argc, char** argv)
{
    try {
        const Arguments args = parseArguments(argc, argv);
        const AuditConfiguration auditConfig = readConfig(args.config);
        const vela::GummelConfig& config = auditConfig.solver;
        vela::JsonMeshReader reader;
        const vela::DeviceMesh mesh = reader.read(args.mesh.string());
        const vela::DopingModel doping =
            readDoping(args.doping, mesh.numNodes(), config.inputScaling);
        const vela::DDSolution state =
            readState(args.state, mesh.numNodes(), config.inputScaling);

        const vela::FixedStateOperatorAuditResult result =
            vela::evaluateFixedStateOperators(
                mesh, doping, state, config, auditConfig.materials,
                vela::FixedStateOperatorAuditOptions{args.generalTri3});
        if (!args.generalTri3 &&
            (result.nodes.size() != 6 || result.edges.size() != 9 ||
             result.triangles.size() != 4)) {
            throw std::runtime_error("fixed-state audit result counts are not exactly 6/9/4");
        }
        if (!args.generalTri3 && args.elementOut &&
            result.elementEdgeGssLauxTriangles.size() != 4) {
            throw std::runtime_error(
                "element-edge GSS/Laux audit result count is not exactly 4");
        }

        const auto& units = config.inputScaling.unitSystem();
        writeNodes(args.nodeOut, result.nodes, units);
        if (args.elementOut)
            writeElementEdgeGssLaux(
                *args.elementOut, result.elementEdgeGssLauxTriangles, mesh, units);
        if (args.processOut)
            writeProcessProbe(*args.processOut, result.processProbe, units);
        writeEdges(args.edgeOut, result.edges, units);
        writeTriangles(args.triangleOut, result.triangles, units);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "pn2d_minimal6_operator_audit: " << error.what() << '\n';
        return 1;
    }
}
