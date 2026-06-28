#include "vela/io/CsvUtils.h"
#include "vela/physics/ImpactIonizationModel.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace vela;

namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::unordered_map<std::string, std::string>> rows;
};

struct NodeRecord {
    std::string biasText;
    double bias = 0.0;
    int nodeId = 0;
    double xUm = 0.0;
    double yUm = 0.0;
    std::optional<double> sentElectronQf;
    std::optional<double> sentHoleQf;
    std::optional<double> sentElectronAlpha_cm_inv;
    std::optional<double> sentHoleAlpha_cm_inv;
};

struct Element {
    int n0 = 0;
    int n1 = 0;
    int n2 = 0;
};

struct ModelSet {
    std::string name;
    ImpactIonizationModelConfig config;
    std::unique_ptr<ImpactIonizationModel> model;
};

struct OutputRow {
    std::string biasText;
    int nodeId = 0;
    double xUm = 0.0;
    double yUm = 0.0;
    std::string regionLabel;
    std::string carrier;
    double field = 0.0;
    double alphaSentaurus = 0.0;
    double alphaDefault = 0.0;
    double alphaFitAOnly = 0.0;
    double alphaFitAB = 0.0;
    double alphaFitABSwitch = 0.0;
    double ratioDefault = 0.0;
    double ratioFitAOnly = 0.0;
    double ratioFitAB = 0.0;
    double ratioFitABSwitch = 0.0;
    std::string branchDefault;
    std::string branchFitABSwitch;
};

std::string requireValue(const std::vector<std::string>& args, int& index, const std::string& option)
{
    if (index + 1 >= static_cast<int>(args.size()))
        throw std::runtime_error("missing value for " + option);
    return args[++index];
}

double parseDouble(const std::string& text, const std::string& context)
{
    try {
        std::size_t consumed = 0;
        const double value = std::stod(text, &consumed);
        if (consumed != text.size() || !std::isfinite(value))
            throw std::runtime_error("");
        return value;
    } catch (...) {
        throw std::runtime_error("invalid numeric value for " + context + ": '" + text + "'");
    }
}

int parseInt(const std::string& text, const std::string& context)
{
    try {
        std::size_t consumed = 0;
        const int value = std::stoi(text, &consumed);
        if (consumed != text.size())
            throw std::runtime_error("");
        return value;
    } catch (...) {
        throw std::runtime_error("invalid integer value for " + context + ": '" + text + "'");
    }
}

std::optional<double> optionalDouble(const std::unordered_map<std::string, std::string>& row,
                                     const std::vector<std::string>& names)
{
    for (const std::string& name : names) {
        const auto it = row.find(name);
        if (it != row.end() && !it->second.empty())
            return parseDouble(it->second, name);
    }
    return std::nullopt;
}

std::string optionalText(const std::unordered_map<std::string, std::string>& row,
                         const std::vector<std::string>& names)
{
    for (const std::string& name : names) {
        const auto it = row.find(name);
        if (it != row.end())
            return it->second;
    }
    return {};
}

CsvTable readCsv(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("failed to open CSV: " + path.string());

    CsvTable table;
    std::string line;
    if (!std::getline(in, line))
        return table;
    table.header = splitCsvLine(line);

    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> cells = splitCsvLine(line);
        std::unordered_map<std::string, std::string> row;
        for (std::size_t i = 0; i < table.header.size() && i < cells.size(); ++i)
            row[table.header[i]] = cells[i];
        table.rows.push_back(std::move(row));
    }
    return table;
}

bool hasColumns(const CsvTable& table, const std::vector<std::string>& names)
{
    const std::set<std::string> header(table.header.begin(), table.header.end());
    for (const auto& name : names) {
        if (!header.contains(name))
            return false;
    }
    return true;
}

std::vector<NodeRecord> loadNodeRecords(const std::filesystem::path& path)
{
    const CsvTable table = readCsv(path);
    std::map<std::pair<double, int>, NodeRecord> records;

    if (hasColumns(table, {"quantity", "sentaurus_value"})) {
        for (const auto& row : table.rows) {
            const std::string biasText = optionalText(row, {"bias_V"});
            const int nodeId = parseInt(optionalText(row, {"node_id"}), "node_id");
            const double bias = parseDouble(biasText, "bias_V");
            NodeRecord& record = records[{bias, nodeId}];
            record.biasText = biasText;
            record.bias = bias;
            record.nodeId = nodeId;
            record.xUm = parseDouble(optionalText(row, {"x_um"}), "x_um");
            record.yUm = parseDouble(optionalText(row, {"y_um"}), "y_um");

            const std::string quantity = optionalText(row, {"quantity"});
            const auto sentaurus = optionalDouble(row, {"sentaurus_value"});
            if (quantity == "electron_qf" || quantity == "eQuasiFermi") {
                record.sentElectronQf = sentaurus;
            } else if (quantity == "hole_qf" || quantity == "hQuasiFermi") {
                record.sentHoleQf = sentaurus;
            } else if (quantity == "electron_alpha_avalanche" ||
                       quantity == "eAlphaAvalanche") {
                record.sentElectronAlpha_cm_inv = sentaurus;
            } else if (quantity == "hole_alpha_avalanche" ||
                       quantity == "hAlphaAvalanche") {
                record.sentHoleAlpha_cm_inv = sentaurus;
            }
        }
    } else {
        for (const auto& row : table.rows) {
            const std::string biasText = optionalText(row, {"bias_V"});
            const int nodeId = parseInt(optionalText(row, {"node_id"}), "node_id");
            const double bias = parseDouble(biasText, "bias_V");
            NodeRecord record;
            record.biasText = biasText;
            record.bias = bias;
            record.nodeId = nodeId;
            record.xUm = parseDouble(optionalText(row, {"x_um"}), "x_um");
            record.yUm = parseDouble(optionalText(row, {"y_um"}), "y_um");
            record.sentElectronQf = optionalDouble(row, {
                "eQuasiFermi_sentaurus",
                "electron_qf_sentaurus",
                "e_qf_sentaurus",
                "eQuasiFermi",
                "electron_qf"});
            record.sentHoleQf = optionalDouble(row, {
                "hQuasiFermi_sentaurus",
                "hole_qf_sentaurus",
                "h_qf_sentaurus",
                "hQuasiFermi",
                "hole_qf"});
            record.sentElectronAlpha_cm_inv = optionalDouble(row, {
                "alpha_n_sentaurus_cm_inv",
                "eAlphaAvalanche",
                "electron_alpha_avalanche"});
            record.sentHoleAlpha_cm_inv = optionalDouble(row, {
                "alpha_p_sentaurus_cm_inv",
                "hAlphaAvalanche",
                "hole_alpha_avalanche"});
            records[{bias, nodeId}] = std::move(record);
        }
    }

    std::vector<NodeRecord> out;
    out.reserve(records.size());
    for (auto& [_, record] : records) {
        if (!record.sentElectronQf || !record.sentHoleQf ||
            !record.sentElectronAlpha_cm_inv || !record.sentHoleAlpha_cm_inv) {
            throw std::runtime_error(
                "input CSV is missing Sentaurus qF or alpha data for bias " +
                record.biasText + ", node " + std::to_string(record.nodeId));
        }
        out.push_back(std::move(record));
    }
    return out;
}

std::vector<Element> loadElements(const std::filesystem::path& path)
{
    const CsvTable table = readCsv(path);
    std::vector<Element> elements;
    for (const auto& row : table.rows) {
        Element element;
        element.n0 = parseInt(optionalText(row, {"node0", "n0"}), "node0");
        element.n1 = parseInt(optionalText(row, {"node1", "n1"}), "node1");
        element.n2 = parseInt(optionalText(row, {"node2", "n2"}), "node2");
        elements.push_back(element);
    }
    return elements;
}

std::unordered_map<int, std::vector<int>> buildNeighbors(const std::vector<Element>& elements)
{
    std::unordered_map<int, std::set<int>> sets;
    for (const Element& e : elements) {
        const int ids[3] = {e.n0, e.n1, e.n2};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i != j)
                    sets[ids[i]].insert(ids[j]);
            }
        }
    }

    std::unordered_map<int, std::vector<int>> neighbors;
    for (auto& [node, values] : sets)
        neighbors[node] = std::vector<int>(values.begin(), values.end());
    return neighbors;
}

double reconstructedGradient_V_per_cm(
    int nodeId,
    const std::unordered_map<int, NodeRecord*>& nodesById,
    const std::unordered_map<int, std::vector<int>>& neighbors,
    const std::unordered_map<int, double>& values)
{
    const auto centerNodeIt = nodesById.find(nodeId);
    const auto centerValueIt = values.find(nodeId);
    const auto neighborsIt = neighbors.find(nodeId);
    if (centerNodeIt == nodesById.end() || centerValueIt == values.end() ||
        neighborsIt == neighbors.end()) {
        return 0.0;
    }

    const NodeRecord& center = *centerNodeIt->second;
    const double centerValue = centerValueIt->second;
    double sxx = 0.0;
    double sxy = 0.0;
    double syy = 0.0;
    double sxv = 0.0;
    double syv = 0.0;
    for (const int neighborId : neighborsIt->second) {
        const auto neighborNodeIt = nodesById.find(neighborId);
        const auto neighborValueIt = values.find(neighborId);
        if (neighborNodeIt == nodesById.end() || neighborValueIt == values.end())
            continue;
        const NodeRecord& neighbor = *neighborNodeIt->second;
        const double dx = neighbor.xUm - center.xUm;
        const double dy = neighbor.yUm - center.yUm;
        const double distance = std::hypot(dx, dy);
        if (distance <= 1.0e-30)
            continue;
        const double weight = 1.0 / distance;
        const double dv = neighborValueIt->second - centerValue;
        sxx += weight * dx * dx;
        sxy += weight * dx * dy;
        syy += weight * dy * dy;
        sxv += weight * dx * dv;
        syv += weight * dy * dv;
    }

    const double det = sxx * syy - sxy * sxy;
    if (std::abs(det) <= 1.0e-60)
        return 0.0;
    const double gradX_V_per_um = (sxv * syy - syv * sxy) / det;
    const double gradY_V_per_um = (sxx * syv - sxy * sxv) / det;
    return std::hypot(gradX_V_per_um, gradY_V_per_um) * 1.0e4;
}

std::string regionLabel(double xUm)
{
    if (xUm >= 0.9 && xUm <= 1.1)
        return "center";
    if (xUm >= 0.7 && xUm <= 0.85)
        return "left_shoulder";
    if (xUm >= 1.15 && xUm <= 1.3)
        return "right_shoulder";
    return "other";
}

double alphaCmInv(const ImpactIonizationModel& model, bool electron, double field_V_per_cm)
{
    const double alpha_m_inv = electron
        ? model.electronCoefficient(field_V_per_cm * 100.0)
        : model.holeCoefficient(field_V_per_cm * 100.0);
    return alpha_m_inv / 100.0;
}

double ratio(double value, double reference)
{
    return reference != 0.0 ? value / reference : 0.0;
}

std::string branchFor(const ImpactIonizationModelConfig& config, double field_V_per_cm)
{
    const ImpactIonizationModelConfig resolved = applyImpactIonizationParameterSet(config);
    return field_V_per_cm * 100.0 < resolved.switchField ? "low" : "high";
}

ModelSet makeSet(const std::string& parameterSet)
{
    ImpactIonizationModelConfig config = impactIonizationModelConfig("van_overstraeten");
    config.parameterSet = parameterSet;
    config.temperature_K = 300.0;
    config.referenceTemperature_K = 300.0;
    config.debugRawVanOverstraeten = true;
    return {parameterSet, config, makeImpactIonizationModel(config)};
}

std::vector<ModelSet> makeModelSets()
{
    std::vector<ModelSet> sets;
    sets.push_back(makeSet("default"));
    sets.push_back(makeSet("sentaurus_fit_A_only"));
    sets.push_back(makeSet("sentaurus_fit_A_B"));
    sets.push_back(makeSet("sentaurus_fit_A_B_switch"));
    return sets;
}

std::vector<OutputRow> buildRows(const std::vector<NodeRecord>& records,
                                 const std::vector<Element>& elements)
{
    const auto neighbors = buildNeighbors(elements);
    std::map<double, std::vector<NodeRecord*>> byBias;
    for (const NodeRecord& record : records)
        byBias[record.bias].push_back(const_cast<NodeRecord*>(&record));

    std::vector<ModelSet> sets = makeModelSets();
    std::vector<OutputRow> rows;
    for (auto& [_, nodes] : byBias) {
        std::unordered_map<int, NodeRecord*> nodesById;
        std::unordered_map<int, double> sentElectronQf;
        std::unordered_map<int, double> sentHoleQf;
        for (NodeRecord* node : nodes) {
            nodesById[node->nodeId] = node;
            sentElectronQf[node->nodeId] = *node->sentElectronQf;
            sentHoleQf[node->nodeId] = *node->sentHoleQf;
        }

        std::sort(nodes.begin(), nodes.end(), [](const NodeRecord* a, const NodeRecord* b) {
            return a->nodeId < b->nodeId;
        });
        for (const NodeRecord* node : nodes) {
            for (const bool electron : {true, false}) {
                OutputRow row;
                row.biasText = node->biasText;
                row.nodeId = node->nodeId;
                row.xUm = node->xUm;
                row.yUm = node->yUm;
                row.regionLabel = regionLabel(node->xUm);
                row.carrier = electron ? "electron" : "hole";
                row.field = reconstructedGradient_V_per_cm(
                    node->nodeId, nodesById, neighbors,
                    electron ? sentElectronQf : sentHoleQf);
                row.alphaSentaurus = electron
                    ? *node->sentElectronAlpha_cm_inv
                    : *node->sentHoleAlpha_cm_inv;
                row.alphaDefault = alphaCmInv(*sets[0].model, electron, row.field);
                row.alphaFitAOnly = alphaCmInv(*sets[1].model, electron, row.field);
                row.alphaFitAB = alphaCmInv(*sets[2].model, electron, row.field);
                row.alphaFitABSwitch = alphaCmInv(*sets[3].model, electron, row.field);
                row.ratioDefault = ratio(row.alphaDefault, row.alphaSentaurus);
                row.ratioFitAOnly = ratio(row.alphaFitAOnly, row.alphaSentaurus);
                row.ratioFitAB = ratio(row.alphaFitAB, row.alphaSentaurus);
                row.ratioFitABSwitch = ratio(row.alphaFitABSwitch, row.alphaSentaurus);
                row.branchDefault = branchFor(sets[0].config, row.field);
                row.branchFitABSwitch = branchFor(sets[3].config, row.field);
                rows.push_back(row);
            }
        }
    }
    return rows;
}

std::string fmt(double value)
{
    if (!std::isfinite(value))
        return "";
    std::ostringstream os;
    os << std::setprecision(12) << value;
    return os.str();
}

void writeOutputCsv(const std::filesystem::path& path, const std::vector<OutputRow>& rows)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to write output CSV: " + path.string());

    out << "bias_V,node_id,x_um,y_um,region_label,carrier,F_sentaurus_V_per_cm,"
           "alpha_sentaurus_cm_inv,alpha_default_cm_inv,alpha_fit_A_only_cm_inv,"
           "alpha_fit_A_B_cm_inv,alpha_fit_A_B_switch_cm_inv,ratio_default,"
           "ratio_fit_A_only,ratio_fit_A_B,ratio_fit_A_B_switch,branch_default,"
           "branch_fit_A_B_switch\n";
    for (const OutputRow& row : rows) {
        out << row.biasText << ',' << row.nodeId << ',' << fmt(row.xUm) << ','
            << fmt(row.yUm) << ',' << row.regionLabel << ',' << row.carrier << ','
            << fmt(row.field) << ',' << fmt(row.alphaSentaurus) << ','
            << fmt(row.alphaDefault) << ',' << fmt(row.alphaFitAOnly) << ','
            << fmt(row.alphaFitAB) << ',' << fmt(row.alphaFitABSwitch) << ','
            << fmt(row.ratioDefault) << ',' << fmt(row.ratioFitAOnly) << ','
            << fmt(row.ratioFitAB) << ',' << fmt(row.ratioFitABSwitch) << ','
            << row.branchDefault << ',' << row.branchFitABSwitch << '\n';
    }
}

std::optional<double> median(std::vector<double> values)
{
    values.erase(std::remove_if(values.begin(), values.end(), [](double value) {
        return !std::isfinite(value) || value <= 0.0;
    }), values.end());
    if (values.empty())
        return std::nullopt;
    std::sort(values.begin(), values.end());
    const std::size_t mid = values.size() / 2;
    if (values.size() % 2 == 1)
        return values[mid];
    return 0.5 * (values[mid - 1] + values[mid]);
}

double ratioBySet(const OutputRow& row, const std::string& parameterSet)
{
    if (parameterSet == "default")
        return row.ratioDefault;
    if (parameterSet == "sentaurus_fit_A_only")
        return row.ratioFitAOnly;
    if (parameterSet == "sentaurus_fit_A_B")
        return row.ratioFitAB;
    if (parameterSet == "sentaurus_fit_A_B_switch")
        return row.ratioFitABSwitch;
    return std::numeric_limits<double>::quiet_NaN();
}

double scoreRows(const std::vector<OutputRow>& rows, const std::string& parameterSet)
{
    std::vector<double> scores;
    for (const OutputRow& row : rows) {
        const double value = ratioBySet(row, parameterSet);
        if (std::isfinite(value) && value > 0.0)
            scores.push_back(std::abs(std::log10(value)));
    }
    const auto med = median(scores);
    return med ? *med : std::numeric_limits<double>::quiet_NaN();
}

std::vector<OutputRow> filterRows(const std::vector<OutputRow>& rows,
                                  const std::string& region,
                                  const std::string& carrier = "")
{
    std::vector<OutputRow> subset;
    std::copy_if(rows.begin(), rows.end(), std::back_inserter(subset),
                 [&](const OutputRow& row) {
                     return row.regionLabel == region &&
                            (carrier.empty() || row.carrier == carrier);
                 });
    return subset;
}

std::string bestParameterSet(const std::vector<OutputRow>& rows)
{
    const std::vector<std::string> parameterSets = {
        "default",
        "sentaurus_fit_A_only",
        "sentaurus_fit_A_B",
        "sentaurus_fit_A_B_switch",
    };
    std::string best = parameterSets.front();
    double bestScore = std::numeric_limits<double>::infinity();
    for (const std::string& parameterSet : parameterSets) {
        const double score = scoreRows(rows, parameterSet);
        if (std::isfinite(score) && score < bestScore) {
            bestScore = score;
            best = parameterSet;
        }
    }
    return best;
}

std::string yesNo(bool value)
{
    return value ? "yes" : "no";
}

void writeSummary(const std::filesystem::path& path, const std::vector<OutputRow>& rows)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to write summary: " + path.string());

    const std::vector<std::string> parameterSets = {
        "default",
        "sentaurus_fit_A_only",
        "sentaurus_fit_A_B",
        "sentaurus_fit_A_B_switch",
    };
    const std::vector<std::string> regions = {
        "center",
        "left_shoulder",
        "right_shoulder",
    };

    out << "# VanOverstraeten Parameter Sweep\n\n";
    out << "- F units: V/cm\n";
    out << "- alpha units: 1/cm\n";
    out << "- Input F: Sentaurus eQuasiFermi/hQuasiFermi reconstructed with 1/d weighted least squares\n";
    out << "- Parameter set under test: van_overstraeten_sentaurus_effective_fit_v1\n\n";

    out << "## Median Alpha Ratio By Region\n\n";
    out << "| region | carrier | default | sentaurus_fit_A_only | sentaurus_fit_A_B | "
           "sentaurus_fit_A_B_switch |\n";
    out << "|---|---|---:|---:|---:|---:|\n";
    for (const std::string& region : regions) {
        for (const std::string& carrier : {"electron", "hole"}) {
            const std::vector<OutputRow> subset = filterRows(rows, region, carrier);
            out << "| " << region << " | " << carrier;
            for (const std::string& parameterSet : parameterSets) {
                std::vector<double> ratios;
                for (const OutputRow& row : subset)
                    ratios.push_back(ratioBySet(row, parameterSet));
                const auto med = median(ratios);
                out << " | " << (med ? fmt(*med) : "");
            }
            out << " |\n";
        }
    }

    const std::string best = bestParameterSet(rows);
    const double defaultScore = scoreRows(rows, "default");
    const double aOnlyScore = scoreRows(rows, "sentaurus_fit_A_only");
    const double abScore = scoreRows(rows, "sentaurus_fit_A_B");
    const double abSwitchScore = scoreRows(rows, "sentaurus_fit_A_B_switch");

    std::vector<OutputRow> shoulders;
    for (const OutputRow& row : rows) {
        if (row.regionLabel == "left_shoulder" || row.regionLabel == "right_shoulder")
            shoulders.push_back(row);
    }
    const double shoulderAOnlyScore = scoreRows(shoulders, "sentaurus_fit_A_only");
    const double shoulderABScore = scoreRows(shoulders, "sentaurus_fit_A_B");

    std::vector<OutputRow> switchBand;
    for (const OutputRow& row : rows) {
        if (row.field >= 2.5e5 && row.field < 4.0e5)
            switchBand.push_back(row);
    }
    const double switchBandABScore = scoreRows(switchBand, "sentaurus_fit_A_B");
    const double switchBandABSwitchScore =
        scoreRows(switchBand, "sentaurus_fit_A_B_switch");

    const std::vector<OutputRow> electronRows = filterRows(rows, "center", "electron");
    const std::vector<OutputRow> holeRows = filterRows(rows, "center", "hole");
    const double electronScore = scoreRows(electronRows, best);
    const double holeScore = scoreRows(holeRows, best);
    const bool carrierConsistent = std::isfinite(electronScore) && std::isfinite(holeScore) &&
        std::abs(electronScore - holeScore) <= 0.5;

    const bool aOnlyImprovesOrder =
        std::isfinite(defaultScore) && std::isfinite(aOnlyScore) &&
        aOnlyScore + 0.25 < defaultScore;
    const bool abImprovesShoulder =
        std::isfinite(shoulderAOnlyScore) && std::isfinite(shoulderABScore) &&
        shoulderABScore + 0.1 < shoulderAOnlyScore;
    const bool switchImprovesBand =
        std::isfinite(switchBandABScore) && std::isfinite(switchBandABSwitchScore) &&
        switchBandABSwitchScore + 0.1 < switchBandABScore;

    out << "\n## Answers\n\n";
    out << "1. Best parameter_set: " << best << "\n";
    out << "2. A_only fixes a broad 1/100 or 1/30 scale error: "
        << yesNo(aOnlyImprovesOrder) << "\n";
    out << "3. A_B improves shoulders: " << yesNo(abImprovesShoulder) << "\n";
    out << "4. switchField=2.5e5 improves the 2.5e5~4e5 V/cm interval: "
        << yesNo(switchImprovesBand) << "\n";
    out << "5. electron/hole behavior is consistent: " << yesNo(carrierConsistent)
        << "\n";
    if (best == "sentaurus_fit_A_B_switch" &&
        std::isfinite(abSwitchScore) && std::isfinite(defaultScore) &&
        abSwitchScore + 0.2 < defaultScore) {
        out << "6. sentaurus_fit_A_B_switch is clearly best; next step: run full BV simulation with it.\n";
    } else {
        out << "6. sentaurus_fit_A_B_switch is not clearly best on this node-output sweep.\n";
    }
    if (std::isfinite(defaultScore) &&
        std::isfinite(aOnlyScore) &&
        std::isfinite(abScore) &&
        std::isfinite(abSwitchScore) &&
        aOnlyScore >= defaultScore &&
        abScore >= defaultScore &&
        abSwitchScore >= defaultScore) {
        out << "7. All schemes fail to improve shoulders; continue checking cutoff/smoothing/RefDens or node/edge/cell alpha mapping.\n";
    } else {
        out << "7. At least one parameter-set scheme improves the node-output ratios; use the CSV to inspect remaining shoulder gaps.\n";
    }
}

void usage()
{
    std::cerr
        << "Usage: vanoverstraeten_parameter_sweep --input-csv FILE --elements-csv FILE "
           "[--output-csv FILE] [--summary-md FILE]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::filesystem::path inputCsv;
        std::filesystem::path elementsCsv;
        std::filesystem::path outputCsv =
            "build/diagnostics/vanoverstraeten_parameter_sweep.csv";
        std::filesystem::path summaryMd =
            "build/diagnostics/vanoverstraeten_parameter_sweep_summary.md";

        std::vector<std::string> args(argv + 1, argv + argc);
        for (int i = 0; i < static_cast<int>(args.size()); ++i) {
            const std::string& arg = args[i];
            if (arg == "--input-csv") {
                inputCsv = requireValue(args, i, arg);
            } else if (arg == "--elements-csv") {
                elementsCsv = requireValue(args, i, arg);
            } else if (arg == "--output-csv") {
                outputCsv = requireValue(args, i, arg);
            } else if (arg == "--summary-md") {
                summaryMd = requireValue(args, i, arg);
            } else if (arg == "--help" || arg == "-h") {
                usage();
                return 0;
            } else {
                throw std::runtime_error("unknown argument: " + arg);
            }
        }
        if (inputCsv.empty() || elementsCsv.empty()) {
            usage();
            return 2;
        }

        const std::vector<NodeRecord> records = loadNodeRecords(inputCsv);
        const std::vector<Element> elements = loadElements(elementsCsv);
        const std::vector<OutputRow> rows = buildRows(records, elements);
        writeOutputCsv(outputCsv, rows);
        writeSummary(summaryMd, rows);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "vanoverstraeten_parameter_sweep: " << ex.what() << '\n';
        return 1;
    }
}
