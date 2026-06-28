#include "vela/io/CsvUtils.h"
#include "vela/physics/ImpactIonizationModel.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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
    std::optional<double> selfElectronQf;
    std::optional<double> selfHoleQf;
    std::optional<double> sentElectronAlpha_cm_inv;
    std::optional<double> sentHoleAlpha_cm_inv;
};

struct Element {
    int n0 = 0;
    int n1 = 0;
    int n2 = 0;
};

struct OutputRow {
    double bias = 0.0;
    std::string biasText;
    int nodeId = 0;
    double xUm = 0.0;
    double yUm = 0.0;
    double fnSent = 0.0;
    double fnSelf = 0.0;
    double fnRatio = 0.0;
    double alphaNSent = 0.0;
    double alphaNFromSentF = 0.0;
    double alphaNFromSelfF = 0.0;
    double alphaNFromSentFRatio = 0.0;
    double alphaNFromSelfFRatio = 0.0;
    double fpSent = 0.0;
    double fpSelf = 0.0;
    double fpRatio = 0.0;
    double alphaPSent = 0.0;
    double alphaPFromSentF = 0.0;
    double alphaPFromSelfF = 0.0;
    double alphaPFromSentFRatio = 0.0;
    double alphaPFromSelfFRatio = 0.0;
};

std::string u8s(const char8_t* text)
{
    return std::string(reinterpret_cast<const char*>(text));
}

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

std::string keyFor(const std::string& bias, int nodeId)
{
    return bias + "|" + std::to_string(nodeId);
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

    if (hasColumns(table, {"quantity", "sentaurus_value", "vela_value_scaled_to_sentaurus_units"})) {
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
            const auto self = optionalDouble(row, {"vela_value_scaled_to_sentaurus_units"});
            if (quantity == "electron_qf") {
                record.sentElectronQf = sentaurus;
                record.selfElectronQf = self;
            } else if (quantity == "hole_qf") {
                record.sentHoleQf = sentaurus;
                record.selfHoleQf = self;
            } else if (quantity == "electron_alpha_avalanche") {
                record.sentElectronAlpha_cm_inv = sentaurus;
            } else if (quantity == "hole_alpha_avalanche") {
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
                "eQuasiFermi_sentaurus", "electron_qf_sentaurus", "e_qf_sentaurus"});
            record.sentHoleQf = optionalDouble(row, {
                "hQuasiFermi_sentaurus", "hole_qf_sentaurus", "h_qf_sentaurus"});
            record.selfElectronQf = optionalDouble(row, {
                "eQuasiFermi_self", "electron_qf_self", "e_qf_self"});
            record.selfHoleQf = optionalDouble(row, {
                "hQuasiFermi_self", "hole_qf_self", "h_qf_self"});
            record.sentElectronAlpha_cm_inv = optionalDouble(row, {
                "alpha_n_sentaurus_cm_inv", "eAlphaAvalanche", "electron_alpha_avalanche"});
            record.sentHoleAlpha_cm_inv = optionalDouble(row, {
                "alpha_p_sentaurus_cm_inv", "hAlphaAvalanche", "hole_alpha_avalanche"});
            records[{bias, nodeId}] = std::move(record);
        }
    }

    std::vector<NodeRecord> out;
    out.reserve(records.size());
    for (auto& [_, record] : records) {
        if (!record.sentElectronQf || !record.sentHoleQf ||
            !record.selfElectronQf || !record.selfHoleQf ||
            !record.sentElectronAlpha_cm_inv || !record.sentHoleAlpha_cm_inv) {
            throw std::runtime_error(
                "input CSV is missing qF or alpha data for bias " + record.biasText +
                ", node " + std::to_string(record.nodeId));
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

std::string fmt(double value)
{
    if (!std::isfinite(value))
        return "";
    std::ostringstream os;
    os << std::setprecision(12) << value;
    return os.str();
}

std::vector<OutputRow> buildRows(const std::vector<NodeRecord>& records,
                                 const std::vector<Element>& elements)
{
    const auto neighbors = buildNeighbors(elements);
    std::map<double, std::vector<NodeRecord*>> byBias;
    for (const NodeRecord& record : records)
        byBias[record.bias].push_back(const_cast<NodeRecord*>(&record));

    ImpactIonizationModelConfig config = impactIonizationModelConfig("van_overstraeten");
    config.temperature_K = 300.0;
    config.referenceTemperature_K = 300.0;
    config.debugRawVanOverstraeten = true;
    const auto model = makeImpactIonizationModel(config);

    std::vector<OutputRow> rows;
    for (auto& [bias, nodes] : byBias) {
        std::unordered_map<int, NodeRecord*> nodesById;
        std::unordered_map<int, double> sentElectronQf;
        std::unordered_map<int, double> sentHoleQf;
        std::unordered_map<int, double> selfElectronQf;
        std::unordered_map<int, double> selfHoleQf;
        for (NodeRecord* node : nodes) {
            nodesById[node->nodeId] = node;
            sentElectronQf[node->nodeId] = *node->sentElectronQf;
            sentHoleQf[node->nodeId] = *node->sentHoleQf;
            selfElectronQf[node->nodeId] = *node->selfElectronQf;
            selfHoleQf[node->nodeId] = *node->selfHoleQf;
        }

        std::sort(nodes.begin(), nodes.end(), [](const NodeRecord* a, const NodeRecord* b) {
            return a->nodeId < b->nodeId;
        });
        for (const NodeRecord* node : nodes) {
            OutputRow row;
            row.bias = bias;
            row.biasText = node->biasText;
            row.nodeId = node->nodeId;
            row.xUm = node->xUm;
            row.yUm = node->yUm;
            row.fnSent = reconstructedGradient_V_per_cm(
                node->nodeId, nodesById, neighbors, sentElectronQf);
            row.fnSelf = reconstructedGradient_V_per_cm(
                node->nodeId, nodesById, neighbors, selfElectronQf);
            row.fpSent = reconstructedGradient_V_per_cm(
                node->nodeId, nodesById, neighbors, sentHoleQf);
            row.fpSelf = reconstructedGradient_V_per_cm(
                node->nodeId, nodesById, neighbors, selfHoleQf);
            row.fnRatio = ratio(row.fnSelf, row.fnSent);
            row.fpRatio = ratio(row.fpSelf, row.fpSent);
            row.alphaNSent = *node->sentElectronAlpha_cm_inv;
            row.alphaPSent = *node->sentHoleAlpha_cm_inv;
            row.alphaNFromSentF = alphaCmInv(*model, true, row.fnSent);
            row.alphaNFromSelfF = alphaCmInv(*model, true, row.fnSelf);
            row.alphaPFromSentF = alphaCmInv(*model, false, row.fpSent);
            row.alphaPFromSelfF = alphaCmInv(*model, false, row.fpSelf);
            row.alphaNFromSentFRatio = ratio(row.alphaNFromSentF, row.alphaNSent);
            row.alphaNFromSelfFRatio = ratio(row.alphaNFromSelfF, row.alphaNSent);
            row.alphaPFromSentFRatio = ratio(row.alphaPFromSentF, row.alphaPSent);
            row.alphaPFromSelfFRatio = ratio(row.alphaPFromSelfF, row.alphaPSent);
            rows.push_back(row);
        }
    }
    return rows;
}

void writeOutputCsv(const std::filesystem::path& path, const std::vector<OutputRow>& rows)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to write output CSV: " + path.string());
    out << "bias_V,node_id,x_um,y_um,Fn_sentaurus_V_per_cm,Fn_self_V_per_cm,Fn_ratio,"
           "alpha_n_sentaurus_cm_inv,alpha_n_self_from_S_F_cm_inv,"
           "alpha_n_self_from_self_F_cm_inv,alpha_n_from_S_F_ratio,"
           "alpha_n_from_self_F_ratio,Fp_sentaurus_V_per_cm,Fp_self_V_per_cm,"
           "Fp_ratio,alpha_p_sentaurus_cm_inv,alpha_p_self_from_S_F_cm_inv,"
           "alpha_p_self_from_self_F_cm_inv,alpha_p_from_S_F_ratio,"
           "alpha_p_from_self_F_ratio\n";
    for (const auto& row : rows) {
        out << row.biasText << ',' << row.nodeId << ',' << fmt(row.xUm) << ','
            << fmt(row.yUm) << ',' << fmt(row.fnSent) << ',' << fmt(row.fnSelf) << ','
            << fmt(row.fnRatio) << ',' << fmt(row.alphaNSent) << ','
            << fmt(row.alphaNFromSentF) << ',' << fmt(row.alphaNFromSelfF) << ','
            << fmt(row.alphaNFromSentFRatio) << ',' << fmt(row.alphaNFromSelfFRatio) << ','
            << fmt(row.fpSent) << ',' << fmt(row.fpSelf) << ',' << fmt(row.fpRatio) << ','
            << fmt(row.alphaPSent) << ',' << fmt(row.alphaPFromSentF) << ','
            << fmt(row.alphaPFromSelfF) << ',' << fmt(row.alphaPFromSentFRatio) << ','
            << fmt(row.alphaPFromSelfFRatio) << '\n';
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

bool inRegion(const OutputRow& row, const std::string& region)
{
    if (region == "center")
        return row.xUm >= 0.9 && row.xUm <= 1.1;
    if (region == "left")
        return row.xUm >= 0.7 && row.xUm <= 0.85;
    if (region == "right")
        return row.xUm >= 1.15 && row.xUm <= 1.3;
    return true;
}

std::string markdownRatioRow(const std::string& label, const std::vector<OutputRow>& rows)
{
    std::vector<double> nFromS;
    std::vector<double> nFromSelf;
    std::vector<double> pFromS;
    std::vector<double> pFromSelf;
    for (const auto& row : rows) {
        nFromS.push_back(row.alphaNFromSentFRatio);
        nFromSelf.push_back(row.alphaNFromSelfFRatio);
        pFromS.push_back(row.alphaPFromSentFRatio);
        pFromSelf.push_back(row.alphaPFromSelfFRatio);
    }
    auto nS = median(nFromS);
    auto nSelf = median(nFromSelf);
    auto pS = median(pFromS);
    auto pSelf = median(pFromSelf);
    std::ostringstream os;
    os << "| " << label << " | " << rows.size() << " | "
       << (nS ? fmt(*nS) : "") << " | " << (nSelf ? fmt(*nSelf) : "") << " | "
       << (pS ? fmt(*pS) : "") << " | " << (pSelf ? fmt(*pSelf) : "") << " |";
    return os.str();
}

bool nearOne(std::optional<double> value)
{
    return value && *value >= 0.5 && *value <= 2.0;
}

bool muchLessThanOne(std::optional<double> value)
{
    return value && *value > 0.0 && *value < 0.25;
}

void writeSummary(const std::filesystem::path& path, const std::vector<OutputRow>& rows)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to write summary: " + path.string());

    out << "# Alpha From Sentaurus F Diagnostic\n\n";
    out << "- F units: V/cm\n";
    out << "- alpha units: 1/cm\n";
    out << "- Gradient recovery: 1/d weighted least squares on node qF values and elements.csv connectivity\n";
    out << "- Alpha model: Vela VanOverstraeten/de Man implementation at T=300 K\n\n";

    out << "## Region Statistics\n\n";
    out << "| region | nodes | median alpha_n S_F/S | median alpha_n self_F/S | "
           "median alpha_p S_F/S | median alpha_p self_F/S |\n";
    out << "|---|---:|---:|---:|---:|---:|\n";
    const std::vector<std::pair<std::string, std::string>> regions = {
        {u8s(u8"\u7ed3\u4e2d\u5fc3"), "center"},
        {u8s(u8"\u5de6\u80a9\u90e8"), "left"},
        {u8s(u8"\u53f3\u80a9\u90e8"), "right"},
    };
    for (const auto& [label, key] : regions) {
        std::vector<OutputRow> subset;
        std::copy_if(rows.begin(), rows.end(), std::back_inserter(subset),
                     [&](const OutputRow& row) { return inRegion(row, key); });
        out << markdownRatioRow(label, subset) << '\n';
    }

    out << "\n## Bias Statistics\n\n";
    out << "| bias | nodes | median alpha_n S_F/S | median alpha_n self_F/S | "
           "median alpha_p S_F/S | median alpha_p self_F/S |\n";
    out << "|---|---:|---:|---:|---:|---:|\n";
    const std::vector<double> requestedBiases = {-5.0, -10.0, -16.0, -18.0, -20.0};
    for (double bias : requestedBiases) {
        std::vector<OutputRow> subset;
        std::copy_if(rows.begin(), rows.end(), std::back_inserter(subset),
                     [&](const OutputRow& row) { return std::abs(row.bias - bias) < 1.0e-9; });
        std::ostringstream label;
        label << std::fixed << std::setprecision(0) << bias << "V";
        out << markdownRatioRow(label.str(), subset) << '\n';
    }

    std::vector<OutputRow> centerRows;
    std::vector<OutputRow> shoulderRows;
    std::copy_if(rows.begin(), rows.end(), std::back_inserter(centerRows),
                 [](const OutputRow& row) { return inRegion(row, "center"); });
    std::copy_if(rows.begin(), rows.end(), std::back_inserter(shoulderRows),
                 [](const OutputRow& row) { return inRegion(row, "left") || inRegion(row, "right"); });
    std::vector<double> centerRatios;
    std::vector<double> shoulderRatios;
    for (const auto& row : centerRows) {
        centerRatios.push_back(row.alphaNFromSentFRatio);
        centerRatios.push_back(row.alphaPFromSentFRatio);
    }
    for (const auto& row : shoulderRows) {
        shoulderRatios.push_back(row.alphaNFromSentFRatio);
        shoulderRatios.push_back(row.alphaPFromSentFRatio);
    }
    const auto centerMedian = median(centerRatios);
    const auto shoulderMedian = median(shoulderRatios);

    out << "\n## Interpretation\n\n";
    if (nearOne(centerMedian) && muchLessThanOne(shoulderMedian)) {
        out << u8s(u8"\u5269\u4f59\u5dee\u5f02\u96c6\u4e2d\u5728\u4e2d\u4f4e\u573a\u80a9\u90e8\uff0c"
                   "\u4f18\u5148\u68c0\u67e5\u4f4e\u573a\u5206\u652f\u3001\u5e73\u6ed1\u548c\u8282\u70b9\u6620\u5c04")
            << "\n\n";
    }
    if (nearOne(centerMedian) && nearOne(shoulderMedian)) {
        out << u8s(u8"\u4f7f\u7528 Sentaurus F \u540e alpha \u660e\u663e\u63a5\u8fd1\uff0c"
                   "\u95ee\u9898\u4f18\u5148\u68c0\u67e5 F \u7684\u6062\u590d\u65b9\u5f0f\u3001"
                   "\u8282\u70b9/\u8fb9/\u5355\u5143\u4f4d\u7f6e\u548c driving force \u7684\u5b9a\u4e49")
            << "\n";
    } else {
        out << u8s(u8"\u4f7f\u7528 Sentaurus F \u540e alpha \u4ecd\u4e0d\u4e00\u81f4\uff0c"
                   "\u95ee\u9898\u4f18\u5148\u68c0\u67e5\u516c\u5f0f\u3001\u53c2\u6570\u3001\u5206\u6bb5\u3001"
                   "cutoff/smoothing \u548c\u5355\u4f4d")
            << "\n";
    }
}

void usage()
{
    std::cerr
        << "Usage: alpha_from_sentaurus_f --input-csv FILE --elements-csv FILE "
           "[--output-csv FILE] [--summary-md FILE]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::filesystem::path inputCsv;
        std::filesystem::path elementsCsv;
        std::filesystem::path outputCsv = "build/diagnostics/alpha_from_sentaurus_F.csv";
        std::filesystem::path summaryMd = "build/diagnostics/alpha_from_sentaurus_F_summary.md";

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
        std::cerr << "alpha_from_sentaurus_f: " << ex.what() << '\n';
        return 1;
    }
}
