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

struct FitSample {
    std::string biasText;
    double bias = 0.0;
    double field_V_per_cm = 0.0;
    double alpha_cm_inv = 0.0;
};

struct FitResult {
    std::string carrier;
    std::string fitScope;
    std::string biasText;
    std::string fieldRegion;
    std::size_t sampleCount = 0;
    double fMin = std::numeric_limits<double>::quiet_NaN();
    double fMax = std::numeric_limits<double>::quiet_NaN();
    double alphaMin = std::numeric_limits<double>::quiet_NaN();
    double alphaMax = std::numeric_limits<double>::quiet_NaN();
    double c = std::numeric_limits<double>::quiet_NaN();
    double d = std::numeric_limits<double>::quiet_NaN();
    double aEff = std::numeric_limits<double>::quiet_NaN();
    double bEff = std::numeric_limits<double>::quiet_NaN();
    double rSquared = std::numeric_limits<double>::quiet_NaN();
    double selfA = std::numeric_limits<double>::quiet_NaN();
    double selfB = std::numeric_limits<double>::quiet_NaN();
    double aRatio = std::numeric_limits<double>::quiet_NaN();
    double bRatio = std::numeric_limits<double>::quiet_NaN();
};

struct SelfParams {
    double a_cm_inv = 0.0;
    double b_V_per_cm = 0.0;
};

struct RegionSpec {
    std::string name;
    double minField = 0.0;
    double maxField = 0.0;
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

std::vector<FitSample> buildSamples(const std::vector<NodeRecord>& records,
                                    const std::vector<Element>& elements,
                                    bool electron)
{
    const auto neighbors = buildNeighbors(elements);
    std::map<double, std::vector<NodeRecord*>> byBias;
    for (const NodeRecord& record : records)
        byBias[record.bias].push_back(const_cast<NodeRecord*>(&record));

    std::vector<FitSample> samples;
    for (auto& [bias, nodes] : byBias) {
        std::unordered_map<int, NodeRecord*> nodesById;
        std::unordered_map<int, double> qfValues;
        for (NodeRecord* node : nodes) {
            nodesById[node->nodeId] = node;
            qfValues[node->nodeId] = electron ? *node->sentElectronQf : *node->sentHoleQf;
        }

        for (const NodeRecord* node : nodes) {
            FitSample sample;
            sample.biasText = node->biasText;
            sample.bias = bias;
            sample.field_V_per_cm = reconstructedGradient_V_per_cm(
                node->nodeId, nodesById, neighbors, qfValues);
            sample.alpha_cm_inv = electron
                ? *node->sentElectronAlpha_cm_inv
                : *node->sentHoleAlpha_cm_inv;
            samples.push_back(sample);
        }
    }
    return samples;
}

SelfParams selfParams(bool electron, const std::string& region)
{
    const ImpactIonizationModelConfig config = impactIonizationModelConfig("van_overstraeten");
    const bool high = region == "high";
    if (electron) {
        return {
            (high ? config.electronAHigh : config.electronALow) / 100.0,
            (high ? config.electronBHigh : config.electronBLow) / 100.0,
        };
    }
    return {
        (high ? config.holeAHigh : config.holeALow) / 100.0,
        (high ? config.holeBHigh : config.holeBLow) / 100.0,
    };
}

double ratio(double value, double reference)
{
    return reference != 0.0 && std::isfinite(reference) ? value / reference
                                                        : std::numeric_limits<double>::quiet_NaN();
}

bool inRegion(const FitSample& sample, const RegionSpec& region, double alphaFloor)
{
    return std::isfinite(sample.field_V_per_cm) &&
           std::isfinite(sample.alpha_cm_inv) &&
           sample.field_V_per_cm >= region.minField &&
           sample.field_V_per_cm <= region.maxField &&
           sample.alpha_cm_inv > alphaFloor;
}

FitResult fitSamples(const std::vector<FitSample>& allSamples,
                     bool electron,
                     const std::string& fitScope,
                     const std::string& biasText,
                     const RegionSpec& region,
                     double alphaFloor)
{
    FitResult result;
    result.carrier = electron ? "electron" : "hole";
    result.fitScope = fitScope;
    result.biasText = biasText;
    result.fieldRegion = region.name;
    const SelfParams self = selfParams(electron, region.name);
    result.selfA = self.a_cm_inv;
    result.selfB = self.b_V_per_cm;

    std::vector<FitSample> samples;
    std::copy_if(allSamples.begin(), allSamples.end(), std::back_inserter(samples),
                 [&](const FitSample& sample) {
                     const bool biasMatches = fitScope == "all_biases" ||
                         std::abs(sample.bias - parseDouble(biasText, "bias_V")) < 1.0e-9;
                     return biasMatches && inRegion(sample, region, alphaFloor);
                 });
    result.sampleCount = samples.size();
    if (samples.empty())
        return result;

    result.fMin = samples.front().field_V_per_cm;
    result.fMax = samples.front().field_V_per_cm;
    result.alphaMin = samples.front().alpha_cm_inv;
    result.alphaMax = samples.front().alpha_cm_inv;
    for (const FitSample& sample : samples) {
        result.fMin = std::min(result.fMin, sample.field_V_per_cm);
        result.fMax = std::max(result.fMax, sample.field_V_per_cm);
        result.alphaMin = std::min(result.alphaMin, sample.alpha_cm_inv);
        result.alphaMax = std::max(result.alphaMax, sample.alpha_cm_inv);
    }

    if (samples.size() < 2)
        return result;

    double sumX = 0.0;
    double sumY = 0.0;
    double sumXX = 0.0;
    double sumXY = 0.0;
    for (const FitSample& sample : samples) {
        const double x = 1.0 / sample.field_V_per_cm;
        const double y = std::log(sample.alpha_cm_inv);
        sumX += x;
        sumY += y;
        sumXX += x * x;
        sumXY += x * y;
    }
    const double n = static_cast<double>(samples.size());
    const double denom = n * sumXX - sumX * sumX;
    if (std::abs(denom) <= 1.0e-60)
        return result;

    const double slope = (n * sumXY - sumX * sumY) / denom;
    result.c = (sumY - slope * sumX) / n;
    result.d = -slope;
    result.aEff = std::exp(result.c);
    result.bEff = result.d;

    const double meanY = sumY / n;
    double sse = 0.0;
    double sst = 0.0;
    for (const FitSample& sample : samples) {
        const double x = 1.0 / sample.field_V_per_cm;
        const double y = std::log(sample.alpha_cm_inv);
        const double predicted = result.c + slope * x;
        sse += (y - predicted) * (y - predicted);
        sst += (y - meanY) * (y - meanY);
    }
    result.rSquared = sst > 0.0 ? 1.0 - sse / sst : (sse <= 1.0e-24 ? 1.0 : 0.0);
    result.aRatio = ratio(result.selfA, result.aEff);
    result.bRatio = ratio(result.selfB, result.bEff);
    return result;
}

std::vector<FitResult> buildFitResults(const std::vector<NodeRecord>& records,
                                       const std::vector<Element>& elements,
                                       double alphaFloor)
{
    const std::vector<RegionSpec> regions = {
        {"low_mid", 5.0e4, 2.5e5},
        {"high", 2.5e5, 5.0e5},
    };
    std::map<double, std::string> biasTexts;
    for (const NodeRecord& record : records)
        biasTexts.emplace(record.bias, record.biasText);

    std::vector<FitResult> results;
    for (const bool electron : {true, false}) {
        const std::vector<FitSample> samples = buildSamples(records, elements, electron);
        for (const RegionSpec& region : regions) {
            results.push_back(fitSamples(
                samples, electron, "all_biases", "all", region, alphaFloor));
            for (const auto& [_, biasText] : biasTexts) {
                results.push_back(fitSamples(
                    samples, electron, "per_bias", biasText, region, alphaFloor));
            }
        }
    }
    return results;
}

std::string fmt(double value)
{
    if (!std::isfinite(value))
        return "";
    std::ostringstream os;
    os << std::setprecision(12) << value;
    return os.str();
}

void writeFitCsv(const std::filesystem::path& path, const std::vector<FitResult>& results)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to write output CSV: " + path.string());

    out << "carrier,fit_scope,bias_V,field_region,sample_count,F_min_V_per_cm,"
           "F_max_V_per_cm,alpha_min_cm_inv,alpha_max_cm_inv,C,D_V_per_cm,"
           "A_eff_cm_inv,B_eff_V_per_cm,r_squared,self_A_cm_inv,self_B_V_per_cm,"
           "A_ratio_self_over_fit,B_ratio_self_over_fit\n";
    for (const FitResult& row : results) {
        out << row.carrier << ',' << row.fitScope << ',' << row.biasText << ','
            << row.fieldRegion << ',' << row.sampleCount << ','
            << fmt(row.fMin) << ',' << fmt(row.fMax) << ','
            << fmt(row.alphaMin) << ',' << fmt(row.alphaMax) << ','
            << fmt(row.c) << ',' << fmt(row.d) << ',' << fmt(row.aEff) << ','
            << fmt(row.bEff) << ',' << fmt(row.rSquared) << ','
            << fmt(row.selfA) << ',' << fmt(row.selfB) << ','
            << fmt(row.aRatio) << ',' << fmt(row.bRatio) << '\n';
    }
}

std::optional<FitResult> findAllBiasResult(const std::vector<FitResult>& results,
                                           const std::string& carrier,
                                           const std::string& region)
{
    for (const FitResult& row : results) {
        if (row.carrier == carrier && row.fitScope == "all_biases" &&
            row.fieldRegion == region && std::isfinite(row.aRatio) &&
            std::isfinite(row.bRatio)) {
            return row;
        }
    }
    return std::nullopt;
}

bool ratioNearOne(double value)
{
    return std::isfinite(value) && value >= 0.5 && value <= 2.0;
}

std::string consistencyText(const std::vector<FitResult>& results, const std::string& region)
{
    const auto electron = findAllBiasResult(results, "electron", region);
    const auto hole = findAllBiasResult(results, "hole", region);
    if (!electron && !hole)
        return u8s(u8"\u6837\u672c\u4e0d\u8db3\uff0c\u65e0\u6cd5\u5224\u65ad");

    const bool electronOk = electron &&
        ratioNearOne(electron->aRatio) && ratioNearOne(electron->bRatio);
    const bool holeOk = hole &&
        ratioNearOne(hole->aRatio) && ratioNearOne(hole->bRatio);
    if (electronOk && holeOk)
        return u8s(u8"\u4e0e\u81ea\u7814\u53c2\u6570\u57fa\u672c\u4e00\u81f4");
    if (electronOk || holeOk)
        return u8s(u8"\u5355\u4e00\u8f7d\u6d41\u5b50\u63a5\u8fd1\uff0c\u53e6\u4e00\u4fa7\u6709\u660e\u663e\u504f\u5dee");
    return u8s(u8"\u4e0e\u81ea\u7814\u53c2\u6570\u4e0d\u4e00\u81f4");
}

double deviationScore(const FitResult& row)
{
    double score = 0.0;
    if (std::isfinite(row.aRatio) && row.aRatio > 0.0)
        score += std::abs(std::log10(row.aRatio));
    if (std::isfinite(row.bRatio) && row.bRatio > 0.0)
        score += std::abs(std::log10(row.bRatio));
    return score;
}

std::string largerCarrierDeviation(const std::vector<FitResult>& results)
{
    std::map<std::string, std::vector<double>> scores;
    for (const FitResult& row : results) {
        if (row.fitScope != "all_biases")
            continue;
        if (std::isfinite(row.aRatio) && std::isfinite(row.bRatio))
            scores[row.carrier].push_back(deviationScore(row));
    }
    const auto meanScore = [&](const std::string& carrier) {
        const auto it = scores.find(carrier);
        if (it == scores.end() || it->second.empty())
            return std::numeric_limits<double>::quiet_NaN();
        double sum = 0.0;
        for (double value : it->second)
            sum += value;
        return sum / static_cast<double>(it->second.size());
    };
    const double electron = meanScore("electron");
    const double hole = meanScore("hole");
    if (!std::isfinite(electron) && !std::isfinite(hole))
        return u8s(u8"\u6837\u672c\u4e0d\u8db3\uff0c\u65e0\u6cd5\u5224\u65ad");
    if (std::isfinite(electron) && (!std::isfinite(hole) || electron > hole * 1.2))
        return "electron";
    if (std::isfinite(hole) && (!std::isfinite(electron) || hole > electron * 1.2))
        return "hole";
    return u8s(u8"\u7535\u5b50\u548c\u7a7a\u7a74\u63a5\u8fd1");
}

std::string segmentSwitchText(const std::vector<FitResult>& results)
{
    bool obvious = false;
    for (const std::string& carrier : {"electron", "hole"}) {
        const auto low = findAllBiasResult(results, carrier, "low_mid");
        const auto high = findAllBiasResult(results, carrier, "high");
        if (!low || !high)
            continue;
        const double aJump = ratio(high->aEff, low->aEff);
        const double bJump = ratio(high->bEff, low->bEff);
        if ((std::isfinite(aJump) && (aJump > 2.0 || aJump < 0.5)) ||
            (std::isfinite(bJump) && (bJump > 1.5 || bJump < 0.67))) {
            obvious = true;
        }
    }
    if (obvious) {
        return u8s(u8"\u5b58\u5728\u660e\u663e\u5206\u6bb5\u8ff9\u8c61\uff1b"
                   "\u672c\u5de5\u5177\u7528 2.5e5 V/cm \u5206\u533a\uff0c"
                   "\u81ea\u7814\u9ed8\u8ba4 switchField=4.0e5 V/cm");
    }
    return u8s(u8"\u672a\u770b\u5230\u5f3a\u5206\u6bb5\u8df3\u53d8\uff1b"
               "\u4ecd\u9700\u7ed3\u5408\u9010\u8282\u70b9 alpha ratio \u5224\u65ad");
}

std::string likelyCauseText(const std::vector<FitResult>& results)
{
    int aTooSmall = 0;
    int bTooLarge = 0;
    int bTooSmall = 0;
    int segmentMismatch = 0;
    for (const FitResult& row : results) {
        if (row.fitScope != "all_biases" || !std::isfinite(row.aRatio) ||
            !std::isfinite(row.bRatio)) {
            continue;
        }
        if (row.aRatio < 0.5)
            ++aTooSmall;
        if (row.bRatio > 1.2)
            ++bTooLarge;
        if (row.bRatio < 0.8)
            ++bTooSmall;
        const bool high = row.fieldRegion == "high";
        if ((high && !ratioNearOne(row.aRatio)) || (high && !ratioNearOne(row.bRatio)))
            ++segmentMismatch;
    }

    if (bTooLarge >= aTooSmall && bTooLarge > 0)
        return u8s(u8"B \u592a\u5927");
    if (aTooSmall > 0)
        return u8s(u8"A \u592a\u5c0f");
    if (segmentMismatch > 0 || bTooSmall > 0)
        return u8s(u8"\u5206\u6bb5\u9519\u8bef");
    return u8s(u8"\u53c2\u6570\u62df\u5408\u672a\u663e\u793a\u660e\u786e A/B \u504f\u79bb\uff0c"
               "\u82e5\u80a9\u90e8 alpha \u4ecd\u504f\u5c0f\uff0c\u66f4\u50cf cutoff \u592a\u5f3a"
               "\u6216\u8282\u70b9\u6620\u5c04\u5bfc\u81f4");
}

void writeSummary(const std::filesystem::path& path,
                  const std::vector<FitResult>& results,
                  double alphaFloor)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("failed to write summary: " + path.string());

    out << "# Sentaurus VanOverstraeten Effective Fit\n\n";
    out << "- F units: V/cm\n";
    out << "- alpha units: 1/cm\n";
    out << "- Fit model: ln(alpha) = C - D/F, A_eff = exp(C), B_eff = D\n";
    out << "- Gradient recovery: 1/d weighted least squares on Sentaurus qF node values and elements.csv connectivity\n";
    out << "- alpha filter: alpha > " << fmt(alphaFloor) << " 1/cm\n";
    out << "- " << u8s(u8"\u8fd9\u662f\u8282\u70b9\u8f93\u51fa\u7b49\u6548\u62df\u5408\uff0c"
                       "\u4e0d\u4e00\u5b9a\u7b49\u4e8e Sentaurus \u5185\u90e8\u53c2\u6570\uff0c"
                       "\u4f46\u8db3\u4ee5\u5224\u65ad\u81ea\u7814\u53c2\u6570\u662f\u5426\u504f\u79bb")
        << "\n\n";

    out << "## All-Bias Fits\n\n";
    out << "| carrier | region | samples | A_eff (1/cm) | B_eff (V/cm) | R2 | "
           "self_A/fit_A | self_B/fit_B |\n";
    out << "|---|---|---:|---:|---:|---:|---:|---:|\n";
    for (const FitResult& row : results) {
        if (row.fitScope != "all_biases")
            continue;
        out << "| " << row.carrier << " | " << row.fieldRegion << " | "
            << row.sampleCount << " | " << fmt(row.aEff) << " | "
            << fmt(row.bEff) << " | " << fmt(row.rSquared) << " | "
            << fmt(row.aRatio) << " | " << fmt(row.bRatio) << " |\n";
    }

    out << "\n## Per-Bias Fits\n\n";
    out << "| carrier | bias | region | samples | A_eff (1/cm) | B_eff (V/cm) | R2 |\n";
    out << "|---|---:|---|---:|---:|---:|---:|\n";
    for (const FitResult& row : results) {
        if (row.fitScope != "per_bias")
            continue;
        out << "| " << row.carrier << " | " << row.biasText << " | "
            << row.fieldRegion << " | " << row.sampleCount << " | "
            << fmt(row.aEff) << " | " << fmt(row.bEff) << " | "
            << fmt(row.rSquared) << " |\n";
    }

    out << "\n## Answers\n\n";
    out << "1. " << u8s(u8"Sentaurus \u7684\u4f4e\u4e2d\u573a\u7b49\u6548 A/B: ")
        << consistencyText(results, "low_mid") << "\n";
    out << "2. " << u8s(u8"Sentaurus \u7684\u9ad8\u573a\u7b49\u6548 A/B: ")
        << consistencyText(results, "high") << "\n";
    out << "3. " << u8s(u8"\u504f\u5dee\u66f4\u5927\u7684\u8f7d\u6d41\u5b50: ")
        << largerCarrierDeviation(results) << "\n";
    out << "4. " << u8s(u8"\u5206\u6bb5\u5207\u6362\u70b9: ")
        << segmentSwitchText(results) << "\n";
    out << "5. " << u8s(u8"\u5f53\u524d\u81ea\u7814 alpha \u504f\u5c0f\u66f4\u50cf: ")
        << likelyCauseText(results) << "\n";
}

void usage()
{
    std::cerr
        << "Usage: sentaurus_vanoverstraeten_fit --input-csv FILE --elements-csv FILE "
           "[--output-csv FILE] [--summary-md FILE] [--alpha-floor VALUE]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::filesystem::path inputCsv;
        std::filesystem::path elementsCsv;
        std::filesystem::path outputCsv =
            "build/diagnostics/sentaurus_vanoverstraeten_fit.csv";
        std::filesystem::path summaryMd =
            "build/diagnostics/sentaurus_vanoverstraeten_fit_summary.md";
        double alphaFloor = 1.0e-20;

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
            } else if (arg == "--alpha-floor") {
                alphaFloor = parseDouble(requireValue(args, i, arg), arg);
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
        if (alphaFloor < 0.0 || !std::isfinite(alphaFloor))
            throw std::runtime_error("--alpha-floor must be finite and non-negative");

        const std::vector<NodeRecord> records = loadNodeRecords(inputCsv);
        const std::vector<Element> elements = loadElements(elementsCsv);
        const std::vector<FitResult> results = buildFitResults(records, elements, alphaFloor);
        writeFitCsv(outputCsv, results);
        writeSummary(summaryMd, results, alphaFloor);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "sentaurus_vanoverstraeten_fit: " << ex.what() << '\n';
        return 1;
    }
}
