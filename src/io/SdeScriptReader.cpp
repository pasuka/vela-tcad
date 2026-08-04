#include "vela/io/SdeScriptReader.h"

#include <cctype>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace vela {
namespace {

struct ConstantProfileDefinition {
    std::string species;
    Real concentration_cm3 = 0.0;
};

struct GaussianProfileDefinition {
    std::string species;
    Real peak_cm3 = 0.0;
    Real background_cm3 = 0.0;
    Real depth_um = 0.0;
};

struct RefinementSizeDefinition {
    Real lateralSizeUm = 0.0;
};

bool isCommentOrEmpty(const std::string& line)
{
    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            continue;
        }
        return c == ';' || c == '#';
    }
    return true;
}

bool speciesIsAcceptor(const std::string& species)
{
    return species.find("Boron") != std::string::npos ||
        species.find("Acceptor") != std::string::npos;
}

std::runtime_error parseError(const std::string& sourceName,
                              std::size_t lineNumber,
                              const std::string& message)
{
    return std::runtime_error(
        "SDE parse error in '" + sourceName + "' line " + std::to_string(lineNumber) +
        ": " + message);
}

RegionIr& ensureRegion(DeviceIr2D& ir, const std::string& regionName, const std::string& material)
{
    if (RegionIr* region = ir.findRegion(regionName)) {
        if (!material.empty()) {
            region->material = material;
        }
        return *region;
    }
    ir.regions.push_back(RegionIr{regionName, material, {}});
    return ir.regions.back();
}

} // namespace

DeviceIr2D SdeScriptReader::parseFile(const std::string& filename) const
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open SDE script: " + filename);
    }
    std::ostringstream text;
    text << in.rdbuf();
    return parseText(text.str(), filename);
}

DeviceIr2D SdeScriptReader::parseText(const std::string& text, const std::string& sourceName) const
{
    // Supported subset is intentionally strict for phase-1.
    const std::regex createRectangle(
        R"regex(\(sdegeo:create-rectangle\s+\(position\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+[-\d.eE+]+\)\s+\(position\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+[-\d.eE+]+\)\s+"([^"]+)"\s+"([^"]+)"\s*\))regex");
    const std::regex createPolygon(
        R"regex(\(sdegeo:create-polygon\s+\(list\s+(.+)\)\s+"([^"]+)"\s+"([^"]+)"\s*\))regex");
    const std::regex position2D(
        R"regex(\(position\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+[-\d.eE+]+\))regex");
    const std::regex defineConstantProfile(
        R"regex(\(sdedr:define-constant-profile\s+"([^"]+)"\s+"([^"]+)"\s+([-\d.eE+]+)\s*\))regex");
    const std::regex defineConstantProfileRegion(
        R"regex(\(sdedr:define-constant-profile-region\s+"([^"]+)"\s+"([^"]+)"\s+"([^"]+)"\s*\))regex");
    const std::regex defineGaussianProfile(
        R"regex(\(sdedr:define-gaussian-profile\s+"([^"]+)"\s+"([^"]+)".*"PeakVal"\s+([-\d.eE+]+).*"ValueAtDepth"\s+([-\d.eE+]+).*"Depth"\s+([-\d.eE+]+).*\))regex");
    const std::regex defineGaussianProfileRegion(
        R"regex(\(sdedr:define-analytical-profile-region\s+"([^"]+)"\s+"([^"]+)"\s+"([^"]+)"\s*\))regex");
    const std::regex defineRefinementSize(
        R"regex(\(sdedr:define-refinement-size\s+"([^"]+)"\s+([-\d.eE+]+).*\))regex");
    const std::regex defineRefinementRegion(
        R"regex(\(sdedr:define-refinement-region\s+"([^"]+)"\s+"([^"]+)"\s+"([^"]+)"\s*\))regex");
    const std::regex defineRefinementFunction(
        R"regex(\(sdedr:define-refinement-function\s+"([^"]+)"\s+"([^"]+)"\s+"([^"]+)"\s+([-\d.eE+]+).*\))regex");

    std::unordered_map<std::string, ConstantProfileDefinition> constantProfiles;
    std::unordered_map<std::string, GaussianProfileDefinition> gaussianProfiles;
    std::unordered_map<std::string, RefinementSizeDefinition> refinementSizes;
    DeviceIr2D ir;
    std::istringstream lines(text);
    std::string line;
    std::size_t lineNumber = 0;
    std::size_t primitiveCounter = 0;

    while (std::getline(lines, line)) {
        ++lineNumber;
        if (isCommentOrEmpty(line)) {
            continue;
        }

        std::smatch match;
        if (std::regex_search(line, match, createRectangle)) {
            const Real x0 = std::stod(match[1].str());
            const Real y0 = std::stod(match[2].str());
            const Real x1 = std::stod(match[3].str());
            const Real y1 = std::stod(match[4].str());
            const std::string material = match[5].str();
            const std::string regionName = match[6].str();

            GeometryPrimitiveIr primitive;
            primitive.kind = GeometryPrimitiveKind::Rectangle;
            primitive.name = "rect_" + std::to_string(primitiveCounter++);
            primitive.region = regionName;
            primitive.material = material;
            primitive.points = {
                Point2D{x0, y0},
                Point2D{x1, y1},
            };
            ir.geometry.push_back(primitive);
            RegionIr& region = ensureRegion(ir, regionName, material);
            region.primitiveIndices.push_back(ir.geometry.size() - 1);
            continue;
        }

        if (std::regex_search(line, match, createPolygon)) {
            const std::string pointsBody = match[1].str();
            const std::string material = match[2].str();
            const std::string regionName = match[3].str();

            GeometryPrimitiveIr primitive;
            primitive.kind = GeometryPrimitiveKind::Polygon;
            primitive.name = "poly_" + std::to_string(primitiveCounter++);
            primitive.region = regionName;
            primitive.material = material;

            for (auto pointIt = std::sregex_iterator(pointsBody.begin(), pointsBody.end(), position2D);
                 pointIt != std::sregex_iterator(); ++pointIt) {
                primitive.points.push_back(Point2D{
                    std::stod((*pointIt)[1].str()),
                    std::stod((*pointIt)[2].str()),
                });
            }
            if (primitive.points.size() < 3) {
                throw parseError(sourceName, lineNumber,
                    "polygon must contain at least 3 vertices");
            }

            ir.geometry.push_back(primitive);
            RegionIr& region = ensureRegion(ir, regionName, material);
            region.primitiveIndices.push_back(ir.geometry.size() - 1);
            continue;
        }

        if (std::regex_search(line, match, defineConstantProfile)) {
            constantProfiles[match[1].str()] = ConstantProfileDefinition{
                match[2].str(),
                std::stod(match[3].str()),
            };
            continue;
        }

        if (std::regex_search(line, match, defineConstantProfileRegion)) {
            const std::string profileName = match[2].str();
            const std::string regionName = match[3].str();
            const auto profileIt = constantProfiles.find(profileName);
            if (profileIt == constantProfiles.end()) {
                throw parseError(sourceName, lineNumber,
                    "constant profile '" + profileName + "' is not defined");
            }
            ensureRegion(ir, regionName, "");

            DopingProfileIr profile;
            profile.name = match[1].str();
            profile.kind = DopingProfileKind::Constant;
            profile.targetRegion = regionName;
            if (speciesIsAcceptor(profileIt->second.species)) {
                profile.acceptors_cm3 = profileIt->second.concentration_cm3;
            } else {
                profile.donors_cm3 = profileIt->second.concentration_cm3;
            }
            ir.dopingProfiles.push_back(profile);
            continue;
        }

        if (std::regex_search(line, match, defineGaussianProfile)) {
            gaussianProfiles[match[1].str()] = GaussianProfileDefinition{
                match[2].str(),
                std::stod(match[3].str()),
                std::stod(match[4].str()),
                std::stod(match[5].str()),
            };
            continue;
        }

        if (std::regex_search(line, match, defineGaussianProfileRegion)) {
            const std::string profileName = match[2].str();
            const std::string regionName = match[3].str();
            const auto profileIt = gaussianProfiles.find(profileName);
            if (profileIt == gaussianProfiles.end()) {
                throw parseError(sourceName, lineNumber,
                    "gaussian profile '" + profileName + "' is not defined");
            }
            ensureRegion(ir, regionName, "");

            DopingProfileIr profile;
            profile.name = match[1].str();
            profile.kind = DopingProfileKind::Gaussian;
            profile.targetRegion = regionName;
            profile.gaussianPeak_cm3 = profileIt->second.peak_cm3;
            profile.gaussianBackground_cm3 = profileIt->second.background_cm3;
            profile.gaussianSigmaXUm = profileIt->second.depth_um;
            profile.gaussianSigmaYUm = profileIt->second.depth_um;
            profile.gaussianActsOnDonors = !speciesIsAcceptor(profileIt->second.species);
            ir.dopingProfiles.push_back(profile);
            continue;
        }

        if (std::regex_search(line, match, defineRefinementSize)) {
            refinementSizes[match[1].str()] = RefinementSizeDefinition{
                std::stod(match[2].str()),
            };
            continue;
        }

        if (std::regex_search(line, match, defineRefinementRegion)) {
            const std::string refinementName = match[2].str();
            const std::string regionName = match[3].str();
            const auto sizeIt = refinementSizes.find(refinementName);
            if (sizeIt == refinementSizes.end()) {
                throw parseError(sourceName, lineNumber,
                    "refinement size '" + refinementName + "' is not defined");
            }
            ensureRegion(ir, regionName, "");
            ir.meshControl.regionTargetSizeUm[regionName] = sizeIt->second.lateralSizeUm;
            continue;
        }

        if (std::regex_search(line, match, defineRefinementFunction)) {
            const std::string quantity = match[2].str();
            const std::string method = match[3].str();
            if (quantity != "DopingConcentration" || method != "MaxTransDiff") {
                throw parseError(sourceName, lineNumber,
                    "only DopingConcentration + MaxTransDiff refinement is supported in phase-1");
            }
            ir.meshControl.refineByDopingGradient = true;
            ir.meshControl.dopingGradientThresholdCm3PerUm = std::stod(match[4].str());
            continue;
        }

        if (line.find("(sde") != std::string::npos || line.find("(sdedr:") != std::string::npos ||
            line.find("(sdegeo:") != std::string::npos) {
            throw parseError(sourceName, lineNumber,
                "unsupported command in phase-1 subset: " + line);
        }
    }

    return ir;
}

} // namespace vela
