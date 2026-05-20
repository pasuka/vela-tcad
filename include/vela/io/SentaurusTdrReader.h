#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace vela {

enum class SentaurusTdrRegionType {
    Material = 0,
    Contact = 1,
    Interface = 2,
    Other = 99,
};

struct SentaurusTdrVertex {
    double x = 0.0;
    double y = 0.0;
};

struct SentaurusTdrRegion {
    int index = -1;
    std::string name;
    std::string material;
    SentaurusTdrRegionType type = SentaurusTdrRegionType::Other;
    std::vector<std::array<std::size_t, 3>> triangles;
    std::vector<std::array<std::size_t, 2>> edges;
};

struct SentaurusTdrField {
    int index = -1;
    std::string name;
    int region_index = -1;
    std::string unit;
    std::size_t value_count = 0;
    std::size_t component_count = 1;
    std::vector<double> values;
};

struct SentaurusTdrInventory {
    std::vector<SentaurusTdrVertex> vertices;
    std::vector<SentaurusTdrRegion> regions;
    std::vector<SentaurusTdrField> fields;

    const SentaurusTdrField* findField(const std::string& name, int regionIndex) const;
};

class SentaurusTdrReader {
public:
    SentaurusTdrInventory readInventory(const std::string& filename) const;
    void exportNeutral(const std::string& filename, const std::string& outputDirectory) const;
};

} // namespace vela
