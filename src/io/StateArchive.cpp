#include "vela/io/StateArchive.h"
#include <highfive/highfive.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace vela {
namespace {
constexpr const char* schema = "vela.state/2";
constexpr std::uint64_t maxBytes = 1024ULL * 1024 * 1024;
const std::map<std::string, std::string> units{
    {"psi", "V"}, {"phin", "V"}, {"phip", "V"},
    {"electrons_m3", "m^-3"}, {"holes_m3", "m^-3"},
    {"electron_quantum_potential_V", "V"}, {"electron_quantum_potential_like_V", "V"},
    {"electron_qf_increment_V", "V"}, {"hole_qf_increment_V", "V"},
    {"electron_qf_reference_V", "V"}, {"hole_qf_reference_V", "V"},
    {"temperature_K", "K"}};
[[noreturn]] void bad(const std::string& message) {
    throw std::runtime_error("HDF5 state: " + message);
}
bool hashValid(const std::string& hash) {
    return hash.size() == 64 && std::all_of(hash.begin(), hash.end(), [](char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    });
}
void checkPath(const std::filesystem::path& path) {
    if (path.extension() != ".h5") bad("state path must end in .h5");
}
void dimensions(std::uint64_t count, std::uint64_t fields) {
    if (count == 0 || fields < 4 || fields > units.size() || count > maxBytes / (8 * fields))
        bad("invalid or oversized dimensions");
}
std::string textAttribute(const HighFive::File& file, const char* name) {
    const auto a = file.getAttribute(name);
    if (a.getSpace().getElementCount() != 1) bad(std::string("non-scalar attribute ") + name);
    return a.read<std::string>();
}
} // namespace

void validateStateArchive(const StateArchive& s) {
    dimensions(s.nodeCount, s.fields.size());
    if (!s.metadata.is_object() || s.metadata.dump().size() > 65536) bad("invalid metadata");
    const auto mode = s.metadata.at("mode").get<std::string>();
    if (mode != "dd" && mode != "electrothermal") bad("unknown equation mode");
    if (!hashValid(s.metadata.at("mesh_sha256").get<std::string>())) bad("invalid mesh identity");
    if (!s.metadata.at("potential_origin_V").is_number() ||
        !std::isfinite(s.metadata.at("potential_origin_V").get<double>())) bad("invalid potential origin");
    for (const auto* name : {"psi", "phin", "phip"})
        if (!s.fields.contains(name)) bad(std::string("missing field ") + name);
    if (mode == "dd" && (!s.fields.contains("electrons_m3") || !s.fields.contains("holes_m3")))
        bad("DD densities required");
    if (s.fields.contains("electrons_m3") != s.fields.contains("holes_m3")) bad("partial density pair");
    if ((mode == "electrothermal") != s.fields.contains("temperature_K")) bad("temperature/mode mismatch");
    if (s.fields.contains("electron_quantum_potential_like_V") &&
        !s.fields.contains("electron_quantum_potential_V")) bad("partial quantum state");
    std::size_t split = 0;
    for (const auto* name : {"electron_qf_increment_V", "hole_qf_increment_V",
                            "electron_qf_reference_V", "hole_qf_reference_V"})
        split += s.fields.contains(name);
    if (split != 0 && split != 4) bad("partial split quasi-Fermi state");
    for (const auto& [name, values] : s.fields) {
        if (!units.contains(name) || values.size() != s.nodeCount) bad("unknown or incomplete field " + name);
        if (!std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }))
            bad("non-finite field " + name);
    }
    if (split) for (std::size_t i = 0; i < s.nodeCount; ++i) {
        for (const auto& carrier : {std::string("electron"), std::string("hole")}) {
            const double value = s.fields.at(carrier == "electron" ? "phin" : "phip")[i];
            const double combined = s.fields.at(carrier + "_qf_reference_V")[i] +
                                    s.fields.at(carrier + "_qf_increment_V")[i];
            if (!std::isfinite(combined) || std::abs(value - combined) >
                32 * std::numeric_limits<double>::epsilon() * std::max({1., std::abs(value), std::abs(combined)}))
                bad("inconsistent split quasi-Fermi state");
        }
    }
}

void writeStateArchive(const std::filesystem::path& path, const StateArchive& state) {
    checkPath(path);
    validateStateArchive(state);
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    static std::atomic<std::uint64_t> sequence{0};
    auto temporary = path;
    temporary += ".tmp." + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
                 "." + std::to_string(sequence.fetch_add(1));
    try {
        {
            HighFive::File f(temporary.string(), HighFive::File::Create | HighFive::File::Excl);
            f.createAttribute("schema", std::string(schema));
            f.createAttribute("node_count", state.nodeCount);
            f.createAttribute("metadata_json", state.metadata.dump());
            auto group = f.createGroup("fields");
            for (const auto& [name, values] : state.fields) {
                HighFive::AtomicType<double> diskType;
                if (H5Tset_order(diskType.getId(), H5T_ORDER_LE) < 0 ||
                    H5Tequal(diskType.getId(), H5T_IEEE_F64LE) <= 0) bad("binary64 representation unavailable");
                auto d = group.createDataSet(name, HighFive::DataSpace({values.size()}), diskType);
                d.createAttribute("unit", units.at(name));
                d.write_raw(values.data());
            }
            f.flush();
        }
#ifdef _WIN32
        if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
            bad("cannot replace checkpoint, Windows error " + std::to_string(GetLastError()));
#else
        std::filesystem::rename(temporary, path);
#endif
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

StateArchive readStateArchive(const std::filesystem::path& path, std::uint64_t expected,
                              const std::string& mesh) {
    checkPath(path);
    if (!hashValid(mesh)) bad("expected mesh identity is required");
    if (std::filesystem::file_size(path) > maxBytes + 1024 * 1024) bad("file exceeds size limit");
    HighFive::File f(path.string(), HighFive::File::ReadOnly);
    if (textAttribute(f, "schema") != schema) bad("unsupported schema");
    StateArchive s;
    const auto countAttribute = f.getAttribute("node_count");
    const auto countType = countAttribute.getDataType();
    if (countAttribute.getSpace().getElementCount() != 1 || H5Tget_class(countType.getId()) != H5T_INTEGER ||
        H5Tget_size(countType.getId()) != 8 || H5Tget_sign(countType.getId()) != H5T_SGN_NONE)
        bad("node count must be scalar uint64");
    s.nodeCount = countAttribute.read<std::uint64_t>();
    if (s.nodeCount != expected) bad("node count mismatch");
    const auto rawMetadata = textAttribute(f, "metadata_json");
    if (rawMetadata.size() > 65536) bad("metadata too large");
    s.metadata = nlohmann::json::parse(rawMetadata);
    if (s.metadata.at("mesh_sha256") != mesh) bad("mesh identity mismatch");
    auto group = f.getGroup("fields");
    const auto names = group.listObjectNames();
    dimensions(s.nodeCount, names.size());
    for (const auto& name : names) {
        if (!units.contains(name)) bad("unknown field " + name);
        auto d = group.getDataSet(name);
        const auto type = d.getDataType();
        if (d.getDimensions() != std::vector<std::size_t>{static_cast<std::size_t>(s.nodeCount)} ||
            H5Tequal(type.getId(), H5T_IEEE_F64LE) <= 0) bad("field type or dimensions mismatch");
        if (d.getAttribute("unit").read<std::string>() != units.at(name)) bad("field unit mismatch");
        auto& values = s.fields[name];
        values.resize(s.nodeCount);
        d.read_raw(values.data());
    }
    validateStateArchive(s);
    return s;
}
} // namespace vela
