#include "vela/io/StateArchive.h"
#include "vela/io/DDSolutionState.h"
#include "vela/io/StateIdentity.h"
#include "vela/io/ElectrothermalState.h"
#include <catch2/catch_test_macros.hpp>
#include <bit>
#include <chrono>
#include <fstream>
#include <limits>

TEST_CASE("Thermal record keeps temperature and independent QF increments across archive references", "[state_archive]") {
    const auto root=std::filesystem::temp_directory_path()/
        ("vela_thermal_record_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    struct Cleanup {std::filesystem::path p;~Cleanup(){std::error_code e;std::filesystem::remove_all(p,e);}} cleanup{root};
    const nlohmann::json record={{"potential_origin_V",28.},
        {"state_interleaved",{.125,28.,0.,327.25}},
        {"referenced_state_interleaved",{.125,1e-18,-0.,327.25}},
        {"electron_qf_reference_V",{28.}},{"hole_qf_reference_V",{0.}},
        {"temperature_K",{327.25}},{"newton_updates",7}};
    const std::string sha(64,'a');
    const nlohmann::json meta={{"mode","electrothermal"},{"mesh_sha256",sha},{"potential_origin_V",28.}};
    const auto path=root/"output.json";
    const auto packed=vela::packElectrothermalRecord(path,record,meta);
    REQUIRE_FALSE(packed.contains("state_interleaved"));
    REQUIRE_FALSE(packed.contains("temperature_K"));
    const auto restored=vela::unpackElectrothermalRecord(path,packed,1,sha,28.);
    REQUIRE(restored==record);
    REQUIRE(std::bit_cast<std::uint64_t>(restored.at("referenced_state_interleaved")[2].get<double>())==
            std::bit_cast<std::uint64_t>(-0.));
    REQUIRE_THROWS(vela::packElectrothermalRecord(path,record,meta));
    REQUIRE_THROWS(vela::unpackElectrothermalRecord(path,packed,1,std::string(64,'b'),28.));
    REQUIRE_THROWS(vela::unpackElectrothermalRecord(path,packed,1,sha,0.));
    REQUIRE_THROWS(vela::unpackElectrothermalRecord(path,record,1,sha,28.));
    auto corrupt=packed;corrupt["state_archive"]["sha256"]=std::string(64,'0');
    REQUIRE_THROWS(vela::unpackElectrothermalRecord(path,corrupt,1,sha,28.));
}

namespace {
struct Workspace {
    std::filesystem::path root = std::filesystem::temp_directory_path() /
        ("vela_state_contract_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    Workspace() { std::filesystem::create_directory(root); }
    ~Workspace() { std::error_code ec; std::filesystem::remove_all(root, ec); }
};
vela::StateArchive example(bool thermal = false) {
    vela::StateArchive s;
    s.nodeCount = 3;
    s.metadata = {{"mode", thermal ? "electrothermal" : "dd"}, {"mesh_sha256", std::string(64, 'a')},
                  {"potential_origin_V", 28.}, {"bias_V", 40.}};
    s.fields = {{"psi", {0., 0.25, -0.5}}, {"phin", {28., 28., 28.}},
        {"phip", {0., -0., 1e-200}}, {"electrons_m3", {1e26, 1e12, 0.}},
        {"holes_m3", {1e-120, 1e23, 1e-30}},
        {"electron_qf_increment_V", {1e-18, -1e-18, 0.}},
        {"electron_qf_reference_V", {28., 28., 28.}},
        {"hole_qf_increment_V", {0., -0., 1e-200}},
        {"hole_qf_reference_V", {0., 0., 0.}}};
    if (thermal) s.fields["temperature_K"] = {300., 327.25, 401.5};
    return s;
}
void equalBits(const vela::StateArchive& a, const vela::StateArchive& b) {
    REQUIRE(a.metadata == b.metadata);
    REQUIRE(a.nodeCount == b.nodeCount);
    REQUIRE(a.fields.size() == b.fields.size());
    for (const auto& [name, values] : a.fields)
        for (std::size_t i = 0; i < values.size(); ++i)
            REQUIRE(std::bit_cast<std::uint64_t>(values[i]) == std::bit_cast<std::uint64_t>(b.fields.at(name)[i]));
}
}

TEST_CASE("HDF5 state retains split QF low bits and temperature", "[state_archive]") {
    Workspace w;
    for (bool thermal : {false, true}) {
        auto s = example(thermal);
        s.fields["electron_quantum_potential_V"] = {0.01, -0.05, 0.};
        vela::writeStateArchive(w.root/"state.h5", s);
        equalBits(s, vela::readStateArchive(w.root/"state.h5", 3, std::string(64, 'a')));
    }
}
TEST_CASE("HDF5 state binds node layout to mesh identity", "[state_archive]") {
    Workspace w;
    vela::writeStateArchive(w.root/"state.h5", example());
    REQUIRE_THROWS(vela::readStateArchive(w.root/"state.h5", 4, std::string(64, 'a')));
    REQUIRE_THROWS(vela::readStateArchive(w.root/"state.h5", 3, std::string(64, 'b')));
    REQUIRE_THROWS(vela::readStateArchive(w.root/"state.h5", 3, ""));
    REQUIRE_THROWS(vela::writeStateArchive(w.root/"state.csv", example()));
}
TEST_CASE("Invalid HDF5 replacement preserves accepted checkpoint", "[state_archive]") {
    Workspace w;
    const auto s = example(true);
    const auto path = w.root/"state.h5";
    vela::writeStateArchive(path, s);
    auto invalid = s;
    invalid.fields.erase("temperature_K");
    REQUIRE_THROWS(vela::writeStateArchive(path, invalid));
    equalBits(s, vela::readStateArchive(path, 3, std::string(64, 'a')));
    invalid = s;
    invalid.fields["psi"][1] = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS(vela::writeStateArchive(path, invalid));
    invalid = s;
    invalid.fields["electron_qf_reference_V"][0] = 0.;
    REQUIRE_THROWS(vela::writeStateArchive(path, invalid));
    equalBits(s, vela::readStateArchive(path, 3, std::string(64, 'a')));
    REQUIRE(std::distance(std::filesystem::directory_iterator(w.root), std::filesystem::directory_iterator{}) == 1);
}
TEST_CASE("Failed HDF5 commit cleans temporary and keeps destination", "[state_archive]") {
    Workspace w;
    const auto target = w.root/"blocked.h5";
    std::filesystem::create_directory(target);
    std::ofstream(target/"marker") << "preserved";
    REQUIRE_THROWS(vela::writeStateArchive(target, example()));
    REQUIRE(std::filesystem::exists(target/"marker"));
    REQUIRE(std::distance(std::filesystem::directory_iterator(w.root), std::filesystem::directory_iterator{}) == 1);
}

TEST_CASE("DD archive preserves solver state without legacy serialization", "[state_archive]") {
    auto source = example();
    source.fields["electron_quantum_potential_V"] = {0.01, -0.05, 0.};
    const auto solution = vela::restoreDDSolution(source);
    equalBits(source, vela::archiveDDSolution(solution, source.metadata));
    auto tiny = solution;
    tiny.phip(0) = std::numeric_limits<double>::denorm_min();
    tiny.phipIncrement(0) = tiny.phip(0);
    REQUIRE(vela::archiveDDSolution(tiny, source.metadata).fields.at("phip")[0] == 0.);
    auto partial = solution;
    partial.holeQfReference.resize(1);
    REQUIRE_THROWS(vela::archiveDDSolution(partial, source.metadata));
}

TEST_CASE("State SHA256 matches known empty and multiblock vectors", "[state_archive]") {
    REQUIRE(vela::stateSha256("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    REQUIRE(vela::stateSha256("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    REQUIRE(vela::stateSha256(std::string(1000000,'a')) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}
TEST_CASE("State mesh identity distinguishes equal sized meshes and units", "[state_archive]") {
    vela::DeviceMesh a,b;
    a.addNode({0,0.,0.});a.addNode({1,1.,0.});
    b.addNode({0,1.,0.});b.addNode({1,0.,0.});
    REQUIRE(vela::stateMeshIdentity(a) != vela::stateMeshIdentity(b));
    REQUIRE(vela::stateMeshIdentity(a) != vela::stateMeshIdentity(a, vela::UnitScalingConfig{vela::UnitScalingMode::UnitScaling}));
    REQUIRE(vela::stateMeshIdentity(a) == vela::stateMeshIdentity(a));
    REQUIRE_THROWS(vela::stateMeshIdentity(a, 0.));
    REQUIRE_THROWS(vela::stateMeshIdentity(a, std::numeric_limits<double>::infinity()));
}
TEST_CASE("State provenance detects changed model bytes at the same path", "[state_archive]") {
    Workspace w;
    std::ofstream(w.root/"material.json") << "first";
    const nlohmann::json config={{"materials_file","material.json"},
        {"output_state_file","not_created.h5"},{"sweep",{{"initial_state_file","not_loaded.h5"}}}};
    const auto before=vela::stateInputProvenance(config,w.root);
    std::ofstream(w.root/"material.json") << "second";
    const auto after=vela::stateInputProvenance(config,w.root);
    REQUIRE(before.size()==1);
    REQUIRE(before.at("materials_file")!=after.at("materials_file"));
    REQUIRE_THROWS(vela::stateInputProvenance(config,w.root/"missing"));
}
