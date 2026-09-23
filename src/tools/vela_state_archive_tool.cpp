// Migration/inspection utility, not a production legacy restart fallback.
#include "vela/io/StateArchive.h"
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc != 4 && argc != 5) throw std::runtime_error("encode JSON H5 | decode H5 N MESH_SHA256");
        if (std::string(argv[1]) == "encode" && argc == 4) {
            std::ifstream stream(argv[2]);
            if (!stream) throw std::runtime_error("Cannot read migration input");
            const auto j = nlohmann::json::parse(stream);
            vela::StateArchive s;
            s.metadata = j.at("metadata");
            s.fields = j.at("fields").get<decltype(s.fields)>();
            s.nodeCount = s.fields.at("psi").size();
            vela::writeStateArchive(argv[3], s);
        } else if (std::string(argv[1]) == "decode" && argc == 5) {
            const auto s = vela::readStateArchive(argv[2], std::stoull(argv[3]), argv[4]);
            std::cout << nlohmann::json({{"fields", s.fields}, {"metadata", s.metadata}}).dump() << '\n';
        } else throw std::runtime_error("Invalid migration utility command");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
