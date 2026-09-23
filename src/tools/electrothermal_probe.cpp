// Compatibility CLI for the shared four-equation simulation service.
#include "vela/simulation/ElectrothermalSimulation.h"
#include "vela/io/ElectrothermalState.h"
#include "vela/io/StateIdentity.h"
#include "vela/io/MeshReader.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
int main(int argc, char** argv) {
    try {
        if(argc!=3) throw std::invalid_argument("Usage: electrothermal_probe INPUT.json OUTPUT.json");
        if(std::filesystem::exists(argv[2])) throw std::invalid_argument("Output already exists");
        std::ifstream input(argv[1]); if(!input) throw std::invalid_argument("Cannot open input");
        nlohmann::json cfg; input>>cfg;
        const auto inputPath=std::filesystem::absolute(argv[1]);
        std::filesystem::path meshPath=cfg.at("mesh_file").get<std::string>();
        if(meshPath.is_relative()) meshPath=inputPath.parent_path()/meshPath;
        cfg["mesh_file"]=meshPath.string();
        const auto mesh=vela::JsonMeshReader{}.read(meshPath.string());
        const auto identity=vela::stateMeshIdentity(mesh,cfg.value("coordinate_to_metres",1.));
        cfg=vela::unpackElectrothermalRecord(inputPath,cfg,mesh.numNodes(),identity,cfg.value("potential_origin_V",0.));
        const auto result=vela::solveElectrothermalPoint(cfg, std::cout);
        auto packed=vela::packElectrothermalRecord(argv[2],result,{{"mode","electrothermal"},
            {"mesh_sha256",identity},{"potential_origin_V",cfg.value("potential_origin_V",0.)},
            {"input_file_sha256",vela::stateInputProvenance(cfg,inputPath.parent_path())},
            {"boundary_values_sha256",vela::stateSha256(cfg.value("boundaries",nlohmann::json::array()).dump())},
            {"source_config_sha256",vela::stateSha256(cfg.dump())}});
        packed["state_archive"]["input_file"]=inputPath.string();
        std::ofstream output(argv[2]); if(!output) throw std::runtime_error("Cannot open output");
        output<<packed.dump(2)<<'\n'; output.flush();
        if(!output) throw std::runtime_error("Failed writing thermal record");return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
