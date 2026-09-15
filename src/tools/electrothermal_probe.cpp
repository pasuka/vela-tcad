// Compatibility CLI for the shared four-equation simulation service.
#include "vela/simulation/ElectrothermalSimulation.h"
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
        const auto result=vela::solveElectrothermalPoint(cfg, std::cout);
        std::ofstream output(argv[2]); if(!output) throw std::runtime_error("Cannot open output");
        output<<result.dump(2)<<'\n'; return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
