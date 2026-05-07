#include "vela/simulation/DCSweep.h"
#include "vela/simulation/PoissonSimulation.h"

#include <exception>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

namespace {

void usage(const char* argv0)
{
    std::cerr << "Usage: " << argv0 << " --config <simulation.json>\n";
}

} // namespace

int main(int argc, char** argv)
{
    std::string configFile;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--config" && i + 1 < argc) {
            configFile = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return 0;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (configFile.empty()) {
        usage(argv[0]);
        return 2;
    }

    try {
        std::ifstream ifs(configFile);
        if (!ifs.is_open()) {
            std::cerr << "Cannot open config file: " << configFile << '\n';
            return 1;
        }

        nlohmann::json cfg;
        ifs >> cfg;
        const std::string type = cfg.value("simulation_type", cfg.contains("sweep") ? "dc_sweep" : "poisson");

        nlohmann::json status;
        status["config"] = configFile;
        status["simulation_type"] = type;
        status["converged"] = true;

        if (type == "dc_sweep") {
            vela::DCSweep sweep;
            const auto points = sweep.run(configFile);
            bool allConverged = !points.empty();
            for (const auto& point : points)
                allConverged = allConverged && point.converged;
            status["converged"] = allConverged;
            status["points"] = points.size();
        } else if (type == "poisson") {
            vela::PoissonSimulation sim;
            const auto psi = sim.run(configFile);
            status["nodes"] = psi.size();
        } else {
            std::cerr << "Unknown simulation_type: " << type << '\n';
            return 2;
        }

        std::cout << status.dump() << '\n';
        return status.value("converged", false) ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cerr << "vela_example_runner failed: " << ex.what() << '\n';
        return 1;
    }
}
