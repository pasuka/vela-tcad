// Diagnostic-only projection/repacking using the production Wide type.
#include "vela/numerics/SplitDDState.h"
#include <iostream>
#include <stdexcept>

int main() {
    using vela::split_dd::Wide;
    nlohmann::json input;
    std::cin >> input;
    const double scale = input.at("scale");
    if (!std::isfinite(scale) || scale <= 0) throw std::runtime_error("Invalid scale");
    nlohmann::json out = nlohmann::json::array();
    for (const auto& row : input.at("rows")) {
        nlohmann::json result;
        for (int block=0; block<3; ++block) {
            const double value=row.at(block), ref=block ? row.at(2+block).get<double>() : 0.;
            if(!std::isfinite(value)||!std::isfinite(ref)) throw std::runtime_error("Nonfinite input");
            const Wide coord=(Wide(value)-Wide(ref))/Wide(scale);
            const double hi=static_cast<double>(coord);
            const double lo=static_cast<double>(coord-Wide(hi));
            const Wide reconstructed=Wide(hi)+Wide(lo);
            result["hi"].push_back(hi);
            result["lo"].push_back(lo);
            result["physical"].push_back(static_cast<double>(reconstructed*Wide(scale)+Wide(ref)));
            result["increment"].push_back(static_cast<double>(reconstructed*Wide(scale)));
        }
        out.push_back(result);
    }
    std::cout << out.dump() << '\n';
}
