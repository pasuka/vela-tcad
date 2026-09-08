#pragma once
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace simplemos_native_geometry {
// Validation-only coefficient interpolation. No change to state or solver gates.
inline void apply(std::vector<double>& values, const char* pathName,
                  const char* alphaName) {
    const char* alphaText=std::getenv(alphaName);
    if (!alphaText) return;
    const double alpha=std::stod(alphaText);
    if (!std::isfinite(alpha)) throw std::runtime_error("Nonfinite geometry alpha");
    if (alpha==0.) return; // The disabled control must remain bit identical.
    const char* path=std::getenv(pathName);
    if (!path) throw std::runtime_error("Missing native geometry ratios");
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot read native geometry ratios");
    std::vector<double> ratios(values.size());
    for (std::size_t expected=0;expected<values.size();++expected) {
        std::size_t id;double ratio;
        if (!(file>>id>>ratio) || id!=expected || !std::isfinite(ratio) || ratio<0.)
            throw std::runtime_error("Invalid native geometry row");
        ratios[id]=ratio;
    }
    std::string extra;
    if (file>>extra) throw std::runtime_error("Unexpected native geometry rows");
    for (std::size_t i=0;i<values.size();++i) {
        const double multiplier=1.+alpha*(ratios[i]-1.);
        if (multiplier<0. || !std::isfinite(multiplier))
            throw std::runtime_error("Negative interpolated geometry");
        values[i]*=multiplier;
    }
}
inline void transport(std::vector<double>& values) {
    apply(values,"VELA_NATIVE_GEOMETRY_EDGE_RATIOS","VELA_NATIVE_GEOMETRY_T");
}
inline void poissonVolume(std::vector<double>& values) {
    apply(values,"VELA_NATIVE_GEOMETRY_NODE_RATIOS","VELA_NATIVE_GEOMETRY_P");
}
}
