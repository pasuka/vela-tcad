#pragma once
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace vela::diagnostics {
inline void writeExperimentJson(const std::filesystem::path& path,const nlohmann::json& value)
{
    auto temporary=path;temporary+=".tmp";
    {
        std::ofstream stream(temporary);
        if(!stream)throw std::runtime_error("cannot open experiment journal temporary file");
        stream<<std::setw(2)<<value<<'\n';stream.flush();
        if(!stream)throw std::runtime_error("cannot write experiment journal");
    }
#ifdef _WIN32
    if(!MoveFileExW(temporary.c_str(),path.c_str(),MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("cannot replace experiment journal");
#else
    std::filesystem::rename(temporary,path);
#endif
}
} // namespace vela::diagnostics
