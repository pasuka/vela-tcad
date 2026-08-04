#pragma once

#include "vela/preprocess/DeviceIr2D.h"
#include <string>

namespace vela {

class SdeScriptReader {
public:
    DeviceIr2D parseFile(const std::string& filename) const;
    DeviceIr2D parseText(const std::string& text, const std::string& sourceName = "<memory>") const;
};

} // namespace vela

