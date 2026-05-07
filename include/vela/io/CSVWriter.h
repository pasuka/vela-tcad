#pragma once

#include <fstream>
#include <string>
#include <vector>

namespace vela {

class CSVWriter {
public:
    explicit CSVWriter(const std::string& filename);

    void writeHeader(const std::vector<std::string>& columns);
    void writeRow(const std::vector<std::string>& values);

private:
    std::ofstream ofs_;
};

} // namespace vela
