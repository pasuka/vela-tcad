#include "vela/io/CSVWriter.h"
#include <stdexcept>

namespace vela {

CSVWriter::CSVWriter(const std::string& filename)
    : ofs_(filename)
{
    if (!ofs_.is_open())
        throw std::runtime_error("CSVWriter: cannot open file: " + filename);
}

void CSVWriter::writeHeader(const std::vector<std::string>& columns)
{
    writeRow(columns);
}

void CSVWriter::writeRow(const std::vector<std::string>& values)
{
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) ofs_ << ',';
        ofs_ << values[i];
    }
    ofs_ << '\n';
}

} // namespace vela
