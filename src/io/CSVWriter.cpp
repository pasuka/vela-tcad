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
        const auto& value = values[i];
        if (value.find_first_of(",\"\r\n") == std::string::npos) {
            ofs_ << value;
        } else {
            ofs_ << '"';
            for (char ch : value) {
                if (ch == '"') ofs_ << '"';
                ofs_ << ch;
            }
            ofs_ << '"';
        }
    }
    ofs_ << '\n';
}

} // namespace vela
