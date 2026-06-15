#include "vela/io/CsvUtils.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace vela {

std::string trimCsvToken(std::string value)
{
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char ch) {
        return std::isspace(ch);
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char ch) {
        return std::isspace(ch);
    }).base();
    if (first >= last)
        return {};
    return std::string(first, last);
}

std::vector<std::string> splitCsvLine(const std::string& line)
{
    return splitCsvLine(line, "CSV parser does not support quoted fields.");
}

std::vector<std::string> splitCsvLine(const std::string& line,
                                      const std::string& quotedFieldsError)
{
    if (line.find('"') != std::string::npos) {
        throw std::runtime_error(quotedFieldsError);
    }

    std::vector<std::string> columns;
    std::stringstream ss(line);
    std::string column;
    while (std::getline(ss, column, ','))
        columns.push_back(trimCsvToken(column));
    return columns;
}

} // namespace vela
