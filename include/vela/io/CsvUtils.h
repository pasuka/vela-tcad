#pragma once

#include <string>
#include <vector>

namespace vela {

std::string trimCsvToken(std::string value);
std::vector<std::string> splitCsvLine(const std::string& line);
std::vector<std::string> splitCsvLine(const std::string& line,
                                      const std::string& quotedFieldsError);

} // namespace vela
