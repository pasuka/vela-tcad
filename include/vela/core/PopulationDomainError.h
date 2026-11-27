#pragma once

#include <stdexcept>

namespace vela {

// A trial state has left the finite, strictly positive population domain.
// Only trial evaluation may recover from this signal; invalid parents and
// unrelated assembly/configuration errors must still propagate.
class PopulationDomainError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

} // namespace vela
