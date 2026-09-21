#pragma once
#include "vela/core/Types.h"
#include <string>
namespace vela::detail {
// Opt-in diagnostic only: captures solver-input coordinates, not unscaled physics.
void captureLinearInput(const std::string& directory,const SparseMatrixd&,const VectorXd&,const VectorXd&);
}
