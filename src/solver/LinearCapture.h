#pragma once
#include "vela/core/Types.h"
#include <string>
namespace vela::detail {
// Opt-in diagnostic only: captures solver-input coordinates, not unscaled physics.
void captureLinearInput(const std::string& directory,const SparseMatrixd&,const VectorXd&,const VectorXd&);
void captureNewtonMatrix(const std::string& path, const SparseMatrixd&);
void captureNewtonVectors(const std::string& path, const VectorXd& state,
    const VectorXd& residual, const VectorXd& direction, const VectorXd& referenceVolts);
}
