#pragma once
#include "vela/core/Types.h"
#include <memory>
#include <string>
#include <vector>

namespace vela::detail {
// Lifecycle is controlled by LinearSolver's exact structure/value cache.
class DirectBackend {
public:
    virtual ~DirectBackend() = default;
    virtual void analyze(const SparseMatrixd&) = 0;
    virtual void factor(const SparseMatrixd&) = 0;
    virtual VectorXd solve(const VectorXd&) = 0;
    virtual void statistics() const {}
};
// Output is new column -> original column. Uses the undirected structural union,
// including explicitly stored zeros; never changes the numerical equations.
std::vector<int> metisPermutation(const SparseMatrixd&);
std::unique_ptr<DirectBackend> makeMetisSparseLU();
std::unique_ptr<DirectBackend> makeMumps(bool metis);
std::unique_ptr<DirectBackend> makeSuperLuMt(bool metis);
std::unique_ptr<DirectBackend> makeStrumpack();
int directBackendThreads();
} // namespace vela::detail
