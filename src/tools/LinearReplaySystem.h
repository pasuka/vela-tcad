#pragma once
// Replay frozen VELALU01 Newton systems through the actual shared solver.
#include "vela/solver/LinearSolver.h"
#include "vela/core/PerformanceProfiler.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using namespace vela;
static_assert(sizeof(double) == 8 && sizeof(SparseMatrixd::StorageIndex) == 4);
static_assert(std::endian::native == std::endian::little && !SparseMatrixd::IsRowMajor);
using Json = nlohmann::json;
using Clock = std::chrono::steady_clock;
static double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

struct Reader {
    std::ifstream input;
    std::uint64_t remaining;
    explicit Reader(const std::filesystem::path& path)
        : input(path, std::ios::binary), remaining(std::filesystem::file_size(path)) {
        require(input.good() && remaining <= 128ULL * 1024 * 1024, "Invalid capture file size");
    }
    void bytes(void* destination, std::uint64_t length) {
        require(length <= remaining, "Truncated capture");
        input.read(static_cast<char*>(destination), static_cast<std::streamsize>(length));
        require(input.good(), "Capture read failed");
        remaining -= length;
    }
    std::uint64_t size() {
        std::uint64_t value; bytes(&value, 8); return value;
    }
    SparseMatrixd matrix() {
        const auto rows = size(), cols = size(), count = size();
        require(rows > 0 && rows == cols && rows <= 100000 && count <= 10000000,
                "Invalid matrix dimensions");
        require(4 * (cols + 1) + 12 * count <= remaining, "Truncated matrix arrays");
        std::vector<int> outer(cols + 1), inner(count);
        std::vector<double> values(count);
        bytes(outer.data(), 4 * outer.size());
        bytes(inner.data(), 4 * count); bytes(values.data(), 8 * count);
        require(outer[0] == 0 && outer[cols] == static_cast<int>(count), "Invalid CSC boundaries");
        for (std::size_t col = 0; col < cols; ++col) {
            require(outer[col] >= 0 && outer[col] <= outer[col + 1] &&
                    outer[col + 1] <= static_cast<int>(count), "Invalid CSC offsets");
            int previous = -1;
            for (int j = outer[col]; j < outer[col + 1]; ++j) {
                require(inner[j] > previous && inner[j] < static_cast<int>(rows) &&
                        std::isfinite(values[j]), "Invalid CSC entries");
                previous = inner[j];
            }
        }
        return Eigen::Map<const SparseMatrixd>(rows, cols, count,
            outer.data(), inner.data(), values.data());
    }
    VectorXd vector(Eigen::Index expected) {
        const auto count = size();
        require(count == static_cast<std::uint64_t>(expected) && 8 * count <= remaining,
                "Invalid vector length");
        VectorXd result(expected); bytes(result.data(), 8 * count);
        require(result.allFinite(), "Nonfinite capture vector"); return result;
    }
};

struct System {
    bool rawAvailable=true;
    std::uint64_t sequence, iteration;
    SparseMatrixd matrix, rawMatrix;
    VectorXd rhs, rawRhs, columns, rows, l2Rows, reference;
    explicit System(const std::filesystem::path& path) {
        Reader reader(path); char magic[8]; reader.bytes(magic, 8);
        rawAvailable=std::memcmp(magic,"VELALU01",8)==0;
        require(rawAvailable || std::memcmp(magic,"VELALU02",8)==0,"Invalid capture magic/version");
        sequence = reader.size(); iteration = reader.size();
        matrix = reader.matrix(); rhs = reader.vector(matrix.rows());
        if(!rawAvailable) {
            require(matrix.rows()%3==0,"Expected three-block solver input");
            reference=reader.vector(matrix.rows());rawMatrix=matrix;rawRhs=rhs;
            columns=rows=l2Rows=VectorXd::Ones(matrix.rows());
            require(reader.remaining==0,"Trailing capture data");return;
        }
        rawMatrix = reader.matrix();
        require(rawMatrix.rows() == matrix.rows() && matrix.rows() % 3 == 0,
                "Capture is not a primary three-block Newton system");
        rawRhs = reader.vector(matrix.rows()); columns = reader.vector(matrix.rows());
        rows = reader.vector(matrix.rows()); l2Rows = reader.vector(matrix.rows());
        reference = reader.vector(matrix.rows());
        require(reader.remaining == 0, "Trailing capture data");
        require((columns.array() > 0).all() && (rows.array() > 0).all() &&
                (l2Rows.array() > 0).all(), "Invalid scaling vectors");
        require(matrix.nonZeros() == rawMatrix.nonZeros(), "Different raw/scaled pattern size");
        require(std::memcmp(matrix.outerIndexPtr(), rawMatrix.outerIndexPtr(),
                            4 * (matrix.cols() + 1)) == 0 &&
                std::memcmp(matrix.innerIndexPtr(), rawMatrix.innerIndexPtr(),
                            4 * matrix.nonZeros()) == 0, "Different raw/scaled pattern");
        SparseMatrixd reconstructed = rawMatrix;
        for (int col = 0; col < reconstructed.outerSize(); ++col) {
            for (SparseMatrixd::InnerIterator entry(reconstructed, col); entry; ++entry) {
                entry.valueRef() *= rows(entry.row());
                entry.valueRef() *= l2Rows(entry.row());
                entry.valueRef() *= columns(entry.col());
            }
        }
        const VectorXd reconstructedRhs = rawRhs.cwiseProduct(rows).cwiseProduct(l2Rows);
        require(std::memcmp(reconstructed.valuePtr(), matrix.valuePtr(), 8 * matrix.nonZeros()) == 0 &&
                std::memcmp(reconstructedRhs.data(), rhs.data(), 8 * rhs.size()) == 0,
                "Scaling reconstruction is not bit-exact");
    }
};

static Json errors(const SparseMatrixd& matrix, const VectorXd& rhs, const VectorXd& x) {
    require(x.allFinite(), "Nonfinite solve result");
    const VectorXd residual = matrix * x - rhs;
    VectorXd denominator = rhs.cwiseAbs();
    for (int col = 0; col < matrix.outerSize(); ++col)
        for (SparseMatrixd::InnerIterator entry(matrix, col); entry; ++entry)
            denominator(entry.row()) += std::abs(entry.value()) * std::abs(x(entry.col()));
    VectorXd ratios(residual.size());
    for (int i = 0; i < residual.size(); ++i) {
        ratios(i) = denominator(i) > 0 ? std::abs(residual(i)) / denominator(i) : 0;
        require(denominator(i) > 0 || residual(i) == 0, "Nonzero residual on empty row");
    }
    require(ratios.allFinite(), "Nonfinite backward error");
    Json blocks = Json::array(); const int n = static_cast<int>(rhs.size()) / 3;
    for (int block = 0; block < 3; ++block)
        blocks.push_back(ratios.segment(block * n, n).maxCoeff());
    return Json{{"componentwise_backward_error", ratios.maxCoeff()},
                {"block_componentwise_backward_error", blocks},
                {"residual_norm", residual.norm()},
                {"normwise_backward_error", residual.norm() /
                    std::max(matrix.norm() * x.norm() + rhs.norm(), std::numeric_limits<double>::min())}};
}
static Json quality(const System& system, const VectorXd& coordinates) {
    const VectorXd step = coordinates.cwiseProduct(system.columns);
    const VectorXd referenceStep = system.reference.cwiseProduct(system.columns);
    Json difference = Json::array(); const int n = static_cast<int>(step.size()) / 3;
    for (int block = 0; block < 3; ++block) {
        const auto expected = referenceStep.segment(block * n, n);
        difference.push_back((step.segment(block * n, n) - expected).norm() /
            std::max(expected.norm(), std::numeric_limits<double>::min()));
    }
    return Json{{"scaled", errors(system.matrix, system.rhs, coordinates)},
                {"raw", errors(system.rawMatrix, system.rawRhs, step)},
                {"relative_raw_step_difference_by_block", difference}};
}
