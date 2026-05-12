#pragma once

#include <cstddef>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace vela {

// Scalar and index types
using Real  = double;
using Index = std::size_t;

// Commonly used Eigen types
using VectorXd = Eigen::VectorXd;
using SparseMatrixd = Eigen::SparseMatrix<double>;

// 2-D point (x, y)
using Point2 = Eigen::Vector2d;

// 3-D point (x, y, z) - reserved for future use
using Point3 = Eigen::Vector3d;

/// Cell element types supported by Vela.
enum class CellType : int {
    Tri3 = 3,   ///< Linear triangle (3 nodes)
};

} // namespace vela
