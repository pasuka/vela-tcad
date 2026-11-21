#pragma once
// Optional high-precision residual accumulation for sparse iterative refinement.
#include "vela/solver/LinearSolver.h"
#include "vela/core/PerformanceProfiler.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <type_traits>

namespace vela::linear_refinement {
using MP = boost::multiprecision::cpp_dec_float_100;

template<class Scalar>
struct BasicDefect {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> values;
    Scalar backwardMax = 0;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> componentwiseBackwardErrors;
};
using Defect = BasicDefect<double>;

// Opt-in forensic export. Never changes a direction or a convergence decision.
// The caller owns an isolated output directory; each system gets a fresh child.
inline void exportFailedSystem(const vela::SparseMatrixd& matrix,
                               const vela::VectorXd& rhs,
                               const vela::VectorXd& direction,
                               const vela::VectorXd& rowScale,
                               const vela::VectorXd& columnScale,
                               const Defect& defect) {
    const char* root = std::getenv("VELA_BORDERED_FAILURE_DUMP");
    if (!root || !*root) return;
    static std::atomic<unsigned long> sequence{0};
    std::filesystem::path folder;
    do {
        folder = std::filesystem::path(root) /
            ("system_" + std::to_string(sequence.fetch_add(1)));
        std::filesystem::create_directories(folder.parent_path());
    } while (!std::filesystem::create_directory(folder));
    std::ofstream a(folder / "A.mtx"), v(folder / "vectors.csv");
    a.exceptions(std::ios::badbit | std::ios::failbit);
    v.exceptions(std::ios::badbit | std::ios::failbit);
    a << "%%MatrixMarket matrix coordinate real general\n"
      << matrix.rows() << ' ' << matrix.cols() << ' ' << matrix.nonZeros() << '\n'
      << std::setprecision(17);
    for (int col=0; col<matrix.outerSize(); ++col)
        for (vela::SparseMatrixd::InnerIterator e(matrix,col); e; ++e)
            a << e.row()+1 << ' ' << e.col()+1 << ' ' << e.value() << '\n';
    v << "index,rhs,direction,row_scale,column_scale,defect,backward_error\n"
      << std::setprecision(17);
    for (Eigen::Index i=0; i<rhs.size(); ++i)
        v << i << ',' << rhs[i] << ',' << direction[i] << ',' << rowScale[i]
          << ',' << columnScale[i] << ',' << defect.values[i] << ','
          << defect.componentwiseBackwardErrors[i] << '\n';
}

namespace detail {
template<class Scalar, class Coefficient>
BasicDefect<Scalar> accumulateDefect(
    const Eigen::SparseMatrix<Scalar>& J,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& F,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& step,
    Coefficient coefficient) {
    if (J.rows() != F.size() || J.cols() != step.size() || !F.allFinite() || !step.allFinite())
        throw std::invalid_argument("Invalid refinement defect dimensions or input");
    ScopedPerformanceTimer timer("linear.mp_defect.total");
    incrementPerformanceCounter("linear.mp_defect.calls");
    std::vector<MP> r(F.size()), denominator(F.size());
    for (int i=0;i<F.size();++i) { r[i]=MP(F(i)); denominator[i]=abs(r[i]); }
    // Reuse the identical 100-digit conversion of each direction coordinate.
    // Products and row accumulation keep their original precision and order.
    std::vector<MP> preciseStep;
    preciseStep.reserve(static_cast<std::size_t>(step.size()));
    for (int col=0;col<step.size();++col) preciseStep.emplace_back(step(col));
    for (int col=0;col<J.outerSize();++col)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(J,col);it;++it) {
            const MP term=coefficient(it)*preciseStep[it.col()];
            r[it.row()]+=term; denominator[it.row()]+=abs(term);
        }
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    BasicDefect<Scalar> answer{Vector(F.size()), 0, Vector::Zero(F.size())};
    for (int i=0;i<F.size();++i) {
        answer.values(i)=static_cast<Scalar>(r[i]);
        if (denominator[i]!=0) {
            answer.componentwiseBackwardErrors(i)=static_cast<Scalar>(abs(r[i])/denominator[i]);
            answer.backwardMax=std::max(answer.backwardMax,answer.componentwiseBackwardErrors(i));
        }
    }
    if (!answer.values.allFinite() || !std::isfinite(answer.backwardMax))
        throw std::runtime_error("Nonfinite linear refinement defect");
    return answer;
}
} // namespace detail

template<class Scalar>
BasicDefect<Scalar> computeHighPrecisionDefect(
    const Eigen::SparseMatrix<Scalar>& J,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& F,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& step) {
    return detail::accumulateDefect(J, F, step,
        [](const typename Eigen::SparseMatrix<Scalar>::InnerIterator& entry) {
            return MP(entry.value());
        });
}

// The matrix and its sparse traversal order must stay fixed while refining one system.
template<class Scalar>
class PreparedMpDefect {
public:
    using Matrix = Eigen::SparseMatrix<Scalar>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    explicit PreparedMpDefect(const Matrix& matrix) : matrix_(matrix) {
        ScopedPerformanceTimer timer("linear.mp_defect.prepare");
        coefficients_.reserve(static_cast<std::size_t>(matrix.nonZeros()));
        for (Eigen::Index col = 0; col < matrix.outerSize(); ++col)
            for (typename Matrix::InnerIterator entry(matrix, col); entry; ++entry)
                coefficients_.emplace_back(entry.value());
    }
    PreparedMpDefect(Matrix&&) = delete;

    BasicDefect<Scalar> evaluate(const Vector& residual, const Vector& step) const {
        std::size_t index = 0;
        return detail::accumulateDefect(matrix_, residual, step,
            [&](const typename Matrix::InnerIterator&) -> const MP& {
                return coefficients_[index++];
            });
    }

private:
    const Matrix& matrix_;
    std::vector<MP> coefficients_;
};

inline Defect highPrecisionDefect(const SparseMatrixd& J, const VectorXd& F,
                                  const VectorXd& step) {
    return computeHighPrecisionDefect<double>(J, F, step);
}

// BVMethods requests this only for the PARDISO bordered path.
bool useBvMethodsEftDefect();
Defect eftDefect(const SparseMatrixd&, const VectorXd&, const VectorXd&);

template<class Scalar, class Solve, class Observer>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> refineWithSolve(
    const Eigen::SparseMatrix<Scalar>& J,
    const std::type_identity_t<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>& F,
    const std::type_identity_t<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>& rowWeights,
    std::type_identity_t<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> step,
    Solve solve, int corrections, Observer observe) {
    if (corrections<0 || J.rows()!=F.size() || J.cols()!=step.size() || rowWeights.size()!=F.size())
        throw std::invalid_argument("Invalid linear refinement dimensions/count");
    if (!rowWeights.allFinite() || (rowWeights.array() <= 0).any())
        throw std::invalid_argument("Refinement row weights must be finite and positive");
    Eigen::SparseMatrix<Scalar> scaled=J;
    for (int col=0;col<scaled.outerSize();++col)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(scaled,col);it;++it)
            it.valueRef()*=rowWeights(it.row());
    for (int k=0;k<=corrections;++k) {
        const auto defect=computeHighPrecisionDefect(J,F,step);
        observe(k,step,defect);
        if (k<corrections) {
            step+=solve(scaled,(-defect.values.cwiseProduct(rowWeights)).eval());
            if (!step.allFinite()) throw std::runtime_error("Nonfinite linear refinement update");
        }
    }
    return step;
}

template<class Observer>
vela::VectorXd refine(const vela::SparseMatrixd& J, const vela::VectorXd& F,
                      const vela::VectorXd& rowWeights, vela::VectorXd step,
                      vela::LinearSolver& solver, int corrections, Observer observe) {
    return refineWithSolve(J,F,rowWeights,std::move(step),
        [&solver](const vela::SparseMatrixd& matrix, const vela::VectorXd& rhs) {
            return solver.solve(matrix,rhs);
        },corrections,observe);
}
} // namespace vela::linear_refinement
