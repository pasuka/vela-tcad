#pragma once
// Optional high-precision residual accumulation for sparse iterative refinement.
#include "vela/solver/LinearSolver.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace vela::linear_refinement {
using MP = boost::multiprecision::cpp_dec_float_100;

struct Defect {
    vela::VectorXd values;
    double backwardMax = 0;
};

inline Defect highPrecisionDefect(const vela::SparseMatrixd& J,
                                 const vela::VectorXd& F,
                                 const vela::VectorXd& step) {
    if (J.rows() != F.size() || J.cols() != step.size() || !F.allFinite() || !step.allFinite())
        throw std::invalid_argument("Invalid refinement defect dimensions or input");
    std::vector<MP> r(F.size()), denominator(F.size());
    for (int i=0;i<F.size();++i) { r[i]=MP(F(i)); denominator[i]=abs(r[i]); }
    for (int col=0;col<J.outerSize();++col)
        for (vela::SparseMatrixd::InnerIterator it(J,col);it;++it) {
            const MP term=MP(it.value())*MP(step(it.col()));
            r[it.row()]+=term; denominator[it.row()]+=abs(term);
        }
    Defect answer{vela::VectorXd(F.size()),0};
    for (int i=0;i<F.size();++i) {
        answer.values(i)=static_cast<double>(r[i]);
        if (denominator[i]!=0)
            answer.backwardMax=std::max(answer.backwardMax,static_cast<double>(abs(r[i])/denominator[i]));
    }
    if (!answer.values.allFinite() || !std::isfinite(answer.backwardMax))
        throw std::runtime_error("Nonfinite linear refinement defect");
    return answer;
}

template<class Observer>
vela::VectorXd refine(const vela::SparseMatrixd& J, const vela::VectorXd& F,
                      const vela::VectorXd& rowWeights, vela::VectorXd step,
                      vela::LinearSolver& solver, int corrections, Observer observe) {
    if (corrections<0 || J.rows()!=F.size() || J.cols()!=step.size() || rowWeights.size()!=F.size())
        throw std::invalid_argument("Invalid linear refinement dimensions/count");
    if (!rowWeights.allFinite() || (rowWeights.array() <= 0).any())
        throw std::invalid_argument("Refinement row weights must be finite and positive");
    vela::SparseMatrixd scaled=J;
    for (int col=0;col<scaled.outerSize();++col)
        for (vela::SparseMatrixd::InnerIterator it(scaled,col);it;++it)
            it.valueRef()*=rowWeights(it.row());
    for (int k=0;k<=corrections;++k) {
        const Defect defect=highPrecisionDefect(J,F,step);
        observe(k,step,defect);
        if (k<corrections) {
            step+=solver.solve(scaled,-defect.values.cwiseProduct(rowWeights));
            if (!step.allFinite()) throw std::runtime_error("Nonfinite linear refinement update");
        }
    }
    return step;
}
} // namespace vela::linear_refinement
