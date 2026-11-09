#pragma once
#include "vela/core/Types.h"
#include <cmath>
#include <stdexcept>
#include <string>
namespace vela {
// A terminal functional depends only on its selected equation rows. Do not
// evaluate 0 * nonfinite in unrelated rows of a full transport Jacobian.
inline VectorXd weightedSparseRows(const SparseMatrixd& matrix,const VectorXd& weights,
                                  Eigen::Index* unrelatedNonfinite=nullptr) {
    if(weights.size()!=matrix.rows() || !weights.allFinite())
        throw std::invalid_argument("weightedSparseRows: invalid row weights");
    if(unrelatedNonfinite)*unrelatedNonfinite=0;
    VectorXd result=VectorXd::Zero(matrix.cols());
    for(Eigen::Index col=0;col<matrix.outerSize();++col)
        for(SparseMatrixd::InnerIterator entry(matrix,col);entry;++entry) {
            const Real weight=weights[entry.row()];
            if(weight==0.) {
                if(unrelatedNonfinite && !std::isfinite(entry.value()))++*unrelatedNonfinite;
                continue;
            }
            if(!std::isfinite(entry.value()))
                throw std::runtime_error("terminal current derivative nonfinite at selected row="+
                    std::to_string(entry.row())+" column="+std::to_string(entry.col()));
            result[entry.col()]+=entry.value()*weight;
        }
    return result;
}
}
