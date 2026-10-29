#pragma once
#include "vela/solver/LinearSolver.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
namespace vela {
struct CommonRelativeLinearResult {
    SparseMatrixd matrix,backendMatrix;
    VectorXd direction,shared,relative,rowScale,columnScale,scaledRhs;
};
// Right-coordinate change: x_i = common + relative_i, with the anchor's
// relative coordinate zero. Caller supplies analytically evaluated common
// columns; these need not equal a cancellation-prone sum of stored doubles.
inline CommonRelativeLinearResult solveCommonRelative(
    const SparseMatrixd& a,const VectorXd& b,const std::vector<std::vector<int>>& groups,
    const std::vector<VectorXd>& commonColumns,const std::string& backend) {
    const int n=a.rows();
    if(a.cols()!=n || b.size()!=n || groups.size()!=commonColumns.size() || !b.allFinite())
        throw std::invalid_argument("Invalid common-relative system");
    std::set<int> used,anchors;
    for(std::size_t k=0;k<groups.size();++k) {
        if(groups[k].empty() || commonColumns[k].size()!=n || !commonColumns[k].allFinite())
            throw std::invalid_argument("Invalid common-relative column");
        for(int i:groups[k])if(i<0 || i>=n || !used.insert(i).second)
            throw std::invalid_argument("Overlapping common-relative coordinates");
        anchors.insert(groups[k].front());
    }
    std::vector<Eigen::Triplet<Real>> entries;
    for(int c=0;c<a.outerSize();++c)if(!anchors.contains(c))
        for(SparseMatrixd::InnerIterator e(a,c);e;++e){if(!std::isfinite(e.value()))throw std::invalid_argument("Nonfinite common-relative matrix");entries.emplace_back(e.row(),c,e.value());}
    for(std::size_t k=0;k<groups.size();++k)for(int r=0;r<n;++r)
        if(commonColumns[k][r]!=0)entries.emplace_back(r,groups[k].front(),commonColumns[k][r]);
    CommonRelativeLinearResult out;out.matrix.resize(n,n);out.matrix.setFromTriplets(entries.begin(),entries.end());
    out.rowScale=VectorXd::Zero(n);out.columnScale=VectorXd::Zero(n);
    for(int c=0;c<n;++c)for(SparseMatrixd::InnerIterator e(out.matrix,c);e;++e)
        out.rowScale[e.row()]=std::max(out.rowScale[e.row()],std::abs(e.value()));
    for(int r=0;r<n;++r){if(!(out.rowScale[r]>0))throw std::runtime_error("Zero common-relative row");out.rowScale[r]=1/std::max(out.rowScale[r],std::abs(b[r]));}
    for(int c=0;c<n;++c)for(SparseMatrixd::InnerIterator e(out.matrix,c);e;++e)
        out.columnScale[c]=std::max(out.columnScale[c],std::abs(e.value()*out.rowScale[e.row()]));
    for(int c=0;c<n;++c){if(!(out.columnScale[c]>0))throw std::runtime_error("Zero common-relative column");out.columnScale[c]=1/out.columnScale[c];}
    out.backendMatrix=out.matrix;
    for(int c=0;c<n;++c)for(SparseMatrixd::InnerIterator e(out.backendMatrix,c);e;++e)
        e.valueRef()*=out.rowScale[e.row()]*out.columnScale[c];
    out.scaledRhs=out.rowScale.array()*b.array();LinearSolver solver(backend);
    out.direction=out.columnScale.array()*solver.solve(out.backendMatrix,out.scaledRhs).array();
    if(!out.direction.allFinite())throw std::runtime_error("Nonfinite common-relative direction");
    out.relative=out.direction;out.shared=VectorXd::Zero(n);
    for(const auto& group:groups){out.relative[group.front()]=0;for(int i:group)out.shared[i]=out.direction[group.front()];}
    return out;
}
}
