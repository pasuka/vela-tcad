#pragma once

#include "vela/solver/LinearSolver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <optional>

namespace vela {

/// Explicit [psi, phin, phip, contact_voltage] layout. Retain every psi row
/// and the circuit row; eliminate exact carrier constraints whose only other
/// coefficient is the shared voltage. Substitute dx_i=c_i+s_i*dV into every
/// retained equation, including the current row.
/// Unlike the 3N keep_psi backend this accepts an actual, non-dummy border.
struct CarrierReducedLoadLine {
    SparseMatrixd full, reduced, equilibrated;
    VectorXd rhs, reducedRhs, scaledRhs, rowScale, columnScale, fixed;
    VectorXd voltageSlope;
    std::vector<int> map;

    VectorXd lift(const VectorXd& freeDirection) const {
        if (freeDirection.size() != reduced.rows())
            throw std::invalid_argument("load-line lift dimension mismatch");
        VectorXd result = fixed;
        for (Eigen::Index i = 0; i < result.size(); ++i)
            if (map[i] >= 0) result[i] = freeDirection[map[i]];
            else result[i] += voltageSlope[i]*freeDirection[map.back()];
        return result;
    }

    VectorXd solve(const std::string& backend) const {
        // Explicit selection: no environment-dependent fallback or second
        // elimination of the already reduced, non-3N system.
        if (backend != "pardiso_mkl" && backend != "sparselu" && backend != "umfpack")
            throw std::invalid_argument("load-line probe needs an unreduced direct backend");
        LinearSolver solver(backend);
        const VectorXd y = solver.solve(equilibrated, scaledRhs);
        return lift(VectorXd(columnScale.array() * y.array()));
    }
};

inline CarrierReducedLoadLine assembleCarrierReducedLoadLine(
    const SparseMatrixd& jacobian, const VectorXd& residual,
    const VectorXd& voltageColumn, const VectorXd& currentGradient,
    Real currentVoltageDerivative, Real current, Real voltage, Real supply,
    Real resistance, int nodeCount,std::optional<Real> circuitResidual = {})
{
    const int n = 3 * nodeCount;
    if (nodeCount <= 0 || jacobian.rows() != n || jacobian.cols() != n ||
        residual.size() != n || voltageColumn.size() != n ||
        currentGradient.size() != n || !(resistance > 0.0) ||
        !residual.allFinite() || !voltageColumn.allFinite() || !currentGradient.allFinite() ||
        !std::isfinite(currentVoltageDerivative) || !std::isfinite(current) ||
        !std::isfinite(voltage) || !std::isfinite(supply) || !std::isfinite(resistance))
        throw std::invalid_argument("invalid carrier-reduced load-line inputs");
    CarrierReducedLoadLine out;
    std::vector<Eigen::Triplet<Real>> entries;
    entries.reserve(jacobian.nonZeros() + 2*n + 1);
    for (int col = 0; col < jacobian.outerSize(); ++col)
        for (SparseMatrixd::InnerIterator e(jacobian, col); e; ++e)
            entries.emplace_back(e.row(), e.col(), e.value());
    for (int i = 0; i < n; ++i) {
        if (voltageColumn[i] != 0) entries.emplace_back(i, n, voltageColumn[i]);
        if (currentGradient[i] != 0) entries.emplace_back(n, i, resistance*currentGradient[i]);
    }
    entries.emplace_back(n, n, 1 + resistance*currentVoltageDerivative);
    out.full.resize(n+1,n+1);
    out.full.setFromTriplets(entries.begin(),entries.end());
    out.full.makeCompressed();
    out.rhs.resize(n+1);
    out.rhs.head(n) = -residual;
    if(circuitResidual && !std::isfinite(*circuitResidual))throw std::invalid_argument("Nonfinite circuit residual");
    out.rhs[n] = -circuitResidual.value_or(voltage + resistance*current - supply);
    std::vector<int> count(n+1,0);
    std::vector<bool> diagonalOrVoltage(n+1,true);
    for (int col=0; col<=n; ++col)
        for (SparseMatrixd::InnerIterator e(out.full,col);e;++e) {
            if (!std::isfinite(e.value())) throw std::runtime_error("nonfinite circuit matrix");
            if (e.value()!=0) {
                ++count[e.row()];
                if (e.col()!=e.row() && e.col()!=n) diagonalOrVoltage[e.row()]=false;
            }
        }
    out.map.assign(n+1,-1);
    out.fixed = VectorXd::Zero(n+1);
    out.voltageSlope=VectorXd::Zero(n+1);
    std::vector<long double> fixed(n+1,0), corrected(n+1);
    int nf=0;
    for (int i=0;i<=n;++i) {
        const Real diagonal=out.full.coeff(i,i);
        if (i>=nodeCount && i<n && count[i]>0 && diagonalOrVoltage[i] && diagonal!=0) {
            out.fixed[i] = static_cast<Real>(fixed[i]=static_cast<long double>(out.rhs[i])/diagonal);
            out.voltageSlope[i]=-out.full.coeff(i,n)/diagonal;
        }
        else out.map[i]=nf++;
        corrected[i]=out.rhs[i];
    }
    entries.clear();
    for (int col=0;col<=n;++col)
        for (SparseMatrixd::InnerIterator e(out.full,col);e;++e) {
            if (out.map[e.row()]<0 || e.value()==0) continue;
            if (out.map[col]<0) {
                corrected[e.row()] -= static_cast<long double>(e.value())*fixed[col];
                if(out.voltageSlope[col]!=0)
                    entries.emplace_back(out.map[e.row()],out.map[n],e.value()*out.voltageSlope[col]);
            }
            else entries.emplace_back(out.map[e.row()],out.map[col],e.value());
        }
    out.reduced.resize(nf,nf);
    out.reduced.setFromTriplets(entries.begin(),entries.end());
    out.reduced.makeCompressed();
    out.reducedRhs.resize(nf);
    for (int i=0;i<=n;++i)
        if (out.map[i]>=0) out.reducedRhs[out.map[i]]=static_cast<Real>(corrected[i]);
    out.rowScale=VectorXd::Zero(nf);
    for (int col=0;col<nf;++col)
        for (SparseMatrixd::InnerIterator e(out.reduced,col);e;++e)
            out.rowScale[e.row()]=std::max(out.rowScale[e.row()],std::abs(e.value()));
    for (int i=0;i<nf;++i) {
        if (!(out.rowScale[i]>0)) throw std::runtime_error("zero reduced circuit row");
        out.rowScale[i]=1/std::max(out.rowScale[i],std::abs(out.reducedRhs[i]));
    }
    out.columnScale=VectorXd::Zero(nf);
    for (int col=0;col<nf;++col)
        for (SparseMatrixd::InnerIterator e(out.reduced,col);e;++e)
            out.columnScale[col]=std::max(out.columnScale[col],std::abs(e.value()*out.rowScale[e.row()]));
    for (int i=0;i<nf;++i) {
        if (!(out.columnScale[i]>0)) throw std::runtime_error("zero reduced circuit column");
        out.columnScale[i]=1/out.columnScale[i];
    }
    out.equilibrated=out.reduced;
    for (int col=0;col<nf;++col)
        for (SparseMatrixd::InnerIterator e(out.equilibrated,col);e;++e)
            e.valueRef()*=out.rowScale[e.row()]*out.columnScale[col];
    out.scaledRhs=out.rowScale.array()*out.reducedRhs.array();
    if (!out.scaledRhs.allFinite()) throw std::runtime_error("nonfinite scaled circuit RHS");
    return out;
}
} // namespace vela
