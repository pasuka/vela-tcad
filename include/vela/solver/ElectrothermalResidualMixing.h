#pragma once
#include "vela/core/Types.h"
#include <Eigen/QR>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace vela::experimental {
struct ElectrothermalResidualMix {
    VectorXd direction;
    Index rank=0;
    Real coefficientL1=0.,predictedNorm=0.;
};
// NGMRES-style residual-space correction after a nonlinear Newton map.
// All history must belong to one bias, model and QF reference system.
inline ElectrothermalResidualMix electrothermalResidualMix(const VectorXd& x,const VectorXd& residual,
        const std::vector<VectorXd>& historyX,const std::vector<VectorXd>& historyR,const VectorXd& rows) {
    if(historyX.size()!=historyR.size() || historyX.empty() || historyX.size()>4 || x.size()!=residual.size() ||
       rows.size()!=x.size() || !x.allFinite() || !residual.allFinite() || !rows.allFinite() || (rows.array()<=0.).any())
        throw std::invalid_argument("Invalid residual mixing history");
    const VectorXd rhs=-(residual.array()*rows.array()).matrix();
    if(!rhs.allFinite() || !std::isfinite(rhs.norm()))throw std::invalid_argument("Invalid scaled mixing residual");
    std::vector<VectorXd> differences,states;std::vector<Real> norms;
    for(std::size_t j=0;j<historyX.size();++j) {
        if(historyX[j].size()!=x.size() || historyR[j].size()!=x.size() || !historyX[j].allFinite() || !historyR[j].allFinite())
            throw std::invalid_argument("Invalid residual mixing sample");
        VectorXd difference=((historyR[j]-residual).array()*rows.array()).matrix();
        const Real norm=difference.norm();
        if(!difference.allFinite() || !std::isfinite(norm))throw std::invalid_argument("Invalid scaled history difference");
        if(norm<=1e-14*std::max(rhs.norm(),1e-100))continue;
        differences.push_back(difference/norm);states.push_back(historyX[j]-x);norms.push_back(norm);
    }
    if(differences.empty())throw std::invalid_argument("Residual mixing has no independent change");
    Eigen::MatrixXd matrix(x.size(),differences.size());
    for(std::size_t j=0;j<differences.size();++j)matrix.col(j)=differences[j];
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(matrix);qr.setThreshold(1e-8);
    const VectorXd normalized=qr.solve(rhs);
    ElectrothermalResidualMix result;result.rank=static_cast<Index>(qr.rank());result.direction=VectorXd::Zero(x.size());
    for(std::size_t j=0;j<states.size();++j) {
        const Real coefficient=normalized[j]/norms[j];result.coefficientL1+=std::abs(coefficient);
        result.direction+=coefficient*states[j];
    }
    result.predictedNorm=(matrix*normalized-rhs).norm();
    // A stagnating map has very small history changes, so useful coefficients
    // can be large. Physical step caps and a fresh residual decide acceptance.
    if(!result.rank || !result.direction.allFinite() || !std::isfinite(result.predictedNorm) || !std::isfinite(result.coefficientL1) || result.coefficientL1>1e8)
        throw std::invalid_argument("Residual mixing is ill-conditioned");
    return result;
}
}
