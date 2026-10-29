#pragma once
#include "vela/equation/experimental/ResidualAudit.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
// This reconstructs the sum of captured binary64 operands. It does not recover
// precision already lost inside a physics kernel or promote a failed Jv gate.
template<class Scalar>
std::vector<Scalar> replayResidual(const ResidualAudit& audit) {
    if(!audit.finished || audit.failed || audit.truncated || audit.rows<0 ||
       audit.final.size()!=audit.rows)
        throw std::invalid_argument("Incomplete residual audit");
    std::vector<Scalar> result(static_cast<std::size_t>(audit.rows),Scalar(0));
    for(const auto& op:audit.operations) {
        if(op.row<0 || op.row>=audit.rows || !std::isfinite(op.operand))
            throw std::invalid_argument("Invalid residual audit operand");
        auto& value=result[static_cast<std::size_t>(op.row)];
        switch(op.operation) {
        case '+':value+=Scalar(op.operand);break;
        case '-':value-=Scalar(op.operand);break;
        case '/':
            if(op.operand==0)throw std::invalid_argument("Zero residual audit divisor");
            value/=Scalar(op.operand);break;
        case '=':value=Scalar(op.operand);break;
        default:throw std::invalid_argument("Unknown residual audit operation");
        }
    }
    return result;
}

inline VectorXd residualAuditIncrement(const ResidualAudit& before,const ResidualAudit& after) {
    if(before.rows!=after.rows)throw std::invalid_argument("Residual audit size mismatch");
    using Wide=boost::multiprecision::cpp_dec_float_100;
    const auto a=replayResidual<Wide>(before),b=replayResidual<Wide>(after);
    VectorXd difference(before.rows);
    for(int i=0;i<before.rows;++i)difference(i)=static_cast<Real>(b[i]-a[i]);
    if(!difference.allFinite())throw std::runtime_error("Nonfinite residual increment");
    return difference;
}
}
