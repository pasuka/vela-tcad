#include "vela/equation/LatticeBandEdgeWork.h"
#include <cmath>
#include <stdexcept>
namespace vela {
LatticeBandEdgeWork latticeBandEdgeWork(Real jn,Real jp,Real ec0,Real ec1,Real ev0,Real ev1){
    for(Real v:{jn,jp,ec0,ec1,ev0,ev1})if(!std::isfinite(v))throw std::invalid_argument("Nonfinite electrical heat input");
    const Real dn=ec1-ec0,dp=ev1-ev0,power=jn*dn+jp*dp;
    if(!std::isfinite(power))throw std::runtime_error("Nonfinite electrical heat output");
    return {power,{dn,dp,-jn,jn,-jp,jp}};
}
} // namespace vela
