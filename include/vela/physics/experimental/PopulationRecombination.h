#pragma once
#include "vela/equation/Tri3LocalForwardAD.h"
#include <cmath>

namespace vela::experimental {
// Population coordinates: no quasi-Fermi reconstruction or potential column.
// Compensated wide products retain cancellation near equilibrium, while the
// wider exponent range avoids overflowing a representable product difference.
inline detail::Tri3LocalForwardDual populationRecombinationExcess(
    const detail::Tri3LocalForwardDual& n,
    const detail::Tri3LocalForwardDual& p,
    const detail::Tri3LocalForwardDual& gammaN,
    const detail::Tri3LocalForwardDual& gammaP,double ni)
{
    using Wide=long double;
    Wide q=Wide(gammaN.value)*gammaP.value;
    Wide qError=std::fma(Wide(gammaN.value),Wide(gammaP.value),-q);
    for(Wide factor:{Wide(ni),Wide(ni)}) {
        const Wide next=q*factor;
        qError=std::fma(q,factor,-next)+qError*factor;
        q=next;
    }
    const Wide product=Wide(n.value)*p.value;
    const Wide difference=product-q;
    const Wide tail=difference-product;
    const Wide subtractionError=(product-(difference-tail))+(-q-tail);
    const Wide productError=std::fma(Wide(n.value),Wide(p.value),-product);
    detail::Tri3LocalForwardDual result(static_cast<double>(
        difference+(subtractionError+productError-qError)));
    for(std::size_t k=0;k<result.derivative.size();++k) {
        const Wide carrier=Wide(p.value)*n.derivative[k]+Wide(n.value)*p.derivative[k];
        const Wide equilibrium=Wide(ni)*ni*(Wide(gammaP.value)*gammaN.derivative[k]
                                            +Wide(gammaN.value)*gammaP.derivative[k]);
        result.derivative[k]=static_cast<double>(carrier-equilibrium);
    }
    return result;
}
} // namespace vela::experimental
