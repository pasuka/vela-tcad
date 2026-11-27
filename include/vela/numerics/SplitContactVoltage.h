#pragma once
#include "vela/numerics/SplitCoordinate.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
namespace vela::numerics {
// Contact voltage and load-line arithmetic share this pair. Never collapse
// hi+lo to double before constructing the boundary target or circuit residual.
struct SplitContactVoltage {
    SplitCoordinate value;
    SplitContactVoltage shifted(double alpha,double direction) const {
        const double product=alpha*direction;
        return {value.shifted(product).shifted(std::fma(alpha,direction,-product))};
    }
    double loadLine(double resistance,double current,double supply) const {
        using W=boost::multiprecision::cpp_bin_float_quad;
        if(!(resistance>0) || !std::isfinite(resistance) || !std::isfinite(current) ||
           !std::isfinite(supply) || !std::isfinite(value.hi) || !std::isfinite(value.lo))
            throw std::invalid_argument("Invalid split contact voltage");
        // Cancel the high voltage first. Adding its tail before subtracting
        // an equal supply would unnecessarily round even a quad expression.
        return static_cast<double>((W(value.hi)-W(supply))+W(resistance)*W(current)+W(value.lo));
    }
};
}
