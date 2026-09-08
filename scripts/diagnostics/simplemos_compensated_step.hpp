#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
namespace simplemos_compensated_step {
struct Sum { double high; double low; };
inline Sum twoSum(double a,double b) {
    const double s=a+b,bv=s-a;
    return {s,(a-(s-bv))+(b-bv)};
}
struct Proposal { double value; double remainder; };
inline Proposal propose(double x,double step,double alpha,double carry) {
    const double product=alpha*step;
    const double productError=std::fma(alpha,step,-product);
    const auto increment=twoSum(product,carry);
    const auto state=twoSum(x,increment.high);
    return {state.high,state.low+increment.low+productError};
}
// Diagnostic executables perform one Newton solve per process. State is reset
// before each coupled solve and committed only after line-search acceptance.
extern bool active;
extern std::vector<double> carry;
}
