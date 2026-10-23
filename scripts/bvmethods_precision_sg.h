// Private arithmetic control: higher precision cannot recover rounded inputs.
#pragma once
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <stdexcept>
#include <cmath>
namespace bvmethods_live {
inline long double decimalSg(double left,double right,double rightPopulation,
    long double qfDrop,double vt,double coefficient,bool electron) {
    using T=boost::multiprecision::cpp_dec_float_100;
    if(!(left>0&&right>0&&rightPopulation>0&&vt>0)||!std::isfinite(qfDrop))throw std::invalid_argument("invalid decimal SG state");
    if(qfDrop==0 || coefficient==0)return 0;
    const T delta=(electron?T(1):T(-1))*T(qfDrop)/T(vt);
    const T arg=log(T(right)/T(left))+delta;
    const T b=arg==0?T(1):T(arg/expm1(arg));
    return (T(coefficient)*b*T(rightPopulation)*expm1(delta)).convert_to<long double>();
}
}
