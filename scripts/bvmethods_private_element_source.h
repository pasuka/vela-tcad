// Fixed-state diagnostic only. Not included by any production CMake target.
// Inputs: cm, cm^-3, V, cm^2/(V s), um^2. Output: G*Measure.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace bvmethods_diagnostic {
template<class T> struct ElementResult {
    std::array<T, 2> current{}, field{};
    std::array<T, 2> referenceDrop{}, incrementDrop{}, qfDrop{};
    T force{}, alpha{}, generation{};
};

template<class T> T bernoulli(T x) {
    using std::abs; using std::exp; using std::expm1;
    if (abs(x)<T(1e-5)) return 1-x/2+x*x/12-x*x*x*x/720;
    if (x>600) return x*exp(-x);
    if (x< -600) return -x;
    return x/expm1(x);
}
template<class T> T alpha(T field, bool electron) {
    using std::cos; using std::exp;
    if (field<=0) return 0;
    T u=(field/T(400000)-T(.95))/T(.1);
    T w=u<=0?T(0):(u>=1?T(1):(1-cos(T(3.141592654)*u))/2);
    const T a0=electron?703000:1582000, a1=electron?703000:671000;
    const T b0=electron?1231000:2036000, b1=electron?1231000:1693000;
    return (1-w)*a0*exp(-b0/field)+w*a1*exp(-b1/field);
}

// The representative is endpoint zero; the two other endpoints are ordered
// by the independently recorded native HFS callback identity, never fitted.
// QF operands are kept separate until subtraction in the selected precision.
template<class T> ElementResult<T> evaluate(
    const std::array<std::array<double,2>,3>& points,
    const std::array<double,3>& density,
    const std::array<double,3>& psi,
    const std::array<double,3>& reference,
    const std::array<double,3>& increment,
    double mobility, bool electron,
    const std::array<std::array<double,2>,2>& normals, int normalCount,
    bool contact)
{
    using std::sqrt; using std::log; using std::expm1;
    ElementResult<T> r;
    std::array<std::array<T,2>,2> edge{}, direction{};
    std::array<T,2> j{}, dp{};
    const T vt=T(.02585199526648486), q=T(1.602192e-19);
    const T sign=electron?1:-1;
    if (mobility<0 || !std::isfinite(mobility) || normalCount<0 || normalCount>2)
        throw std::invalid_argument("invalid diagnostic mobility/normals");
    for(int k=0;k<3;++k)
        if (!(density[k]>0) || !std::isfinite(density[k]) ||
            !std::isfinite(psi[k]) || !std::isfinite(reference[k]) ||
            !std::isfinite(increment[k]))
            throw std::invalid_argument("invalid diagnostic primary state");
    for(int k=0;k<2;++k) {
        for(int d=0;d<2;++d) edge[k][d]=T(points[k+1][d])-T(points[0][d]);
        const T h=sqrt(edge[k][0]*edge[k][0]+edge[k][1]*edge[k][1]);
        if (!(h>0)) throw std::invalid_argument("zero diagnostic edge");
        for(int d=0;d<2;++d) direction[k][d]=edge[k][d]/h;
        r.referenceDrop[k]=T(reference[k+1])-T(reference[0]);
        r.incrementDrop[k]=T(increment[k+1])-T(increment[0]);
        r.qfDrop[k]=r.referenceDrop[k]+r.incrementDrop[k];
        const T delta=sign*r.qfDrop[k]/vt;
        const T arg=log(T(density[k+1])/T(density[0]))+delta;
        const T value=delta>0 ? T(density[0])*bernoulli(-arg)*(-expm1(-delta))
                              : T(density[k+1])*bernoulli(arg)*expm1(delta);
        j[k]=-sign*q*vt*value/h;
        dp[k]=T(psi[k+1])-T(psi[0]);
    }
    const auto solve=[](const auto& a, const auto& b) {
        const T det=a[0][0]*a[1][1]-a[0][1]*a[1][0];
        if (det==0) throw std::invalid_argument("degenerate diagnostic triangle");
        return std::array<T,2>{(b[0]*a[1][1]-b[1]*a[0][1])/det,
                               (a[0][0]*b[1]-a[1][0]*b[0])/det};
    };
    r.current=solve(direction,j); r.field=solve(edge,dp);
    for(auto& f:r.field) f=-f;
    if (!contact && normalCount) {
        const T nx=T(normals[0][0]), ny=T(normals[0][1]);
        const T norm=nx*nx+ny*ny;
        if (!(norm>0)) throw std::invalid_argument("zero diagnostic normal");
        if (normalCount==2 && std::abs(nx*T(normals[1][1])-ny*T(normals[1][0]))>T(1e-12)) {
            r.current={0,0}; r.field={0,0};
        } else {
            for(auto* v:{&r.current,&r.field}) {
                const T dot=((*v)[0]*nx+(*v)[1]*ny)/norm;
                (*v)[0]-=dot*nx; (*v)[1]-=dot*ny;
            }
        }
    }
    const T mag=sqrt(r.current[0]*r.current[0]+r.current[1]*r.current[1]);
    if (mag>0) r.force=std::max(T(0),(r.field[0]*r.current[0]+r.field[1]*r.current[1])/mag);
    r.alpha=alpha(r.force,electron);
    r.generation=r.alpha*mag*T(mobility)/q;
    return r;
}
} // namespace bvmethods_diagnostic
