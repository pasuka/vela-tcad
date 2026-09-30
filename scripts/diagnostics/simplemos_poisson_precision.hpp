#pragma once
// Diagnostic only: preserve the frozen double geometry/material coefficients.
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <vector>
#include <stdexcept>
namespace simplemos_poisson_precision {
using Q = boost::multiprecision::cpp_bin_float_quad;
struct Node {
    double psi,n,p,xpsi,xn,xp,potentialScale,eref,href,ni,vt;
    double doping,voln,volp,vold,q,area,scale,interfaceRhs;
};
struct Edge { int i,j; double g; };
inline Q clippedExp(Q z) {
    if(z < -500)z=-500;
    if(z > 500)z=500;
    return exp(z);
}
inline std::vector<Q> evaluate(const std::vector<Node>& nodes,
                               const std::vector<Edge>& edges, bool packed) {
    std::vector<Q> psi(nodes.size()),n(nodes.size()),p(nodes.size()),r(nodes.size());
    for(std::size_t i=0;i<nodes.size();++i) {
        const auto& a=nodes[i];
        psi[i]=packed?Q(a.xpsi)*Q(a.potentialScale):Q(a.psi);
        n[i]=a.n;p[i]=a.p;
        if(packed) {
            n[i]=a.ni>0 ? Q(a.ni)*clippedExp((psi[i]-Q(a.eref)-Q(a.xn)*Q(a.potentialScale))/Q(a.vt)):Q(0);
            p[i]=a.ni>0 ? Q(a.ni)*clippedExp((Q(a.href)+Q(a.xp)*Q(a.potentialScale)-psi[i])/Q(a.vt)):Q(0);
        }
    }
    for(const auto& e:edges) {
        const Q f=Q(e.g)*(psi.at(e.i)-psi.at(e.j));
        r.at(e.i)+=f;r.at(e.j)-=f;
    }
    for(std::size_t i=0;i<nodes.size();++i) {
        const auto& a=nodes[i];
        r[i]=(r[i]+Q(a.q)*Q(a.area)*(n[i]*Q(a.voln)-p[i]*Q(a.volp)-Q(a.doping)*Q(a.vold))-Q(a.interfaceRhs))/Q(a.scale);
    }
    return r;
}
}
