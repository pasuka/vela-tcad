#pragma once
// Extended arithmetic for Poisson; geometry and material inputs retain their stored precision.
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <vector>
#include <stdexcept>
#include <utility>
namespace vela::poisson_precision {
using Q = boost::multiprecision::cpp_bin_float_quad;
struct Node {
    // Match the production carrier mapping, including material reference offsets.
    // These are separate from the intrinsic density used by recombination.
    double psi,n,p,xpsi,xn,xp,potentialScale,eref,href,niElectron,niHole,vt;
    double doping,voln,volp,vold,q,area,scale,interfaceRhs;
};
struct Edge { int i,j; double g; };
// A read-only snapshot of the production Poisson arithmetic. Carrier QFs stay
// fixed; this is not a complete split drift-diffusion state or acceptance mode.
struct Model {
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<std::pair<int,double>> dirichlet;
};
inline Q clippedExp(Q z) {
    if(z < -500)z=-500;
    if(z > 500)z=500;
    return exp(z);
}
inline std::vector<Q> evaluate(const std::vector<Node>& nodes,
                               const std::vector<Edge>& edges, bool packed,
                               const std::vector<double>& packedPsiLow = {}) {
    if (!packedPsiLow.empty() && (!packed || packedPsiLow.size()!=nodes.size()))
        throw std::invalid_argument("Poisson low coordinates require one packed value per node");
    std::vector<Q> psi(nodes.size()),n(nodes.size()),p(nodes.size()),r(nodes.size());
    for(std::size_t i=0;i<nodes.size();++i) {
        const auto& a=nodes[i];
        psi[i]=packed?(Q(a.xpsi)+(packedPsiLow.empty()?Q(0):Q(packedPsiLow[i])))*Q(a.potentialScale):Q(a.psi);
        n[i]=a.n;p[i]=a.p;
        if(packed) {
            n[i]=a.niElectron>0 ? Q(a.niElectron)*clippedExp((psi[i]-Q(a.eref)-Q(a.xn)*Q(a.potentialScale))/Q(a.vt)):Q(0);
            p[i]=a.niHole>0 ? Q(a.niHole)*clippedExp((Q(a.href)+Q(a.xp)*Q(a.potentialScale)-psi[i])/Q(a.vt)):Q(0);
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
inline std::vector<Q> evaluate(const Model& model,
                              const std::vector<double>& packedPsiLow = {}) {
    auto r=evaluate(model.nodes,model.edges,true,packedPsiLow);
    for(const auto& [i,value]:model.dirichlet) {
        if(i<0 || static_cast<std::size_t>(i)>=model.nodes.size())
            throw std::invalid_argument("Poisson boundary node is out of range");
        r[i]=Q(model.nodes[i].xpsi)+(packedPsiLow.empty()?Q(0):Q(packedPsiLow[i]))-Q(value);
    }
    return r;
}
}
