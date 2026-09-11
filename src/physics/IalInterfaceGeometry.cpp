#include "vela/physics/IalInterfaceGeometry.h"
#include "vela/physics/IalMobility.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace vela {
std::vector<IalInterfaceNode> buildIalInterfaceGeometry(
    const std::vector<std::array<Real,2>>& coordinates_m,
    const std::vector<IalInterfaceSegment>& segments,
    const std::array<Real,3>& crystalX,
    const std::array<Real,3>& crystalY)
{
    if (coordinates_m.empty()||segments.empty())
        throw std::invalid_argument("IALMob geometry requires nodes and interface segments");
    Real xx=0.,yy=0.,xy=0.;
    for (int k=0;k<3;++k) {
        if (!std::isfinite(crystalX[k])||!std::isfinite(crystalY[k]))
            throw std::invalid_argument("Invalid IALMob crystal axes");
        xx+=crystalX[k]*crystalX[k];yy+=crystalY[k]*crystalY[k];xy+=crystalX[k]*crystalY[k];
    }
    if (std::abs(xx-1.)>1e-12||std::abs(yy-1.)>1e-12||std::abs(xy)>1e-12)
        throw std::invalid_argument("IALMob crystal axes must be orthonormal");
    for (const auto& p:coordinates_m)
        for (Real x:p)
            if (!std::isfinite(x)) throw std::invalid_argument("Invalid IALMob mesh coordinate");
    std::map<std::size_t,std::array<Real,2>> normal;
    std::set<std::pair<std::size_t,std::size_t>> seen;
    for (const auto& s:segments) {
        if (s.node0>=coordinates_m.size()||s.node1>=coordinates_m.size()||s.node0==s.node1)
            throw std::invalid_argument("Invalid IALMob interface endpoint");
        if (!seen.insert(std::minmax(s.node0,s.node1)).second)
            throw std::invalid_argument("Duplicate IALMob interface segment");
        const auto& a=coordinates_m[s.node0];const auto& b=coordinates_m[s.node1];
        const Real dx=b[0]-a[0],dy=b[1]-a[1],length=std::hypot(dx,dy);
        if (!std::isfinite(length)||length==0.) throw std::invalid_argument("Degenerate IALMob interface");
        for (Real x:s.semiconductorPoint_m)
            if (!std::isfinite(x)) throw std::invalid_argument("Invalid IALMob semiconductor point");
        Real nx=-dy/length,ny=dx/length;
        const Real side=nx*(s.semiconductorPoint_m[0]-a[0])+ny*(s.semiconductorPoint_m[1]-a[1]);
        if (side==0.) throw std::invalid_argument("IALMob semiconductor point must select a side");
        if (side<0.) {nx=-nx;ny=-ny;}
        for (std::size_t node:{s.node0,s.node1}) {normal[node][0]+=nx;normal[node][1]+=ny;}
    }
    std::map<std::size_t,int> family;
    for (const auto& [node,n]:normal) {
        std::array<Real,3> crystal{};
        for (int k=0;k<3;++k) crystal[k]=n[0]*crystalX[k]+n[1]*crystalY[k];
        family[node]=IalMobility::orientationFamily(crystal);
    }
    std::vector<IalInterfaceNode> result;
    result.reserve(coordinates_m.size());
    for (std::size_t node=0;node<coordinates_m.size();++node) {
        const auto& p=coordinates_m[node];
        Real distance=std::numeric_limits<Real>::infinity();
        for (const auto& s:segments) {
            const auto& a=coordinates_m[s.node0];const auto& b=coordinates_m[s.node1];
            const Real dx=b[0]-a[0],dy=b[1]-a[1],length=std::hypot(dx,dy);
            const Real tx=dx/length,ty=dy/length;
            const Real t=std::clamp((p[0]-a[0])*tx+(p[1]-a[1])*ty,Real(0.),length);
            distance=std::min(distance,std::hypot(p[0]-a[0]-t*tx,p[1]-a[1]-t*ty));
        }
        Real nearest=std::numeric_limits<Real>::infinity();std::size_t selected=0;
        for (const auto& [candidate,f]:family) {
            const auto& q=coordinates_m[candidate];
            const Real d=std::hypot(p[0]-q[0],p[1]-q[1]);
            if (d<nearest) {nearest=d;selected=candidate;}
        }
        result.push_back({distance,family.at(selected),selected,family.count(node)!=0});
    }
    return result;
}
} // namespace vela
