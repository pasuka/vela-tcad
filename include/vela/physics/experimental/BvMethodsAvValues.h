#pragma once

#include "vela/physics/experimental/BvMethodsTri3Transport.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {

// Ordinary double coefficient arithmetic is identical while the dual
// exponential and its tangent are finite. Exclude parameter and field
// ranges where the dual can take its extended-precision recovery.
inline bool avCoefficientsSafe(const AvCoefficient& p)
{
    return p.lowA>=0. && p.lowA<=1e30 && p.highA>=0. && p.highA<=1e30 &&
           p.lowB>0. && p.lowB<=1e30 && p.highB>0. && p.highB<=1e30 &&
           p.switchField>=1e-30 && p.switchField<=1e30;
}

inline bool avFieldSafe(double field)
{
    return field==0. || (field>=1e-60 && field<=1e30);
}

inline double avCoefficientValue(double field,const AvCoefficient& p)
{
    if(field<=0.)return 0.;
    const auto exponential=[&](double a,double b) {
        return a==0.?0.:a*std::exp(-b/field);
    };
    const double low=exponential(p.lowA,p.lowB);
    if(field<=.95*p.switchField)return low;
    const double high=exponential(p.highA,p.highB);
    if(field>=1.05*p.switchField)return high;
    const double t=(field/p.switchField-.95)/.1;
    const double weight=.5*(1-std::cos(std::acos(-1.)*t));
    return (1-weight)*low+weight*high;
}

inline std::array<double,2> avParticleCurrentValues(
    const AvTri3Geometry& g,const std::array<double,3>& density,
    const std::array<double,3>& qf,const std::array<double,3>& tangent,
    double thermalVoltage,bool electron)
{
    if(g.representative>=3 || !(thermalVoltage>0.))
        throw std::invalid_argument("Invalid AV representative or thermal voltage");
    for(double n:density)if(!(n>0.) || !std::isfinite(n))
        throw std::domain_error("AV requires finite positive populations");
    const auto r=g.representative;
    std::array<std::array<double,2>,2> t{};
    std::array<double,2> j{};
    std::size_t edge=0;
    const double sign=electron?1.:-1.;
    for(std::size_t i=0;i<3;++i)if(i!=r) {
        const double dx=g.points_cm[i][0]-g.points_cm[r][0];
        const double dy=g.points_cm[i][1]-g.points_cm[r][1];
        const double length=std::hypot(dx,dy);
        if(!(length>0.))throw std::invalid_argument("Degenerate AV edge");
        t[edge]={dx/length,dy/length};
        const double delta=(sign/thermalVoltage)*(qf[i]-qf[r]);
        const double argument=populationTransportLogRatio(density[i],density[r],tangent[i],tangent[r])+delta;
        const double coefficient=-sign*thermalVoltage/length;
        j[edge] = detail::sgBalancedEdgeFlux(coefficient, density[r], density[i],
            delta, delta, argument, tri3BernoulliValue, [](double x) { return std::expm1(x); });
        ++edge;
    }
    const double determinant=t[0][0]*t[1][1]-t[0][1]*t[1][0];
    if(std::abs(determinant)<1e-14)throw std::invalid_argument("Degenerate AV triangle");
    const double x=(j[0]*t[1][1]-t[0][1]*j[1])/determinant;
    const double y=(t[0][0]*j[1]-j[0]*t[1][0])/determinant;
    return {g.currentProjection[0]*x+g.currentProjection[1]*y,
            g.currentProjection[2]*x+g.currentProjection[3]*y};
}

struct AvGenerationValues {
    std::array<double,3> node{};
    double electron=0.,hole=0.;
    bool electronZero=false,holeZero=false;
    bool electronActive=false,holeActive=false;
};

inline AvGenerationValues avGenerationValues(
    const AvTri3Geometry& geometry,const std::array<double,2>& jn,
    const std::array<double,2>& jp,double mn,double mp,double en,double hp,
    const AvCoefficient& electron,const AvCoefficient& hole)
{
    AvGenerationValues out;
    const double n=std::hypot(jn[0],jn[1]);
    const double p=std::hypot(jp[0],jp[1]);
    const double an=avCoefficientValue(en,electron),ap=avCoefficientValue(hp,hole);
    out.electron=an*mn*n;
    out.hole=ap*mp*p;
    out.electronZero=n==0.;
    out.holeZero=p==0.;
    out.electronActive=out.electronZero&&an>0.;
    out.holeActive=out.holeZero&&ap>0.;
    const double total=out.electron+out.hole;
    for(std::size_t k=0;k<3;++k) {
        if(!(geometry.boxAreas_cm2[k]>=0.))throw std::invalid_argument("Negative AV box measure");
        out.node[k]=geometry.boxAreas_cm2[k]*total;
    }
    return out;
}

} // namespace vela::experimental
