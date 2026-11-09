#pragma once

#include "vela/equation/Tri3LocalForwardAD.h"
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela::experimental {

// Local contract for the explicit bvmethods_tri3_qf profile and diagnostics.
// Units: cm, cm^-3, V, cm^2/(V s); outputs are particle rates, without q.
// Local derivative order is psi[0:3], fn[0:3], fp[0:3]. Density, mobility
// and driving-field jets must come from the caller's same material/state law.
using AvJet = detail::Tri3LocalForwardDual;

inline AvJet avLog(const AvJet& x)
{
    if (!(x.value > 0.)) throw std::domain_error("AV log requires positive density");
    AvJet y(std::log(x.value));
    for (std::size_t k=0;k<y.derivative.size();++k) y.derivative[k]=x.derivative[k]/x.value;
    return y;
}

inline AvJet avLogRatio(const AvJet& numerator,const AvJet& denominator)
{
    // Keep the qualified near-equilibrium arithmetic. Outside the safe
    // ratio range, logs and their derivatives remain representable even
    // when forming the density quotient (or its derivative) does not.
    const double ratio=numerator.value/denominator.value;
    if(ratio>=std::sqrt(std::numeric_limits<double>::min()) &&
       ratio<=std::sqrt(std::numeric_limits<double>::max()))
    {
        const auto quotient=numerator/denominator;
        bool finite=true;
        for(double d:quotient.derivative)finite=finite&&std::isfinite(d);
        if(finite)return avLog(quotient);
        // The quotient can be ordinary while its intermediate tangent
        // overflows. d log(a/b)=da/a-db/b remains representable.
    }
    return avLog(numerator)-avLog(denominator);
}

inline AvJet avNorm(const std::array<AvJet,2>& v)
{
    AvJet y(std::hypot(v[0].value,v[1].value));
    // Zero is an explicit generalized-derivative convention, NOT a claim
    // that the Euclidean norm has a Frechet derivative at the origin.
    if (y.value>0.) for(std::size_t k=0;k<y.derivative.size();++k)
        y.derivative[k]=(v[0].value/y.value)*v[0].derivative[k]
                       +(v[1].value/y.value)*v[1].derivative[k];
    return y;
}

inline double avNormDirectional(const std::array<double,2>& v,
                                const std::array<double,2>& direction)
{
    const double n=std::hypot(v[0],v[1]);
    return n>0. ? v[0]/n*direction[0]+v[1]/n*direction[1]
                : std::hypot(direction[0],direction[1]);
}

inline AvJet avBernoulli(const AvJet& x)
{
    if(std::abs(x.value)<1e-5) {
        const auto z=x*x;
        return AvJet(1)-x/AvJet(2)+z/AvJet(12)-z*z/AvJet(720)+z*z*z/AvJet(30240);
    }
    if(x.value>50.) {
        const auto e=detail::dualExp(-x);
        return x*e/(AvJet(1)-e);
    }
    if(x.value<-50.) return -x/(AvJet(1)-detail::dualExp(x));
    return x/detail::dualExpm1(x);
}

inline AvJet avSgFlux(const AvJet& coefficient,const AvJet& from,const AvJet& to,
                     const AvJet& delta,const AvJet& argument)
{
    // B(-a)=exp(a)B(a), a=log(to/from)+delta. The two Bernoulli
    // terms cannot nearly cancel for |delta|>50, and avoid both an
    // overflowing expm1(delta) and loss of the tiny endpoint derivative.
    if(std::abs(delta.value)>50.)
        return coefficient*(from*avBernoulli(-argument)-to*avBernoulli(argument));
    return coefficient*to*avBernoulli(argument)*detail::dualExpm1(delta);
}

struct AvCoefficient {
    double lowA=0.,lowB=1.,highA=0.,highB=1.,switchField=4e5;
};

inline AvJet avExponentialCoefficient(const AvJet& field,double a,double b)
{
    if(a==0.)return {};
    auto result=AvJet(a)*detail::dualExp(-AvJet(b)/field);
    bool finite=std::isfinite(result.value);
    for(double d:result.derivative)finite=finite&&std::isfinite(d);
    if(finite)return result;
    // A exp(-B/E) can have a finite zero value AND zero tangent although
    // differentiating -B/E first overflows. Keep the combined exponent.
    const long double exponent=-(long double)b/field.value;
    const long double value=(long double)a*std::exp(exponent);
    result.value=static_cast<double>(value);
    const long double slope=std::exp(std::log((long double)a)+
        std::log((long double)b)-2*std::log((long double)field.value)+exponent);
    for(std::size_t k=0;k<result.derivative.size();++k)
        result.derivative[k]=static_cast<double>(slope*field.derivative[k]);
    return result;
}

inline AvJet avCoefficient(const AvJet& field,const AvCoefficient& p)
{
    if(field.value<=0.) return {};
    const auto low=avExponentialCoefficient(field,p.lowA,p.lowB);
    if(field.value<=.95*p.switchField) return low;
    const auto high=avExponentialCoefficient(field,p.highA,p.highB);
    if(field.value>=1.05*p.switchField) return high;
    const auto t=(field/AvJet(p.switchField)-AvJet(.95))/AvJet(.1);
    AvJet weight(.5*(1-std::cos(std::acos(-1.)*t.value)));
    for(std::size_t k=0;k<weight.derivative.size();++k)
        weight.derivative[k]=.5*std::acos(-1.)*std::sin(std::acos(-1.)*t.value)*t.derivative[k];
    return (AvJet(1)-weight)*low+weight*high;
}

struct AvTri3Geometry {
    std::array<std::array<double,2>,3> points_cm{};
    std::array<double,3> boxAreas_cm2{};
    std::size_t representative=0;
    std::array<double,4> currentProjection{1.,0.,0.,1.};
    std::array<std::size_t,3> nodeOrder{0,1,2};
};

inline std::array<AvJet,2> avParticleCurrentWithoutMobility(
    const AvTri3Geometry& g,const std::array<AvJet,3>& density,
    const std::array<AvJet,3>& qf,double thermalVoltage,bool electron)
{
    if(g.representative>=3 || !(thermalVoltage>0.))
        throw std::invalid_argument("Invalid AV representative or thermal voltage");
    for(const auto& n:density) if(!(n.value>0.) || !std::isfinite(n.value))
        throw std::domain_error("AV requires finite positive populations");
    const auto r=g.representative;
    std::array<std::array<double,2>,2> t{};
    std::array<AvJet,2> j{};std::size_t edge=0;
    const double sign=electron?1.:-1.;
    for(std::size_t i=0;i<3;++i)if(i!=r){
        const double dx=g.points_cm[i][0]-g.points_cm[r][0];
        const double dy=g.points_cm[i][1]-g.points_cm[r][1];
        const double length=std::hypot(dx,dy);
        if(!(length>0.))throw std::invalid_argument("Degenerate AV edge");
        t[edge]={dx/length,dy/length};
        const auto delta=AvJet(sign/thermalVoltage)*(qf[i]-qf[r]);
        const auto argument=avLogRatio(density[i],density[r])+delta;
        j[edge]=avSgFlux(AvJet(-sign*thermalVoltage/length),density[r],density[i],delta,argument);
        ++edge;
    }
    const double determinant=t[0][0]*t[1][1]-t[0][1]*t[1][0];
    if(std::abs(determinant)<1e-14)throw std::invalid_argument("Degenerate AV triangle");
    const auto x=(j[0]*AvJet(t[1][1])-AvJet(t[0][1])*j[1])/AvJet(determinant);
    const auto y=(AvJet(t[0][0])*j[1]-j[0]*AvJet(t[1][0]))/AvJet(determinant);
    return {AvJet(g.currentProjection[0])*x+AvJet(g.currentProjection[1])*y,
            AvJet(g.currentProjection[2])*x+AvJet(g.currentProjection[3])*y};
}

struct AvTri3Result {
    AvJet generation_cm3_s;
    std::array<AvJet,3> nodeGeneration_per_cm_s;
    std::array<AvJet,2> electronParticleCurrentWithoutMobility;
    std::array<AvJet,2> holeParticleCurrentWithoutMobility;
    bool electronZeroNorm=false,holeZeroNorm=false;
};

inline AvTri3Result evaluateBvMethodsTri3Avalanche(
    const AvTri3Geometry& geometry,const std::array<AvJet,3>& n,
    const std::array<AvJet,3>& p,const std::array<AvJet,3>& fn,
    const std::array<AvJet,3>& fp,const AvJet& mun,const AvJet& mup,
    const AvJet& electronField,const AvJet& holeField,double thermalVoltage,
    const AvCoefficient& electron,const AvCoefficient& hole,
    const std::array<AvJet,2>* directElectronCurrent=nullptr,
    const std::array<AvJet,2>* directHoleCurrent=nullptr)
{
    AvTri3Result out;
    out.electronParticleCurrentWithoutMobility=
        directElectronCurrent?*directElectronCurrent:avParticleCurrentWithoutMobility(geometry,n,fn,thermalVoltage,true);
    out.holeParticleCurrentWithoutMobility=
        directHoleCurrent?*directHoleCurrent:avParticleCurrentWithoutMobility(geometry,p,fp,thermalVoltage,false);
    const auto jn=avNorm(out.electronParticleCurrentWithoutMobility);
    const auto jp=avNorm(out.holeParticleCurrentWithoutMobility);
    out.electronZeroNorm=jn.value==0.;out.holeZeroNorm=jp.value==0.;
    out.generation_cm3_s=avCoefficient(electronField,electron)*mun*jn
                       +avCoefficient(holeField,hole)*mup*jp;
    for(std::size_t k=0;k<3;++k){
        if(!(geometry.boxAreas_cm2[k]>=0.))throw std::invalid_argument("Negative AV box measure");
        out.nodeGeneration_per_cm_s[k]=AvJet(geometry.boxAreas_cm2[k])*out.generation_cm3_s;
    }
    return out;
}

} // namespace vela::experimental
