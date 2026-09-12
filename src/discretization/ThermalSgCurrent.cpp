#include "vela/discretization/ThermalSgCurrent.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/core/PhysicalConstants.h"
#include <cmath>
#include <stdexcept>
namespace vela {
namespace {
struct D {Real v=0.;std::array<Real,8> d{};D()=default;D(Real x):v(x){}};
D operator+(const D&a,const D&b){D r(a.v+b.v);for(int i=0;i<8;++i)r.d[i]=a.d[i]+b.d[i];return r;}
D operator-(const D&a,const D&b){D r(a.v-b.v);for(int i=0;i<8;++i)r.d[i]=a.d[i]-b.d[i];return r;}
D operator*(const D&a,const D&b){D r(a.v*b.v);for(int i=0;i<8;++i)r.d[i]=a.d[i]*b.v+a.v*b.d[i];return r;}
D operator/(const D&a,const D&b){D r(a.v/b.v);for(int i=0;i<8;++i)r.d[i]=(a.d[i]-r.v*b.d[i])/b.v;return r;}
D exp(const D&a){D r(std::exp(a.v));for(int i=0;i<8;++i)r.d[i]=r.v*a.d[i];return r;}
D expm1(const D&a){D r(std::expm1(a.v));for(int i=0;i<8;++i)r.d[i]=std::exp(a.v)*a.d[i];return r;}
D log(const D&a){D r(std::log(a.v));for(int i=0;i<8;++i)r.d[i]=a.d[i]/a.v;return r;}
D variable(Real x,int i){D r(x);r.d[i]=1.;return r;}
D quantity(const ThermalQuantity&q,int offset){D r(q.value);for(int i=0;i<4;++i)r.d[offset+i]=q.derivative[i];return r;}
D logB(const D&x){D r;Real slope;
    if(std::abs(x.v)<1e-4){const Real z=x.v*x.v;r.v=-x.v/2.-z/24.+z*z/2880.;slope=-.5-x.v/12.+x.v*z/720.;}
    else if(x.v>50.){r.v=std::log(x.v)-x.v-std::log1p(-std::exp(-x.v));slope=1./x.v-1./(1.-std::exp(-x.v));}
    else if(x.v< -50.){r.v=std::log(-x.v)-std::log1p(-std::exp(x.v));slope=1./x.v+std::exp(x.v)/(1.-std::exp(x.v));}
    else {r.v=std::log(x.v/std::expm1(x.v));slope=1./x.v-std::exp(x.v)/std::expm1(x.v);}
    for(int i=0;i<8;++i)r.d[i]=slope*x.d[i];return r;
}
D einsteinLimit(const D&eta){Real f=fermiDiracHalf(eta.v),df=fermiDiracHalfDerivative(eta.v),ddf=fermiDiracHalfSecondDerivative(eta.v);D r(df>0.?f/df:1.);
    if(r.v<1.)return D(1.);if(df>0.)for(int i=0;i<8;++i)r.d[i]=(1.-(f/df)*(ddf/df))*eta.d[i];return r;
}
}
ThermalSgCurrent thermalSgCurrent(const SiliconThermalState&a,const SiliconThermalResult&pa,
    const SiliconThermalState&b,const SiliconThermalResult&pb,Real mu,Real weight,bool electron){
    if(!std::isfinite(mu)||mu<=0.||!std::isfinite(weight)||weight<0.)throw std::invalid_argument("Invalid thermal SG mobility/weight");
    if(weight==0.)return {};
    const D ta=variable(a.temperature_K,3),tb=variable(b.temperature_K,7),vt=D(constants::kb/constants::q)*(ta+tb)/D(2.);
    const D qa=variable(electron?a.electronQf_V:a.holeQf_V,electron?1:2),qb=variable(electron?b.electronQf_V:b.holeQf_V,electron?5:6);
    const D n0=quantity(electron?pa.electrons_m3:pa.holes_m3,0),n1=quantity(electron?pb.electrons_m3:pb.holes_m3,4);
    const D dos0=quantity(electron?pa.Nc_m3:pa.Nv_m3,0),dos1=quantity(electron?pb.Nc_m3:pb.Nv_m3,4);
    if(!(ta.v>=50.)||!(tb.v>=50.)||!(n0.v>0.)||!(n1.v>0.))throw std::invalid_argument("Invalid thermal SG temperature/density");
    const D eta0=quantity(electron?pa.electronEta:pa.holeEta,0);
    const D eta1=quantity(electron?pb.electronEta:pb.holeEta,4);
    const D logF=(log(n1)-log(dos1))-(log(n0)-log(dos0));
    D g=std::abs(logF.v)>1e-8?(eta1-eta0)/logF:einsteinLimit((eta0+eta1)/D(2.));
    if(!(g.v>0.)||!std::isfinite(g.v))g=einsteinLimit((eta0+eta1)/D(2.));
    const D qfDifference=D(electron?b.electronQfReference_V-a.electronQfReference_V:b.holeQfReference_V-a.holeQfReference_V)+qb-qa;
    const D y=D(electron?1.:-1.)*qfDifference/(vt*g),arg=log(n1)-log(n0)+y;
    D current;
    if(std::abs(y.v)<1e-5)current=exp(log(vt*g*n1)+logB(arg))*expm1(y)*D((electron?-1.:1.)*constants::q*mu*weight);
    else {
        D logRelative;
        logRelative.v=y.v>50.?y.v+std::log1p(-std::exp(-y.v)):std::log(std::abs(std::expm1(y.v)));
        const Real slope=y.v>50.?1./(1.-std::exp(-y.v)):std::exp(y.v)/std::expm1(y.v);
        for(int i=0;i<8;++i)logRelative.d[i]=slope*y.d[i];
        current=exp(log(vt*g*n1)+logB(arg)+logRelative)*D((electron?-1.:1.)*(y.v>0.?1.:-1.)*constants::q*mu*weight);
    }
    ThermalSgCurrent result{current.v,current.d,current.v/mu};
    if(!std::isfinite(current.v))throw std::runtime_error("Nonfinite thermal SG current");
    for(Real d:current.d)if(!std::isfinite(d))throw std::runtime_error("Nonfinite thermal SG derivative");
    return result;
}
} // namespace vela
