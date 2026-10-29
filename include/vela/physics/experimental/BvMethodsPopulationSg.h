#pragma once

#include "vela/physics/experimental/BvMethodsTri3Avalanche.h"
#include "vela/physics/experimental/JoyceDixonStatistics.h"

namespace vela::experimental {

// Experimental population-coordinate arithmetic. The primary density is never
// reconstructed through QF. Pairing B(-a)=B(a)+a keeps edge reversal consistent.
inline std::array<AvJet,2> populationBernoulliPair(const AvJet& x)
{
    const auto a=x.value<0.?-x:x;
    AvJet positive;
    if(a.value>10.) {
        const auto e=detail::dualExp(-a);
        positive=a*e/(AvJet(1)-e);
    } else if(a.value>.5) {
        positive=AvJet(1)/((detail::dualExp(a)-AvJet(1))/a);
    } else {
        AvJet term(1),sum(1);double d=-.5;
        for(int k=0;k<1000;++k) {
            term=term*(a/AvJet(k+1));
            const auto next=sum+term/AvJet(k+2);
            const double nd=d-term.value/(k+3);
            const bool stop=next.value==sum.value && nd==d;
            sum=next;d=nd;if(stop)break;
        }
        positive=AvJet(1)/sum;
    }
    const auto negative=positive+a;
    return x.value<0.?std::array<AvJet,2>{negative,positive}
                     :std::array<AvJet,2>{positive,negative};
}

inline AvJet populationLogGamma(const AvJet& density,double dos)
{
    AvJet eta(JoyceDixonStatistics::inverseDensity(density.value,dos));
    const double slope=JoyceDixonStatistics::densitySlopeFromPopulation(density.value,dos);
    for(std::size_t k=0;k<eta.derivative.size();++k)eta.derivative[k]=density.derivative[k]/slope;
    return avLogRatio(density,AvJet(dos))-eta;
}

inline AvJet populationSgArgument(const AvJet& from,const AvJet& to,
    const AvJet& psiFrom,const AvJet& psiTo,double niFrom,double niTo,
    double dosFrom,double dosTo,double vt,bool electron,const AvJet& einstein=AvJet(1))
{
    if(!(vt>0.) || !(einstein.value>0.))throw std::invalid_argument("Invalid population SG scale");
    const auto band=AvJet(electron?1.:-1.)*(psiTo-psiFrom)/AvJet(vt)
                   +avLogRatio(AvJet(niTo),AvJet(niFrom));
    const auto gamma=populationLogGamma(to,dosTo)-populationLogGamma(from,dosFrom);
    return (band+gamma+(einstein-AvJet(1))*avLogRatio(to,from))/einstein;
}

inline AvJet populationSgFlux(const AvJet& coefficient,const AvJet& from,
                              const AvJet& to,const AvJet& argument)
{
    const auto b=populationBernoulliPair(argument);
    return coefficient*(from*b[1]-to*b[0]);
}

inline std::array<AvJet,2> populationAvCurrent(const AvTri3Geometry& g,
    const std::array<AvJet,3>& density,const std::array<AvJet,3>& psi,
    const std::array<double,3>& ni,const std::array<double,3>& dos,double vt,bool electron)
{
    if(g.representative>=3)throw std::invalid_argument("Invalid population AV representative");
    std::array<std::array<double,2>,2> t{};std::array<AvJet,2> j{};
    std::size_t edge=0;const auto r=g.representative;
    for(std::size_t i=0;i<3;++i)if(i!=r) {
        const double dx=g.points_cm[i][0]-g.points_cm[r][0],dy=g.points_cm[i][1]-g.points_cm[r][1];
        const double length=std::hypot(dx,dy);
        if(!(length>0.))throw std::invalid_argument("Degenerate population AV edge");
        t[edge]={dx/length,dy/length};
        auto u=r,v=i;const bool reverse=g.nodeOrder[u]>g.nodeOrder[v];if(reverse)std::swap(u,v);
        const auto a=populationSgArgument(density[u],density[v],psi[u],psi[v],ni[u],ni[v],dos[u],dos[v],vt,electron);
        j[edge]=populationSgFlux(AvJet((electron?-1.:1.)*vt/length),density[u],density[v],a);
        if(reverse)j[edge]=-j[edge];++edge;
    }
    const double determinant=t[0][0]*t[1][1]-t[0][1]*t[1][0];
    if(std::abs(determinant)<1e-14)throw std::invalid_argument("Degenerate population AV triangle");
    const auto x=(j[0]*AvJet(t[1][1])-AvJet(t[0][1])*j[1])/AvJet(determinant);
    const auto y=(AvJet(t[0][0])*j[1]-j[0]*AvJet(t[1][0]))/AvJet(determinant);
    return {AvJet(g.currentProjection[0])*x+AvJet(g.currentProjection[1])*y,
            AvJet(g.currentProjection[2])*x+AvJet(g.currentProjection[3])*y};
}
} // namespace vela::experimental
