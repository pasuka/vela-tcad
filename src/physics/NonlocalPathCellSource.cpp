#include "vela/physics/NonlocalPathCellSource.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace vela::experimental {
namespace {
constexpr Real tol=128*std::numeric_limits<Real>::epsilon();
void check(bool ok,const char* text) {if(!ok) throw std::invalid_argument(text);}
bool positive(Real v) {return std::isfinite(v) && v>0;}
Real softplus(Real v) {return std::max(v,Real(0))+std::log1p(std::exp(-std::abs(v)));}
Real sample(const TunnelingEndpoint& e,const std::vector<Real>& v) {
    long double out=0;for(const auto& [n,w]:e.weights) out+=static_cast<long double>(w)*v.at(n);
    return static_cast<Real>(out);
}
struct Profile {
    std::vector<TunnelingBandKnot> bands;
    std::vector<Real> ec,ev,coefficient,efn,efp;
};
Real interpolate(const Profile& p,const std::vector<Real>& v,Real x) {
    for(Index i=1;i<v.size();++i) if(p.bands[i-1].distance_m<=x && x<=p.bands[i].distance_m) {
        const Real w=(x-p.bands[i-1].distance_m)/(p.bands[i].distance_m-p.bands[i-1].distance_m);
        return v[i-1]+w*(v[i]-v[i-1]);
    }
    throw std::invalid_argument("turning point outside profile");
}
Real turning(const Profile& p,const std::vector<Real>& v,Real energy) {
    for(Index i=1;i<v.size();++i) if(v[i-1]>=energy && energy>=v[i] && v[i-1]!=v[i])
        return p.bands[i-1].distance_m+(p.bands[i].distance_m-p.bands[i-1].distance_m)*(v[i-1]-energy)/(v[i-1]-v[i]);
    throw std::invalid_argument("missing downward energy turning point");
}
// Physical scalar window integral [m^-2 s^-1]. Spatial quadrature is applied
// afterwards using corrected source-cell area / source chord.
Real scalar(const Profile& p,Real lower,Real upper,Real globalEv,Real globalEc,
    const NonlocalPathBandToBandConfig& c) {
    const auto constants=tunnelingPhysicalConstants(c.constants);const auto mass=resolveIndirectParameters(c);
    const Real mc=mass.conductionMass_kg,mv=mass.valenceMass_kg,mr=mc/(1+mc/mv);
    const Real h=constants.planck_J_s,hb=h/(2*std::numbers::pi_v<Real>);
    const Real kt=constants.boltzmann_J_per_K*p.bands[0].temperature_K;
    const Real energy=lower/2+upper/2;
    const Real xs=turning(p,p.ev,energy),xe=turning(p,p.ec,energy),length=xe-xs;
    check(positive(length),"nonpositive physical tunneling length");
    Index end=0;while(end<p.ec.size() && p.ec[end]>energy) ++end;
    check(end>0 && end<p.ec.size(),"missing first conduction turn");
    const auto a=tunnelingFractionIntegrals(std::vector<TunnelingBandKnot>(p.bands.begin(),p.bands.begin()+end+1),energy,mc,mv,hb);
    long double gap=0;
    for(Index i=1;i<p.bands.size();++i) {
        const Real lo=std::max(xs,p.bands[i-1].distance_m),hi=std::min(xe,p.bands[i].distance_m);
        if(lo<hi) gap+=static_cast<long double>(hi-lo)*(interpolate(p,p.coefficient,lo)+interpolate(p,p.coefficient,hi))/2;
    }
    const Real en=interpolate(p,p.efn,xe),ep=interpolate(p,p.efp,xs),delta=(ep-en)/kt;
    if(delta==0) return 0;
    const Real logOccupation=std::log(-std::expm1(-std::abs(delta)))-softplus((energy-std::max(ep,en))/kt)-softplus((std::min(ep,en)-energy)/kt);
    const Real vw=-std::expm1(-2*mv*std::max(Real(0),globalEv-energy)*a.valenceInverse_m2/(hb*hb));
    const Real cw=-std::expm1(-2*mc*std::max(Real(0),energy-globalEc)*a.conductionInverse_m2/(hb*hb));
    if(vw==0 || cw==0 || gap==0) return 0;
    const Real prefactor=mass.coupling_J2_m_per_kg*(1+2/std::expm1(c.phononEnergy_J/kt))/(64*std::numbers::pi_v<Real>*std::numbers::pi_v<Real>*c.phononEnergy_J)*
        std::sqrt(mv*mc/(h*length*std::sqrt(2*mr)))*static_cast<Real>(gap)/a.valenceInverse_m2/a.conductionInverse_m2;
    check(positive(prefactor),"unrepresentable cell prefactor");
    const Real value=std::copysign(std::exp(std::log(prefactor)-2*(a.valenceAction+a.conductionAction)+
        std::log(vw)+std::log(cw)+logOccupation+std::log(upper-lower)),delta);
    check(std::isfinite(value),"nonfinite cell scalar integral");return value;
}
}

TunnelingFractionIntegrals tunnelingFractionIntegrals(const std::vector<TunnelingBandKnot>& p,
    Real energy,Real mc,Real mv,Real hb) {
    check(p.size()>=2 && std::isfinite(energy) && positive(mc) && positive(mv) && positive(hb),"invalid fraction kernel input");
    std::vector<Real> z,gap;
    for(const auto& n:p) {
        check(std::isfinite(n.distance_m) && std::isfinite(n.conduction_J) && std::isfinite(n.valence_J) && positive(n.conduction_J-n.valence_J),"invalid fraction band knot");
        gap.push_back(n.conduction_J-n.valence_J);z.push_back((energy-n.valence_J)/gap.back());
    }
    check(z.front()<=0 && z.back()>=1,"fraction profile must bracket both turning points");
    TunnelingFractionIntegrals out;const Real cross=mc/(mc+mv);
    for(Index i=1;i<p.size();++i) {
        const Real span=p[i].distance_m-p[i-1].distance_m;check(positive(span),"unordered fraction path");
        const Real delta=std::abs(z[i]-z[i-1]);
        for(bool valence:{true,false}) {
            const Real mass=valence?mv:mc;
            const Real k=(std::sqrt(2*mass*gap[i-1])+std::sqrt(2*mass*gap[i]))/(2*hb);
            Real action=0,inverse=0;
            if(delta==0) {
                check(z[i]!=0 && z[i]!=1,"singular turning plateau");
                if(z[i]<=0 || z[i]>=1 || (z[i]<cross)!=valence) continue;
                const Real f=valence?z[i]:1-z[i];action=span*k*std::sqrt(f);inverse=span/k/std::sqrt(f);
            } else {
                Real lo=std::max(valence?Real(0):cross,std::min(z[i-1],z[i]));
                Real hi=std::min(valence?cross:Real(1),std::max(z[i-1],z[i]));
                if(lo>=hi) continue;
                if(!valence) {const Real old=lo;lo=1-hi;hi=1-old;}
                const Real a=std::sqrt(lo),b=std::sqrt(hi),scale=span/delta;
                action=k*scale*(2./3)*(hi-lo)*(lo+a*b+hi)/(a+b);
                inverse=scale/k*2*(hi-lo)/(a+b);
            }
            (valence?out.valenceAction:out.conductionAction)+=action;
            (valence?out.valenceInverse_m2:out.conductionInverse_m2)+=inverse;
        }
    }
    check(positive(out.valenceAction) && positive(out.conductionAction) && positive(out.valenceInverse_m2) && positive(out.conductionInverse_m2),"incomplete fraction integrals");return out;
}

std::vector<Real> tunnelingFirstCrossingWeights(const std::vector<Real>& v,Real lower,Real upper) {
    check(!v.empty() && std::isfinite(lower) && std::isfinite(upper) && positive(upper-lower) &&
        std::all_of(v.begin(),v.end(),[](Real x){return std::isfinite(x);}),"invalid first crossing window");
    const Real width=upper-lower;std::vector<Real> weights(v.size(),0);
    weights[0]=std::max(Real(0),upper-std::max(lower,v[0]))/width;
    Real record=v[0];Index lowest=0;
    for(Index i=1;i<v.size();++i) {
        const Real lo=std::max(lower,v[i]),hi=std::min({upper,v[i-1],record}),mid=v[i-1]/2+v[i]/2;
        if(lo<hi) {
            weights[i-1]+=std::max(Real(0),hi-std::max(lo,mid))/width;
            weights[i]+=std::max(Real(0),std::min(hi,mid)-lo)/width;
        }
        if(v[i]<record) {record=v[i];lowest=i;}
    }
    weights[lowest]+=std::max(Real(0),std::min(upper,record)-lower)/width;
    long double sum=0;for(Real w:weights) sum+=w;
    check(std::abs(sum-1)<=tol,"first crossing collection does not cover window");return weights;
}

TunnelingSourceAssembly assembleNonlocalPathCellSources(const TunnelingGeometryMesh& m,
    const TunnelingBandState& s,const NonlocalPathBandToBandConfig& c,
    const TunnelingCellIntegrationPolicy& policy,const TunnelingGeometryOptions& options) {
    s.validate();c.validate();
    check(positive(policy.minimumWidth_kT) && positive(policy.maximumWidth_kT) && policy.maximumWidth_kT>=policy.minimumWidth_kT &&
        positive(policy.variationTolerance) && std::isfinite(policy.variationReference_per_m2_s) && policy.variationReference_per_m2_s>=0 && policy.maximumDepth<=20,"invalid private cell integration policy");
    TunnelingSourceAssembly out;out.constants=c.constants;out.sourceMesh=m;out.cellSeeded=true;out.selectedRegions=c.regions;
    out.boxes=tunnelingSourceBoxes(m,c,TunnelingBoxPolicy::Delaunay);out.resolved=true;
    const auto paths=NonlocalPathGeometry(m).traceFromCells(s,c,options);
    std::vector<Real> ev,ec,efn,efp,coefficient;
    Real globalEv=-std::numeric_limits<Real>::infinity(),globalEc=-globalEv,T=0;
    for(const auto& n:s.nodes) {
        ev.push_back(n.valence_J);ec.push_back(n.conduction_J);efn.push_back(n.electronFermi_J);efp.push_back(n.holeFermi_J);
        if(n.semiconductor && n.material==m.siliconMaterial) {
            if(T==0) T=n.temperature_K;
            check(std::abs(T-n.temperature_K)<=tol*T,"cell source requires isothermal Si");
            check(positive(n.conduction_J-n.valence_J),"nonpositive node gap");
            coefficient.push_back(std::pow(n.conduction_J-n.valence_J,-1.25));
            globalEv=std::max(globalEv,n.valence_J);globalEc=std::min(globalEc,n.conduction_J);
        } else coefficient.push_back(0);
    }
    const Real kt=tunnelingPhysicalConstants(c.constants).boltzmann_J_per_K*T;
    std::vector<NonlocalPairEvent> all;std::uint64_t eventId=0;
    for(const auto& geometry:paths) {
        TunnelingSourceRecord record;record.geometry=geometry;record.seedArea_m2=out.boxes.cellSeedAreas_m2.at(geometry.path.id);
        const auto termination=geometry.path.termination;
        if(termination==TunnelingTermination::ZeroDirection) {record.resolved=true;record.disposition="explicit_geometric_zero";}
        else if(termination==TunnelingTermination::NoBarrier && geometry.path.points.size()==1) {
            record.resolved=true;record.disposition="length_excluded_energy_window";
        }
        else if(termination!=TunnelingTermination::Complete && termination!=TunnelingTermination::ContactBoundary && termination!=TunnelingTermination::LengthLimit && termination!=TunnelingTermination::NoBarrier)
            record.disposition="unresolved_cell_path";
        else try {
            check(geometry.path.points.size()>=2,"incomplete source chord at path length limit");
            const Real lower=geometry.path.energyLower_J,upper=geometry.path.energyUpper_J;
            if(lower>=upper) {record.resolved=true;record.disposition="empty_energy_window";}
            else {
                check(geometry.path.points.size()>=2 && positive(geometry.sourceChord_m),"incomplete source chord");
                Profile p;Real distance=0;
                for(Index i=0;i<geometry.path.points.size();++i) {
                    const auto& point=geometry.path.points[i];if(i) distance+=(point.position_m-geometry.path.points[i-1].position_m).norm();
                    p.ev.push_back(sample(point.support,ev));p.ec.push_back(sample(point.support,ec));
                    p.efn.push_back(sample(point.support,efn));p.efp.push_back(sample(point.support,efp));
                    p.coefficient.push_back(sample(point.support,coefficient));p.bands.push_back({distance,p.ev.back(),p.ec.back(),T});
                }
                struct Sample {Real lower,upper,value;};std::vector<Sample> samples;
                Index evaluations=0;
                auto evaluate=[&](Real a,Real b) {check(++evaluations<=c.maximumIntegrationIntervals,"cell energy quadrature resource limit");return scalar(p,a,b,globalEv,globalEc,c);};
                std::function<void(Real,Real,Real,Index)> visit=[&](Real a,Real b,Real parent,Index depth) {
                    const Real mid=a/2+b/2;check(a<mid && mid<b,"unrepresentable cell energy interval");
                    const Real left=evaluate(a,mid),right=evaluate(mid,b),half=parent/2;
                    const bool split=b-a>policy.maximumWidth_kT*kt || std::max(std::abs(half-left),std::abs(half-right))>
                        policy.variationTolerance*(std::abs(half)+policy.variationReference_per_m2_s);
                    if(split && depth<policy.maximumDepth && b-a>policy.minimumWidth_kT*kt) {visit(a,mid,left,depth+1);visit(mid,b,right,depth+1);return;}
                    if(split) {if(depth>=policy.maximumDepth) ++record.depthPolicyStops;else ++record.minimumPolicyStops;}
                    samples.push_back({a,mid,left});samples.push_back({mid,b,right});
                };
                const Real root=evaluate(lower,upper);
                if(upper-lower<=policy.minimumWidth_kT*kt) samples.push_back({lower,upper,root});else visit(lower,upper,root,0);
                record.energySamples=samples.size();long double signedTotal=0;
                for(const auto& a:samples) {
                    const Real rate=static_cast<Real>(static_cast<long double>(a.value)*record.seedArea_m2/geometry.sourceChord_m);signedTotal+=rate;
                    const auto ew=tunnelingFirstCrossingWeights(p.ec,a.lower,a.upper),hw=tunnelingFirstCrossingWeights(p.ev,a.lower,a.upper);
                    // Cartesian pairing preserves each marginal without an
                    // empirical renormalization, source floor, or lost contact.
                    for(Index h=0;h<hw.size();++h) if(hw[h]>0) for(Index e=0;e<ew.size();++e) if(ew[e]>0)
                        record.collectionEvents.push_back({eventId++,rate*hw[h]*ew[e],geometry.path.points[h].support,geometry.path.points[e].support});
                }
                record.signedEvent_per_m_s=static_cast<Real>(signedTotal);record.resolved=true;record.disposition="cell_evaluated";
            }
        } catch(const std::invalid_argument& e) {record.disposition=e.what();record.collectionEvents.clear();record.resolved=false;}
        out.resolved &= record.resolved;all.insert(all.end(),record.collectionEvents.begin(),record.collectionEvents.end());out.records.push_back(std::move(record));
    }
    out.sources=distributePairEvents(s,all);
    for(Index i=0;i<s.nodes.size();++i) {
        (s.nodes[i].contact?out.contacts.contactElectron_per_m_s:out.contacts.freeElectron_per_m_s)+=out.sources.electronIntegrals_per_m_s[i];
        (s.nodes[i].contact?out.contacts.contactHole_per_m_s:out.contacts.freeHole_per_m_s)+=out.sources.holeIntegrals_per_m_s[i];
    }
    return out;
}
} // namespace vela::experimental
