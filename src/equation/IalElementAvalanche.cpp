#include "vela/equation/IalElementAvalanche.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace vela {
namespace {
using D = detail::Tri3LocalForwardDual;
using Vec = std::array<D,2>;

void nonnegative(Real x) {
    if (!std::isfinite(x) || x < 0.)
        throw std::invalid_argument("Element avalanche requires finite nonnegative parameters");
}

void validate(const ElementIonizationCurve& c) {
    for (Real x : {c.lowA_per_m,c.lowB_V_per_m,c.highA_per_m,c.highB_V_per_m})
        nonnegative(x);
    if (!std::isfinite(c.switchField_V_per_m) || c.switchField_V_per_m <= 0.)
        throw std::invalid_argument("Element avalanche requires a positive switch field");
}

D coefficient(const D& field, const ElementIonizationCurve& c, Real cutoff) {
    if (field.value <= cutoff || field.value == 0.) return D{};
    const auto branch = [&](Real a, Real b) {
        return D(a)*detail::dualExp(-D(b)/field);
    };
    const D low = branch(c.lowA_per_m,c.lowB_V_per_m);
    if (field.value <= .95*c.switchField_V_per_m) return low;
    const D high = branch(c.highA_per_m,c.highB_V_per_m);
    if (field.value >= 1.05*c.switchField_V_per_m) return high;
    const D z = (field/D(c.switchField_V_per_m)-D(.95))/D(.1);
    D weight((1.-std::cos(std::numbers::pi*z.value))/2.);
    for (std::size_t k=0; k<9; ++k)
        weight.derivative[k]=.5*std::numbers::pi*std::sin(std::numbers::pi*z.value)*z.derivative[k];
    return (D(1.)-weight)*low+weight*high;
}

D extremaDensity(const std::array<D,3>& density) {
    const auto less=[](const D& a,const D& b){return a.value<b.value;};
    const D& low=*std::min_element(density.begin(),density.end(),less);
    const D& high=*std::max_element(density.begin(),density.end(),less);
    if (high.value==0.) return D{};
    // Avoid the product of two large densities in the primal expression.
    return D(2.)*low/(D(1.)+low/high);
}

D norm(const Vec& v) {return detail::dualSqrt(v[0]*v[0]+v[1]*v[1]);}
void check(const D& x) {
    if (!std::isfinite(x.value)) throw std::runtime_error("Nonfinite element avalanche value");
    for (Real d:x.derivative)
        if (!std::isfinite(d)) throw std::runtime_error("Nonfinite element avalanche derivative");
}
} // namespace

IalElementAvalancheResult evaluateIalElementAvalanche(
    const IalElementGeometry& geometry,
    const std::array<IalElementVertexState,3>& state,
    const std::array<const IalMobility*,3>& electronModels,
    const std::array<const IalMobility*,3>& holeModels,
    const IalElementAvalancheOptions& options)
{
    if (options.mobilityVertex>=3)
        throw std::invalid_argument("Element avalanche mobility vertex must be 0, 1 or 2");
    validate(options.electron); validate(options.hole);
    nonnegative(options.minimumField_V_per_m);
    Real total=0.;
    for (Real v:options.sourceMeasure_m2) {nonnegative(v);total+=v;}
    if (!std::isfinite(total) || total<=0.)
        throw std::invalid_argument("Element avalanche requires positive total source measure");
    for (const auto& s:state)
        if (s.temperature_K!=300.)
            throw std::invalid_argument("Element avalanche identification is limited to 300 K");

    auto selectedGeometry=geometry;
    selectedGeometry.vertexMeasure_m2={0.,0.,0.};
    selectedGeometry.vertexMeasure_m2[options.mobilityVertex]=1.;
    IalElementMobilityOptions mobilityOptions;
    mobilityOptions.referenceDensity_m3=0.;
    mobilityOptions.highField=options.highFieldMobility;
    const auto mobility=evaluateIalElementMobility(selectedGeometry,state,
        electronModels,holeModels,mobilityOptions);

    const auto& xy=geometry.coordinates_m;
    const Real ax=xy[1][0]-xy[0][0],ay=xy[1][1]-xy[0][1];
    const Real bx=xy[2][0]-xy[0][0],by=xy[2][1]-xy[0][1];
    const Real det=ax*by-ay*bx; // Already validated by the mobility kernel.
    const auto gradient=[&](const std::array<D,3>& v)->Vec {
        const D a=v[1]-v[0],b=v[2]-v[0];
        return {(a*D(by)-b*D(ay))/D(det),(b*D(ax)-a*D(bx))/D(det)};
    };
    std::array<D,3> psi,fn,fp,n,p;
    for (std::size_t i=0;i<3;++i) {
        const auto& s=state[i];
        psi[i]=D::variable(s.potential_V,3*i);
        fn[i]=D::variable(s.electronQf_V,3*i+1);
        fp[i]=D::variable(s.holeQf_V,3*i+2);
        n[i]=D(s.electrons_m3);p[i]=D(s.holes_m3);
        n[i].derivative[3*i]=s.electronResponse_m3_per_V;
        n[i].derivative[3*i+1]=-s.electronResponse_m3_per_V;
        p[i].derivative[3*i]=-s.holeResponse_m3_per_V;
        p[i].derivative[3*i+2]=s.holeResponse_m3_per_V;
    }
    const D electric=norm(gradient(psi));
    const auto carrier=[&](const std::array<D,3>& density,const std::array<D,3>& qf,
                           const D& mu,const ElementIonizationCurve& curve) {
        IalElementAvalancheCarrier r;
        const Vec g=gradient(qf);
        r.currentField_V_per_m=geometry.partialBoundaryLayer
            ?detail::dualAbs(g[0]*D(geometry.boundaryTangent[0])+g[1]*D(geometry.boundaryTangent[1]))
            :norm(g);
        r.ionizationField_V_per_m=geometry.touchesEffectiveElectrode?electric:r.currentField_V_per_m;
        r.density_m3=extremaDensity(density);
        r.mobility_m2_per_Vs=mu;
        r.alpha_per_m=coefficient(r.ionizationField_V_per_m,curve,options.minimumField_V_per_m);
        r.generation_m3_per_s=r.alpha_per_m*r.density_m3*mu*r.currentField_V_per_m;
        for (std::size_t i=0;i<3;++i) {
            r.source_per_m_s[i]=r.generation_m3_per_s*D(options.sourceMeasure_m2[i]);
            check(r.source_per_m_s[i]);
        }
        for (const auto& x:{r.currentField_V_per_m,r.ionizationField_V_per_m,
             r.density_m3,r.mobility_m2_per_Vs,r.alpha_per_m,r.generation_m3_per_s}) check(x);
        return r;
    };
    return {carrier(n,fn,mobility.electron,options.electron),
            carrier(p,fp,mobility.hole,options.hole)};
}
} // namespace vela
