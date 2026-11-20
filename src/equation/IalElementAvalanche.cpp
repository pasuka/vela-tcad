#include "vela/equation/IalElementAvalanche.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace vela {
namespace {
using D = detail::Tri3LocalForwardDual;


void nonnegative(long double x) {
    if (!std::isfinite(x) || x < 0.)
        throw std::invalid_argument("Element avalanche requires finite nonnegative parameters");
}

void validate(const auto& c) {
    for (auto x : {c.lowA_per_m,c.lowB_V_per_m,c.highA_per_m,c.highB_V_per_m})
        nonnegative(x);
    if (!std::isfinite(c.switchField_V_per_m) || c.switchField_V_per_m <= 0.)
        throw std::invalid_argument("Element avalanche requires a positive switch field");
}

template<class Scalar, class Curve>
detail::BasicTri3LocalForwardDual<Scalar> coefficient(
    const detail::BasicTri3LocalForwardDual<Scalar>& field, const Curve& c, Scalar cutoff) {
    using D = detail::BasicTri3LocalForwardDual<Scalar>;
    if (field.value <= cutoff || field.value == 0.) return D{};
    const auto branch = [&](Scalar a, Scalar b) {
        return D(a)*detail::dualExp(-D(b)/field);
    };
    const D low = branch(c.lowA_per_m,c.lowB_V_per_m);
    const auto lower = static_cast<Scalar>(0.95L);
    const auto upper = static_cast<Scalar>(1.05L);
    if (field.value <= lower*c.switchField_V_per_m) return low;
    const D high = branch(c.highA_per_m,c.highB_V_per_m);
    if (field.value >= upper*c.switchField_V_per_m) return high;
    const D z = (field/D(c.switchField_V_per_m)-D(lower))/D(static_cast<Scalar>(0.1L));
    constexpr auto pi = std::numbers::pi_v<Scalar>;
    D weight((1.-std::cos(pi*z.value))/2.);
    for (std::size_t k=0; k<9; ++k)
        weight.derivative[k]=.5*pi*std::sin(pi*z.value)*z.derivative[k];
    return (D(1.)-weight)*low+weight*high;
}

template<class D>
D extremaDensity(const std::array<D,3>& density) {
    const auto less=[](const D& a,const D& b){return a.value<b.value;};
    const D& low=*std::min_element(density.begin(),density.end(),less);
    const D& high=*std::max_element(density.begin(),density.end(),less);
    if (high.value==0.) return D{};
    // Avoid the product of two large densities in the primal expression.
    return D(2.)*low/(D(1.)+low/high);
}

template<class D>
D norm(const std::array<D,2>& v) {return detail::dualSqrt(v[0]*v[0]+v[1]*v[1]);}

template<class D>
void check(const D& x) {
    if (!std::isfinite(x.value)) throw std::runtime_error("Nonfinite element avalanche value");
    for (const auto d:x.derivative)
        if (!std::isfinite(d)) throw std::runtime_error("Nonfinite element avalanche derivative");
}
template<class Scalar>
void validateSource(const std::array<Scalar, 3>& shares) {
    Scalar total = 0;
    for (auto value : shares) {
        if (!std::isfinite(value) || value < 0)
            throw std::invalid_argument("Element avalanche requires finite nonnegative source shares");
        total += value;
    }
    if (!std::isfinite(total) || total<=0.)
        throw std::invalid_argument("Element avalanche requires positive total source measure");
}

template<class Scalar, class Curve>
BasicIalElementAvalancheResult<Scalar> evaluatePrepared(
    const BasicIalElementGeometry<Scalar>& geometry,
    const std::array<detail::BasicTri3LocalForwardDual<Scalar>, 3>& psi,
    const std::array<detail::BasicTri3LocalForwardDual<Scalar>, 3>& fn,
    const std::array<detail::BasicTri3LocalForwardDual<Scalar>, 3>& fp,
    const std::array<detail::BasicTri3LocalForwardDual<Scalar>, 3>& n,
    const std::array<detail::BasicTri3LocalForwardDual<Scalar>, 3>& p,
    const std::array<detail::BasicTri3LocalForwardDual<Scalar>, 2>& mobilities,
    GradQfContactField contactField, const std::array<Scalar, 3>& shares,
    const Curve& electron, const Curve& hole, Scalar minimumField) {
    using D = detail::BasicTri3LocalForwardDual<Scalar>;
    using Vec = std::array<D, 2>;
    validate(electron);
    validate(hole);
    nonnegative(minimumField);
    validateSource(shares);
    const detail::IalTri3Gradient gradient(geometry.coordinates_m);
    const D electric=norm(gradient(psi));
    const auto carrier=[&](const std::array<D,3>& density,const std::array<D,3>& qf,
                           const D& mu,const Curve& curve) {
        BasicIalElementAvalancheCarrier<Scalar> r;
        const Vec g=gradient(qf);
        // At an electrode, native PMI retains the full QF current even when
        // the element also belongs to an interface layer. Only alpha falls
        // back to the electric field; it must not replace the current drive.
        const bool tangentialCurrent=geometry.partialBoundaryLayer
            && !geometry.touchesEffectiveElectrode;
        r.currentField_V_per_m=tangentialCurrent
            ?detail::dualAbs(g[0]*D(geometry.boundaryTangent[0])+g[1]*D(geometry.boundaryTangent[1]))
            :norm(g);
        const bool electricAtContact=geometry.touchesEffectiveElectrode &&
            contactField==GradQfContactField::ElectricField;
        r.ionizationField_V_per_m=electricAtContact?electric:r.currentField_V_per_m;
        r.density_m3=extremaDensity(density);
        r.mobility_m2_per_Vs=mu;
        r.alpha_per_m=coefficient(r.ionizationField_V_per_m,curve,minimumField);
        r.generation_m3_per_s=r.alpha_per_m*r.density_m3*mu*r.currentField_V_per_m;
        for (std::size_t i=0;i<3;++i) {
            r.source_per_m_s[i]=r.generation_m3_per_s*D(shares[i]);
            check(r.source_per_m_s[i]);
        }
        for (const auto& x:{r.currentField_V_per_m,r.ionizationField_V_per_m,
             r.density_m3,r.mobility_m2_per_Vs,r.alpha_per_m,r.generation_m3_per_s}) check(x);
        return r;
    };
    return {carrier(n,fn,mobilities[0],electron),
            carrier(p,fp,mobilities[1],hole)};
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
    for (const auto& s:state)
        if (s.temperature_K!=300.)
            throw std::invalid_argument("Element avalanche identification is limited to 300 K");

    auto selectedGeometry=geometry;
    selectedGeometry.vertexMeasure_m2={0.,0.,0.};
    selectedGeometry.vertexMeasure_m2[options.mobilityVertex]=1.;
    IalElementMobilityOptions mobilityOptions;
    mobilityOptions.referenceDensity_m3=0.;
    mobilityOptions.highField=options.highFieldMobility;
    mobilityOptions.contactField=options.mobilityContactField;
    mobilityOptions.populationThermalVoltage_V=options.populationThermalVoltage_V;
    const auto mobility=evaluateIalElementMobility(selectedGeometry,state,
        electronModels,holeModels,mobilityOptions);

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
        if (options.populationThermalVoltage_V>0.)
            detail::seedBoltzmannPopulationVertex(
                i,options.populationThermalVoltage_V,psi[i],fn[i],fp[i],n[i],p[i]);
    }
    return evaluatePrepared(geometry, psi, fn, fp, n, p,
                            std::array{mobility.electron, mobility.hole}, options.contactField, options.sourceMeasure_m2,
                            options.electron, options.hole, options.minimumField_V_per_m);
}

BasicIalElementAvalancheResult<ArithmeticReal> evaluateIalPopulationAvalancheFromMobility(
    const BasicIalElementGeometry<ArithmeticReal>& geometry,
    const std::array<IalAvalanchePopulationVertex, 3>& state,
    ArithmeticReal potentialScale_V, ArithmeticReal thermalVoltage_V,
    const std::array<detail::BasicTri3LocalForwardDual<ArithmeticReal>, 2>& mobility,
    const IalElementAvalancheOptions& options,
    const std::array<long double, 3>* arithmeticSourceMeasures,
    const ArithmeticIalIonizationCoefficients* arithmeticCoefficients)
{
    std::array<ArithmeticReal, 3> shares;
    for (std::size_t i = 0; i < shares.size(); ++i)
        shares[i] = arithmeticSourceMeasures ? (*arithmeticSourceMeasures)[i] : options.sourceMeasure_m2[i];
    for (const auto& value : mobility) {
        check(value);
        if (!(value.value > 0))
            throw std::invalid_argument("Population avalanche requires positive mobility");
    }
    const auto fields = detail::buildIalPopulationFields(state, potentialScale_V, thermalVoltage_V);
    if (arithmeticCoefficients)
        return evaluatePrepared(geometry, fields.psi, fields.electronQf, fields.holeQf,
                                fields.electrons, fields.holes, mobility, options.contactField, shares,
                                arithmeticCoefficients->electron, arithmeticCoefficients->hole,
                                arithmeticCoefficients->minimumField_V_per_m);
    return evaluatePrepared(geometry, fields.psi, fields.electronQf, fields.holeQf,
                            fields.electrons, fields.holes, mobility, options.contactField, shares,
                            options.electron, options.hole, ArithmeticReal(options.minimumField_V_per_m));
}

} // namespace vela
