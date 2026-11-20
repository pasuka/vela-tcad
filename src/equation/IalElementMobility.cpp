#include "vela/equation/IalElementMobility.h"
#include "vela/physics/IalHighFieldMobility.h"
#include "vela/core/IalKernelProfiling.h"
#include "vela/discretization/ScharfetterGummel.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace vela {
namespace {
using D=detail::Tri3LocalForwardDual;
struct PassTimer {
    using Clock=std::chrono::steady_clock;
    int phase;bool enabled;Clock::time_point start;
    explicit PassTimer(int p):phase(p),enabled(ialKernelProfile.timingEnabled) {
        ++ialKernelProfile.passes[phase];if(enabled)start=Clock::now();
    }
    ~PassTimer(){if(enabled)ialKernelProfile.seconds[phase]+=std::chrono::duration<double>(Clock::now()-start).count();}
};
template<class Dual> Dual norm(const std::array<Dual,2>& v) {
    return detail::dualSqrt(v[0]*v[0]+v[1]*v[1]);
}
template<class Dual, class Scalar>
Dual dot(const std::array<Dual,2>& v, const std::array<Scalar,2>& g) {
    return v[0]*Dual(g[0])+v[1]*Dual(g[1]);
}
template<class Scalar, class Dual>
std::array<Dual,4> drives(const BasicIalElementGeometry<Scalar>& geometry,
    const std::array<Dual,3>& psi, const std::array<Dual,3>& fn,
    const std::array<Dual,3>& fp, GradQfContactField contactField) {
    const detail::IalTri3Gradient gradient(geometry.coordinates_m);
    const std::array<Dual,3> distance{Dual(geometry.interfaceDistance_m[0]),
        Dual(geometry.interfaceDistance_m[1]), Dual(geometry.interfaceDistance_m[2])};
    const auto e=gradient(psi), gn=gradient(fn), gp=gradient(fp), gd=gradient(distance);
    const std::array<Scalar,2> normal{gd[0].value,gd[1].value};
    const auto normalComponent=dot(e,normal);
    const auto parallel=norm(std::array{e[0]-normalComponent*Dual(normal[0]),
                                       e[1]-normalComponent*Dual(normal[1])});
    const bool tangential=geometry.partialBoundaryLayer && !geometry.touchesEffectiveElectrode;
    auto driveN=tangential?detail::dualAbs(dot(gn,geometry.boundaryTangent)):norm(gn);
    auto driveP=tangential?detail::dualAbs(dot(gp,geometry.boundaryTangent)):norm(gp);
    if (geometry.touchesEffectiveElectrode && contactField==GradQfContactField::ElectricField)
        driveN=driveP=norm(e);
    return {detail::dualAbs(normalComponent),parallel,driveN,driveP};
}
template<class Scalar> void nonnegative(Scalar v) {
    if (!std::isfinite(v) || v<0.) throw std::invalid_argument("IALMob element expects finite nonnegative SI inputs");
}
template<class Scalar> void finite(Scalar v) {
    if (!std::isfinite(v)) throw std::invalid_argument("IALMob element expects finite potentials and geometry");
}
D limited(const D& mu, const D& drive, const FieldMobilityParameters& p) {
    if (drive.value==0.) return mu;
    return mu/detail::dualPow(D(1.)+detail::dualPow(mu*drive/D(p.saturationVelocity),p.beta),1./p.beta);
}
template<class Scalar>
detail::BasicTri3LocalForwardDual<Scalar> limited(
    const detail::BasicTri3LocalForwardDual<Scalar>& mu,
    const detail::BasicTri3LocalForwardDual<Scalar>& drive,
    const IalHighFieldDecimal& parameters)
{
    if (drive.value == 0.)
        return mu;
    const auto field = static_cast<ArithmeticReal>(drive.value);
    const auto mobility = static_cast<ArithmeticReal>(mu.value);
    const auto power = std::pow(mobility * field / parameters.velocity_m_s, parameters.beta);
    const auto denominator = VELA_ARITHMETIC_LITERAL(1.0) + power;
    const auto scale = std::pow(denominator, -VELA_ARITHMETIC_LITERAL(1.0) / parameters.beta);
    const auto value = mobility * scale;
    const auto dMobility = scale / denominator;
    const auto dField = -value * (power / denominator) / field;
    detail::BasicTri3LocalForwardDual<Scalar> result(static_cast<Scalar>(value));
    for (std::size_t i = 0; i < result.derivative.size(); ++i)
        result.derivative[i] = static_cast<Scalar>(dMobility * mu.derivative[i] + dField * drive.derivative[i]);
    return result;
}
}

static IalElementMobilityResult evaluateImpl(
    const IalElementGeometry& geometry,
    const std::array<IalElementVertexState,3>& state,
    const std::array<const IalMobility*,3>& electronModels,
    const std::array<const IalMobility*,3>& holeModels,
    const IalElementMobilityOptions& options, bool thermalDirections,
    std::array<IalMobilityDifferential,6>& localDifferentials,
    std::array<IalHighFieldResult,6>& highFieldDifferentials)
{
    PassTimer timer(thermalDirections?2:options.spatialDerivatives?1:0);
    nonnegative(options.referenceDensity_m3);
    if(options.temperatureDerivatives && !options.spatialDerivatives)
        throw std::invalid_argument("IALMob values-only evaluation cannot request temperature columns");
    if (options.decimalHighField && options.temperatureDependentHighField)
        throw std::invalid_argument("IALMob decimal high-field coefficients require the isothermal SLOT evaluation");
    for (const auto& p:{options.electronField,options.holeField})
        if (!std::isfinite(p.saturationVelocity)||p.saturationVelocity<=0.||
            !std::isfinite(p.beta)||p.beta<=0.)
            throw std::invalid_argument("IALMob HFS requires positive velocity and beta");
    for (Real v:{options.electronVelocityTemperatureExponent,options.holeVelocityTemperatureExponent,
        options.electronBetaTemperatureExponent,options.holeBetaTemperatureExponent}) finite(v);
    Real measure=0.;
    for (int i=0;i<3;++i) {
        for (Real x:geometry.coordinates_m[i]) finite(x);
        nonnegative(geometry.interfaceDistance_m[i]);
        nonnegative(geometry.vertexMeasure_m2[i]);
        measure+=geometry.vertexMeasure_m2[i];
        if (!electronModels[i]||!holeModels[i]) throw std::invalid_argument("IALMob missing orientation model");
        for (Real x:{state[i].potential_V,state[i].electronQf_V,state[i].holeQf_V}) finite(x);
        for (Real x:{state[i].donors_m3,state[i].acceptors_m3,state[i].electrons_m3,
            state[i].holes_m3,state[i].electronResponse_m3_per_V,state[i].holeResponse_m3_per_V}) nonnegative(x);
    }
    if (!std::isfinite(measure)||measure<=0.) throw std::invalid_argument("IALMob requires positive total vertex measure");
    if (geometry.partialBoundaryLayer) {
        for (Real x:geometry.boundaryTangent) finite(x);
        if (std::abs(std::hypot(geometry.boundaryTangent[0],geometry.boundaryTangent[1])-1.)>1e-12)
            throw std::invalid_argument("IALMob boundary tangent must be unit length");
    }
    for (const auto& v:state) {
        finite(v.temperature_K);finite(v.electronTemperatureResponse_m3_per_K);finite(v.holeTemperatureResponse_m3_per_K);
        if(v.temperature_K<50.)throw std::invalid_argument("IALMob element requires temperature >=50 K");
        if(options.highField && !options.temperatureDependentHighField && v.temperature_K!=300.)
            throw std::invalid_argument("Hot IALMob element requires explicit temperature-dependent HFS");
    }
    std::array<D,3> psi,qfn,qfp,n,p,temperature;
    for (int i=0;i<3;++i) {
        psi[i]=D::variable(state[i].potential_V,3*i);
        qfn[i]=D::variable(state[i].electronQf_V,3*i+1);
        qfp[i]=D::variable(state[i].holeQf_V,3*i+2);
        n[i]=D(state[i].electrons_m3);
        n[i].derivative[3*i]=state[i].electronResponse_m3_per_V;
        n[i].derivative[3*i+1]=-state[i].electronResponse_m3_per_V;
        p[i]=D(state[i].holes_m3);
        p[i].derivative[3*i]=-state[i].holeResponse_m3_per_V;
        p[i].derivative[3*i+2]=state[i].holeResponse_m3_per_V;
        if (options.populationThermalVoltage_V>0.)
            detail::seedBoltzmannPopulationVertex(
                i,options.populationThermalVoltage_V,psi[i],qfn[i],qfp[i],n[i],p[i]);
        temperature[i]=D(state[i].temperature_K);
        if(!options.spatialDerivatives){
            psi[i]=D(state[i].potential_V);qfn[i]=D(state[i].electronQf_V);qfp[i]=D(state[i].holeQf_V);
            n[i]=D(state[i].electrons_m3);p[i]=D(state[i].holes_m3);
        }
        if(thermalDirections){
            psi[i]=D(state[i].potential_V);qfn[i]=D(state[i].electronQf_V);qfp[i]=D(state[i].holeQf_V);
            n[i]=D(state[i].electrons_m3);p[i]=D(state[i].holes_m3);
            n[i].derivative[i]=state[i].electronTemperatureResponse_m3_per_K;
            p[i].derivative[i]=state[i].holeTemperatureResponse_m3_per_K;
            temperature[i].derivative[i]=1.;
        }
    }
    const auto [en,parallel,driveN,driveP] = drives(geometry,psi,qfn,qfp,options.contactField);
    IalElementMobilityResult result;
    for (int i=0;i<3;++i) {
        const IalMobilityState local{state[i].donors_m3,state[i].acceptors_m3,n[i].value,p[i].value,
            en.value,geometry.interfaceDistance_m[i],state[i].temperature_K};
        const auto low=[&](const IalMobility& model,int carrier) {
            if(!options.spatialDerivatives)return D(model.evaluate(local,options.screeningCache,
                options.preparationCache,options.generatedLowFieldPartials).mobility_m2_per_Vs);
            auto& saved=localDifferentials[2*i+carrier];
            if(!thermalDirections || !options.reuseThermalLocalDifferentials)
                saved=model.evaluateWithDerivatives(local,options.screeningCache,options.preparationCache,options.generatedLowFieldPartials);
            const auto& r=saved;
            D mu(r.result.mobility_m2_per_Vs);
            for (int k=0;k<9;++k)
                mu.derivative[k]=r.derivative_SI[2]*n[i].derivative[k]+
                    r.derivative_SI[3]*p[i].derivative[k]+r.derivative_SI[4]*en.derivative[k]+
                    r.temperatureDerivative_m2_per_Vs_K*temperature[i].derivative[k];
            return mu;
        };
        const D lowN=low(*electronModels[i],0), lowP=low(*holeModels[i],1);
        D finalN=lowN,finalP=lowP;
        if (options.highField) {
            const D fractionN=options.referenceDensity_m3==0.?D(1.):n[i]/(n[i]+D(options.referenceDensity_m3));
            const D fractionP=options.referenceDensity_m3==0.?D(1.):p[i]/(p[i]+D(options.referenceDensity_m3));
            const D fieldN=fractionN*driveN+(D(1.)-fractionN)*parallel;
            const D fieldP=fractionP*driveP+(D(1.)-fractionP)*parallel;
            if(options.temperatureDependentHighField){
                const auto hot=[&](const D& mu,const D& field,const FieldMobilityParameters& base,Real ve,Real be,int carrier){
                    const IalHighFieldParameters parameters{base.saturationVelocity,base.beta,ve,be};
                    if(options.explicitHighFieldPartials && !options.spatialDerivatives)
                        return D(evaluateIalHighFieldMobilityValue(mu.value,field.value,state[i].temperature_K,parameters));
                    auto& h=highFieldDifferentials[2*i+carrier];
                    if(!thermalDirections || !options.reuseThermalHighField)
                        h=(options.explicitHighFieldPartials?evaluateIalHighFieldMobilityExplicit:evaluateIalHighFieldMobility)
                            (mu.value,field.value,state[i].temperature_K,parameters);
                    else ++ialKernelProfile.highFieldReuses;
                    D result(h.mobility_m2_per_Vs);
                    for(int k=0;k<9;++k)result.derivative[k]=h.lowFieldDerivative*mu.derivative[k]+
                        h.drivingFieldDerivative_m3_per_V2s*field.derivative[k]+h.temperatureDerivative_m2_per_Vs_K*temperature[i].derivative[k];
                    return result;
                };
                finalN=hot(lowN,fieldN,options.electronField,options.electronVelocityTemperatureExponent,options.electronBetaTemperatureExponent,0);
                finalP=hot(lowP,fieldP,options.holeField,options.holeVelocityTemperatureExponent,options.holeBetaTemperatureExponent,1);
            } else if (options.decimalHighField) {
                finalN = limited(lowN, fieldN, (*options.decimalHighField)[0]);
                finalP = limited(lowP, fieldP, (*options.decimalHighField)[1]);
            }else{
                finalN=limited(lowN,fieldN,options.electronField);
                finalP=limited(lowP,fieldP,options.holeField);
            }
        }
        const D w(geometry.vertexMeasure_m2[i]/measure);
        result.electronLowField=result.electronLowField+w*lowN;
        result.holeLowField=result.holeLowField+w*lowP;
        result.electron=result.electron+w*finalN;
        result.hole=result.hole+w*finalP;
    }
    for (const auto& r:{result.electron,result.hole,result.electronLowField,result.holeLowField}) {
        if (!std::isfinite(r.value)||r.value<=0.) throw std::runtime_error("Invalid IALMob element mobility");
        for (Real d:r.derivative)
            if (!std::isfinite(d)) throw std::runtime_error("Nonfinite IALMob element derivative");
    }
    return result;
}
IalElementMobilityResult evaluateIalElementMobility(
    const IalElementGeometry& geometry,const std::array<IalElementVertexState,3>& state,
    const std::array<const IalMobility*,3>& electrons,const std::array<const IalMobility*,3>& holes,
    const IalElementMobilityOptions& options)
{
    std::array<IalMobilityDifferential,6> localDifferentials;
    std::array<IalHighFieldResult,6> highFieldDifferentials;
    auto result=evaluateImpl(geometry,state,electrons,holes,options,false,localDifferentials,highFieldDifferentials);
    if(options.temperatureDerivatives){
        const auto t=evaluateImpl(geometry,state,electrons,holes,options,true,localDifferentials,highFieldDifferentials);
        for(int i=0;i<3;++i){
            result.electronTemperatureDerivative[i]=t.electron.derivative[i];
            result.holeTemperatureDerivative[i]=t.hole.derivative[i];
            result.electronLowTemperatureDerivative[i]=t.electronLowField.derivative[i];
            result.holeLowTemperatureDerivative[i]=t.holeLowField.derivative[i];
        }
    }
    return result;
}

detail::IalPopulationFields detail::buildIalPopulationFields(
    const std::array<IalPopulationVertex,3>& state,
    ArithmeticReal potentialScale_V, ArithmeticReal thermalVoltage_V)
{
    using Dual = IalPopulationDual;
    if (!(potentialScale_V > 0) || !std::isfinite(potentialScale_V) ||
        !(thermalVoltage_V > 0) || !std::isfinite(thermalVoltage_V))
        throw std::invalid_argument("IALMob population basis requires positive potential and thermal scales");
    for (const auto& vertex : state) {
        finite(vertex.potential);
        for (const auto value : {vertex.electrons_m3, vertex.holes_m3,
                                 vertex.electronIntrinsic_m3, vertex.holeIntrinsic_m3})
            if (!(value > 0) || !std::isfinite(value))
                throw std::invalid_argument("IALMob population basis requires finite positive populations");
    }
    IalPopulationFields out;
    const auto& ref = state.front();
    for (std::size_t i = 0; i < state.size(); ++i) {
        const auto& vertex = state[i];
        // An element gradient only depends on differences. Subtract the
        // reference in the solver coordinate before converting to volts.
        const auto delta = (vertex.potential - ref.potential) * potentialScale_V;
        out.psi[i] = Dual::variable(delta, 3*i);
        out.electronQf[i] = Dual(delta - thermalVoltage_V *
            (arithmeticPopulationLogRatio(vertex.electrons_m3, ref.electrons_m3) -
             arithmeticPopulationLogRatio(vertex.electronIntrinsic_m3, ref.electronIntrinsic_m3)));
        out.holeQf[i] = Dual(delta + thermalVoltage_V *
            (arithmeticPopulationLogRatio(vertex.holes_m3, ref.holes_m3) -
             arithmeticPopulationLogRatio(vertex.holeIntrinsic_m3, ref.holeIntrinsic_m3)));
        out.electrons[i] = Dual(vertex.electrons_m3);
        out.holes[i] = Dual(vertex.holes_m3);
        seedBoltzmannPopulationVertex(i, thermalVoltage_V, out.psi[i],
                                      out.electronQf[i], out.holeQf[i], out.electrons[i], out.holes[i]);
        for (auto* value : {&out.psi[i], &out.electronQf[i], &out.holeQf[i]})
            value->derivative[3*i] *= potentialScale_V;
    }
    return out;
}

IalPopulationMobility evaluateIalPopulationMobility(
    const BasicIalElementGeometry<ArithmeticReal>& geometry,
    const detail::IalPopulationFields& fields,
    const std::array<std::array<ArithmeticReal,2>,3>& doping_m3,
    const std::array<const BasicIalMobilityParameters<long double>*,3>& electrons,
    const std::array<const BasicIalMobilityParameters<long double>*,3>& holes,
    ArithmeticReal temperature_K, ArithmeticReal referenceDensity_m3,
    const std::optional<std::array<IalHighFieldDecimal,2>>& highField,
    GradQfContactField contactField, IalScreeningMethod method)
{
    using Dual = detail::IalPopulationDual;
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80)
        throw std::logic_error("SLOT population mobility requires the extended80 runner");
    nonnegative(referenceDensity_m3);
    if (highField && temperature_K != 300.L)
        throw std::invalid_argument("SLOT population HFS requires 300 K");
    if (highField)
        for (const auto& value : *highField)
            if (!std::isfinite(value.velocity_m_s) || value.velocity_m_s <= 0 ||
                !std::isfinite(value.beta) || value.beta <= 0)
                throw std::invalid_argument("SLOT population HFS requires positive velocity and beta");
    const auto [normal, parallel, driveN, driveP] = drives(
        geometry, fields.psi, fields.electronQf, fields.holeQf, contactField);
    IalPopulationMobility result;
    for (std::size_t i = 0; i < doping_m3.size(); ++i) {
        if (!electrons[i] || !holes[i])
            throw std::invalid_argument("SLOT population mobility lacks a bound orientation model");
        const auto& n = fields.electrons[i];
        const auto& p = fields.holes[i];
        const BasicIalMobilityState<long double> local{
            doping_m3[i][0], doping_m3[i][1], n.value, p.value,
            normal.value, geometry.interfaceDistance_m[i], temperature_K};
        const auto low = [&](const auto& parameters, bool electron) {
            const auto value = evaluateExtendedSlotIalMobility(local, parameters, electron, method);
            Dual out(static_cast<ArithmeticReal>(value.result.mobility_m2_per_Vs));
            for (std::size_t k = 0; k < out.derivative.size(); ++k)
                out.derivative[k] = static_cast<ArithmeticReal>(
                    value.derivative_SI[2]*n.derivative[k] + value.derivative_SI[3]*p.derivative[k]
                    + value.derivative_SI[4]*normal.derivative[k]);
            return out;
        };
        result.electronLowField[i] = low(*electrons[i], true);
        result.holeLowField[i] = low(*holes[i], false);
        const auto applyHighField = [&](const Dual& mu, const Dual& density, const Dual& drive, std::size_t carrier) {
            if (!highField) return mu;
            const auto fraction = referenceDensity_m3 == 0
                ? Dual(1) : density/(density+Dual(referenceDensity_m3));
            return limited(mu, fraction*drive + (Dual(1)-fraction)*parallel, (*highField)[carrier]);
        };
        result.electron[i] = applyHighField(result.electronLowField[i], n, driveN, 0);
        result.hole[i] = applyHighField(result.holeLowField[i], p, driveP, 1);
        for (const auto& value : {result.electron[i], result.hole[i]}) {
            if (!std::isfinite(value.value) || value.value <= 0 ||
                !std::ranges::all_of(value.derivative, [](auto derivative) { return std::isfinite(derivative); }))
                throw std::runtime_error("Invalid SLOT population mobility or derivative");
        }
    }
    return result;
}
} // namespace vela
