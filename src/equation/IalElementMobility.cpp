#include "vela/equation/IalElementMobility.h"
#include "vela/physics/IalHighFieldMobility.h"
#include "vela/core/IalKernelProfiling.h"
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace vela {
namespace {
using D=detail::Tri3LocalForwardDual;
using Vec=std::array<D,2>;
struct PassTimer {
    using Clock=std::chrono::steady_clock;
    int phase;bool enabled;Clock::time_point start;
    explicit PassTimer(int p):phase(p),enabled(ialKernelProfile.timingEnabled) {
        ++ialKernelProfile.passes[phase];if(enabled)start=Clock::now();
    }
    ~PassTimer(){if(enabled)ialKernelProfile.seconds[phase]+=std::chrono::duration<double>(Clock::now()-start).count();}
};
D norm(const Vec& v) { return detail::dualSqrt(v[0]*v[0]+v[1]*v[1]); }
D dot(const Vec& v, const std::array<Real,2>& g) { return v[0]*D(g[0])+v[1]*D(g[1]); }
void nonnegative(Real v) {
    if (!std::isfinite(v) || v<0.) throw std::invalid_argument("IALMob element expects finite nonnegative SI inputs");
}
void finite(Real v) {
    if (!std::isfinite(v)) throw std::invalid_argument("IALMob element expects finite potentials and geometry");
}
D limited(const D& mu, const D& drive, const FieldMobilityParameters& p) {
    if (drive.value==0.) return mu;
    return mu/detail::dualPow(D(1.)+detail::dualPow(mu*drive/D(p.saturationVelocity),p.beta),1./p.beta);
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
    const auto& xy=geometry.coordinates_m;
    const Real ax=xy[1][0]-xy[0][0], ay=xy[1][1]-xy[0][1];
    const Real bx=xy[2][0]-xy[0][0], by=xy[2][1]-xy[0][1];
    const Real det=ax*by-ay*bx;
    if (!std::isfinite(det)||det==0.) throw std::invalid_argument("IALMob requires a nondegenerate triangle");
    const auto gradient=[&](const std::array<D,3>& v)->Vec {
        const D a=v[1]-v[0], b=v[2]-v[0];
        return {(a*D(by)-b*D(ay))/D(det),(b*D(ax)-a*D(bx))/D(det)};
    };
    std::array<D,3> psi,qfn,qfp,n,p,distance,temperature;
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
        distance[i]=D(geometry.interfaceDistance_m[i]);
    }
    const Vec e=gradient(psi), gn=gradient(qfn), gp=gradient(qfp), gd=gradient(distance);
    const std::array<Real,2> normal{gd[0].value,gd[1].value};
    const D en=detail::dualAbs(dot(e,normal));
    const D normalComponent=dot(e,normal);
    const D parallel=norm({e[0]-normalComponent*D(normal[0]),e[1]-normalComponent*D(normal[1])});
    D driveN=geometry.partialBoundaryLayer?detail::dualAbs(dot(gn,geometry.boundaryTangent)):norm(gn);
    D driveP=geometry.partialBoundaryLayer?detail::dualAbs(dot(gp,geometry.boundaryTangent)):norm(gp);
    if (geometry.touchesEffectiveElectrode) driveN=driveP=norm(e);
    IalElementMobilityResult result;
    for (int i=0;i<3;++i) {
        const IalMobilityState local{state[i].donors_m3,state[i].acceptors_m3,n[i].value,p[i].value,
            en.value,geometry.interfaceDistance_m[i],state[i].temperature_K};
        const auto low=[&](const IalMobility& model,int carrier) {
            if(!options.spatialDerivatives)return D(model.evaluate(local,options.screeningCache,
                options.preparationCache).mobility_m2_per_Vs);
            auto& saved=localDifferentials[2*i+carrier];
            if(!thermalDirections || !options.reuseThermalLocalDifferentials)
                saved=model.evaluateWithDerivatives(local,options.screeningCache,options.preparationCache);
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
                    auto& h=highFieldDifferentials[2*i+carrier];
                    if(!thermalDirections || !options.reuseThermalHighField)
                        h=evaluateIalHighFieldMobility(mu.value,field.value,state[i].temperature_K,
                            {base.saturationVelocity,base.beta,ve,be});
                    else ++ialKernelProfile.highFieldReuses;
                    D result(h.mobility_m2_per_Vs);
                    for(int k=0;k<9;++k)result.derivative[k]=h.lowFieldDerivative*mu.derivative[k]+
                        h.drivingFieldDerivative_m3_per_V2s*field.derivative[k]+h.temperatureDerivative_m2_per_Vs_K*temperature[i].derivative[k];
                    return result;
                };
                finalN=hot(lowN,fieldN,options.electronField,options.electronVelocityTemperatureExponent,options.electronBetaTemperatureExponent,0);
                finalP=hot(lowP,fieldP,options.holeField,options.holeVelocityTemperatureExponent,options.holeBetaTemperatureExponent,1);
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
} // namespace vela
