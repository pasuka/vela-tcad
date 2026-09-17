// Controlled formula microbenchmark; no screening-root or cache lookup timing.
// A synthetic, reproducible coverage set is not a device workload distribution.
#include "vela/physics/detail/IalMobilityGenerated.h"
#include "vela/physics/IalHighFieldMobility.h"
#include <nlohmann/json.hpp>
#include <chrono>
#include <ctime>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <type_traits>

using namespace vela;
using namespace vela::ial_detail;
using json=nlohmann::json;
namespace {
volatile Real sink=0.;
struct Sample {
    std::array<Real,7> x;
    std::array<Dual,7> ad;
    IalMobilityParameters p;
    bool electron;
    Real minimum;
    Preparation<Real> scalar,generatedScalar;
    Preparation<Dual> partial,generatedPartial;
};
template<class T> Real checksum(const Preparation<T>& q) {
    Real sum=0.;
    for(const auto& v:std::array<T,10>{q.mu3,q.mu2,q.transitionScale,q.coulombDamping,q.damping,
            q.phonon3d,q.phononNumerator,q.phononTemperaturePower,q.roughnessExponent,q.roughnessNumerator}) {
        sum+=value(v);
        if constexpr(std::is_same_v<T,Dual>)for(int k=0;k<7;++k)sum+=(k+1)*v.derivative[k];
    }
    return sum;
}
template<class F> json measure(const std::vector<Sample>& samples,int rounds,const char* name,F f) {
    // Warm this path outside the measured region. Keep a checksum observable.
    Real sum=0.;for(const auto& s:samples)sum+=f(s);sink=sum;
    const auto begin=std::chrono::steady_clock::now();const auto cpu=std::clock();sum=0.;
    for(int r=0;r<rounds;++r)for(const auto& s:samples)sum+=f(s);
    sink=sum;
    if(!std::isfinite(sum))throw std::runtime_error("Nonfinite benchmark checksum");
    return {{"name",name},{"calls",rounds*samples.size()},{"checksum",sum},
        {"wall_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()},
        {"cpu_seconds",double(std::clock()-cpu)/CLOCKS_PER_SEC}};
}
}
int main(int argc,char** argv) {
    try {
        const int rounds=argc>1?std::stoi(argv[1]):100;
        const bool reverse=argc>2 && std::string(argv[2])=="reverse";
        if(rounds<1 || rounds>100000)throw std::invalid_argument("rounds outside [1,100000]");
        std::vector<Sample> samples;IalScreeningCache roots;
        for(bool electron:{true,false})for(int k=0;k<256;++k) {
            Sample s;s.electron=electron;s.p=IalMobility::siliconDefaults(electron);
            s.x={std::pow(10.,18.+(k%8)),std::pow(10.,18.+((k/8)%8)),
                 std::pow(10.,15.+(k%11)),std::pow(10.,15.+((k*7)%11)),
                 std::pow(10.,2.+(k%7)),(k%8)*2e-9,200.+(k%9)*50.};
            for(int j=0;j<7;++j){s.ad[j]=Dual(s.x[j]);s.ad[j].derivative[j]=1.;}
            s.minimum=roots.minimum(s.p.mass,s.x[6]);
            s.scalar=prepare(s.x,s.p,s.electron,s.minimum);
            s.partial=prepare(s.ad,s.p,s.electron,s.minimum);
            s.generatedScalar=generatedPrepareValue(s.x,s.p,s.electron,s.minimum);
            s.generatedPartial=generatedPreparePartials(s.x,s.p,s.electron,s.minimum);
            const auto a=evaluatePrepared(s.partial,s.ad[4],s.p);
            const auto b=generatedEvaluatePartials(s.generatedPartial,s.x[4],s.p);
            if(std::abs(a[0].value-b[0].value)>2e-12*a[0].value)
                throw std::runtime_error("Generated value mismatch in benchmark set");
            for(int j=0;j<7;++j){
                const Real scale=std::max(s.x[j],j==5?1e-10:1.)/a[0].value;
                if(std::abs(a[0].derivative[j]-b[0].derivative[j])*scale>
                   2e-10*std::abs(a[0].derivative[j])*scale+2e-12)
                    throw std::runtime_error("Generated derivative mismatch in benchmark set");
            }
            samples.push_back(s);
        }
        using Task=std::pair<const char*,Real(*)(const Sample&)>;
        std::vector<Task> tasks{
            {"prepare_value_reference",[](const Sample& s){return checksum(prepare(s.x,s.p,s.electron,s.minimum));}},
            {"prepare_value_generated",[](const Sample& s){return checksum(generatedPrepareValue(s.x,s.p,s.electron,s.minimum));}},
            {"prepare_partials_reference",[](const Sample& s){const auto p=prepare(s.ad,s.p,s.electron,s.minimum);return checksum(p)+p.mu3.derivative[6];}},
            {"prepare_partials_generated",[](const Sample& s){const auto p=generatedPreparePartials(s.x,s.p,s.electron,s.minimum);return checksum(p)+p.mu3.derivative[6];}},
            {"field_value_reference",[](const Sample& s){return evaluatePrepared(s.scalar,s.x[4],s.p)[0];}},
            {"field_value_generated",[](const Sample& s){return generatedEvaluateValue(s.generatedScalar,s.x[4],s.p)[0];}},
            {"field_partials_reference",[](const Sample& s){const auto p=evaluatePrepared(s.partial,s.ad[4],s.p)[0];Real v=p.value;for(int k=0;k<7;++k)v+=p.derivative[k]*s.x[k];return v;}},
            {"field_partials_generated",[](const Sample& s){const auto p=generatedEvaluatePartials(s.generatedPartial,s.x[4],s.p)[0];Real v=p.value;for(int k=0;k<7;++k)v+=p.derivative[k]*s.x[k];return v;}},
            {"high_partials_reference",[](const Sample& s){const auto p=evaluateIalHighFieldMobility(.03,s.x[4],s.x[6],{});return p.mobility_m2_per_Vs+p.lowFieldDerivative+p.drivingFieldDerivative_m3_per_V2s+p.temperatureDerivative_m2_per_Vs_K;}},
            {"high_partials_explicit",[](const Sample& s){const auto p=evaluateIalHighFieldMobilityExplicit(.03,s.x[4],s.x[6],{});return p.mobility_m2_per_Vs+p.lowFieldDerivative+p.drivingFieldDerivative_m3_per_V2s+p.temperatureDerivative_m2_per_Vs_K;}},
            {"high_value_explicit",[](const Sample& s){return evaluateIalHighFieldMobilityValue(.03,s.x[4],s.x[6],{});}}};
        if(reverse)std::reverse(tasks.begin(),tasks.end());
        json result={{"samples",samples.size()},{"rounds",rounds},{"reverse",reverse},
            {"scope","Synthetic kernel CPU/wall; excludes screening solve, caches, transport chains and sparse solve"},{"measurements",json::array()}};
        for(const auto& task:tasks)result["measurements"].push_back(measure(samples,rounds,task.first,task.second));
        std::cout<<result.dump(2)<<'\n';
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
