#pragma once
#include <cmath>
#include <cstddef>
#include <nlohmann/json.hpp>

namespace vela::update_diagnostics {
// Read-only observer. Subtract partitions before summing: forming ref+inc
// first would hide precisely the sub-ULP writeback this observer measures.
inline long double coordinateChange(double ref0,double inc0,double ref1,double inc1) {
    return (static_cast<long double>(ref1)-ref0)+
           (static_cast<long double>(inc1)-inc0);
}
struct Rms {
    long double rawSum=0,actualSum=0;
    std::size_t count=0;
    double rawMax=0,actualMax=0;
    int rawNode=-1,actualNode=-1;
    void add(double raw,long double actual,double denominator,int node) {
        const double r=raw/denominator*1e5;
        const double a=static_cast<double>(actual/denominator*1e5);
        rawSum+=static_cast<long double>(r)*r;
        actualSum+=static_cast<long double>(a)*a;++count;
        if(std::abs(r)>rawMax){rawMax=std::abs(r);rawNode=node;}
        if(std::abs(a)>actualMax){actualMax=std::abs(a);actualNode=node;}
    }
    void merge(const Rms& x){
        rawSum+=x.rawSum;actualSum+=x.actualSum;count+=x.count;
        if(x.rawMax>rawMax){rawMax=x.rawMax;rawNode=x.rawNode;}
        if(x.actualMax>actualMax){actualMax=x.actualMax;actualNode=x.actualNode;}
    }
    nlohmann::json json() const {
        return {{"count",count},{"raw_rms",count?std::sqrt(static_cast<double>(rawSum/count)):0.},
            {"actual_rms",count?std::sqrt(static_cast<double>(actualSum/count)):0.},
            {"raw_max",rawMax},{"raw_max_node",rawNode},
            {"actual_max",actualMax},{"actual_max_node",actualNode}};
    }
};
inline nlohmann::json stopLayers(double rhs,double rawU,bool originalPolicy) {
    return {{"rhs",rhs},{"rhs_satisfied",std::isfinite(rhs)&&rhs>=0.&&rhs<=1e-10},
        {"raw_U",rawU},{"raw_U_satisfied",std::isfinite(rawU)&&rawU>=0.&&rawU<1.},
        {"original_policy_satisfied",originalPolicy},{"actual_U_used_for_acceptance",false}};
}
}
