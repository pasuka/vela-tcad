#pragma once
#include "vela/numerics/SplitCoordinate.h"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <nlohmann/json.hpp>
#include <vector>
#include <string>

namespace vela::split_dd {
// References can be O(1 V) while their increments are O(1e-25 V). Quad
// precision alone would lose an increment's sub-ULP response upon reference
// addition. Keep 100 decimal digits throughout this qualification operator.
using Wide=boost::multiprecision::cpp_bin_float_100;
// Explicit experimental checkpoint contract. Ordinary DDSolution readers must
// not consume this schema by dropping its low components.
class State {
    std::vector<numerics::SplitCoordinate> coordinates_;
    std::vector<double> electronReference_,holeReference_;
    double scale_;
    std::string mesh_;
public:
    State(std::vector<numerics::SplitCoordinate> x,std::vector<double> en,
          std::vector<double> hp,double scale,std::string mesh)
        :coordinates_(std::move(x)),electronReference_(std::move(en)),
         holeReference_(std::move(hp)),scale_(scale),mesh_(std::move(mesh)) {
        if(electronReference_.empty() || holeReference_.size()!=electronReference_.size() ||
           coordinates_.size()!=3*electronReference_.size() || !std::isfinite(scale_) ||
           scale_<=0 || mesh_.empty())throw std::invalid_argument("Invalid split DD state layout");
        for(const auto& p:coordinates_)if(!std::isfinite(p.hi)||!std::isfinite(p.lo))
            throw std::invalid_argument("Nonfinite split coordinate");
        for(double r:electronReference_)if(!std::isfinite(r))throw std::invalid_argument("Nonfinite QF reference");
        for(double r:holeReference_)if(!std::isfinite(r))throw std::invalid_argument("Nonfinite QF reference");
    }
    std::size_t size() const{return electronReference_.size();}
    double potentialScale() const{return scale_;}
    double reference(std::size_t block,std::size_t node) const {
        if(block==1)return electronReference_.at(node);
        if(block==2)return holeReference_.at(node);
        throw std::out_of_range("QF reference block");
    }
    const std::string& meshFingerprint() const{return mesh_;}
    Wide coordinate(std::size_t i) const {const auto& p=coordinates_.at(i);return Wide(p.hi)+Wide(p.lo);}
    Wide potential(std::size_t block,std::size_t node) const {
        if(block>2 || node>=size())throw std::out_of_range("Split potential index");
        return coordinate(block*size()+node)*Wide(scale_)+
            (block==1?Wide(electronReference_[node]):block==2?Wide(holeReference_[node]):Wide(0));
    }
    State shifted(const std::vector<double>& step,double alpha) const {
        if(step.size()!=coordinates_.size() || !std::isfinite(alpha))
            throw std::invalid_argument("Invalid split DD update");
        auto result=*this;
        for(std::size_t i=0;i<step.size();++i) {
            const double product=alpha*step[i];
            const double tail=std::fma(alpha,step[i],-product);
            result.coordinates_[i]=coordinates_[i].shifted(product).shifted(tail);
        }
        return result;
    }
    nlohmann::json checkpoint() const {
        nlohmann::json pairs=nlohmann::json::array();
        for(const auto& p:coordinates_)pairs.push_back({p.hi,p.lo});
        return {{"schema","vela.split-dd-state.v1"},{"mesh_fingerprint",mesh_},
            {"potential_scale_V",scale_},{"electron_reference_V",electronReference_},
            {"hole_reference_V",holeReference_},{"coordinates",pairs}};
    }
    static State restore(const nlohmann::json& j,const std::string& mesh) {
        if(j.at("schema")!="vela.split-dd-state.v1" || j.at("mesh_fingerprint")!=mesh)
            throw std::invalid_argument("Split checkpoint schema or mesh mismatch");
        std::vector<numerics::SplitCoordinate> x;
        for(const auto& p:j.at("coordinates")) {
            if(!p.is_array() || p.size()!=2)throw std::invalid_argument("Incomplete split pair");
            x.push_back({p[0].get<double>(),p[1].get<double>()});
        }
        return State(x,j.at("electron_reference_V").get<std::vector<double>>(),
            j.at("hole_reference_V").get<std::vector<double>>(),j.at("potential_scale_V"),mesh);
    }
};
}
