// Private fixed-state audit only; no production CMake target includes this.
// Geometry is independently generated from mesh coordinates. No native E/G.
#pragma once
#include <array>
#include <cmath>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <set>
#include <nlohmann/json.hpp>

namespace bvmethods_live {
struct Entry { int row,column; double ex,ey; }; // V/cm per physical V
struct Geometry {
    std::vector<double> volume_um2;
    std::vector<Entry> entries;
    explicit Geometry(const nlohmann::json& data) {
        volume_um2=data.at("volume_um2").get<std::vector<double>>();
        if(volume_um2.empty())throw std::invalid_argument("empty geometry");
        for(double v:volume_um2)if(!std::isfinite(v)||v<0)throw std::invalid_argument("invalid volume");
        std::set<std::pair<int,int>> seen;
        for(const auto& r:data.at("stencil")) {
            Entry e{r.at("row"),r.at("column"),r.at("Ex"),r.at("Ey")};
            if(e.row<0||e.column<0||e.row>=int(volume_um2.size())||e.column>=int(volume_um2.size())||
               !std::isfinite(e.ex)||!std::isfinite(e.ey)||!seen.insert({e.row,e.column}).second)
                throw std::invalid_argument("invalid field stencil");
            entries.push_back(e);
        }
    }
    std::vector<std::array<double,2>> field(const std::vector<double>& psi) const {
        if(psi.size()!=volume_um2.size())throw std::invalid_argument("state size");
        for(double v:psi)if(!std::isfinite(v))throw std::invalid_argument("nonfinite potential");
        std::vector<std::array<long double,2>> sum(psi.size());
        for(const auto& e:entries) {
            if(e.row==e.column)continue;
            const long double drop=(long double)psi[e.column]-(long double)psi[e.row];
            sum[e.row][0]+=e.ex*drop;sum[e.row][1]+=e.ey*drop;
        }
        std::vector<std::array<double,2>> out(psi.size());
        for(std::size_t i=0;i<out.size();++i)out[i]={double(sum[i][0]),double(sum[i][1])};
        return out;
    }
    template<class Rate> std::vector<double> integrals(const std::vector<double>& psi,Rate rate) const {
        const auto E=field(psi);std::vector<double> out(E.size());
        for(std::size_t i=0;i<E.size();++i)out[i]=rate(100*std::hypot(E[i][0],E[i][1]))*volume_um2[i];
        return out; // rate(SI field)*um^2; caller owns concentration/area units
    }
    template<class Partial> std::vector<std::tuple<int,int,double>> chain(const std::vector<double>& psi,Partial partial) const {
        const auto E=field(psi);std::vector<std::tuple<int,int,double>> out;
        for(const auto& e:entries) {
            const auto f=E[e.row];const double norm=std::hypot(f[0],f[1]);
            const double value=norm>0 ? partial(norm*100)*volume_um2[e.row]*100*(f[0]*e.ex+f[1]*e.ey)/norm : 0.;
            out.emplace_back(e.row,e.column,value);
        }
        return out; // d(rate*um^2)/d(physical psi V)
    }
};
}
