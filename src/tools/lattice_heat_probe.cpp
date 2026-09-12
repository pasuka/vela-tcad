// Prescribed-source thermal diagnostic. Does not claim electrothermal closure.
#include "vela/equation/LatticeHeatAssembler.h"
#include "vela/io/MeshReader.h"
#include <Eigen/SparseLU>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <numeric>
#include <stdexcept>

using namespace vela;
using json=nlohmann::json;
namespace {
VectorXd values(const json& a) {
    const auto data=a.get<std::vector<Real>>();
    return Eigen::Map<const VectorXd>(data.data(),data.size());
}
std::vector<Real> list(const VectorXd& a) {return {a.data(),a.data()+a.size()};}
}
int main(int argc,char** argv) {
    try {
        if (argc!=3) throw std::invalid_argument("Usage: lattice_heat_probe INPUT.json OUTPUT.json");
        if (std::filesystem::exists(argv[2])) throw std::invalid_argument("Output already exists");
        std::ifstream input(argv[1]);if (!input) throw std::invalid_argument("Cannot open input");
        json cfg;input>>cfg;
        JsonMeshReader reader;
        auto mesh=reader.read(cfg.at("mesh_file").get<std::string>());
        std::map<Index,LatticeConductivity> laws;
        for (const auto& r:cfg.at("region_conductivity")) {
            LatticeConductivity law;
            const auto model=r.at("model").get<std::string>();
            if (model=="constant") law.constant_W_per_m_K=r.at("value_W_per_m_K");
            else if (model=="inverse_quadratic") {
                law.model=LatticeConductivity::Model::InverseQuadratic;
                law.numerator=r.at("numerator");law.denominator=r.at("denominator").get<std::array<Real,3>>();
            } else throw std::invalid_argument("Unknown conductivity model");
            if (!laws.emplace(r.at("region_id").get<Index>(),law).second)
                throw std::invalid_argument("Duplicate region conductivity");
        }
        std::vector<LatticeThermodeEdge> edges;
        for (const auto& b:cfg.at("thermodes"))
            edges.push_back({b.at("nodes").get<std::array<Index,2>>(),b.at("ambient_K"),b.at("conductance_W_per_m2_K")});
        LatticeHeatAssembler assembler(mesh,cfg.at("coordinate_to_metres").get<Real>(),laws,edges);
        VectorXd t=values(cfg.at("temperature_K")),q=values(cfg.at("cell_source_W_per_m3"));
        const Real tolerance=cfg.value("absolute_residual_W_per_m",1e-8);
        if (!std::isfinite(tolerance) || tolerance<=0) throw std::invalid_argument("Invalid tolerance");
        json history=json::array();
        auto a=assembler.assemble(t,q);
        const bool solve=cfg.value("solve_prescribed_source",false);
        bool converged=a.residual_W_per_m.lpNorm<Eigen::Infinity>()<=tolerance;
        for (unsigned iteration=0; solve && !converged && iteration<50; ++iteration) {
            Eigen::SparseLU<SparseMatrixd> lu;lu.compute(a.jacobian_W_per_m_K);
            if (lu.info()!=Eigen::Success) throw std::runtime_error("Thermal factorization failed");
            VectorXd update=lu.solve(-a.residual_W_per_m);
            if (lu.info()!=Eigen::Success || !update.allFinite()) throw std::runtime_error("Thermal solve failed");
            bool accepted=false;Real alpha=1.;
            for (unsigned backtrack=0;backtrack<30;++backtrack,alpha*=.5) {
                const VectorXd trial=t+alpha*update;
                if (trial.minCoeff()<=0) continue;
                auto next=assembler.assemble(trial,q);
                if (next.residual_W_per_m.norm()<a.residual_W_per_m.norm()) {
                    t=trial;a=std::move(next);accepted=true;break;
                }
            }
            if (!accepted) throw std::runtime_error("Thermal line search failed");
            history.push_back({{"iteration",iteration+1},{"alpha",alpha},{"residual_inf_W_per_m",a.residual_W_per_m.lpNorm<Eigen::Infinity>()}});
            converged=a.residual_W_per_m.lpNorm<Eigen::Infinity>()<=tolerance;
        }
        std::vector<Index> nodeIds(mesh.numNodes());
        std::iota(nodeIds.begin(),nodeIds.end(),Index{0});
        json result={
            {"schema","vela.prescribed_lattice_heat_probe.v1"},
            {"scope","prescribed heat sources only; electrical state is not solved"},
            {"solved",solve},{"thermal_residual_converged",converged},
            {"node_id",nodeIds},
            {"temperature_K",list(t)},{"nodal_area_m2",list(assembler.nodalAreas_m2())},
            {"residual_W_per_m",list(a.residual_W_per_m)},
            {"integrated_source_W_per_m",a.integrated_source_W_per_m},
            {"outward_boundary_heat_W_per_m",a.outward_boundary_heat_W_per_m},
            {"history",history}};
        std::ofstream output(argv[2]);if (!output) throw std::runtime_error("Cannot open output");
        output<<result.dump(2)<<'\n';
        return solve && !converged ? 2:0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
