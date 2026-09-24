#include "vela/equation/LatticeHeatAssembler.h"
#include <algorithm>
#include <bit>
#include <cmath>
#include <set>
#include <stdexcept>

namespace vela {
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
using EdgeKey = std::array<Index, 2>;
EdgeKey edgeKey(Index i, Index j) { return {std::min(i,j), std::max(i,j)}; }
}

std::pair<Real, Real> LatticeConductivity::valueAndDerivative(Real t) const {
    require(std::isfinite(t) && t > 0, "Lattice temperature must be finite and positive");
    if (model == Model::Constant) {
        require(std::isfinite(constant_W_per_m_K) && constant_W_per_m_K > 0,
                "Lattice conductivity must be finite and positive");
        return {constant_W_per_m_K, 0};
    }
    require(model == Model::InverseQuadratic, "Unknown lattice conductivity model");
    require(std::isfinite(numerator) && numerator > 0, "Invalid conductivity numerator");
    for (Real a : denominator) require(std::isfinite(a), "Invalid conductivity coefficient");
    const Real d = denominator[0] + t * (denominator[1] + t * denominator[2]);
    require(std::isfinite(d) && d > 0, "Conductivity denominator must be finite and positive");
    const Real k = numerator / d;
    const Real dk = -k * (denominator[1] + 2 * t * denominator[2]) / d;
    require(std::isfinite(k) && k > 0 && std::isfinite(dk), "Non-finite lattice conductivity");
    return {k, dk};
}

LatticeHeatAssembler::LatticeHeatAssembler(
    const DeviceMesh& mesh, Real scale,
    const std::map<Index, LatticeConductivity>& laws,
    const std::vector<LatticeThermodeEdge>& thermodes)
    : nodalAreas_(VectorXd::Zero(mesh.numNodes())) {
    require(std::isfinite(scale) && scale > 0, "Thermal coordinate scale must be positive");
    require(mesh.numNodes() > 0 && mesh.numCells() > 0, "Thermal mesh must not be empty");
    std::map<EdgeKey, unsigned> edgeCount;
    for (const auto& cell : mesh.cells()) {
        require(cell.type == CellType::Tri3 && cell.node_ids.size() == 3,
                "Lattice heat requires Tri3 cells");
        require(cell.id == elements_.size(), "Thermal cell IDs must be contiguous");
        const auto it = laws.find(cell.region_id);
        require(it != laws.end(), "Every thermal region requires an explicit conductivity law");
        Element e{{cell.node_ids[0],cell.node_ids[1],cell.node_ids[2]},0,{},it->second};
        std::array<Point2,3> p;
        for (int a=0; a<3; ++a) {
            require(e.nodes[a] < mesh.numNodes(), "Invalid thermal node ID");
            const auto& n = mesh.getNode(e.nodes[a]);
            p[a] = Point2(n.x*scale,n.y*scale);
            require(p[a].allFinite(), "Non-finite thermal mesh coordinates");
        }
        const auto u = p[1]-p[0], v = p[2]-p[0];
        const Real det = u.x()*v.y()-u.y()*v.x();
        require(std::isfinite(det) && det != 0, "Degenerate thermal triangle");
        e.area = std::abs(det)*0.5;
        Eigen::Matrix<Real,2,3> grad;
        for (int a=0; a<3; ++a) {
            const int j=(a+1)%3, k=(a+2)%3;
            grad.col(a)=Point2(p[j].y()-p[k].y(),p[k].x()-p[j].x())/det;
            nodalAreas_[e.nodes[a]] += e.area/3;
            require(++edgeCount[edgeKey(e.nodes[a],e.nodes[j])] <= 2,
                    "Non-manifold thermal edge");
        }
        e.gradientIntegral=e.area*grad.transpose()*grad;
        require(e.gradientIntegral.allFinite(), "Non-finite thermal geometry coefficients");
        elements_.push_back(e);
    }
    require((nodalAreas_.array()>0).all(), "Unconnected thermal node");
    std::set<EdgeKey> assigned;
    for (const auto& b : thermodes) {
        const auto key=edgeKey(b.nodes[0],b.nodes[1]);
        require(edgeCount.contains(key) && edgeCount.at(key)==1,
                "Thermode must reference an exterior mesh edge");
        require(assigned.insert(key).second, "Duplicate thermal boundary edge");
        require(std::isfinite(b.ambient_K) && b.ambient_K>0,
                "Thermode ambient must be finite and positive");
        require(std::isfinite(b.conductance_W_per_m2_K) && b.conductance_W_per_m2_K>0,
                "Thermode conductance must be finite and positive");
        const auto& a=mesh.getNode(key[0]); const auto& z=mesh.getNode(key[1]);
        boundaries_.push_back({b,std::hypot(a.x-z.x,a.y-z.y)*scale});
    }
}

std::vector<std::uint64_t> LatticeHeatAssembler::structureIdentity() const {
    std::vector<std::uint64_t> key{static_cast<std::uint64_t>(nodalAreas_.size()),elements_.size(),boundaries_.size()};
    const auto real=[&](Real x){key.push_back(std::bit_cast<std::uint64_t>(x));};
    for(const auto& e:elements_){
        for(auto i:e.nodes)key.push_back(i);
        real(e.area);for(int i=0;i<9;++i)real(e.gradientIntegral.data()[i]);
        key.push_back(static_cast<std::uint64_t>(e.law.model));real(e.law.constant_W_per_m_K);
        real(e.law.numerator);for(auto v:e.law.denominator)real(v);
    }
    for(const auto& b:boundaries_){for(auto i:b.input.nodes)key.push_back(i);
        real(b.length);real(b.input.ambient_K);real(b.input.conductance_W_per_m2_K);}
    return key;
}
void LatticeHeatAssembler::appendStructure(std::vector<Eigen::Triplet<Real>>& entries,Index stride,Index offset) const {
    for(const auto& e:elements_)for(auto i:e.nodes)for(auto j:e.nodes)entries.emplace_back(stride*i+offset,stride*j+offset,0.);
    for(const auto& b:boundaries_)for(auto i:b.input.nodes)for(auto j:b.input.nodes)entries.emplace_back(stride*i+offset,stride*j+offset,0.);
}
LatticeHeatAssembly LatticeHeatAssembler::assemble(const VectorXd& t, const VectorXd& q) const {
    require(t.size()==nodalAreas_.size() && t.allFinite() && (t.array()>0).all(),
            "Invalid nodal lattice temperature vector");
    require(q.size()==static_cast<Eigen::Index>(elements_.size()) && q.allFinite(),
            "Invalid prescribed cell heat source vector");
    LatticeHeatAssembly result;
    result.residual_W_per_m=VectorXd::Zero(t.size());
    std::vector<Eigen::Triplet<Real>> entries;
    std::size_t cursor=0;
    if(structure_){
        auto key=structureIdentity();
        if(!structure_->matches(key)){
            appendStructure(entries);structure_->build(std::move(key),t.size(),entries);entries.clear();
        }else ++structure_->hits;
        result.jacobian_W_per_m_K=structure_->zeroMatrix();
    }else entries.reserve(elements_.size()*9+boundaries_.size()*4);
    const auto add=[&](Index i,Index j,Real value){
        if(structure_)structure_->add(result.jacobian_W_per_m_K,cursor,i,j,value);
        else entries.emplace_back(i,j,value);
    };
    for (Index c=0; c<elements_.size(); ++c) {
        const auto& e=elements_[c];
        Eigen::Vector3d local(t[e.nodes[0]],t[e.nodes[1]],t[e.nodes[2]]);
        const auto [k,dk]=e.law.valueAndDerivative(local.mean());
        // Subtract a local reference to preserve exact constant-temperature nullspace.
        const Eigen::Vector3d gradient=e.gradientIntegral*(local.array()-local[0]).matrix();
        result.integrated_source_W_per_m += q[c]*e.area;
        for (int a=0; a<3; ++a) {
            result.residual_W_per_m[e.nodes[a]] += k*gradient[a]-q[c]*e.area/3;
            for (int b=0; b<3; ++b)
                add(e.nodes[a],e.nodes[b],k*e.gradientIntegral(a,b)+dk*gradient[a]/3);
        }
    }
    for (const auto& boundary : boundaries_) {
        const auto& b=boundary.input;
        const Real factor=b.conductance_W_per_m2_K*boundary.length/6;
        const Real d0=t[b.nodes[0]]-b.ambient_K, d1=t[b.nodes[1]]-b.ambient_K;
        result.residual_W_per_m[b.nodes[0]] += factor*(2*d0+d1);
        result.residual_W_per_m[b.nodes[1]] += factor*(d0+2*d1);
        result.outward_boundary_heat_W_per_m += 3*factor*(d0+d1);
        for (int i=0;i<2;++i) for (int j=0;j<2;++j)
            add(b.nodes[i],b.nodes[j],factor*(i==j?2:1));
    }
    if(!structure_){result.jacobian_W_per_m_K.resize(t.size(),t.size());
        result.jacobian_W_per_m_K.setFromTriplets(entries.begin(),entries.end());}
    return result;
}
} // namespace vela
