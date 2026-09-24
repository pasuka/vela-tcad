#pragma once

#include "vela/mesh/DeviceMesh.h"
#include "vela/equation/SparseAssemblyStructure.h"
#include <map>
#include <memory>

namespace vela {

/// Explicit SI conductivity. No material-name-based defaults are inferred.
struct LatticeConductivity {
    enum class Model { Constant, InverseQuadratic };
    Model model = Model::Constant;
    Real constant_W_per_m_K = 1.0;
    /// k(T) = numerator / (a + b*T + c*T*T), in W/(m K).
    Real numerator = 1.0;
    std::array<Real, 3> denominator{1.0, 0.0, 0.0};
    std::pair<Real, Real> valueAndDerivative(Real temperature_K) const;
};

/// Explicit exterior edge, independent of electrical contacts. Others insulate.
struct LatticeThermodeEdge {
    std::array<Index, 2> nodes{};
    Real ambient_K = 300.0;
    Real conductance_W_per_m2_K = 0.0;
};

struct LatticeHeatAssembly {
    VectorXd residual_W_per_m;
    SparseMatrixd jacobian_W_per_m_K;
    Real integrated_source_W_per_m = 0.0;
    Real outward_boundary_heat_W_per_m = 0.0;
};

/// Steady Tri3 heat conduction with prescribed cell heat sources, per unit width.
/// This building block does not implement self-consistent electrical heating.
/// Coordinates are converted by the explicit coordinateToMetres factor; all
/// temperature, conductivity, source and boundary inputs/outputs are SI.
class LatticeHeatAssembler {
public:
    LatticeHeatAssembler(const DeviceMesh& mesh, Real coordinateToMetres,
                         const std::map<Index, LatticeConductivity>& regionLaws,
                         const std::vector<LatticeThermodeEdge>& thermodes);
    LatticeHeatAssembly assemble(const VectorXd& temperature_K,
                                const VectorXd& cell_source_W_per_m3) const;
    const VectorXd& nodalAreas_m2() const { return nodalAreas_; }
    void setStructureCache(std::shared_ptr<SparseAssemblyStructure> cache) { structure_=std::move(cache); }
    std::vector<std::uint64_t> structureIdentity() const;
    void appendStructure(std::vector<Eigen::Triplet<Real>>& entries,Index stride=1,Index offset=0) const;

private:
    struct Element {
        std::array<Index, 3> nodes;
        Real area;
        Eigen::Matrix3d gradientIntegral;
        LatticeConductivity law;
    };
    struct Boundary { LatticeThermodeEdge input; Real length; };
    std::vector<Element> elements_;
    std::vector<Boundary> boundaries_;
    VectorXd nodalAreas_;
    std::shared_ptr<SparseAssemblyStructure> structure_;
};

} // namespace vela
