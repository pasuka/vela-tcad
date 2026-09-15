#pragma once
#include "vela/equation/IalTransport.h"
#include "vela/equation/LatticeHeatAssembler.h"
#include "vela/physics/SiliconThermalPhysics.h"
#include <optional>

namespace vela {
/// Explicit SI operator inputs. Geometry/material preparation remains external.
/// No legacy solver configuration or acceptance criterion is inferred.
struct ElectrothermalGeometry {
    VectorXd poissonEdge_F_per_m; // permittivity times dual/edge length
    VectorXd siliconArea_m2;      // zero at oxide-only nodes
    /// Explicit source volume; empty uses siliconArea_m2. The LDMOS adapter
    /// supplies its qualified global-node source volumes at silicon nodes.
    VectorXd recombinationArea_m2;
    VectorXd fixedCharge_C_per_m;
    VectorXd transportWeight;     // AverageBox dual/edge length, zero in oxide
};
struct ElectrothermalBoundary {
    struct HoleRecombination {
        Real velocity_m_per_s=0.;
        Real boundaryLength_m=0.; // half of each incident contact edge, per unit width
    };
    std::map<Index,Real> potential_V, electronQf_V, holeQf_V, temperature_K;
    /// Ideal neutral silicon contacts. Replaces all three electrical rows;
    /// includes the temperature derivative of the neutrality potential.
    std::map<Index,Real> neutralContactBias_V;
    /// Finite hole exchange at a neutral-potential Ohmic contact. The electron
    /// constraint remains ideal; the hole continuity row remains active.
    /// Velocity is constant: electrode-specific RecVel(TempDep) is not enabled.
    std::map<Index,HoleRecombination> holeRecombination;
};
struct ElectrothermalAssembly {
    /// Interleaved psi, fn, fp, T. Interior units: C/m, A/m, A/m, W/m.
    /// Dirichlet rows have V, V, V, K units respectively.
    VectorXd residual;
    SparseMatrixd jacobian;
    VectorXd electronOutflow_A_per_m, holeOutflow_A_per_m;
    VectorXd electronFluxAbs_A_per_m, holeFluxAbs_A_per_m, recombination_A_per_m;
    Real latticeSource_W_per_m=0., boundaryHeat_W_per_m=0.;
};
/// Experimental four-equation DD/lattice operator, with explicit silicon
/// statistics and the same oriented edge current in continuity and heating.
/// This does not change the qualified isothermal production solver.
class ElectrothermalAssembler {
public:
    ElectrothermalAssembler(const DeviceMesh& mesh, const DopingModel& doping_SI,
        ElectrothermalGeometry geometry, LatticeHeatAssembler heat,
        MobilityModelConfig mobility_SI, SiliconThermalPhysics physics=SiliconThermalPhysics{},
        Real constantElectronMobility=.1, Real constantHoleMobility=.04,
        bool reusePhysicsPreparation=false,bool reuseIalScreening=false);
    ElectrothermalAssembly assemble(const VectorXd& state,
                                   const ElectrothermalBoundary& boundary,
        const VectorXd& electronQfReference_V={}, const VectorXd& holeQfReference_V={},
        bool buildJacobian=true,bool skipEquilibriumTransport=false) const;
    /// psi and dpsi/dT for n-p=Nd-Na, fn=fp=bias.
    std::pair<Real,Real> neutralPotential(Index node, Real bias, Real temperature) const;
    /// Doping preparations, temperature preparations, exact preparation hits.
    std::array<std::size_t,3> preparationCounts() const {return preparationCounts_;}
private:
    const SiliconThermalPhysics::TemperaturePreparation& preparedAt(Index,Real temperature) const;
    const DeviceMesh& mesh_;
    const DopingModel& doping_;
    ElectrothermalGeometry geometry_;
    LatticeHeatAssembler heat_;
    mutable MobilityModelConfig mobility_;
    SiliconThermalPhysics physics_;
    Real muE_,muH_;
    bool reusePreparation_;
    bool reuseIalScreening_;
    mutable std::vector<std::optional<SiliconThermalPhysics::DopingPreparation>> dopingPreparation_;
    mutable std::vector<std::optional<SiliconThermalPhysics::TemperaturePreparation>> temperaturePreparation_;
    mutable std::array<std::size_t,3> preparationCounts_{};
};
} // namespace vela
