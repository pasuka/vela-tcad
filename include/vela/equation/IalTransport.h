#pragma once
#include "vela/equation/IalElementMobility.h"
#include "vela/equation/IalElementAvalanche.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include <map>
#include <optional>

namespace vela {
struct IalTransportOptions {
    std::string geometryFile;
    /// Optional, mesh-specific observed selector and independent source shares.
    std::string avalancheGeometryFile;
    std::vector<std::string> effectiveElectrodes;
    std::array<Real,3> crystalX{1.,0.,0.},crystalY{0.,1.,0.};
    std::map<int,IalMobility> electrons,holes;
    std::map<int,BasicIalMobilityParameters<long double>> extendedElectrons,extendedHoles;
    IalElementMobilityOptions element;
    // Explicit non-Si transport controls. Values are cm^2/(V s); no HFS or
    // element avalanche is applied to these materials.
    std::map<std::string,std::array<Real,2>> constantMaterialMobility;
};
struct IalTransportGeometry {
    struct CellSupport {
        Index cellId;
        IalElementGeometry geometry;
        std::array<int,3> family;
        std::optional<IalElementAvalancheOptions> avalanche;
    };
    struct EdgeContribution { std::size_t support; Real weight; };
    std::vector<CellSupport> cells;
    // -1 for cells outside the silicon transport support.
    std::vector<int> cellSupportIndex;
    std::vector<std::vector<EdgeContribution>> edges;
    std::vector<Real> constantElectronEdges,constantHoleEdges;
};
// Measures follow the cell and edge support order of validated IalTransportGeometry.
struct ExtendedIalTransportMeasures {
    std::vector<std::array<long double, 3>> vertex, avalanche;
    std::vector<std::vector<long double>> edgeWeights;
    std::optional<ArithmeticIalIonizationCoefficients> ionization;
    struct SurfaceCell {
        std::array<long double, 3> distances;
        std::array<long double, 2> tangent;
        std::array<int, 3> family;
    };
    std::vector<SurfaceCell> surface;
    std::vector<std::array<long double, 2>> constantMobility_m2_per_Vs;
};
ExtendedIalTransportMeasures readExtendedIalTransportMeasures(
    const MobilityModelConfig& config, const DeviceMesh& mesh,
    std::optional<long double> lengthScale_m = {},
    const std::map<std::string, std::array<long double, 2>>& constantMobility_cm2_per_Vs = {});

struct IalTransportState {
    bool hasTemperatureDerivatives=true;
    bool hasSpatialDerivatives=true;
    VectorXd psi,n,p,phin,phip,dn,dp,temperature,dn_dT,dp_dT;
    std::vector<IalElementMobilityResult> cells;
    std::vector<Real> electronEdges,holeEdges;
};

std::shared_ptr<const IalTransportOptions> ialTransportOptionsFromJson(
    const nlohmann::json& input, const nlohmann::json* mobility = nullptr,
    UnitScalingConfig scaling = {});
// Prepare once for an immutable mesh/options pair, then share the read-only
// geometry through config copies. Live carrier states are not shared here.
void prepareIalTransportGeometry(MobilityModelConfig& config,const DeviceMesh& mesh);
// Evaluate the supplied local state, never the cached global state. Inputs are
// SI; returned mobilities use the configured internal units, partials per volt.
// Used by cell-local avalanche residuals and independent perturbed-state audits.
IalElementMobilityResult evaluateIalTransportCell(
    const MobilityModelConfig& config,Index cellId,
    const std::array<IalElementVertexState,3>& states,bool spatialDerivatives,
    Real populationThermalVoltage_V=0.);
/// Controlled density-gradient source, SI input and output. No cached state.
IalElementAvalancheResult evaluateIalTransportAvalancheCell(
    const MobilityModelConfig& config,Index cellId,
    const std::array<IalElementVertexState,3>& states,
    IalElementAvalancheOptions options);
struct ExtendedIalTransportCell {
    std::array<detail::IalPopulationDual,2> mobility_m2_per_Vs;
    std::optional<BasicIalElementAvalancheResult<ArithmeticReal>> avalanche;
};
/// Reassemble one silicon cell from the bound PAR, mesh, and population state.
/// The caller supplies the scales of the authoritative input coordinates.
ExtendedIalTransportCell evaluateExtendedIalTransportCell(
    const MobilityModelConfig& config, const DeviceMesh& mesh, const DopingModel& doping,
    Index cellId, const std::array<IalPopulationVertex,3>& state,
    ArithmeticReal lengthScale_m, ArithmeticReal dopingScale_m3,
    ArithmeticReal potentialScale_V, ArithmeticReal thermalVoltage_V,
    ArithmeticReal temperature_K, std::optional<IalElementAvalancheOptions> avalanche = {},
    const ExtendedIalTransportMeasures* measures = nullptr);
void updateIalTransportState(MobilityModelConfig& config,const DeviceMesh& mesh,
    const DopingModel& doping,const VectorXd& psi_V,const VectorXd& n,const VectorXd& p,
    const VectorXd& phin_V,const VectorXd& phip_V,
    const VectorXd& dn_per_V = {},const VectorXd& dp_per_V = {},
    const VectorXd& temperature_K = {},const VectorXd& dn_per_K = {},const VectorXd& dp_per_K = {},
    bool reuseScreening=false,bool temperatureDerivatives=true);
Real ialEdgeMobility(const MobilityModelConfig& config,Index edge,CarrierType carrier,bool lowField=false);
} // namespace vela
