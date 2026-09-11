#pragma once
#include "vela/equation/IalElementMobility.h"
#include "vela/mesh/DeviceMesh.h"
#include "vela/physics/DopingModel.h"
#include <map>

namespace vela {
struct IalTransportOptions {
    std::string geometryFile;
    std::vector<std::string> effectiveElectrodes;
    std::array<Real,3> crystalX{1.,0.,0.},crystalY{0.,1.,0.};
    std::map<int,IalMobility> electrons,holes;
    IalElementMobilityOptions element;
};
struct IalTransportGeometry {
    struct CellSupport {
        Index cellId;
        IalElementGeometry geometry;
        std::array<int,3> family;
    };
    struct EdgeContribution { std::size_t support; Real weight; };
    std::vector<CellSupport> cells;
    std::vector<std::vector<EdgeContribution>> edges;
};
struct IalTransportState {
    VectorXd psi,n,p,phin,phip,dn,dp;
    std::vector<IalElementMobilityResult> cells;
    std::vector<Real> electronEdges,holeEdges;
};

std::shared_ptr<const IalTransportOptions> ialTransportOptionsFromJson(const nlohmann::json& input);
void updateIalTransportState(MobilityModelConfig& config,const DeviceMesh& mesh,
    const DopingModel& doping,const VectorXd& psi_V,const VectorXd& n,const VectorXd& p,
    const VectorXd& phin_V,const VectorXd& phip_V,
    const VectorXd& dn_per_V = {},const VectorXd& dp_per_V = {});
Real ialEdgeMobility(const MobilityModelConfig& config,Index edge,CarrierType carrier,bool lowField=false);
} // namespace vela
