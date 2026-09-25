#include "vela/equation/IalTransport.h"
#include "vela/physics/IalInterfaceGeometry.h"
#include "vela/physics/IalMobilityJson.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/IalKernelProfiling.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <set>
#include <numeric>

namespace vela {
namespace {
using json=nlohmann::json;
using Pair=std::pair<Index,Index>;
Pair pair(Index a,Index b) {return std::minmax(a,b);}
bool silicon(const DeviceMesh& mesh,const Cell& c) {return mesh.getRegion(c.region_id).material=="Si";}
void keys(const json& input,std::initializer_list<std::string> allowed) {
    if (!input.is_object()) throw std::invalid_argument("IALMob expects an object");
    for (const auto& [key,v]:input.items())
        if (std::find(allowed.begin(),allowed.end(),key)==allowed.end())
            throw std::invalid_argument("Unsupported IALMob key: "+key);
}
std::shared_ptr<const IalTransportGeometry> geometry(const MobilityModelConfig& config,const DeviceMesh& mesh) {
    ScopedPerformanceTimer timer("ialmob.geometry");
    const auto& options=*config.ialmob;
    std::ifstream stream(options.geometryFile);
    if (!stream) throw std::invalid_argument("Cannot read IALMob transport geometry: "+options.geometryFile);
    const json input=json::parse(stream);
    keys(input,{"schema","node_count","cell_count","cells"});
    if (input.at("schema")!="vela.ialmob.transport_geometry.v1"||
        input.at("node_count").get<Index>()!=mesh.numNodes()||
        input.at("cell_count").get<Index>()!=mesh.numCells())
        throw std::invalid_argument("IALMob geometry does not match the mesh");
    auto result=std::make_shared<IalTransportGeometry>();result->edges.resize(mesh.numEdges());
    std::vector<std::array<Real,2>> xy;
    for (const auto& n:mesh.nodes()) xy.push_back({n.x*config.internalLengthToM,n.y*config.internalLengthToM});
    std::map<Pair,std::vector<Index>> adjacent;std::map<Pair,Index> edgeIds;
    for (const auto& e:mesh.edges()) edgeIds[pair(e.n0,e.n1)]=e.id;
    for (const auto& c:mesh.cells()) {
        if (c.node_ids.size()!=3) throw std::invalid_argument("IALMob transport requires a Tri3 mesh");
        for (int k=0;k<3;++k) adjacent[pair(c.node_ids[k],c.node_ids[(k+1)%3])].push_back(c.id);
    }
    std::set<Index> effectiveNodes;std::set<Pair> contactEdges;std::set<std::string> foundContacts;
    for (const auto& contact:mesh.contacts()) {
        const bool effective=std::find(options.effectiveElectrodes.begin(),options.effectiveElectrodes.end(),contact.name)!=options.effectiveElectrodes.end();
        if (effective) {effectiveNodes.insert(contact.node_ids.begin(),contact.node_ids.end());foundContacts.insert(contact.name);}
        for (const auto& edge:contact.edge_node_ids) contactEdges.insert(pair(edge[0],edge[1]));
        if (contact.edge_node_ids.empty()) {
            const std::set<Index> nodes(contact.node_ids.begin(),contact.node_ids.end());
            for (const auto& [edge,cells]:adjacent)
                if (cells.size()==1&&nodes.count(edge.first)&&nodes.count(edge.second)) contactEdges.insert(edge);
        }
    }
    for (const auto& name:options.effectiveElectrodes)
        if (!foundContacts.count(name)) throw std::invalid_argument("Unknown IALMob effective electrode: "+name);
    std::vector<IalInterfaceSegment> segments;std::map<Index,std::array<Real,2>> tangents;
    std::map<Index,bool> tangentIsMaterialInterface;
    for (const auto& [edge,cells]:adjacent) {
        bool si=false,oxide=false,other=false;Index siCell=0;
        for (Index id:cells) {
            const auto& c=mesh.getCell(id);const auto& material=mesh.getRegion(c.region_id).material;
            if (material=="Si") {si=true;siCell=id;} else other=true;
            if (material=="SiO2") oxide=true;
        }
        if (!si) continue;
        const auto& a=xy[edge.first];const auto& b=xy[edge.second];
        const Real dx=b[0]-a[0],dy=b[1]-a[1],length=std::hypot(dx,dy);
        if (oxide) {
            std::array<Real,2> center{};
            for (Index node:mesh.getCell(siCell).node_ids) for (int d=0;d<2;++d) center[d]+=xy[node][d]/3.;
            segments.push_back({edge.first,edge.second,center});
        }
        if (other||(cells.size()==1&&!contactEdges.count(edge)))
            for (Index id:cells) if (silicon(mesh,mesh.getCell(id))) {
                // At an interface termination, use the material interface for
                // ParallelToInterface, ahead of an ordinary exterior face.
                // Multiple nonparallel faces of equal priority remain unsupported.
                if (tangents.count(id)) {
                    const bool previous=tangentIsMaterialInterface.at(id);
                    if (previous&&!other) continue;
                    if (previous==other) {
                        const auto& t=tangents.at(id);
                        if (std::abs(t[0]*dy/length-t[1]*dx/length)>1e-12)
                            throw std::invalid_argument("IALMob PartialLayer cell has ambiguous boundary interfaces");
                        continue;
                    }
                }
                tangents[id]={dx/length,dy/length};tangentIsMaterialInterface[id]=other;
            }
    }
    const auto nodes=buildIalInterfaceGeometry(xy,segments,options.crystalX,options.crystalY);
    std::set<Index> seenCells;
    std::vector<Real> coefficientSum(mesh.numEdges(),0.);
    for (const auto& row:input.at("cells")) {
        keys(row,{"cell_id","node_ids","vertex_measure_m2","edge_coefficients"});
        const Index id=row.at("cell_id");const auto& c=mesh.getCell(id);
        if (!silicon(mesh,c)||!seenCells.insert(id).second||row.at("node_ids").get<std::vector<Index>>()!=c.node_ids)
            throw std::invalid_argument("IALMob geometry cell topology mismatch");
        IalTransportGeometry::CellSupport support;support.cellId=id;
        support.geometry.vertexMeasure_m2=row.at("vertex_measure_m2");
        Real total=0.;for (Real m:support.geometry.vertexMeasure_m2) {
            if (!std::isfinite(m)||m<0.) throw std::invalid_argument("Invalid IALMob vertex measure");total+=m;
        }
        if (!std::isfinite(total)||total<=0.) throw std::invalid_argument("Invalid IALMob total measure");
        for (int k=0;k<3;++k) {
            const Index node=c.node_ids[k];support.geometry.coordinates_m[k]=xy[node];
            support.geometry.interfaceDistance_m[k]=nodes[node].distance_m;
            support.family[k]=nodes[node].orientationFamily;
            options.electrons.at(support.family[k]);options.holes.at(support.family[k]);
            support.geometry.touchesEffectiveElectrode|=effectiveNodes.count(node)!=0;
        }
        support.geometry.partialBoundaryLayer=tangents.count(id)!=0;
        if (support.geometry.partialBoundaryLayer) support.geometry.boundaryTangent=tangents.at(id);
        const auto coeff=row.at("edge_coefficients").get<std::array<Real,3>>();
        const std::size_t index=result->cells.size();result->cells.push_back(support);
        for (int k=0;k<3;++k) {
            if (!std::isfinite(coeff[k])||coeff[k]<0.) throw std::invalid_argument("Invalid IALMob transport coefficient");
            const Index edge=edgeIds.at(pair(c.node_ids[k],c.node_ids[(k+1)%3]));
            coefficientSum[edge]+=coeff[k];result->edges[edge].push_back({index,coeff[k]});
        }
    }
    for (const auto& c:mesh.cells())
        if (silicon(mesh,c)&&!seenCells.count(c.id)) throw std::invalid_argument("Missing IALMob Si cell measure");
    for (const auto& edge:mesh.edges()) {
        const Real expected=edge.transport_couple>=0.?edge.transport_couple:edge.couple;
        const Real actual=coefficientSum[edge.id]*edge.length;
        if (!result->edges[edge.id].empty()&&std::abs(expected-actual)>1e-10*std::max(std::abs(expected),std::abs(actual)))
            throw std::invalid_argument("IALMob cell coefficients disagree with transport edge coupling");
        if (coefficientSum[edge.id]>0.)
            for (auto& c:result->edges[edge.id]) c.weight/=coefficientSum[edge.id];
    }
    return result;
}
}

std::shared_ptr<const IalTransportOptions> ialTransportOptionsFromJson(const json& input) {
    keys(input,{"geometry_file","effective_electrodes","crystal_x","crystal_y","electron_parameters_cm","hole_parameters_cm","high_field","reference_density_m3","screening_method"});
    const auto method=ial_json::screeningMethod(input);
    auto options=std::make_shared<IalTransportOptions>();
    options->geometryFile=input.at("geometry_file");
    if (!std::filesystem::path(options->geometryFile).is_absolute())
        throw std::invalid_argument("IALMob geometry_file must be an explicit absolute path");
    options->effectiveElectrodes=input.at("effective_electrodes").get<std::vector<std::string>>();
    options->crystalX=input.at("crystal_x");options->crystalY=input.at("crystal_y");
    options->element.highField=input.value("high_field",true);
    options->element.referenceDensity_m3=input.value("reference_density_m3",1e18);
    for (const auto& [name,p]:input.at("electron_parameters_cm").items())
        options->electrons.emplace(std::stoi(name),IalMobility(ial_json::parameters(p,true),true,method));
    for (const auto& [name,p]:input.at("hole_parameters_cm").items())
        options->holes.emplace(std::stoi(name),IalMobility(ial_json::parameters(p,false),false,method));
    return options;
}

void prepareIalTransportGeometry(MobilityModelConfig& config,const DeviceMesh& mesh) {
    if (config.model!="ialmob") return;
    if (!config.ialmob) throw std::invalid_argument("IALMob requires an explicit parameter/geometry contract");
    if (!config.ialmobGeometry) config.ialmobGeometry=geometry(config,mesh);
}

void updateIalTransportState(MobilityModelConfig& config,const DeviceMesh& mesh,
    const DopingModel& doping,const VectorXd& psi,const VectorXd& n,const VectorXd& p,
    const VectorXd& phin,const VectorXd& phip,const VectorXd& dn,const VectorXd& dp,
    const VectorXd& temperature,const VectorXd& dn_dT,const VectorXd& dp_dT,bool reuseScreening,bool temperatureDerivatives) {
    if (config.model!="ialmob") return;
    ScopedPerformanceTimer timer("ialmob.update");
    if (!config.ialmob) throw std::invalid_argument("IALMob requires an explicit parameter/geometry contract");
    for (const auto* v:{&psi,&n,&p,&phin,&phip})
        if (v->size()!=static_cast<int>(mesh.numNodes())||!v->allFinite()) throw std::invalid_argument("IALMob state size or finiteness mismatch");
    if ((dn.size()!=0&&dn.size()!=n.size())||(dp.size()!=0&&dp.size()!=p.size())) throw std::invalid_argument("IALMob response size mismatch");
    if(temperature.size()!=0){
        for(const auto* v:{&temperature,&dn_dT,&dp_dT})
            if(v->size()!=n.size()||!v->allFinite())throw std::invalid_argument("IALMob thermal state size/finiteness mismatch");
        if((temperature.array()<50.).any())throw std::invalid_argument("IALMob requires temperature >=50 K");
    }else if(dn_dT.size()!=0||dp_dT.size()!=0)throw std::invalid_argument("IALMob thermal responses require temperature");
    const auto same=[](const VectorXd& a,const VectorXd& b) {return a.size()==b.size()&&(a.array()==b.array()).all();};
    const bool spatialDerivatives=temperatureDerivatives || !config.ialmob->element.residualValuesOnly;
    if (config.ialmobState) {
        const auto& s=*config.ialmobState;
        if ((!temperatureDerivatives || s.hasTemperatureDerivatives) && (!spatialDerivatives || s.hasSpatialDerivatives) && same(s.psi,psi)&&same(s.n,n)&&same(s.p,p)&&same(s.phin,phin)&&same(s.phip,phip)&&same(s.dn,dn)&&same(s.dp,dp)&&same(s.temperature,temperature)&&same(s.dn_dT,dn_dT)&&same(s.dp_dT,dp_dT)) {
            incrementPerformanceCounter("ialmob.state_hits");
            return;
        }
    }
    prepareIalTransportGeometry(config,mesh);
    ScopedPerformanceTimer stateTimer("ialmob.state_build");
    auto s=std::make_shared<IalTransportState>();s->psi=psi;s->n=n;s->p=p;s->phin=phin;s->phip=phip;s->dn=dn;s->dp=dp;s->temperature=temperature;s->dn_dT=dn_dT;s->dp_dT=dp_dT;
    s->hasTemperatureDerivatives=temperatureDerivatives;
    s->hasSpatialDerivatives=spatialDerivatives;
    s->electronEdges.resize(mesh.numEdges(),0.);s->holeEdges.resize(mesh.numEdges(),0.);
    auto options=config.ialmob->element;
    IalScreeningCache screening;
    IalMobilityPreparationCache preparation;
    options.screeningCache=reuseScreening?&screening:nullptr;
    options.preparationCache=options.reuseLocalPreparation?&preparation:nullptr;
    options.spatialDerivatives=spatialDerivatives;
    if(temperature.size()!=0){options.temperatureDependentHighField=true;options.temperatureDerivatives=true;}
    if(!temperatureDerivatives)options.temperatureDerivatives=false;
    options.electronField=config.electronField;options.holeField=config.holeField;
    const Real velocityFactor=config.internalMobilityToM2PerVS*config.internalFieldToVPerM;
    options.electronField.saturationVelocity*=velocityFactor;options.holeField.saturationVelocity*=velocityFactor;
    for (const auto& support:config.ialmobGeometry->cells) {
        std::array<IalElementVertexState,3> states;std::array<const IalMobility*,3> em,hm;
        const auto& c=mesh.getCell(support.cellId);
        for (int k=0;k<3;++k) {
            const Index node=c.node_ids[k];const Real factor=config.internalConcentrationToM3;
            states[k]={psi[node],phin[node],phip[node],doping.donors(node)*factor,doping.acceptors(node)*factor,
                n[node]*factor,p[node]*factor,dn.size()?dn[node]*factor:0.,dp.size()?dp[node]*factor:0.,
                temperature.size()?temperature[node]:300.,dn_dT.size()?dn_dT[node]*factor:0.,dp_dT.size()?dp_dT[node]*factor:0.};
            em[k]=&config.ialmob->electrons.at(support.family[k]);hm[k]=&config.ialmob->holes.at(support.family[k]);
        }
        auto r=evaluateIalElementMobility(support.geometry,states,em,hm,options);
        for (auto* v:{&r.electron,&r.hole,&r.electronLowField,&r.holeLowField}) {
            v->value/=config.internalMobilityToM2PerVS;
            for (auto& d:v->derivative) d/=config.internalMobilityToM2PerVS;
        }
        for(auto* v:{&r.electronTemperatureDerivative,&r.holeTemperatureDerivative,
            &r.electronLowTemperatureDerivative,&r.holeLowTemperatureDerivative})
            for(auto& d:*v)d/=config.internalMobilityToM2PerVS;
        s->cells.push_back(r);
    }
    for (Index e=0;e<mesh.numEdges();++e)
        for (const auto& c:config.ialmobGeometry->edges[e]) {
            s->electronEdges[e]+=c.weight*s->cells[c.support].electron.value;
            s->holeEdges[e]+=c.weight*s->cells[c.support].hole.value;
        }
    config.ialmobState=s;
    incrementPerformanceCounter("ialmob.local_preparation_hits",preparation.hits());
    incrementPerformanceCounter("ialmob.local_preparation_builds",preparation.size());
    ialKernelProfile.localPreparationHits+=preparation.hits();
    ialKernelProfile.localPreparationBuilds+=preparation.size();
}
Real ialEdgeMobility(const MobilityModelConfig& config,Index edge,CarrierType carrier,bool lowField) {
    if (!config.ialmobState) throw std::logic_error("IALMob edge requested without a live state");
    if (lowField) {
        Real result=0.;
        for (const auto& c:config.ialmobGeometry->edges.at(edge)) {
            const auto& r=config.ialmobState->cells[c.support];
            result+=c.weight*(carrier==CarrierType::Electron?r.electronLowField.value:r.holeLowField.value);
        }
        return result;
    }
    return carrier==CarrierType::Electron?config.ialmobState->electronEdges.at(edge):config.ialmobState->holeEdges.at(edge);
}
} // namespace vela
