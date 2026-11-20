#include "vela/equation/IalTransport.h"
#include "vela/mesh/ContactCell.h"
#include "vela/physics/IalInterfaceGeometry.h"
#include "vela/physics/IalMobilityJson.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/core/IalKernelProfiling.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <set>
#include <numeric>
#include <sstream>
#include <iomanip>

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
template<class Scalar, class Segment>
auto buildSurfaceGeometry(const DeviceMesh& mesh, const std::vector<std::array<Scalar,2>>& xy,
                          const std::map<Pair,std::vector<Index>>& adjacent,
                          const std::set<Pair>& contactEdges,
                          const std::array<Scalar,3>& crystalX, const std::array<Scalar,3>& crystalY)
{
    std::vector<Segment> segments;std::map<Index,std::array<Scalar,2>> tangents;
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
        const Scalar dx=b[0]-a[0],dy=b[1]-a[1],length=std::hypot(dx,dy);
        if (oxide) {
            std::array<Scalar,2> center{};
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
    return std::pair{buildIalInterfaceGeometry(xy, segments, crystalX, crystalY), std::move(tangents)};
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
    result->cellSupportIndex.assign(mesh.numCells(),-1);
    std::vector<std::array<Real,2>> xy;
    for (const auto& n:mesh.nodes()) xy.push_back({n.x*config.internalLengthToM,n.y*config.internalLengthToM});
    std::map<Pair,std::vector<Index>> adjacent;std::map<Pair,Index> edgeIds;
    for (const auto& e:mesh.edges()) edgeIds[pair(e.n0,e.n1)]=e.id;
    for (const auto& c:mesh.cells()) {
        if (c.node_ids.size()!=3) throw std::invalid_argument("IALMob transport requires a Tri3 mesh");
        for (int k=0;k<3;++k) adjacent[pair(c.node_ids[k],c.node_ids[(k+1)%3])].push_back(c.id);
    }
    std::vector<bool> effectiveNodes(mesh.numNodes(),false);std::set<Pair> contactEdges;std::set<std::string> foundContacts;
    for (const auto& contact:mesh.contacts()) {
        const bool effective=std::find(options.effectiveElectrodes.begin(),options.effectiveElectrodes.end(),contact.name)!=options.effectiveElectrodes.end();
        if (effective) {
            for (Index node:contact.node_ids) effectiveNodes.at(node)=true;
            foundContacts.insert(contact.name);
        }
        for (const auto& edge:contact.edge_node_ids) contactEdges.insert(pair(edge[0],edge[1]));
        if (contact.edge_node_ids.empty()) {
            const std::set<Index> nodes(contact.node_ids.begin(),contact.node_ids.end());
            for (const auto& [edge,cells]:adjacent)
                if (cells.size()==1&&nodes.count(edge.first)&&nodes.count(edge.second)) contactEdges.insert(edge);
        }
    }
    for (const auto& name:options.effectiveElectrodes)
        if (!foundContacts.count(name)) throw std::invalid_argument("Unknown IALMob effective electrode: "+name);
    const auto [nodes,tangents]=buildSurfaceGeometry<Real,IalInterfaceSegment>(mesh,xy,adjacent,contactEdges,options.crystalX,options.crystalY);
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
        }
        support.geometry.touchesEffectiveElectrode=detail::cellTouchesContact(c,effectiveNodes);
        support.geometry.partialBoundaryLayer=tangents.count(id)!=0;
        if (support.geometry.partialBoundaryLayer) support.geometry.boundaryTangent=tangents.at(id);
        const auto coeff=row.at("edge_coefficients").get<std::array<Real,3>>();
        const std::size_t index=result->cells.size();result->cells.push_back(support);
        result->cellSupportIndex[id]=static_cast<int>(index);
        for (int k=0;k<3;++k) {
            if (!std::isfinite(coeff[k])||coeff[k]<0.) throw std::invalid_argument("Invalid IALMob transport coefficient");
            const Index edge=edgeIds.at(pair(c.node_ids[k],c.node_ids[(k+1)%3]));
            coefficientSum[edge]+=coeff[k];result->edges[edge].push_back({index,coeff[k]});
        }
    }
    for (const auto& c:mesh.cells())
        if (silicon(mesh,c)&&!seenCells.count(c.id)) throw std::invalid_argument("Missing IALMob Si cell measure");
    if (!options.avalancheGeometryFile.empty()) {
        std::ifstream sourceStream(options.avalancheGeometryFile);
        if (!sourceStream) throw std::invalid_argument("Cannot read IALMob avalanche geometry");
        const json source=json::parse(sourceStream);
        keys(source,{"schema","scope","node_count","cell_count","cells"});
        const auto indexValue=[](const json& value) -> Index {
            if (!value.is_number_integer() || value.get<long double>()<0 ||
                value.get<long double>()>std::numeric_limits<Index>::max())
                throw std::invalid_argument("IALMob avalanche IDs/counts must be nonnegative integers");
            return value.get<Index>();
        };
        if (source.at("schema")!="vela.ialmob.observed_avalanche_geometry.v1" ||
            indexValue(source.at("node_count"))!=mesh.numNodes() ||
            indexValue(source.at("cell_count"))!=mesh.numCells() ||
            !source.at("cells").is_array())
            throw std::invalid_argument("IALMob avalanche geometry does not match the mesh");
        std::set<Index> sourceSeen;
        for (const auto& row:source.at("cells")) {
            keys(row,{"cell_id","node_ids","coordinates_m","mobility_vertex","source_measure_m2"});
            const Index id=indexValue(row.at("cell_id"));
            for (const auto& key:{"node_ids","coordinates_m","source_measure_m2"})
                if (!row.at(key).is_array() || row.at(key).size()!=3)
                    throw std::invalid_argument("IALMob avalanche requires exactly three vertex entries");
            for (const auto& node:row.at("node_ids")) indexValue(node);
            for (const auto& coordinate:row.at("coordinates_m"))
                if (!coordinate.is_array() || coordinate.size()!=2)
                    throw std::invalid_argument("IALMob avalanche coordinates must be 2-D");
            if (id>=mesh.numCells() || result->cellSupportIndex.at(id)<0 ||
                !sourceSeen.insert(id).second ||
                row.at("node_ids").get<std::vector<Index>>()!=mesh.getCell(id).node_ids)
                throw std::invalid_argument("IALMob avalanche cell topology mismatch");
            auto& support=result->cells.at(result->cellSupportIndex.at(id));
            const auto coordinates=row.at("coordinates_m").get<std::array<std::array<Real,2>,3>>();
            for (int i=0;i<3;++i) for (int d=0;d<2;++d) {
                const Real expected=support.geometry.coordinates_m[i][d],actual=coordinates[i][d];
                if (!std::isfinite(actual) || std::abs(actual-expected)>
                    16.*std::numeric_limits<Real>::epsilon()*std::max(std::abs(expected),Real{1e-30}))
                    throw std::invalid_argument("IALMob avalanche coordinate mismatch");
            }
            IalElementAvalancheOptions selected;
            selected.mobilityVertex=indexValue(row.at("mobility_vertex"));
            selected.sourceMeasure_m2=row.at("source_measure_m2");
            Real total=0.;
            for (Real v:selected.sourceMeasure_m2) {
                if (!std::isfinite(v)||v<0.) throw std::invalid_argument("Invalid IALMob avalanche source share");
                total+=v;
            }
            if (selected.mobilityVertex>=3 || !std::isfinite(total) || total<=0.)
                throw std::invalid_argument("Invalid IALMob avalanche selection or source measure");
            const auto& xy=support.geometry.coordinates_m;
            const Real area=.5*std::abs((xy[1][0]-xy[0][0])*(xy[2][1]-xy[0][1])-
                (xy[1][1]-xy[0][1])*(xy[2][0]-xy[0][0]));
            // Native truncated shares may sum to less than the triangle area.
            // The relative allowance covers coordinate/export roundoff only.
            if (total>area*(1.+1e-10))
                throw std::invalid_argument("IALMob avalanche source shares exceed cell area");
            support.avalanche=selected;
        }
        if (sourceSeen!=seenCells) throw std::invalid_argument("Missing IALMob avalanche Si cell");
    }
    for (const auto& edge:mesh.edges()) {
        const Real expected=edge.transport_couple>=0.?edge.transport_couple:edge.couple;
        const Real actual=coefficientSum[edge.id]*edge.length;
        if (!result->edges[edge.id].empty()&&std::abs(expected-actual)>1e-10*std::max(std::abs(expected),std::abs(actual))) {
            std::ostringstream message;
            message << std::setprecision(17) << "IALMob cell coefficients disagree with transport edge coupling: edge "
                    << edge.id << " (" << edge.n0 << ',' << edge.n1 << "), expected " << expected << ", actual " << actual;
            throw std::invalid_argument(message.str());
        }
        if (coefficientSum[edge.id]>0.)
            for (auto& c:result->edges[edge.id]) c.weight/=coefficientSum[edge.id];
    }
    result->constantElectronEdges.assign(mesh.numEdges(),0.);
    result->constantHoleEdges.assign(mesh.numEdges(),0.);
    std::set<std::string> foundConstantMaterials;
    for (const auto& edge:mesh.edges()) {
        std::array<Real,2> sum{};int count=0;
        for (Index id:adjacent.at(pair(edge.n0,edge.n1))) {
            const auto& name=mesh.getRegion(mesh.getCell(id).region_id).material;
            const auto entry=options.constantMaterialMobility.find(name);
            if (entry==options.constantMaterialMobility.end()) continue;
            foundConstantMaterials.insert(name);++count;
            for(int k=0;k<2;++k) sum[k]+=entry->second[k]*1e-4/config.internalMobilityToM2PerVS;
        }
        if (!count) continue;
        if (!result->edges[edge.id].empty())
            throw std::invalid_argument("IALMob constant material must not share a transport edge with Si");
        result->constantElectronEdges[edge.id]=sum[0]/count;
        result->constantHoleEdges[edge.id]=sum[1]/count;
    }
    for (const auto& [name,value]:options.constantMaterialMobility)
        if (!foundConstantMaterials.count(name))
            throw std::invalid_argument("Unknown IALMob constant material: "+name);
    return result;
}
}

std::shared_ptr<const IalTransportOptions> ialTransportOptionsFromJson(
    const json& input, const json* mobility, UnitScalingConfig scaling)
{
    keys(input,{"screening_method","geometry_file","avalanche_geometry_file","effective_electrodes","crystal_x","crystal_y","electron_parameters_cm","hole_parameters_cm","high_field","reference_density_m3","constant_material_mobility_cm2_per_Vs","decimal_coefficients"});
    const auto method=ial_json::screeningMethod(input);
    auto options=std::make_shared<IalTransportOptions>();
    options->geometryFile=input.at("geometry_file");
    options->avalancheGeometryFile=input.value("avalanche_geometry_file",std::string{});
    if (input.contains("avalanche_geometry_file") &&
        (options->avalancheGeometryFile.empty() || !std::filesystem::path(options->avalancheGeometryFile).is_absolute()))
        throw std::invalid_argument("IALMob avalanche_geometry_file must be an explicit absolute path");
    if (!std::filesystem::path(options->geometryFile).is_absolute())
        throw std::invalid_argument("IALMob geometry_file must be an explicit absolute path");
    options->effectiveElectrodes=input.at("effective_electrodes").get<std::vector<std::string>>();
    options->crystalX=input.at("crystal_x");options->crystalY=input.at("crystal_y");
    options->element.highField=input.value("high_field",true);
    options->element.referenceDensity_m3=input.value("reference_density_m3",1e18);
    if (input.contains("constant_material_mobility_cm2_per_Vs")) {
        const auto& entries=input.at("constant_material_mobility_cm2_per_Vs");
        if (!entries.is_object()) throw std::invalid_argument("IALMob constant materials must be an object");
        for (const auto& [name,row]:entries.items()) {
            if (name.empty()||name=="Si") throw std::invalid_argument("IALMob constant material must be explicit non-Si");
            keys(row,{"electron","hole"});
            const std::array<Real,2> values{row.at("electron").get<Real>(),row.at("hole").get<Real>()};
            for(Real v:values) if (!std::isfinite(v)||v<=0.) throw std::invalid_argument("IALMob constant mobility must be positive finite");
            options->constantMaterialMobility.emplace(name,values);
        }
    }
    const auto* decimals = input.contains("decimal_coefficients") ? &input.at("decimal_coefficients") : nullptr;
    if (decimals) {
        if (!mobility)
            throw std::invalid_argument("IALMob decimal coefficients require the owning mobility configuration");
        keys(*decimals, {"high_field", "electron_parameters_cm", "hole_parameters_cm"});
        const auto& high = decimals->at("high_field");
        keys(high, {"electron_field_beta", "hole_field_beta", "electron_saturation_velocity_m_s",
                    "hole_saturation_velocity_m_s"});
        const auto number = [&](const std::string& key) {
            const auto value = decimalCompatibilityNumber(high.at(key), "IALMob decimal_coefficients." + key,
                                                           mobility->at(key).get<Real>());
            if (!(value > 0))
                throw std::invalid_argument("IALMob decimal high-field coefficients must be positive");
            return value;
        };
        const ArithmeticReal velocityDivisor = scaling.isUnitScaling() ? 100 : 1;
        options->element.decimalHighField = std::array<IalHighFieldDecimal, 2>{{
            {number("electron_saturation_velocity_m_s") / velocityDivisor, number("electron_field_beta")},
            {number("hole_saturation_velocity_m_s") / velocityDivisor, number("hole_field_beta")}}};
    }
    const auto loadFamilies = [&](const char* key, bool electron, auto& models,
                                  auto& extendedModels) {
        const auto& families = input.at(key);
        const auto* decimalFamilies = decimals ? &decimals->at(key) : nullptr;
        if (decimalFamilies && (!decimalFamilies->is_object() || decimalFamilies->size() != families.size()))
            throw std::invalid_argument("IALMob decimal coefficients must cover the compatibility families");
        for (const auto& [name, coefficients] : families.items()) {
            const auto* family = decimalFamilies ? &decimalFamilies->at(name) : nullptr;
            const int orientation = std::stoi(name);
            models.emplace(orientation, IalMobility(ial_json::parameters(coefficients, electron, family),
                                                    electron, method));
            if (family)
                extendedModels.emplace(orientation, ial_json::extendedParameters(coefficients, electron, *family));
        }
    };
    loadFamilies("electron_parameters_cm", true, options->electrons, options->extendedElectrons);
    loadFamilies("hole_parameters_cm", false, options->holes, options->extendedHoles);
    return options;
}

void prepareIalTransportGeometry(MobilityModelConfig& config,const DeviceMesh& mesh) {
    if (config.model!="ialmob") return;
    if (!config.ialmob) throw std::invalid_argument("IALMob requires an explicit parameter/geometry contract");
    if (!config.ialmobGeometry) config.ialmobGeometry=geometry(config,mesh);
}

ExtendedIalTransportMeasures readExtendedIalTransportMeasures(
    const MobilityModelConfig& config, const DeviceMesh& mesh,
    std::optional<long double> lengthScale_m,
    const std::map<std::string, std::array<long double,2>>& constantMobility_cm2_per_Vs)
{
    using WideJson = nlohmann::basic_json<std::map, std::vector, std::string, bool,
                                         std::int64_t, std::uint64_t, long double>;
    if (!config.ialmob || !config.ialmobGeometry)
        throw std::invalid_argument("Extended IALMob measures require validated geometry");
    const auto& geometry = *config.ialmobGeometry;
    const auto read = [](const std::string& file) {
        std::ifstream input(file);
        if (!input)
            throw std::invalid_argument("Cannot read extended IALMob geometry: " + file);
        return WideJson::parse(input);
    };
    const auto matches = [](const auto& exact, const auto& compatibility) {
        for (std::size_t i = 0; i < exact.size(); ++i)
            if (!std::isfinite(exact[i]) || exact[i] < 0 || static_cast<Real>(exact[i]) != compatibility[i])
                throw std::invalid_argument("Extended IALMob measures disagree with validated geometry");
    };
    const auto cellMatches = [&](const auto& row, Index id) {
        return row.at("cell_id") == id &&
            row.at("node_ids").template get<std::vector<Index>>() == mesh.getCell(id).node_ids;
    };

    ExtendedIalTransportMeasures out;
    out.vertex.resize(geometry.cells.size());
    out.edgeWeights.resize(mesh.numEdges());
    std::vector<std::array<long double, 3>> coefficients(geometry.cells.size());
    const auto transport = read(config.ialmob->geometryFile);
    const auto& rows = transport.at("cells");
    if (rows.size() != geometry.cells.size())
        throw std::invalid_argument("Extended IALMob cell count mismatch");
    for (std::size_t i = 0; i < geometry.cells.size(); ++i) {
        const auto& bound = geometry.cells[i];
        const auto& row = rows.at(i);
        if (!cellMatches(row, bound.cellId))
            throw std::invalid_argument("Extended IALMob cell topology mismatch");
        out.vertex[i] = row.at("vertex_measure_m2").get<std::array<long double, 3>>();
        matches(out.vertex[i], bound.geometry.vertexMeasure_m2);
        coefficients[i] = row.at("edge_coefficients").get<std::array<long double, 3>>();
        for (const auto value : coefficients[i])
            if (!std::isfinite(value) || value < 0)
                throw std::invalid_argument("Invalid extended IALMob edge coefficient");
    }
    for (const auto& edge : mesh.edges()) {
        const auto& support = geometry.edges.at(edge.id);
        auto& weights = out.edgeWeights.at(edge.id);
        long double sum = 0;
        for (const auto& contribution : support) {
            const auto& cell = mesh.getCell(geometry.cells.at(contribution.support).cellId);
            for (int k = 0; k < 3; ++k)
                if (pair(cell.node_ids[k], cell.node_ids[(k + 1) % 3]) == pair(edge.n0, edge.n1)) {
                    const auto value = coefficients[contribution.support][k];
                    weights.push_back(value);
                    sum += value;
                    break;
                }
        }
        if (weights.size() != support.size() || !std::isfinite(sum))
            throw std::invalid_argument("Invalid extended IALMob edge support");
        for (std::size_t k = 0; k < weights.size(); ++k) {
            if (sum > 0)
                weights[k] /= sum;
            if (std::abs(weights[k] - support[k].weight) > 8 * std::numeric_limits<Real>::epsilon())
                throw std::invalid_argument("Extended IALMob edge coefficient mismatch");
        }
    }
    if (!config.ialmob->avalancheGeometryFile.empty()) {
        const auto source = read(config.ialmob->avalancheGeometryFile);
        const auto& sourceRows = source.at("cells");
        if (sourceRows.size() != geometry.cells.size())
            throw std::invalid_argument("Extended IALMob avalanche count mismatch");
        out.avalanche.resize(geometry.cells.size());
        std::set<int> seen;
        for (const auto& row : sourceRows) {
            const auto id = row.at("cell_id").get<Index>();
            const auto i = geometry.cellSupportIndex.at(id);
            if (i < 0 || !seen.insert(i).second || !cellMatches(row, id))
                throw std::invalid_argument("Extended IALMob avalanche topology mismatch");
            const auto& bound = geometry.cells.at(i).avalanche;
            if (!bound || row.at("mobility_vertex") != bound->mobilityVertex)
                throw std::invalid_argument("Extended IALMob avalanche mobility vertex mismatch");
            out.avalanche[i] = row.at("source_measure_m2").get<std::array<long double, 3>>();
            matches(out.avalanche[i], bound->sourceMeasure_m2);
        }
    }
    if (lengthScale_m) {
        if (constantMobility_cm2_per_Vs.size() != config.ialmob->constantMaterialMobility.size())
            throw std::invalid_argument("Incomplete extended constant mobility binding");
        std::vector<std::array<long double,2>> xy;
        xy.reserve(mesh.numNodes());
        for (const auto& n : mesh.nodes()) {
            const auto point = mesh.arithmeticCoordinates(n.id);
            xy.push_back({point[0] * *lengthScale_m, point[1] * *lengthScale_m});
        }
        std::map<Pair,std::vector<Index>> adjacent;
        for (const auto& c : mesh.cells())
            for (int k = 0; k < 3; ++k)
                adjacent[pair(c.node_ids[k],c.node_ids[(k+1)%3])].push_back(c.id);
        std::set<Pair> contactEdges;
        for (const auto& electrode : mesh.contacts()) {
            for (const auto& edge : electrode.edge_node_ids)
                contactEdges.insert(pair(edge[0], edge[1]));
            if (electrode.edge_node_ids.empty()) {
                const std::set<Index> nodes(electrode.node_ids.begin(), electrode.node_ids.end());
                for (const auto& [edge, cells] : adjacent)
                    if (cells.size() == 1 && nodes.contains(edge.first) && nodes.contains(edge.second))
                        contactEdges.insert(edge);
            }
        }
        std::array<long double,3> crystalX, crystalY;
        std::ranges::copy(config.ialmob->crystalX, crystalX.begin());
        std::ranges::copy(config.ialmob->crystalY, crystalY.begin());
        const auto [nodes, tangents] = buildSurfaceGeometry<long double,ExtendedIalInterfaceSegment>(
            mesh, xy, adjacent, contactEdges, crystalX, crystalY);
        out.surface.resize(geometry.cells.size());
        for (std::size_t i = 0; i < geometry.cells.size(); ++i) {
            const auto cellId = geometry.cells[i].cellId;
            const auto& c = mesh.getCell(cellId);
            auto& support = out.surface[i];
            if (const auto tangent = tangents.find(cellId); tangent != tangents.end())
                support.tangent = tangent->second;
            for (int k = 0; k < 3; ++k) {
                const auto& node = nodes.at(c.node_ids[k]);
                support.distances[k] = node.distance_m;
                support.family[k] = node.orientationFamily;
            }
        }
        out.constantMobility_m2_per_Vs.resize(mesh.numEdges());
        for (const auto& edge : mesh.edges()) {
            auto& sum = out.constantMobility_m2_per_Vs[edge.id];
            int count = 0;
            for (const auto id : adjacent.at(pair(edge.n0,edge.n1))) {
                const auto& name = mesh.getRegion(mesh.getCell(id).region_id).material;
                const auto entry = constantMobility_cm2_per_Vs.find(name);
                if (entry == constantMobility_cm2_per_Vs.end())
                    continue;
                ++count;
                for (int k = 0; k < 2; ++k)
                    sum[k] += entry->second[k] * 1e-4L;
            }
            if (count)
                for (auto& carrier : sum) carrier /= count;
        }
    }
    return out;
}

IalElementMobilityResult evaluateIalTransportCell(
    const MobilityModelConfig& config,Index cellId,
    const std::array<IalElementVertexState,3>& states,bool spatialDerivatives,
    Real populationThermalVoltage_V) {
    if (!config.ialmob || !config.ialmobGeometry)
        throw std::logic_error("IALMob cell evaluation requires prepared geometry");
    const int index=config.ialmobGeometry->cellSupportIndex.at(cellId);
    if (index<0) return {};
    const auto& support=config.ialmobGeometry->cells.at(index);
    std::array<const IalMobility*,3> em,hm;
    for (int k=0;k<3;++k) {
        em[k]=&config.ialmob->electrons.at(support.family[k]);
        hm[k]=&config.ialmob->holes.at(support.family[k]);
    }
    auto options=config.ialmob->element;
    options.contactField=config.gradQfContactField.value_or(GradQfContactField::ElectricField);
    options.spatialDerivatives=spatialDerivatives;
    options.populationThermalVoltage_V=populationThermalVoltage_V;
    options.temperatureDerivatives=false;
    options.screeningCache=nullptr;options.preparationCache=nullptr;
    options.electronField=config.electronField;options.holeField=config.holeField;
    const Real velocityFactor=config.internalMobilityToM2PerVS*config.internalFieldToVPerM;
    options.electronField.saturationVelocity*=velocityFactor;
    options.holeField.saturationVelocity*=velocityFactor;
    auto r=evaluateIalElementMobility(support.geometry,states,em,hm,options);
    for (auto* value:{&r.electron,&r.hole,&r.electronLowField,&r.holeLowField}) {
        value->value/=config.internalMobilityToM2PerVS;
        for (auto& d:value->derivative) d/=config.internalMobilityToM2PerVS;
    }
    return r;
}

IalElementAvalancheResult evaluateIalTransportAvalancheCell(
    const MobilityModelConfig& config,Index cellId,
    const std::array<IalElementVertexState,3>& states,
    IalElementAvalancheOptions options) {
    if (config.model!="ialmob" || !config.ialmob || !config.ialmobGeometry)
        throw std::invalid_argument("Element density-gradient source requires prepared IALMob geometry");
    const int index=config.ialmobGeometry->cellSupportIndex.at(cellId);
    if (index<0) return {}; // Source is restricted to the original Si material.
    const auto& support=config.ialmobGeometry->cells.at(index);
    if (!support.avalanche)
        throw std::invalid_argument("Element density-gradient source requires observed avalanche geometry");
    if (config.ialmob->element.referenceDensity_m3!=0.)
        throw std::invalid_argument("Element density-gradient source requires zero HFS reference density");
    const Real velocityFactor=config.internalMobilityToM2PerVS*config.internalFieldToVPerM;
    const auto same=[](Real a,Real b){return std::abs(a-b)<=1e-12*std::max(std::abs(a),std::abs(b));};
    if (!same(config.electronField.saturationVelocity*velocityFactor,1.07e5) ||
        !same(config.holeField.saturationVelocity*velocityFactor,8.37e4) ||
        config.electronField.beta!=1.109 || config.holeField.beta!=1.213)
        throw std::invalid_argument("Element density-gradient source requires identified 300 K HFS parameters");
    options.mobilityVertex=support.avalanche->mobilityVertex;
    options.sourceMeasure_m2=support.avalanche->sourceMeasure_m2;
    options.highFieldMobility=config.ialmob->element.highField;
    options.mobilityContactField=config.gradQfContactField.value_or(GradQfContactField::ElectricField);
    std::array<const IalMobility*,3> em,hm;
    for (int k=0;k<3;++k) {
        em[k]=&config.ialmob->electrons.at(support.family[k]);
        hm[k]=&config.ialmob->holes.at(support.family[k]);
    }
    return evaluateIalElementAvalanche(support.geometry,states,em,hm,options);
}

ExtendedIalTransportCell evaluateExtendedIalTransportCell(
    const MobilityModelConfig& config, const DeviceMesh& mesh, const DopingModel& doping,
    Index cellId, const std::array<IalPopulationVertex,3>& state,
    ArithmeticReal lengthScale_m, ArithmeticReal dopingScale_m3,
    ArithmeticReal potentialScale_V, ArithmeticReal thermalVoltage_V,
    ArithmeticReal temperature_K, std::optional<IalElementAvalancheOptions> avalanche,
    const ExtendedIalTransportMeasures* measures)
{
    if (!config.ialmob || !config.ialmobGeometry || doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument("Extended IALMob requires bound mesh, doping, and geometry");
    for (auto scale : {lengthScale_m, dopingScale_m3})
        if (!(scale > 0) || !std::isfinite(scale))
            throw std::invalid_argument("Extended IALMob requires finite positive coordinate scales");
    const auto supportIndex = config.ialmobGeometry->cellSupportIndex.at(cellId);
    if (supportIndex < 0)
        throw std::invalid_argument("Extended IALMob cell is outside the silicon transport support");
    const auto& support = config.ialmobGeometry->cells.at(supportIndex);
    const auto& options = *config.ialmob;
    if (options.element.highField && !options.element.decimalHighField)
        throw std::invalid_argument("Extended IALMob high-field coefficients are unavailable");
    if (options.element.generatedLowFieldPartials || options.element.temperatureDependentHighField)
        throw std::invalid_argument("Extended IALMob cell uses the isothermal reference kernel");
    const auto fields = detail::buildIalPopulationFields(state, potentialScale_V, thermalVoltage_V);
    BasicIalElementGeometry<ArithmeticReal> geometry;
    geometry.partialBoundaryLayer = support.geometry.partialBoundaryLayer;
    geometry.touchesEffectiveElectrode = support.geometry.touchesEffectiveElectrode;
    const auto* surface = measures && !measures->surface.empty() ? &measures->surface.at(supportIndex) : nullptr;
    if (surface)
        std::ranges::copy(surface->tangent, geometry.boundaryTangent.begin());
    else
        std::ranges::copy(support.geometry.boundaryTangent, geometry.boundaryTangent.begin());
    const auto& cell = mesh.getCell(cellId);
    std::array<std::array<ArithmeticReal,2>,3> doping_m3;
    using Parameters = BasicIalMobilityParameters<long double>;
    std::array<const Parameters*,3> electrons, holes;
    for (std::size_t i = 0; i < cell.node_ids.size(); ++i) {
        const auto node = cell.node_ids[i];
        auto point = mesh.arithmeticCoordinates(node);
        for (auto& coordinate : point) coordinate *= lengthScale_m;
        geometry.coordinates_m[i] = point;
        geometry.interfaceDistance_m[i] = surface ? surface->distances[i] : support.geometry.interfaceDistance_m[i];
        geometry.vertexMeasure_m2[i] = measures ? measures->vertex.at(supportIndex)[i]
                                               : support.geometry.vertexMeasure_m2[i];
        auto operands = doping.dopingOperandsForAssembly(node);
        for (auto& concentration : operands) concentration *= dopingScale_m3;
        doping_m3[i] = operands;
        const auto family = surface ? surface->family[i] : support.family[i];
        electrons[i] = &options.extendedElectrons.at(family);
        holes[i] = &options.extendedHoles.at(family);
    }
    const auto contactField = config.gradQfContactField.value_or(GradQfContactField::ElectricField);
    const auto screeningMethod = options.electrons.at(support.family.front()).screeningMethod();
    const auto highField = options.element.highField ? options.element.decimalHighField : std::nullopt;
    const auto referenceDensity = static_cast<ArithmeticReal>(options.element.referenceDensity_m3);
    const auto mobility = evaluateIalPopulationMobility(geometry, fields, doping_m3, electrons, holes,
        temperature_K, referenceDensity, highField, contactField, screeningMethod);
    const auto measure = std::accumulate(geometry.vertexMeasure_m2.begin(),
                                         geometry.vertexMeasure_m2.end(), ArithmeticReal{0});
    if (!(measure > 0) || !std::isfinite(measure))
        throw std::invalid_argument("Extended IALMob requires positive vertex measure");
    ExtendedIalTransportCell out;
    for (std::size_t i = 0; i < mobility.electron.size(); ++i) {
        const detail::IalPopulationDual weight(geometry.vertexMeasure_m2[i]/measure);
        out.mobility_m2_per_Vs[0] = out.mobility_m2_per_Vs[0] + weight*mobility.electron[i];
        out.mobility_m2_per_Vs[1] = out.mobility_m2_per_Vs[1] + weight*mobility.hole[i];
    }
    if (avalanche) {
        if (!support.avalanche || referenceDensity != 0 || temperature_K != 300.L)
            throw std::invalid_argument("Extended SLOT avalanche requires bound shares, zero HFS reference density and 300 K");
        const auto chosen = support.avalanche->mobilityVertex;
        avalanche->sourceMeasure_m2 = support.avalanche->sourceMeasure_m2;
        out.avalanche = evaluateIalPopulationAvalancheFromMobility(geometry, state,
            potentialScale_V, thermalVoltage_V, {mobility.electron.at(chosen), mobility.hole.at(chosen)}, *avalanche,
            measures ? &measures->avalanche.at(supportIndex) : nullptr,
            measures && measures->ionization ? &*measures->ionization : nullptr);
    }
    return out;
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
    s->electronEdges=config.ialmobGeometry->constantElectronEdges;
    s->holeEdges=config.ialmobGeometry->constantHoleEdges;
    auto options=config.ialmob->element;
    options.contactField=config.gradQfContactField.value_or(GradQfContactField::ElectricField);
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
        Real result=carrier==CarrierType::Electron?config.ialmobGeometry->constantElectronEdges.at(edge):config.ialmobGeometry->constantHoleEdges.at(edge);
        for (const auto& c:config.ialmobGeometry->edges.at(edge)) {
            const auto& r=config.ialmobState->cells[c.support];
            result+=c.weight*(carrier==CarrierType::Electron?r.electronLowField.value:r.holeLowField.value);
        }
        return result;
    }
    return carrier==CarrierType::Electron?config.ialmobState->electronEdges.at(edge):config.ialmobState->holeEdges.at(edge);
}
} // namespace vela
