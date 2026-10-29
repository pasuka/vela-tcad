#include "vela/equation/BVProcessProbe.h"

#include "vela/core/PhysicalConstants.h"
#include "vela/core/PerformanceProfiler.h"
#include "vela/equation/AssemblerUtils.h"
#include "vela/physics/IalMobilityJson.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace vela {
namespace {

constexpr Index InvalidIndex = std::numeric_limits<Index>::max();

std::string fnv1a64(std::string_view text)
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char value : text) {
        hash ^= value;
        hash *= 1099511628211ULL;
    }
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(16) << hash;
    return output.str();
}

Real safeRatio(Real numerator, Real denominator)
{
    return denominator != 0.0 ? numerator / denominator : 0.0;
}

Point2 edgeUnit(const DeviceMesh& mesh, Index node0, Index node1)
{
    const Node& first = mesh.getNode(node0);
    const Node& second = mesh.getNode(node1);
    const Point2 delta{second.x - first.x, second.y - first.y};
    const Real length = delta.norm();
    if (length > 0.0)
        return delta / length;
    return Point2::Zero();
}

Point2 cellGradient(
    const DeviceMesh& mesh,
    Index cellId,
    const VectorXd& values)
{
    bool valid = false;
    Real doubleArea = 0.0;
    const Point2 gradient = detail::cellScalarGradient(
        mesh,
        mesh.getCell(cellId),
        [&](Index node) { return values(static_cast<Eigen::Index>(node)); },
        valid,
        doubleArea);
    return valid ? gradient : Point2::Zero();
}

bool edgeTouchesContact(
    const std::vector<bool>& contactNodes,
    Index node0,
    Index node1)
{
    return contactNodes.at(static_cast<std::size_t>(node0)) ||
        contactNodes.at(static_cast<std::size_t>(node1));
}

void assignEndpointState(
    BVProcessProbeRecord& record,
    const DDSolution& state)
{
    const Eigen::Index i = static_cast<Eigen::Index>(record.node0);
    const Eigen::Index j = static_cast<Eigen::Index>(record.node1);
    record.psi0 = state.psi(i);
    record.psi1 = state.psi(j);
    if (record.carrier == "electron") {
        record.quasiFermi0 = state.phin(i);
        record.quasiFermi1 = state.phin(j);
        record.density0 = state.n(i);
        record.density1 = state.n(j);
    } else {
        record.quasiFermi0 = state.phip(i);
        record.quasiFermi1 = state.phip(j);
        record.density0 = state.p(i);
        record.density1 = state.p(j);
    }
}

void finalizeRecord(
    BVProcessProbeRecord& record,
    const ImpactIonizationModelConfig& impact,
    const std::string& configurationFingerprint)
{
    record.solverCoupled =
        impact.model != "none" && impact.couplingMode == "self_consistent";
    record.zeroMeasure = record.sourceMeasure <= 0.0;
    record.zeroMobility = record.finalMobility <= 0.0;
    record.zeroAlpha = record.alpha <= 0.0;
    record.mobilityLimiter =
        safeRatio(record.finalMobility, record.lowFieldMobility);
    record.generationRate =
        safeRatio(record.sourceIntegral, record.sourceMeasure);
    record.qGContribution = constants::q * record.sourceIntegral;
    record.configurationFingerprint = configurationFingerprint;

    std::ostringstream branches;
    branches << "support=" << record.supportKind
             << ";carrier=" << record.carrier
             << ";coupling=" << impact.couplingMode
             << ";contact=" << (record.contactAdjacent ? 1 : 0)
             << ";zero_measure=" << (record.zeroMeasure ? 1 : 0)
             << ";zero_mobility=" << (record.zeroMobility ? 1 : 0)
             << ";zero_alpha=" << (record.zeroAlpha ? 1 : 0)
             << ";reconstructed_current="
             << (record.reconstructedCurrent ? 1 : 0)
             << ";directional_partition="
             << (record.directionalPartition ? 1 : 0)
             << ";source_mapping=" << impact.sourceMappingMode
             << ";source_volume=" << impact.sourceVolumePolicy
             << ";qf_discretization="
             << impact.quasiFermiGradientDiscretization
             << ";contact_fallback_mode="
             << impact.contactElectricFieldFallbackMode
             << ";qf_carrier_floor="
             << (impact.quasiFermiCarrierTruncation > 0.0 ? 1 : 0)
             << ";minimum_field_cutoff="
             << (impact.minimumField > 0.0 ? 1 : 0)
             << ";refdens_interpolation="
             << ((impact.electronDrivingForceRefDensity > 0.0 ||
                  impact.holeDrivingForceRefDensity > 0.0)
                     ? 1
                     : 0)
             << ";flux_sign="
             << ((record.supportKind=="observed_ialmob_element_vertex" ||
                  record.supportKind=="phumob_element_vertex") ? "not_applicable" : record.directedSgFlux > 0.0
                     ? "positive"
                     : (record.directedSgFlux < 0.0 ? "negative" : "zero"));
    if (!record.productionBranchDetails.empty())
        branches << ';' << record.productionBranchDetails;
    record.activeBranches = branches.str();
    record.activeBranchFingerprint =
        fnv1a64(configurationFingerprint + ";" + record.activeBranches);

    for (std::size_t i = 0; i < record.scatterCount; ++i) {
        const Real source = record.sourceWeights[i] * record.sourceIntegral;
        const Real residual = record.solverCoupled ? -source : 0.0;
        record.electronResidualContributions[i] = residual;
        record.holeResidualContributions[i] = residual;
    }
}

void appendRecord(
    BVProcessProbeResult& result,
    BVProcessProbeRecord record,
    const ImpactIonizationModelConfig& impact)
{
    finalizeRecord(record, impact, result.configurationFingerprint);
    result.totalSourceIntegral += record.sourceIntegral;
    result.totalQGContribution += record.qGContribution;
    for (std::size_t i = 0; i < record.scatterCount; ++i) {
        result.electronResidualContribution +=
            record.electronResidualContributions[i];
        result.holeResidualContribution +=
            record.holeResidualContributions[i];
    }
    result.records.push_back(std::move(record));
}

void assignScatterFromNodeIntegrals(
    BVProcessProbeRecord& record,
    const std::vector<Real>& nodeIntegrals)
{
    for (Index node = 0; node < static_cast<Index>(nodeIntegrals.size()); ++node) {
        const Real contribution =
            nodeIntegrals.at(static_cast<std::size_t>(node));
        if (contribution == 0.0)
            continue;
        if (record.scatterCount >= record.scatterNodes.size()) {
            throw std::runtime_error(
                "BV process probe scatter support exceeds Tri3 capacity");
        }
        const std::size_t slot = record.scatterCount++;
        record.scatterNodes[slot] = node;
        record.sourceWeights[slot] =
            safeRatio(contribution, record.sourceIntegral);
    }
}

void requireState(const DDSolution& state, Index nodeCount)
{
    const auto valid = [nodeCount](const VectorXd& values) {
        return values.size() == static_cast<Eigen::Index>(nodeCount) &&
            values.allFinite();
    };
    if (!valid(state.psi) || !valid(state.phin) || !valid(state.phip) ||
        !valid(state.n) || !valid(state.p)) {
        throw std::invalid_argument(
            "BV process probe requires finite state fields matching the mesh");
    }
}

} // namespace

std::string bvProcessConfigurationFingerprint(
    const MobilityModelConfig& mobility,
    const ImpactIonizationModelConfig& impact)
{
    std::ostringstream canonical;
    canonical << std::setprecision(17)
              << "mobility.model=" << mobility.model
              << ";mobility.drive=" << mobility.highFieldDrivingForce
              << ";mobility.doping_basis=" << mobility.dopingConcentrationBasis
              << ";mobility.jacobian_field_derivatives="
              << (mobility.jacobianFieldDerivatives ? 1 : 0)
              << ";mobility.electron_ct=" << mobility.electronCT.muMin << ':'
              << mobility.electronCT.nRef << ':' << mobility.electronCT.alpha
              << ";mobility.hole_ct=" << mobility.holeCT.muMin << ':'
              << mobility.holeCT.nRef << ':' << mobility.holeCT.alpha
              << ";mobility.electron_masetti="
              << mobility.electronMasetti.muConst << ':'
              << mobility.electronMasetti.muMin1 << ':'
              << mobility.electronMasetti.muMin2 << ':'
              << mobility.electronMasetti.mu1 << ':'
              << mobility.electronMasetti.pc << ':'
              << mobility.electronMasetti.cr << ':'
              << mobility.electronMasetti.cs << ':'
              << mobility.electronMasetti.alpha << ':'
              << mobility.electronMasetti.beta
              << ";mobility.hole_masetti="
              << mobility.holeMasetti.muConst << ':'
              << mobility.holeMasetti.muMin1 << ':'
              << mobility.holeMasetti.muMin2 << ':'
              << mobility.holeMasetti.mu1 << ':'
              << mobility.holeMasetti.pc << ':'
              << mobility.holeMasetti.cr << ':'
              << mobility.holeMasetti.cs << ':'
              << mobility.holeMasetti.alpha << ':'
              << mobility.holeMasetti.beta
              << ";mobility.electron_field="
              << mobility.electronField.saturationVelocity << ':'
              << mobility.electronField.beta
              << ";mobility.hole_field="
              << mobility.holeField.saturationVelocity << ':'
              << mobility.holeField.beta
              << ";mobility.surface=" << mobility.surface.thetaElectron << ':'
              << mobility.surface.thetaHole << ':' << mobility.surface.beta
              << ':' << mobility.surface.referenceField << ':'
              << mobility.surface.minFactor << ':'
              << mobility.surface.maxFactor << ':'
              << mobility.surface.surfaceRegion
              << ";impact.model=" << impact.model
              << ";impact.coupling=" << impact.couplingMode
              << ";impact.parameter_set=" << impact.parameterSet
              << ";impact.drive=" << impact.drivingForce
              << ";impact.generation=" << impact.generation
              << ";impact.current=" << impact.currentApproximation
              << ";impact.current_magnitude=" << impact.currentMagnitudeMode
              << ";impact.qf_discretization="
              << impact.quasiFermiGradientDiscretization
              << ";impact.source_mapping=" << impact.sourceMappingMode
              << ";impact.source_partition=" << impact.edgeSourcePartition
              << ";impact.source_volume=" << impact.sourceVolumePolicy
              << ";impact.midpoint_density="
              << impact.cellReconstructedMidpointDensity
              << ";impact.drive_interpolation="
              << impact.drivingForceInterpolation
              << ";impact.source_geometry_scale=" << impact.sourceGeometryScale
              << ";impact.source_volume_factor=" << impact.sourceVolumeFactor
              << ";impact.contact_fallback="
              << (impact.contactElectricFieldFallback ? 1 : 0)
              << ";impact.ialmob_boundary_projection="
              << (impact.ialmobBoundaryFieldProjection ? 1 : 0)
              << ";impact.contact_fallback_scope="
              << impact.contactElectricFieldFallbackScope
              << ";impact.contact_fallback_mode="
              << impact.contactElectricFieldFallbackMode
              << ";impact.refdens_n=" << impact.electronDrivingForceRefDensity
              << ";impact.refdens_p=" << impact.holeDrivingForceRefDensity
              << ";impact.qf_truncation=" << impact.quasiFermiCarrierTruncation
              << ";impact.minimum_field=" << impact.minimumField
              << ";impact.debug_raw="
              << (impact.debugRawVanOverstraeten ? 1 : 0)
              << ";impact.a_scale=" << impact.aScale
              << ";impact.b_scale=" << impact.bScale
              << ";impact.selberherr=" << impact.electronA << ':'
              << impact.electronB << ':' << impact.holeA << ':'
              << impact.holeB << ':' << impact.carrierVelocity
              << ";impact.van_overstraeten=" << impact.electronALow << ':'
              << impact.electronAHigh << ':' << impact.electronBLow << ':'
              << impact.electronBHigh << ':' << impact.holeALow << ':'
              << impact.holeAHigh << ':' << impact.holeBLow << ':'
              << impact.holeBHigh << ':' << impact.switchField << ':'
              << impact.phononEnergy << ':' << impact.referenceTemperature_K
              << ':' << impact.temperature_K;
    for (const std::string& selector : mobility.surface.surfaceInterface)
        canonical << ";mobility.surface_interface=" << selector;
    if (mobility.model=="ialmob" && mobility.ialmob) {
        const auto& options=*mobility.ialmob;
        canonical << ";ialmob.high_field=" << options.element.highField
                  << ";ialmob.refdens=" << options.element.referenceDensity_m3;
        for (const auto* models:{&options.electrons,&options.holes}) {
            canonical << ";carrier_parameters=";
            for (const auto& [family,model]:*models) {
                canonical << family << ':';
                for (const auto& [name,member]:ial_json::fields)
                    canonical << name << '=' << model.parameters().*member << ',';
            }
        }
        // Hash the prepared values, not file names: replacing an input file at
        // the same path must not make distinct source configurations identical.
        if (mobility.ialmobGeometry) for (const auto& cell:mobility.ialmobGeometry->cells) {
            canonical << ";ialmob.cell=" << cell.cellId << ':';
            for (int i=0;i<3;++i) {
                for (Real coordinate:cell.geometry.coordinates_m[i]) canonical << coordinate << ',';
                canonical << cell.geometry.interfaceDistance_m[i] << ',' << cell.family[i] << ','
                          << cell.geometry.vertexMeasure_m2[i] << ';';
            }
            canonical << cell.geometry.touchesEffectiveElectrode << ':'
                      << cell.geometry.partialBoundaryLayer << ':'
                      << cell.geometry.boundaryTangent[0] << ':' << cell.geometry.boundaryTangent[1];
            if (cell.avalanche) {
                canonical << ";source_vertex=" << cell.avalanche->mobilityVertex;
                for (Real measure:cell.avalanche->sourceMeasure_m2) canonical << ':' << measure;
            }
        }
    }
    if (detail::usesPhuMobElementAvalancheSource(impact)) {
        canonical << ";element_box=" << mobility.edgeAveraging
                  << ";surface_discretization=" << mobility.surface.discretization
                  << ";hfs_discretization=" << mobility.highFieldGradientDiscretization;
        for(const auto& material:mobility.constantMaterials)
            canonical << ";constant_material=" << material;
    }
    return fnv1a64(canonical.str());
}

BVProcessProbeResult evaluateBVProcessProbe(
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const DDSolution& state,
    const MobilityModelConfig& inputMobilityConfig,
    const ImpactIonizationModelConfig& impactConfig,
    const BandgapNarrowingConfig& bandgapNarrowing,
    const MaterialDatabase& materials,
    Real temperature_K,
    Real fieldFactor)
{
    ScopedPerformanceTimer timer("bv.process_probe");
    incrementPerformanceCounter("bv.process_probe_calls");
    requireState(state, mesh.numNodes());
    if (doping.numNodes() != mesh.numNodes())
        throw std::invalid_argument("BV process probe doping size mismatch");

    auto mobilityConfig=inputMobilityConfig;
    if (detail::usesIalElementAvalancheSource(impactConfig) && mobilityConfig.model!="ialmob")
        throw std::invalid_argument("Observed element source requires IALMob transport");
    if (mobilityConfig.model=="ialmob" && impactConfig.model!="none") {
        if (!detail::usesElementLocalAvalancheSource(impactConfig) || std::abs(temperature_K-300.)>1e-10)
            throw std::invalid_argument("IALMob BV probe requires a supported element source at 300 K");
        // The standalone probe and branch fingerprint can precede any coupled
        // residual. Cell mobility is evaluated from the supplied state below.
        prepareIalTransportGeometry(mobilityConfig,mesh);
    }

    BVProcessProbeResult result;
    result.configurationFingerprint =
        bvProcessConfigurationFingerprint(mobilityConfig, impactConfig);
    if (impactConfig.model == "none")
        return result;
    if (!detail::usesEdgeCurrentAvalancheSource(impactConfig)) {
        throw std::invalid_argument(
            "BV process probe currently requires current-density avalanche");
    }

    const Real thermalVoltage =
        constants::kb * temperature_K / constants::q;
    const std::vector<Material> cellMaterials =
        detail::buildCellMaterials(mesh, materials, temperature_K);
    const std::vector<Real> ni = detail::buildValidatedEffectiveNodeNi(
        "BVProcessProbe",
        mesh,
        materials,
        doping,
        bandgapNarrowing,
        thermalVoltage);
    const auto mobility = makeMobilityModel(mobilityConfig);
    const auto impact = makeImpactIonizationModel(impactConfig);
    const auto edgeCells = detail::buildEdgeCellMap(mesh);
    const std::vector<bool> contactNodes = detail::contactNodeMask(mesh);

    if (detail::usesPhuMobElementAvalancheSource(impactConfig)) {
        // A standalone diagnostic receives configuration, not the assembler's
        // populated geometry cache. Build the same material/contact support
        // before evaluating HFS and recording contact adjacency.
        detail::updateElementHighFieldGeometry(mobilityConfig, mesh, edgeCells, cellMaterials);
        detail::updateSurfaceMobilityCellGeometry(mobilityConfig,mesh,edgeCells,state.psi,fieldFactor,&cellMaterials);
        for(const auto& cell:mesh.cells()) {
            if(mesh.getRegion(cell.region_id).material!="Si") continue;
            std::array<Real,3> psi,fn,fp,n,p;
            for(int k=0;k<3;++k) {
                const auto node=cell.node_ids[k];psi[k]=state.psi[node];fn[k]=state.phin[node];fp[k]=state.phip[node];
                n[k]=state.n[node];p[k]=state.p[node];
            }
            detail::PhuMobElementAvalancheDiagnostics diagnostic;
            const auto evaluated=detail::phumobElementAvalancheSourceIntegralsLocal(
                impactConfig,mobilityConfig,mesh,doping,cellMaterials,cell.id,psi,fn,fp,n,p,&diagnostic);
            const auto measures=detail::cellBoxNodeMeasures(mesh,cell.id);
            for(int carrier=0;carrier<2;++carrier) for(int k=0;k<3;++k) {
                const bool electron=carrier==0;
                BVProcessProbeRecord record;
                record.supportKind="phumob_element_vertex";record.carrier=electron?"electron":"hole";
                record.cellId=cell.id;record.edgeId=InvalidIndex;record.node0=record.node1=cell.node_ids[k];
                assignEndpointState(record,state);
                record.midpointDensity=electron?n[k]:p[k];
                record.electricFieldVector=-cellGradient(mesh,cell.id,state.psi)*fieldFactor;
                record.quasiFermiGradientVector=cellGradient(mesh,cell.id,electron?state.phin:state.phip)*fieldFactor;
                record.finalMobility=electron?diagnostic.electronMobility[k]:diagnostic.holeMobility[k];
                record.lowFieldMobility=electron?diagnostic.electronLowFieldMobility[k]:diagnostic.holeLowFieldMobility[k];
                record.impactField=electron?diagnostic.electronDrive:diagnostic.holeDrive;
                record.highFieldDrive=electron?diagnostic.electronHighFieldDrive:diagnostic.holeHighFieldDrive;
                record.alpha=electron?diagnostic.electronAlpha:diagnostic.holeAlpha;
                record.currentVector=-record.midpointDensity*record.finalMobility*record.quasiFermiGradientVector;
                record.currentVectorProvenance="vertex_density_phumob_gradqf";
                if(impactConfig.currentApproximation=="phumob_element_sg_best_vertex") {
                    const auto& j=electron?diagnostic.electronSgCurrent:diagnostic.holeSgCurrent;
                    record.currentVector=Point2(j[0],j[1]);
                    record.currentVectorProvenance="sg_best_vertex_effective_band";
                    record.productionBranchDetails="source_vertex="+std::to_string(cell.node_ids[diagnostic.sourceVertex]);
                }
                record.selectedFluxMagnitude=record.currentVector.norm();
                record.sourceMeasure=measures[k];
                record.sourceIntegral=electron?evaluated.electron[k]:evaluated.hole[k];
                record.scatterCount=1;record.scatterNodes[0]=record.node0;record.sourceWeights[0]=1.;
                record.contactAdjacent=mobilityConfig.highFieldCells.at(cell.id).contact;
                record.reconstructedCurrent=true;
                record.productionBranchDetails+=";material=Si;source=phumob_vertex;hfs=element_vertex_partial_layer";
                appendRecord(result,std::move(record),impactConfig);
            }
        }
        return result;
    }

    if (detail::usesIalElementAvalancheSource(impactConfig)) {
        for (const auto& support:mobilityConfig.ialmobGeometry->cells) {
            const auto& cell=mesh.getCell(support.cellId);
            std::array<IalElementVertexState,3> local;
            const Real concentration=mobilityConfig.internalConcentrationToM3;
            const Real length=mobilityConfig.internalLengthToM;
            const Real measure=concentration*length*length;
            for (int i=0;i<3;++i) {
                const auto n=cell.node_ids[i];
                local[i]={state.psi[n],state.phin[n],state.phip[n],doping.donors(n)*concentration,
                    doping.acceptors(n)*concentration,state.n[n]*concentration,state.p[n]*concentration};
            }
            const auto evaluated=evaluateIalTransportAvalancheCell(mobilityConfig,cell.id,local,
                detail::ialElementAvalancheOptionsFromConfig(impactConfig));
            // A one-hot evaluation supplies the matching low-field mobility.
            auto geometry=support.geometry;
            geometry.vertexMeasure_m2={0.,0.,0.};
            geometry.vertexMeasure_m2[support.avalanche->mobilityVertex]=1.;
            std::array<const IalMobility*,3> em,hm;
            for(int i=0;i<3;++i) {
                em[i]=&mobilityConfig.ialmob->electrons.at(support.family[i]);
                hm[i]=&mobilityConfig.ialmob->holes.at(support.family[i]);
            }
            IalElementMobilityOptions lowOptions;lowOptions.highField=false;lowOptions.spatialDerivatives=false;
            const auto low=evaluateIalElementMobility(geometry,local,em,hm,lowOptions);
            for (bool electron:{true,false}) {
                const auto& value=electron?evaluated.electron:evaluated.hole;
                Point2 gradient=cellGradient(mesh,cell.id,electron?state.phin:state.phip)*fieldFactor;
                if (geometry.partialBoundaryLayer&&!geometry.touchesEffectiveElectrode) {
                    const Point2 tangent(geometry.boundaryTangent[0],geometry.boundaryTangent[1]);
                    gradient=tangent*gradient.dot(tangent);
                }
                const auto density=[&](int i){return electron?local[i].electrons_m3:local[i].holes_m3;};
                int minimum=0,maximum=0;
                for(int i=1;i<3;++i){if(density(i)<density(minimum))minimum=i;if(density(i)>density(maximum))maximum=i;}
                for (int i=0;i<3;++i) {
                    BVProcessProbeRecord record;
                    record.supportKind="observed_ialmob_element_vertex";
                    record.carrier=electron?"electron":"hole";
                    record.cellId=cell.id;record.edgeId=InvalidIndex;
                    record.node0=record.node1=cell.node_ids[i];assignEndpointState(record,state);
                    record.midpointDensity=value.density_m3.value/concentration;
                    record.electricFieldVector=-cellGradient(mesh,cell.id,state.psi)*fieldFactor;
                    record.quasiFermiGradientVector=gradient;
                    record.lowFieldMobility=(electron?low.electronLowField.value:low.holeLowField.value)/mobilityConfig.internalMobilityToM2PerVS;
                    record.finalMobility=value.mobility_m2_per_Vs.value/mobilityConfig.internalMobilityToM2PerVS;
                    record.highFieldDrive=value.ionizationField_V_per_m.value/mobilityConfig.internalFieldToVPerM;
                    record.currentVector=-record.midpointDensity*record.finalMobility*gradient;
                    record.currentVectorProvenance="element_extrema_density_gradqf";
                    record.selectedFluxMagnitude=record.currentVector.norm();
                    record.impactField=record.highFieldDrive;
                    record.alpha=impactConfig.unitSystem.mInvToInternalInverseLength(value.alpha_per_m.value);
                    record.sourceMeasure=support.avalanche->sourceMeasure_m2[i]/(length*length);
                    record.sourceIntegral=value.source_per_m_s[i].value/measure;
                    record.scatterCount=1;record.scatterNodes[0]=record.node0;record.sourceWeights[0]=1.;
                    record.contactAdjacent=geometry.touchesEffectiveElectrode;record.reconstructedCurrent=true;
                    record.productionBranchDetails="density_min="+std::to_string(minimum)+";density_max="+std::to_string(maximum)+
                        ";mobility_vertex="+std::to_string(support.avalanche->mobilityVertex)+
                        ";tangent="+std::to_string(geometry.partialBoundaryLayer&&!geometry.touchesEffectiveElectrode);
                    appendRecord(result,std::move(record),impactConfig);
                }
            }
        }
        return result;
    }

    if (detail::usesElementEdgeGssLauxAvalancheSource(impactConfig)) {
        const auto production =
            detail::elementEdgeGssLauxAvalancheSourceRecords(
                impactConfig,
                *impact,
                mobilityConfig,
                *mobility,
                edgeCells,
                mesh,
                doping,
                cellMaterials,
                state.psi,
                state.phin,
                state.phip,
                state.n,
                state.p,
                ni,
                thermalVoltage,
                fieldFactor);
        for (const auto& cellRecord : production) {
            const Cell& cell = mesh.getCell(cellRecord.cellId);
            const Point2 electricVector =
                -cellGradient(mesh, cellRecord.cellId, state.psi) * fieldFactor;
            const Point2 electronQfGradient =
                cellGradient(mesh, cellRecord.cellId, state.phin) * fieldFactor;
            const Point2 holeQfGradient =
                cellGradient(mesh, cellRecord.cellId, state.phip) * fieldFactor;
            for (const std::string carrier : {"electron", "hole"}) {
                for (std::size_t local = 0; local < 3; ++local) {
                    const Index node0 = cell.node_ids[local];
                    const Index node1 =
                        cell.node_ids[(local + 1) % 3];
                    BVProcessProbeRecord edgeRecord;
                    edgeRecord.supportKind = "element_edge_gss_laux";
                    edgeRecord.carrier = carrier;
                    edgeRecord.cellId = cellRecord.cellId;
                    edgeRecord.localEdge = static_cast<int>(local);
                    edgeRecord.edgeId = cellRecord.edgeIds[local];
                    edgeRecord.node0 = node0;
                    edgeRecord.node1 = node1;
                    assignEndpointState(edgeRecord, state);
                    edgeRecord.midpointDensity =
                        0.5 * (edgeRecord.density0 + edgeRecord.density1);
                    edgeRecord.electricFieldVector = electricVector;
                    edgeRecord.quasiFermiGradientVector =
                        carrier == "electron"
                        ? electronQfGradient
                        : holeQfGradient;
                    edgeRecord.lowFieldMobility = carrier == "electron"
                        ? cellRecord.electronLowFieldMobilities[local]
                        : cellRecord.holeLowFieldMobilities[local];
                    edgeRecord.finalMobility = carrier == "electron"
                        ? cellRecord.electronMobilities[local]
                        : cellRecord.holeMobilities[local];
                    edgeRecord.highFieldDrive = carrier == "electron"
                        ? cellRecord.electronMobilityDrivingFields[local]
                        : cellRecord.holeMobilityDrivingFields[local];
                    edgeRecord.directedSgFlux = carrier == "electron"
                        ? cellRecord.electronSignedEdgeFlux[local]
                        : cellRecord.holeSignedEdgeFlux[local];
                    edgeRecord.selectedFluxMagnitude =
                        std::abs(edgeRecord.directedSgFlux);
                    edgeRecord.currentVector =
                        edgeRecord.directedSgFlux *
                        edgeUnit(mesh, node0, node1);
                    edgeRecord.currentVectorProvenance =
                        "solver_native_element_edge_sg_flux";
                    edgeRecord.impactField = carrier == "electron"
                        ? cellRecord.electronImpactField
                        : cellRecord.holeImpactField;
                    edgeRecord.alpha = carrier == "electron"
                        ? cellRecord.electronAlpha
                        : cellRecord.holeAlpha;
                    edgeRecord.sourceMeasure =
                        cellRecord.edgePartialVolumes[local];
                    edgeRecord.contactAdjacent =
                        edgeTouchesContact(contactNodes, node0, node1);
                    edgeRecord.reconstructedCurrent = false;
                    edgeRecord.productionBranchDetails =
                        "source_scatter=element_vertex_records";
                    appendRecord(
                        result, std::move(edgeRecord), impactConfig);

                    BVProcessProbeRecord vertexRecord;
                    vertexRecord.supportKind = "element_vertex_gss_laux";
                    vertexRecord.carrier = carrier;
                    vertexRecord.cellId = cellRecord.cellId;
                    vertexRecord.edgeId = InvalidIndex;
                    vertexRecord.node0 = node0;
                    vertexRecord.node1 = node0;
                    assignEndpointState(vertexRecord, state);
                    vertexRecord.midpointDensity = vertexRecord.density0;
                    vertexRecord.electricFieldVector = electricVector;
                    vertexRecord.quasiFermiGradientVector =
                        carrier == "electron"
                        ? electronQfGradient
                        : holeQfGradient;
                    vertexRecord.currentVector = carrier == "electron"
                        ? cellRecord.electronCurrentVector
                        : cellRecord.holeCurrentVector;
                    vertexRecord.currentVectorProvenance =
                        "solver_native_gss_laux_cell_vector";
                    vertexRecord.lowFieldMobility = 0.5 * (
                        (carrier == "electron"
                             ? cellRecord.electronLowFieldMobilities[local]
                             : cellRecord.holeLowFieldMobilities[local]) +
                        (carrier == "electron"
                             ? cellRecord.electronLowFieldMobilities[
                                   (local + 2) % 3]
                             : cellRecord.holeLowFieldMobilities[
                                   (local + 2) % 3]));
                    vertexRecord.finalMobility = 0.5 * (
                        (carrier == "electron"
                             ? cellRecord.electronMobilities[local]
                             : cellRecord.holeMobilities[local]) +
                        (carrier == "electron"
                             ? cellRecord.electronMobilities[(local + 2) % 3]
                             : cellRecord.holeMobilities[(local + 2) % 3]));
                    vertexRecord.highFieldDrive = 0.5 * (
                        (carrier == "electron"
                             ? cellRecord.electronMobilityDrivingFields[local]
                             : cellRecord.holeMobilityDrivingFields[local]) +
                        (carrier == "electron"
                             ? cellRecord.electronMobilityDrivingFields[
                                   (local + 2) % 3]
                             : cellRecord.holeMobilityDrivingFields[
                                   (local + 2) % 3]));
                    vertexRecord.selectedFluxMagnitude =
                        vertexRecord.currentVector.norm();
                    vertexRecord.impactField = carrier == "electron"
                        ? cellRecord.electronImpactField
                        : cellRecord.holeImpactField;
                    vertexRecord.alpha = carrier == "electron"
                        ? cellRecord.electronAlpha
                        : cellRecord.holeAlpha;
                    vertexRecord.sourceMeasure =
                        cellRecord.vertexMeasures[local];
                    vertexRecord.sourceIntegral = carrier == "electron"
                        ? cellRecord.electronSourceIntegrals[local]
                        : cellRecord.holeSourceIntegrals[local];
                    vertexRecord.scatterCount = 1;
                    vertexRecord.scatterNodes[0] = node0;
                    vertexRecord.sourceWeights[0] = 1.0;
                    vertexRecord.contactAdjacent =
                        contactNodes.at(static_cast<std::size_t>(node0));
                    vertexRecord.reconstructedCurrent = true;
                    appendRecord(
                        result, std::move(vertexRecord), impactConfig);
                }
            }
        }
        return result;
    }

    if (detail::usesTriangleGssAvalancheSource(impactConfig)) {
        const auto production = detail::triangleGssAvalancheSourceRecords(
            impactConfig,
            *impact,
            mobilityConfig,
            *mobility,
            edgeCells,
            mesh,
            doping,
            cellMaterials,
            state.psi,
            state.phin,
            state.phip,
            state.n,
            state.p,
            ni,
            thermalVoltage,
            fieldFactor);
        for (const auto& source : production) {
            const Point2 unit = edgeUnit(mesh, source.node0, source.node1);
            const Point2 electricVector =
                -safeRatio(
                    state.psi(static_cast<Eigen::Index>(source.node1)) -
                        state.psi(static_cast<Eigen::Index>(source.node0)),
                    source.edgeLength) *
                fieldFactor * unit;
            for (const std::string carrier : {"electron", "hole"}) {
                BVProcessProbeRecord record;
                record.supportKind = "triangle_gss_local_edge";
                record.carrier = carrier;
                record.cellId = source.cellId;
                record.localEdge = source.localEdge;
                record.edgeId = source.edgeId;
                record.node0 = source.node0;
                record.node1 = source.node1;
                assignEndpointState(record, state);
                record.midpointDensity = carrier == "electron"
                    ? source.electronMidpointDensity
                    : source.holeMidpointDensity;
                record.electricFieldVector = electricVector;
                record.quasiFermiGradientVector =
                    cellGradient(
                        mesh,
                        source.cellId,
                        carrier == "electron" ? state.phin : state.phip) *
                    fieldFactor;
                record.lowFieldMobility = carrier == "electron"
                    ? source.electronLowFieldMobility
                    : source.holeLowFieldMobility;
                record.highFieldDrive = carrier == "electron"
                    ? source.electronMobilityDrivingField
                    : source.holeMobilityDrivingField;
                record.finalMobility = carrier == "electron"
                    ? source.electronMobility
                    : source.holeMobility;
                record.selectedFluxMagnitude = carrier == "electron"
                    ? source.electronFluxProxy
                    : source.holeFluxProxy;
                const Real qfNorm = record.quasiFermiGradientVector.norm();
                record.currentVector = Point2::Zero();
                if (qfNorm > 0.0) {
                    record.currentVector =
                        -record.selectedFluxMagnitude *
                        record.quasiFermiGradientVector / qfNorm;
                }
                record.currentVectorProvenance =
                    "solver_reconstructed_triangle_gss_gradqf_vector";
                record.impactField = carrier == "electron"
                    ? source.electronImpactField
                    : source.holeImpactField;
                record.alpha = carrier == "electron"
                    ? source.electronAlpha
                    : source.holeAlpha;
                record.sourceMeasure = source.truncatedPartialVolume;
                record.sourceIntegral = carrier == "electron"
                    ? source.electronSourceIntegral
                    : source.holeSourceIntegral;
                record.scatterCount = 2;
                record.scatterNodes[0] = source.node0;
                record.scatterNodes[1] = source.node1;
                record.sourceWeights[0] = 0.5;
                record.sourceWeights[1] = 0.5;
                record.contactAdjacent =
                    edgeTouchesContact(contactNodes, source.node0, source.node1);
                record.reconstructedCurrent = true;
                appendRecord(result, std::move(record), impactConfig);
            }
        }
        return result;
    }

    const auto production = detail::sgEdgeCurrentAvalancheSourceRecords(
        impactConfig,
        *impact,
        mobilityConfig,
        *mobility,
        edgeCells,
        mesh,
        doping,
        cellMaterials,
        state.psi,
        state.phin,
        state.phip,
        state.n,
        state.p,
        ni,
        thermalVoltage,
        fieldFactor,
        true);
    for (const auto& source : production) {
        const Point2 unit = edgeUnit(mesh, source.node0, source.node1);
        const Point2 electricVector =
            -safeRatio(
                state.psi(static_cast<Eigen::Index>(source.node1)) -
                    state.psi(static_cast<Eigen::Index>(source.node0)),
                source.edgeLength) *
            fieldFactor * unit;
        for (const std::string carrier : {"electron", "hole"}) {
            BVProcessProbeRecord record;
            record.supportKind = "sg_edge";
            record.carrier = carrier;
            record.cellId = edgeCells.at(source.edgeId).empty()
                ? InvalidIndex
                : edgeCells.at(source.edgeId).front();
            record.edgeId = source.edgeId;
            record.node0 = source.node0;
            record.node1 = source.node1;
            assignEndpointState(record, state);
            record.midpointDensity =
                0.5 * (record.density0 + record.density1);
            record.electricFieldVector = electricVector;
            record.quasiFermiGradientVector =
                safeRatio(
                    record.quasiFermi1 - record.quasiFermi0,
                    source.edgeLength) *
                fieldFactor * unit;
            record.lowFieldMobility = carrier == "electron"
                ? source.electronLowFieldMobility
                : source.holeLowFieldMobility;
            record.highFieldDrive = carrier == "electron"
                ? source.electronMobilityDrivingField
                : source.holeMobilityDrivingField;
            record.finalMobility = carrier == "electron"
                ? source.electronMobility
                : source.holeMobility;
            record.directedSgFlux = carrier == "electron"
                ? source.electronRawSignedFluxProxy
                : source.holeRawSignedFluxProxy;
            record.selectedFluxMagnitude = carrier == "electron"
                ? source.electronFluxProxy
                : source.holeFluxProxy;
            record.currentVector = record.directedSgFlux * unit;
            record.currentVectorProvenance = "solver_native_directed_sg_flux";
            record.impactField = carrier == "electron"
                ? source.electronImpactField
                : source.holeImpactField;
            record.alpha = carrier == "electron"
                ? source.electronAlpha
                : source.holeAlpha;
            record.sourceMeasure = source.edgeAreaProxy;
            record.sourceIntegral = carrier == "electron"
                ? source.electronSourceIntegral
                : source.holeSourceIntegral;
            std::vector<Real> scatter(
                static_cast<std::size_t>(mesh.numNodes()), 0.0);
            detail::addMappedEdgeSourceToNodes(
                impactConfig,
                scatter,
                edgeCells,
                mesh,
                source,
                carrier == "electron"
                    ? source.electronNode0SourceIntegral
                    : source.holeNode0SourceIntegral,
                carrier == "electron"
                    ? source.electronNode1SourceIntegral
                    : source.holeNode1SourceIntegral,
                record.sourceIntegral);
            assignScatterFromNodeIntegrals(record, scatter);
            record.contactAdjacent =
                edgeTouchesContact(contactNodes, source.node0, source.node1);
            record.reconstructedCurrent =
                impactConfig.currentMagnitudeMode != "edge_scalar_abs" ||
                detail::usesCellCurrentReconstructedAvalancheCurrent(impactConfig) ||
                detail::usesCellVectorCurrentReconstructedAvalancheCurrent(
                    impactConfig);
            record.directionalPartition =
                detail::usesDirectionalEdgeAvalancheSourcePartition(impactConfig);
            if (carrier == "electron" &&
                source.electronSgDiagnosticsCollected) {
                const auto& sg = source.electronSgFluxDecomposition;
                std::ostringstream details;
                details
                    << "sg_flat_qf_short_circuit="
                    << (sg.flatQuasiFermiShortCircuit ? 1 : 0)
                    << ";sg_node0_clamped_low="
                    << (sg.node0ExponentClampedLow ? 1 : 0)
                    << ";sg_node0_clamped_high="
                    << (sg.node0ExponentClampedHigh ? 1 : 0)
                    << ";sg_node1_clamped_low="
                    << (sg.node1ExponentClampedLow ? 1 : 0)
                    << ";sg_node1_clamped_high="
                    << (sg.node1ExponentClampedHigh ? 1 : 0)
                    << ";sg_ni_gradient_drift="
                    << (sg.includeNiGradientDrift ? 1 : 0);
                record.productionBranchDetails = details.str();
            }
            appendRecord(result, std::move(record), impactConfig);
        }
    }
    return result;
}

} // namespace vela
