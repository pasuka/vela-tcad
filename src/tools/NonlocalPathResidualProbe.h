#pragma once
// Uninstalled fixed-state diagnostic, included only by the N4 source probe.
#include "vela/equation/CoupledDDAssembler.h"
#include "vela/solver/NewtonSolver.h"
#include "vela/core/PhysicalConstants.h"
#include <map>

inline nlohmann::json nonlocalResidualProbe(const nlohmann::json& j,
    const vela::experimental::TunnelingGeometryMesh& sm,
    const vela::experimental::TunnelingBandState& bands,
    const vela::experimental::TunnelingSourceAssembly& source) {
    using namespace vela;
    const auto& input=j.at("actual_residual");
    UnitScalingConfig units(UnitScalingMode::UnitScaling);const auto u=units.unitSystem();
    const auto cfg=newtonConfigFromJson(input.at("solver"),units);
    DeviceMesh mesh;
    for(Index i=0;i<sm.positions_m.size();++i) {
        Node n;n.id=i;n.x=u.metersToInternalLength(sm.positions_m[i].x());n.y=u.metersToInternalLength(sm.positions_m[i].y());mesh.addNode(n);
    }
    std::map<Index,Region> regions;
    for(Index i=0;i<sm.cells.size();++i) {
        const auto& s=sm.cells[i];Cell c;c.id=i;c.type=CellType::Tri3;c.region_id=s.region;c.node_ids.assign(s.nodes.begin(),s.nodes.end());mesh.addCell(c);
        auto& r=regions[s.region];r.id=s.region;r.name="region_"+std::to_string(s.region);
        r.material=input.at("region_materials").at(std::to_string(s.region)).get<std::string>();r.cell_ids.push_back(i);
    }
    for(auto& [id,r]:regions) mesh.addRegion(r);
    mesh.buildEdges();BoxGeometryBuilder::Options geometry;geometry.nodeVolumePolicy=BoxGeometryBuilder::NodeVolumePolicy::MixedVoronoi;
    geometry.poissonPermittivityPolicy=BoxGeometryBuilder::PoissonPermittivityPolicy::CellMaterial;
    geometry.cellBoxPolicy=BoxGeometryBuilder::CellBoxPolicy::DelaunayTransfer;
    mesh.buildBoxGeometry(geometry);
    MaterialDatabase materials(input.at("materials_file").get<std::string>(),units);
    DopingModel doping(mesh.numNodes());
    for(Index i=0;i<mesh.numNodes();++i) {
        doping.setNodeDoping(i,input.at("donors_cm3").at(i).get<Real>(),input.at("acceptors_cm3").at(i).get<Real>());
    }
    const auto constants=experimental::tunnelingPhysicalConstants(source.constants);
    DDScalingSpec scale;scale.enabled=true;scale.unitSystem=u;
    scale.V0=constants.boltzmann_J_per_K*cfg.temperature_K/constants.elementaryCharge_C;
    scale.C0=input.at("concentration_reference_cm3").get<Real>();scale.D0=1;scale.mu0=1/scale.V0;
    scale.L0=u.metersToInternalLength(input.at("length_reference_m").get<Real>());
    scale.permittivityReference_F_per_m=constants.elementaryCharge_C*scale.C0*1e6*
        std::pow(input.at("length_reference_m").get<Real>(),2)/scale.V0;
    scale.chargeAreaFactor=u.chargeAreaFactor();scale.chargeLineFactor=u.chargeLineFactor();
    scale.fieldFromCoordinateDeltaFactor=u.fieldFromCoordinateDeltaFactor();
    scale.currentDensityLineIntegralFactor=u.currentDensityAM2PerInternal()*u.lengthMPerInternal();
    scale.regionResolvedInterfaceAssembly=cfg.regionResolvedInterfaceAssembly;
    auto recombination=recombinationModelConfig(cfg.recombination,cfg.taun,cfg.taup,cfg.srhDopingDependence);
    // The residual guard below rejects a second local BTBT source.
    recombination.bandToBand=cfg.bandToBand;
    // V0 is a coordinate normalization, not the core's physical thermal
    // voltage. Keep the core temperature at the explicitly configured T.
    const Real coreVt=vela::constants::kb*cfg.temperature_K/vela::constants::q;
    CoupledDDAssembler dd(mesh,materials,doping,coreVt,cfg.mobility,recombination,
        cfg.bandgapNarrowing,cfg.impactIonization,{},{},scale,{},cfg.carrierStatistics);
    const Index N=mesh.numNodes();VectorXd x=VectorXd::Zero(3*N),en(N),hp(N),refE(N),refH(N);
    for(Index i=0;i<N;++i) {
        x(i)=input.at("potential_V").at(i).get<Real>()/scale.V0;
        refE(i)=-bands.nodes[i].electronFermi_J/constants.elementaryCharge_C;
        refH(i)=-bands.nodes[i].holeFermi_J/constants.elementaryCharge_C;
        en(i)=bands.nodes[i].semiconductor?bands.nodes[i].electrons_per_m3/1e6:1;
        hp(i)=bands.nodes[i].semiconductor?bands.nodes[i].holes_per_m3/1e6:1;
    }
    dd.setQuasiFermiReferenceFields(refE,refH);
    CoupledDDBoundaryConditions bc;
    for(Index i=0;i<N;++i) if(bands.nodes[i].contact) {
        bc.psi[i]=x(i);if(bands.nodes[i].semiconductor) {bc.phin[i]=refE(i)/scale.V0;bc.phip[i]=refH(i)/scale.V0;}
    }
    nlohmann::json result{{"population_frozen",true},{"constitutive_mapping_qualified",false},
        {"temperature_K",cfg.temperature_K},{"core_thermal_voltage_V",coreVt},{"coordinate_scale_V",scale.V0},
        {"core_constants","CODATA2018"},{"source_constants","explicit source assembly policy"},
        {"QF_representation","exported physical QF as reference, zero increment; internal low part unqualified"},
        {"boundary_values","frozen state; not an independently coupled resistor load row"},
        {"backend","residual assembly; no linear solve"},{"rows",nlohmann::json::array()}};
    try {const auto r=dd.experimentalNonlocalResidual(x,bc,source,bands);result["strict_state_accepted"]=true;}
    catch(const std::invalid_argument& e) {result["strict_state_accepted"]=false;result["strict_state_reason"]=e.what();}
    const auto rn=dd.electronDensity(x),rp=dd.holeDensity(x);
    Real nerr=0,perr=0;
    for(Index i=0;i<N;++i) if(bands.nodes[i].semiconductor) {
        nerr=std::max(nerr,std::abs(std::log(rn(i)/en(i))));perr=std::max(perr,std::abs(std::log(rp(i)/hp(i))));
    }
    result["reconstructed_population_max_abs_log_ratio"]={{"electron",nerr},{"hole",perr}};
    CoupledDDFeedbackStateSubstitution frozen;frozen.replaceElectronDensity=frozen.replaceHoleDensity=true;
    frozen.electronDensity=en;frozen.holeDensity=hp;
    const auto base=dd.feedbackSubstitutionResidual(x,bc,frozen);
    const auto actual=dd.experimentalNonlocalFrozenResidual(x,bc,source,bands,frozen);
    const auto terms=dd.feedbackSubstitutionCarrierContinuityTermDiagnostics(x,bc,frozen);
    const auto poisson=dd.feedbackSubstitutionPoissonTermDiagnostics(x,bc,frozen);
    const Real factor=1/(u.currentDensityAM2PerInternal()*u.lengthMPerInternal()*dd.continuityResidualScale());
    for(Index i=0;i<N;++i) {
        const auto& t=terms[i];nlohmann::json row{{"node",i},{"semiconductor",bands.nodes[i].semiconductor},{"contact",bands.nodes[i].contact},
            {"psi",actual(i)},{"electron",actual(N+i)},{"hole",actual(2*N+i)}};
        for(bool e:{true,false}) {
            const Real sourceTerm=-(e?source.sources.electronIntegrals_per_m_s[i]:source.sources.holeIntegrals_per_m_s[i])*factor;
            const Real flux=e?t.electronFlux:t.holeFlux,fluxAbs=e?t.electronFluxAbsSum:t.holeFluxAbsSum;
            const Real rec=e?t.electronRecombination:t.holeRecombination,impact=e?t.electronImpact:t.holeImpact;
            const Real gauge=e?t.electronGauge:t.holeGauge,boundary=e?t.electronBoundary:t.holeBoundary;
            const bool constrained=e?bc.phin.contains(i):bc.phip.contains(i);
            const Real expected=(e?base(N+i):base(2*N+i))+(constrained?0:sourceTerm);
            row[e?"electron_terms":"hole_terms"]={{"transport",flux},{"absolute_transport",fluxAbs},
                {"SRH",rec},{"avalanche",impact},{"nonlocal",sourceTerm},{"gauge",gauge},{"boundary",boundary},
                {"actual_source_addition_error",(e?actual(N+i):actual(2*N+i))-expected},
                {"absolute_Vela_assembly_terms",fluxAbs+std::abs(rec)+std::abs(impact)+std::abs(sourceTerm)+std::abs(gauge)+std::abs(boundary)}};
        }
        const auto& p=poisson[i];row["poisson_terms"]={{"dielectric",p.dielectricFlux},{"electron_charge",p.electronCharge},
            {"absolute_dielectric_edges",p.dielectricFluxAbsSum},
            {"hole_charge",p.holeCharge},{"dopant_charge",p.dopingCharge},{"fixed_charge",p.fixedInterfaceCharge},
            {"mapping_error_e_V",p.electronQuasiFermiMappingError},{"mapping_error_h_V",p.holeQuasiFermiMappingError},
            {"uses_reconstructed_population",false},{"ledger_closure_error",p.closureError}};
        result["rows"].push_back(std::move(row));
    }
    result["full_rhs_qualified"]=false;return result;
}
