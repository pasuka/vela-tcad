#include "vela/equation/PhuMobSgSourceCurrent.inl"

// Included inside detail after the Tri3 local AD helpers. This explicit
// Si-only profile uses live Fermi or Boltzmann populations supplied by assembly.
struct PhuMobElementAvalancheDiagnostics {
    std::array<Real,3> electronMobility{},holeMobility{};
    std::array<Real,3> electronLowFieldMobility{},holeLowFieldMobility{};
    Real electronAlpha=0.,holeAlpha=0.,electronDrive=0.,holeDrive=0.;
    Real electronHighFieldDrive=0.,holeHighFieldDrive=0.;
    std::array<Real,2> electronSgCurrent{},holeSgCurrent{};
    std::size_t sourceVertex=0;
};
template<class S>
inline ElementEdgeGssLauxLocalSourceIntegrals<S> phumobElementAvalancheSourceIntegralsLocal(
    const ImpactIonizationModelConfig& impact,const MobilityModelConfig& mobility,
    const DeviceMesh& mesh,const DopingModel& doping,const std::vector<Material>& materials,
    Index cid,const std::array<S,3>& psi,const std::array<S,3>& fn,const std::array<S,3>& fp,
    const std::array<S,3>& n,const std::array<S,3>& p,
    PhuMobElementAvalancheDiagnostics* diagnostics=nullptr,
    bool diagnosticFreezeMobilityQfChain=false)
{
    ElementEdgeGssLauxLocalSourceIntegrals<S> out;
    if constexpr (!std::is_same_v<S,Real> && !std::is_same_v<S,Tri3LocalForwardDual>) {
        throw std::invalid_argument("PhuMob element source has no multiprecision audit implementation");
    } else {
    const auto& cell=mesh.getCell(cid);
    if(mesh.getRegion(cell.region_id).material!="Si") return out;
    if(mobility.edgeAveraging!="element_box_phumob" || !usesElementHighField(mobility))
        throw std::invalid_argument("PhuMob element avalanche requires full element Lombardi/GradQF transport");
    const bool sg=impact.currentApproximation=="phumob_element_sg_best_vertex";
    const auto gp=localAdTri3Gradient(mesh,cell,psi);
    const auto gn=localAdTri3Gradient(mesh,cell,fn);
    const auto gh=localAdTri3Gradient(mesh,cell,fp);
    const Real toField=1./mobility.internalLengthToM/mobility.internalFieldToVPerM;
    const auto norm=[&](const auto& g){return localAdNorm2(g[0],g[1])*S(toField);};
    const S en=norm(gn),eh=norm(gh),electric=norm(gp);
    const bool contact=mobility.highFieldCells.at(cid).contact;
    const auto highField=[&](const std::array<S,3>& qf) {
        const auto& z=contact && mobility.useElectricFieldAtGradQfContacts(true)?psi:qf;
        const auto& s=mobility.highFieldCells.at(cid).stencil;
        return localAdNorm2((z[1]-z[0])*S(s[0])+(z[2]-z[0])*S(s[2]),
                            (z[1]-z[0])*S(s[1])+(z[2]-z[0])*S(s[3]));
    };
    const S driveN=highField(fn),driveP=highField(fp);
    const auto projected=[&](const std::array<S,3>& z) {
        const auto& st=mobility.highFieldCells.at(cid).stencil;
        return localAdNorm2((z[1]-z[0])*S(st[0])+(z[2]-z[0])*S(st[2]),
                            (z[1]-z[0])*S(st[1])+(z[2]-z[0])*S(st[3]));
    };
    const S forceN=impact.useElectricFieldAtGradQfContacts()&&contact?electric:(sg?projected(fn):en);
    const S forceP=impact.useElectricFieldAtGradQfContacts()&&contact?electric:(sg?projected(fp):eh);
    const S alphaN=localAdImpactCoefficient(impact,CarrierType::Electron,
        forceN);
    const S alphaP=localAdImpactCoefficient(impact,CarrierType::Hole,
        forceP);
    if(diagnostics) {
        diagnostics->electronAlpha=localAdValue(alphaN);diagnostics->holeAlpha=localAdValue(alphaP);
        diagnostics->electronDrive=localAdValue(forceN);
        diagnostics->holeDrive=localAdValue(forceP);
        diagnostics->electronHighFieldDrive=localAdValue(driveN);diagnostics->holeHighFieldDrive=localAdValue(driveP);
    }
    S normal(0.);
    if(std::isfinite(mobility.surface.cellNormalX.at(cid)) && std::isfinite(mobility.surface.cellNormalY.at(cid)))
        normal=localAdAbs(gp[0]*S(mobility.surface.cellNormalX[cid])+gp[1]*S(mobility.surface.cellNormalY[cid]))*S(toField);
    const auto measures=cellBoxNodeMeasures(mesh,cid);
    const auto& material=materials.at(cid);
    std::size_t best=0;
    for(std::size_t k=1;k<3;++k)
        if(measures[k]>measures[best] ||
           (measures[k]==measures[best] && cell.node_ids[k]<cell.node_ids[best])) best=k;
    const Real vt=constants::kb*material.temperature_K.value_or(300.)/constants::q;
    auto jn=sg?phumobSourceSgVector(mesh,cell,best,n,fn,vt,true):std::array<S,2>{};
    auto jp=sg?phumobSourceSgVector(mesh,cell,best,p,fp,vt,false):std::array<S,2>{};
    if(sg && !contact) {
        // Native partial-layer source projects the SG vector as well as alpha's
        // force. Apply the same geometric map to an affine scalar with gradient J.
        const auto projectCurrent=[&](const std::array<S,2>& j) {
            const auto& a=mesh.getNode(cell.node_ids[0]);
            const auto& b=mesh.getNode(cell.node_ids[1]);
            const auto& c=mesh.getNode(cell.node_ids[2]);
            const S v1=j[0]*S(b.x-a.x)+j[1]*S(b.y-a.y);
            const S v2=j[0]*S(c.x-a.x)+j[1]*S(c.y-a.y);
            const auto& st=mobility.highFieldCells.at(cid).stencil;
            return std::array<S,2>{(v1*S(st[0])+v2*S(st[2]))/S(toField),
                                   (v1*S(st[1])+v2*S(st[3]))/S(toField)};
        };
        jn=projectCurrent(jn);jp=projectCurrent(jp);
    }
    const Real toSource=mobility.internalMobilityToM2PerVS*mobility.internalFieldToVPerM*
        impact.unitSystem.internalInverseLengthToMInv(1.);
    for(std::size_t k=0;k<3;++k) {
        const std::size_t v=sg?best:k;
        const Index node=cell.node_ids[v];
        auto mu=[&](CarrierType carrier,const S& field) {
            const Real nv=localAdValue(n[v]),pv=localAdValue(p[v]);
            PhuMobScalarState state{doping.donors(node),doping.acceptors(node),nv,pv,material.temperature_K.value_or(300.)};
            const Real value=evaluatePhuMobScalar(carrier,state,mobility.phuMob).mobility;
            const auto d=evaluatePhuMobLogDensityDerivatives(carrier,state,mobility.phuMob);
            S bulk(value);
            if(nv>0.) bulk=bulk+S(d.electrons/nv)*(n[v]-S(nv));
            if(pv>0.) bulk=bulk+S(d.holes/pv)*(p[v]-S(pv));
            const auto surface=evaluateElementLombardi(value,localAdValue(normal),mobility.surface.nodeDistances.at(node),
                doping.totalImpurity(node),state.temperature_K,carrier==CarrierType::Electron?mobility.electronLombardi:mobility.holeLombardi,mobility);
            if(diagnostics) (carrier==CarrierType::Electron?diagnostics->electronLowFieldMobility:diagnostics->holeLowFieldMobility)[k]=surface.mobility;
            S low=S(surface.mobility)+S(surface.bulkDerivative)*(bulk-S(value))+
                S(surface.fieldDerivative)*(normal-S(localAdValue(normal)));
            const auto high=evaluateElementCanali(surface.mobility,localAdValue(field),
                carrier==CarrierType::Electron?mobility.electronField:mobility.holeField,100./mobility.internalFieldToVPerM);
            // Contact HFS is driven by psi, so retain that chain. Also retain
            // population, low-field and normal-field derivatives everywhere.
            const Real fieldDerivative = diagnosticFreezeMobilityQfChain && !contact
                ? 0. : high.fieldDerivative;
            return S(high.mobility)+S(high.lowDerivative)*(low-S(surface.mobility))+
                S(fieldDerivative)*(field-S(localAdValue(field)));
        };
        // mu*n*|grad phi| is particle current. Multiplication by alpha and
        // the cell vertex measure yields particle generation per width/time.
        const S mun=mu(CarrierType::Electron,driveN),mup=mu(CarrierType::Hole,driveP);
        if(diagnostics) {diagnostics->electronMobility[k]=localAdValue(mun);diagnostics->holeMobility[k]=localAdValue(mup);}
        if(sg) {
            out.electron[k]=alphaN*mun*norm(jn)*S(measures[k]*toSource);
            out.hole[k]=alphaP*mup*norm(jp)*S(measures[k]*toSource);
        } else {
            out.electron[k]=alphaN*mun*n[k]*en*S(measures[k]*toSource);
            out.hole[k]=alphaP*mup*p[k]*eh*S(measures[k]*toSource);
        }
        out.combined[k]=out.electron[k]+out.hole[k];
        if(sg && diagnostics) {
            diagnostics->sourceVertex=best;
            for(int d=0;d<2;++d) {
                diagnostics->electronSgCurrent[d]=localAdValue(mun*jn[d])*toField;
                diagnostics->holeSgCurrent[d]=localAdValue(mup*jp[d])*toField;
            }
        }
    }
    return out;
    }
}
