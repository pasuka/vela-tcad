// Included inside vela::detail after the local source integral types.
inline IalElementAvalancheOptions ialElementAvalancheOptionsFromConfig(
    const ImpactIonizationModelConfig& configuredImpact)
{
    const auto impact=applyImpactIonizationParameterSet(configuredImpact);
    const auto& units=impact.unitSystem;
    const auto a=[&](Real v){return units.internalInverseLengthToMInv(v);};
    const auto b=[&](Real v){return units.internalElectricFieldToVPerM(v);};
    IalElementAvalancheOptions options;
    options.electron={a(impact.electronALow),b(impact.electronBLow),a(impact.electronAHigh),b(impact.electronBHigh),b(impact.switchField)};
    options.hole={a(impact.holeALow),b(impact.holeBLow),a(impact.holeAHigh),b(impact.holeBHigh),b(impact.switchField)};
    options.minimumField_V_per_m=b(impact.minimumField);
    return options;
}

template <typename Scalar>
inline ElementEdgeGssLauxLocalSourceIntegrals<Scalar>
ialElementAvalancheSourceIntegralsLocal(
    const ImpactIonizationModelConfig& impact,
    const MobilityModelConfig& mobility,const DeviceMesh& mesh,
    const DopingModel& doping,Index cellId,
    const std::array<Scalar,3>& psi,const std::array<Scalar,3>& fn,
    const std::array<Scalar,3>& fp,const std::array<Scalar,3>& n,
    const std::array<Scalar,3>& p)
{
    ElementEdgeGssLauxLocalSourceIntegrals<Scalar> result;
    if constexpr (!std::is_same_v<Scalar,Real> && !std::is_same_v<Scalar,Tri3LocalForwardDual>) {
        throw std::invalid_argument("Observed IALMob source has no multiprecision audit implementation");
    } else {
        const auto& cell=mesh.getCell(cellId);
        if (mesh.getRegion(cell.region_id).material!="Si") return result;
        if (cell.node_ids.size()!=3) throw std::invalid_argument("Observed IALMob source requires Tri3");
        std::array<IalElementVertexState,3> state;
        const Real concentration=mobility.internalConcentrationToM3;
        for (std::size_t i=0;i<3;++i) {
            const auto node=cell.node_ids[i];
            state[i]={localAdValue(psi[i]),localAdValue(fn[i]),localAdValue(fp[i]),
                doping.donors(node)*concentration,doping.acceptors(node)*concentration,
                localAdValue(n[i])*concentration,localAdValue(p[i])*concentration};
            if constexpr (std::is_same_v<Scalar,Tri3LocalForwardDual>) {
                // The coupled assembler seeds block-ordered physical potentials.
                state[i].electronResponse_m3_per_V=n[i].derivative[i]*concentration;
                state[i].holeResponse_m3_per_V=p[i].derivative[6+i]*concentration;
            }
        }
        const auto options=ialElementAvalancheOptionsFromConfig(impact);
        const auto evaluated=evaluateIalTransportAvalancheCell(mobility,cellId,state,options);
        const Real measure=concentration*mobility.internalLengthToM*mobility.internalLengthToM;
        const auto convert=[&](const Tri3LocalForwardDual& v) {
            Scalar out(v.value/measure);
            if constexpr (std::is_same_v<Scalar,Tri3LocalForwardDual>)
                for (std::size_t i=0;i<3;++i) for (std::size_t k=0;k<3;++k)
                    out.derivative[3*k+i]=v.derivative[3*i+k]/measure;
            return out;
        };
        for (std::size_t i=0;i<3;++i) {
            result.electron[i]=convert(evaluated.electron.source_per_m_s[i]);
            result.hole[i]=convert(evaluated.hole.source_per_m_s[i]);
            result.combined[i]=result.electron[i]+result.hole[i];
        }
        return result;
    }
}

#include "vela/equation/PhuMobElementAvalancheAssembly.inl"

template <typename Scalar>
inline ElementEdgeGssLauxLocalSourceIntegrals<Scalar>
elementLocalAvalancheSourceIntegrals(
    const ImpactIonizationModelConfig& impact,const MobilityModelConfig& mobility,
    const MobilityModel& model,const std::vector<Index>& edges,
    const DeviceMesh& mesh,const DopingModel& doping,const std::vector<Material>& materials,
    Index cell,const std::array<Scalar,3>& psi,const std::array<Scalar,3>& fn,
    const std::array<Scalar,3>& fp,const std::array<Scalar,3>& n,const std::array<Scalar,3>& p,
    const std::array<Real,3>& ni,Real vt,Real fieldFactor)
{
    if (usesPhuMobElementAvalancheSource(impact))
        return phumobElementAvalancheSourceIntegralsLocal(impact,mobility,mesh,doping,materials,cell,psi,fn,fp,n,p);
    if (usesIalElementAvalancheSource(impact))
        return ialElementAvalancheSourceIntegralsLocal(impact,mobility,mesh,doping,cell,psi,fn,fp,n,p);
    return elementEdgeGssLauxAvalancheSourceIntegralsLocal<Scalar>(
        impact,mobility,model,edges,mesh,doping,materials,cell,psi,fn,fp,n,p,ni,vt,fieldFactor);
}
