#pragma once
#include "vela/core/Types.h"
#include <array>
#include <memory>

namespace vela {
/// Explicit SI local model for the audited Siliconc100.par Formula-1 DOS,
/// OldSlotboom, Fermi, Scharfetter TempDependence and Auger H/N0 combination.
/// This interface does not change Material::atTemperature or solver defaults.
struct SiliconThermalParameters {
    Real bandgap300_eV=1.12416, affinity300_eV=4.0727;
    Real bandgapAlpha_eV_per_K=4.73e-4, bandgapBeta_K=636.;
    Real bgnToAffinity=.5, bgnCoefficient_eV=.009, bgnReference_m3=1e23;
    bool fermiBgnCorrection=true;
    Real electronTransverseMass=.1905, electronLongitudinalMass=.9163, electronMassOffset=0.;
    std::array<Real,5> holeMassNumerator{.443587,.003609528,.0001173515,1.263218e-6,3.025581e-9};
    std::array<Real,5> holeMassDenominator{1.,.004683382,.0002286895,7.469271e-7,1.727481e-9};
    Real holeMassOffset=0.;
    std::array<Real,2> srhTau300_s{1e-5,3e-6}, srhReference_m3{1e22,1e22};
    std::array<Real,2> srhTemperatureExponent{-1.5,-1.5};
    std::array<Real,3> augerElectron_m6_per_s{6.7e-44,2.45e-43,-2.2e-44};
    std::array<Real,3> augerHole_m6_per_s{7.2e-44,4.5e-45,2.63e-44};
    std::array<Real,2> augerEnhancement{3.46667,8.25688}, augerReference_m3{1e24,1e24};
    // Original D0 specifies Auger without WithGeneration (UG p.489).
    bool augerWithGeneration=false;
};
struct ThermalQuantity {
    Real value=0.;
    /// Partial order: psi [V], electron QF [V], hole QF [V], T [K].
    std::array<Real,4> derivative{};
};
struct SiliconThermalState {
    Real potential_V=0.,electronQf_V=0.,holeQf_V=0.,temperature_K=300.;
    Real donors_m3=0.,acceptors_m3=0.;
    /// Optional QF reference plus the corresponding increment above.
    Real electronQfReference_V=0.,holeQfReference_V=0.;
};
struct SiliconThermalResult {
    ThermalQuantity bandgap_eV,affinity_eV,Nc_m3,Nv_m3,ni_m3,effectiveNi_m3;
    ThermalQuantity conductionBand_eV,valenceBand_eV,electrons_m3,holes_m3;
    ThermalQuantity electronEta,holeEta;
    ThermalQuantity srhRate_m3_per_s,augerRate_m3_per_s;
    ThermalQuantity electronLifetime_s,holeLifetime_s,augerElectron_m6_per_s,augerHole_m6_per_s;
    Real bandgapNarrowing_eV=0.;
};
class SiliconThermalPhysics {
public:
    class DopingPreparation {
    public:
        bool matches(Real donors,Real acceptors) const {return donors_==donors && acceptors_==acceptors;}
    private:
        friend class SiliconThermalPhysics;
        Real donors_=0.,acceptors_=0.,delta_=0.;
        std::shared_ptr<const unsigned char> owner_;
    };
    class TemperaturePreparation {
    public:
        bool matches(Real temperature,Real donors,Real acceptors) const {
            return temperature_==temperature && donors_==donors && acceptors_==acceptors;
        }
    private:
        friend class SiliconThermalPhysics;
        Real temperature_=0.,donors_=0.,acceptors_=0.;
        std::shared_ptr<const unsigned char> owner_;
        SiliconThermalResult base_;
        ThermalQuantity vt_,effectiveEg_;
    };
    explicit SiliconThermalPhysics(SiliconThermalParameters parameters={});
    /// Reference fixed at intrinsic silicon at 300 K, independent of local T.
    Real referencePotential_V() const {return reference_;}
    SiliconThermalResult evaluate(const SiliconThermalState&) const;
    /// Preparations retain the model identity; stale temperature/doping or a
    /// preparation from another parameter instance is rejected, never reused.
    DopingPreparation prepareDoping(Real donors_m3,Real acceptors_m3) const;
    TemperaturePreparation prepareTemperature(Real temperature_K,const DopingPreparation&) const;
    SiliconThermalResult evaluate(const SiliconThermalState&,const TemperaturePreparation&) const;
    /// Neutrality needs densities/partials only, without recombination work.
    std::array<ThermalQuantity,2> carrierDensities(const SiliconThermalState&,const TemperaturePreparation&) const;
private:
    struct CarrierEvaluation {SiliconThermalResult result;ThermalQuantity fN,fP;};
    CarrierEvaluation carriers(const SiliconThermalState&,const TemperaturePreparation&) const;
    SiliconThermalParameters p_;
    Real reference_,Nc300_,Nv300_;
    std::shared_ptr<const unsigned char> identity_=std::make_shared<const unsigned char>(0);
};
} // namespace vela
