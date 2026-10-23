#pragma once

#include "vela/core/Types.h"
#include <nlohmann/json_fwd.hpp>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Experimental contracts. Geometry/kernel/source diagnostics are separate;
// no public solver dispatch.
namespace vela::experimental {

enum class TunnelingConstants { Codata2018, SentaurusT202203 };
enum class TunnelingDerivativeMode { NativeDefault, NativeExplicit, FullMathematical };
enum class TunnelingSolveMode { CoupledDC, Gummel, Transient, AC, Noise };

struct IndirectABParameters {
    Real A_per_m3_s = 0; // Normalized to F0=100 V/m; NOT the local E2 A.
    Real B_V_per_m = 0;
    Real referenceGap_J = 0; // Includes constant gap shift, excludes density BGN.
    Real massRatio = 0; // mv/mc (native Rpath); zero selects automatic masses.
};
struct IndirectMicroscopicParameters {
    Real conductionMass_kg = 0;
    Real valenceMass_kg = 0;
    Real coupling_J2_m_per_kg = 0;
};
struct NonlocalPathBandToBandConfig {
    std::string material = "Si";
    std::vector<Index> regions; // Explicit region IDs, no implicit all-material fallback.
    std::optional<IndirectABParameters> ab;
    std::optional<IndirectMicroscopicParameters> microscopic;
    Real phononEnergy_J = 0;
    Real conductionOffset_J = 0;
    Real valenceOffset_J = 0;
    Real maximumPathLength_m = 0;
    Real relativeIntegrationTolerance = 1e-8;
    Index maximumIntegrationIntervals = 1024;
    TunnelingConstants constants = TunnelingConstants::Codata2018;
    TunnelingDerivativeMode derivatives = TunnelingDerivativeMode::FullMathematical;
    void validate(TunnelingSolveMode mode = TunnelingSolveMode::CoupledDC) const;
};

// Audit metadata is deliberately not part of runtime numerical configuration.
struct TunnelingParameterProvenance {
    std::string sourceDigest, sourceVersion, scope;
    std::vector<std::pair<std::string, std::string>> parameterOrigins;
};

// Strict private SI serialization. Never called by public simulation config.
NonlocalPathBandToBandConfig nonlocalPathConfigFromJson(const nlohmann::json&);
nlohmann::json nonlocalPathConfigToJson(const NonlocalPathBandToBandConfig&);
IndirectMicroscopicParameters resolveIndirectParameters(const NonlocalPathBandToBandConfig&);
Real tunnelingElementaryCharge(TunnelingConstants);
struct TunnelingPhysicalConstants {
    Real elementaryCharge_C, planck_J_s, boltzmann_J_per_K, electronMass_kg;
};
TunnelingPhysicalConstants tunnelingPhysicalConstants(TunnelingConstants);

// All values and derivatives use physical, unscaled SI coordinates.
enum class TunnelingVariable { Potential_V, ElectronQF_J, HoleQF_J,
                               ElectronDensity_per_m3, HoleDensity_per_m3, Temperature_K };
struct TunnelingDependency {
    Index node = 0;
    TunnelingVariable variable = TunnelingVariable::Potential_V;
    bool operator==(const TunnelingDependency&) const = default;
};
struct TunnelingBandDerivative {
    TunnelingDependency column;
    Real conduction_J = 0, valence_J = 0, electronFermi_J = 0, holeFermi_J = 0, narrowing_J = 0;
    Point2 conductionGradient_J_per_m = Point2::Zero();
    Point2 valenceGradient_J_per_m = Point2::Zero();
    Real electrons_per_m3 = 0, holes_per_m3 = 0;
    Real conductionDOS_per_m3 = 0, valenceDOS_per_m3 = 0, temperature_K = 0;
};
struct TunnelingBandNode {
    Index material = 0, region = 0;
    bool semiconductor = true;
    bool contact = false;
    Real conduction_J = 0, valence_J = 0;
    Real electronFermi_J = 0, holeFermi_J = 0;
    Real narrowing_J = 0, temperature_K = 300;
    Real electrons_per_m3 = 0, holes_per_m3 = 0;
    Real conductionDOS_per_m3 = 0, valenceDOS_per_m3 = 0;
    Point2 conductionGradient_J_per_m = Point2::Zero();
    Point2 valenceGradient_J_per_m = Point2::Zero();
    std::vector<TunnelingBandDerivative> derivatives;
    // Shared vertices may belong to multiple same-material regions.
    std::vector<Index> additionalRegions;
};
struct TunnelingStateIdentity {
    std::uint64_t mesh = 0, material = 0, parameters = 0, contacts = 0;
    std::uint64_t population = 0, temperature = 0, trial = 0;
    bool operator==(const TunnelingStateIdentity&) const = default;
};
struct TunnelingBandState {
    TunnelingStateIdentity identity;
    std::vector<TunnelingBandNode> nodes;
    void validate() const;
};

struct TunnelingEndpoint {
    Index cell = 0, material = 0, region = 0;
    std::vector<std::pair<Index, Real>> weights; // Barycentric, sum=1, nonnegative.
    void validate(const TunnelingBandState&) const;
};
enum class TunnelingTermination { Complete, LengthLimit, InsulatorBoundary,
    ZeroDirection, ReflectionLimit, NumericalFailure, DomainBoundary, ContactBoundary,
    RegionExcluded, MaterialMismatch, TraversalLimit, AmbiguousCorner, NoBarrier };
struct TunnelingPathPoint {
    Point2 position_m = Point2::Zero();
    TunnelingEndpoint support;
    bool reflection = false;
};
struct TunnelingPath {
    std::uint64_t id = 0;
    TunnelingStateIdentity identity;
    std::vector<TunnelingPathPoint> points;
    Real energyLower_J = 0, energyUpper_J = 0;
    std::vector<std::pair<Real, Real>> integrationIntervals_m;
    TunnelingTermination termination = TunnelingTermination::NumericalFailure;
    void validate(const TunnelingBandState&, Real maximumLength_m) const;
};

struct NonlocalPairEvent {
    std::uint64_t pathId = 0;
    Real signedRate_per_m_s = 0; // Area integral per metre depth; negative is recombination.
    TunnelingEndpoint holeStart, electronEnd;
};
struct NonlocalSourceDerivative {
    bool electron = true;
    Index rowNode = 0;
    TunnelingDependency column;
    Real integralDerivative = 0; // [m^-1 s^-1] per physical column unit.
};
struct NonlocalPairSources {
    TunnelingStateIdentity identity;
    std::vector<Real> electronIntegrals_per_m_s, holeIntegrals_per_m_s;
    std::vector<NonlocalPairEvent> ledger;
    // Empty derivatives without an evaluated mode mean NOT EVALUATED, not zero.
    std::optional<TunnelingDerivativeMode> evaluatedDerivativeMode;
    std::vector<NonlocalSourceDerivative> derivatives;
    void validate(const TunnelingBandState&) const;
    Real balanceError() const; // |sum(e)-sum(h)| / sum(abs(events)); no unit floor.
    std::vector<Real> volumeRates(bool electron, const std::vector<Real>& areas_m2) const;
    Real currentPerUm(bool electron, Real areaFactor, TunnelingConstants) const;
};
// Conservative transfer contract only, not a physical generation kernel.
NonlocalPairSources distributePairEvents(const TunnelingBandState&,
                                        const std::vector<NonlocalPairEvent>&);

class NonlocalJacobianSupport {
public:
    using Entry = std::pair<Index, Index>; // Actual assembler row/column indices.
    // Returns true iff symbolic support changed; shrinking also changes revision.
    bool replace(Index rows, Index columns, std::vector<Entry> entries);
    void require(Index row, Index column) const; // Missing support must rebuild, never drop.
    std::uint64_t revision() const { return revision_; }
    const std::vector<Entry>& entries() const { return entries_; }
private:
    Index rows_ = 0, columns_ = 0;
    std::uint64_t revision_ = 0;
    std::vector<Entry> entries_;
};
} // namespace vela::experimental
