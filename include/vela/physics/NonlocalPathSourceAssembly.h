#pragma once
#include "vela/physics/NonlocalPathPhononKernel.h"

namespace vela::experimental {
enum class TunnelingBoxPolicy { Barycentric, Delaunay };
struct TunnelingSourceBoxes {
    std::vector<Real> seedAreas_m2, fieldAreas_m2;
    std::vector<Real> cellSeedAreas_m2;
};
// Geometry-only construction; no native rate, fitted weight or source table.
TunnelingSourceBoxes tunnelingSourceBoxes(const TunnelingGeometryMesh&,
    const NonlocalPathBandToBandConfig&, TunnelingBoxPolicy);
struct TunnelingSeedWindow {
    Real valenceAvailableEnergy_J = 0, conductionAvailableEnergy_J = 0;
};
struct TunnelingWindowInput {
    TunnelingStateIdentity identity;
    std::vector<TunnelingSeedWindow> nodes;
};
struct TunnelingSourceRecord {
    TunnelingGeometryResult geometry;
    Real seedArea_m2 = 0;
    std::optional<PhononTunnelingResult> kernel;
    Real signedEvent_per_m_s = 0;
    bool resolved = false;
    std::string disposition;
    std::vector<NonlocalPairEvent> collectionEvents;
    Index energySamples = 0, minimumPolicyStops = 0, depthPolicyStops = 0;
};
struct TunnelingContactLedger {
    Real freeElectron_per_m_s = 0, contactElectron_per_m_s = 0;
    Real freeHole_per_m_s = 0, contactHole_per_m_s = 0;
};
struct TunnelingSourceAssembly {
    TunnelingConstants constants = TunnelingConstants::Codata2018;
    TunnelingGeometryMesh sourceMesh;
    NonlocalPairSources sources;
    TunnelingSourceBoxes boxes;
    std::vector<TunnelingSourceRecord> records;
    TunnelingContactLedger contacts;
    bool resolved = false;
    bool cellSeeded = false;
    std::vector<Index> selectedRegions;
    // Recheck all paths and the ledger. Partial diagnostics cannot enter RHS.
    void requireResolved(const TunnelingBandState&) const;
};
// One autonomous N2 path per selected Si vertex, with geometry-derived box
// quadrature and actual endpoint Fermi interpolation. Explicit window input
// remains necessary: this API does NOT claim native window sampling identity.
TunnelingSourceAssembly assembleNonlocalPathSources(const TunnelingGeometryMesh&,
    const TunnelingBandState&, const NonlocalPathBandToBandConfig&,
    const TunnelingWindowInput&, TunnelingBoxPolicy,
    const TunnelingGeometryOptions& = {});
nlohmann::json tunnelingSourceAudit(const TunnelingSourceAssembly&,
    const TunnelingBandState&);
} // namespace vela::experimental
