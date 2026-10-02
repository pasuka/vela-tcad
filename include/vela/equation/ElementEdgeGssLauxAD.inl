// Included from AssemblerUtils.h inside namespace vela::detail after the
// Tri3 geometry/current-reconstruction helpers have been declared.

template <typename Scalar>
inline auto localAdValue(const Scalar& value)
{
    if constexpr (std::is_same_v<Scalar, Tri3LocalForwardDual>)
        return value.value;
    else
        return value;
}

template <typename Scalar>
inline Scalar localAdAbs(const Scalar& value)
{
    if constexpr (std::is_same_v<Scalar, Tri3LocalForwardDual>)
        return dualAbs(value);
    else {
        using std::abs;
        return abs(value);
    }
}

template <typename Scalar>
inline Scalar localAdExp(const Scalar& value)
{
    if constexpr (std::is_same_v<Scalar, Tri3LocalForwardDual>)
        return dualExp(value);
    else {
        using std::exp;
        return exp(value);
    }
}

template <typename Scalar>
inline Scalar localAdExpm1(const Scalar& value)
{
    if constexpr (std::is_same_v<Scalar, Tri3LocalForwardDual>)
        return dualExpm1(value);
    else {
        using std::expm1;
        return expm1(value);
    }
}

template <typename Scalar>
inline Scalar localAdSqrt(const Scalar& value)
{
    if constexpr (std::is_same_v<Scalar, Tri3LocalForwardDual>)
        return dualSqrt(value);
    else {
        using std::sqrt;
        return sqrt(value);
    }
}

template <typename Scalar>
inline Scalar localAdPow(const Scalar& value, Real exponent)
{
    if constexpr (std::is_same_v<Scalar, Tri3LocalForwardDual>)
        return dualPow(value, exponent);
    else {
        using std::pow;
        return pow(value, exponent);
    }
}

template <typename Scalar>
inline Scalar localAdLimitedExp(const Scalar& exponent)
{
    if (localAdValue(exponent) <= -500.0)
        return Scalar(std::exp(-500.0));
    if (localAdValue(exponent) >= 500.0)
        return Scalar(std::exp(500.0));
    return localAdExp(exponent);
}

template <typename Scalar>
inline Scalar localAdBernoulli(const Scalar& argument)
{
    const auto value = localAdValue(argument);
    if (localAdValue(localAdAbs(argument)) < 1.0e-10)
        return Scalar(1.0) - argument * Scalar(0.5) +
               argument * argument / Scalar(12.0);
    if (value > 500.0)
        return argument * localAdExp(-argument);
    if (value < -500.0)
        return -argument;
    return argument / localAdExpm1(argument);
}

template <typename Scalar>
inline Scalar localAdNorm2(const Scalar& x, const Scalar& y)
{
    return localAdSqrt(x * x + y * y);
}

template <typename Scalar>
inline std::array<Scalar, 2> localAdTri3Gradient(
    const DeviceMesh& mesh,
    const Cell& cell,
    const std::array<Scalar, 3>& values)
{
    const Node& p0 = mesh.getNode(cell.node_ids[0]);
    const Node& p1 = mesh.getNode(cell.node_ids[1]);
    const Node& p2 = mesh.getNode(cell.node_ids[2]);
    const Real dx10 = p1.x - p0.x;
    const Real dy10 = p1.y - p0.y;
    const Real dx20 = p2.x - p0.x;
    const Real dy20 = p2.y - p0.y;
    const Real determinant = dx10 * dy20 - dy10 * dx20;
    if (std::abs(determinant) <= 1.0e-300) {
        throw std::invalid_argument(
            "degenerate triangle cannot evaluate avalanche driving fields");
    }
    const Scalar dv10 = values[1] - values[0];
    const Scalar dv20 = values[2] - values[0];
    return {
        (dv10 * Scalar(dy20) - dv20 * Scalar(dy10)) /
            Scalar(determinant),
        (Scalar(dx10) * dv20 - Scalar(dx20) * dv10) /
            Scalar(determinant)};
}

template <typename Scalar>
inline Scalar localAdContactElectricFallbackImpactField(
    const ImpactIonizationModelConfig& config,
    const DeviceMesh& mesh,
    const Cell& cell,
    const std::array<Scalar, 2>& electricGradient,
    const Scalar& quasiFermiField,
    const std::array<Scalar, 3>& potential,
    Real fieldFactor)
{
    const Scalar electricMagnitude = localAdNorm2(
        electricGradient[0], electricGradient[1]) * Scalar(fieldFactor);
    if (config.contactElectricFieldFallbackMode == "cell_gradient_magnitude")
        return electricMagnitude;

    ContactBoundaryFaceGeometry geometry;
    if (!contactBoundaryFaceGeometry(mesh, cell, geometry))
        return quasiFermiField;
    const Scalar faceNormalField = localAdAbs(
        electricGradient[0] * Scalar(geometry.unitNormal.x()) +
        electricGradient[1] * Scalar(geometry.unitNormal.y())) *
        Scalar(fieldFactor);
    if (config.contactElectricFieldFallbackMode == "face_normal")
        return faceNormalField;

    const std::size_t face0 = geometry.faceLocalNodes[0];
    const std::size_t face1 = geometry.faceLocalNodes[1];
    const std::size_t opposite = geometry.oppositeLocalNode;
    const Scalar oneSidedField = localAdAbs(
        potential[opposite] -
        Scalar(0.5) * (potential[face0] + potential[face1])) *
        Scalar(fieldFactor / geometry.oppositeToFaceMidpointDistance);
    if (config.contactElectricFieldFallbackMode == "one_sided")
        return oneSidedField;

    return Scalar(geometry.centroidElectricWeight) * faceNormalField +
           Scalar(1.0 - geometry.centroidElectricWeight) * quasiFermiField;
}

template <typename Scalar>
inline Scalar localAdFieldLimitedMobility(
    Real lowFieldMobility,
    const Scalar& drivingField,
    const FieldMobilityParameters& parameters)
{
    if (lowFieldMobility <= 0.0)
        return Scalar(0.0);
    if (parameters.saturationVelocity <= 0.0 || parameters.beta <= 0.0) {
        throw std::invalid_argument(
            "field saturation velocity and beta must be positive");
    }
    const Scalar field = localAdAbs(drivingField);
    if (localAdValue(field) <= 0.0)
        return Scalar(lowFieldMobility);
    const Scalar ratio =
        Scalar(lowFieldMobility / parameters.saturationVelocity) * field;
    return Scalar(lowFieldMobility) /
        localAdPow(
            Scalar(1.0) + localAdPow(ratio, parameters.beta),
            1.0 / parameters.beta);
}

inline bool localAdUsesFieldMobility(const std::string& model)
{
    return model == "constant_field" ||
           model == "constant_field_lombardi" ||
           model == "caughey_thomas_field" ||
           model == "caughey_thomas_field_surface" ||
           model == "masetti_field" ||
           model == "masetti_field_lombardi";
}

inline bool localAdSupportedMobility(const std::string& model)
{
    return model == "constant" ||
           model == "constant_field" ||
           model == "constant_lombardi" ||
           model == "constant_field_lombardi" ||
           model == "caughey_thomas" ||
           model == "caughey_thomas_field" ||
           model == "caughey_thomas_field_surface" ||
           model == "masetti" ||
           model == "masetti_field" ||
           model == "masetti_field_lombardi";
}

template <typename Scalar>
inline Scalar localAdEndpointAveragedMobility(
    const MobilityModelConfig& mobilityConfig,
    const MobilityModel& mobility,
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const std::vector<Material>& cellMaterials,
    Index cellId,
    Index node0,
    Index node1,
    int localNode0,
    int localNode1,
    CarrierType carrier,
    const Scalar& drivingField,
    const std::array<Scalar, 3>& n,
    const std::array<Scalar, 3>& p)
{
    if (!localAdSupportedMobility(mobilityConfig.model)) {
        throw std::invalid_argument(
            "element-edge local AD does not support mobility model '" +
            mobilityConfig.model + "'");
    }
    const Material& material =
        cellMaterials.at(static_cast<std::size_t>(cellId));
    const Real surfaceField =
        isSurfaceMobilityModel(mobilityConfig) &&
            cellId < mobilityConfig.surface.cellNormalFields.size()
        ? mobilityConfig.surface.cellNormalFields[cellId]
        : 0.0;
    const Real surfaceDistance =
        isSurfaceMobilityModel(mobilityConfig) &&
            cellId < mobilityConfig.surface.cellDistances.size()
        ? mobilityConfig.surface.cellDistances[cellId]
        : 0.0;
    const auto lowFieldAt = [&](Index node, int localNode) {
        const Real mobilityDoping = nodeMobilityDopingConcentration(
            mesh, doping, node, cellId, &mobilityConfig);
        return carrier == CarrierType::Electron
            ? mobility.electronMobility(
                material,
                mobilityDoping,
                static_cast<Real>(
                    localAdValue(n[static_cast<std::size_t>(localNode)])),
                static_cast<Real>(
                    localAdValue(p[static_cast<std::size_t>(localNode)])),
                0.0,
                surfaceField,
                surfaceDistance)
            : mobility.holeMobility(
                material,
                mobilityDoping,
                static_cast<Real>(
                    localAdValue(n[static_cast<std::size_t>(localNode)])),
                static_cast<Real>(
                    localAdValue(p[static_cast<std::size_t>(localNode)])),
                0.0,
                surfaceField,
                surfaceDistance);
    };
    const Real low0 = lowFieldAt(node0, localNode0);
    const Real low1 = lowFieldAt(node1, localNode1);
    if (!localAdUsesFieldMobility(mobilityConfig.model))
        return Scalar(0.5 * (low0 + low1));
    const FieldMobilityParameters& parameters =
        carrier == CarrierType::Electron
        ? mobilityConfig.electronField
        : mobilityConfig.holeField;
    return Scalar(0.5) * (
        localAdFieldLimitedMobility(low0, drivingField, parameters) +
        localAdFieldLimitedMobility(low1, drivingField, parameters));
}

template <typename Scalar>
inline Scalar localAdImpactCoefficient(
    const ImpactIonizationModelConfig& inputConfig,
    CarrierType carrier,
    const Scalar& inputField)
{
    if (inputConfig.model == "none")
        return Scalar(0.0);
    const Scalar field = localAdAbs(inputField);
    const auto fieldValue = localAdValue(field);
    if (inputConfig.model == "selberherr") {
        if (fieldValue < inputConfig.minimumField || fieldValue <= 0.0)
            return Scalar(0.0);
        const Real prefactor = carrier == CarrierType::Electron
            ? inputConfig.electronA : inputConfig.holeA;
        const Real critical = carrier == CarrierType::Electron
            ? inputConfig.electronB : inputConfig.holeB;
        if (prefactor <= 0.0)
            return Scalar(0.0);
        const Scalar rawExponent = -Scalar(critical) / field;
        if (localAdValue(rawExponent) <= -700.0)
            return Scalar(prefactor * std::exp(-700.0));
        if (localAdValue(rawExponent) >= 0.0)
            return Scalar(prefactor);
        return Scalar(prefactor) * localAdExp(rawExponent);
    }
    if (inputConfig.model != "van_overstraeten") {
        throw std::invalid_argument(
            "element-edge local AD does not support impact model '" +
            inputConfig.model + "'");
    }

    const ImpactIonizationModelConfig config =
        applyImpactIonizationParameterSet(inputConfig);
    if (!config.debugRawVanOverstraeten &&
        fieldValue < config.minimumField) {
        return Scalar(0.0);
    }
    if (fieldValue <= 0.0)
        return Scalar(0.0);

    constexpr Real kBoltzmannEvPerK = constants::kb / constants::q;
    const Real referenceArgument =
        config.phononEnergy /
        (2.0 * kBoltzmannEvPerK * config.referenceTemperature_K);
    const Real argument =
        config.phononEnergy /
        (2.0 * kBoltzmannEvPerK * config.temperature_K);
    const Real denominator = std::tanh(argument);
    const Real gamma = std::abs(denominator) <= 0.0
        ? 1.0
        : std::tanh(referenceArgument) / denominator;
    if (gamma <= 0.0)
        return Scalar(0.0);

    const bool lowField = fieldValue < config.switchField;
    const Real prefactor = carrier == CarrierType::Electron
        ? (lowField ? config.electronALow : config.electronAHigh)
        : (lowField ? config.holeALow : config.holeAHigh);
    const Real critical = carrier == CarrierType::Electron
        ? (lowField ? config.electronBLow : config.electronBHigh)
        : (lowField ? config.holeBLow : config.holeBHigh);
    if (prefactor <= 0.0)
        return Scalar(0.0);
    const Scalar rawExponent = -Scalar(critical * gamma) / field;
    if (localAdValue(rawExponent) <= -700.0)
        return Scalar(gamma * prefactor * std::exp(-700.0));
    if (localAdValue(rawExponent) >= 0.0)
        return Scalar(gamma * prefactor);
    return Scalar(gamma * prefactor) * localAdExp(rawExponent);
}

template <typename Scalar>
inline Scalar localAdElectronSgFlux(
    Real ni0,
    Real ni1,
    const Scalar& psi0,
    const Scalar& psi1,
    const Scalar& phin0,
    const Scalar& phin1,
    Real Vt,
    const Scalar& coefficient)
{
    if constexpr (!std::is_same_v<Scalar, Tri3LocalForwardDual>) {
        if (phin0 == phin1)
            return Scalar(0.0);
    }
    const Scalar exponent0 = (psi0 - phin0) / Scalar(Vt);
    const Scalar exponent1 = (psi1 - phin1) / Scalar(Vt);
    const Scalar n0 = Scalar(ni0) * localAdLimitedExp(exponent0);
    const Scalar n1 = Scalar(ni1) * localAdLimitedExp(exponent1);
    if (ni0 <= 0.0 || ni1 <= 0.0) {
        const Scalar electrostaticEta = (psi1 - psi0) / Scalar(Vt);
        return coefficient *
            (localAdBernoulli(-electrostaticEta) * n0 -
             localAdBernoulli(electrostaticEta) * n1);
    }
    const Scalar eta =
        (psi1 - psi0) / Scalar(Vt) + Scalar(std::log(ni1 / ni0));

    const Scalar qfDifference = (phin1 - phin0) / Scalar(Vt);
    if (localAdValue(localAdAbs(qfDifference)) < 50.0 &&
        localAdValue(exponent0) > -500.0 &&
        localAdValue(exponent0) < 500.0 &&
        localAdValue(exponent1) > -500.0 &&
        localAdValue(exponent1) < 500.0) {
        return coefficient * localAdBernoulli(eta) * n1 *
            localAdExpm1(qfDifference);
    }
    return coefficient *
        (localAdBernoulli(-eta) * n0 - localAdBernoulli(eta) * n1);
}

template <typename Scalar>
inline Scalar localAdHoleSgFlux(
    Real ni0,
    Real ni1,
    const Scalar& psi0,
    const Scalar& psi1,
    const Scalar& phip0,
    const Scalar& phip1,
    Real Vt,
    const Scalar& coefficient)
{
    if constexpr (!std::is_same_v<Scalar, Tri3LocalForwardDual>) {
        if (phip0 == phip1)
            return Scalar(0.0);
    }
    const Scalar exponent0 = (phip0 - psi0) / Scalar(Vt);
    const Scalar exponent1 = (phip1 - psi1) / Scalar(Vt);
    const Scalar p0 = Scalar(ni0) * localAdLimitedExp(exponent0);
    const Scalar p1 = Scalar(ni1) * localAdLimitedExp(exponent1);
    if (ni0 <= 0.0 || ni1 <= 0.0) {
        const Scalar electrostaticEta = (psi1 - psi0) / Scalar(Vt);
        return coefficient *
            (localAdBernoulli(electrostaticEta) * p0 -
             localAdBernoulli(-electrostaticEta) * p1);
    }
    const Scalar eta =
        (psi1 - psi0) / Scalar(Vt) + Scalar(std::log(ni0 / ni1));

    const Scalar qfDifference = (phip1 - phip0) / Scalar(Vt);
    if (localAdValue(localAdAbs(qfDifference)) < 50.0 &&
        localAdValue(exponent0) > -500.0 &&
        localAdValue(exponent0) < 500.0 &&
        localAdValue(exponent1) > -500.0 &&
        localAdValue(exponent1) < 500.0) {
        return -coefficient * localAdBernoulli(eta) * p0 *
            localAdExpm1(qfDifference);
    }
    return coefficient *
        (localAdBernoulli(eta) * p0 - localAdBernoulli(-eta) * p1);
}

template <typename Scalar>
inline std::array<Scalar, 2> localAdGssLauxCurrentVector(
    const DeviceMesh& mesh,
    const Cell& cell,
    const std::array<Scalar, 3>& signedEdgeCurrent)
{
    std::array<Scalar, 2> result{Scalar(0.0), Scalar(0.0)};
    for (std::size_t edge = 0; edge < 3; ++edge) {
        std::array<Real, 3> basis{};
        basis[edge] = 1.0;
        const Point2 column =
            gssLauxTri3CurrentVector(mesh, cell, basis);
        result[0] = result[0] + signedEdgeCurrent[edge] * Scalar(column.x());
        result[1] = result[1] + signedEdgeCurrent[edge] * Scalar(column.y());
    }
    return result;
}

template <typename Scalar>
struct ElementEdgeGssLauxLocalSourceIntegrals {
    std::array<Scalar, 3> electron{};
    std::array<Scalar, 3> hole{};
    std::array<Scalar, 3> combined{};
};

template <typename Scalar>
inline Scalar localAdAvalancheMidpointAux2(const Scalar& argument)
{
    return Scalar(1.0) /
        (Scalar(1.0) + localAdLimitedExp(argument));
}

/// Complete nine-column local derivative of the legacy triangle GSS source.
///
/// This is deliberately kept separate from the element-edge GSS/Laux bundle:
/// the legacy source uses cell GradQF for alpha, edge GradQF for |J|, the GSS
/// logistic midpoint density, Genius truncated edge volumes, and a symmetric
/// endpoint split.  All state-dependent factors are carried by Scalar so a
/// Tri3LocalForwardDual evaluation differentiates alpha, density, high-field
/// mobility, edge field, and source partition in one pass.
template <typename Scalar>
inline ElementEdgeGssLauxLocalSourceIntegrals<Scalar>
triangleGssAvalancheSourceIntegralsLocal(
    const ImpactIonizationModelConfig& impactConfig,
    const MobilityModelConfig& mobilityConfig,
    const MobilityModel& mobility,
    const std::vector<Index>& cellEdgeIds,
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const std::vector<Material>& cellMaterials,
    Index cellId,
    const std::array<Scalar, 3>& psi,
    const std::array<Scalar, 3>& phin,
    const std::array<Scalar, 3>& phip,
    const std::array<Scalar, 3>& n,
    const std::array<Scalar, 3>& p,
    Real Vt,
    Real fieldFactor)
{
    if (impactConfig.sourceMappingMode !=
            "triangle_gss_gradqf_truncated" ||
        impactConfig.generation != "current_density" ||
        impactConfig.currentApproximation != "cell_reconstructed" ||
        impactConfig.cellReconstructedMidpointDensity != "gss_logistic" ||
        impactConfig.quasiFermiGradientDiscretization != "cell_gradient" ||
        impactConfig.sourceVolumePolicy != "genius_truncated" ||
        impactConfig.edgeSourcePartition != "symmetric") {
        throw std::invalid_argument(
            "triangle GSS local AD requires the canonical legacy profile");
    }
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "triangle GSS local AD requires Tri3 cells");
    }

    const auto electricGradient = localAdTri3Gradient(mesh, cell, psi);
    const auto electronGradient = localAdTri3Gradient(mesh, cell, phin);
    const auto holeGradient = localAdTri3Gradient(mesh, cell, phip);
    const Scalar electronCellField =
        localAdNorm2(electronGradient[0], electronGradient[1]) *
        Scalar(fieldFactor);
    const Scalar holeCellField =
        localAdNorm2(holeGradient[0], holeGradient[1]) *
        Scalar(fieldFactor);
    Scalar electronImpactField = electronCellField;
    Scalar holeImpactField = holeCellField;
    if (cellUsesContactElectricFieldFallback(
            impactConfig, mesh, cell)) {
        electronImpactField = localAdContactElectricFallbackImpactField(
            impactConfig, mesh, cell, electricGradient,
            electronCellField, psi, fieldFactor);
        holeImpactField = localAdContactElectricFallbackImpactField(
            impactConfig, mesh, cell, electricGradient,
            holeCellField, psi, fieldFactor);
    }
    const Scalar electronAlpha = localAdImpactCoefficient(
        impactConfig, CarrierType::Electron, electronImpactField);
    const Scalar holeAlpha = localAdImpactCoefficient(
        impactConfig, CarrierType::Hole, holeImpactField);

    ElementEdgeGssLauxLocalSourceIntegrals<Scalar> source{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const std::size_t local = static_cast<std::size_t>(localEdge);
        const std::size_t next =
            static_cast<std::size_t>((localEdge + 1) % 3);
        const Index node0 = cell.node_ids[local];
        const Index node1 = cell.node_ids[next];
        const Index edgeId =
            edgeIdForNodePair(mesh, cellEdgeIds, node0, node1);
        if (edgeId >= mesh.numEdges()) {
            throw std::runtime_error(
                "triangle GSS local AD could not map a cell edge");
        }
        const Real edgeLength = mesh.getEdge(edgeId).length;
        if (edgeLength <= 1.0e-30) {
            throw std::invalid_argument(
                "degenerate triangle edge cannot evaluate a GSS source");
        }

        const Scalar electronEdgeField =
            localAdAbs(phin[next] - phin[local]) *
            Scalar(fieldFactor / edgeLength);
        const Scalar holeEdgeField =
            localAdAbs(phip[next] - phip[local]) *
            Scalar(fieldFactor / edgeLength);
        const Scalar electronMobility = localAdEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, cellId,
            node0, node1, localEdge, (localEdge + 1) % 3,
            CarrierType::Electron, electronEdgeField, n, p);
        const Scalar holeMobility = localAdEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, cellId,
            node0, node1, localEdge, (localEdge + 1) % 3,
            CarrierType::Hole, holeEdgeField, n, p);

        const Scalar electronArgument =
            (psi[next] - psi[local]) / Scalar(2.0 * Vt);
        const Scalar holeArgument = -electronArgument;
        const Scalar electronMidpointDensity =
            n[local] * localAdAvalancheMidpointAux2(electronArgument) +
            n[next] * localAdAvalancheMidpointAux2(-electronArgument);
        const Scalar holeMidpointDensity =
            p[local] * localAdAvalancheMidpointAux2(holeArgument) +
            p[next] * localAdAvalancheMidpointAux2(-holeArgument);
        const Scalar partialVolume(
            geniusTri3TruncatedPartialVolumeWithEdge(
                mesh, cell, node0, node1) *
            impactConfig.sourceGeometryScale);
        const Scalar electronIntegral = electronAlpha *
            electronMobility * electronMidpointDensity *
            electronEdgeField * partialVolume;
        const Scalar holeIntegral = holeAlpha *
            holeMobility * holeMidpointDensity *
            holeEdgeField * partialVolume;
        const Scalar combined = electronIntegral + holeIntegral;
        source.electron[local] = source.electron[local] +
            Scalar(0.5) * electronIntegral;
        source.electron[next] = source.electron[next] +
            Scalar(0.5) * electronIntegral;
        source.hole[local] = source.hole[local] +
            Scalar(0.5) * holeIntegral;
        source.hole[next] = source.hole[next] +
            Scalar(0.5) * holeIntegral;
        source.combined[local] = source.combined[local] +
            Scalar(0.5) * combined;
        source.combined[next] = source.combined[next] +
            Scalar(0.5) * combined;
    }
    return source;
}

template <typename Scalar>
inline ElementEdgeGssLauxLocalSourceIntegrals<Scalar>
elementEdgeGssLauxAvalancheSourceIntegralsLocal(
    const ImpactIonizationModelConfig& impactConfig,
    const MobilityModelConfig& mobilityConfig,
    const MobilityModel& mobility,
    const std::vector<Index>& cellEdgeIds,
    const DeviceMesh& mesh,
    const DopingModel& doping,
    const std::vector<Material>& cellMaterials,
    Index cellId,
    const std::array<Scalar, 3>& psi,
    const std::array<Scalar, 3>& phin,
    const std::array<Scalar, 3>& phip,
    const std::array<Scalar, 3>& n,
    const std::array<Scalar, 3>& p,
    const std::array<Real, 3>& ni,
    Real Vt,
    Real fieldFactor)
{
    if (!usesElementEdgeGssLauxAvalancheSource(impactConfig)) {
        throw std::invalid_argument(
            "element-edge local AD requires its canonical configuration");
    }
    const Cell& cell = mesh.getCell(cellId);
    if (cell.type != CellType::Tri3 || cell.node_ids.size() != 3) {
        throw std::invalid_argument(
            "element-edge local AD requires Tri3 cells");
    }

    const auto electricGradient = localAdTri3Gradient(mesh, cell, psi);
    const auto electronGradient = localAdTri3Gradient(mesh, cell, phin);
    const auto holeGradient = localAdTri3Gradient(mesh, cell, phip);
    const Scalar electronQfField =
        localAdNorm2(electronGradient[0], electronGradient[1]) *
        Scalar(fieldFactor);
    const Scalar holeQfField =
        localAdNorm2(holeGradient[0], holeGradient[1]) *
        Scalar(fieldFactor);
    Scalar electronImpactField = electronQfField;
    Scalar holeImpactField = holeQfField;
    if (impactConfig.drivingForce == "electric_field") {
        electronImpactField = localAdNorm2(
            electricGradient[0], electricGradient[1]) * Scalar(fieldFactor);
        holeImpactField = electronImpactField;
    } else if (cellUsesContactElectricFieldFallback(
                   impactConfig, mesh, cell)) {
        electronImpactField = localAdContactElectricFallbackImpactField(
            impactConfig, mesh, cell, electricGradient, electronQfField,
            psi, fieldFactor);
        holeImpactField = localAdContactElectricFallbackImpactField(
            impactConfig, mesh, cell, electricGradient, holeQfField,
            psi, fieldFactor);
    }
    const Scalar electronAlpha = localAdImpactCoefficient(
        impactConfig, CarrierType::Electron, electronImpactField);
    const Scalar holeAlpha = localAdImpactCoefficient(
        impactConfig, CarrierType::Hole, holeImpactField);

    std::array<Scalar, 3> electronFlux{};
    std::array<Scalar, 3> holeFlux{};
    for (int localEdge = 0; localEdge < 3; ++localEdge) {
        const std::size_t local = static_cast<std::size_t>(localEdge);
        const std::size_t next =
            static_cast<std::size_t>((localEdge + 1) % 3);
        const Index node0 = cell.node_ids[local];
        const Index node1 = cell.node_ids[next];
        const Index edgeId =
            edgeIdForNodePair(mesh, cellEdgeIds, node0, node1);
        if (edgeId >= mesh.numEdges()) {
            throw std::runtime_error(
                "element-edge local AD could not map a cell edge");
        }
        const Real edgeLength =
            (meshPoint(mesh, node1) - meshPoint(mesh, node0)).norm();
        if (edgeLength <= 1.0e-30) {
            throw std::invalid_argument(
                "degenerate triangle edge cannot evaluate an SG current");
        }
        const Scalar electricEdgeField =
            localAdAbs(psi[next] - psi[local]) *
            Scalar(fieldFactor / edgeLength);
        const Scalar electronEdgeField =
            localAdAbs(phin[next] - phin[local]) *
            Scalar(fieldFactor / edgeLength);
        const Scalar holeEdgeField =
            localAdAbs(phip[next] - phip[local]) *
            Scalar(fieldFactor / edgeLength);
        const bool qfMobility =
            mobilityConfig.highFieldDrivingForce == "quasi_fermi_gradient";
        const Scalar electronMobility = localAdEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, cellId,
            node0, node1, localEdge, (localEdge + 1) % 3,
            CarrierType::Electron,
            qfMobility ? electronEdgeField : electricEdgeField,
            n, p);
        const Scalar holeMobility = localAdEndpointAveragedMobility(
            mobilityConfig, mobility, mesh, doping, cellMaterials, cellId,
            node0, node1, localEdge, (localEdge + 1) % 3,
            CarrierType::Hole,
            qfMobility ? holeEdgeField : electricEdgeField,
            n, p);
        electronFlux[local] =
            localAdValue(electronMobility) > 0.0
            ? localAdElectronSgFlux(
                ni[local], ni[next],
                psi[local], psi[next], phin[local], phin[next], Vt,
                electronMobility * Scalar(Vt * fieldFactor / edgeLength))
            : Scalar(0.0);
        holeFlux[local] =
            localAdValue(holeMobility) > 0.0
            ? localAdHoleSgFlux(
                ni[local], ni[next],
                psi[local], psi[next], phip[local], phip[next], Vt,
                holeMobility * Scalar(Vt * fieldFactor / edgeLength))
            : Scalar(0.0);
    }

    const auto electronCurrent =
        localAdGssLauxCurrentVector(mesh, cell, electronFlux);
    const auto holeCurrent =
        localAdGssLauxCurrentVector(mesh, cell, holeFlux);
    const Scalar electronGeneration =
        electronAlpha *
        localAdNorm2(electronCurrent[0], electronCurrent[1]);
    const Scalar holeGeneration =
        holeAlpha * localAdNorm2(holeCurrent[0], holeCurrent[1]);
    const auto vertexMeasures = tri3ElementVertexBoxMeasures(mesh, cell);
    ElementEdgeGssLauxLocalSourceIntegrals<Scalar> source{};
    for (std::size_t localNode = 0; localNode < 3; ++localNode) {
        const Scalar measure =
            Scalar(vertexMeasures[localNode] *
                   impactConfig.sourceGeometryScale);
        source.electron[localNode] = electronGeneration * measure;
        source.hole[localNode] = holeGeneration * measure;
        source.combined[localNode] =
            source.electron[localNode] + source.hole[localNode];
    }
    return source;
}
