#include "vela/solver/ExtendedBoltzmannColdStart.h"

#include "vela/solver/ConstraintElimination.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace vela {

ExtendedOhmicTargets buildExtendedZeroBiasOhmicTargets(
    const DeviceMesh& mesh, const DopingModel& doping,
    const ExtendedBoltzmannMaterials& materials,
    const std::vector<ContactBoundarySpec>& contacts,
    long double thermalVoltage_V, bool dominantSignedContactMean)
{
    const auto count = mesh.numNodes();
    const auto dimension = static_cast<Eigen::Index>(count);
    if (!count || doping.numNodes() != count || materials.intrinsic.size() != dimension ||
        materials.intrinsicPotential_V.size() != dimension ||
        !std::isfinite(thermalVoltage_V) || thermalVoltage_V <= 0)
        throw std::invalid_argument("Inconsistent extended Ohmic inputs");

    ExtendedOhmicTargets targets{
        std::vector<std::optional<long double>>(count),
        ExtendedVector::Zero(count), ExtendedVector::Zero(count)};
    std::unordered_set<std::string> bound;
    for (const auto& terminal : mesh.contacts()) {
        const auto spec = std::ranges::find(contacts, terminal.name, &ContactBoundarySpec::name);
        if (spec == contacts.end())
            continue;
        if (terminal.node_ids.empty() || !std::isfinite(spec->bias) ||
            spec->type != ContactType::Ohmic || spec->bias != 0 ||
            (spec->decimalInitialBias_V && *spec->decimalInitialBias_V != 0) ||
            !bound.insert(terminal.name).second)
            throw std::invalid_argument("Extended cold start requires distinct zero-bias Ohmic contacts");

        long double mean = 0;
        for (const auto node : terminal.node_ids)
            mean += doping.netDopingForAssembly(node);
        mean /= terminal.node_ids.size();
        for (const auto node : terminal.node_ids) {
            auto& potential = targets.potential.at(node);
            if (potential)
                throw std::invalid_argument("Multiple extended Ohmic targets on one node");
            auto net = doping.netDopingForAssembly(node);
            if (dominantSignedContactMean && mean != 0 &&
                (net == 0 || std::signbit(net) != std::signbit(mean)))
                net = mean;
            const auto target = materials.intrinsic[node] > 0
                ? extendedBoltzmannOhmicTarget(net, materials.intrinsic[node],
                    materials.intrinsicPotential_V[node], thermalVoltage_V, 0, thermalVoltage_V, 0)
                : ExtendedOhmicPopulation{};
            potential = target.potential;
            targets.electrons[node] = target.electrons;
            targets.holes[node] = target.holes;
        }
    }
    if (contacts.empty() || bound.size() != contacts.size())
        throw std::invalid_argument("Extended Ohmic contacts do not match the mesh");
    return targets;
}

ExtendedBoltzmannColdStartResult solveExtendedBoltzmannColdStart(
    BasicPoissonAssembler<long double>& poisson, const DopingModel& doping,
    const ExtendedBoltzmannMaterials& materials, const ExtendedOhmicTargets& targets,
    long double thermalVoltage_V, long double concentrationScale_m3,
    int maxUpdates, long double absoluteResidualTolerance,
    const ExtendedAveragePoissonGeometry* geometry)
{
    if constexpr (compiledArithmeticPrecision != ArithmeticPrecision::Extended80) {
        throw std::logic_error("Extended cold start requires an extended80 build");
    } else {
        const auto count = static_cast<Eigen::Index>(doping.numNodes());
        if (count <= 0 || materials.intrinsic.size() != count ||
            materials.electronIntrinsic.size() != count || materials.holeIntrinsic.size() != count ||
            materials.intrinsicPotential_V.size() != count || std::ssize(targets.potential) != count ||
            targets.electrons.size() != count || targets.holes.size() != count || maxUpdates < 0 ||
            !std::isfinite(absoluteResidualTolerance) || absoluteResidualTolerance <= 0 ||
            !std::isfinite(thermalVoltage_V) || thermalVoltage_V <= 0 ||
            thermalVoltage_V != poisson.potentialScale_V() ||
            !std::isfinite(concentrationScale_m3) || concentrationScale_m3 <= 0)
            throw std::invalid_argument("Inconsistent extended Poisson cold-start inputs");

        std::unordered_map<Index, long double> boundary;
        std::vector<int> prescribed;
        ExtendedVector potential = ExtendedVector::Zero(count);
        for (Index node = 0; node < doping.numNodes(); ++node) {
            if (targets.potential[node]) {
                potential[node] = *targets.potential[node];
                boundary.emplace(node, potential[node]);
                prescribed.push_back(static_cast<int>(node));
            } else if (materials.intrinsic[node] > 0) {
                potential[node] = extendedBoltzmannOhmicTarget(doping.netDopingForAssembly(node),
                    materials.intrinsic[node], materials.intrinsicPotential_V[node],
                    thermalVoltage_V, 0, thermalVoltage_V, 0).potential;
            }
        }
        if (prescribed.empty())
            throw std::invalid_argument("Extended Poisson cold start requires a potential reference");
        const ExtendedVector qf = ExtendedVector::Zero(count);
        const auto assemble = [&](const ExtendedVector& state) {
            return poisson.linearizeBoltzmannEquilibrium(state, qf, qf,
                materials.electronIntrinsic, materials.holeIntrinsic, thermalVoltage_V, boundary,
                ArithmeticPoissonAssembler::DirichletCoordinate::SolverPotential, geometry);
        };
        const auto norm = [](const auto& frame) {
            if (!frame.rhs.allFinite() || !frame.electrons.allFinite() || !frame.holes.allFinite())
                throw std::runtime_error("Nonfinite extended Poisson equilibrium");
            return frame.rhs.cwiseAbs().maxCoeff();
        };
        auto current = assemble(potential);
        auto residual = norm(current);
        int updates = 0;
        while (updates < maxUpdates && residual > absoluteResidualTolerance) {
            const auto reduced = reduceLinearConstraints(current.jacobian, current.rhs, prescribed);
            const auto direction = reduced.recover(solveExtendedSparseLu(reduced.matrix, reduced.rhs));
            bool progressed = false;
            for (int backtrack = 0; backtrack <= 10; ++backtrack) {
                ExtendedVector trial = potential + std::ldexp(1.L, -backtrack) * direction;
                auto next = assemble(trial);
                const auto nextResidual = norm(next);
                if (nextResidual >= residual)
                    continue;
                potential = std::move(trial);
                current = std::move(next);
                residual = nextResidual;
                ++updates;
                progressed = true;
                break;
            }
            if (!progressed)
                break;
        }
        ExtendedPopulationState state;
        state.potentialScale_V = thermalVoltage_V;
        state.concentrationScale_m3 = concentrationScale_m3;
        state.coordinates.resize(3 * count);
        state.coordinates << potential, current.electrons, current.holes;
        return {std::move(state), residual, updates, residual <= absoluteResidualTolerance};
    }
}

} // namespace vela
