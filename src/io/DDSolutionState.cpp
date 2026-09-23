#include "vela/io/DDSolutionState.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vela {
namespace { thread_local const DDStateArchiveScope* activeArchive = nullptr; }
DDStateArchiveScope::DDStateArchiveScope(nlohmann::json value)
    : metadata(std::move(value)), previous_(activeArchive) { activeArchive = this; }
DDStateArchiveScope::~DDStateArchiveScope() { activeArchive = previous_; }
const nlohmann::json* activeDDStateArchiveMetadata() {
    return activeArchive ? &activeArchive->metadata : nullptr;
}
StateArchive archiveDDSolution(const DDSolution& solution, nlohmann::json metadata,
                               UnitScalingConfig scaling) {
    StateArchive s;
    s.nodeCount = solution.psi.size();
    s.metadata = std::move(metadata);
    if (s.metadata.at("mode") != "dd") throw std::invalid_argument("DDSolution requires dd archive mode");
    auto field = [&](const char* name, const VectorXd& values) {
        if (values.size() != solution.psi.size()) throw std::invalid_argument("Incomplete DD field");
        s.fields[name] = std::vector<double>(values.data(), values.data()+values.size());
    };
    field("psi", solution.psi); field("phin", solution.phin); field("phip", solution.phip);
    field("electrons_m3", solution.n); field("holes_m3", solution.p);
    const auto& units = scaling.unitSystem();
    for (const char* name : {"electrons_m3", "holes_m3"})
        for (auto& v : s.fields.at(name)) v = units.internalConcentrationToM3(v);
    if (solution.electronQuantumPotential.size()) field("electron_quantum_potential_V", solution.electronQuantumPotential);
    if (solution.electronQuantumPotentialLike.size()) field("electron_quantum_potential_like_V", solution.electronQuantumPotentialLike);
    if (solution.phinIncrement.size() || solution.phipIncrement.size()) {
        field("electron_qf_increment_V", solution.phinIncrement);
        field("hole_qf_increment_V", solution.phipIncrement);
        if (solution.electronQfReference.size()) field("electron_qf_reference_V", solution.electronQfReference);
        else s.fields["electron_qf_reference_V"] = std::vector<double>(s.nodeCount, solution.electronQfReference_V);
        if (solution.holeQfReference.size()) field("hole_qf_reference_V", solution.holeQfReference);
        else s.fields["hole_qf_reference_V"] = std::vector<double>(s.nodeCount, solution.holeQfReference_V);
    } else if (solution.electronQfReference.size() || solution.holeQfReference.size()) {
        throw std::invalid_argument("DD reference without increment");
    }
    // Preserve the pre-existing solver restart writer's subnormal policy.
    // The general archive codec itself is lossless and never normalizes values.
    for (auto& [name, values] : s.fields)
        for (auto& v : values)
            if (v != 0. && std::abs(v) < std::numeric_limits<double>::min()) v = 0.;
    validateStateArchive(s);
    return s;
}

DDSolution restoreDDSolution(const StateArchive& archive, UnitScalingConfig scaling) {
    validateStateArchive(archive);
    if (archive.metadata.at("mode") != "dd") throw std::invalid_argument("DDSolution requires dd archive mode");
    DDSolution s;
    auto field = [&](const char* name, VectorXd& target) {
        const auto& values = archive.fields.at(name);
        target = Eigen::Map<const VectorXd>(values.data(), static_cast<Eigen::Index>(values.size()));
    };
    field("psi", s.psi); field("phin", s.phin); field("phip", s.phip);
    field("electrons_m3", s.n); field("holes_m3", s.p);
    for (auto* density : {&s.n, &s.p}) for (auto& v : *density) {
        v = scaling.unitSystem().m3ToInternalConcentration(v);
        if (!std::isfinite(v)) throw std::invalid_argument("DD density conversion overflow");
    }
    s.electronQuantumPotential = VectorXd::Zero(s.psi.size());
    if (archive.fields.contains("electron_quantum_potential_V")) field("electron_quantum_potential_V", s.electronQuantumPotential);
    if (archive.fields.contains("electron_quantum_potential_like_V")) field("electron_quantum_potential_like_V", s.electronQuantumPotentialLike);
    if (archive.fields.contains("electron_qf_increment_V")) {
        field("electron_qf_increment_V", s.phinIncrement);
        field("hole_qf_increment_V", s.phipIncrement);
        field("electron_qf_reference_V", s.electronQfReference);
        field("hole_qf_reference_V", s.holeQfReference);
        s.electronQfReference_V = s.electronQfReference(0);
        s.holeQfReference_V = s.holeQfReference(0);
    }
    s.converged = true; s.iters = 0;
    return s;
}
} // namespace vela
