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
    if (solution.packedLow.size() && !solution.hasConsistentSplitPackedState())
        throw std::invalid_argument("Cannot archive stale or incomplete split DD coordinates");
    if (solution.hasConsistentPackedState()) {
        const int count = static_cast<int>(s.nodeCount);
        const char* names[] = {"packed_psi", "packed_electron_qf_increment", "packed_hole_qf_increment"};
        for (int block = 0; block < 3; ++block) {
            field(names[block], solution.packedState.segment(block * count, count));
            if (solution.packedLow.size()) {
                const std::string low = std::string(names[block]) + "_low";
                field(low.c_str(), solution.packedLow.segment(block * count, count));
            }
        }
        s.metadata["packed_potential_scale_V"] = solution.packedPotentialScale_V;
        if (solution.packedLow.size()) {
            s.metadata["split_state_schema"] = "vela.split-dd-state.v1";
            s.metadata["split_mesh_fingerprint"] = solution.packedMeshFingerprint;
        }
    }
    // Preserve the pre-existing solver restart writer's subnormal policy.
    // The general archive codec itself is lossless and never normalizes values.
    for (auto& [name, values] : s.fields) {
        // Low coordinate tails are solver state, including subnormal tails.
        // They must not inherit physical-output normalization.
        if (name.starts_with("packed_") && name.ends_with("_low")) continue;
        for (auto& v : values)
            if (v != 0. && std::abs(v) < std::numeric_limits<double>::min()) v = 0.;
    }
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
    if (archive.fields.contains("packed_psi")) {
        const int count = static_cast<int>(archive.nodeCount);
        s.packedState.resize(3 * count);
        const bool low = archive.fields.contains("packed_psi_low");
        if (low) s.packedLow.resize(3 * count);
        const char* names[] = {"packed_psi", "packed_electron_qf_increment", "packed_hole_qf_increment"};
        for (int block = 0; block < 3; ++block) {
            VectorXd values;
            field(names[block], values);
            s.packedState.segment(block * count, count) = values;
            if (low) {
                const std::string name = std::string(names[block]) + "_low";
                field(name.c_str(), values);
                s.packedLow.segment(block * count, count) = values;
            }
        }
        s.packedPotentialScale_V = archive.metadata.at("packed_potential_scale_V");
        if (low) s.packedMeshFingerprint = archive.metadata.at("split_mesh_fingerprint");
        if (!s.hasConsistentPackedState())
            throw std::invalid_argument("Archived DD coordinates disagree with physical state");
    }
    s.converged = true; s.iters = 0;
    return s;
}
} // namespace vela
