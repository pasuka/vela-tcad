#include "vela/io/DDSolutionState.h"
#include "vela/numerics/ExpansionProductSum.h"
#include <tuple>
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
    if(solution.stateVariables=="potential_density") {
        if(solution.packedState.size() || solution.packedLow.size() || solution.phinIncrement.size() || solution.phipIncrement.size())
            throw std::invalid_argument("population archive cannot carry QF packed coordinates");
        s.metadata["state_variables"]="potential_density";
        s.metadata["authoritative_fields"]={"psi","electrons_m3","holes_m3"};
    }
    if (s.metadata.at("mode") != "dd") throw std::invalid_argument("DDSolution requires dd archive mode");
    auto field = [&](const char* name, const VectorXd& values) {
        if (values.size() != solution.psi.size()) throw std::invalid_argument("Incomplete DD field");
        s.fields[name] = std::vector<double>(values.data(), values.data()+values.size());
    };
    field("psi", solution.psi); field("phin", solution.phin); field("phip", solution.phip);
    field("electrons_m3", solution.n); field("holes_m3", solution.p);
    const auto& units = scaling.unitSystem();
    if(solution.stateVariables=="potential_density") {
        field("population_electrons_internal",solution.n);field("population_holes_internal",solution.p);
        s.metadata["population_coordinate_schema"]="vela.population-coordinates.v1";
        s.metadata["population_concentration_m3_per_internal"]=units.internalConcentrationToM3(1.);
        s.metadata["authoritative_fields"]={"psi","population_electrons_internal","population_holes_internal"};
        if(solution.populationPsi.size()) {
            if(!solution.hasConsistentPopulationPotential())
                throw std::invalid_argument("stale population potential coordinates");
            const double origin=s.metadata.at("potential_origin_V");
            if(solution.populationPotentialOrigin_V && *solution.populationPotentialOrigin_V!=origin)
                throw std::invalid_argument("population potential reference changed without translation");
            field("population_psi_normalized",solution.populationPsi);
            s.metadata["population_potential_schema"]="vela.population-potential.v1";
            s.metadata["population_potential_scale_V"]=solution.populationPotentialScale_V;
            s.metadata["population_potential_origin_V"]=origin;
            s.metadata["population_potential_unit"]="V";
            s.metadata["authoritative_fields"]={"population_psi_normalized","population_electrons_internal","population_holes_internal"};
        }
    }
    if (solution.populationLow.size()) {
        const int count = solution.psi.size();
        if (solution.stateVariables != "potential_density" ||
            solution.populationPsi.size() || solution.populationLow.size() != 3*count)
            throw std::invalid_argument("Incomplete compensated population state");
        field("population_psi_low_V",solution.populationLow.head(count));
        field("population_electrons_low_internal",solution.populationLow.segment(count,count));
        field("population_holes_low_internal",solution.populationLow.tail(count));
    }
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
            if (solution.stateVariables!="potential_density" && v != 0. && std::abs(v) < std::numeric_limits<double>::min()) v = 0.;
    }
    validateStateArchive(s);
    return s;
}

DDSolution restoreDDSolution(const StateArchive& archive, UnitScalingConfig scaling) {
    validateStateArchive(archive);
    if (archive.metadata.at("mode") != "dd") throw std::invalid_argument("DDSolution requires dd archive mode");
    DDSolution s;
    s.stateVariables=archive.metadata.value("state_variables",std::string{"potential_quasi_fermi"});
    if(s.stateVariables!="potential_density" && s.stateVariables!="potential_quasi_fermi")
        throw std::invalid_argument("unknown archived state_variables");
    if(s.stateVariables=="potential_density" && (archive.fields.contains("packed_psi") || archive.fields.contains("electron_qf_increment_V")))
        throw std::invalid_argument("population archive has incompatible QF coordinates");
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
    if(s.stateVariables=="potential_density") {
        s.populationRestoreMode="legacy_si";
        if(archive.fields.contains("population_electrons_internal")) {
            const double stored=archive.metadata.at("population_concentration_m3_per_internal");
            const double current=scaling.unitSystem().internalConcentrationToM3(1.);
            if(stored==current) {
                field("population_electrons_internal",s.n);field("population_holes_internal",s.p);
                s.populationRestoreMode="exact_internal";
            } else s.populationRestoreMode="converted_units";
        }
    }
    if(archive.fields.contains("population_psi_normalized")) {
        field("population_psi_normalized",s.populationPsi);
        s.populationPotentialScale_V=archive.metadata.at("population_potential_scale_V");
        s.populationPotentialOrigin_V=archive.metadata.at("population_potential_origin_V").get<double>();
        if(!s.hasConsistentPopulationPotential())throw std::invalid_argument("inconsistent population potential coordinates");
    }
    if (archive.fields.contains("population_psi_low_V")) {
        const int count = archive.nodeCount;
        const double stored = archive.metadata.at("population_concentration_m3_per_internal");
        const double current = scaling.unitSystem().internalConcentrationToM3(1.);
        s.populationLow.resize(3*count);
        VectorXd values;
        field("population_psi_low_V",values);
        s.populationLow.head(count) = values;
        for (const auto& [name,density,offset] :
             {std::tuple{"population_electrons_internal",&s.n,count},
              std::tuple{"population_holes_internal",&s.p,2*count}}) {
            const char* lowField = offset == count ? "population_electrons_low_internal" :
                "population_holes_low_internal";
            field(lowField,values);
            if (stored == current) {
                s.populationLow.segment(offset,count) = values;
                continue;
            }
            const auto& high = archive.fields.at(name);
            const double ratio = stored/current;
            for (int i = 0; i < count; ++i) {
                experimental::ExpansionProductSum product;
                product.addProduct({high[i],values[i]},ratio);
                const auto converted = product.parts();
                (*density)[i] = converted.high;
                s.populationLow[offset+i] = converted.low;
            }
        }
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
