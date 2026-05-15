#include "vela/solver/SolutionValidation.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace vela {

namespace {

std::string formatDiagnosticReal(Real value)
{
    std::ostringstream oss;
    oss.precision(17);
    oss << value;
    return oss.str();
}

void addDiagnostic(DDSolutionValidationResult& result, std::string message)
{
    result.valid = false;
    result.diagnostics.push_back(std::move(message));
}

DDSolutionFieldStats computeStats(const VectorXd& values)
{
    DDSolutionFieldStats stats;
    if (values.size() == 0)
        return stats;

    stats.min = values(0);
    stats.max = values(0);
    for (int i = 1; i < values.size(); ++i) {
        stats.min = std::min(stats.min, values(i));
        stats.max = std::max(stats.max, values(i));
    }
    return stats;
}

void validateFiniteField(DDSolutionValidationResult& result,
                         const VectorXd& values,
                         const char* fieldName)
{
    for (int i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values(i))) {
            addDiagnostic(result,
                          std::string(fieldName) + "[" + std::to_string(i) + "] is not finite");
            return;
        }
    }
}

void validateCarrierFloor(DDSolutionValidationResult& result,
                          const VectorXd& values,
                          const char* fieldName,
                          Real carrierFloor)
{
    const Real lowerBound = -std::abs(carrierFloor);
    for (int i = 0; i < values.size(); ++i) {
        if (values(i) < lowerBound) {
            addDiagnostic(result,
                          std::string(fieldName) + "[" + std::to_string(i) + "]=" +
                          formatDiagnosticReal(values(i)) + " is below -carrier_floor=" +
                          formatDiagnosticReal(lowerBound));
            return;
        }
    }
}

void validateMinimumCarrier(DDSolutionValidationResult& result,
                            const VectorXd& values,
                            const char* fieldName,
                            Real minimumCarrierDensity,
                            Real carrierFloor)
{
    const Real lowerBound = minimumCarrierDensity - std::abs(carrierFloor);
    for (int i = 0; i < values.size(); ++i) {
        if (values(i) < lowerBound) {
            addDiagnostic(result,
                          std::string(fieldName) + "[" + std::to_string(i) + "]=" +
                          formatDiagnosticReal(values(i)) + " is below minimum=" +
                          formatDiagnosticReal(minimumCarrierDensity));
            return;
        }
    }
}

bool nearlyEqual(Real actual, Real expected, Real absTol, Real relTol)
{
    return std::abs(actual - expected) <=
           std::max(absTol, relTol * std::max(std::abs(actual), std::abs(expected)));
}

void validateContactQuasiFermi(DDSolutionValidationResult& result,
                               const DDSolution& sol,
                               const DeviceMesh& mesh,
                               const std::unordered_map<std::string, Real>& contactBiases,
                               const DDSolutionValidationOptions& options)
{
    for (const Contact& contact : mesh.contacts()) {
        const auto it = contactBiases.find(contact.name);
        if (it == contactBiases.end())
            continue;

        const Real bias = it->second;
        for (Index node : contact.node_ids) {
            if (node >= mesh.numNodes()) {
                addDiagnostic(result,
                              "contact '" + contact.name + "' references invalid node " +
                              std::to_string(node));
                continue;
            }
            const int i = static_cast<int>(node);
            if (!nearlyEqual(sol.phin(i), bias,
                             options.contactPotentialAbsTolerance,
                             options.contactPotentialRelTolerance)) {
                addDiagnostic(result,
                              "contact '" + contact.name + "' node " + std::to_string(node) +
                              " phin=" + formatDiagnosticReal(sol.phin(i)) +
                              " does not match bias=" + formatDiagnosticReal(bias));
                return;
            }
            if (!nearlyEqual(sol.phip(i), bias,
                             options.contactPotentialAbsTolerance,
                             options.contactPotentialRelTolerance)) {
                addDiagnostic(result,
                              "contact '" + contact.name + "' node " + std::to_string(node) +
                              " phip=" + formatDiagnosticReal(sol.phip(i)) +
                              " does not match bias=" + formatDiagnosticReal(bias));
                return;
            }
        }
    }
}

} // namespace

std::string DDSolutionValidationResult::diagnosticsString() const
{
    std::ostringstream oss;
    for (std::size_t i = 0; i < diagnostics.size(); ++i) {
        if (i > 0)
            oss << ';';
        oss << diagnostics[i];
    }
    return oss.str();
}

DDSolutionValidationResult validateDDSolution(
    const DDSolution& sol,
    const DeviceMesh& mesh,
    const std::unordered_map<std::string, Real>& contactBiases,
    const DDSolutionValidationOptions& options)
{
    DDSolutionValidationResult result;
    const int expectedSize = static_cast<int>(mesh.numNodes());
    const auto validateSize = [&](const VectorXd& values, const char* fieldName) {
        if (values.size() != expectedSize) {
            addDiagnostic(result,
                          std::string(fieldName) + " size=" + std::to_string(values.size()) +
                          " does not match mesh node count=" + std::to_string(expectedSize));
        }
    };

    validateSize(sol.psi, "psi");
    validateSize(sol.phin, "phin");
    validateSize(sol.phip, "phip");
    validateSize(sol.n, "n");
    validateSize(sol.p, "p");
    if (!result.valid)
        return result;

    result.psi = computeStats(sol.psi);
    result.phin = computeStats(sol.phin);
    result.phip = computeStats(sol.phip);
    result.n = computeStats(sol.n);
    result.p = computeStats(sol.p);

    validateFiniteField(result, sol.psi, "psi");
    validateFiniteField(result, sol.phin, "phin");
    validateFiniteField(result, sol.phip, "phip");
    validateFiniteField(result, sol.n, "n");
    validateFiniteField(result, sol.p, "p");
    if (!result.valid)
        return result;

    validateCarrierFloor(result, sol.n, "n", options.carrierFloor);
    validateCarrierFloor(result, sol.p, "p", options.carrierFloor);

    if (options.enforceMinimumCarrierDensity) {
        validateMinimumCarrier(result, sol.n, "n", options.minimumCarrierDensity, options.carrierFloor);
        validateMinimumCarrier(result, sol.p, "p", options.minimumCarrierDensity, options.carrierFloor);
    }

    if (options.checkContactQuasiFermiBias)
        validateContactQuasiFermi(result, sol, mesh, contactBiases, options);

    return result;
}

} // namespace vela
