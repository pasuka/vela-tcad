#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"

#include <map>
#include <vector>

namespace vela {

struct RecoveredField2 {
    Point2 vector = Point2::Zero();
    Real magnitude = 0.0;
    bool valid = false;
};

struct CellField2 {
    Index cellId = 0;
    Index regionId = 0;
    Point2 vector = Point2::Zero();
    Real magnitude = 0.0;
    Real area = 0.0;
    bool valid = false;
};

struct EdgeField2 {
    Index edgeId = 0;
    Index node0 = 0;
    Index node1 = 0;
    Point2 vector = Point2::Zero();
    Real projected = 0.0;
    Real magnitude = 0.0;
    bool valid = false;
};

struct NodeField2 {
    Index nodeId = 0;
    Index regionId = 0;
    Point2 vector = Point2::Zero();
    Real magnitude = 0.0;
    bool valid = false;
    std::map<Index, RecoveredField2> regionSamples;
};

enum class ElectricFieldLeastSquaresWeight {
    InverseDistance,
    InverseDistanceSquared,
};

/**
 * @brief Return the maximum edge-projected electric field magnitude from nodal potential.
 *
 * This is an engineering diagnostic computed as max(|delta psi| / edge length)
 * over valid mesh edges. It is intentionally edge-based and should not be
 * interpreted as a higher-order reconstructed electric-field solution.
 * Degenerate edges with non-positive or non-finite length are ignored.
 */
Real maxEdgeElectricFieldMagnitude(const DeviceMesh& mesh, const VectorXd& potential_V);

std::vector<CellField2> computeCellElectricField(const DeviceMesh& mesh,
                                                 const VectorXd& potential_V);

std::vector<CellField2> computeCellGradElectronQuasiFermi(const DeviceMesh& mesh,
                                                          const VectorXd& electronQf_V);

std::vector<CellField2> computeCellGradHoleQuasiFermi(const DeviceMesh& mesh,
                                                      const VectorXd& holeQf_V);

std::vector<EdgeField2> computeEdgeElectricField(const DeviceMesh& mesh,
                                                 const VectorXd& potential_V);

std::vector<NodeField2> computeNodeElectricFieldAreaAverage(const DeviceMesh& mesh,
                                                            const VectorXd& potential_V);

std::vector<NodeField2> computeNodeElectricFieldLeastSquares(
    const DeviceMesh& mesh,
    const VectorXd& potential_V,
    ElectricFieldLeastSquaresWeight weight);

std::vector<NodeField2> computeNodeElectricFieldSPR(const DeviceMesh& mesh,
                                                    const VectorXd& potential_V);

} // namespace vela
