#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"
#include <string>
#include <vector>

namespace vela {

/**
 * @brief Top-level driver for a linear electrostatic Poisson simulation.
 *
 * Reads a JSON configuration file, builds the mesh, assigns doping,
 * assembles the Poisson equation, applies Dirichlet contacts, solves the
 * linear system, and writes the electrostatic potential to a VTK file.
 *
 * JSON configuration schema
 * -------------------------
 * @code
 * {
 *   "mesh_file":   "path/to/mesh.json",
 *   "output_vtk":  "path/to/output.vtk",
 *   "doping": [
 *     { "region": "n_region", "donors": 1e23, "acceptors": 0.0 },
 *     { "region": "p_region", "donors": 0.0,  "acceptors": 1e23 }
 *   ],
 *   "regions": [
 *     { "name": "oxide", "fixed_charge_m3": 1e21 }
 *   ],
 *   "interfaces": [
 *     { "regions": ["silicon", "oxide"], "sheet_charge_m2": 1e15 }
 *   ],
 *   "contacts": [
 *     { "name": "cathode", "bias": 0.0, "flatband_voltage": 0.1 },
 *     { "name": "anode",   "bias": 0.0, "work_function_eV": 0.0 }
 *   ]
 * }
 * @endcode
 *
 * Charge and contact conventions
 * ------------------------------
 * fixed_charge_m3 and sheet_charge_m2 are signed elementary-charge number
 * densities. The assembler multiplies them by q. Sheet charge on a shared-node
 * interface is distributed edge-by-edge: q * sheet_charge_m2 * edge_length / 2
 * is added to each endpoint of every edge whose adjacent cells match the
 * configured region pair. Contact flatband_voltage or work_function_eV shifts
 * the Dirichlet electrostatic potential as psi = bias - offset.
 *
 * Path resolution
 * ---------------
 * Relative paths in "mesh_file" and "output_vtk" are resolved relative to
 * the directory containing the config JSON.
 */
struct PoissonResult {
    DeviceMesh mesh;
    VectorXd potential;
    std::vector<Real> netDoping;
};

class PoissonSimulation {
public:
    /**
     * @brief Run the simulation from a config JSON file.
     *
     * @param configFile  Path to the JSON configuration file.
     * @return  Solved electrostatic potential vector (one value per node).
     * @throws std::runtime_error on I/O or solver failure.
     */
    VectorXd run(const std::string& configFile);

    /**
     * @brief Run the simulation and return reusable artifacts for integrations.
     *
     * @param configFile  Path to the JSON configuration file.
     * @return Mesh, solved electrostatic potential, and net doping used by the run.
     * @throws std::runtime_error on I/O or solver failure.
     */
    PoissonResult runWithResult(const std::string& configFile);
};

} // namespace vela
