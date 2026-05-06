#pragma once

#include "vela/core/Types.h"
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
 *   "contacts": [
 *     { "name": "cathode", "bias": 0.0 },
 *     { "name": "anode",   "bias": 0.0 }
 *   ]
 * }
 * @endcode
 *
 * Path resolution
 * ---------------
 * Relative paths in "mesh_file" and "output_vtk" are resolved relative to
 * the directory containing the config JSON.
 */
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
};

} // namespace vela
