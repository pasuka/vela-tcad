#pragma once

#include "vela/core/Types.h"
#include <string>
#include <vector>

namespace vela {

// ------------------------------------------------------------------
// Basic mesh entities
// ------------------------------------------------------------------

/// A mesh vertex with its geometric coordinates and dual-cell volume.
struct Node {
    Index  id     = 0;
    Real   x      = 0.0;
    Real   y      = 0.0;
    Real   volume = 0.0;  ///< Voronoi control-volume area [m^2] (computed later)
};

/// A mesh edge connecting two nodes.
struct Edge {
    Index  id     = 0;
    Index  n0     = 0;  ///< First node id
    Index  n1     = 0;  ///< Second node id
    Real   length = 0.0; ///< Euclidean length [m]
    Real   couple = 0.0; ///< Voronoi coupling length [m] (computed later)
};

/// A triangular (or higher-order) mesh cell.
struct Cell {
    Index             id        = 0;
    CellType          type      = CellType::Tri3;
    Index             region_id = 0;
    std::vector<Index> node_ids;
};

/// A named material/doping region composed of cells.
struct Region {
    Index              id       = 0;
    std::string        name;
    std::string        material;
    std::vector<Index> cell_ids;
};

/// A device contact composed of boundary nodes.
struct Contact {
    Index              id        = 0;
    std::string        name;
    Index              region_id = 0;
    std::vector<Index> node_ids;
};

} // namespace vela
