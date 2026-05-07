#pragma once

#include "vela/mesh/MeshEntity.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <stdexcept>

namespace vela {

class BoxGeometryBuilder;

/**
 * @brief Container for the 2-D device mesh.
 *
 * Stores nodes, cells, edges, regions and contacts.
 * Edges are not stored explicitly in the input; call buildEdges() after
 * all nodes and cells have been added to generate them automatically.
 */
class DeviceMesh {
public:
    DeviceMesh() = default;

    // ------------------------------------------------------------------
    // Population helpers
    // ------------------------------------------------------------------
    void addNode   (const Node&    node);
    void addCell   (const Cell&    cell);
    void addRegion (const Region&  region);
    void addContact(const Contact& contact);

    /**
     * @brief Generate unique edges, compute edge lengths, and build box geometry.
     *
     * Must be called after all nodes and cells have been added.
     */
    void buildEdges();

    /**
     * @brief Compute node control volumes and edge coupling lengths.
     *
     * Uses 2-D Tri3 barycentric node volumes and cotangent box edge
     * couplings. buildEdges() also invokes this after generating edges.
     */
    void buildBoxGeometry();

    // ------------------------------------------------------------------
    // Queries
    // ------------------------------------------------------------------
    Index numNodes()    const { return nodes_.size();    }
    Index numEdges()    const { return edges_.size();    }
    Index numCells()    const { return cells_.size();    }
    Index numRegions()  const { return regions_.size();  }
    Index numContacts() const { return contacts_.size(); }

    const Node&    getNode   (Index id) const;
    const Edge&    getEdge   (Index id) const;
    const Cell&    getCell   (Index id) const;
    const Region&  getRegion (Index id) const;
    const Contact& getContact(Index id) const;

    // Accessors to full collections (read-only)
    const std::vector<Node>&    nodes()    const { return nodes_;    }
    const std::vector<Edge>&    edges()    const { return edges_;    }
    const std::vector<Cell>&    cells()    const { return cells_;    }
    const std::vector<Region>&  regions()  const { return regions_;  }
    const std::vector<Contact>& contacts() const { return contacts_; }

private:
    friend class BoxGeometryBuilder;

    void buildEdgesOnly();

    std::vector<Node>    nodes_;
    std::vector<Edge>    edges_;
    std::vector<Cell>    cells_;
    std::vector<Region>  regions_;
    std::vector<Contact> contacts_;
};

} // namespace vela
