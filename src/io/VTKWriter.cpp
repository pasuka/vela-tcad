#include "vela/io/VTKWriter.h"
#include <fstream>
#include <stdexcept>

namespace vela {

VTKWriter::VTKWriter(const std::string& filename, const DeviceMesh& mesh)
    : filename_(filename), mesh_(mesh)
{}

void VTKWriter::write()
{
    std::ofstream ofs(filename_);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open VTK file for writing: " + filename_);

    const auto& nodes = mesh_.nodes();
    const auto& cells = mesh_.cells();

    // ------------------------------------------------------------------
    // VTK Legacy ASCII header
    // ------------------------------------------------------------------
    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Vela TCAD output\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";

    // ------------------------------------------------------------------
    // Points
    // ------------------------------------------------------------------
    ofs << "POINTS " << nodes.size() << " double\n";
    for (const auto& n : nodes) {
        ofs << n.x << " " << n.y << " 0.0\n";
    }

    // ------------------------------------------------------------------
    // Cells (only Tri3 supported in this stage)
    // ------------------------------------------------------------------
    // Each Tri3 cell: 1 count + 3 node ids -> 4 integers per cell
    const Index numCells = cells.size();
    ofs << "CELLS " << numCells << " " << numCells * 4 << "\n";
    for (const auto& c : cells) {
        ofs << "3";
        for (Index nid : c.node_ids)
            ofs << " " << nid;
        ofs << "\n";
    }

    // Cell types (VTK_TRIANGLE = 5)
    ofs << "CELL_TYPES " << numCells << "\n";
    for (Index i = 0; i < numCells; ++i)
        ofs << "5\n";

    // ------------------------------------------------------------------
    // Cell data: region id
    // ------------------------------------------------------------------
    ofs << "CELL_DATA " << numCells << "\n";
    ofs << "SCALARS region_id int 1\n";
    ofs << "LOOKUP_TABLE default\n";
    // Check stream state after all writes
    if (!ofs)
        throw std::runtime_error("Error writing VTK file: " + filename_);
}

void VTKWriter::addNodeScalar(const std::string& fieldName,
                              const std::vector<Real>& values)
{
    if (values.size() != mesh_.numNodes())
        throw std::invalid_argument("addNodeScalar: values size mismatch.");

    // Open in append mode so previous write() output is preserved
    std::ofstream ofs(filename_, std::ios::app);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open VTK file for appending: " + filename_);

    if (!pointDataHeaderWritten_) {
        ofs << "POINT_DATA " << mesh_.numNodes() << "\n";
        pointDataHeaderWritten_ = true;
    }

    ofs << "SCALARS " << fieldName << " double 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (const auto& v : values)
        ofs << v << "\n";
}

} // namespace vela
