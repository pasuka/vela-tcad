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

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Vela TCAD output\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";

    ofs << "POINTS " << nodes.size() << " double\n";
    for (const auto& n : nodes)
        ofs << n.x << " " << n.y << " 0.0\n";

    const Index numCells = cells.size();
    ofs << "CELLS " << numCells << " " << numCells * 4 << "\n";
    for (const auto& c : cells) {
        ofs << "3";
        for (Index nid : c.node_ids)
            ofs << " " << nid;
        ofs << "\n";
    }

    ofs << "CELL_TYPES " << numCells << "\n";
    for (Index i = 0; i < numCells; ++i)
        ofs << "5\n";

    ofs << "CELL_DATA " << numCells << "\n";
    ofs << "SCALARS region_id int 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (const auto& c : cells)
        ofs << c.region_id << "\n";

    if (!ofs)
        throw std::runtime_error("Error writing VTK file: " + filename_);
}

void VTKWriter::addNodeScalar(const std::string& fieldName,
                              const std::vector<Real>& values)
{
    if (values.size() != mesh_.numNodes())
        throw std::invalid_argument("addNodeScalar: values size mismatch.");

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

void VTKWriter::addNodeVector(const std::string& fieldName,
                              const std::vector<Point3>& values)
{
    if (values.size() != mesh_.numNodes())
        throw std::invalid_argument("addNodeVector: values size mismatch.");

    std::ofstream ofs(filename_, std::ios::app);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open VTK file for appending: " + filename_);

    if (!pointDataHeaderWritten_) {
        ofs << "POINT_DATA " << mesh_.numNodes() << "\n";
        pointDataHeaderWritten_ = true;
    }

    ofs << "VECTORS " << fieldName << " double\n";
    for (const auto& v : values)
        ofs << v.x() << " " << v.y() << " " << v.z() << "\n";
}

void VTKWriter::addCellScalar(const std::string& fieldName,
                              const std::vector<Real>& values)
{
    if (values.size() != mesh_.numCells())
        throw std::invalid_argument("addCellScalar: values size mismatch.");

    std::ofstream ofs(filename_, std::ios::app);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open VTK file for appending: " + filename_);

    ofs << "SCALARS " << fieldName << " double 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (const auto& v : values)
        ofs << v << "\n";
}

void VTKWriter::addCellVector(const std::string& fieldName,
                              const std::vector<Point3>& values)
{
    if (values.size() != mesh_.numCells())
        throw std::invalid_argument("addCellVector: values size mismatch.");

    std::ofstream ofs(filename_, std::ios::app);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open VTK file for appending: " + filename_);

    ofs << "VECTORS " << fieldName << " double\n";
    for (const auto& v : values)
        ofs << v.x() << " " << v.y() << " " << v.z() << "\n";
}

} // namespace vela
