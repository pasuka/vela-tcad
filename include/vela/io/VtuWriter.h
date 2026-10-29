#pragma once

#include "vela/core/Types.h"
#include "vela/mesh/DeviceMesh.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace vela {

/**
 * @brief Writer for VTK XML UnstructuredGrid files (``.vtu``).
 *
 * Output contract (docs/sentaurus_sdevice_workflow.md section 3.3): Float64
 * point/cell arrays named with Sentaurus TDR data names, vectors as
 * three-component arrays of the same name, ``region_id`` as Int32 CellData,
 * per-array units and run metadata in ``FieldData``.  ``ascii`` is the
 * default encoding; ``binary`` writes base64 inline arrays with a UInt32
 * byte-count header, which any VTK XML reader accepts.
 */
class VtuWriter {
public:
    enum class Encoding { Ascii, Binary };

    VtuWriter(const DeviceMesh& mesh, Encoding encoding = Encoding::Ascii);

    /// Node-centred scalar (one value per node, indexed by node id).
    void addPointScalar(const std::string& name, const std::string& unit,
                        const std::vector<Real>& values);
    /// Node-centred vector (three components per node).
    void addPointVector(const std::string& name, const std::string& unit,
                        const std::vector<Point3>& values);
    /// Cell-centred scalar (one value per cell).
    void addCellScalar(const std::string& name, const std::string& unit,
                       const std::vector<Real>& values);
    /// Cell-centred vector (three components per cell).
    void addCellVector(const std::string& name, const std::string& unit,
                       const std::vector<Point3>& values);
    /// Global metadata (one Float64 tuple).
    void addFieldValue(const std::string& name, Real value);
    /// Global metadata (one String tuple, written as VTK character codes).
    void addFieldString(const std::string& name, const std::string& value);

    /// Names of the point arrays added so far, in insertion order.
    std::vector<std::string> pointArrayNames() const;
    /// Names of the cell arrays added so far, in insertion order (excluding region_id).
    std::vector<std::string> cellArrayNames() const;

    /// Write the file; throws std::runtime_error on I/O failure.
    void write(const std::string& filename) const;
    /// Render the document to a string (used by write and by tests).
    std::string render() const;

    /// FieldData array name carrying the unit of @p arrayName.
    static std::string unitArrayName(const std::string& arrayName);

private:
    struct Array {
        std::string name;
        std::string unit;
        int components = 1;
        std::vector<Real> data;
    };
    struct FieldEntry {
        std::string name;
        bool isString = false;
        Real value = 0.0;
        std::string text;
    };

    void appendDataArray(std::string& out, const std::string& indent, const std::string& type,
                         const std::string& name, int components,
                         const std::vector<Real>& values) const;
    void appendIntArray(std::string& out, const std::string& indent, const std::string& name,
                        const std::vector<std::int32_t>& values) const;
    std::string encodeBinary(const std::vector<unsigned char>& bytes) const;

    const DeviceMesh& mesh_;
    Encoding encoding_;
    std::vector<Array> pointArrays_;
    std::vector<Array> cellArrays_;
    std::vector<FieldEntry> fieldData_;
};

/**
 * @brief One exportable Sentaurus TDR data set and the Vela quantity behind it.
 *
 * ``velaSource`` names the array of ``writeDDSolutionVTK`` (legacy ``.vtk``
 * name) the value is derived from; ``unit`` follows the Sentaurus manual
 * Appendix F.  ``components`` is 1 for scalars and 3 for vectors.
 */
struct TdrFieldSpec {
    std::string tdrName;
    std::string velaSource;
    std::string unit;
    int components = 1;
};

/// The catalogue of TDR data names the DD exporter can provide.
const std::vector<TdrFieldSpec>& tdrFieldCatalogue();

/// Look up a catalogue entry by TDR name (nullptr when not provided).
const TdrFieldSpec* tdrFieldSpec(const std::string& tdrName);

/**
 * @brief Selection of TDR fields requested for a ``.vtu`` export.
 *
 * Requests use Sentaurus ``Plot{}`` output names: ``Name`` selects the scalar
 * form (magnitude for vector quantities), ``Name/Vector`` the vector form.
 * An empty request list selects every catalogue entry (vectors as vectors).
 */
class TdrFieldSelection {
public:
    explicit TdrFieldSelection(const std::vector<std::string>& requests = {});

    /// True when the scalar form of @p tdrName was requested (or everything).
    bool wantsScalar(const std::string& tdrName) const;
    /// True when the vector form of @p tdrName was requested (or everything).
    bool wantsVector(const std::string& tdrName) const;
    /// True when either form was requested.
    bool wants(const std::string& tdrName) const;
    /// Requested names that no catalogue entry provides (``Name`` or ``Name/Vector``).
    std::vector<std::string> unavailable() const;
    /// Mark @p tdrName as delivered so it no longer counts as missing.
    void markProvided(const std::string& tdrName, bool vector);
    /// Requested and available names the exporter did not deliver.
    std::vector<std::string> undelivered() const;
    bool selectsAll() const { return all_; }

private:
    bool all_ = false;
    std::vector<std::pair<std::string, bool>> requests_; // (name, vector form)
    std::vector<std::pair<std::string, bool>> provided_;
};

} // namespace vela
