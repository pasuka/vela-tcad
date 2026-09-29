#include "vela/io/VtuWriter.h"

#include <algorithm>
#include <charconv>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace vela {

namespace {

std::string xmlEscape(const std::string& text)
{
    std::string out;
    out.reserve(text.size());
    for (const char c : text) {
        switch (c) {
        case '&': out += "&amp;"; break;
        case '<': out += "&lt;"; break;
        case '>': out += "&gt;"; break;
        case '"': out += "&quot;"; break;
        default: out += c; break;
        }
    }
    return out;
}

std::string formatReal(Real value)
{
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<Real>::max_digits10) << value;
    return oss.str();
}

const char kBase64[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string base64(const std::vector<unsigned char>& bytes)
{
    std::string out;
    out.reserve((bytes.size() + 2) / 3 * 4);
    std::size_t i = 0;
    for (; i + 2 < bytes.size(); i += 3) {
        const unsigned int triple = (bytes[i] << 16) | (bytes[i + 1] << 8) | bytes[i + 2];
        out += kBase64[(triple >> 18) & 0x3F];
        out += kBase64[(triple >> 12) & 0x3F];
        out += kBase64[(triple >> 6) & 0x3F];
        out += kBase64[triple & 0x3F];
    }
    if (i < bytes.size()) {
        unsigned int triple = bytes[i] << 16;
        if (i + 1 < bytes.size())
            triple |= bytes[i + 1] << 8;
        out += kBase64[(triple >> 18) & 0x3F];
        out += kBase64[(triple >> 12) & 0x3F];
        out += (i + 1 < bytes.size()) ? kBase64[(triple >> 6) & 0x3F] : '=';
        out += '=';
    }
    return out;
}

template <typename T>
void appendBytes(std::vector<unsigned char>& out, const T& value)
{
    unsigned char buffer[sizeof(T)];
    std::memcpy(buffer, &value, sizeof(T));
    out.insert(out.end(), buffer, buffer + sizeof(T));
}

} // namespace

VtuWriter::VtuWriter(const DeviceMesh& mesh, Encoding encoding)
    : mesh_(mesh), encoding_(encoding)
{}

std::string VtuWriter::unitArrayName(const std::string& arrayName)
{
    return "Unit:" + arrayName;
}

void VtuWriter::addPointScalar(const std::string& name, const std::string& unit,
                               const std::vector<Real>& values)
{
    if (values.size() != mesh_.numNodes())
        throw std::invalid_argument("VtuWriter::addPointScalar: values size mismatch for " + name);
    pointArrays_.push_back({name, unit, 1, values});
}

void VtuWriter::addPointVector(const std::string& name, const std::string& unit,
                               const std::vector<Point3>& values)
{
    if (values.size() != mesh_.numNodes())
        throw std::invalid_argument("VtuWriter::addPointVector: values size mismatch for " + name);
    Array array{name, unit, 3, {}};
    array.data.reserve(values.size() * 3);
    for (const Point3& v : values) {
        array.data.push_back(v.x());
        array.data.push_back(v.y());
        array.data.push_back(v.z());
    }
    pointArrays_.push_back(std::move(array));
}

void VtuWriter::addCellScalar(const std::string& name, const std::string& unit,
                              const std::vector<Real>& values)
{
    if (values.size() != mesh_.numCells())
        throw std::invalid_argument("VtuWriter::addCellScalar: values size mismatch for " + name);
    cellArrays_.push_back({name, unit, 1, values});
}

void VtuWriter::addCellVector(const std::string& name, const std::string& unit,
                              const std::vector<Point3>& values)
{
    if (values.size() != mesh_.numCells())
        throw std::invalid_argument("VtuWriter::addCellVector: values size mismatch for " + name);
    Array array{name, unit, 3, {}};
    array.data.reserve(values.size() * 3);
    for (const Point3& v : values) {
        array.data.push_back(v.x());
        array.data.push_back(v.y());
        array.data.push_back(v.z());
    }
    cellArrays_.push_back(std::move(array));
}

void VtuWriter::addFieldValue(const std::string& name, Real value)
{
    fieldData_.push_back({name, false, value, {}});
}

void VtuWriter::addFieldString(const std::string& name, const std::string& value)
{
    fieldData_.push_back({name, true, 0.0, value});
}

std::vector<std::string> VtuWriter::pointArrayNames() const
{
    std::vector<std::string> names;
    for (const Array& array : pointArrays_)
        names.push_back(array.name);
    return names;
}

std::vector<std::string> VtuWriter::cellArrayNames() const
{
    std::vector<std::string> names;
    for (const Array& array : cellArrays_)
        names.push_back(array.name);
    return names;
}

std::string VtuWriter::encodeBinary(const std::vector<unsigned char>& bytes) const
{
    std::vector<unsigned char> framed;
    framed.reserve(bytes.size() + 4);
    appendBytes(framed, static_cast<std::uint32_t>(bytes.size()));
    framed.insert(framed.end(), bytes.begin(), bytes.end());
    return base64(framed);
}

void VtuWriter::appendDataArray(std::string& out, const std::string& indent, const std::string& type,
                                const std::string& name, int components,
                                const std::vector<Real>& values) const
{
    out += indent + "<DataArray type=\"" + type + "\" Name=\"" + xmlEscape(name) + "\"";
    if (components != 1)
        out += " NumberOfComponents=\"" + std::to_string(components) + "\"";
    out += " format=\"" + std::string(encoding_ == Encoding::Ascii ? "ascii" : "binary") + "\">\n";
    if (encoding_ == Encoding::Ascii) {
        std::size_t column = 0;
        out += indent + "  ";
        for (std::size_t i = 0; i < values.size(); ++i) {
            out += formatReal(values[i]);
            ++column;
            if (i + 1 < values.size()) {
                if (column % static_cast<std::size_t>(std::max(components, 1) * 2) == 0)
                    out += "\n" + indent + "  ";
                else
                    out += ' ';
            }
        }
        out += "\n";
    } else {
        std::vector<unsigned char> bytes;
        bytes.reserve(values.size() * sizeof(double));
        for (const Real value : values)
            appendBytes(bytes, static_cast<double>(value));
        out += indent + "  " + encodeBinary(bytes) + "\n";
    }
    out += indent + "</DataArray>\n";
}

void VtuWriter::appendIntArray(std::string& out, const std::string& indent, const std::string& name,
                               const std::vector<std::int32_t>& values) const
{
    out += indent + "<DataArray type=\"Int32\" Name=\"" + xmlEscape(name) + "\" format=\"" +
        std::string(encoding_ == Encoding::Ascii ? "ascii" : "binary") + "\">\n";
    if (encoding_ == Encoding::Ascii) {
        out += indent + "  ";
        for (std::size_t i = 0; i < values.size(); ++i) {
            out += std::to_string(values[i]);
            if (i + 1 < values.size())
                out += ((i + 1) % 12 == 0) ? "\n" + indent + "  " : " ";
        }
        out += "\n";
    } else {
        std::vector<unsigned char> bytes;
        bytes.reserve(values.size() * sizeof(std::int32_t));
        for (const std::int32_t value : values)
            appendBytes(bytes, value);
        out += indent + "  " + encodeBinary(bytes) + "\n";
    }
    out += indent + "</DataArray>\n";
}

std::string VtuWriter::render() const
{
    const auto& nodes = mesh_.nodes();
    const auto& cells = mesh_.cells();
    std::string out;
    out += "<?xml version=\"1.0\"?>\n";
    out += "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
           "header_type=\"UInt32\">\n";
    out += "  <UnstructuredGrid>\n";

    // FieldData: run metadata plus one unit string per data array.
    out += "    <FieldData>\n";
    auto appendString = [&](const std::string& name, const std::string& text) {
        out += "      <DataArray type=\"String\" Name=\"" + xmlEscape(name) +
            "\" NumberOfTuples=\"1\" format=\"ascii\">\n        ";
        for (const unsigned char c : text)
            out += std::to_string(static_cast<int>(c)) + " ";
        out += "0\n      </DataArray>\n";
    };
    for (const FieldEntry& entry : fieldData_) {
        if (entry.isString) {
            appendString(entry.name, entry.text);
        } else {
            out += "      <DataArray type=\"Float64\" Name=\"" + xmlEscape(entry.name) +
                "\" NumberOfTuples=\"1\" format=\"ascii\">\n        " + formatReal(entry.value) +
                "\n      </DataArray>\n";
        }
    }
    for (const Array& array : pointArrays_)
        appendString(unitArrayName(array.name), array.unit);
    for (const Array& array : cellArrays_)
        appendString(unitArrayName(array.name), array.unit);
    out += "    </FieldData>\n";

    out += "    <Piece NumberOfPoints=\"" + std::to_string(nodes.size()) + "\" NumberOfCells=\"" +
        std::to_string(cells.size()) + "\">\n";

    out += "      <Points>\n";
    std::vector<Real> coordinates;
    coordinates.reserve(nodes.size() * 3);
    for (const auto& node : nodes) {
        coordinates.push_back(node.x);
        coordinates.push_back(node.y);
        coordinates.push_back(0.0);
    }
    appendDataArray(out, "        ", "Float64", "Points", 3, coordinates);
    out += "      </Points>\n";

    out += "      <Cells>\n";
    std::vector<std::int32_t> connectivity;
    std::vector<std::int32_t> offsets;
    std::vector<std::int32_t> types;
    connectivity.reserve(cells.size() * 3);
    offsets.reserve(cells.size());
    types.reserve(cells.size());
    std::int32_t offset = 0;
    for (const auto& cell : cells) {
        for (const Index nid : cell.node_ids)
            connectivity.push_back(static_cast<std::int32_t>(nid));
        offset += static_cast<std::int32_t>(cell.node_ids.size());
        offsets.push_back(offset);
        types.push_back(cell.node_ids.size() == 3 ? 5 : (cell.node_ids.size() == 4 ? 9 : 7));
    }
    appendIntArray(out, "        ", "connectivity", connectivity);
    appendIntArray(out, "        ", "offsets", offsets);
    appendIntArray(out, "        ", "types", types);
    out += "      </Cells>\n";

    out += "      <PointData>\n";
    for (const Array& array : pointArrays_)
        appendDataArray(out, "        ", "Float64", array.name, array.components, array.data);
    out += "      </PointData>\n";

    out += "      <CellData>\n";
    std::vector<std::int32_t> regionIds;
    regionIds.reserve(cells.size());
    for (const auto& cell : cells)
        regionIds.push_back(static_cast<std::int32_t>(cell.region_id));
    appendIntArray(out, "        ", "region_id", regionIds);
    for (const Array& array : cellArrays_)
        appendDataArray(out, "        ", "Float64", array.name, array.components, array.data);
    out += "      </CellData>\n";

    out += "    </Piece>\n";
    out += "  </UnstructuredGrid>\n";
    out += "</VTKFile>\n";
    return out;
}

void VtuWriter::write(const std::string& filename) const
{
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open())
        throw std::runtime_error("Cannot open VTU file for writing: " + filename);
    const std::string text = render();
    ofs.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!ofs)
        throw std::runtime_error("Error writing VTU file: " + filename);
}

// ---------------------------------------------------------------------------
// TDR field catalogue
// ---------------------------------------------------------------------------

const std::vector<TdrFieldSpec>& tdrFieldCatalogue()
{
    // Units follow the Sentaurus manual Appendix F (Tables 195/196); the Vela
    // source names are the legacy writeDDSolutionVTK arrays each value is
    // derived from (converted to the listed unit by the exporter).
    static const std::vector<TdrFieldSpec> catalogue = {
        {"ElectrostaticPotential", "Potential", "V", 1},
        {"eQuasiFermiPotential", "ElectronQuasiFermi", "V", 1},
        {"hQuasiFermiPotential", "HoleQuasiFermi", "V", 1},
        {"eDensity", "Electrons", "cm^-3", 1},
        {"hDensity", "Holes", "cm^-3", 1},
        {"DopingConcentration", "NetDoping", "cm^-3", 1},
        {"DonorConcentration", "DopingModel::donors", "cm^-3", 1},
        {"AcceptorConcentration", "DopingModel::acceptors", "cm^-3", 1},
        {"SpaceCharge", "SpaceCharge", "cm^-3", 1},
        {"ElectricField", "ElectricFieldVector", "V/cm", 3},
        {"eCurrentDensity", "SentaurusElectronCurrentDensityVector", "A/cm^2", 3},
        {"hCurrentDensity", "SentaurusHoleCurrentDensityVector", "A/cm^2", 3},
        {"TotalCurrentDensity", "SentaurusTotalCurrentDensityVector", "A/cm^2", 3},
        {"ConductionCurrentDensity", "SentaurusTotalCurrentDensityVector", "A/cm^2", 3},
        {"eMobility", "ElectronMobilityCm2PerVs", "cm^2/(V*s)", 1},
        {"hMobility", "HoleMobilityCm2PerVs", "cm^2/(V*s)", 1},
        {"eVelocity", "ElectronMobilityCm2PerVs*ElectronImpactIonizationDrive", "cm/s", 1},
        {"hVelocity", "HoleMobilityCm2PerVs*HoleImpactIonizationDrive", "cm/s", 1},
        {"srhRecombination", "SRHRecombinationCm3PerS", "cm^-3*s^-1", 1},
        {"AugerRecombination", "AugerRecombinationCm3PerS", "cm^-3*s^-1", 1},
        {"Band2BandGeneration", "Band2BandGeneration", "cm^-3*s^-1", 1},
        {"AvalancheGeneration", "AvalancheGeneration", "cm^-3*s^-1", 1},
        {"TotalRecombination", "SRHRecombination+AugerRecombination-Band2BandGeneration-AvalancheGeneration",
         "cm^-3*s^-1", 1},
        {"eAlphaAvalanche", "ElectronAlphaAvalanche", "cm^-1", 1},
        {"hAlphaAvalanche", "HoleAlphaAvalanche", "cm^-1", 1},
        {"eIonIntegral", "ElectronIonIntegral", "1", 1},
        {"hIonIntegral", "HoleIonIntegral", "1", 1},
        {"MeanIonIntegral", "MeanIonIntegral", "1", 1},
        {"BandGap", "BandGap", "eV", 1},
        {"BandgapNarrowing", "BandgapNarrowing", "eV", 1},
        {"EffectiveBandGap", "BandGap-BandgapNarrowing", "eV", 1},
        {"ElectronAffinity", "ElectronAffinity", "eV", 1},
        {"ConductionBandEnergy", "ConductionBandEnergy", "eV", 1},
        {"ValenceBandEnergy", "ValenceBandEnergy", "eV", 1},
        {"EffectiveIntrinsicDensity", "EffectiveIntrinsicDensity", "cm^-3", 1},
        {"eEffectiveStateDensity", "densityOfStatesNc", "cm^-3", 1},
        {"hEffectiveStateDensity", "densityOfStatesNv", "cm^-3", 1},
        {"eEparallel", "ElectronEparallel", "V/cm", 1},
        {"hEparallel", "HoleEparallel", "V/cm", 1},
        {"eEnormal", "ElectronEnormal", "V/cm", 1},
        {"hEnormal", "HoleEnormal", "V/cm", 1},
        {"eGradQuasiFermi", "ElectronGradQuasiFermiVector", "V/cm", 3},
        {"hGradQuasiFermi", "HoleGradQuasiFermiVector", "V/cm", 3},
    };
    return catalogue;
}

const TdrFieldSpec* tdrFieldSpec(const std::string& tdrName)
{
    for (const TdrFieldSpec& spec : tdrFieldCatalogue()) {
        if (spec.tdrName == tdrName)
            return &spec;
    }
    return nullptr;
}

// ---------------------------------------------------------------------------
// TdrFieldSelection
// ---------------------------------------------------------------------------

TdrFieldSelection::TdrFieldSelection(const std::vector<std::string>& requests)
    : all_(requests.empty())
{
    for (const std::string& request : requests) {
        const std::size_t slash = request.find('/');
        const std::string name = request.substr(0, slash);
        bool vector = false;
        if (slash != std::string::npos) {
            const std::string suffix = request.substr(slash + 1);
            vector = suffix == "Vector" || suffix == "vector";
        }
        if (name.empty())
            continue;
        const std::pair<std::string, bool> entry{name, vector};
        if (std::find(requests_.begin(), requests_.end(), entry) == requests_.end())
            requests_.push_back(entry);
    }
}

bool TdrFieldSelection::wantsScalar(const std::string& tdrName) const
{
    if (all_) {
        const TdrFieldSpec* spec = tdrFieldSpec(tdrName);
        return spec != nullptr && spec->components == 1;
    }
    return std::find(requests_.begin(), requests_.end(), std::make_pair(tdrName, false)) != requests_.end();
}

bool TdrFieldSelection::wantsVector(const std::string& tdrName) const
{
    if (all_) {
        const TdrFieldSpec* spec = tdrFieldSpec(tdrName);
        return spec != nullptr && spec->components == 3;
    }
    return std::find(requests_.begin(), requests_.end(), std::make_pair(tdrName, true)) != requests_.end();
}

bool TdrFieldSelection::wants(const std::string& tdrName) const
{
    return wantsScalar(tdrName) || wantsVector(tdrName);
}

std::vector<std::string> TdrFieldSelection::unavailable() const
{
    std::vector<std::string> missing;
    for (const auto& [name, vector] : requests_) {
        const TdrFieldSpec* spec = tdrFieldSpec(name);
        // Scalar magnitude of a vector quantity is always derivable; the
        // vector form of a scalar quantity is not.
        if (spec == nullptr || (vector && spec->components != 3))
            missing.push_back(vector ? name + "/Vector" : name);
    }
    return missing;
}

void TdrFieldSelection::markProvided(const std::string& tdrName, bool vector)
{
    provided_.emplace_back(tdrName, vector);
}

std::vector<std::string> TdrFieldSelection::undelivered() const
{
    std::vector<std::string> missing;
    for (const auto& request : requests_) {
        const TdrFieldSpec* spec = tdrFieldSpec(request.first);
        if (spec == nullptr || (request.second && spec->components != 3))
            continue;
        if (std::find(provided_.begin(), provided_.end(), request) == provided_.end())
            missing.push_back(request.second ? request.first + "/Vector" : request.first);
    }
    return missing;
}

} // namespace vela
