#include "vela/material/MaterialDatabase.h"

namespace vela {

MaterialDatabase::MaterialDatabase()
{
    // Silicon
    Material si;
    si.name  = "Si";
    si.eps_r = 11.7;
    si.ni    = 1.0e16;   // [m^-3]
    si.mun   = 0.135;    // [m^2/V/s]
    si.mup   = 0.048;    // [m^2/V/s]
    db_["Si"] = si;

    // Silicon dioxide (insulator – ni, mun, mup remain 0)
    Material sio2;
    sio2.name  = "SiO2";
    sio2.eps_r = 3.9;
    db_["SiO2"] = sio2;
}

void MaterialDatabase::addMaterial(const Material& mat)
{
    db_[mat.name] = mat;
}

const Material& MaterialDatabase::getMaterial(const std::string& name) const
{
    auto it = db_.find(name);
    if (it == db_.end())
        throw std::out_of_range("MaterialDatabase: unknown material '" + name + "'");
    return it->second;
}

bool MaterialDatabase::hasMaterial(const std::string& name) const
{
    return db_.count(name) > 0;
}

} // namespace vela
