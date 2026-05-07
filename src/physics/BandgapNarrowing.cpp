#include "vela/physics/BandgapNarrowing.h"

namespace vela {

Real BandgapNarrowing::deltaEg(Real, Real, Real) const
{
    return 0.0;
}

Real NoBandgapNarrowing::deltaEg(Real, Real, Real) const
{
    return 0.0;
}

} // namespace vela
