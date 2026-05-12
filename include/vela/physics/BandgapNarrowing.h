#pragma once

#include "vela/core/Types.h"

namespace vela {

class BandgapNarrowing {
public:
    virtual ~BandgapNarrowing() = default;

    /// Return the effective bandgap narrowing DeltaEg [eV] at a node.
    virtual Real deltaEg(Real netDoping, Real n, Real p) const;
};

class NoBandgapNarrowing final : public BandgapNarrowing {
public:
    Real deltaEg(Real netDoping, Real n, Real p) const override;
};

} // namespace vela
