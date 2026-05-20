# PN diode reference fixture

This directory is the first complete reference_tcad validation chain for a
unit_scaling PN diode deck. The source data is an explicit CSV/text export; no
proprietary binary formats are parsed.

## Structure

- Device: 2D silicon PN diode.
- Regions: p_region and n_region.
- Contacts: anode on p_region, cathode on n_region.
- Junction: abrupt junction at the p_region / n_region interface.
- Input units: length in um, concentration in cm^-3, voltage in V,
  capacitance in F/m, and high-field diagnostic in V/cm.
- Vela decks: generated with `"scaling": {"mode": "unit_scaling"}`.

## Validation coverage

- forward IV: anode bias sweep verifies monotonic current increase.
- reverse quasi-static CV: reverse anode sweep verifies finite capacitance.
- reverse BV: reverse anode sweep verifies max field non-decreasing with bias.

The first checked-in dataset is intentionally small and trend oriented. It is
used to verify signs, trends, and key orders of magnitude, not calibration.
