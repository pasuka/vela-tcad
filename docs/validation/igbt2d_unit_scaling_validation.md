# IGBT Unit Scaling Validation

This page records the checked-in IGBT reference_tcad chain for Vela
`unit_scaling`. It is engineering trend validation only and makes no calibration
claim.

## Scope

- low-current collector IV
- high-injection IV with stored charge proxy
- stored charge and quasi-static CV trend
- breakdown diagnostic max-field trend
- impact-ionization BV smoke diagnostic

## Result

The checked-in reports show monotone high-injection current, non-negative stored
charge trend data, and non-decreasing max field over the sampled reverse-bias
points. These are engineering trend validation fixtures, not calibrated IGBT
or calibrated breakdown diagnostic results. This document makes no calibration claim.
