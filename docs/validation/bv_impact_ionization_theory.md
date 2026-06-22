# BV Impact-Ionization Theory And Vela Mapping

## Scope

This document describes the PN2D Sentaurus2018 BV validation target. The checked-in
BV deck uses `Avalanche(VanOverstraeten)`, not Okuto-Crowell, and sweeps the
Anode to `-20.0 V`. Okuto-Crowell is included only as contrast unless a fresh
Sentaurus run proves a different active model.

## Chynoweth Form

The common local ionization coefficient form is:

```math
\alpha(E) = A \exp(-B / |E|)
```

`A` is a prefactor in inverse length and `B` is a critical field. Vela evaluates
coefficient functions through `ImpactIonizationModel::electronCoefficient` and
`ImpactIonizationModel::holeCoefficient`.

## Van Overstraeten-de Man

Vela's `VanOverstraetenImpactIonization` uses low-field and high-field
coefficient sets selected by `switchField`, plus a temperature factor:

```math
\gamma(T) =
\frac{\tanh(\hbar\omega / (2 k_B T_\mathrm{ref}))}
{\tanh(\hbar\omega / (2 k_B T))}
```

```math
\alpha(E, T) = \gamma(T) A_\mathrm{region} \exp(-\gamma(T) B_\mathrm{region} / |E|)
```

The code mapping is
`ImpactIonizationModelConfig::{electronALow,electronAHigh,electronBLow,electronBHigh,holeALow,holeAHigh,holeBLow,holeBHigh,switchField,phononEnergy,referenceTemperature_K,temperature_K}`.

## Okuto-Crowell Contrast

Okuto-Crowell is commonly written as:

```math
\alpha(E) = a E^2 \exp(-(b/E)^2)
```

This is not the current PN2D BV fixture target. Do not implement
`OkutoCrowellImpactIonization` as part of this validation unless the freshly run
deck or exported parameter provenance contradicts the checked-in
`Avalanche(VanOverstraeten)` source.

## Driving Force

Sentaurus isothermal avalanche defaults are interpreted as quasi-Fermi-gradient
driven unless the deck enables a different option. Vela maps this through
`impact_ionization.driving_force = "quasi_fermi_gradient"` and
`current_approximation = "density_gradient"` for the Sentaurus-default SG
edge-current avalanche path.

The Sentaurus Device user guide describes `GradQuasiFermi` as the default
driving-force model for drift-diffusion simulations. For mesh elements touching
a contact, Sentaurus defaults to replacing the quasi-Fermi gradient with the
electric field unless `ComputeGradQuasiFermiAtContacts = UseQuasiFermi` is
requested. The same manual documents optional interpolation of avalanche
driving forces back to the electric field at low carrier concentration through
`RefDens_eGradQuasiFermi_ElectricField` and
`RefDens_hGradQuasiFermi_ElectricField`; the checked-in PN2D deck does not set
those keywords.

Sentaurus also enables boundary-layer parallel-field correction by default for
mobility and avalanche driving forces near interfaces and external boundaries.
Vela does not yet claim parity for that full boundary-layer model; the current
implementation mirrors only the documented contact-element electric-field
fallback for avalanche coefficient driving fields.

## SG Edge-Current Source

The same-dimension comparison target is not `alpha(E)` versus
`AvalancheGeneration`. Compare coefficient to coefficient, or compare
generation/source-integral to generation/source-integral. Vela's SG path records
`electronAlpha`, `holeAlpha`, electron/hole flux proxies, and edge/node source
integrals through `sgEdgeCurrentAvalancheSourceRecords`.

## Ionization Integral And Multiplication

The one-dimensional intuition is:

```math
M = \frac{1}{1 - \int \alpha \, dl}
```

This is useful as a diagnostic for field-line breakdown propensity. It is not
the production acceptance path for this PN2D Sentaurus-default SG edge-current
validation unless a later task explicitly adds a post-processing criterion.

## Current Validation Boundary

This validation promotes reproducibility, model identity, documented
field/current trends, and windowed current diagnostics. It does not claim full
absolute-current parity over the entire `0..-20 V` sweep and does not promote
hidden scalar calibration knobs such as `source_geometry_scale`.
