# Sentaurus 2018 BVmethods NMOS reference

This directory preserves the text inputs used by the Sentaurus Training
`BVmethods` NMOS example and the method mapping used by Vela validation. Binary
TDR/PLT/log outputs remain generated artifacts under `build*/reference_tcad/`
and are not checked in.

## Source inventory

| Checked-in file | Training workbench node | Method |
|---|---:|---|
| `source/bvmethods_nmos_sde.cmd` | 1 | symmetric NMOS structure and mesh |
| `source/aba_poisson_sdevice.cmd` | 3 | ABA, Poisson-only sweep |
| `source/aba_coupled_sdevice.cmd` | 4 | ABA postprocessing on coupled DD |
| `source/external_resistor_sdevice.cmd` | 5 | `1e7 ohm*um` series-resistor load line |
| `source/voltage_to_current_sdevice.cmd` | 6 | voltage ramp to 6 V, then current control |
| `source/continuation_sdevice.cmd` | 7 | Sentaurus continuation curve tracing |
| `source/transient_sdevice.cmd` | 8 | transient voltage ramp with resistor |
| `source/models.par` | 3--8 | shared SRH lifetime parameters |

The SDevice files retain the original contact names, physical models, target
values, step controls and current direction. File/output names were made
descriptive so the decks can coexist outside Sentaurus Workbench.

Run the structure deck first. It builds `bvmethods_nmos_half_msh.tdr`, mirrors
it about `x=0`, renames the mirrored drain contact to source, and writes
`bvmethods_nmos_msh.tdr`. Each SDevice deck expects that final mesh in its
working directory.

## Vela template mapping

- `vela/materials_sentaurus2018.json` preserves the Si/Nitride material
  overrides used by the boundary-method and IIC replay scripts. The values
  are unchanged; these inputs are now stored with the BVmethods case.
- `configs/templates/bvmethods_nmos_external_resistor.template.json` maps the
  series-resistor load-line method.
- `configs/templates/bvmethods_nmos_voltage_to_current.template.json` maps the
  voltage-to-current boundary switch.
- IIC/ABA validation remains represented by the existing BV/IIC configuration
  workflow and the checked-in validation scripts.
- Sentaurus `Continuation` is mapped to Vela's existing pseudo-arclength BV
  mode. The five non-transient references and the current Vela closure status
  are sealed in `bvmethods_nontransient_validation_20260817.json`.
- Sentaurus `Transient` remains archived as a reference only and is explicitly
  outside the non-transient closure scope.

The validated Sentaurus threshold voltages are `6.379791636 V` for the
external-resistor method and `6.383184201 V` for voltage-to-current control at
the example drain-current criterion. No mobility or avalanche coefficient
scaling is used by the Vela templates.

`scripts/freeze_bvmethods_nontransient.py` regenerates the five compact
reference curves from the archived PLT files. The Continuation deck is
prepared by `scripts/prepare_bvmethods_nmos_continuation.py`, which copies an
accepted base deck's physics verbatim and changes only its sweep controller.
The preparer requires two adjacent accepted states: the earlier state supplies
the initial secant direction and the later state is the restart anchor. This
is a numerical initialization only; it does not modify the device equations.

The current NMOS continuation status remains **pending**. Analytic-tangent
trials with and without Enormal mobility advanced only about `1.5e-10 V`
before reaching the minimum arclength step. A bounded two-state secant trial
from 6.000000 V to 6.056459 V produced no accepted point after the anchor.
These results are sealed in `continuation_diagnostic_20260818.json`; they do
not block the already passed IIC, external-resistor, or voltage-to-current
branches and do not justify adding another physical model.
