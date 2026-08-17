# Sentaurus `.par` parameter import contract

This document defines what may be imported from a Sentaurus parameter file
(`.par`) into a Vela deck, and with what loss of fidelity. It is the written
form of the machine-readable matrix in
[`scripts/sentaurus_parameter_map.py`](../scripts/sentaurus_parameter_map.py);
the code is authoritative, and
[`tests/regression/test_sentaurus_parameter_map.py`](../tests/regression/test_sentaurus_parameter_map.py)
keeps the two in step.

## Why a contract is needed

The tempting way to size this work is to count Vela JSON fields that look like
they could receive a `.par` value. That number is misleading. A shared field
name is not a shared equation. Concretely, in this repository:

- Vela's Auger recombination has constant `Cn`/`Cp`
  (`include/vela/physics/RecombinationModel.h`), while Sentaurus adds a
  carrier-density enhancement `(1 + H exp(-n/N0))`. There is nowhere for
  `H` and `N0` to go.
- `Material::atTemperature` (`src/material/Material.cpp`) scales both carrier
  mobilities by the same hardcoded `(T/T0)^-2.2`, while `.par` carries
  separate electron and hole exponents (`2.5, 2.2` for silicon).
- `ImpactIonizationModel` (`src/physics/ImpactIonizationModel.cpp`) implements
  `selberherr` and `van_overstraeten` only. The `Okuto`, `Lackner`, and
  `UniBo` sections present in every `.par` file have no counterpart.
- Vela has an electron quantum potential only
  (`include/vela/material/Material.h`); the hole parameters in
  `QuantumPotentialParameters` cannot be imported.

So the scope of this work is not "import the Sentaurus material library". It
is:

> Convert, traceably, the parameters of explicitly selected models that Vela
> already implements with equivalent semantics. Mark everything else.

## Status vocabulary

Every mapping row carries exactly one status.

| Status | Meaning | Importable |
| --- | --- | --- |
| `exact` | Vela evaluates the same formula from the same parameters. | yes |
| `frozen_at_temperature` | Numerically correct at one temperature only, because Vela stores a number where Sentaurus has a law. | yes, with the temperature recorded |
| `approximated` | Vela implements a reduced form; importing loses physics. | only with `--allow-lossy` |
| `unsupported_formula` | The section exists in Vela but this `Formula` or named sub-model does not. | never |
| `unsupported_model` | Vela has no counterpart. | never |

A sixth outcome, `unmapped`, is produced at classification time rather than
stored in the matrix. It means the parameter was not reached: either no row
describes it, or its model is not active. `unmapped` is **not** an error — a
`.par` file is a candidate library, and most of it is inert in any given run.

`--allow-lossy` widens the gate to `approximated`. It never admits
`unsupported_formula` or `unsupported_model`: those represent physics Vela
cannot reproduce at any temperature, and accepting them would produce a deck
that silently solves a different problem.

## The activated-model rule

A `.par` file lists candidates, not selections. `models.par` in this
repository defines `Slotboom`, `OldSlotboom`, `Bennett`, and `delAlamo`
bandgap narrowing simultaneously, along with five impact-ionization models and
several mobility formulas. The file alone cannot say which are used.

The rule is therefore:

> Importing parameters must never switch on a physics model. Parameters are
> generated only for models that were independently shown to be active.

The activation context comes from the SDevice execution IR
(`scripts/sentaurus_execution_ir.py`), which classifies the models a `.cmd`
deck actually selects. Each matrix row names the model it depends on via
`requires_model`; a row whose model is absent from the active set classifies
as `unmapped` and contributes nothing.

Rows for models Vela does not implement are gated on *their own* name — for
example the `Lackner` rows require the `Lackner` model. Because that name is
not in the execution IR's supported list, those rows can never be reached
through a legitimate import, but if the gate is ever widened they immediately
report `unsupported_model` instead of being silently forgotten.

When no activation context is supplied at all, every model-gated row is
`unmapped`. Only material constants such as permittivity import without one.

## Inert values

Strictness becomes useless if it blocks every real file. Several parameters
that Vela cannot represent ship at values that switch them off:
`HighFieldDependence.alpha = 0`, `ku = kv = 1`,
`QuantumPotentialParameters.nu = 0`. Dropping a term that Sentaurus is also
not using loses nothing.

A matrix row may therefore declare a `neutral_value`. When the observed value
matches it, the status is demoted to `exact` with a reason recording why.
Both carriers must be neutral: a term switched off for electrons but active
for holes still changes the solve. A file that moves such a parameter off its
neutral value fails closed, which is the behaviour that makes the rule safe.

## Temperature and the intrinsic density

Two consequences of Vela's data model drive most of the `frozen_at_temperature`
rows.

**Vela stores numbers where Sentaurus stores laws.** `Material` holds scalar
`bandgap_eV`, `Nc_m3`, `Nv_m3`, `mun`, and `mup`. There is no Varshni
evaluator and no DOS-mass temperature law. The emitter must evaluate those
formulas once, at a declared temperature, and record which temperature it
used. Note that `Bandgap.Tpar` is `0` in this repository's `models.par`, so
`Eg0 = 1.16964 eV` is a 0 K reference, not the 300 K gap.

**Vela never derives `ni`.** `Material::atTemperature` only rescales an `ni`
that is already present, and `MaterialDatabase` supplies a built-in silicon
default. A materials JSON that carries `bandgap_eV`, `Nc_m3`, and `Nv_m3` but
no `ni` therefore silently keeps the built-in value, and every carrier density
in the solve is wrong while every field looks populated.

Any material emitter must consequently write `ni` explicitly, alongside
`bandgap_eV`, `Nc_m3`, `Nv_m3`, the temperature it froze at, and the formula
and statistics it assumed.

**`dEg0` must not be counted twice.** `docs/config_schema.md` documents
`dEg0(OldSlotboom)` as being expressed through the material `ni` override,
while `bandgap_eV` carries the base Varshni gap. Writing the narrowing into
both would double-count it. The matrix pins this: the `dEg0(OldSlotboom)` row
targets `materials[].ni` and the `Eg0`/`alpha`/`beta`/`Tpar` rows target
`materials[].bandgap_eV`, and a test asserts the two never overlap.

## Units

The matrix is a semantic contract and deliberately says nothing about units.
Unit handling belongs to the emitter, and follows one rule: the syntax IR
preserves the original value and the original unit annotation exactly as
written, and no conversion happens before the emitter has decided which Vela
input a value feeds.

This matters because Vela's own input units are not uniform. Deck input has
two modes — legacy SI, and `scaling.mode: "unit_scaling"` for cm/µm TCAD units
— and band-to-band tunnelling is a documented special case whose canonical
keys are SI regardless of the active mode, with `*_per_cm` aliases provided for
direct Sentaurus transcription. A converter that normalises early cannot
express that.

## Pipeline position

```
.par ──▶ sentaurus_parameter_ir.py ──▶ syntax IR ──┐
                                                    ├─▶ mapping matrix ──▶ emitter ──▶ Vela JSON
.cmd ──▶ sentaurus_execution_ir.py ─▶ active models ┘
```

The syntax IR is lossless and performs no physics interpretation. The mapping
matrix decides what is allowed. Only the emitter converts units and evaluates
temperature-dependent formulas. Keeping these three concerns apart is what
makes the coverage numbers auditable: syntax coverage, activated-model
coverage, and parameter coverage are measured separately and never conflated.

## Usage

Dump the matrix for review:

```bash
python3 scripts/sentaurus_parameter_map.py --dump-matrix
```

Classify a parsed `.par` against an activation context:

```bash
python3 scripts/sentaurus_parameter_ir.py path/to/models.par --out models.ir.json
python3 scripts/sentaurus_parameter_map.py --classify models.ir.json \
    --active-model SRH --active-model TempDependence
```

The command exits non-zero when any active parameter cannot be imported, and
prints the offending rows with the reason for each.

## Current coverage

Against `reference_tcad/pn2d_sentaurus2018/source/models.par` with the PN2D
breakdown model set active, the matrix reports a single honest blocker:
`Auger.H` and `Auger.N0`, the carrier-density enhancement Vela does not
implement. Six parameters classify as `approximated` — the two mobility
temperature exponents, and the three Auger polynomial coefficients plus the
constant-mobility exponent. Everything else is either importable or inert.

`reference_tcad/bvmethods_sentaurus2018/source/models.par`, which contains only
the SRH `Scharfetter` section, imports with no blockers and no approximations
under an `SRH` + `TempDependence` context. That file is the reference case for
what a clean import looks like.
