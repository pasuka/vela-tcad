# SimpleMOS M39 BGN-chain first-divergence audit

## Scope

M39 is a read-only audit of the frozen n23, Vd = 0.05 V, Vg = 0.05 V
states used by M27 and M29.  It does not rerun Sentaurus or Vela and does not
change a production default.  The audit order is deliberately fail-closed:

1. C-1: identify the effective Sentaurus BGN model from the solver log;
2. C0: verify that donor and acceptor inputs are independently observable;
3. C1: compare the exported narrowing against Vela OldSlotboom using candidate
   impurity-concentration conventions;
4. C2/C3: audit band-gap convention and frozen Boltzmann density identities;
5. C4/C5: inspect previously localized nodes/edges and the conduction/valence
   band split.

## First divergence

The first divergence occurs before the numeric BGN chain, at C-1 model
selection.  M29 generated its nominal `bgn_off` deck by omitting
`EffectiveIntrinsicDensity(OldSlotboom)`.  The corresponding sdevice
2022.03-SP2 solver log states:

```text
Bandgap narrowing model: Bennett/Wilson with bandgap narrowing (no Fermi)
```

The control is therefore Bennett/Wilson versus OldSlotboom, not BGN off versus
BGN on.  M29 terminal currents remain measured results, but the former
BGN-on/off causal label and its effect estimate are invalid.

M29 also omitted `DonorConcentration` and `AcceptorConcentration` from its
plot list.  Its generated `doping.csv` contains zeros and cannot support C0.
M39 records C0 as `not_scored_missing_fields`; it never interprets those zeros
as physical doping.  Downstream diagnostics use the M27-full input doping only
after exact node-id and coordinate identity is established for all nodes.

## Qualified OldSlotboom chain

The two M27 states provide complete, independently exported doping fields and
both solver logs select OldSlotboom.  Their results are:

| Audit | Result | Interpretation |
|---|---:|---|
| C0 donor/acceptor identity | max relative difference `0` | input doping is identical |
| C1 total impurity | max absolute difference `2.776e-17 eV` | Sentaurus and Vela OldSlotboom formulas close |
| C1 absolute net doping | max absolute difference `2.540e-2 eV` | net doping is not the selected convention |
| C3 density joint-fit residual | max RMS `1.362e-14` | frozen density identity closes after constant/convention fit |
| C5 conduction share | `0.5` | symmetric narrowing split |
| C5 valence share | `0.5` | symmetric narrowing split |
| C5 gap slope | `-1.0` | exported band gap narrows by exported delta Eg |

The M29 OldSlotboom/SRH-off state gives the same C1 total-impurity closure when
the qualified input-doping oracle is used.  The nominal off state differs from
Vela OldSlotboom by up to `7.447e-2 eV`, but this is explicitly a cross-model
Bennett/Wilson-versus-OldSlotboom diagnostic and not a parity score.

## Conclusion and next control

M39 does not find a discrepancy in Vela's OldSlotboom delta-Eg formula, the
total-impurity concentration convention, or the symmetric band split.  It
instead invalidates the M29 no-BGN control at model selection.  A subsequent
causal experiment must first obtain a documented Sentaurus configuration that
actually disables BGN (or construct a parameter-equivalent no-narrowing
control) and must export donor and acceptor fields.  Until then, BGN cannot be
claimed as the closed cause of the deep-off current gap from M29 alone.
