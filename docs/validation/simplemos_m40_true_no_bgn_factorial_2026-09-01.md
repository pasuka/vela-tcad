# SimpleMOS M40 explicit no-BGN factorial

## Manual and execution contract

M39 showed that omitting `EffectiveIntrinsicDensity(OldSlotboom)` selects the
Sentaurus default Bennett/Wilson model and does not disable BGN.  M40 therefore
uses the explicit syntax documented on page 305 of the T-2022.03 Sentaurus
Device User Guide:

```text
EffectiveIntrinsicDensity(NoBandGapNarrowing)
```

The paired OldSlotboom cells use
`EffectiveIntrinsicDensity(BandGapNarrowing(OldSlotboom))`.  The manual PDF,
release banner, decks, solver logs, TDRs, and current files are fingerprinted.
The manual PDF SHA-256 is
`f33068f5040c1c8c7993bdc1c86a0d777925792780ab4405752accd682706e88`.

Four new Sentaurus states and four new Vela states were solved at n23,
Vd = 0.05 V, Vg = 0.05 V.  HFS is disabled in all cells; BGN and SRH are the
only factors.  All Sentaurus decks export donor concentration, acceptor
concentration, effective intrinsic density, effective band gap, and BGN.

## Observability gates

Both no-BGN logs report `Bennett/Wilson without bandgap narrowing`; the model
name identifies the base formulation while the final phrase confirms that
narrowing is disabled.  In both states the exported `BandgapNarrowing` field is
exactly zero.  Across all four states, 942 silicon nodes have independent donor
and acceptor fields, and the maximum relative difference against `doping.csv`
is zero.

## Terminal-current result

| BGN | SRH | Sentaurus (A/um) | Vela (A/um) | Gap (dex) |
|---|---|---:|---:|---:|
| OldSlotboom | on | `1.12874e-16` | `1.38994e-16` | `0.090399` |
| OldSlotboom | off | `9.58830e-17` | `1.22229e-16` | `0.105432` |
| none | on | `2.63634e-17` | `2.70316e-17` | `0.010871` |
| none | off | `1.93877e-17` | `4.16936e-16` | `1.332544` |

With SRH held on, the BGN response is `0.631594 dex` in Sentaurus and
`0.711122 dex` in Vela.  Their difference is `0.079528 dex`, or 87.97% of the
production-cell `0.090399 dex` gap.  The true no-BGN/SRH-on residual is
`0.010871 dex` (2.535% relative current error).

The previous M29 Bennett/Wilson cell had a smaller `0.002310 dex` gap.  This
does not contradict M40: it was a third physical model, not a no-BGN control.
True no-BGN leaves `0.008560 dex` more residual than that Bennett/Wilson cell.

## Causal interpretation

The SRH-on conditional contrast confirms that BGN-dependent self-consistent
response is the dominant measured contribution at this deep-off point, but it
does not identify a mismatch in the OldSlotboom delta-Eg formula: M39 already
closed that formula and its symmetric band split.  The remaining difference
must arise downstream of delta-Eg, through its coupling to carrier state,
continuity, SRH, or current extraction.

The no-BGN/SRH-off Vela corner has a `1.332544 dex` parity gap and produces a
large BGN-by-SRH interaction.  Therefore the unconditional 2x2 main effects are
not used as a single-cause estimate.  Follow-up work should compare the
OldSlotboom and no-BGN SRH-on states node by node, then localize the conditional
response through electron quasi-Fermi state and continuity residuals.
