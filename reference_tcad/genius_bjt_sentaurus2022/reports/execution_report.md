# Genius BJT Sentaurus WP0-WP2 execution report

Run date: 2026-09-02 (Asia/Shanghai)

## Outcome

WP0, WP1, and WP2 pass. The accepted structure is an exact SDE evaluation of
the Genius BJT analytical doping formula, not a fitted Gaussian approximation.
Both SDevice model levels completed the VBE=0.70 V preparation and the
VCE=0..3 V sweep at 0.1 V reporting intervals.

## Environment

- Local worktree branch: `codex/genius-bjt-sentaurus2022`
- Local base commit: `270f5258a3c6354b2b3564af98ae4ab9434de2ca`
- Genius source commit: `543da8452d5dfd33e6f8c457f962f6f670f0fce7`
- Remote host/user: `tcad` / `root`
- SDE and SDevice: T-2022.03-SP2
- Accepted remote run: `/root/sentaurus_runs/vela_oracle_2022/genius_bjt_wp0_wp2_20260902/run_v2`

The actual remote account is `root`; this is intentionally recorded because an
older workflow note expected `tcad` as the login user.

## WP1 structure gate

- Mesh bounds: x=[0,6] um, y=[0,2] um
- Nodes/elements: 5611 / 10940
- Base contact: x=[1.25,2.00] um, 18 nodes
- Emitter contact: x=[2.75,4.25] um, 34 nodes
- Collector contact: x=[0,6] um at y=2 um, 65 nodes
- Maximum contact coordinate error: 0 um
- Maximum donor relative error against the Genius formula: 3.60e-15
- Maximum acceptor relative error against the Genius formula: 3.62e-15

The first diagnostic mesh was rejected before use because coordinate-dependent
terms had incorrectly been placed in the Sentaurus `General(init=...)` string.
Sentaurus evaluates `init` once, which made the resulting doping spatially
constant and unphysical. In the accepted run, `init` contains constants only
and every x/y-dependent expression is evaluated in `function`.

## WP2 device gate

Both models completed without `does not converge`, `Linear system cannot`,
`NAN`, or syntax-error messages. PLT files contain exactly 31 requested rows.
Currents below are direct Sentaurus `TotalCurrent` values in A/um.

| VCE (V) | M0 Ic | M0 Ib | M0 beta | M1 Ic | M1 Ib | M1 beta |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | -8.4941e-8 | 1.0135e-7 | 0.838 | -2.9925e-7 | 3.6010e-7 | 0.831 |
| 0.1 | 1.6548e-6 | 1.7971e-8 | 92.080 | 2.4776e-6 | 6.7090e-8 | 36.929 |
| 1.0 | 1.7607e-6 | 1.6085e-8 | 109.462 | 2.6589e-6 | 6.0378e-8 | 44.038 |
| 3.0 | 1.8433e-6 | 1.6039e-8 | 114.928 | 2.8075e-6 | 6.0188e-8 | 46.646 |

- M0 maximum absolute/relative KCL residual: 1.61e-20 A/um / 8.93e-15
- M1 maximum absolute/relative KCL residual: 4.63e-20 A/um / 1.69e-14
- Maximum requested-voltage selection error: 0 V for both models

## Scope boundary

No Vela solver, importer, or BJT configuration was changed in WP0-WP2. The two
normalized reference curves and four representative state TDRs per model are
the inputs for the later WP3+ Vela conversion and cross-solver comparison.

