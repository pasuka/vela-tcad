# Charon-derived Sentaurus 2018 Schottky reference

WP13 (2026-10-05): the historical voltage/arclength decks explicitly load
`vela/materials_legacy_builtin_300k.json`; the partial mobility ablations load
their corresponding `_legacy_300k.json` snapshots. These retain pre-alignment
effective parameters. Original partial material inputs and archived results
are unchanged. Regenerate/check frozen files with
`python scripts/generate_legacy_material_baselines.py --check`.

This fixture translates Charon's two-dimensional `n_diode` Schottky-contact
test into Sentaurus SDE/SDevice and then compares Vela on the exact imported
Sentaurus mesh. The source benchmark is a 1 um by 1 um, uniformly doped
`1e16 cm^-3` n-silicon slab with a 4.75 eV Schottky anode and an Ohmic
cathode at 300 K.

The Sentaurus O-2018.06-SP2 reference completed the full 0--1 V sweep. The
small committed curve contains the converged equilibrium point and all forward
points; the TDR, PLT, and logs remain under the ignored
`build-release/reference_tcad/schottky_charon_sentaurus2018/raw` tree and are
sealed by hashes in `schottky_charon_sentaurus2018_reference.json`.

Vela's pre-existing `dirichlet_barrier` Gummel model did not converge at 0 V
on the imported 697-node mesh after 150 iterations. The implemented minimum
parity feature is therefore the boundary equation used by the source deck:
electrostatic Schottky barrier pinning plus independent electron/hole
thermionic Robin fluxes in coupled Newton. No image-force lowering, tunnelling,
series resistance, AC, or high-field model is enabled.

The Vela acceptance now covers the full 0.01--1 V Sentaurus range. The
comparison uses all 24 nonzero Sentaurus points, has matching monotonic trend,
and stays within the unchanged 0.5-decade current-shape tolerance. The maximum
error is 0.478652 decade. At 1 V the log-interpolated Vela current is
`1.08262e-4 A/um`, versus `1.18824e-4 A/um` in Sentaurus (8.89% relative
error).

The Vela solve is intentionally split into two numerical stages without
changing the physics. Stage A uses voltage continuation through 0.82 V and
writes restart states. Stage B resumes the accepted 0.82 V state and uses the
existing pseudo-arclength continuation through the high-injection branch to
1 V:

```text
build/vela_example_runner.exe --config reference_tcad/schottky_charon_sentaurus2018/vela/simulation_iv.json
build/vela_example_runner.exe --config reference_tcad/schottky_charon_sentaurus2018/vela/simulation_iv_arclength.json
```

A fresh two-command rerun produces 42 stage-A points and 114 stage-B points,
ending at `1.0002943642565851 V` without voltage backsteps or current
decreases. Stage B explicitly sets `abstol = 2e-9`: the imported 0.82 V state
has already reached the attainable normalized residual floor, while restarting
a purely relative solve would reset that small residual to a relative norm of
one and incorrectly demand another five orders of magnitude reduction. The
freshly merged curve agrees with the checked-in curve to `1.46e-14 V` in bias
and `6.77e-13` relative current.

The clean workflow is automated by
`scripts/run_schottky_reference_workflow.py`. It materializes both stages in
an isolated output directory, verifies the exact 0.82 V restart dependency,
merges the curve, applies the 24-point Sentaurus gate, and audits terminal KCL
at 0, 0.4, and approximately 1 V. The compact passed ledger is
`schottky_workflow_validation_20260817.json`.

The earlier 0.563125 V failure was caused by the deck's
`quasi_fermi_update_limit_V = 0.05` hard update cap. Both analytic and finite-
difference Jacobians failed with that cap, while the same restart transition
converged when the cap was disabled. No new physical model or solver feature
was required. Secant prediction alone was rejected because it could converge
to a non-monotonic branch near 0.835 V.

The decks use `scaling.mode = unit_scaling`, so their thermionic velocity
values are in cm/s despite the stable historical `_m_per_s` field spelling.
The checked-in `vela_schottky_iv_combined.csv` is the small two-column merge of
the two generated stage curves; large restart states and nonlinear histories
remain ignored build artifacts.
