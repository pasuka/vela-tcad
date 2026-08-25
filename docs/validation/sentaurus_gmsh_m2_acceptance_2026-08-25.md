# Sentaurus SDE to Gmsh M2 acceptance

## Scope

M2 extends the public, license-free Gmsh preprocessing path without changing
the M1 frozen subset. The default `--capability-level M2` supports:

- Rectangle bodies with the creation-time ABA or BAB mode replayed over an
  audited OpenCASCADE fragment map;
- conformal region interfaces exported as named physical curves;
- `insert-vertex` on one unambiguous straight authored edge and pick-point
  contacts bound to the resulting complete curve;
- source-ordered NoReplace, LocalReplace and Replace doping semantics;
- line-baseline Gaussian profiles with Positive, Negative or Both direction,
  Eval/NoEval, lateral Gauss/Erf, and zero-decay region/material/eval-window
  restrictions;
- the M1 MSH 4.1, VTU, VTK, generic DF-ISE GRD/DAT, Vela JSON/CSV and strict
  cross-format round-trip contract.

M2 fails closed for AB/ABiA/ABiB/XX, Polygon bodies, body/entity-list/remove/
imprint contacts, refinement controls, nonzero constant-profile decay and
nonzero analytical-restriction decay. Those gaps are not silently treated as
ABA, unrefined, or zero-decay cases.

## Public command

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\sentaurus_gmsh_pipeline.py `
  --sde tests\fixtures\sentaurus_gmsh_m2\m2_core_sde.cmd `
  --output-dir build-release\gmsh-m2 `
  --runtime-dir build-release\gmsh-m2-runtime `
  --mesh-size 0.2 `
  --capability-level M2
```

The input is valid Sentaurus Structure Editor T-2022.03 syntax. In particular,
constant NoReplace placements omit the replacement argument, while a no-decay
LocalReplace placement supplies `0.0 "LocalReplace"`. An analytical restriction
uses `eval-window 0.0 "region"`, in the official positional order.

## Geometry, contact and field invariants

The fixture authors two overlapping 2 x 1 um rectangles under ABA. The expected
owned areas are exactly 1 um2 for `R.P` and 2 um2 for `R.N`. The fragment audit
requires every authored body area to equal the sum of its source-mapped OCC
fragments and independently verifies every membership by boundary-vertex
containment. Ownership is then replayed from body creation order; no OCC tag or
centroid ordering is used.

The inserted vertex at `(0, 0.5)` splits the left boundary. `Anode`, picked at
`(0, 0.25)`, must span exactly `y=0..0.5`; `Cathode` owns the complete selected
right curve. The `R.P+R.N` boundary is an interface physical group, using
official SDE's creation-order naming, and cannot be selected as an exterior
contact because it has two adjacent device surfaces.

At `(2.0, 1.0)`, the east-going baseline's Positive side is upward (its left
normal) and the Gaussian is at `PeakPos=0.2`: donor concentration is
`1.02e18 cm^-3`. The lower side receives only the `2e16 cm^-3` region constant.
The baseline itself also receives only that constant because `NoEval` excludes
it. This orientation is frozen from T-2022.03-SP2 TDR values, correcting the
earlier unverified right-normal assumption. Donor and acceptor remain distinct
arrays; signed net doping is derived only afterward.

At each of the two shared `x=1` interface vertices, official TDR contains one
region-local field record for each side. The public single-valued Gmsh/VTK view
sums those records, matching the existing TDR importer: donor is `2.01e16` and
acceptor is `1e16 cm^-3`. Region-local summation is explicit provenance; no
body tag, centroid, last-owner, or adjacent-cell majority heuristic is used.

## Registered regression gates

CTest `sentaurus_gmsh_m2` contains backend-independent analytic tests and real
Gmsh 4.15.2 integrations. It covers partial BAB overlap, ABA and BAB complete
containment, one ABA body replacing fragments of multiple old bodies, edge
touch, point touch, inserted-edge contact extent, endpoint ambiguity,
NoReplace/LocalReplace/Replace, Gaussian Positive/Negative/Both and NoEval,
strict JSON manifests, and MSH/VTU/VTK/GRD/DAT readback.

The test asserts geometry areas and semantic invariants, not complete Gmsh
connectivity, because OCC tags and triangulation may change across Gmsh
versions. Runtime M0 still pins the accepted implementation to Gmsh 4.15.2.

The final UCRT64 Release build succeeded. The full CTest run passed 708 of 709;
the sole failure was the older PN2D oracle IR golden correctly detecting newly
retained source locations. After updating that hash-sealed golden, its complete
CTest passed on rerun. Thus every one of the 709 registered tests has Release
pass evidence, including the real `sentaurus_gmsh_m2` integration.

## Sentaurus 2022 internal oracle

The licensed-VM runner executes this same owned fixture under Structure Editor
T-2022.03-SP2, enables the evaluated command log, invokes `sde:build-mesh`, and
hash-seals the input, log, version, stdout/stderr, boundary TDR, mesh command and
mesh TDR. Commercial binary artifacts remain under ignored `build-release` and
are not required by public CI. The reviewed text/hash summary is recorded only
after the official run succeeds; a license wait or SDE failure is reported as
such and is never converted into a passing oracle.

Run `20260825_m2_v4` completed with exit code 0 and empty stderr. Official SDE
reported 7 points, 5 triangles and 4 mesh regions (two device regions plus
interface/contact categories). The mesh TDR SHA-256 is
`afc7f365c198bd536290072ed6c5b33f375f92d5586276797b6ac2b2200fe752`;
the evaluated log SHA-256 is
`a141dcd53038778d59586ba500a2300952bf269a6dfb6dcf1c60d90b030b50dd`.
The public, machine-readable golden is
`tests/fixtures/sentaurus_gmsh_m2/sentaurus2022_oracle_summary.json`. It retains
the source and artifact hashes, region areas, contact/interface coordinates,
region-local species values, shared-node aggregate values and all rules derived
from them without requiring a proprietary runtime in CI.
