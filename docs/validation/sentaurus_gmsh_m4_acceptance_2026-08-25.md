# Sentaurus SDE to Gmsh M4 acceptance - 2026-08-25

## Outcome

The public M4 Gmsh-to-Vela execution path is implemented and exercised with
the complete PN2D IV and BV decks. IV is accepted. BV reaches the contracted
-20 V endpoint but remains failed because its mixed terminal-current
conservation gate is 4.3% over the limit. The internal TDX gate is implemented
and its hash-sealed dry-run passes, but a live run is pending explicit approval
to transfer the generated GRD/DAT device data to the configured `sentaurus`
SSH target. Overall M0-M4 status therefore remains incomplete.

## Public orchestration contract

- `sentaurus_import.py device` defaults to `offline-static + gmsh`.
- The M3 Device IR is passed directly to the Gmsh backend; an evaluated frontend
  is never silently reparsed as offline source.
- `tensor-legacy` is explicit compatibility mode.
- `official-tdr` fails closed unless paired with `sentaurus2022-evaluated` and
  an actual `_msh.tdr` result.
- The unified manifest retains the Gmsh runtime, cross-format, mesh-quality and
  source hashes and distinguishes generated, solved, accepted and failed.
- A Gmsh BAMG mesh is qualified as an explicitly obtuse-capable barycentric
  mesh. It is not mislabeled as a non-obtuse mixed-Voronoi mesh.

## Complete PN2D IV evidence

The complete checked-in PN2D SDE deck and IV SDevice deck were run with Gmsh
4.15.2 and the UCRT64 Release solver:

- mesh: 1062 nodes, 1990 Tri3 cells and 11 nodes on each terminal;
- maximum triangle angle: 134.64 degrees;
- solver policy: `pn2d_sentaurus2018_iv_qualified` with barycentric volumes;
- curve: 16 accepted points, 0.0 through 0.3 V;
- maximum relative terminal imbalance: 0.001656;
- maximum absolute imbalance: 1.103e-17 A/um;
- maximum mixed tolerance ratio: 0.008119, below the required 1.0;
- final status: `accepted`.

## Complete PN2D BV evidence

The same public command path was run with the complete BV SDevice deck:

- mesh: 1038 nodes and 1942 Tri3 cells;
- maximum triangle angle: 126.01 degrees;
- curve: 406 converged points and final accepted bias -20 V;
- maximum relative terminal imbalance: 0.08018;
- maximum absolute imbalance: 5.966e-18 A/um;
- maximum mixed tolerance ratio: 1.04299 versus the required maximum 1.0;
- final status: `failed`, stage `acceptance`;
- exact failure: terminal-current imbalance exceeds the mixed absolute/relative
  contract maximum.

The runner did not fail to converge and did not truncate the requested voltage.
The remaining BV discrepancy is a terminal-current closure issue on the Gmsh
barycentric mesh and is not hidden by relaxing the contract.

## Internal TDX gate

`scripts/sentaurus_gmsh_tdx_gate.py` validates the preprocessing manifest and
GRD/DAT hashes before any transfer. A live run is designed to:

1. capture and hash `tdx -h`;
2. stage only `device.grd/device.dat` in an isolated VM directory;
3. execute `tdx -d device.grd device.dat device.tdr`;
4. retrieve the TDR and logs;
5. import the TDR with the C++ importer;
6. compare nodes, Tri3 topology, regions, contacts and datasets in a strict JSON
   report.

The dry-run manifest passed. Live execution was blocked before any network
transfer because the `sentaurus` SSH alias was not explicitly approved as a
trusted destination for device geometry and doping data. This is an external
authorization gate, not a public-product dependency or a successful TDX claim.

## Verification status

- M4 unit/integration regression: 5/5 passed.
- M4-affected Sentaurus/Gmsh Python regression: 114/114 passed.
- UCRT64 Release configure/build: passed; no C++ relink was required.
- UCRT64 Release CTest: 711/711 passed in 157.26 seconds.
- Complete Gmsh PN2D IV: accepted.
- Complete Gmsh PN2D BV: reached -20 V, acceptance failed only on terminal
  current conservation.
- Internal TDX dry-run: passed.
- Internal TDX live round-trip: pending explicit destination/data-transfer
  approval.

No live TDX result may be added without the returned TDR and machine-readable
comparison report.
