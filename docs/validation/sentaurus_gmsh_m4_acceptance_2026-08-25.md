# Sentaurus SDE to Gmsh M4 acceptance - 2026-08-25

> Historical snapshot: PN2D BV status was superseded by the M5 acceptance record
> (`docs/validation/sentaurus_gmsh_m5_acceptance_2026-08-26.md`). New main sync and
> reference-script coverage are tracked in the 2026-09-05 reference SDE plan.


## Outcome

The public M4 Gmsh-to-Vela execution path is implemented and exercised with
the complete PN2D IV and BV decks. IV is accepted. BV reaches the contracted
-20 V endpoint but remains failed because its mixed terminal-current
conservation gate is 4.3% over the limit. The internal TDX gate is implemented
and its hash-sealed dry-run and licensed live round-trip both pass. The returned
TDR preserves the complete Gmsh topology, contacts and all seven node datasets.
Overall M0-M4 implementation status is complete; the PN2D BV conservation
failure remains an explicitly reported solver acceptance issue.

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

The dry-run manifest passed. After explicit authorization, the live gate ran on
2026-08-26 against Sentaurus Data Explorer T-2022.03-SP2. TDX cold start took
several minutes, then conversion completed with exit code 0 and empty stderr.

The final sealed run is
`build-release/m4-tdx-live-final-20260826/tdx_gate_manifest.json` with remote
directory
`~/sentaurus_runs/vela_gmsh_tdx/pn2d_gmsh_m4_20260826_final`. Its evidence is:

- input GRD SHA-256:
  `f9be73280a55b104a96e099aeb101f6033882002a770c1d01dea9ed42641b230`;
- input DAT SHA-256:
  `b2c96cb0b95f5533f95e2f467e5284cc1c2cbe24be93607ce1d84c458ef7dde5`;
- returned TDR SHA-256:
  `6cf588ff06919f295cfbb010f4d3589009aa634848c54c598aa2b197ea43d7d4`;
- frozen `tdx -h` SHA-256:
  `f50182d90e3120999661a6588880189348e20386bab7fcb36dd00b8f919bb815`;
- 1062 vertices, 2010 total elements and 1990 bulk Tri3 elements match;
- region names `R.Si`, `Anode` and `Cathode` match;
- Anode and Cathode each retain 10 contact edges;
- all seven datasets compare at all 1062 nodes with zero maximum absolute and
  relative difference.

## Verification status

- M4 unit/integration regression: 5/5 passed.
- M4-affected Sentaurus/Gmsh Python regression: 114/114 passed.
- UCRT64 Release configure/build: passed; no C++ relink was required.
- UCRT64 Release CTest: 711/711 passed in 157.26 seconds.
- Complete Gmsh PN2D IV: accepted.
- Complete Gmsh PN2D BV: reached -20 V, acceptance failed only on terminal
  current conservation.
- Internal TDX dry-run: passed.
- Internal TDX live round-trip: passed on T-2022.03-SP2; conversion exit code 0,
  strict topology/contact/dataset comparison passed.

The returned TDR, importer export, logs and machine-readable comparison report
are retained under `build-release/m4-tdx-live-final-20260826/`.
