# PN2D CI reference inputs

`pn2d-ci-reference-v1.tar.gz` contains the GRD/DAT inputs, CMD decks, parameter
files, numerical sidecars and numerical references used in the Vela VM
qualification of commit `86ec6199`. The archive includes these inputs:

- `base/{0v,iv,bv}` holds the original 1,943-node mesh and reference data
  derived from the [v2 registry](../artifact_manifest_v2.json).
- `variants/<id>/{forward,reverse}` holds the 19,511-node M3 cases from the
  [variants campaign](../variants/README.md). Each branch includes its native
  PLT and exported field CSVs plus the independently solved Vela JSON terminal
  reference.

The Sentaurus numerical reference is T-2022.03-SP2. Fields retained in the
archive are electrostatic potential, electron and hole density, and electric
field. The archive stores text inputs and references so the regression can run
without a Sentaurus installation or access to either VM. Generated candidate
results stay in the build directory and CI report artifacts.

The runner checks the archive SHA-256 before extraction and checks CMD tokens
against the checked-in source or the M3 deck generator. If the baseline changes,
refresh the archive from qualified reference outputs and update its pinned
SHA-256 in [the runner](../../../scripts/run_pn2d_ci.py) in the same change.
