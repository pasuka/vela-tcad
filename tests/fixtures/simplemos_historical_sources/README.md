# Historical SimpleMOS source provenance

This archive verifies the implementation records attached to M9/M10/M12 and
M30-M46 historical evidence. It does **not** qualify the current solver or
assert that the old numerical results are reproduced by current code.

The original evidence files, numerical assertions, failed hypotheses and
artifact hashes remain unchanged. `manifest.json` pins the evidence bytes and
each original source mapping. `sources.zip` contains 132 recovered path/hash
pairs, stored under their original SHA-256 digests. Exact bytes are recovered
from the recorded Git blobs; `materialization` records whether the historical
checkout used blob bytes, LF or CRLF. Tests hash archive bytes without newline
normalization and need neither a Git checkout nor the original source paths.

One source revision was not recovered: `CoupledDDAssembler.cpp`, SHA-256
`e525d828d0199dafdea947377dc6faed7cea419b5cb76fc1f51ae05d3dfe64ed`, referenced by
M33/M35/M37. Its original metadata is preserved and checked, but its source bytes
remain unavailable. The dedicated regression asserts this exact gap; unknown
gaps fail. Passing archive-integrity tests must not be described as recovery of
that source or complete historical reproducibility.

Current implementation qualification is separate: C++ numerical/property tests,
state-interface and workflow regressions, and the independently sealed 816-point
engineering campaign described in
`docs/validation/simplemos_engineering_joined_closeout_2026-09-30.md`.
