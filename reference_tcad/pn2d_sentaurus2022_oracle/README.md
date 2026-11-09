# PN2D Sentaurus 2022 SDE oracle

The owned source and wrapper, frozen Device IR, generated mesh command and
hash manifest describe the historical T-2022.03-SP2 validation run (1943 points,
3680 elements, 3 regions). These small reviewable fixtures remain in Git.

TDR files, evaluated logs, stdout/stderr, version output and exit-code files
are local validation evidence only. They are ignored and untracked. The
manifest retains their original SHA-256 and byte counts in `files`; entries
listed in `local_validation_files` are not required by a public checkout.
The public manifest contains only a concise version identifier, not raw output.

Public tests always compare the owned SDE source with the frozen IR semantics
and verify hashes of distributed files. Raw evidence is hash-checked whenever
present; log replay is explicitly skipped when absent. For a required internal
oracle audit, restore the original bytes under `artifacts/` and set
`VELA_REQUIRE_LOCAL_SDE_ORACLE=1`; missing evidence then fails the test.
Absence of local evidence never establishes a fresh live Sentaurus/TDX pass.

During the 2026-09-05 cleanup, original bytes and the original manifest were
also preserved under ignored `build-release/private-sde-oracle-20260905/`.
Historical Git objects were not rewritten; a separate history-cleaning decision
is needed if publishing must exclude those objects as well as the current tree.
Do not use `git add -f` for validation evidence or bulk-renormalize sealed files.
