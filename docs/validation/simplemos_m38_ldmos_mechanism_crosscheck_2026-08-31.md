# SimpleMOS M38 LDMOS mechanism crosscheck

## Evidence status

M38 reads five existing summary files from the dirty
`codex/templates-ldmos-phase-a` worktree without modifying that worktree or
rerunning Sentaurus.  The source files are recorded by absolute path and
SHA-256.  Because the LDMOS worktree is not committed, this is a qualified
read-only snapshot rather than frozen final LDMOS evidence.

## Cross-device result

| Mechanism | SimpleMOS n23 | Templates/LDMOS G3 | Interpretation |
|---|---:|---:|---|
| Boundary/contact signed measure | 3.52% of all-node measure response | AverageBox reduces frozen electron residual L2 to 1.73% | Geometry matters much more globally in LDMOS, but is not sufficient |
| Direct contact HFS/SG | 0 A/um electron-cut response | Fixed replay improves by 1.0248 dex | Device-specific; not the shared mechanism |
| QF self-consistent feedback | 94.49% adjoint relaxation | 0.4097 of 0.4293 dex total error | Shared dominant pathway |
| Electron-QF specificity | 99.9997% signed electron-equation share | `phin` recovers 0.4082 dex; `psi` only 0.00145 dex | Shared electronic continuity/QF control |

For LDMOS at Vg = 0.5 V, correcting the contact HFS treatment reduces the
fixed-state median error from `1.04443 dex` to `0.0196168 dex`.  The
self-consistent error remains `0.429281 dex`, however, and the state-feedback
increment is `0.409664 dex`.  Replacing only `phin` recovers `0.408183 dex`, or
99.64% of that feedback increment; `psi` and `phip` are negligible by
comparison.

Switching the LDMOS coupling to the signed-cotangent proxy changes operator,
self-consistent, and feedback attribution by less than `9.7e-5 dex`.  A
separate AverageBox audit reduces the frozen electron residual L2 from its
baseline to `1.73%`, yet the candidate same-bias reclose still has
`0.29297 dex` current error.  Thus geometry is a necessary LDMOS correction,
but does not explain the remaining self-consistent current gap.

## Conclusion

The two devices do not share the same direct contact-HFS defect: it is absent
on the SimpleMOS drain cut and material in the current LDMOS implementation.
They do share the more important remaining mechanism after contact/geometry
correction: an internal electron-continuity disturbance is amplified and
transmitted to terminal current through the electron quasi-Fermi state.  The
next LDMOS work should therefore localize the `phin` feedback source and its
Jacobian pathway, while retaining the already-qualified contact and
AverageBox corrections as separate factors.
